import importlib.util
import json
from pathlib import Path

import numpy as np
import pytest


SCRIPT = Path(__file__).parents[1] / "tools" / "potato_invariants_to_simple.py"
SPEC = importlib.util.spec_from_file_location("potato_invariants_to_simple", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


HEADER = """# potato-simple-invariants-v1
# id R_cm Z_cm p xi sigma H0 J_perp psi_star psi_axis psi_edge phi_elec v0_cm_s status
"""


def write_handoff(path: Path, *, h0=1.0, phi=0.0, v0=2.0):
    path.write_text(
        HEADER
        + f"1 180 0 1 0.3 1 {h0} 0.00005 4558211 0 15693266 {phi} {v0} 0\n"
    )


def fake_resolver(h0, j_perp, psi_star, **kwargs):
    assert (h0, j_perp, psi_star) == (1.0, 0.00005, 4558211.0)
    assert kwargs["v0_ratio"] == 0.5
    return {
        "states": np.array([[0.2, 0.0, 0.0, 0.5, 0.3]]),
        "cylindrical": np.array([[180.0, 0.0, 0.0]]),
        "sigma": np.array([1]),
        "section_branch": np.array([2]),
        "section_kind_code": np.array([1]),
        "section_kind": np.array(["B_min"]),
        "residual": np.array([1.0e-12]),
        "n_solutions": 1,
        "status": 0,
        "ordering": "POTATO R_c (HFS to LFS)",
        "simple_invariants": np.array([0.25, 0.0000125, 4558211.0]),
    }


def fake_inspector(state):
    assert np.array_equal(state, np.array([0.2, 0.0, 0.0, 0.5, 0.3]))
    return {"h0": 0.25, "j_perp": 0.0000125, "p_phi": 4558211.0, "status": 0}


def test_parse_and_explicit_velocity_rescaling(tmp_path):
    handoff = tmp_path / "potato_invariants.dat"
    write_handoff(handoff)

    rows = MODULE.read_handoff(handoff)
    result = MODULE.convert_rows(
        rows,
        fake_resolver,
        simple_v0_cm_s=4.0,
        allow_v0_rescale=True,
        inspector=fake_inspector,
    )

    assert result["schema"] == "potato-simple-candidates-v1"
    assert result["time_convention"] == "tau=v0*t"
    assert result["rows"][0]["source"]["p_phi_convention"] == (
        "raw POTATO psi_star=(c/q)P_phi"
    )
    assert result["rows"][0]["candidates"][0]["state"][3:] == [0.5, 0.3]
    assert result["rows"][0]["candidates"][0]["reconstructed_invariants"]["h0"] == 0.25
    json.dumps(result)


def test_velocity_mismatch_is_rejected_by_default(tmp_path):
    handoff = tmp_path / "potato_invariants.dat"
    write_handoff(handoff)
    rows = MODULE.read_handoff(handoff)

    with pytest.raises(ValueError, match="reference velocity"):
        MODULE.convert_rows(rows, fake_resolver, simple_v0_cm_s=4.0)


@pytest.mark.parametrize(
    ("changes", "message"),
    [
        ({"h0": 1.1}, "H0"),
        ({"phi": 0.01}, "potential"),
        ({"v0": -1.0}, "v0"),
    ],
)
def test_normalization_errors_are_rejected(tmp_path, changes, message):
    handoff = tmp_path / "potato_invariants.dat"
    write_handoff(handoff, **changes)

    with pytest.raises(ValueError, match=message):
        MODULE.read_handoff(handoff)
