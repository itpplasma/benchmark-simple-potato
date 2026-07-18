import importlib.util
import math
import re
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pytest


REPO = Path(__file__).parents[1]
POTATO = REPO.parent / "NEO-RT" / "POTATO" / "build" / "potato.x"
SCRIPT = REPO / "tools" / "potato_invariants_to_simple.py"
SPEC = importlib.util.spec_from_file_location("potato_invariants_to_simple", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


@pytest.mark.skipif(not POTATO.exists(), reason="sibling POTATO Release build is absent")
def test_reactor_case_round_trip(tmp_path):
    pysimple = pytest.importorskip("pysimple")
    source = REPO / "rung0" / "potato_run"
    for name in (
        "circ.eqdsk",
        "convexwall.dat",
        "field_divB0.inp",
        "potato.in",
        "profile_poly.in",
    ):
        shutil.copy(source / name, tmp_path / name)

    potato_run = subprocess.run(
        [POTATO],
        cwd=tmp_path,
        check=True,
        timeout=120,
        text=True,
        capture_output=True,
    )
    rows = MODULE.read_handoff(tmp_path / "potato_invariants.dat")
    pysimple.init(
        str(REPO / "rung0" / "circ_chartmap_simple.nc"),
        deterministic=True,
        ntestpart=1,
        npoiper2=1024,
        trace_time=1.0e-3,
        n_e=1,
        n_d=2,
        facE_al=700.0,
    )
    converted = MODULE.convert_rows(
        rows,
        pysimple.states_from_potato_invariants,
        simple_v0_cm_s=float(pysimple.params.v0),
        inspector=pysimple.invariants_from_state,
    )

    result = converted["rows"][0]
    assert np.isclose(result["v0_ratio_potato_over_simple"], 1.0, rtol=1e-12)
    assert len(result["candidates"]) >= 1
    source_rz = np.array([rows[0]["R_cm"], rows[0]["Z_cm"]])
    nearest = min(
        result["candidates"],
        key=lambda item: np.linalg.norm(
            np.array(item["cylindrical_native"])[[0, 2]] - source_rz
        ),
    )
    nearest_rz = np.array(nearest["cylindrical_native"])[[0, 2]]
    assert np.linalg.norm(nearest_rz - source_rz) < 0.1
    assert abs(nearest["state"][4] - rows[0]["xi"]) < 1.0e-3

    match = re.search(
        r"single orbit: taub=\s*([+\-0-9.E]+)\s+delphi=\s*([+\-0-9.E]+)\s+ierr=0",
        potato_run.stdout,
    )
    assert match is not None
    taub, delta_phi = map(float, match.groups())
    potato_omega_b = 2.0 * math.pi * rows[0]["v0_cm_s"] / taub
    potato_omega_phi = delta_phi * rows[0]["v0_cm_s"] / taub
    simple_frequency = pysimple.compute_canonical_frequencies(
        np.array(nearest["state"]), integrator="midpoint", n_periods=1
    )
    assert simple_frequency["status"] == 0
    assert math.isclose(simple_frequency["omega_b"], potato_omega_b, rel_tol=0.02)
    assert math.isclose(
        simple_frequency["omega_phi"], potato_omega_phi, rel_tol=0.02
    )
