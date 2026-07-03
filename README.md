# benchmark-simple-potato

Benchmark of collisionless guiding-centre orbits in axisymmetric tokamak
fields: [SIMPLE](https://github.com/itpplasma/SIMPLE) (symplectic, canonical
flux coordinates) against [POTATO](https://github.com/itpplasma/NEO-RT)
(adaptive RK, cylindrical coordinates, finite orbit width), with NEO-RT
thin-orbit frequencies as third reference.

- `doc/benchmark.tex`: physics, definitions, metrics, literature
  (build: `latexmk -pdf benchmark.tex` in `doc/`).
- `WORKPLAN.md`: execution ladder, current status, task split.
- `rung0/`: common circular-tokamak field inputs for both codes.
