# 5D thermal spherical-particle breakup

This standalone pair compares CQMOM with direct Monte Carlo for the homogeneous
phase space `(V,T,u,v,w)`. It shares the physical theme of the 4D recoil example:
Schiller-Naumann drag, a constant carrier velocity, and smooth equal-volume
binary breakup. This case adds Ranz-Marshall lumped-capacitance heat transfer.

Its intentionally different closure is the key feature. The CQMOM breakup
source does not visit terminal nodes. It evaluates the analytical
`(2^(1-p)-1)` breakup expression with calls such as
`context.generalized_moment({p + 1.0/3.0, q, r, s, t})`. These are derived
stage-local fractional moments required because `d(V)` is proportional to
`V^(1/3)`. Daughters have half the parent volume and inherit temperature and
velocity exactly; there is no recoil.

`configuration.qmom` controls the fixed `(2,2,2,2,2)` hierarchy, conditioning
order, solver policy, Strang splitting, node-advection tolerances, and output.
The short `0.02 s` horizon and `0.001 s` macro step are application-local.
The two-node volume direction is a deliberate stability compromise for the
fixed-branch sampled demonstration.

Both executables regenerate the same 16,000-particle ensemble using seed
`2026082311`; each prints the difference between its CQMOM input tensor and
the MC moments recomputed from that ensemble. The MC breakup stream uses seed
`2026082312`. Diagnostics include number density, volume, volume-weighted
temperature, volume-weighted momentum, and volume-weighted `u^2`.

Build against installed QMOMKit and QMOMKitMonteCarlo:

```bash
cmake -S . -B build -DCMAKE_PREFIX_PATH=/path/to/install
cmake --build build
cd build
./cqmom_thermal_spherical_breakup_example
./cqmom_thermal_spherical_breakup_mc_example
```

The CQMOM moment/quadrature histories and the MC diagnostic history are written
under `output/` when enabled in the local configuration.
