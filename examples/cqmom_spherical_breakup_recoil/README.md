# 4D spherical-particle breakup with recoil

This standalone pair compares CQMOM with a direct Monte Carlo reference for a
homogeneous spherical-particle population in phase space `(d,u,v,w)`. Both use
Schiller-Naumann drag in a constant carrier velocity, equal-volume binary
breakup, and the same deterministic sampled initial population.

The CQMOM source deliberately evaluates every parent-loss/two-daughter-birth
contribution through `CqmomSourceContext::for_each_terminal_node(...)`. The
daughter pair has equal and opposite Gaussian recoil components, so its
conditional expected raw moments are used in the terminal-quadrature source.
The MC companion samples the same symmetric recoil pair at each event.

`configuration.qmom` controls the fixed hierarchy `(2,2,2,2)`, conditioning
order, solver policy, Strang splitting, node-advection tolerances, and output.
The application constants set a short `0.02 s` horizon and `0.001 s` macro
step. The smaller-than-initially-targeted diameter count keeps the fixed-branch
source evolution realizable for this lightweight sampled case.

Both executables regenerate the same 16,000-particle ensemble using seed
`2026082301`; they print the maximum difference between the CQMOM input tensor
and moments recomputed from that ensemble. The MC breakup stream uses seed
`2026082302`. Reported diagnostics are number density, total particle volume,
volume-weighted momentum, and volume-weighted `u^2`, which exposes recoil.

Build against installed QMOMKit and QMOMKitMonteCarlo:

```bash
cmake -S . -B build -DCMAKE_PREFIX_PATH=/path/to/install
cmake --build build
cd build
./cqmom_spherical_breakup_recoil_example
./cqmom_spherical_breakup_recoil_mc_example
```

The CQMOM moment/quadrature histories and the MC diagnostic history are written
under `output/` when enabled in the local configuration.
