# Rust MD Workstation

A Rust-native, CPU-focused molecular dynamics workstation for education,
prototyping, small Lennard-Jones simulations, trajectory analysis, and
reproducible reports.

This project starts small on purpose. The first milestone is a CLI-driven
Lennard-Jones engine with Velocity Verlet integration, steepest-descent energy
minimization, XYZ trajectory output, CSV energy logging, simple XYZ/PDB
molecular input, harmonic bond/angle/Coulomb force-field terms, non-bonded
exclusions, run manifests, and Markdown reports.

## What This Is

- A native Rust workspace with a built-in CPU simulation engine.
- A CLI-first tool for small-scale molecular dynamics experiments.
- A clean base for future force fields, analysis tools, and reports.
- An educational and extensible codebase.

## What This Is Not

- Not a replacement for GROMACS, AMBER, NAMD, OpenMM, or Desmond.
- Not a clinically validated drug design engine.
- Not a production pharmaceutical simulation package.
- Not dependent on Python, Conda, Docker, WSL2, GPU acceleration, or external MD engines.

## Quick Start

```bash
cargo run -p md-cli -- run examples/lj-fluid.toml
```

The command writes a run directory like:

```text
runs/lj-fluid-001/
  config.toml
  trajectory.xyz
  energy.csv
  summary.json
  run-manifest.json
  run-report.md
  checkpoint.json
```

Validate a config without running a simulation:

```bash
cargo run -p md-cli -- validate examples/lj-fluid.toml
```

Run tests:

```bash
cargo test
```

## Workspace Layout

```text
crates/
  md-core/        Core state, units, energies, and validation
  md-force/       Lennard-Jones force and potential calculation
  md-integrator/  Velocity Verlet integration
  md-neighbor/    Cell-list neighbor search and rebuild checks
  md-analysis/    Energy, RMSD, and distance analysis
  md-report/      Markdown run report generation
  md-io/          XYZ, CSV, and JSON output helpers
  md-cli/         Command-line workflow
  md-ui/          Optional desktop UI prototype
examples/
  lj-fluid.toml              First runnable Lennard-Jones example
  molecules/argon-cluster.xyz
  molecules/argon-dimer.xyz
  molecules/argon-angle.xyz
  molecules/argon-force-field.xyz
  molecules/argon-dihedral.xyz
  molecules/argon-minimize.xyz
  molecules/mixed-argon-neon.xyz
  molecules/alanine-fragment.pdb
  topology/argon-dimer.toml
  topology/argon-angle.toml
  topology/argon-force-field.toml
  topology/argon-dihedral.toml
  validation/small-nve.toml  Small deterministic NVE check
  validation/long-nve-drift.toml
  validation/cold-lattice.toml
  validation/pbc-small-nve.toml
  validation/neighbor-compare.toml
  validation/larger-neighbor-compare.toml
  validation/topology-neighbor-compare.toml
  validation/parallel-open.toml
  validation/xyz-input.toml
  validation/bonded-dimer.toml
  validation/angle-trimer.toml
  validation/force-field-report.toml
  validation/dihedral-chain.toml
  validation/berendsen-thermostat.toml
  validation/minimize-lj.toml
  validation/pdb-input.toml
  validation/binary-checkpoint.toml
  validation/force-field-types.toml
```

## Desktop UI Prototype

The CLI remains the primary workflow. An optional egui/eframe desktop prototype
is isolated in `md-ui` and built only when the `desktop` feature is enabled:

```bash
cargo run -p md-ui --features desktop
```

The prototype can open and edit validation TOML configs, start a simulation via
the Rust CLI workflow, list run directories, preview `energy.csv`, and display
`run-report.md`. It does not call external MD engines.

## Units

The current engine uses reduced Lennard-Jones units:

- mass is dimensionless
- sigma defines the length scale
- epsilon defines the energy scale
- Boltzmann constant is treated as 1.0 for temperature estimates
- the current temperature estimate uses 3N degrees of freedom and does not yet
  subtract center-of-mass or constraint degrees of freedom

These assumptions are intentionally simple and are documented here so the first
engine remains honest and reproducible.

## Validation

The validation suite adds deterministic checks around the first engine:

- Lennard-Jones force symmetry, cutoff behavior, overlap rejection, and known
  potential values.
- Velocity Verlet stability for simple NVE cases.
- A 10k-step NVE energy-drift regression with a shifted LJ potential and small
  timestep.
- Deterministic seed handling for initial velocities.
- Regression coverage for repeated CLI runs with identical energy logs.
- Config validation for unsupported force types and invalid initial lattices.
- Periodic minimum image behavior, cutoff constraints, and shifted-potential
  behavior.
- Neighbor-list pair generation, rebuild checks, and naive-vs-neighbor force
  equivalence.
- Rayon-based parallel force evaluation checked against serial force
  evaluation.
- Basic Berendsen thermostat temperature coupling with deterministic metadata.
- Steepest-descent energy minimization with finite-state rejection and
  before/after energy history.
- Energy summaries, RMSD series, and pair-distance analysis over run outputs.
- Standard XYZ molecule input, a small PDB coordinate subset, config-relative
  input paths, atom element labels in generated trajectories, and PDB atom
  metadata preservation in run outputs.
- Simple topology files with per-atom charges, harmonic bonds, harmonic angles,
  bonded non-bonded exclusions, optional Coulomb interactions, and combined
  LJ + bonded force integration.
- Per-element force-field defaults for mass, charge, sigma, and epsilon, with
  Lorentz-Berthelot LJ mixing when multiple LJ types are present.
- Product-layer run manifests, reproducibility metadata, project metadata,
  Markdown run reports, and JSON or compact binary checkpoint/resume support.

Run the validation examples:

```bash
cargo run -p md-cli -- validate-suite
cargo run -p md-cli -- validate-suite --long
cargo run -p md-cli -- run examples/validation/small-nve.toml
cargo run -p md-cli -- run examples/validation/long-nve-drift.toml
cargo run -p md-cli -- run examples/validation/cold-lattice.toml
cargo run -p md-cli -- run examples/validation/pbc-small-nve.toml
cargo run -p md-cli -- run examples/validation/neighbor-compare.toml
cargo run -p md-cli -- run examples/validation/topology-neighbor-compare.toml
cargo run -p md-cli -- run examples/validation/parallel-open.toml
cargo run -p md-cli -- run examples/validation/xyz-input.toml
cargo run -p md-cli -- run examples/validation/bonded-dimer.toml
cargo run -p md-cli -- run examples/validation/angle-trimer.toml
cargo run -p md-cli -- run examples/validation/force-field-report.toml
cargo run -p md-cli -- run examples/validation/dihedral-chain.toml
cargo run -p md-cli -- run examples/validation/berendsen-thermostat.toml
cargo run -p md-cli -- minimize examples/validation/minimize-lj.toml
cargo run -p md-cli -- run examples/validation/pdb-input.toml
cargo run -p md-cli -- run examples/validation/binary-checkpoint.toml
cargo run -p md-cli -- run examples/validation/force-field-types.toml
```

`validate-suite` performs a compact config, input, topology, neighbor, and force
setup smoke check for the curated validation configs. The default suite skips
longer stress configs; pass `--long` to include `long-nve-drift.toml` and the
larger neighbor benchmark config. The manual `long-nve-drift-100k.toml` file is
kept out of the suite.

For a manual longer stress check, copy `examples/validation/long-nve-drift.toml`
and raise `simulation.steps` to `100000`; keep `dt = 0.0001` unless you are
intentionally exploring timestep sensitivity.

## Periodic Boundaries

Phase 3 supports orthorhombic periodic boundary conditions:

```toml
[box]
x = 10.0
y = 10.0
z = 10.0
periodic = true

[force]
type = "lennard-jones"
epsilon = 1.0
sigma = 1.0
cutoff = 2.5
shift_potential = true
```

When `periodic = true`, pair distances use the minimum image convention and
positions are wrapped back into the simulation box after each integration step.
The cutoff must be no larger than half the shortest box length. `shift_potential`
subtracts the LJ potential at the cutoff so the reported potential energy is
continuous at the cutoff; forces are not force-shifted yet.

## Neighbor Lists

Phase 4 adds a cell-list based Verlet neighbor list:

```toml
[neighbor]
enabled = true
skin = 0.3
rebuild_interval = 10
```

The neighbor list is built using `cutoff + skin`. During integration it rebuilds
when any particle has moved more than half the skin distance since the last
build, or when `rebuild_interval` requests a forced rebuild. For periodic boxes,
`cutoff + skin` must fit within half the shortest box length.

Compare naive and neighbor-list force evaluation:

```bash
cargo run -p md-cli -- bench-neighbor examples/validation/neighbor-compare.toml --repeats 50
```

For a larger 512-particle workload:

```bash
cargo run -p md-cli -- bench-neighbor examples/validation/larger-neighbor-compare.toml --repeats 10
```

The benchmark path also honors topology, Coulomb, and non-bonded exclusions:

```bash
cargo run -p md-cli -- bench-neighbor examples/validation/topology-neighbor-compare.toml --repeats 50
```

## Parallel CPU Execution

Phase 5 adds Rayon-based CPU parallel force evaluation:

```toml
[execution]
parallel = true
```

The parallel force path uses bounded local force buffers and then reduces them
back into the system state, avoiding shared mutable force accumulation. Naive
parallel pair evaluation streams over the upper-triangular pair space instead of
allocating the full pair list. It works with both naive pair evaluation and
neighbor-list evaluation. Rayon chooses the thread count by default; set
`RAYON_NUM_THREADS` in the environment to control it for experiments.

The neighbor benchmark also reports serial and parallel timings:

```bash
RAYON_NUM_THREADS=4 cargo run -p md-cli -- bench-neighbor examples/validation/larger-neighbor-compare.toml --repeats 10
```

## Thermostats

Runs default to NVE with no thermostat:

```toml
[thermostat]
type = "none"
```

A simple Berendsen thermostat can be enabled for deterministic NVT-style
temperature coupling:

```toml
[thermostat]
type = "berendsen"
target_temperature = 0.2
tau = 0.02
```

If `target_temperature` is omitted, `[simulation].temperature` is used as the
target. The implementation rescales velocities after each integration step with
the standard Berendsen first-order coupling factor. This is useful for gentle
temperature control in small examples, but it is not a rigorous canonical
ensemble sampler and should not be used for production thermodynamic sampling.

Run the thermostat validation example:

```bash
cargo run -p md-cli -- run examples/validation/berendsen-thermostat.toml
```

## Energy Minimization

The CLI can relax an initial structure with a conservative steepest-descent
workflow:

```bash
cargo run -p md-cli -- minimize examples/validation/minimize-lj.toml
```

Minimization uses the same configured Lennard-Jones, Coulomb, topology, boundary,
parallel force, and shifted-potential options as a dynamics run. Velocities are
set to zero, trial steps move along the force direction, `max_displacement`
limits the largest per-step atomic displacement, and `max_backtracks` halves a
trial step until it lowers the potential energy.

```toml
[minimization]
steps = 80
output_interval = 10
step_size = 0.001
max_displacement = 0.02
force_tolerance = 0.000001
max_backtracks = 12
```

The command writes a normal run directory plus minimization-specific artifacts:

```text
trajectory.xyz
energy.csv
minimization.csv
minimized.xyz
summary.json
run-manifest.json
run-report.md
```

`minimization.csv` records `step,potential,max_force,step_scale` for every
accepted step. This is a simple relaxation tool, not a conjugate-gradient or
production-quality optimizer.

## Analysis

Phase 6 adds run-directory analysis:

```bash
cargo run -p md-cli -- analyze runs/lj-fluid-001 --reference-frame 0 --pair 0,1
```

The command reads `energy.csv` and `trajectory.xyz`, then writes:

```text
analysis-summary.json
analysis-rmsd.csv
analysis-distance.csv
analysis-report.md
```

Current analysis tools include energy drift and temperature summaries, RMSD
against a reference frame, and distance series for one atom pair. RMSD is
computed after centering and Kabsch-style rigid alignment, so global translation
and rotation do not dominate the reported value.

When `config.toml` is present in the run directory and `[box].periodic = true`,
analysis first unwraps XYZ frames with frame-to-frame minimum-image
displacements using the configured orthorhombic box. This reduces artificial
coordinate jumps when particles cross a periodic boundary before RMSD or pair
distance calculations are made.

Generate or refresh the product report for a run directory:

```bash
cargo run -p md-cli -- report runs/lj-fluid-001
```

Inspect and verify reproducibility metadata without rerunning the simulation:

```bash
cargo run -p md-cli -- inspect runs/lj-fluid-001
cargo run -p md-cli -- reproduce runs/lj-fluid-001
```

`run-manifest.json` records `fnv1a64:<hex>` hashes for `config.toml` plus
copied input/topology files when present, engine version/git metadata, platform
metadata, Rayon thread count, and the command line captured at run creation.
`reproduce` validates those recorded files still match the manifest and fails
fast if a file is missing, metadata is absent, or a hash differs.

Resume a completed or interrupted run from the recorded checkpoint:

```bash
cargo run -p md-cli -- resume runs/lj-fluid-001 --additional-steps 100
```

## Molecular File Input

Phase 7 adds simple molecular coordinate input from standard XYZ files:

```toml
[input]
path = "../molecules/argon-cluster.xyz"
format = "xyz"
frame = 0
```

Input paths are resolved relative to the config file location. The selected XYZ
frame provides atom count, element labels, and initial coordinates. The current
runner still uses a single mass from `[system]` for all atoms and initializes
velocities from `[system].seed` and `[simulation].temperature`.

Run the XYZ-input example:

```bash
cargo run -p md-cli -- run examples/validation/xyz-input.toml
```

For reproducibility, runs started from an input file copy that source coordinate
file into the run directory as `input.xyz` and record input metadata in
`summary.json`.

The runner also accepts a deliberately small PDB coordinate subset:

```toml
[input]
path = "../molecules/alanine-fragment.pdb"
format = "pdb"
frame = 0
```

Supported PDB records are `ATOM`, `HETATM`, `MODEL`, `ENDMDL`, and `END`.
Coordinates are read from fixed-width PDB columns. Atom name, residue name,
residue sequence ID, chain ID, and element are preserved when available. PDB
runs copy the source file as `input.pdb`, write parsed metadata to
`atom-metadata.csv`, and record both files in `summary.json`,
`run-manifest.json`, and `run-report.md`.

This is intentionally not full PDB compliance: alternate locations, insertion
codes, occupancies, B-factors, connectivity records, unit cells, and residue
semantics are not interpreted yet.

Run the PDB-input example:

```bash
cargo run -p md-cli -- run examples/validation/pdb-input.toml
```

## Simple Topology And Force Fields

The molecular-mechanics layer uses a small TOML topology file:

```toml
[topology]
path = "../topology/argon-dimer.toml"

[coulomb]
enabled = true
constant = 0.1
cutoff = 3.0
shift_potential = false
```

Topology files can define per-atom charges, harmonic bonds, harmonic angles,
periodic dihedrals, and explicit non-bonded exclusions:

```toml
charges = [0.0, 0.0, 0.0]

[[bonds]]
i = 0
j = 1
k = 25.0
r0 = 1.25

[[bonds]]
i = 1
j = 2
k = 25.0
r0 = 1.25

[[angles]]
i = 0
j = 1
k = 2
force_constant = 10.0
theta0 = 1.5707963267948966

[[dihedrals]]
i = 0
j = 1
k = 2
l = 3
force_constant = 0.4
multiplicity = 3
phase = 0.0

[[exclusions]]
i = 0
j = 2
```

The runner computes Lennard-Jones forces first, then adds optional Coulomb and
harmonic bonded forces into the same Velocity Verlet step. Bonded atom pairs are
automatically excluded from LJ/Coulomb pair loops; explicit `[[exclusions]]`
can add more skipped non-bonded pairs. Angle `theta0` values are in radians.
Dihedral `phase` values are in radians, and the implemented periodic torsion
form is `force_constant * (1 + cos(multiplicity * phi - phase))`.

Configs can also provide lightweight per-element force-field defaults:

```toml
[force_field]
mixing_rule = "lorentz-berthelot"

[force_field.types.Ar]
mass = 2.0
sigma = 1.0
epsilon = 1.0
charge = 0.5

[force_field.types.Ne]
mass = 4.0
sigma = 2.0
epsilon = 4.0
charge = -0.25
```

Type names currently match the atom element labels loaded from the generated
system, XYZ input, or PDB input. Missing `sigma` and `epsilon` values fall back
to the global `[force]` parameters, so existing configs keep their old uniform
LJ behavior. Force-field charge defaults are applied before topology files are
loaded; topology `charges = [...]` still override them when present.

The topology file is copied into the run directory as `topology.toml`, and
`summary.json` records bond, angle, dihedral, exclusion, and Coulomb settings.

Run the force-field examples:

```bash
cargo run -p md-cli -- run examples/validation/bonded-dimer.toml
cargo run -p md-cli -- run examples/validation/angle-trimer.toml
cargo run -p md-cli -- run examples/validation/dihedral-chain.toml
cargo run -p md-cli -- run examples/validation/force-field-report.toml
cargo run -p md-cli -- run examples/validation/force-field-types.toml
```

The force-field report example is the most compact way to inspect all currently
implemented molecular mechanics metadata in one generated `run-report.md`: it
uses harmonic bonds, harmonic angles, one periodic dihedral, enabled Coulomb,
bonded exclusions, and one explicit exclusion.

## Product Run Directories

Phase 9 adds a lightweight product layer around run outputs. Configs may include
project metadata:

```toml
[project]
name = "Bonded dimer product run"
description = "Small XYZ + topology run used to validate bonded exclusions, reports, and run manifests."
```

Each completed run now writes:

```text
config.toml
trajectory.xyz
energy.csv
summary.json
run-manifest.json
run-report.md
checkpoint.json
```

Runs that use coordinate or topology inputs also copy those files as `input.xyz`
and `topology.toml`. Runs configured with `[checkpoint] format = "binary"` write
`checkpoint.bin` in place of `checkpoint.json`. The manifest records the stable
file names, checkpoint format, run status, project name, and report/analysis
artifacts. The report can be regenerated after analysis so `run-report.md`
includes the analysis summary.

## Checkpoint And Resume

Each dynamics run writes a checkpoint at step 0 and whenever trajectory/energy
output is emitted. The checkpoint stores the full `SystemState`, current step,
time, and potential energy so it can be used for continuation.

JSON checkpoints remain the default and are easy to inspect:

```toml
[checkpoint]
format = "json"
```

For larger systems, a compact binary checkpoint can be selected:

```toml
[checkpoint]
format = "binary"
```

Binary runs write `checkpoint.bin` instead of `checkpoint.json`. The manifest
records both `checkpoint_file` and `checkpoint_format`, and `resume` reads the
recorded format. Older JSON run directories remain resumable.

Resume a run to the original configured step count when a checkpoint is behind:

```bash
cargo run -p md-cli -- resume runs/validation-bonded-dimer-001
```

Extend a completed run by additional steps:

```bash
cargo run -p md-cli -- resume runs/validation-bonded-dimer-001 --additional-steps 50
```

Resume appends to `trajectory.xyz` and `energy.csv`, refreshes `summary.json`,
updates `run-manifest.json`, and regenerates `run-report.md`.

## Limitations

This project is not a replacement for GROMACS, AMBER, NAMD, OpenMM, or Desmond.

The current engine is intended for education, prototyping, and small-scale
simulations.

Results should not be used for clinical, pharmaceutical, or high-stakes
scientific decisions without independent validation.

Current technical limitations:

- force-field support is still minimal: Lennard-Jones, optional Coulomb,
  harmonic bonds, harmonic angles, periodic dihedrals, and non-bonded exclusions
  are implemented
- thermostat support is limited to simple Berendsen velocity coupling; no
  barostat is implemented yet
- energy minimization is limited to conservative steepest descent with
  backtracking; no conjugate-gradient, L-BFGS, or constraints are implemented
  yet
- only orthorhombic periodic boxes are supported
- shifted potential is available, but force shifting is not implemented yet
- parallel force accumulation is currently CPU-only and still reduces local
  force buffers on the host; no GPU or domain-decomposition backend exists yet
- analysis assumes stable atom ordering and atom count; periodic unwrap assumes
  particles move less than half the box length between saved frames and uses
  the first frame as the initial continuous image
- molecular input currently supports XYZ plus a deliberately small PDB subset;
  PDB residue and chain metadata are preserved for reporting, but not yet used
  for topology inference, unit cells, or connectivity
- topology support currently covers per-atom charges, harmonic bonds, harmonic
  angles, periodic dihedrals, and non-bonded exclusions only; constraints,
  bonded parameter typing, residue-based assignment, and mixing rules beyond
  Lorentz-Berthelot are not implemented yet
- compact binary checkpoints are supported, but they are a project-local format
  rather than a stable cross-tool interchange format

## Roadmap

1. Project skeleton and minimal LJ engine.
2. Correctness tests and deterministic validation examples.
3. Periodic boundary conditions and minimum image convention.
4. Neighbor list and cell list acceleration.
5. Careful CPU parallelism.
6. Analysis tools and reproducible reports.
7. Simple molecular file support.
8. Basic molecular mechanics terms.
9. Product run directories, manifests, and reports.
