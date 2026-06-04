# MD Workstation

[![CI](https://github.com/sftrkr/md-workstation/actions/workflows/ci.yml/badge.svg)](https://github.com/sftrkr/md-workstation/actions/workflows/ci.yml)

MD Workstation is a Rust-native, CPU-focused molecular dynamics workstation for
education, prototyping, force-field experiments, trajectory analysis, and
reproducible small-scale reports.

The project is CLI-first and ships with its own built-in MD engine. It does not
wrap GROMACS, AMBER, NAMD, OpenMM, Desmond, RDKit, Open Babel, Python, Conda,
Docker, WSL2, or GPU tooling.

## Current Scope

- Velocity Verlet dynamics in reduced Lennard-Jones units.
- Orthorhombic open or periodic boxes with minimum-image pair distances.
- Lennard-Jones, optional Coulomb, harmonic bonds, harmonic angles, periodic
  dihedrals, bonded exclusions, and explicit non-bonded exclusions.
- Per-element force-field defaults for mass, sigma, epsilon, and charge with
  Lorentz-Berthelot LJ mixing.
- XYZ input, a deliberately small PDB coordinate subset, and generated lattice
  systems.
- Cell-list neighbor search, optional Rayon CPU parallel force evaluation, and
  benchmark comparisons.
- NVE dynamics by default plus a simple deterministic Berendsen thermostat.
- Conservative steepest-descent minimization with backtracking.
- Energy, RMSD, pair-distance, and periodic-unwrap analysis.
- Product run directories with summaries, manifests, reports, checkpoints, and
  reproducibility metadata.
- Optional egui/eframe desktop UI prototype behind the `desktop` feature.

## What This Is Not

- Not a replacement for production MD engines.
- Not a clinically validated drug-design or pharmaceutical simulation package.
- Not a general-purpose molecular file conversion toolkit.
- Not a high-performance GPU or domain-decomposition engine.

## Quick Start

Run the default example:

```bash
cargo run -p md-cli -- run examples/lj-fluid.toml
```

Validate a config without running it:

```bash
cargo run -p md-cli -- validate examples/lj-fluid.toml
```

Run the compact validation suite:

```bash
cargo run -p md-cli -- validate-suite
```

Run tests:

```bash
cargo test
```

A dynamics run writes a directory like:

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

## CLI Commands

```bash
# Dynamics
cargo run -p md-cli -- run examples/validation/small-nve.toml

# Energy minimization
cargo run -p md-cli -- minimize examples/validation/minimize-lj.toml

# Config and suite validation
cargo run -p md-cli -- validate examples/validation/force-field-report.toml
cargo run -p md-cli -- validate-suite
cargo run -p md-cli -- validate-suite --long

# Force-path benchmarks
cargo run -p md-cli -- bench-neighbor examples/validation/larger-neighbor-compare.toml --repeats 10

# Run-directory analysis and report refresh
cargo run -p md-cli -- analyze runs/lj-fluid-001 --reference-frame 0 --pair 0,1
cargo run -p md-cli -- report runs/lj-fluid-001

# Reproducibility inspection
cargo run -p md-cli -- inspect runs/lj-fluid-001
cargo run -p md-cli -- reproduce runs/lj-fluid-001

# Resume from checkpoint
cargo run -p md-cli -- resume runs/lj-fluid-001 --additional-steps 100
```

## Desktop UI

The CLI remains the primary workflow. The optional desktop prototype lives in
`md-ui` and is only built when requested:

```bash
cargo run -p md-ui --features desktop
```

The UI can list validation configs, edit TOML, invoke the Rust CLI run workflow,
list run directories, preview `trajectory.xyz`, plot total energy from
`energy.csv`, and display `run-report.md`. It does not call external MD
engines.

Running `md-ui` without the feature prints a short launch hint:

```bash
cargo run -p md-ui
```

## Workspace Layout

```text
crates/
  md-core/        System state, boxes, validation, and shared energy samples
  md-force/       LJ, Coulomb, bonded terms, exclusions, and force reports
  md-integrator/  Velocity Verlet integration
  md-neighbor/    Cell-list neighbor search and rebuild checks
  md-analysis/    Energy drift, RMSD, pair distances, and periodic unwrap
  md-report/      Markdown run report rendering
  md-io/          XYZ, PDB subset, CSV, JSON, and checkpoint I/O
  md-cli/         CLI workflow
  md-ui/          Optional desktop UI prototype
examples/
  lj-fluid.toml
  molecules/
  topology/
  validation/
```

## Units And Assumptions

The current engine uses reduced Lennard-Jones units:

- mass is dimensionless
- sigma defines the length scale
- epsilon defines the energy scale
- Boltzmann constant is treated as 1.0 for temperature estimates
- the current temperature estimate uses 3N degrees of freedom and does not yet
  subtract center-of-mass or constraint degrees of freedom

These assumptions are intentionally simple so the engine remains transparent and
reproducible while the project grows.

## Configuration Guide

Minimal generated-system dynamics config:

```toml
[simulation]
steps = 100
dt = 0.001
output_interval = 10
temperature = 0.2

[system]
particles = 8
mass = 1.0
element = "Ar"
seed = 123
lattice_spacing = 1.4

[box]
x = 6.0
y = 6.0
z = 6.0
periodic = false

[force]
type = "lennard-jones"
epsilon = 1.0
sigma = 1.0
cutoff = 2.5
shift_potential = false

[output]
directory = "runs"
name = "small-nve"
```

Optional coordinate input:

```toml
[input]
path = "../molecules/argon-cluster.xyz"
format = "xyz"
frame = 0
```

The runner supports XYZ plus a deliberately small PDB subset: `ATOM`, `HETATM`,
`MODEL`, `ENDMDL`, and `END`. PDB atom name, residue name, residue ID, chain ID,
and element metadata are preserved for reports when available.

Periodic boundaries:

```toml
[box]
x = 10.0
y = 10.0
z = 10.0
periodic = true
```

When periodic boundaries are enabled, the cutoff must be no larger than half the
shortest box length. Positions are wrapped after integration steps. Shifted LJ
potential is supported; force shifting is not implemented yet.

Neighbor list and parallel execution:

```toml
[neighbor]
enabled = true
skin = 0.3
rebuild_interval = 10

[execution]
parallel = true
```

Rayon chooses the thread count by default. Set `RAYON_NUM_THREADS` for manual
experiments.

Thermostat:

```toml
[thermostat]
type = "berendsen"
target_temperature = 0.2
tau = 0.02
```

The default is `type = "none"` for NVE. The Berendsen thermostat is useful for
gentle deterministic temperature coupling in examples, but it is not a rigorous
canonical ensemble sampler.

Minimization:

```toml
[minimization]
steps = 80
output_interval = 10
step_size = 0.001
max_displacement = 0.02
force_tolerance = 0.000001
max_backtracks = 12
```

Checkpoint format:

```toml
[checkpoint]
format = "json"   # or "binary"
```

## Topology And Force Fields

Configs can reference a small TOML topology file:

```toml
[topology]
path = "../topology/argon-force-field.toml"

[coulomb]
enabled = true
constant = 0.1
cutoff = 3.0
shift_potential = false
```

Topology files can define charges, bonds, angles, periodic dihedrals, and
explicit exclusions:

```toml
charges = [0.5, 0.0, -0.5]

[[bonds]]
i = 0
j = 1
k = 25.0
r0 = 1.2

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

Bonded pairs are excluded from non-bonded LJ/Coulomb loops automatically.
Explicit exclusions add more skipped non-bonded pairs. Angle and dihedral phases
are in radians.

Per-element force-field defaults:

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

Type names currently match element labels from generated systems, XYZ input, or
PDB input. Missing `sigma` and `epsilon` values fall back to global `[force]`
parameters. Topology `charges = [...]` override force-field charge defaults.

## Run Outputs And Reproducibility

Every run records product metadata in `summary.json`, `run-manifest.json`, and
`run-report.md`.

`run-manifest.json` includes:

- stable output file names
- run status and project metadata
- checkpoint file and checkpoint format
- config, input, and topology hashes in `fnv1a64:<hex>` format
- engine version and git commit when available
- Rust target, platform, Rayon thread count, and captured command line

Runs that use coordinate or topology inputs copy those files into the run
directory as `input.xyz`, `input.pdb`, and/or `topology.toml`.

Use:

```bash
cargo run -p md-cli -- inspect runs/lj-fluid-001
cargo run -p md-cli -- reproduce runs/lj-fluid-001
```

`inspect` prints the reproducibility metadata. `reproduce` verifies that the
recorded config/input/topology files still match the manifest without rerunning
the simulation.

## Analysis

Analyze a completed run directory:

```bash
cargo run -p md-cli -- analyze runs/lj-fluid-001 --reference-frame 0 --pair 0,1
```

The analysis command reads `energy.csv` and `trajectory.xyz`, then writes:

```text
analysis-summary.json
analysis-rmsd.csv
analysis-distance.csv
analysis-report.md
```

Current analysis includes energy drift, temperature summaries, centered
Kabsch-style RMSD alignment, one atom-pair distance series, and periodic unwrap
when the run config contains an orthorhombic periodic box.

Refresh the run report after analysis:

```bash
cargo run -p md-cli -- report runs/lj-fluid-001
```

## Validation Examples

Useful examples:

```bash
cargo run -p md-cli -- run examples/validation/small-nve.toml
cargo run -p md-cli -- run examples/validation/long-nve-drift.toml
cargo run -p md-cli -- run examples/validation/pbc-small-nve.toml
cargo run -p md-cli -- run examples/validation/neighbor-compare.toml
cargo run -p md-cli -- run examples/validation/parallel-open.toml
cargo run -p md-cli -- run examples/validation/xyz-input.toml
cargo run -p md-cli -- run examples/validation/pdb-input.toml
cargo run -p md-cli -- minimize examples/validation/minimize-lj.toml
cargo run -p md-cli -- run examples/validation/force-field-report.toml
cargo run -p md-cli -- run examples/validation/force-field-types.toml
cargo run -p md-cli -- run examples/validation/binary-checkpoint.toml
```

`validate-suite` performs a compact config/input/topology/force setup smoke
check. The default suite skips longer stress configs. `--long` includes
`long-nve-drift.toml` and the larger neighbor benchmark config. The manual
`long-nve-drift-100k.toml` file is intentionally outside the suite.

For a manual longer stress check, copy or edit the long NVE config and use a
larger `simulation.steps` value such as `100000`. Keep `dt = 0.0001` unless you
are intentionally exploring timestep sensitivity.

## Development

Common checks:

```bash
cargo fmt --all -- --check
cargo test
cargo check -p md-ui --features desktop
cargo clippy --workspace --all-targets -- -D warnings
cargo run -p md-cli -- validate-suite
```

The desktop UI dependencies are optional and isolated behind the `desktop`
feature. The normal CLI/test workflow does not require launching a GUI.

GitHub Actions runs the same checks on pushes to `main`, pull requests, and
manual workflow dispatches.

## Releases

Linux, Windows, and macOS binary artifacts are published from existing `v*` git
tags:

```bash
git tag -a v0.1.0 -m "v0.1.0"
git push origin v0.1.0
```

The release workflow tests the workspace, validates curated configs, builds
release binaries for `md` and `md-ui`, packages README/Cargo metadata and
examples, writes SHA-256 checksums, and publishes a GitHub Release with
generated notes. Linux and macOS packages are `.tar.gz` archives; Windows
packages are `.zip` archives. The workflow can also be started manually for an
existing tag from GitHub Actions.

## Limitations

- Reduced Lennard-Jones units only; no real unit conversion layer yet.
- Only orthorhombic boxes are supported.
- Force shifting is not implemented; shifted potential is available.
- Thermostat support is limited to Berendsen velocity coupling.
- Minimization is conservative steepest descent, not conjugate-gradient,
  L-BFGS, constraints, or SHAKE/RATTLE.
- PDB support is intentionally a small coordinate subset; unit cells,
  connectivity records, alternate locations, insertion codes, occupancies, and
  B-factors are not interpreted.
- Topology support is minimal and project-local; there is no residue template
  assignment or broad force-field file format import.
- Binary checkpoints are project-local, not a stable cross-tool interchange
  format.
- Parallelism is CPU-only and host-local; there is no GPU or domain
  decomposition backend.
- Analysis assumes stable atom ordering and atom count. Periodic unwrap assumes
  particles move less than half the box length between saved frames.

## Roadmap

The completed foundation now covers the original CLI-first phases. Sensible next
directions are:

- richer trajectory storage for larger runs
- stronger minimizers and constraints
- broader molecular file support
- residue/template-based force-field assignment
- better UI ergonomics around config editing and run monitoring
- packaging and release automation
