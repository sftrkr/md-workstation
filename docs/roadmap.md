# Roadmap

This roadmap is trimmed to the current codebase. It lists near-term product
areas and research backlog items after the current v0.1 foundation.

## Current Baseline

The project now has a credible v0.1 foundation:

- CLI-first Rust workspace
- Lennard-Jones dynamics with Velocity Verlet
- NVE plus deterministic Berendsen thermostat support
- orthorhombic open and periodic boxes
- minimum-image distances and periodic wrapping
- cutoff and shifted-potential support
- cell-list neighbor search
- Rayon CPU parallel force paths
- XYZ input and trajectory output
- small PDB coordinate subset with atom metadata preservation
- JSON and binary checkpoints
- checkpoint resume
- run directories, summaries, manifests, reports, and reproducibility metadata
- RMSD, pair-distance, energy drift, and periodic unwrap analysis
- harmonic bonds, harmonic angles, periodic dihedrals
- cutoff Coulomb
- non-bonded exclusions
- per-element force-field defaults and LJ mixing
- energy minimization
- curated validation suite
- optional desktop UI prototype
- GitHub CI and multi-platform release artifacts

## Recommended Next Work

### 1. Diagnostics And Error Polish

Goal: Make invalid or unstable runs easier to understand.

Suggested work:

- overlap diagnostics with atom indices and distance
- cutoff/box/skin recovery suggestions
- invalid topology index hints
- energy explosion detection during dynamics
- timestep warning heuristics
- clearer periodic box sizing and lattice spacing messages

Acceptance criteria:

- Errors remain structured internally.
- CLI messages include a useful recovery hint where possible.
- Tests cover cutoff/box, overlap, and topology-index diagnostics.

### 2. Public API And Crate Documentation

Goal: Make the Rust crates easier to understand before the API surface grows.

Suggested work:

- crate-level docs for core crates
- examples for intended public APIs
- `cargo doc --workspace --no-deps`
- public API review for `md-core`, `md-force`, `md-integrator`, `md-io`, and
  `md-analysis`
- note that the API is pre-1.0 unless explicitly stabilized

Acceptance criteria:

- Important public structs and functions have useful docs.
- `cargo doc --workspace --no-deps` succeeds.
- README or docs state the API stability expectation.

### 3. Analysis Expansion

Goal: Add useful trajectory analyses without pretending to be a full MD suite.

Recommended first slice:

- radius of gyration
- RMSF
- radial distribution function for one simple atom or element selection

Acceptance criteria:

- Outputs deterministic CSV files.
- Updates `analysis-summary.json`.
- Records assumptions around reduced units, stable atom order, finite boxes, and
  RDF normalization.
- Includes tests with tiny trajectories.

### 4. Benchmarking V2

Goal: Keep performance work honest and repeatable.

Suggested work:

- Criterion benchmarks for core force kernels
- benchmark scenarios for 128, 512, 1024, and 4096 particles
- integration-step benchmark
- XYZ I/O benchmark
- optional CLI benchmark command only if it avoids large accidental outputs

Acceptance criteria:

- `cargo bench` works without extra system dependencies.
- Benchmarks avoid writing large run directories unless explicitly requested.
- Docs explain how to interpret noisy local benchmark results.

### 5. HTML Report Export

Goal: Add a static self-contained HTML report while keeping Markdown as the
simple canonical format.

Suggested CLI:

```bash
md report runs/run-001 --format markdown
md report runs/run-001 --format html
```

Acceptance criteria:

- Existing Markdown behavior remains unchanged.
- HTML includes run metadata, force-field settings, final energy, analysis
  summary, output files, and limitations.
- No web service or external renderer is required.

### 6. Binary Trajectory Format

Goal: Add a Rust-native binary trajectory format only after text workflows
remain stable.

Suggested work:

- `.rtraj` magic header
- format version
- atom count
- frame metadata
- reduced-unit metadata
- byte-order documentation
- conversion from existing XYZ trajectories

Acceptance criteria:

- Format is documented.
- Reader rejects wrong magic/version cleanly.
- Round-trip tests cover small trajectories.
- XYZ remains the easiest human-readable default.

### 7. Units System V2

Goal: Clarify units and avoid accidental overclaiming.

Suggested work:

- keep reduced LJ units as the default
- add explicit unit metadata in summaries, manifests, and reports
- document every parameter in reduced units
- prepare, but do not rush, a real-units mode

Acceptance criteria:

- Reports and summaries consistently state units.
- Config validation rejects ambiguous unit-mode combinations.
- Docs explain that real biomolecular units are not implemented yet.

### 8. Release Polish

Goal: Make public releases easier to consume.

Suggested work:

- changelog
- release notes template
- artifact smoke checks after upload
- optional detached signatures after checksums are stable

Acceptance criteria:

- Release notes distinguish features, fixes, validation, and limitations.
- Artifacts remain available for Linux, Windows, and macOS.
- Checksums are documented in the release process.

## Research Backlog

These areas need design notes and validation strategies before implementation:

- Langevin thermostat
- constraints
- pressure, virial, NPT, and barostats
- Ewald or Particle Mesh Ewald
- residue templates and automatic atom typing
- water and ion models
- ligand parameterization
- mmCIF and broader molecular I/O
- enhanced sampling
- free-energy methods
- GPU backends
- SIMD-specific kernels
- domain decomposition
