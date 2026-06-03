use std::collections::BTreeMap;
use std::fs::{self, File, OpenOptions};
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::time::{Duration, Instant, SystemTime, UNIX_EPOCH};

use anyhow::{bail, Context, Result};
use clap::{Parser, Subcommand};
use md_analysis::{
    pair_distance_series, rmsd_series, summarize_analysis_with_metadata, unwrap_periodic_frames,
    AnalysisSummary, PairDistanceSample, RmsdSample,
};
use md_core::{EnergySample, SimulationBox, SystemState};
use md_force::{
    add_coulomb_forces_with_options_and_exclusions, add_harmonic_angle_forces,
    add_harmonic_bond_forces, add_periodic_dihedral_forces,
    compute_lennard_jones_forces_parallel_with_mixed_neighbor_list_and_exclusions,
    compute_lennard_jones_forces_parallel_with_mixed_options_and_exclusions,
    compute_lennard_jones_forces_parallel_with_neighbor_list_and_exclusions,
    compute_lennard_jones_forces_parallel_with_options_and_exclusions,
    compute_lennard_jones_forces_with_mixed_neighbor_list_and_exclusions,
    compute_lennard_jones_forces_with_mixed_options_and_exclusions,
    compute_lennard_jones_forces_with_neighbor_list_and_exclusions,
    compute_lennard_jones_forces_with_options_and_exclusions, normalize_excluded_pairs,
    BoundaryCondition, CoulombOptions, ExcludedPair, ForceError, HarmonicAngle, HarmonicBond,
    LennardJonesMixingRule, LennardJonesOptions, LennardJonesParams, LennardJonesParticleParams,
    MixedLennardJonesOptions, PeriodicDihedral,
};
use md_integrator::{IntegratorError, VelocityVerlet};
use md_io::{
    read_checkpoint_binary, read_checkpoint_json, read_energy_csv, read_manifest_json,
    read_pdb_coordinate_frames, read_summary_json, read_xyz_coordinate_frames, read_xyz_trajectory,
    write_atom_metadata_csv, write_checkpoint_binary, write_checkpoint_json, write_manifest_json,
    write_summary_json, write_xyz_frame, EnergyCsvWriter, PdbFrame, RunCheckpoint, RunManifest,
    RunSummary,
};
use md_neighbor::{NeighborBoundary, NeighborList, NeighborListConfig};
use md_report::{render_markdown_report, MarkdownReportInput};
use rand::{Rng, SeedableRng};
use rand_chacha::ChaCha8Rng;
use serde::{Deserialize, Serialize};

#[derive(Debug, Parser)]
#[command(name = "md")]
#[command(about = "Rust-native molecular dynamics workstation")]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Debug, Subcommand)]
enum Commands {
    /// Run a simulation from a TOML config.
    Run { config: PathBuf },
    /// Minimize an initial structure from a TOML config.
    Minimize { config: PathBuf },
    /// Validate a simulation config.
    Validate { config: PathBuf },
    /// Compare naive and neighbor-list force evaluation for a config.
    BenchNeighbor {
        config: PathBuf,
        #[arg(long, default_value_t = 20)]
        repeats: usize,
    },
    /// Analyze a completed run directory.
    Analyze {
        run_dir: PathBuf,
        #[arg(long, default_value_t = 0)]
        reference_frame: usize,
        #[arg(long)]
        pair: Option<String>,
    },
    /// Generate a Markdown report for a completed run directory.
    Report { run_dir: PathBuf },
    /// Continue a run directory from checkpoint.json.
    Resume {
        run_dir: PathBuf,
        #[arg(long)]
        additional_steps: Option<usize>,
    },
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    match cli.command {
        Commands::Run { config } => run_command(&config),
        Commands::Minimize { config } => minimize_command(&config),
        Commands::Validate { config } => {
            let run_config = load_config(&config)?;
            run_config.validate()?;
            let (mut state, topology) =
                prepare_initial_state_and_topology(&run_config, config.parent())?;
            apply_topology_to_state(&mut state, topology.as_ref())?;
            let force_options = run_config.force_options()?;
            run_config.mixed_lennard_jones_options(&state, &force_options)?;
            println!("Config is valid: {}", config.display());
            Ok(())
        }
        Commands::BenchNeighbor { config, repeats } => benchmark_neighbor_command(&config, repeats),
        Commands::Analyze {
            run_dir,
            reference_frame,
            pair,
        } => analyze_command(&run_dir, reference_frame, pair.as_deref()),
        Commands::Report { run_dir } => report_command(&run_dir),
        Commands::Resume {
            run_dir,
            additional_steps,
        } => resume_command(&run_dir, additional_steps),
    }
}

fn run_command(config_path: &Path) -> Result<()> {
    let config = load_config(config_path)?;
    config.validate()?;
    let config_dir = config_path.parent();

    let input_source_path = config
        .input
        .as_ref()
        .map(|input| resolve_config_relative_path(&input.path, config_dir));
    let (mut state, topology) = prepare_initial_state_and_topology(&config, config_dir)?;
    apply_topology_to_state(&mut state, topology.as_ref())?;
    let force_options = config.force_options()?;
    let mixed_force_options = config.mixed_lennard_jones_options(&state, &force_options)?;
    let coulomb_options = config.coulomb_options()?;
    let thermostat_options = config.thermostat_options();
    let mut neighbor_list = if config.neighbor.enabled {
        Some(NeighborList::build(
            &state,
            config.neighbor_list_config(&force_options)?,
        )?)
    } else {
        None
    };
    let initial_potential_energy = compute_configured_forces(
        &mut state,
        &force_options,
        mixed_force_options.as_ref(),
        neighbor_list.as_ref(),
        topology.as_ref(),
        coulomb_options.as_ref(),
        config.execution.parallel,
    )?;
    let integrator = VelocityVerlet::new(config.simulation.dt)?;

    let run_name = config
        .output
        .name
        .clone()
        .unwrap_or_else(|| config_stem(config_path));
    let run_name = sanitize_run_name(&run_name);
    let run_dir = next_run_dir(&config.output.directory, &run_name)?;
    let checkpoint_file = config.checkpoint.file_name().to_string();
    let checkpoint_format = config.checkpoint.normalized_format();
    fs::create_dir_all(&run_dir)
        .with_context(|| format!("failed to create run directory {}", run_dir.display()))?;

    fs::copy(config_path, run_dir.join("config.toml")).with_context(|| {
        format!(
            "failed to copy config into run directory {}",
            run_dir.display()
        )
    })?;

    let input_file =
        if let (Some(input), Some(input_source_path)) = (&config.input, &input_source_path) {
            let input_file = copied_input_file_name(input);
            fs::copy(input_source_path, run_dir.join(&input_file)).with_context(|| {
                format!(
                    "failed to copy input file {} into run directory {}",
                    input_source_path.display(),
                    run_dir.display()
                )
            })?;
            Some(input_file)
        } else {
            None
        };

    let topology_file = if let Some(topology) = &topology {
        fs::copy(&topology.source_path, run_dir.join("topology.toml")).with_context(|| {
            format!(
                "failed to copy topology file {} into run directory {}",
                topology.source_path.display(),
                run_dir.display()
            )
        })?;
        Some("topology.toml".to_string())
    } else {
        None
    };
    let atom_metadata_file = write_optional_atom_metadata_file(&run_dir, &state)?;

    let trajectory_file = File::create(run_dir.join("trajectory.xyz")).with_context(|| {
        format!(
            "failed to create {}",
            run_dir.join("trajectory.xyz").display()
        )
    })?;
    let mut trajectory = BufWriter::new(trajectory_file);

    let energy_file = File::create(run_dir.join("energy.csv"))
        .with_context(|| format!("failed to create {}", run_dir.join("energy.csv").display()))?;
    let mut energy_writer = EnergyCsvWriter::new(BufWriter::new(energy_file))?;

    println!("Running simulation...");
    println!("Particles: {}", state.particle_count());
    println!("Steps: {}", config.simulation.steps);
    println!("dt: {}", integrator.dt());
    println!("Output interval: {}", config.simulation.output_interval);
    println!("Boundary: {}", config.boundary_label());
    if let Some(input_source_path) = &input_source_path {
        println!(
            "Input: {} frame {}",
            input_source_path.display(),
            config.input.as_ref().map_or(0, |input| input.frame)
        );
    }
    println!("Shifted potential: {}", config.force.shift_potential);
    println!("Neighbor list: {}", config.neighbor.enabled);
    println!("Ensemble: {}", config.ensemble_label());
    println!("Thermostat: {}", thermostat_options.label());
    if let Some(target_temperature) = thermostat_options.target_temperature() {
        println!("Thermostat target temperature: {}", target_temperature);
    }
    if let Some(tau) = thermostat_options.tau() {
        println!("Thermostat tau: {}", tau);
    }
    if let Some(topology) = &topology {
        println!(
            "Topology: {} (bonds: {}, angles: {}, dihedrals: {}, exclusions: {})",
            topology.source_path.display(),
            topology.bonds.len(),
            topology.angles.len(),
            topology.dihedrals.len(),
            topology.exclusions.len()
        );
    }
    println!("Coulomb: {}", coulomb_options.is_some());
    println!("Parallel force: {}", config.execution.parallel);
    if config.execution.parallel {
        println!("Rayon threads: {}", rayon::current_num_threads());
    }
    println!("Output: {}", run_dir.display());

    let mut final_sample = sample_energy(0, integrator.dt(), &state, initial_potential_energy);
    write_xyz_frame(&mut trajectory, &state, 0)?;
    energy_writer.write_sample(&final_sample)?;
    write_run_checkpoint(
        &run_dir,
        &run_name,
        final_sample.step,
        final_sample.time,
        initial_potential_energy,
        &state,
        &config.checkpoint,
    )?;

    let progress_interval = (config.simulation.steps / 10).max(1);

    for step in 1..=config.simulation.steps {
        let force_rebuild = step % config.neighbor.rebuild_interval == 0;
        let potential_energy =
            integrator.step_with_force(&mut state, force_options.boundary, |state| {
                if let Some(neighbor_list) = neighbor_list.as_mut() {
                    if force_rebuild || neighbor_list.needs_rebuild(state)? {
                        neighbor_list.rebuild(state)?;
                    }
                }
                compute_configured_forces(
                    state,
                    &force_options,
                    mixed_force_options.as_ref(),
                    neighbor_list.as_ref(),
                    topology.as_ref(),
                    coulomb_options.as_ref(),
                    config.execution.parallel,
                )
                .map_err(IntegratorError::from)
            })?;
        apply_thermostat(&mut state, thermostat_options, integrator.dt())?;

        if should_output(
            step,
            config.simulation.steps,
            config.simulation.output_interval,
        ) {
            final_sample = sample_energy(step, integrator.dt(), &state, potential_energy);
            ensure_finite_energy(&final_sample)?;
            write_xyz_frame(&mut trajectory, &state, step)?;
            energy_writer.write_sample(&final_sample)?;
            write_run_checkpoint(
                &run_dir,
                &run_name,
                final_sample.step,
                final_sample.time,
                potential_energy,
                &state,
                &config.checkpoint,
            )?;
        }

        if step % progress_interval == 0 || step == config.simulation.steps {
            let sample = sample_energy(step, integrator.dt(), &state, potential_energy);
            println!(
                "Step {}/{}  E_kin={:.6}  E_pot={:.6}  E_total={:.6}  T={:.6}",
                step,
                config.simulation.steps,
                sample.kinetic,
                sample.potential,
                sample.total,
                sample.temperature
            );
        }
    }

    trajectory.flush()?;
    energy_writer.flush()?;

    let mut outputs = vec![
        "config.toml".to_string(),
        "trajectory.xyz".to_string(),
        "energy.csv".to_string(),
        "summary.json".to_string(),
    ];
    if let Some(input_file) = &input_file {
        outputs.push(input_file.clone());
    }
    if topology_file.is_some() {
        outputs.push("topology.toml".to_string());
    }
    if atom_metadata_file.is_some() {
        outputs.push("atom-metadata.csv".to_string());
    }
    outputs.push("run-report.md".to_string());
    outputs.push("run-manifest.json".to_string());
    outputs.push(checkpoint_file.clone());

    let summary = RunSummary {
        run_name: run_name.clone(),
        project_name: config.project.name.clone(),
        project_description: config.project.description.clone(),
        particle_count: state.particle_count(),
        steps: config.simulation.steps,
        dt: integrator.dt(),
        output_interval: config.simulation.output_interval,
        seed: config.system.seed,
        workflow: "dynamics".to_string(),
        ensemble: config.ensemble_label().to_string(),
        thermostat: thermostat_options.label().to_string(),
        thermostat_target_temperature: thermostat_options.target_temperature(),
        thermostat_tau: thermostat_options.tau(),
        minimization_step_size: None,
        minimization_force_tolerance: None,
        minimization_final_max_force: None,
        input_file,
        input_format: config.input.as_ref().map(|input| input.normalized_format()),
        input_frame: config.input.as_ref().map(|input| input.frame),
        atom_metadata_file,
        topology_file,
        bond_count: topology.as_ref().map_or(0, |topology| topology.bonds.len()),
        angle_count: topology
            .as_ref()
            .map_or(0, |topology| topology.angles.len()),
        dihedral_count: topology
            .as_ref()
            .map_or(0, |topology| topology.dihedrals.len()),
        excluded_pair_count: topology
            .as_ref()
            .map_or(0, |topology| topology.exclusions.len()),
        coulomb: coulomb_options.is_some(),
        coulomb_cutoff: coulomb_options.as_ref().map(|options| options.cutoff),
        force: "lennard-jones".to_string(),
        boundary: config.boundary_label().to_string(),
        shifted_potential: config.force.shift_potential,
        neighbor_list: config.neighbor.enabled,
        neighbor_skin: config.neighbor.enabled.then_some(config.neighbor.skin),
        neighbor_rebuild_interval: config
            .neighbor
            .enabled
            .then_some(config.neighbor.rebuild_interval),
        parallel: config.execution.parallel,
        rayon_threads: config
            .execution
            .parallel
            .then_some(rayon::current_num_threads()),
        units: "reduced-lennard-jones".to_string(),
        final_energy: final_sample,
        outputs,
    };
    write_summary_json(run_dir.join("summary.json"), &summary)?;

    let completed_at = current_unix_seconds();
    let manifest = RunManifest {
        schema_version: 1,
        application: "md-workstation".to_string(),
        run_name: run_name.clone(),
        project_name: config.project.name.clone(),
        status: "completed".to_string(),
        created_at_unix_seconds: completed_at,
        completed_at_unix_seconds: Some(completed_at),
        config_file: "config.toml".to_string(),
        summary_file: "summary.json".to_string(),
        trajectory_file: "trajectory.xyz".to_string(),
        energy_file: "energy.csv".to_string(),
        input_file: summary.input_file.clone(),
        topology_file: summary.topology_file.clone(),
        checkpoint_file: Some(checkpoint_file),
        checkpoint_format: Some(checkpoint_format),
        analysis_summary_file: None,
        report_file: Some("run-report.md".to_string()),
        outputs: summary.outputs.clone(),
    };
    let report = render_markdown_report(&MarkdownReportInput {
        run_dir: &run_dir.display().to_string(),
        summary: &summary,
        manifest: Some(&manifest),
        analysis: None,
    });
    fs::write(run_dir.join("run-report.md"), report).with_context(|| {
        format!(
            "failed to write {}",
            run_dir.join("run-report.md").display()
        )
    })?;
    write_manifest_json(run_dir.join("run-manifest.json"), &manifest)?;

    println!("Simulation completed.");
    println!("Output: {}", run_dir.display());

    Ok(())
}

fn minimize_command(config_path: &Path) -> Result<()> {
    let config = load_config(config_path)?;
    config.validate()?;
    let config_dir = config_path.parent();

    let input_source_path = config
        .input
        .as_ref()
        .map(|input| resolve_config_relative_path(&input.path, config_dir));
    let (mut state, topology) = prepare_initial_state_and_topology(&config, config_dir)?;
    apply_topology_to_state(&mut state, topology.as_ref())?;
    state.vx.fill(0.0);
    state.vy.fill(0.0);
    state.vz.fill(0.0);

    let force_options = config.force_options()?;
    let mixed_force_options = config.mixed_lennard_jones_options(&state, &force_options)?;
    let coulomb_options = config.coulomb_options()?;
    let mut potential_energy = compute_configured_forces(
        &mut state,
        &force_options,
        mixed_force_options.as_ref(),
        None,
        topology.as_ref(),
        coulomb_options.as_ref(),
        config.execution.parallel,
    )?;
    let mut max_force = max_force_norm(&state);
    ensure_finite_minimization_state(0, potential_energy, max_force, &state)?;

    let run_name = config
        .output
        .name
        .clone()
        .unwrap_or_else(|| config_stem(config_path));
    let run_name = sanitize_run_name(&run_name);
    let run_dir = next_run_dir(&config.output.directory, &run_name)?;
    fs::create_dir_all(&run_dir)
        .with_context(|| format!("failed to create run directory {}", run_dir.display()))?;

    fs::copy(config_path, run_dir.join("config.toml")).with_context(|| {
        format!(
            "failed to copy config into run directory {}",
            run_dir.display()
        )
    })?;

    let input_file =
        if let (Some(input), Some(input_source_path)) = (&config.input, &input_source_path) {
            let input_file = copied_input_file_name(input);
            fs::copy(input_source_path, run_dir.join(&input_file)).with_context(|| {
                format!(
                    "failed to copy input file {} into run directory {}",
                    input_source_path.display(),
                    run_dir.display()
                )
            })?;
            Some(input_file)
        } else {
            None
        };

    let topology_file = if let Some(topology) = &topology {
        fs::copy(&topology.source_path, run_dir.join("topology.toml")).with_context(|| {
            format!(
                "failed to copy topology file {} into run directory {}",
                topology.source_path.display(),
                run_dir.display()
            )
        })?;
        Some("topology.toml".to_string())
    } else {
        None
    };
    let atom_metadata_file = write_optional_atom_metadata_file(&run_dir, &state)?;

    let trajectory_file = File::create(run_dir.join("trajectory.xyz")).with_context(|| {
        format!(
            "failed to create {}",
            run_dir.join("trajectory.xyz").display()
        )
    })?;
    let mut trajectory = BufWriter::new(trajectory_file);

    let energy_file = File::create(run_dir.join("energy.csv"))
        .with_context(|| format!("failed to create {}", run_dir.join("energy.csv").display()))?;
    let mut energy_writer = EnergyCsvWriter::new(BufWriter::new(energy_file))?;

    let history_file = File::create(run_dir.join("minimization.csv")).with_context(|| {
        format!(
            "failed to create {}",
            run_dir.join("minimization.csv").display()
        )
    })?;
    let mut history_writer = BufWriter::new(history_file);
    write_minimization_history_header(&mut history_writer)?;

    println!("Minimizing structure...");
    println!("Particles: {}", state.particle_count());
    println!("Steps: {}", config.minimization.steps);
    println!("Step size: {}", config.minimization.step_size);
    println!("Max displacement: {}", config.minimization.max_displacement);
    println!("Force tolerance: {}", config.minimization.force_tolerance);
    println!("Output interval: {}", config.minimization.output_interval);
    println!("Boundary: {}", config.boundary_label());
    if let Some(input_source_path) = &input_source_path {
        println!(
            "Input: {} frame {}",
            input_source_path.display(),
            config.input.as_ref().map_or(0, |input| input.frame)
        );
    }
    println!("Shifted potential: {}", config.force.shift_potential);
    println!("Neighbor list: false (direct minimization force path)");
    if let Some(topology) = &topology {
        println!(
            "Topology: {} (bonds: {}, angles: {}, dihedrals: {}, exclusions: {})",
            topology.source_path.display(),
            topology.bonds.len(),
            topology.angles.len(),
            topology.dihedrals.len(),
            topology.exclusions.len()
        );
    }
    println!("Coulomb: {}", coulomb_options.is_some());
    println!("Parallel force: {}", config.execution.parallel);
    if config.execution.parallel {
        println!("Rayon threads: {}", rayon::current_num_threads());
    }
    println!("Initial potential: {:.10}", potential_energy);
    println!("Initial max force: {:.10}", max_force);
    println!("Output: {}", run_dir.display());

    let mut final_step = 0;
    let mut final_sample =
        minimization_energy_sample(final_step, config.minimization.step_size, potential_energy);
    let mut last_output_step = 0;
    write_xyz_frame(&mut trajectory, &state, 0)?;
    energy_writer.write_sample(&final_sample)?;
    write_minimization_history_sample(&mut history_writer, 0, potential_energy, max_force, 0.0)?;

    let progress_interval = (config.minimization.steps / 10).max(1);
    for step in 1..=config.minimization.steps {
        if max_force <= config.minimization.force_tolerance {
            break;
        }

        let mut step_scale = minimization_step_scale(max_force, &config.minimization);
        let mut accepted = None;
        let mut last_rejection = None;

        for _ in 0..=config.minimization.max_backtracks {
            let mut trial = state.clone();
            if let Err(error) =
                displace_along_forces(&mut trial, step_scale, force_options.boundary)
            {
                last_rejection = Some(error.to_string());
                step_scale *= 0.5;
                continue;
            }

            match compute_configured_forces(
                &mut trial,
                &force_options,
                mixed_force_options.as_ref(),
                None,
                topology.as_ref(),
                coulomb_options.as_ref(),
                config.execution.parallel,
            ) {
                Ok(trial_potential) => {
                    let trial_max_force = max_force_norm(&trial);
                    if trial_potential.is_finite()
                        && trial_max_force.is_finite()
                        && trial_potential <= potential_energy + 1.0e-12
                    {
                        accepted = Some((trial, trial_potential, trial_max_force, step_scale));
                        break;
                    }
                    last_rejection = Some(format!(
                        "trial potential {:.10} did not improve current potential {:.10}",
                        trial_potential, potential_energy
                    ));
                }
                Err(error) => {
                    last_rejection = Some(error.to_string());
                }
            }
            step_scale *= 0.5;
        }

        let Some((accepted_state, accepted_potential, accepted_max_force, accepted_scale)) =
            accepted
        else {
            let detail = last_rejection
                .map(|message| format!("; last rejection: {message}"))
                .unwrap_or_default();
            bail!("minimization could not find a downhill step at step {step}{detail}");
        };

        state = accepted_state;
        potential_energy = accepted_potential;
        max_force = accepted_max_force;
        final_step = step;
        ensure_finite_minimization_state(final_step, potential_energy, max_force, &state)?;
        write_minimization_history_sample(
            &mut history_writer,
            final_step,
            potential_energy,
            max_force,
            accepted_scale,
        )?;

        if should_output(
            final_step,
            config.minimization.steps,
            config.minimization.output_interval,
        ) {
            final_sample = minimization_energy_sample(
                final_step,
                config.minimization.step_size,
                potential_energy,
            );
            energy_writer.write_sample(&final_sample)?;
            write_xyz_frame(&mut trajectory, &state, final_step)?;
            last_output_step = final_step;
        }

        if final_step % progress_interval == 0
            || final_step == config.minimization.steps
            || max_force <= config.minimization.force_tolerance
        {
            println!(
                "Step {}/{}  E_pot={:.6}  max_force={:.6}  step_scale={:.6e}",
                final_step, config.minimization.steps, potential_energy, max_force, accepted_scale
            );
        }
    }

    final_sample =
        minimization_energy_sample(final_step, config.minimization.step_size, potential_energy);
    if final_step != last_output_step {
        energy_writer.write_sample(&final_sample)?;
        write_xyz_frame(&mut trajectory, &state, final_step)?;
    }

    let minimized_file = File::create(run_dir.join("minimized.xyz")).with_context(|| {
        format!(
            "failed to create {}",
            run_dir.join("minimized.xyz").display()
        )
    })?;
    let mut minimized = BufWriter::new(minimized_file);
    write_xyz_frame(&mut minimized, &state, final_step)?;

    trajectory.flush()?;
    energy_writer.flush()?;
    history_writer.flush()?;
    minimized.flush()?;

    let mut outputs = vec![
        "config.toml".to_string(),
        "trajectory.xyz".to_string(),
        "energy.csv".to_string(),
        "minimization.csv".to_string(),
        "minimized.xyz".to_string(),
        "summary.json".to_string(),
    ];
    if let Some(input_file) = &input_file {
        outputs.push(input_file.clone());
    }
    if topology_file.is_some() {
        outputs.push("topology.toml".to_string());
    }
    if atom_metadata_file.is_some() {
        outputs.push("atom-metadata.csv".to_string());
    }
    outputs.push("run-report.md".to_string());
    outputs.push("run-manifest.json".to_string());

    let summary = RunSummary {
        run_name: run_name.clone(),
        project_name: config.project.name.clone(),
        project_description: config.project.description.clone(),
        particle_count: state.particle_count(),
        steps: final_step,
        dt: config.minimization.step_size,
        output_interval: config.minimization.output_interval,
        seed: config.system.seed,
        workflow: "minimization".to_string(),
        ensemble: "minimization".to_string(),
        thermostat: "none".to_string(),
        thermostat_target_temperature: None,
        thermostat_tau: None,
        minimization_step_size: Some(config.minimization.step_size),
        minimization_force_tolerance: Some(config.minimization.force_tolerance),
        minimization_final_max_force: Some(max_force),
        input_file,
        input_format: config.input.as_ref().map(|input| input.normalized_format()),
        input_frame: config.input.as_ref().map(|input| input.frame),
        atom_metadata_file,
        topology_file,
        bond_count: topology.as_ref().map_or(0, |topology| topology.bonds.len()),
        angle_count: topology
            .as_ref()
            .map_or(0, |topology| topology.angles.len()),
        dihedral_count: topology
            .as_ref()
            .map_or(0, |topology| topology.dihedrals.len()),
        excluded_pair_count: topology
            .as_ref()
            .map_or(0, |topology| topology.exclusions.len()),
        coulomb: coulomb_options.is_some(),
        coulomb_cutoff: coulomb_options.as_ref().map(|options| options.cutoff),
        force: "lennard-jones".to_string(),
        boundary: config.boundary_label().to_string(),
        shifted_potential: config.force.shift_potential,
        neighbor_list: false,
        neighbor_skin: None,
        neighbor_rebuild_interval: None,
        parallel: config.execution.parallel,
        rayon_threads: config
            .execution
            .parallel
            .then_some(rayon::current_num_threads()),
        units: "reduced-lennard-jones".to_string(),
        final_energy: final_sample,
        outputs,
    };
    write_summary_json(run_dir.join("summary.json"), &summary)?;

    let completed_at = current_unix_seconds();
    let manifest = RunManifest {
        schema_version: 1,
        application: "md-workstation".to_string(),
        run_name: run_name.clone(),
        project_name: config.project.name.clone(),
        status: "completed".to_string(),
        created_at_unix_seconds: completed_at,
        completed_at_unix_seconds: Some(completed_at),
        config_file: "config.toml".to_string(),
        summary_file: "summary.json".to_string(),
        trajectory_file: "trajectory.xyz".to_string(),
        energy_file: "energy.csv".to_string(),
        input_file: summary.input_file.clone(),
        topology_file: summary.topology_file.clone(),
        checkpoint_file: None,
        checkpoint_format: None,
        analysis_summary_file: None,
        report_file: Some("run-report.md".to_string()),
        outputs: summary.outputs.clone(),
    };
    let report = render_markdown_report(&MarkdownReportInput {
        run_dir: &run_dir.display().to_string(),
        summary: &summary,
        manifest: Some(&manifest),
        analysis: None,
    });
    fs::write(run_dir.join("run-report.md"), report).with_context(|| {
        format!(
            "failed to write {}",
            run_dir.join("run-report.md").display()
        )
    })?;
    write_manifest_json(run_dir.join("run-manifest.json"), &manifest)?;

    println!("Minimization completed.");
    println!("Final step: {}", final_step);
    println!("Final potential: {:.10}", potential_energy);
    println!("Final max force: {:.10}", max_force);
    println!("Output: {}", run_dir.display());

    Ok(())
}

fn benchmark_neighbor_command(config_path: &Path, repeats: usize) -> Result<()> {
    if repeats == 0 {
        bail!("--repeats must be greater than zero");
    }

    let config = load_config(config_path)?;
    config.validate()?;
    let (mut state, topology) = prepare_initial_state_and_topology(&config, config_path.parent())?;
    apply_topology_to_state(&mut state, topology.as_ref())?;
    let force_options = config.force_options()?;
    let mixed_force_options = config.mixed_lennard_jones_options(&state, &force_options)?;
    let coulomb_options = config.coulomb_options()?;
    let neighbor_config = config.neighbor_list_config(&force_options)?;

    let mut naive_state = state.clone();
    let naive_report = compute_configured_force_report(
        &mut naive_state,
        &force_options,
        mixed_force_options.as_ref(),
        None,
        topology.as_ref(),
        coulomb_options.as_ref(),
        false,
    )?;

    let build_start = Instant::now();
    let neighbor_list = NeighborList::build(&state, neighbor_config)?;
    let build_elapsed = build_start.elapsed();

    let mut neighbor_state = state.clone();
    let neighbor_report = compute_configured_force_report(
        &mut neighbor_state,
        &force_options,
        mixed_force_options.as_ref(),
        Some(&neighbor_list),
        topology.as_ref(),
        coulomb_options.as_ref(),
        false,
    )?;
    let mut parallel_naive_state = state.clone();
    let parallel_naive_report = compute_configured_force_report(
        &mut parallel_naive_state,
        &force_options,
        mixed_force_options.as_ref(),
        None,
        topology.as_ref(),
        coulomb_options.as_ref(),
        true,
    )?;
    let mut parallel_neighbor_state = state.clone();
    let parallel_neighbor_report = compute_configured_force_report(
        &mut parallel_neighbor_state,
        &force_options,
        mixed_force_options.as_ref(),
        Some(&neighbor_list),
        topology.as_ref(),
        coulomb_options.as_ref(),
        true,
    )?;

    let energy_delta = (neighbor_report.potential_energy - naive_report.potential_energy).abs();
    let serial_max_force_delta = max_force_delta(&naive_state, &neighbor_state);
    let parallel_energy_delta =
        (parallel_neighbor_report.potential_energy - naive_report.potential_energy).abs();
    let parallel_max_force_delta = max_force_delta(&naive_state, &parallel_neighbor_state);
    let parallel_naive_energy_delta =
        (parallel_naive_report.potential_energy - naive_report.potential_energy).abs();
    let parallel_naive_max_force_delta = max_force_delta(&naive_state, &parallel_naive_state);

    let naive_start = Instant::now();
    for _ in 0..repeats {
        let mut repeated = state.clone();
        compute_configured_force_report(
            &mut repeated,
            &force_options,
            mixed_force_options.as_ref(),
            None,
            topology.as_ref(),
            coulomb_options.as_ref(),
            false,
        )?;
    }
    let naive_elapsed = naive_start.elapsed();

    let neighbor_start = Instant::now();
    for _ in 0..repeats {
        let mut repeated = state.clone();
        compute_configured_force_report(
            &mut repeated,
            &force_options,
            mixed_force_options.as_ref(),
            Some(&neighbor_list),
            topology.as_ref(),
            coulomb_options.as_ref(),
            false,
        )?;
    }
    let neighbor_elapsed = neighbor_start.elapsed();

    let parallel_naive_start = Instant::now();
    for _ in 0..repeats {
        let mut repeated = state.clone();
        compute_configured_force_report(
            &mut repeated,
            &force_options,
            mixed_force_options.as_ref(),
            None,
            topology.as_ref(),
            coulomb_options.as_ref(),
            true,
        )?;
    }
    let parallel_naive_elapsed = parallel_naive_start.elapsed();

    let parallel_neighbor_start = Instant::now();
    for _ in 0..repeats {
        let mut repeated = state.clone();
        compute_configured_force_report(
            &mut repeated,
            &force_options,
            mixed_force_options.as_ref(),
            Some(&neighbor_list),
            topology.as_ref(),
            coulomb_options.as_ref(),
            true,
        )?;
    }
    let parallel_neighbor_elapsed = parallel_neighbor_start.elapsed();

    let bond_count = topology.as_ref().map_or(0, |topology| topology.bonds.len());
    let angle_count = topology
        .as_ref()
        .map_or(0, |topology| topology.angles.len());
    let dihedral_count = topology
        .as_ref()
        .map_or(0, |topology| topology.dihedrals.len());
    let exclusion_count = topology
        .as_ref()
        .map_or(0, |topology| topology.exclusions.len());

    println!("Neighbor benchmark: {}", config_path.display());
    println!("Particles: {}", state.particle_count());
    println!("Boundary: {}", config.boundary_label());
    println!("Cutoff: {}", config.force.cutoff);
    println!("Skin: {}", neighbor_config.skin);
    match &topology {
        Some(topology) => println!("Topology: {}", topology.source_path.display()),
        None => println!("Topology: none"),
    }
    println!("Harmonic bonds: {bond_count}");
    println!("Harmonic angles: {angle_count}");
    println!("Periodic dihedrals: {dihedral_count}");
    println!("Non-bonded exclusions: {exclusion_count}");
    println!("Coulomb enabled: {}", coulomb_options.is_some());
    if let Some(coulomb_options) = &coulomb_options {
        println!("Coulomb cutoff: {}", coulomb_options.cutoff);
    }
    println!("Rayon threads: {}", rayon::current_num_threads());
    println!("Neighbor pairs: {}", neighbor_list.pairs().len());
    println!("Naive LJ interacting pairs: {}", naive_report.lj_pair_count);
    println!(
        "Neighbor LJ interacting pairs: {}",
        neighbor_report.lj_pair_count
    );
    println!(
        "Naive Coulomb interacting pairs: {}",
        naive_report.coulomb_pair_count
    );
    println!("Neighbor build: {:.6?}", build_elapsed);
    println!("Naive force x{repeats}: {:.6?}", naive_elapsed);
    println!("Neighbor force x{repeats}: {:.6?}", neighbor_elapsed);
    println!(
        "Parallel naive force x{repeats}: {:.6?}",
        parallel_naive_elapsed
    );
    println!(
        "Parallel neighbor force x{repeats}: {:.6?}",
        parallel_neighbor_elapsed
    );
    println!(
        "Neighbor speedup vs naive: {:.3}x",
        speedup_ratio(naive_elapsed, neighbor_elapsed)
    );
    println!(
        "Parallel naive speedup vs serial naive: {:.3}x",
        speedup_ratio(naive_elapsed, parallel_naive_elapsed)
    );
    println!(
        "Parallel neighbor speedup vs serial neighbor: {:.3}x",
        speedup_ratio(neighbor_elapsed, parallel_neighbor_elapsed)
    );
    println!("Neighbor energy delta: {:.6e}", energy_delta);
    println!("Neighbor max force delta: {:.6e}", serial_max_force_delta);
    println!(
        "Parallel naive energy delta: {:.6e}",
        parallel_naive_energy_delta
    );
    println!(
        "Parallel naive max force delta: {:.6e}",
        parallel_naive_max_force_delta
    );
    println!(
        "Parallel neighbor energy delta: {:.6e}",
        parallel_energy_delta
    );
    println!(
        "Parallel neighbor max force delta: {:.6e}",
        parallel_max_force_delta
    );

    Ok(())
}

fn analyze_command(run_dir: &Path, reference_frame: usize, pair: Option<&str>) -> Result<()> {
    let energy_path = run_dir.join("energy.csv");
    let trajectory_path = run_dir.join("trajectory.xyz");

    let energy_samples = read_energy_csv(&energy_path)
        .with_context(|| format!("failed to read {}", energy_path.display()))?;
    let raw_frames = read_xyz_trajectory(&trajectory_path)
        .with_context(|| format!("failed to read {}", trajectory_path.display()))?;
    let periodic_box = analysis_periodic_box(run_dir)?;
    let frames = match periodic_box {
        Some(simulation_box) => unwrap_periodic_frames(&raw_frames, simulation_box)
            .context("failed to unwrap periodic trajectory for analysis")?,
        None => raw_frames,
    };

    let rmsd_samples = rmsd_series(&frames, reference_frame)?;
    let pair_indices = match pair {
        Some(pair) => Some(parse_pair_indices(pair)?),
        None if frames.first().map_or(0, |frame| frame.atoms.len()) >= 2 => Some((0, 1)),
        None => None,
    };
    let pair_samples = pair_indices
        .map(|(atom_i, atom_j)| pair_distance_series(&frames, atom_i, atom_j))
        .transpose()?;

    let summary = summarize_analysis_with_metadata(
        &energy_samples,
        &frames,
        reference_frame,
        pair_samples.as_deref(),
        periodic_box.is_some(),
        periodic_box,
    )?;

    let summary_path = run_dir.join("analysis-summary.json");
    let rmsd_path = run_dir.join("analysis-rmsd.csv");
    let distance_path = run_dir.join("analysis-distance.csv");
    let report_path = run_dir.join("analysis-report.md");

    write_json(&summary_path, &summary)?;
    write_rmsd_csv(&rmsd_path, &rmsd_samples)?;
    if let Some(pair_samples) = &pair_samples {
        write_pair_distance_csv(&distance_path, pair_samples)?;
    }
    write_analysis_report(&report_path, &summary)?;
    update_manifest_after_analysis(run_dir)?;

    println!("Analysis completed.");
    println!("Run: {}", run_dir.display());
    println!("Frames: {}", summary.frame_count);
    println!("Atoms: {}", summary.atom_count);
    println!("Energy drift: {:.6e}", summary.energy.total_drift);
    println!("Final RMSD: {:.6}", summary.rmsd.final_value);
    println!("RMSD alignment: {}", summary.rmsd_alignment);
    println!("Periodic unwrap: {}", summary.periodic_unwrap);
    if let Some(pair_distance) = &summary.pair_distance {
        println!(
            "Pair distance {}-{} final: {:.6}",
            pair_distance.atom_i, pair_distance.atom_j, pair_distance.final_value
        );
    }
    println!("Output: {}", summary_path.display());
    println!("Output: {}", rmsd_path.display());
    if pair_samples.is_some() {
        println!("Output: {}", distance_path.display());
    }
    println!("Output: {}", report_path.display());

    Ok(())
}

fn analysis_periodic_box(run_dir: &Path) -> Result<Option<SimulationBox>> {
    let config_path = run_dir.join("config.toml");
    if !config_path.exists() {
        return Ok(None);
    }

    let config = load_config(&config_path)?;
    if config.simulation_box.periodic {
        Ok(Some(config.simulation_box.to_simulation_box()?))
    } else {
        Ok(None)
    }
}

fn report_command(run_dir: &Path) -> Result<()> {
    let summary_path = run_dir.join("summary.json");
    let manifest_path = run_dir.join("run-manifest.json");
    let report_path = run_dir.join("run-report.md");

    let summary = read_summary_json(&summary_path)
        .with_context(|| format!("failed to read {}", summary_path.display()))?;
    let analysis = read_optional_analysis_summary(run_dir)?;
    let mut manifest = if manifest_path.exists() {
        read_manifest_json(&manifest_path)
            .with_context(|| format!("failed to read {}", manifest_path.display()))?
    } else {
        manifest_from_summary(&summary)
    };
    manifest.report_file = Some("run-report.md".to_string());
    if let Some(checkpoint_file) = existing_checkpoint_file(run_dir) {
        manifest.checkpoint_file = Some(checkpoint_file.to_string());
        manifest.checkpoint_format = Some(checkpoint_format_from_file(checkpoint_file).to_string());
        push_unique(&mut manifest.outputs, checkpoint_file);
    }
    if analysis.is_some() {
        manifest.analysis_summary_file = Some("analysis-summary.json".to_string());
        push_unique(&mut manifest.outputs, "analysis-summary.json");
        push_unique(&mut manifest.outputs, "analysis-rmsd.csv");
        if run_dir.join("analysis-distance.csv").exists() {
            push_unique(&mut manifest.outputs, "analysis-distance.csv");
        }
        push_unique(&mut manifest.outputs, "analysis-report.md");
    }
    push_unique(&mut manifest.outputs, "run-report.md");
    push_unique(&mut manifest.outputs, "run-manifest.json");

    let report = render_markdown_report(&MarkdownReportInput {
        run_dir: &run_dir.display().to_string(),
        summary: &summary,
        manifest: Some(&manifest),
        analysis: analysis.as_ref(),
    });
    fs::write(&report_path, report)
        .with_context(|| format!("failed to write {}", report_path.display()))?;
    write_manifest_json(&manifest_path, &manifest)
        .with_context(|| format!("failed to write {}", manifest_path.display()))?;

    println!("Report completed.");
    println!("Run: {}", run_dir.display());
    println!("Output: {}", report_path.display());
    println!("Output: {}", manifest_path.display());

    Ok(())
}

fn resume_command(run_dir: &Path, additional_steps: Option<usize>) -> Result<()> {
    let config_path = run_dir.join("config.toml");
    let summary_path = run_dir.join("summary.json");

    let config = load_config(&config_path)?;
    config.validate()?;
    let mut summary = read_summary_json(&summary_path)
        .with_context(|| format!("failed to read {}", summary_path.display()))?;
    let (checkpoint_file, checkpoint_format) = checkpoint_info_for_resume(run_dir, &config)?;
    let resume_checkpoint = CheckpointSection {
        format: checkpoint_format.clone(),
    };
    resume_checkpoint.validate()?;
    let checkpoint_path = run_dir.join(&checkpoint_file);
    let checkpoint = read_checkpoint_by_format(&checkpoint_path, &checkpoint_format)
        .with_context(|| format!("failed to read {}", checkpoint_path.display()))?;
    if checkpoint.run_name != summary.run_name {
        bail!(
            "checkpoint run name {:?} does not match summary run name {:?}",
            checkpoint.run_name,
            summary.run_name
        );
    }

    let mut state = checkpoint.state;
    state.validate()?;
    let force_options = config.force_options()?;
    let mixed_force_options = config.mixed_lennard_jones_options(&state, &force_options)?;
    let coulomb_options = config.coulomb_options()?;
    let thermostat_options = config.thermostat_options();
    let topology = load_topology_for_resume(&summary, run_dir, state.particle_count())?;
    let mut neighbor_list = if config.neighbor.enabled {
        Some(NeighborList::build(
            &state,
            config.neighbor_list_config(&force_options)?,
        )?)
    } else {
        None
    };
    let current_potential_energy = compute_configured_forces(
        &mut state,
        &force_options,
        mixed_force_options.as_ref(),
        neighbor_list.as_ref(),
        topology.as_ref(),
        coulomb_options.as_ref(),
        config.execution.parallel,
    )?;
    let integrator = VelocityVerlet::new(config.simulation.dt)?;
    let target_step = match additional_steps {
        Some(steps) if steps > 0 => checkpoint.step + steps,
        Some(_) => bail!("--additional-steps must be greater than zero"),
        None if checkpoint.step < config.simulation.steps => config.simulation.steps,
        None => bail!(
            "checkpoint is already at step {}; use --additional-steps to extend the run",
            checkpoint.step
        ),
    };

    let trajectory_file = OpenOptions::new()
        .append(true)
        .open(run_dir.join("trajectory.xyz"))
        .with_context(|| {
            format!(
                "failed to open {}",
                run_dir.join("trajectory.xyz").display()
            )
        })?;
    let mut trajectory = BufWriter::new(trajectory_file);
    let energy_file = OpenOptions::new()
        .append(true)
        .open(run_dir.join("energy.csv"))
        .with_context(|| format!("failed to open {}", run_dir.join("energy.csv").display()))?;
    let mut energy_writer = EnergyCsvWriter::append(BufWriter::new(energy_file));

    println!("Resuming simulation...");
    println!("Run: {}", run_dir.display());
    println!("Checkpoint step: {}", checkpoint.step);
    println!("Target step: {}", target_step);
    println!("Particles: {}", state.particle_count());
    println!("Boundary: {}", config.boundary_label());
    println!("Ensemble: {}", config.ensemble_label());
    println!("Thermostat: {}", thermostat_options.label());
    if let Some(target_temperature) = thermostat_options.target_temperature() {
        println!("Thermostat target temperature: {}", target_temperature);
    }
    if let Some(tau) = thermostat_options.tau() {
        println!("Thermostat tau: {}", tau);
    }
    println!("Coulomb: {}", coulomb_options.is_some());
    println!("Parallel force: {}", config.execution.parallel);

    let mut final_sample = sample_energy(
        checkpoint.step,
        integrator.dt(),
        &state,
        current_potential_energy,
    );
    let progress_interval = ((target_step - checkpoint.step) / 10).max(1);

    for step in (checkpoint.step + 1)..=target_step {
        let force_rebuild = step % config.neighbor.rebuild_interval == 0;
        let potential_energy =
            integrator.step_with_force(&mut state, force_options.boundary, |state| {
                if let Some(neighbor_list) = neighbor_list.as_mut() {
                    if force_rebuild || neighbor_list.needs_rebuild(state)? {
                        neighbor_list.rebuild(state)?;
                    }
                }
                compute_configured_forces(
                    state,
                    &force_options,
                    mixed_force_options.as_ref(),
                    neighbor_list.as_ref(),
                    topology.as_ref(),
                    coulomb_options.as_ref(),
                    config.execution.parallel,
                )
                .map_err(IntegratorError::from)
            })?;
        apply_thermostat(&mut state, thermostat_options, integrator.dt())?;

        if should_output(step, target_step, config.simulation.output_interval) {
            final_sample = sample_energy(step, integrator.dt(), &state, potential_energy);
            ensure_finite_energy(&final_sample)?;
            write_xyz_frame(&mut trajectory, &state, step)?;
            energy_writer.write_sample(&final_sample)?;
            write_run_checkpoint(
                run_dir,
                &summary.run_name,
                final_sample.step,
                final_sample.time,
                potential_energy,
                &state,
                &resume_checkpoint,
            )?;
        }

        let elapsed = step - checkpoint.step;
        if elapsed % progress_interval == 0 || step == target_step {
            let sample = sample_energy(step, integrator.dt(), &state, potential_energy);
            println!(
                "Step {}/{}  E_kin={:.6}  E_pot={:.6}  E_total={:.6}  T={:.6}",
                step,
                target_step,
                sample.kinetic,
                sample.potential,
                sample.total,
                sample.temperature
            );
        }
    }

    trajectory.flush()?;
    energy_writer.flush()?;

    summary.steps = target_step;
    summary.final_energy = final_sample;
    summary.workflow = "dynamics".to_string();
    summary.ensemble = config.ensemble_label().to_string();
    summary.thermostat = thermostat_options.label().to_string();
    summary.thermostat_target_temperature = thermostat_options.target_temperature();
    summary.thermostat_tau = thermostat_options.tau();
    summary.minimization_step_size = None;
    summary.minimization_force_tolerance = None;
    summary.minimization_final_max_force = None;
    summary.bond_count = topology.as_ref().map_or(0, |topology| topology.bonds.len());
    summary.angle_count = topology
        .as_ref()
        .map_or(0, |topology| topology.angles.len());
    summary.dihedral_count = topology
        .as_ref()
        .map_or(0, |topology| topology.dihedrals.len());
    summary.excluded_pair_count = topology
        .as_ref()
        .map_or(0, |topology| topology.exclusions.len());
    push_unique(&mut summary.outputs, &checkpoint_file);
    write_summary_json(&summary_path, &summary)?;

    let manifest_path = run_dir.join("run-manifest.json");
    let mut manifest = if manifest_path.exists() {
        read_manifest_json(&manifest_path)
            .with_context(|| format!("failed to read {}", manifest_path.display()))?
    } else {
        manifest_from_summary(&summary)
    };
    manifest.status = "completed".to_string();
    manifest.completed_at_unix_seconds = Some(current_unix_seconds());
    manifest.checkpoint_file = Some(checkpoint_file.clone());
    manifest.checkpoint_format = Some(checkpoint_format);
    push_unique(&mut manifest.outputs, &checkpoint_file);
    write_manifest_json(&manifest_path, &manifest)?;
    report_command(run_dir)?;

    println!("Resume completed.");
    println!("Output: {}", run_dir.display());

    Ok(())
}

#[derive(Debug, Deserialize)]
struct RunConfig {
    #[serde(default)]
    project: ProjectSection,
    simulation: SimulationSection,
    #[serde(default)]
    input: Option<InputSection>,
    #[serde(default)]
    topology: Option<TopologySection>,
    #[serde(default)]
    system: SystemSection,
    #[serde(rename = "box")]
    simulation_box: BoxSection,
    force: ForceSection,
    #[serde(default)]
    force_field: ForceFieldSection,
    #[serde(default)]
    coulomb: CoulombSection,
    #[serde(default)]
    neighbor: NeighborSection,
    #[serde(default)]
    thermostat: ThermostatSection,
    #[serde(default)]
    minimization: MinimizationSection,
    #[serde(default)]
    checkpoint: CheckpointSection,
    #[serde(default)]
    execution: ExecutionSection,
    #[serde(default)]
    output: OutputSection,
}

#[derive(Debug, Default, Deserialize)]
struct ProjectSection {
    name: Option<String>,
    description: Option<String>,
}

impl ProjectSection {
    fn validate(&self) -> Result<()> {
        if self
            .name
            .as_deref()
            .is_some_and(|name| name.trim().is_empty())
        {
            bail!("project.name must not be empty when provided");
        }
        if self
            .description
            .as_deref()
            .is_some_and(|description| description.trim().is_empty())
        {
            bail!("project.description must not be empty when provided");
        }
        Ok(())
    }
}

impl RunConfig {
    fn validate(&self) -> Result<()> {
        self.project.validate()?;
        self.simulation.validate()?;
        self.system.validate(self.input.is_none())?;
        if let Some(input) = &self.input {
            input.validate()?;
        }
        if let Some(topology) = &self.topology {
            topology.validate()?;
        }
        self.coulomb.validate()?;
        self.force_field.validate()?;
        self.neighbor.validate()?;
        self.thermostat.validate(self.simulation.dt)?;
        self.minimization.validate()?;
        self.checkpoint.validate()?;
        self.execution.validate()?;
        let force_options = self.force_options()?;
        self.coulomb_options()?;
        if self.neighbor.enabled {
            self.neighbor_list_config(&force_options)?;
        }

        if self.input.is_none() {
            let grid = lattice_grid_size(self.system.particles);
            let spacing = self
                .system
                .lattice_spacing
                .unwrap_or(self.force.sigma * 1.35);
            let extent = spacing * (grid.saturating_sub(1) as f64);
            let simulation_box = self.simulation_box.to_simulation_box()?;

            if extent >= simulation_box.x
                || extent >= simulation_box.y
                || extent >= simulation_box.z
            {
                bail!(
                    "initial lattice extent {extent} does not fit in box ({}, {}, {})",
                    simulation_box.x,
                    simulation_box.y,
                    simulation_box.z
                );
            }
        }

        Ok(())
    }

    fn force_options(&self) -> Result<LennardJonesOptions> {
        let params = self.force.to_lennard_jones_params()?;
        let simulation_box = self.simulation_box.to_simulation_box()?;

        let options = if self.simulation_box.periodic {
            LennardJonesOptions::periodic(params, simulation_box, self.force.shift_potential)
        } else {
            LennardJonesOptions {
                params,
                boundary: BoundaryCondition::Open,
                shift_potential: self.force.shift_potential,
            }
        };

        options.validate()?;
        Ok(options)
    }

    fn mixed_lennard_jones_options(
        &self,
        state: &SystemState,
        force_options: &LennardJonesOptions,
    ) -> Result<Option<MixedLennardJonesOptions>> {
        self.force_field
            .mixed_lennard_jones_options(state, &self.force, force_options)
    }

    fn neighbor_list_config(
        &self,
        force_options: &LennardJonesOptions,
    ) -> Result<NeighborListConfig> {
        let boundary = match force_options.boundary {
            BoundaryCondition::Open => NeighborBoundary::Open,
            BoundaryCondition::Periodic(simulation_box) => {
                NeighborBoundary::Periodic(simulation_box)
            }
        };
        let config = NeighborListConfig {
            cutoff: force_options.params.cutoff,
            skin: self.neighbor.skin,
            boundary,
        };
        config.validate()?;
        Ok(config)
    }

    fn boundary_label(&self) -> &'static str {
        if self.simulation_box.periodic {
            "periodic"
        } else {
            "open"
        }
    }

    fn coulomb_options(&self) -> Result<Option<CoulombOptions>> {
        if !self.coulomb.enabled {
            return Ok(None);
        }

        let cutoff = self.coulomb.cutoff.unwrap_or(self.force.cutoff);
        let simulation_box = self.simulation_box.to_simulation_box()?;
        let options = if self.simulation_box.periodic {
            CoulombOptions::periodic(
                self.coulomb.constant,
                cutoff,
                simulation_box,
                self.coulomb.shift_potential,
            )
        } else {
            CoulombOptions::open(self.coulomb.constant, cutoff, self.coulomb.shift_potential)
        };
        options.validate()?;
        Ok(Some(options))
    }

    fn thermostat_options(&self) -> ThermostatOptions {
        self.thermostat.to_options(self.simulation.temperature)
    }

    fn ensemble_label(&self) -> &'static str {
        match self.thermostat_options() {
            ThermostatOptions::None => "NVE",
            ThermostatOptions::Berendsen { .. } => "NVT",
        }
    }
}

#[derive(Debug, Deserialize)]
struct InputSection {
    path: PathBuf,
    #[serde(default = "default_input_format")]
    format: String,
    #[serde(default)]
    frame: usize,
}

impl InputSection {
    fn validate(&self) -> Result<()> {
        if self.path.as_os_str().is_empty() {
            bail!("input.path must not be empty");
        }
        let format = self.normalized_format();
        if format != "xyz" && format != "pdb" {
            bail!(
                "unsupported input.format {:?}; expected \"xyz\" or \"pdb\"",
                self.format
            );
        }
        Ok(())
    }

    fn normalized_format(&self) -> String {
        self.format.trim().to_ascii_lowercase()
    }
}

#[derive(Debug, Deserialize)]
struct TopologySection {
    path: PathBuf,
}

impl TopologySection {
    fn validate(&self) -> Result<()> {
        if self.path.as_os_str().is_empty() {
            bail!("topology.path must not be empty");
        }
        Ok(())
    }
}

#[derive(Debug, Deserialize)]
struct CoulombSection {
    #[serde(default)]
    enabled: bool,
    #[serde(default = "default_coulomb_constant")]
    constant: f64,
    cutoff: Option<f64>,
    #[serde(default)]
    shift_potential: bool,
}

impl Default for CoulombSection {
    fn default() -> Self {
        Self {
            enabled: false,
            constant: default_coulomb_constant(),
            cutoff: None,
            shift_potential: false,
        }
    }
}

impl CoulombSection {
    fn validate(&self) -> Result<()> {
        if !self.constant.is_finite() {
            bail!("coulomb.constant must be finite");
        }
        if let Some(cutoff) = self.cutoff {
            if !cutoff.is_finite() || cutoff <= 0.0 {
                bail!("coulomb.cutoff must be positive and finite");
            }
        }
        Ok(())
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
enum ThermostatOptions {
    None,
    Berendsen { target_temperature: f64, tau: f64 },
}

impl ThermostatOptions {
    fn label(self) -> &'static str {
        match self {
            ThermostatOptions::None => "none",
            ThermostatOptions::Berendsen { .. } => "berendsen",
        }
    }

    fn target_temperature(self) -> Option<f64> {
        match self {
            ThermostatOptions::None => None,
            ThermostatOptions::Berendsen {
                target_temperature, ..
            } => Some(target_temperature),
        }
    }

    fn tau(self) -> Option<f64> {
        match self {
            ThermostatOptions::None => None,
            ThermostatOptions::Berendsen { tau, .. } => Some(tau),
        }
    }
}

#[derive(Debug, Deserialize)]
struct ThermostatSection {
    #[serde(rename = "type", default = "default_thermostat_type")]
    kind: String,
    target_temperature: Option<f64>,
    #[serde(default = "default_thermostat_tau")]
    tau: f64,
}

impl Default for ThermostatSection {
    fn default() -> Self {
        Self {
            kind: default_thermostat_type(),
            target_temperature: None,
            tau: default_thermostat_tau(),
        }
    }
}

impl ThermostatSection {
    fn validate(&self, dt: f64) -> Result<()> {
        let kind = self.normalized_kind();
        if kind != "none" && kind != "berendsen" {
            bail!(
                "unsupported thermostat.type {:?}; expected \"none\" or \"berendsen\"",
                self.kind
            );
        }
        if let Some(target_temperature) = self.target_temperature {
            if !target_temperature.is_finite() || target_temperature < 0.0 {
                bail!("thermostat.target_temperature must be finite and non-negative");
            }
        }
        if kind == "berendsen" {
            if !self.tau.is_finite() || self.tau <= 0.0 {
                bail!("thermostat.tau must be positive and finite");
            }
            if self.tau < dt {
                bail!("thermostat.tau must be greater than or equal to simulation.dt");
            }
        }
        Ok(())
    }

    fn to_options(&self, simulation_temperature: f64) -> ThermostatOptions {
        match self.normalized_kind().as_str() {
            "berendsen" => ThermostatOptions::Berendsen {
                target_temperature: self.target_temperature.unwrap_or(simulation_temperature),
                tau: self.tau,
            },
            _ => ThermostatOptions::None,
        }
    }

    fn normalized_kind(&self) -> String {
        self.kind.trim().to_ascii_lowercase()
    }
}

#[derive(Debug, Deserialize)]
struct MinimizationSection {
    #[serde(default = "default_minimization_steps")]
    steps: usize,
    #[serde(default = "default_minimization_output_interval")]
    output_interval: usize,
    #[serde(default = "default_minimization_step_size")]
    step_size: f64,
    #[serde(default = "default_minimization_max_displacement")]
    max_displacement: f64,
    #[serde(default = "default_minimization_force_tolerance")]
    force_tolerance: f64,
    #[serde(default = "default_minimization_max_backtracks")]
    max_backtracks: usize,
}

impl Default for MinimizationSection {
    fn default() -> Self {
        Self {
            steps: default_minimization_steps(),
            output_interval: default_minimization_output_interval(),
            step_size: default_minimization_step_size(),
            max_displacement: default_minimization_max_displacement(),
            force_tolerance: default_minimization_force_tolerance(),
            max_backtracks: default_minimization_max_backtracks(),
        }
    }
}

impl MinimizationSection {
    fn validate(&self) -> Result<()> {
        if self.steps == 0 {
            bail!("minimization.steps must be greater than zero");
        }
        if self.output_interval == 0 {
            bail!("minimization.output_interval must be greater than zero");
        }
        if !self.step_size.is_finite() || self.step_size <= 0.0 {
            bail!("minimization.step_size must be positive and finite");
        }
        if !self.max_displacement.is_finite() || self.max_displacement <= 0.0 {
            bail!("minimization.max_displacement must be positive and finite");
        }
        if !self.force_tolerance.is_finite() || self.force_tolerance < 0.0 {
            bail!("minimization.force_tolerance must be finite and non-negative");
        }
        if self.max_backtracks == 0 {
            bail!("minimization.max_backtracks must be greater than zero");
        }
        Ok(())
    }
}

#[derive(Debug, Deserialize)]
struct CheckpointSection {
    #[serde(default = "default_checkpoint_format")]
    format: String,
}

impl Default for CheckpointSection {
    fn default() -> Self {
        Self {
            format: default_checkpoint_format(),
        }
    }
}

impl CheckpointSection {
    fn validate(&self) -> Result<()> {
        let format = self.normalized_format();
        if format != "json" && format != "binary" {
            bail!(
                "unsupported checkpoint.format {:?}; expected \"json\" or \"binary\"",
                self.format
            );
        }
        Ok(())
    }

    fn normalized_format(&self) -> String {
        self.format.trim().to_ascii_lowercase()
    }

    fn file_name(&self) -> &'static str {
        match self.normalized_format().as_str() {
            "binary" => "checkpoint.bin",
            _ => "checkpoint.json",
        }
    }
}

#[derive(Debug, Deserialize)]
struct ExecutionSection {
    #[serde(default)]
    parallel: bool,
}

impl Default for ExecutionSection {
    fn default() -> Self {
        Self { parallel: false }
    }
}

impl ExecutionSection {
    fn validate(&self) -> Result<()> {
        Ok(())
    }
}

#[derive(Debug, Deserialize)]
struct NeighborSection {
    #[serde(default)]
    enabled: bool,
    #[serde(default = "default_neighbor_skin")]
    skin: f64,
    #[serde(default = "default_neighbor_rebuild_interval")]
    rebuild_interval: usize,
}

impl Default for NeighborSection {
    fn default() -> Self {
        Self {
            enabled: false,
            skin: default_neighbor_skin(),
            rebuild_interval: default_neighbor_rebuild_interval(),
        }
    }
}

impl NeighborSection {
    fn validate(&self) -> Result<()> {
        if !self.skin.is_finite() || self.skin < 0.0 {
            bail!("neighbor.skin must be finite and non-negative");
        }
        if self.rebuild_interval == 0 {
            bail!("neighbor.rebuild_interval must be greater than zero");
        }
        Ok(())
    }
}

#[derive(Debug, Deserialize)]
struct SimulationSection {
    steps: usize,
    dt: f64,
    output_interval: usize,
    temperature: f64,
}

impl SimulationSection {
    fn validate(&self) -> Result<()> {
        if self.steps == 0 {
            bail!("simulation.steps must be greater than zero");
        }
        if self.output_interval == 0 {
            bail!("simulation.output_interval must be greater than zero");
        }
        if !self.dt.is_finite() || self.dt <= 0.0 {
            bail!("simulation.dt must be positive and finite");
        }
        if !self.temperature.is_finite() || self.temperature < 0.0 {
            bail!("simulation.temperature must be finite and non-negative");
        }
        if self.dt > 0.01 {
            bail!("simulation.dt is too large for the current conservative Phase 1 runner");
        }
        Ok(())
    }
}

#[derive(Debug, Deserialize)]
struct SystemSection {
    #[serde(default = "default_particles")]
    particles: usize,
    #[serde(default = "default_mass")]
    mass: f64,
    #[serde(default = "default_element")]
    element: String,
    #[serde(default = "default_seed")]
    seed: u64,
    lattice_spacing: Option<f64>,
}

impl Default for SystemSection {
    fn default() -> Self {
        Self {
            particles: default_particles(),
            mass: default_mass(),
            element: default_element(),
            seed: default_seed(),
            lattice_spacing: None,
        }
    }
}

impl SystemSection {
    fn validate(&self, require_generated_particles: bool) -> Result<()> {
        if require_generated_particles && self.particles < 2 {
            bail!("system.particles must be at least 2");
        }
        if !self.mass.is_finite() || self.mass <= 0.0 {
            bail!("system.mass must be positive and finite");
        }
        if require_generated_particles && self.element.trim().is_empty() {
            bail!("system.element must not be empty");
        }
        if let Some(spacing) = self.lattice_spacing {
            if !spacing.is_finite() || spacing <= 0.0 {
                bail!("system.lattice_spacing must be positive and finite");
            }
        }
        Ok(())
    }
}

#[derive(Debug, Deserialize)]
struct BoxSection {
    x: f64,
    y: f64,
    z: f64,
    #[serde(default)]
    periodic: bool,
}

impl BoxSection {
    fn to_simulation_box(&self) -> Result<SimulationBox> {
        SimulationBox::new(self.x, self.y, self.z).context("invalid [box] section")
    }
}

#[derive(Debug, Deserialize)]
struct ForceSection {
    #[serde(rename = "type")]
    kind: String,
    epsilon: f64,
    sigma: f64,
    cutoff: f64,
    #[serde(default)]
    shift_potential: bool,
}

impl ForceSection {
    fn to_lennard_jones_params(&self) -> Result<LennardJonesParams> {
        if self.kind != "lennard-jones" {
            bail!(
                "unsupported force.type {:?}; only \"lennard-jones\" is implemented",
                self.kind
            );
        }

        let params = LennardJonesParams {
            epsilon: self.epsilon,
            sigma: self.sigma,
            cutoff: self.cutoff,
        };
        params.validate()?;
        Ok(params)
    }
}

#[derive(Debug, Deserialize)]
struct ForceFieldSection {
    #[serde(default = "default_mixing_rule")]
    mixing_rule: String,
    #[serde(default)]
    types: BTreeMap<String, ForceFieldTypeSection>,
}

impl Default for ForceFieldSection {
    fn default() -> Self {
        Self {
            mixing_rule: default_mixing_rule(),
            types: BTreeMap::new(),
        }
    }
}

impl ForceFieldSection {
    fn validate(&self) -> Result<()> {
        if self.normalized_mixing_rule() != "lorentz-berthelot" {
            bail!(
                "unsupported force_field.mixing_rule {:?}; expected \"lorentz-berthelot\"",
                self.mixing_rule
            );
        }
        for (name, params) in &self.types {
            if name.trim().is_empty() {
                bail!("force_field type names must not be empty");
            }
            params
                .validate()
                .with_context(|| format!("invalid force_field.types.{name}"))?;
        }
        Ok(())
    }

    fn type_for_element(&self, element: &str) -> Option<&ForceFieldTypeSection> {
        self.types.get(element.trim())
    }

    fn mixed_lennard_jones_options(
        &self,
        state: &SystemState,
        force: &ForceSection,
        force_options: &LennardJonesOptions,
    ) -> Result<Option<MixedLennardJonesOptions>> {
        if !self.has_lj_parameters() {
            return Ok(None);
        }

        let particle_params = state
            .element
            .iter()
            .map(|element| {
                let params = self.type_for_element(element);
                LennardJonesParticleParams {
                    epsilon: params
                        .and_then(|params| params.epsilon)
                        .unwrap_or(force.epsilon),
                    sigma: params
                        .and_then(|params| params.sigma)
                        .unwrap_or(force.sigma),
                }
            })
            .collect();
        let options = MixedLennardJonesOptions {
            particle_params,
            cutoff: force.cutoff,
            boundary: force_options.boundary,
            shift_potential: force_options.shift_potential,
            mixing_rule: LennardJonesMixingRule::LorentzBerthelot,
        };
        options.validate(state.particle_count())?;
        Ok(Some(options))
    }

    fn has_lj_parameters(&self) -> bool {
        self.types
            .values()
            .any(|params| params.sigma.is_some() || params.epsilon.is_some())
    }

    fn normalized_mixing_rule(&self) -> String {
        self.mixing_rule.trim().to_ascii_lowercase()
    }
}

#[derive(Debug, Default, Deserialize)]
struct ForceFieldTypeSection {
    mass: Option<f64>,
    sigma: Option<f64>,
    epsilon: Option<f64>,
    charge: Option<f64>,
}

impl ForceFieldTypeSection {
    fn validate(&self) -> Result<()> {
        if let Some(mass) = self.mass {
            if !mass.is_finite() || mass <= 0.0 {
                bail!("mass must be positive and finite");
            }
        }
        if let Some(sigma) = self.sigma {
            if !sigma.is_finite() || sigma <= 0.0 {
                bail!("sigma must be positive and finite");
            }
        }
        if let Some(epsilon) = self.epsilon {
            if !epsilon.is_finite() || epsilon <= 0.0 {
                bail!("epsilon must be positive and finite");
            }
        }
        if let Some(charge) = self.charge {
            if !charge.is_finite() {
                bail!("charge must be finite");
            }
        }
        Ok(())
    }
}

#[derive(Debug, Deserialize)]
struct OutputSection {
    #[serde(default = "default_output_directory")]
    directory: PathBuf,
    name: Option<String>,
}

impl Default for OutputSection {
    fn default() -> Self {
        Self {
            directory: default_output_directory(),
            name: None,
        }
    }
}

#[derive(Debug, Deserialize)]
struct TopologyFile {
    charges: Option<Vec<f64>>,
    #[serde(default)]
    bonds: Vec<HarmonicBond>,
    #[serde(default)]
    angles: Vec<HarmonicAngle>,
    #[serde(default)]
    dihedrals: Vec<PeriodicDihedral>,
    #[serde(default)]
    exclusions: Vec<ExcludedPair>,
}

#[derive(Debug)]
struct LoadedTopology {
    source_path: PathBuf,
    charges: Option<Vec<f64>>,
    bonds: Vec<HarmonicBond>,
    angles: Vec<HarmonicAngle>,
    dihedrals: Vec<PeriodicDihedral>,
    exclusions: Vec<(usize, usize)>,
}

fn load_config(path: &Path) -> Result<RunConfig> {
    let source = fs::read_to_string(path)
        .with_context(|| format!("failed to read config {}", path.display()))?;
    toml::from_str(&source).with_context(|| format!("failed to parse config {}", path.display()))
}

fn prepare_initial_state_and_topology(
    config: &RunConfig,
    config_dir: Option<&Path>,
) -> Result<(SystemState, Option<LoadedTopology>)> {
    let mut state = build_initial_system(config, config_dir)?;
    apply_force_field_to_state(config, &mut state)?;
    let topology = load_topology(config, config_dir, state.particle_count())?;
    Ok((state, topology))
}

fn load_topology(
    config: &RunConfig,
    config_dir: Option<&Path>,
    particle_count: usize,
) -> Result<Option<LoadedTopology>> {
    let Some(topology_section) = &config.topology else {
        return Ok(None);
    };

    let source_path = resolve_config_relative_path(&topology_section.path, config_dir);
    load_topology_from_path(source_path, particle_count).map(Some)
}

fn load_topology_from_path(
    source_path: impl Into<PathBuf>,
    particle_count: usize,
) -> Result<LoadedTopology> {
    let source_path = source_path.into();
    let source = fs::read_to_string(&source_path)
        .with_context(|| format!("failed to read topology {}", source_path.display()))?;
    let topology: TopologyFile = toml::from_str(&source)
        .with_context(|| format!("failed to parse topology {}", source_path.display()))?;

    if let Some(charges) = &topology.charges {
        if charges.len() != particle_count {
            bail!(
                "topology charges length {} does not match particle count {particle_count}",
                charges.len()
            );
        }
        for (index, charge) in charges.iter().copied().enumerate() {
            if !charge.is_finite() {
                bail!("topology charge[{index}] must be finite");
            }
        }
    }

    for bond in &topology.bonds {
        bond.validate(particle_count)
            .with_context(|| format!("invalid topology {}", source_path.display()))?;
    }
    for angle in &topology.angles {
        angle
            .validate(particle_count)
            .with_context(|| format!("invalid topology {}", source_path.display()))?;
    }
    for dihedral in &topology.dihedrals {
        dihedral
            .validate(particle_count)
            .with_context(|| format!("invalid topology {}", source_path.display()))?;
    }

    let mut exclusions = topology.exclusions;
    exclusions.extend(topology.bonds.iter().map(|bond| ExcludedPair {
        i: bond.i,
        j: bond.j,
    }));
    let exclusions = normalize_excluded_pairs(&exclusions, particle_count)
        .with_context(|| format!("invalid topology {}", source_path.display()))?;

    Ok(LoadedTopology {
        source_path,
        charges: topology.charges,
        bonds: topology.bonds,
        angles: topology.angles,
        dihedrals: topology.dihedrals,
        exclusions,
    })
}

fn apply_topology_to_state(
    state: &mut SystemState,
    topology: Option<&LoadedTopology>,
) -> Result<()> {
    let Some(topology) = topology else {
        return Ok(());
    };

    if let Some(charges) = &topology.charges {
        state.charge.clone_from(charges);
    }
    state.validate()?;
    Ok(())
}

fn apply_force_field_to_state(config: &RunConfig, state: &mut SystemState) -> Result<()> {
    for index in 0..state.particle_count() {
        if let Some(params) = config.force_field.type_for_element(&state.element[index]) {
            if let Some(mass) = params.mass {
                state.mass[index] = mass;
            }
            if let Some(charge) = params.charge {
                state.charge[index] = charge;
            }
        }
    }
    state.validate()?;
    Ok(())
}

fn copied_input_file_name(input: &InputSection) -> String {
    format!("input.{}", input.normalized_format())
}

fn write_optional_atom_metadata_file(
    run_dir: &Path,
    state: &SystemState,
) -> Result<Option<String>> {
    if !state.has_molecular_metadata() {
        return Ok(None);
    }

    let file_name = "atom-metadata.csv";
    let path = run_dir.join(file_name);
    let file =
        File::create(&path).with_context(|| format!("failed to create {}", path.display()))?;
    let mut writer = BufWriter::new(file);
    write_atom_metadata_csv(&mut writer, state)
        .with_context(|| format!("failed to write {}", path.display()))?;
    writer.flush()?;
    Ok(Some(file_name.to_string()))
}

fn build_initial_system(config: &RunConfig, config_dir: Option<&Path>) -> Result<SystemState> {
    let simulation_box = config.simulation_box.to_simulation_box()?;

    if let Some(input) = &config.input {
        return build_initial_system_from_input(config, input, config_dir, &simulation_box);
    }

    build_initial_lattice_system(config, &simulation_box)
}

fn build_initial_lattice_system(
    config: &RunConfig,
    simulation_box: &SimulationBox,
) -> Result<SystemState> {
    let mut state = SystemState::new(config.system.particles);
    let spacing = config
        .system
        .lattice_spacing
        .unwrap_or(config.force.sigma * 1.35);
    let grid = lattice_grid_size(config.system.particles);
    let extent = spacing * (grid.saturating_sub(1) as f64);
    let origin_x = 0.5 * (simulation_box.x - extent);
    let origin_y = 0.5 * (simulation_box.y - extent);
    let origin_z = 0.5 * (simulation_box.z - extent);

    for i in 0..state.particle_count() {
        let ix = i % grid;
        let iy = (i / grid) % grid;
        let iz = i / (grid * grid);

        state.x[i] = origin_x + ix as f64 * spacing;
        state.y[i] = origin_y + iy as f64 * spacing;
        state.z[i] = origin_z + iz as f64 * spacing;
        state.mass[i] = config.system.mass;
        state.element[i] = config.system.element.clone();
    }

    initialize_velocities(
        &mut state,
        config.system.seed,
        config.simulation.temperature,
    )?;
    state.validate()?;
    Ok(state)
}

fn build_initial_system_from_input(
    config: &RunConfig,
    input: &InputSection,
    config_dir: Option<&Path>,
    simulation_box: &SimulationBox,
) -> Result<SystemState> {
    let input_path = resolve_config_relative_path(&input.path, config_dir);
    match input.normalized_format().as_str() {
        "xyz" => build_initial_system_from_xyz_input(config, input, &input_path, simulation_box),
        "pdb" => build_initial_system_from_pdb_input(config, input, &input_path, simulation_box),
        _ => unreachable!("input format is validated before building the initial system"),
    }
}

fn build_initial_system_from_xyz_input(
    config: &RunConfig,
    input: &InputSection,
    input_path: &Path,
    simulation_box: &SimulationBox,
) -> Result<SystemState> {
    let frames = read_xyz_coordinate_frames(input_path)
        .with_context(|| format!("failed to read input XYZ {}", input_path.display()))?;

    let frame = frames.get(input.frame).ok_or_else(|| {
        anyhow::anyhow!(
            "input XYZ frame {} is out of range; file {} has {} frame(s)",
            input.frame,
            input_path.display(),
            frames.len()
        )
    })?;
    if frame.atoms.len() < 2 {
        bail!(
            "input XYZ frame {} must contain at least 2 atoms",
            input.frame
        );
    }

    let mut state = SystemState::new(frame.atoms.len());
    for (index, atom) in frame.atoms.iter().enumerate() {
        if atom.element.trim().is_empty() {
            bail!("input XYZ atom {index} has an empty element");
        }
        state.element[index] = atom.element.clone();
        state.x[index] = atom.x;
        state.y[index] = atom.y;
        state.z[index] = atom.z;
        state.mass[index] = config.system.mass;
    }

    if config.simulation_box.periodic {
        state.wrap_positions(simulation_box);
    }

    initialize_velocities(
        &mut state,
        config.system.seed,
        config.simulation.temperature,
    )?;
    state.validate()?;
    Ok(state)
}

fn build_initial_system_from_pdb_input(
    config: &RunConfig,
    input: &InputSection,
    input_path: &Path,
    simulation_box: &SimulationBox,
) -> Result<SystemState> {
    let frames = read_pdb_coordinate_frames(input_path)
        .with_context(|| format!("failed to read input PDB {}", input_path.display()))?;
    let frame = frames.get(input.frame).ok_or_else(|| {
        anyhow::anyhow!(
            "input PDB frame {} is out of range; file {} has {} frame(s)",
            input.frame,
            input_path.display(),
            frames.len()
        )
    })?;
    pdb_frame_to_state(config, input.frame, frame, simulation_box)
}

fn pdb_frame_to_state(
    config: &RunConfig,
    frame_index: usize,
    frame: &PdbFrame,
    simulation_box: &SimulationBox,
) -> Result<SystemState> {
    if frame.atoms.len() < 2 {
        bail!("input PDB frame {frame_index} must contain at least 2 atoms");
    }

    let mut state = SystemState::new(frame.atoms.len());
    for (index, atom) in frame.atoms.iter().enumerate() {
        if atom.element.trim().is_empty() {
            bail!("input PDB atom {index} has an empty element");
        }
        state.element[index] = atom.element.clone();
        state.atom_name[index] = Some(atom.atom_name.clone());
        state.residue_name[index] = atom.residue_name.clone();
        state.residue_id[index] = atom.residue_id;
        state.chain_id[index] = atom.chain_id.clone();
        state.x[index] = atom.x;
        state.y[index] = atom.y;
        state.z[index] = atom.z;
        state.mass[index] = config.system.mass;
    }

    if config.simulation_box.periodic {
        state.wrap_positions(simulation_box);
    }

    initialize_velocities(
        &mut state,
        config.system.seed,
        config.simulation.temperature,
    )?;
    state.validate()?;
    Ok(state)
}

fn initialize_velocities(
    state: &mut SystemState,
    seed: u64,
    target_temperature: f64,
) -> Result<()> {
    let mut rng = ChaCha8Rng::seed_from_u64(seed);
    for i in 0..state.particle_count() {
        state.vx[i] = rng.gen_range(-1.0..1.0);
        state.vy[i] = rng.gen_range(-1.0..1.0);
        state.vz[i] = rng.gen_range(-1.0..1.0);
    }

    state.remove_center_of_mass_velocity();
    state.rescale_temperature(target_temperature)?;
    Ok(())
}

fn sample_energy(step: usize, dt: f64, state: &SystemState, potential: f64) -> EnergySample {
    EnergySample::new(
        step,
        step as f64 * dt,
        state.kinetic_energy(),
        potential,
        state.temperature(),
    )
}

fn minimization_energy_sample(step: usize, step_size: f64, potential: f64) -> EnergySample {
    EnergySample::new(step, step as f64 * step_size, 0.0, potential, 0.0)
}

fn max_force_norm(state: &SystemState) -> f64 {
    let mut max_force: f64 = 0.0;
    for i in 0..state.particle_count() {
        let force = state.fx[i].hypot(state.fy[i]).hypot(state.fz[i]);
        max_force = max_force.max(force);
    }
    max_force
}

fn minimization_step_scale(max_force: f64, config: &MinimizationSection) -> f64 {
    if max_force > 0.0 {
        config.step_size.min(config.max_displacement / max_force)
    } else {
        config.step_size
    }
}

fn displace_along_forces(
    state: &mut SystemState,
    step_scale: f64,
    boundary: BoundaryCondition,
) -> Result<()> {
    if !step_scale.is_finite() || step_scale <= 0.0 {
        bail!("minimization step scale must be positive and finite");
    }

    for i in 0..state.particle_count() {
        state.x[i] += step_scale * state.fx[i];
        state.y[i] += step_scale * state.fy[i];
        state.z[i] += step_scale * state.fz[i];
    }
    if let BoundaryCondition::Periodic(simulation_box) = boundary {
        state.wrap_positions(&simulation_box);
    }
    state.vx.fill(0.0);
    state.vy.fill(0.0);
    state.vz.fill(0.0);
    state.validate()?;
    Ok(())
}

fn ensure_finite_minimization_state(
    step: usize,
    potential: f64,
    max_force: f64,
    state: &SystemState,
) -> Result<()> {
    if !potential.is_finite() {
        bail!("minimization produced non-finite potential at step {step}");
    }
    if !max_force.is_finite() {
        bail!("minimization produced non-finite max force at step {step}");
    }
    state
        .validate()
        .with_context(|| format!("minimization produced invalid state at step {step}"))?;
    Ok(())
}

fn write_minimization_history_header<W: Write>(writer: &mut W) -> Result<()> {
    writeln!(writer, "step,potential,max_force,step_scale")?;
    Ok(())
}

fn write_minimization_history_sample<W: Write>(
    writer: &mut W,
    step: usize,
    potential: f64,
    max_force: f64,
    step_scale: f64,
) -> Result<()> {
    writeln!(
        writer,
        "{},{:.10},{:.10},{:.10}",
        step, potential, max_force, step_scale
    )?;
    Ok(())
}

fn apply_thermostat(state: &mut SystemState, options: ThermostatOptions, dt: f64) -> Result<()> {
    let ThermostatOptions::Berendsen {
        target_temperature,
        tau,
    } = options
    else {
        return Ok(());
    };

    if target_temperature == 0.0 {
        state.rescale_temperature(0.0)?;
        return Ok(());
    }

    let current_temperature = state.temperature();
    if current_temperature <= 0.0 || !current_temperature.is_finite() {
        bail!("berendsen thermostat cannot heat a zero-temperature state");
    }

    let coupling = dt / tau;
    let scale_sq = 1.0 + coupling * (target_temperature / current_temperature - 1.0);
    if !scale_sq.is_finite() || scale_sq < 0.0 {
        bail!("berendsen thermostat produced an invalid velocity scale");
    }
    scale_velocities(state, scale_sq.sqrt());
    state.validate()?;
    Ok(())
}

fn scale_velocities(state: &mut SystemState, scale: f64) {
    for i in 0..state.particle_count() {
        state.vx[i] *= scale;
        state.vy[i] *= scale;
        state.vz[i] *= scale;
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
struct ConfiguredForceReport {
    potential_energy: f64,
    lj_pair_count: usize,
    coulomb_pair_count: usize,
}

fn compute_configured_forces(
    state: &mut SystemState,
    force_options: &LennardJonesOptions,
    mixed_force_options: Option<&MixedLennardJonesOptions>,
    neighbor_list: Option<&NeighborList>,
    topology: Option<&LoadedTopology>,
    coulomb_options: Option<&CoulombOptions>,
    parallel: bool,
) -> std::result::Result<f64, ForceError> {
    Ok(compute_configured_force_report(
        state,
        force_options,
        mixed_force_options,
        neighbor_list,
        topology,
        coulomb_options,
        parallel,
    )?
    .potential_energy)
}

fn compute_configured_force_report(
    state: &mut SystemState,
    force_options: &LennardJonesOptions,
    mixed_force_options: Option<&MixedLennardJonesOptions>,
    neighbor_list: Option<&NeighborList>,
    topology: Option<&LoadedTopology>,
    coulomb_options: Option<&CoulombOptions>,
    parallel: bool,
) -> std::result::Result<ConfiguredForceReport, ForceError> {
    let excluded_pairs = topology.map_or(&[][..], |topology| topology.exclusions.as_slice());
    let lj_report = if let Some(mixed_force_options) = mixed_force_options {
        if let Some(neighbor_list) = neighbor_list {
            if parallel {
                compute_lennard_jones_forces_parallel_with_mixed_neighbor_list_and_exclusions(
                    state,
                    mixed_force_options,
                    neighbor_list,
                    excluded_pairs,
                )?
            } else {
                compute_lennard_jones_forces_with_mixed_neighbor_list_and_exclusions(
                    state,
                    mixed_force_options,
                    neighbor_list,
                    excluded_pairs,
                )?
            }
        } else if parallel {
            compute_lennard_jones_forces_parallel_with_mixed_options_and_exclusions(
                state,
                mixed_force_options,
                excluded_pairs,
            )?
        } else {
            compute_lennard_jones_forces_with_mixed_options_and_exclusions(
                state,
                mixed_force_options,
                excluded_pairs,
            )?
        }
    } else if let Some(neighbor_list) = neighbor_list {
        if parallel {
            compute_lennard_jones_forces_parallel_with_neighbor_list_and_exclusions(
                state,
                force_options,
                neighbor_list,
                excluded_pairs,
            )?
        } else {
            compute_lennard_jones_forces_with_neighbor_list_and_exclusions(
                state,
                force_options,
                neighbor_list,
                excluded_pairs,
            )?
        }
    } else if parallel {
        compute_lennard_jones_forces_parallel_with_options_and_exclusions(
            state,
            force_options,
            excluded_pairs,
        )?
    } else {
        compute_lennard_jones_forces_with_options_and_exclusions(
            state,
            force_options,
            excluded_pairs,
        )?
    };

    let mut potential_energy = lj_report.potential_energy;
    let mut coulomb_pair_count = 0;
    if let Some(coulomb_options) = coulomb_options {
        let coulomb_report =
            add_coulomb_forces_with_options_and_exclusions(state, coulomb_options, excluded_pairs)?;
        potential_energy += coulomb_report.potential_energy;
        coulomb_pair_count = coulomb_report.pair_count;
    }
    if let Some(topology) = topology {
        potential_energy +=
            add_harmonic_bond_forces(state, &topology.bonds, force_options.boundary)?
                .potential_energy;
        potential_energy +=
            add_harmonic_angle_forces(state, &topology.angles, force_options.boundary)?
                .potential_energy;
        potential_energy +=
            add_periodic_dihedral_forces(state, &topology.dihedrals, force_options.boundary)?
                .potential_energy;
    }

    Ok(ConfiguredForceReport {
        potential_energy,
        lj_pair_count: lj_report.pair_count,
        coulomb_pair_count,
    })
}

fn ensure_finite_energy(sample: &EnergySample) -> Result<()> {
    if sample.kinetic.is_finite()
        && sample.potential.is_finite()
        && sample.total.is_finite()
        && sample.temperature.is_finite()
    {
        Ok(())
    } else {
        bail!(
            "simulation produced non-finite energy at step {}",
            sample.step
        );
    }
}

fn write_json<T: Serialize>(path: &Path, value: &T) -> Result<()> {
    let file =
        File::create(path).with_context(|| format!("failed to create {}", path.display()))?;
    let writer = BufWriter::new(file);
    serde_json::to_writer_pretty(writer, value)
        .with_context(|| format!("failed to write {}", path.display()))?;
    Ok(())
}

fn write_rmsd_csv(path: &Path, samples: &[RmsdSample]) -> Result<()> {
    let file =
        File::create(path).with_context(|| format!("failed to create {}", path.display()))?;
    let mut writer = BufWriter::new(file);
    writeln!(writer, "frame_index,step,rmsd")?;
    for sample in samples {
        writeln!(
            writer,
            "{},{},{:.10}",
            sample.frame_index, sample.step, sample.rmsd
        )?;
    }
    Ok(())
}

fn write_pair_distance_csv(path: &Path, samples: &[PairDistanceSample]) -> Result<()> {
    let file =
        File::create(path).with_context(|| format!("failed to create {}", path.display()))?;
    let mut writer = BufWriter::new(file);
    writeln!(writer, "frame_index,step,atom_i,atom_j,distance")?;
    for sample in samples {
        writeln!(
            writer,
            "{},{},{},{},{:.10}",
            sample.frame_index, sample.step, sample.atom_i, sample.atom_j, sample.distance
        )?;
    }
    Ok(())
}

fn write_analysis_report(path: &Path, summary: &md_analysis::AnalysisSummary) -> Result<()> {
    let file =
        File::create(path).with_context(|| format!("failed to create {}", path.display()))?;
    let mut writer = BufWriter::new(file);

    writeln!(writer, "# MD Analysis Report")?;
    writeln!(writer)?;
    writeln!(writer, "## System")?;
    writeln!(writer)?;
    writeln!(writer, "- Frames: {}", summary.frame_count)?;
    writeln!(writer, "- Atoms: {}", summary.atom_count)?;
    writeln!(writer)?;
    writeln!(writer, "## Methods")?;
    writeln!(writer)?;
    writeln!(writer, "- RMSD alignment: {}", summary.rmsd_alignment)?;
    writeln!(writer, "- Periodic unwrap: {}", summary.periodic_unwrap)?;
    if let Some([x, y, z]) = summary.periodic_box {
        writeln!(writer, "- Periodic box: {:.10} x {:.10} x {:.10}", x, y, z)?;
    }
    writeln!(
        writer,
        "- Assumptions: stable atom ordering, unchanged atom count, and one continuous trajectory image seeded from the first frame"
    )?;
    writeln!(writer)?;
    writeln!(writer, "## Energy")?;
    writeln!(writer)?;
    writeln!(writer, "- Samples: {}", summary.energy.sample_count)?;
    writeln!(
        writer,
        "- Initial total: {:.10}",
        summary.energy.initial_total
    )?;
    writeln!(writer, "- Final total: {:.10}", summary.energy.final_total)?;
    writeln!(writer, "- Total drift: {:.10}", summary.energy.total_drift)?;
    if let Some(drift_per_time) = summary.energy.total_drift_per_time {
        writeln!(writer, "- Drift per time: {:.10}", drift_per_time)?;
    }
    writeln!(writer, "- Mean total: {:.10}", summary.energy.mean_total)?;
    writeln!(
        writer,
        "- Temperature min/mean/max: {:.10} / {:.10} / {:.10}",
        summary.energy.min_temperature,
        summary.energy.mean_temperature,
        summary.energy.max_temperature
    )?;
    writeln!(writer)?;
    writeln!(writer, "## RMSD")?;
    writeln!(writer)?;
    writeln!(
        writer,
        "- Reference frame: {}",
        summary.rmsd.reference_frame
    )?;
    writeln!(writer, "- Final RMSD: {:.10}", summary.rmsd.final_value)?;
    writeln!(
        writer,
        "- RMSD min/mean/max: {:.10} / {:.10} / {:.10}",
        summary.rmsd.min, summary.rmsd.mean, summary.rmsd.max
    )?;

    if let Some(pair_distance) = &summary.pair_distance {
        writeln!(writer)?;
        writeln!(writer, "## Pair Distance")?;
        writeln!(writer)?;
        writeln!(
            writer,
            "- Atoms: {} and {}",
            pair_distance.atom_i, pair_distance.atom_j
        )?;
        writeln!(
            writer,
            "- Final distance: {:.10}",
            pair_distance.final_value
        )?;
        writeln!(
            writer,
            "- Distance min/mean/max: {:.10} / {:.10} / {:.10}",
            pair_distance.min, pair_distance.mean, pair_distance.max
        )?;
    }

    Ok(())
}

fn update_manifest_after_analysis(run_dir: &Path) -> Result<()> {
    let manifest_path = run_dir.join("run-manifest.json");
    let summary_path = run_dir.join("summary.json");

    let mut manifest = if manifest_path.exists() {
        read_manifest_json(&manifest_path)
            .with_context(|| format!("failed to read {}", manifest_path.display()))?
    } else if summary_path.exists() {
        let summary = read_summary_json(&summary_path)
            .with_context(|| format!("failed to read {}", summary_path.display()))?;
        manifest_from_summary(&summary)
    } else {
        return Ok(());
    };

    manifest.analysis_summary_file = Some("analysis-summary.json".to_string());
    push_unique(&mut manifest.outputs, "analysis-summary.json");
    push_unique(&mut manifest.outputs, "analysis-rmsd.csv");
    if run_dir.join("analysis-distance.csv").exists() {
        push_unique(&mut manifest.outputs, "analysis-distance.csv");
    }
    push_unique(&mut manifest.outputs, "analysis-report.md");
    write_manifest_json(&manifest_path, &manifest)
        .with_context(|| format!("failed to write {}", manifest_path.display()))?;

    Ok(())
}

fn read_optional_analysis_summary(run_dir: &Path) -> Result<Option<AnalysisSummary>> {
    let path = run_dir.join("analysis-summary.json");
    if !path.exists() {
        return Ok(None);
    }
    let file = File::open(&path).with_context(|| format!("failed to read {}", path.display()))?;
    let reader = std::io::BufReader::new(file);
    let summary = serde_json::from_reader(reader)
        .with_context(|| format!("failed to parse {}", path.display()))?;
    Ok(Some(summary))
}

fn manifest_from_summary(summary: &RunSummary) -> RunManifest {
    let mut outputs = summary.outputs.clone();
    push_unique(&mut outputs, "run-manifest.json");

    RunManifest {
        schema_version: 1,
        application: "md-workstation".to_string(),
        run_name: summary.run_name.clone(),
        project_name: summary.project_name.clone(),
        status: "completed".to_string(),
        created_at_unix_seconds: current_unix_seconds(),
        completed_at_unix_seconds: Some(current_unix_seconds()),
        config_file: "config.toml".to_string(),
        summary_file: "summary.json".to_string(),
        trajectory_file: "trajectory.xyz".to_string(),
        energy_file: "energy.csv".to_string(),
        input_file: summary.input_file.clone(),
        topology_file: summary.topology_file.clone(),
        checkpoint_file: summary
            .outputs
            .iter()
            .find(|output| {
                output.as_str() == "checkpoint.json" || output.as_str() == "checkpoint.bin"
            })
            .cloned(),
        checkpoint_format: summary
            .outputs
            .iter()
            .find(|output| {
                output.as_str() == "checkpoint.json" || output.as_str() == "checkpoint.bin"
            })
            .map(|output| checkpoint_format_from_file(output).to_string()),
        analysis_summary_file: None,
        report_file: None,
        outputs,
    }
}

fn push_unique(outputs: &mut Vec<String>, value: &str) {
    if !outputs.iter().any(|output| output == value) {
        outputs.push(value.to_string());
    }
}

fn checkpoint_info_for_resume(run_dir: &Path, config: &RunConfig) -> Result<(String, String)> {
    let manifest_path = run_dir.join("run-manifest.json");
    if manifest_path.exists() {
        let manifest = read_manifest_json(&manifest_path)
            .with_context(|| format!("failed to read {}", manifest_path.display()))?;
        if let Some(checkpoint_file) = manifest.checkpoint_file {
            let checkpoint_format = manifest
                .checkpoint_format
                .unwrap_or_else(|| checkpoint_format_from_file(&checkpoint_file).to_string());
            return Ok((checkpoint_file, checkpoint_format));
        }
    }

    let configured_file = config.checkpoint.file_name().to_string();
    if run_dir.join(&configured_file).exists() {
        return Ok((configured_file, config.checkpoint.normalized_format()));
    }
    if run_dir.join("checkpoint.json").exists() {
        return Ok(("checkpoint.json".to_string(), "json".to_string()));
    }
    if run_dir.join("checkpoint.bin").exists() {
        return Ok(("checkpoint.bin".to_string(), "binary".to_string()));
    }

    bail!(
        "no checkpoint file found in {}; expected checkpoint.json or checkpoint.bin",
        run_dir.display()
    )
}

fn write_run_checkpoint(
    run_dir: &Path,
    run_name: &str,
    step: usize,
    time: f64,
    potential_energy: f64,
    state: &SystemState,
    checkpoint: &CheckpointSection,
) -> Result<()> {
    let checkpoint_data = RunCheckpoint {
        schema_version: 1,
        run_name: run_name.to_string(),
        step,
        time,
        potential_energy,
        state: state.clone(),
    };
    let checkpoint_path = run_dir.join(checkpoint.file_name());
    write_checkpoint_by_format(
        &checkpoint_path,
        &checkpoint_data,
        &checkpoint.normalized_format(),
    )
    .with_context(|| format!("failed to write {}", checkpoint_path.display()))
}

fn write_checkpoint_by_format(path: &Path, checkpoint: &RunCheckpoint, format: &str) -> Result<()> {
    match format {
        "binary" => Ok(write_checkpoint_binary(path, checkpoint)?),
        "json" => Ok(write_checkpoint_json(path, checkpoint)?),
        other => bail!("unsupported checkpoint format {other:?}"),
    }
}

fn read_checkpoint_by_format(path: &Path, format: &str) -> Result<RunCheckpoint> {
    match format {
        "binary" => Ok(read_checkpoint_binary(path)?),
        "json" => Ok(read_checkpoint_json(path)?),
        other => bail!("unsupported checkpoint format {other:?}"),
    }
}

fn checkpoint_format_from_file(file_name: &str) -> &'static str {
    if file_name.ends_with(".bin") {
        "binary"
    } else {
        "json"
    }
}

fn existing_checkpoint_file(run_dir: &Path) -> Option<&'static str> {
    if run_dir.join("checkpoint.bin").exists() {
        Some("checkpoint.bin")
    } else if run_dir.join("checkpoint.json").exists() {
        Some("checkpoint.json")
    } else {
        None
    }
}

fn load_topology_for_resume(
    summary: &RunSummary,
    run_dir: &Path,
    particle_count: usize,
) -> Result<Option<LoadedTopology>> {
    let Some(topology_file) = &summary.topology_file else {
        return Ok(None);
    };
    load_topology_from_path(run_dir.join(topology_file), particle_count).map(Some)
}

fn should_output(step: usize, steps: usize, output_interval: usize) -> bool {
    step % output_interval == 0 || step == steps
}

fn parse_pair_indices(pair: &str) -> Result<(usize, usize)> {
    let (left, right) = pair
        .split_once(',')
        .ok_or_else(|| anyhow::anyhow!("pair must be formatted as atom_i,atom_j"))?;
    let atom_i = left
        .trim()
        .parse()
        .with_context(|| format!("invalid atom index {:?}", left.trim()))?;
    let atom_j = right
        .trim()
        .parse()
        .with_context(|| format!("invalid atom index {:?}", right.trim()))?;
    if atom_i == atom_j {
        bail!("pair atom indices must be different");
    }
    Ok((atom_i, atom_j))
}

fn next_run_dir(base: &Path, name: &str) -> Result<PathBuf> {
    for index in 1..=999 {
        let candidate = base.join(format!("{name}-{index:03}"));
        if !candidate.exists() {
            return Ok(candidate);
        }
    }

    bail!(
        "could not find an available run directory for {name} under {}",
        base.display()
    )
}

fn lattice_grid_size(particle_count: usize) -> usize {
    let mut grid = 1;
    while grid * grid * grid < particle_count {
        grid += 1;
    }
    grid
}

fn config_stem(config_path: &Path) -> String {
    config_path
        .file_stem()
        .and_then(|stem| stem.to_str())
        .unwrap_or("run")
        .to_string()
}

fn sanitize_run_name(name: &str) -> String {
    let sanitized: String = name
        .chars()
        .map(|character| {
            if character.is_ascii_alphanumeric() || character == '-' || character == '_' {
                character
            } else {
                '-'
            }
        })
        .collect();

    let trimmed = sanitized.trim_matches('-');
    if trimmed.is_empty() {
        "run".to_string()
    } else {
        trimmed.to_string()
    }
}

fn resolve_config_relative_path(path: &Path, config_dir: Option<&Path>) -> PathBuf {
    if path.is_absolute() {
        path.to_path_buf()
    } else {
        config_dir.unwrap_or_else(|| Path::new(".")).join(path)
    }
}

fn current_unix_seconds() -> u64 {
    SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap_or_default()
        .as_secs()
}

fn max_force_delta(left: &SystemState, right: &SystemState) -> f64 {
    let mut max_delta: f64 = 0.0;
    for i in 0..left.particle_count().min(right.particle_count()) {
        max_delta = max_delta.max((left.fx[i] - right.fx[i]).abs());
        max_delta = max_delta.max((left.fy[i] - right.fy[i]).abs());
        max_delta = max_delta.max((left.fz[i] - right.fz[i]).abs());
    }
    max_delta
}

fn speedup_ratio(baseline: Duration, candidate: Duration) -> f64 {
    let candidate_seconds = candidate.as_secs_f64();
    if candidate_seconds == 0.0 {
        f64::INFINITY
    } else {
        baseline.as_secs_f64() / candidate_seconds
    }
}

fn default_particles() -> usize {
    32
}

fn default_mass() -> f64 {
    1.0
}

fn default_element() -> String {
    "Ar".to_string()
}

fn default_seed() -> u64 {
    42
}

fn default_input_format() -> String {
    "xyz".to_string()
}

fn default_coulomb_constant() -> f64 {
    1.0
}

fn default_neighbor_skin() -> f64 {
    0.3
}

fn default_neighbor_rebuild_interval() -> usize {
    10
}

fn default_thermostat_type() -> String {
    "none".to_string()
}

fn default_thermostat_tau() -> f64 {
    0.1
}

fn default_minimization_steps() -> usize {
    200
}

fn default_minimization_output_interval() -> usize {
    10
}

fn default_minimization_step_size() -> f64 {
    0.001
}

fn default_minimization_max_displacement() -> f64 {
    0.02
}

fn default_minimization_force_tolerance() -> f64 {
    1.0e-4
}

fn default_minimization_max_backtracks() -> usize {
    12
}

fn default_checkpoint_format() -> String {
    "json".to_string()
}

fn default_mixing_rule() -> String {
    "lorentz-berthelot".to_string()
}

fn default_output_directory() -> PathBuf {
    PathBuf::from("runs")
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::time::{SystemTime, UNIX_EPOCH};

    #[test]
    fn lattice_grid_fits_particle_count() {
        assert_eq!(lattice_grid_size(1), 1);
        assert_eq!(lattice_grid_size(8), 2);
        assert_eq!(lattice_grid_size(9), 3);
        assert_eq!(lattice_grid_size(32), 4);
    }

    #[test]
    fn sanitize_run_name_keeps_cross_platform_safe_subset() {
        assert_eq!(sanitize_run_name("lj fluid/01"), "lj-fluid-01");
        assert_eq!(sanitize_run_name("..."), "run");
    }

    #[test]
    fn parses_pair_indices() {
        assert_eq!(parse_pair_indices("0,1").unwrap(), (0, 1));
        assert_eq!(parse_pair_indices(" 2, 5 ").unwrap(), (2, 5));
        assert!(parse_pair_indices("1,1").is_err());
        assert!(parse_pair_indices("1:2").is_err());
    }

    #[test]
    fn parses_and_validates_minimal_lj_config() {
        let config: RunConfig = toml::from_str(
            r#"
            [simulation]
            steps = 10
            dt = 0.001
            output_interval = 5
            temperature = 0.2

            [system]
            particles = 8
            mass = 1.0
            element = "Ar"
            seed = 7
            lattice_spacing = 1.4

            [box]
            x = 6.0
            y = 6.0
            z = 6.0

            [force]
            type = "lennard-jones"
            epsilon = 1.0
            sigma = 1.0
            cutoff = 2.5
            "#,
        )
        .unwrap();

        config.validate().unwrap();
    }

    #[test]
    fn rejects_invalid_force_type() {
        let config: RunConfig = toml::from_str(
            r#"
            [simulation]
            steps = 10
            dt = 0.001
            output_interval = 5
            temperature = 0.2

            [system]
            particles = 8

            [box]
            x = 6.0
            y = 6.0
            z = 6.0

            [force]
            type = "unknown"
            epsilon = 1.0
            sigma = 1.0
            cutoff = 2.5
            "#,
        )
        .unwrap();

        let error = config.validate().unwrap_err().to_string();

        assert!(error.contains("unsupported force.type"));
    }

    #[test]
    fn rejects_lattice_that_does_not_fit_box() {
        let config: RunConfig = toml::from_str(
            r#"
            [simulation]
            steps = 10
            dt = 0.001
            output_interval = 5
            temperature = 0.2

            [system]
            particles = 27
            lattice_spacing = 2.0

            [box]
            x = 3.0
            y = 3.0
            z = 3.0

            [force]
            type = "lennard-jones"
            epsilon = 1.0
            sigma = 1.0
            cutoff = 2.5
            "#,
        )
        .unwrap();

        let error = config.validate().unwrap_err().to_string();

        assert!(error.contains("does not fit in box"));
    }

    #[test]
    fn rejects_periodic_cutoff_larger_than_half_shortest_box_length() {
        let config: RunConfig = toml::from_str(
            r#"
            [simulation]
            steps = 10
            dt = 0.001
            output_interval = 5
            temperature = 0.2

            [system]
            particles = 8
            lattice_spacing = 1.0

            [box]
            x = 4.0
            y = 4.0
            z = 4.0
            periodic = true

            [force]
            type = "lennard-jones"
            epsilon = 1.0
            sigma = 1.0
            cutoff = 2.1
            "#,
        )
        .unwrap();

        let error = config.validate().unwrap_err().to_string();

        assert!(error.contains("periodic cutoff"));
    }

    #[test]
    fn periodic_config_builds_periodic_shifted_force_options() {
        let config: RunConfig = toml::from_str(
            r#"
            [simulation]
            steps = 10
            dt = 0.001
            output_interval = 5
            temperature = 0.2

            [system]
            particles = 8
            lattice_spacing = 1.0

            [box]
            x = 6.0
            y = 6.0
            z = 6.0
            periodic = true

            [force]
            type = "lennard-jones"
            epsilon = 1.0
            sigma = 1.0
            cutoff = 2.5
            shift_potential = true
            "#,
        )
        .unwrap();

        let options = config.force_options().unwrap();

        assert!(matches!(options.boundary, BoundaryCondition::Periodic(_)));
        assert!(options.shift_potential);
        assert_eq!(config.boundary_label(), "periodic");
    }

    #[test]
    fn parallel_execution_config_is_parsed() {
        let config: RunConfig = toml::from_str(
            r#"
            [simulation]
            steps = 10
            dt = 0.001
            output_interval = 5
            temperature = 0.2

            [system]
            particles = 8
            lattice_spacing = 1.0

            [box]
            x = 6.0
            y = 6.0
            z = 6.0
            periodic = true

            [force]
            type = "lennard-jones"
            epsilon = 1.0
            sigma = 1.0
            cutoff = 2.5
            shift_potential = true

            [execution]
            parallel = true
            "#,
        )
        .unwrap();

        config.validate().unwrap();
        assert!(config.execution.parallel);
    }

    #[test]
    fn rejects_invalid_neighbor_settings() {
        let config: RunConfig = toml::from_str(
            r#"
            [simulation]
            steps = 10
            dt = 0.001
            output_interval = 5
            temperature = 0.2

            [system]
            particles = 8
            lattice_spacing = 1.0

            [box]
            x = 6.0
            y = 6.0
            z = 6.0
            periodic = true

            [force]
            type = "lennard-jones"
            epsilon = 1.0
            sigma = 1.0
            cutoff = 2.5
            shift_potential = true

            [neighbor]
            enabled = true
            skin = -0.1
            rebuild_interval = 10
            "#,
        )
        .unwrap();

        let error = config.validate().unwrap_err().to_string();

        assert!(error.contains("neighbor.skin"));
    }

    #[test]
    fn rejects_invalid_thermostat_settings() {
        let config: RunConfig = toml::from_str(
            r#"
            [simulation]
            steps = 10
            dt = 0.001
            output_interval = 5
            temperature = 0.2

            [system]
            particles = 8
            lattice_spacing = 1.0

            [box]
            x = 6.0
            y = 6.0
            z = 6.0

            [force]
            type = "lennard-jones"
            epsilon = 1.0
            sigma = 1.0
            cutoff = 2.5

            [thermostat]
            type = "berendsen"
            tau = 0.0005
            "#,
        )
        .unwrap();

        let error = config.validate().unwrap_err().to_string();

        assert!(error.contains("thermostat.tau"));
    }

    #[test]
    fn berendsen_thermostat_moves_temperature_toward_target() {
        let mut state = SystemState::new(2);
        state.vx[0] = 1.0;
        state.vx[1] = -1.0;
        let initial_temperature = state.temperature();
        let target_temperature = 0.2;
        let dt = 0.001;
        let tau = 0.01;

        apply_thermostat(
            &mut state,
            ThermostatOptions::Berendsen {
                target_temperature,
                tau,
            },
            dt,
        )
        .unwrap();

        let expected_temperature =
            initial_temperature + (dt / tau) * (target_temperature - initial_temperature);
        assert!((state.temperature() - expected_temperature).abs() < 1.0e-12);
    }

    #[test]
    fn deterministic_seed_recreates_initial_velocities() {
        let config = minimal_test_config(123, std::env::temp_dir());

        let first = build_initial_system(&config, None).unwrap();
        let second = build_initial_system(&config, None).unwrap();

        assert_eq!(first.vx, second.vx);
        assert_eq!(first.vy, second.vy);
        assert_eq!(first.vz, second.vz);
        assert!((first.temperature() - config.simulation.temperature).abs() < 1.0e-12);
    }

    #[test]
    fn different_seeds_change_initial_velocities() {
        let first =
            build_initial_system(&minimal_test_config(123, std::env::temp_dir()), None).unwrap();
        let second =
            build_initial_system(&minimal_test_config(456, std::env::temp_dir()), None).unwrap();

        assert_ne!(first.vx, second.vx);
        assert_ne!(first.vy, second.vy);
        assert_ne!(first.vz, second.vz);
    }

    #[test]
    fn xyz_input_builds_system_from_config_relative_path() {
        let root = unique_temp_dir("md-cli-xyz-build");
        let config_dir = root.join("configs");
        let molecule_dir = root.join("molecules");
        std::fs::create_dir_all(&config_dir).unwrap();
        std::fs::create_dir_all(&molecule_dir).unwrap();
        std::fs::write(
            molecule_dir.join("dimer.xyz"),
            "2\nplain molecule comment\nAr 1.0 2.0 3.0\nNe 2.5 2.0 3.0\n",
        )
        .unwrap();

        let config: RunConfig = toml::from_str(
            r#"
            [simulation]
            steps = 2
            dt = 0.001
            output_interval = 1
            temperature = 0.0

            [input]
            path = "../molecules/dimer.xyz"
            format = "xyz"
            frame = 0

            [system]
            mass = 2.0
            seed = 123

            [box]
            x = 4.0
            y = 4.0
            z = 4.0

            [force]
            type = "lennard-jones"
            epsilon = 1.0
            sigma = 1.0
            cutoff = 1.8
            "#,
        )
        .unwrap();

        config.validate().unwrap();
        let state = build_initial_system(&config, Some(&config_dir)).unwrap();

        assert_eq!(state.particle_count(), 2);
        assert_eq!(state.element, vec!["Ar".to_string(), "Ne".to_string()]);
        assert_eq!(state.mass, vec![2.0, 2.0]);
        assert_eq!(state.x, vec![1.0, 2.5]);
        assert_eq!(state.vx, vec![0.0, 0.0]);

        std::fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn force_field_types_apply_defaults_and_mixed_lj_parameters() {
        let root = unique_temp_dir("md-cli-force-field-types");
        let config_dir = root.join("configs");
        let molecule_dir = root.join("molecules");
        std::fs::create_dir_all(&config_dir).unwrap();
        std::fs::create_dir_all(&molecule_dir).unwrap();
        let mixed_sigma = 1.5_f64;
        let equilibrium_distance = mixed_sigma * 2.0_f64.powf(1.0 / 6.0);
        std::fs::write(
            molecule_dir.join("argon-neon.xyz"),
            format!(
                "2\nmixed LJ pair\nAr 1.0 2.0 3.0\nNe {} 2.0 3.0\n",
                1.0 + equilibrium_distance
            ),
        )
        .unwrap();

        let config: RunConfig = toml::from_str(
            r#"
            [simulation]
            steps = 1
            dt = 0.001
            output_interval = 1
            temperature = 0.0

            [input]
            path = "../molecules/argon-neon.xyz"
            format = "xyz"
            frame = 0

            [system]
            mass = 1.0
            seed = 123

            [box]
            x = 6.0
            y = 6.0
            z = 6.0

            [force]
            type = "lennard-jones"
            epsilon = 1.0
            sigma = 1.0
            cutoff = 4.0

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
            "#,
        )
        .unwrap();

        config.validate().unwrap();
        let (mut state, topology) =
            prepare_initial_state_and_topology(&config, Some(&config_dir)).unwrap();
        apply_topology_to_state(&mut state, topology.as_ref()).unwrap();
        let force_options = config.force_options().unwrap();
        let mixed_force_options = config
            .mixed_lennard_jones_options(&state, &force_options)
            .unwrap()
            .unwrap();

        let potential = compute_configured_forces(
            &mut state,
            &force_options,
            Some(&mixed_force_options),
            None,
            topology.as_ref(),
            None,
            false,
        )
        .unwrap();

        assert_eq!(state.mass, vec![2.0, 4.0]);
        assert_eq!(state.charge, vec![0.5, -0.25]);
        assert!((potential - -2.0).abs() < 1.0e-12);

        std::fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn run_with_xyz_input_copies_input_and_writes_metadata() {
        let root = unique_temp_dir("md-cli-xyz-run");
        let config_dir = root.join("configs");
        let molecule_dir = root.join("molecules");
        let output_dir = root.join("runs");
        std::fs::create_dir_all(&config_dir).unwrap();
        std::fs::create_dir_all(&molecule_dir).unwrap();
        std::fs::write(
            molecule_dir.join("dimer.xyz"),
            "2\nplain molecule comment\nAr 1.0 2.0 3.0\nAr 2.5 2.0 3.0\n",
        )
        .unwrap();
        let config_path = config_dir.join("xyz-input.toml");
        std::fs::write(
            &config_path,
            format!(
                r#"
                [simulation]
                steps = 2
                dt = 0.001
                output_interval = 1
                temperature = 0.0

                [input]
                path = "../molecules/dimer.xyz"
                format = "xyz"
                frame = 0

                [system]
                mass = 1.0
                seed = 123

                [box]
                x = 4.0
                y = 4.0
                z = 4.0

                [force]
                type = "lennard-jones"
                epsilon = 1.0
                sigma = 1.0
                cutoff = 1.8

                [output]
                directory = "{}"
                name = "xyz-input-test"
                "#,
                toml_string(&output_dir)
            ),
        )
        .unwrap();

        run_command(&config_path).unwrap();

        let run_dir = output_dir.join("xyz-input-test-001");
        let summary = std::fs::read_to_string(run_dir.join("summary.json")).unwrap();
        let manifest = std::fs::read_to_string(run_dir.join("run-manifest.json")).unwrap();
        let report = std::fs::read_to_string(run_dir.join("run-report.md")).unwrap();
        let copied_input = std::fs::read_to_string(run_dir.join("input.xyz")).unwrap();
        let trajectory = std::fs::read_to_string(run_dir.join("trajectory.xyz")).unwrap();

        assert!(summary.contains(r#""input_file": "input.xyz""#));
        assert!(summary.contains(r#""input_format": "xyz""#));
        assert!(summary.contains(r#""input_frame": 0"#));
        assert!(manifest.contains(r#""report_file": "run-report.md""#));
        assert!(report.contains("# MD Run Report"));
        assert!(copied_input.contains("plain molecule comment"));
        assert!(trajectory.starts_with("2\nstep=0\nAr"));

        std::fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn run_with_pdb_input_preserves_atom_metadata() {
        let root = unique_temp_dir("md-cli-pdb-run");
        let config_dir = root.join("configs");
        let molecule_dir = root.join("molecules");
        let output_dir = root.join("runs");
        std::fs::create_dir_all(&config_dir).unwrap();
        std::fs::create_dir_all(&molecule_dir).unwrap();
        std::fs::write(
            molecule_dir.join("alanine-fragment.pdb"),
            "ATOM      1  N   ALA A   7       1.000   2.000   3.000  1.00 10.00           N  \nATOM      2  CA  ALA A   7       2.000   2.000   3.000  1.00 10.00           C  \nATOM      3  C   ALA A   7       2.500   3.200   3.100  1.00 10.00           C  \nEND\n",
        )
        .unwrap();
        let config_path = config_dir.join("pdb-input.toml");
        std::fs::write(
            &config_path,
            format!(
                r#"
                [simulation]
                steps = 2
                dt = 0.001
                output_interval = 1
                temperature = 0.0

                [input]
                path = "../molecules/alanine-fragment.pdb"
                format = "pdb"
                frame = 0

                [system]
                mass = 1.0
                seed = 123

                [box]
                x = 6.0
                y = 6.0
                z = 6.0

                [force]
                type = "lennard-jones"
                epsilon = 0.1
                sigma = 1.0
                cutoff = 2.5

                [output]
                directory = "{}"
                name = "pdb-input-test"
                "#,
                toml_string(&output_dir)
            ),
        )
        .unwrap();

        run_command(&config_path).unwrap();

        let run_dir = output_dir.join("pdb-input-test-001");
        let summary = std::fs::read_to_string(run_dir.join("summary.json")).unwrap();
        let manifest = std::fs::read_to_string(run_dir.join("run-manifest.json")).unwrap();
        let report = std::fs::read_to_string(run_dir.join("run-report.md")).unwrap();
        let metadata = std::fs::read_to_string(run_dir.join("atom-metadata.csv")).unwrap();
        let copied_input = std::fs::read_to_string(run_dir.join("input.pdb")).unwrap();

        assert!(summary.contains(r#""input_file": "input.pdb""#));
        assert!(summary.contains(r#""input_format": "pdb""#));
        assert!(summary.contains(r#""atom_metadata_file": "atom-metadata.csv""#));
        assert!(manifest.contains(r#""input.pdb""#));
        assert!(manifest.contains(r#""atom-metadata.csv""#));
        assert!(report.contains("- Atom metadata: atom-metadata.csv"));
        assert!(metadata.contains("0,N,N,ALA,7,A"));
        assert!(metadata.contains("1,C,CA,ALA,7,A"));
        assert!(copied_input.contains("ALA A   7"));

        std::fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn minimize_command_lowers_potential_and_writes_history() {
        let root = unique_temp_dir("md-cli-minimize");
        let config_dir = root.join("configs");
        let molecule_dir = root.join("molecules");
        let output_dir = root.join("runs");
        std::fs::create_dir_all(&config_dir).unwrap();
        std::fs::create_dir_all(&molecule_dir).unwrap();
        std::fs::write(
            molecule_dir.join("close-dimer.xyz"),
            "2\nclose LJ dimer\nAr 2.0 2.0 2.0\nAr 2.9 2.0 2.0\n",
        )
        .unwrap();
        let config_path = config_dir.join("minimize.toml");
        std::fs::write(
            &config_path,
            format!(
                r#"
                [simulation]
                steps = 1
                dt = 0.001
                output_interval = 1
                temperature = 0.0

                [input]
                path = "../molecules/close-dimer.xyz"
                format = "xyz"
                frame = 0

                [system]
                mass = 1.0
                seed = 123

                [box]
                x = 6.0
                y = 6.0
                z = 6.0

                [force]
                type = "lennard-jones"
                epsilon = 1.0
                sigma = 1.0
                cutoff = 2.5

                [minimization]
                steps = 40
                output_interval = 10
                step_size = 0.001
                max_displacement = 0.02
                force_tolerance = 0.000001

                [output]
                directory = "{}"
                name = "minimize-test"
                "#,
                toml_string(&output_dir)
            ),
        )
        .unwrap();

        minimize_command(&config_path).unwrap();

        let run_dir = output_dir.join("minimize-test-001");
        let history = std::fs::read_to_string(run_dir.join("minimization.csv")).unwrap();
        let report = std::fs::read_to_string(run_dir.join("run-report.md")).unwrap();
        let summary = std::fs::read_to_string(run_dir.join("summary.json")).unwrap();
        let potentials: Vec<f64> = history
            .lines()
            .skip(1)
            .map(|line| line.split(',').nth(1).unwrap().parse::<f64>().unwrap())
            .collect();

        assert!(potentials.last().unwrap() < potentials.first().unwrap());
        assert!(run_dir.join("minimized.xyz").exists());
        assert!(report.contains("- Workflow: minimization"));
        assert!(report.contains("- Minimization final max force:"));
        assert!(summary.contains(r#""workflow": "minimization""#));

        std::fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn run_with_topology_copies_topology_and_writes_force_field_metadata() {
        let root = unique_temp_dir("md-cli-topology-run");
        let config_dir = root.join("configs");
        let molecule_dir = root.join("molecules");
        let topology_dir = root.join("topology");
        let output_dir = root.join("runs");
        std::fs::create_dir_all(&config_dir).unwrap();
        std::fs::create_dir_all(&molecule_dir).unwrap();
        std::fs::create_dir_all(&topology_dir).unwrap();
        std::fs::write(
            molecule_dir.join("force-field.xyz"),
            "4\nforce-field report\nAr 1.8 3.0 3.0\nAr 3.0 3.0 3.0\nAr 3.0 4.2 3.0\nAr 4.0 4.6 3.8\n",
        )
        .unwrap();
        std::fs::write(
            topology_dir.join("force-field.toml"),
            r#"
            charges = [0.5, 0.0, -0.5, 0.25]

            [[bonds]]
            i = 0
            j = 1
            k = 25.0
            r0 = 1.2

            [[bonds]]
            i = 1
            j = 2
            k = 25.0
            r0 = 1.2

            [[bonds]]
            i = 2
            j = 3
            k = 20.0
            r0 = 1.3416407865

            [[angles]]
            i = 0
            j = 1
            k = 2
            force_constant = 10.0
            theta0 = 1.5707963267948966

            [[angles]]
            i = 1
            j = 2
            k = 3
            force_constant = 8.0
            theta0 = 1.8736811952

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
            "#,
        )
        .unwrap();
        let config_path = config_dir.join("force-field-report.toml");
        std::fs::write(
            &config_path,
            format!(
                r#"
                [simulation]
                steps = 2
                dt = 0.001
                output_interval = 1
                temperature = 0.0

                [input]
                path = "../molecules/force-field.xyz"
                format = "xyz"
                frame = 0

                [topology]
                path = "../topology/force-field.toml"

                [system]
                mass = 1.0
                seed = 123

                [box]
                x = 6.0
                y = 6.0
                z = 6.0

                [force]
                type = "lennard-jones"
                epsilon = 0.02
                sigma = 1.0
                cutoff = 2.5

                [coulomb]
                enabled = true
                constant = 0.1
                cutoff = 3.0

                [output]
                directory = "{}"
                name = "force-field-report-test"
                "#,
                toml_string(&output_dir)
            ),
        )
        .unwrap();

        run_command(&config_path).unwrap();

        let run_dir = output_dir.join("force-field-report-test-001");
        let summary = std::fs::read_to_string(run_dir.join("summary.json")).unwrap();
        let manifest = std::fs::read_to_string(run_dir.join("run-manifest.json")).unwrap();
        let report = std::fs::read_to_string(run_dir.join("run-report.md")).unwrap();
        let topology = std::fs::read_to_string(run_dir.join("topology.toml")).unwrap();

        assert!(summary.contains(r#""topology_file": "topology.toml""#));
        assert!(summary.contains(r#""bond_count": 3"#));
        assert!(summary.contains(r#""angle_count": 2"#));
        assert!(summary.contains(r#""dihedral_count": 1"#));
        assert!(summary.contains(r#""excluded_pair_count": 4"#));
        assert!(summary.contains(r#""coulomb": true"#));
        assert!(manifest.contains(r#""topology_file": "topology.toml""#));
        assert!(report.contains("- Harmonic bonds: 3"));
        assert!(report.contains("- Harmonic angles: 2"));
        assert!(report.contains("- Periodic dihedrals: 1"));
        assert!(report.contains("- Non-bonded exclusions: 4"));
        assert!(report.contains("- Coulomb enabled: true"));
        assert!(topology.contains("[[bonds]]"));
        assert!(topology.contains("[[angles]]"));
        assert!(topology.contains("[[dihedrals]]"));
        assert!(topology.contains("[[exclusions]]"));

        std::fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn benchmark_neighbor_accepts_topology_coulomb_and_exclusions() {
        let root = unique_temp_dir("md-cli-topology-benchmark");
        let config_dir = root.join("configs");
        let molecule_dir = root.join("molecules");
        let topology_dir = root.join("topology");
        std::fs::create_dir_all(&config_dir).unwrap();
        std::fs::create_dir_all(&molecule_dir).unwrap();
        std::fs::create_dir_all(&topology_dir).unwrap();
        std::fs::write(
            molecule_dir.join("trimer.xyz"),
            "3\ntopology benchmark\nAr 1.8 3.0 3.0\nAr 3.0 3.0 3.0\nAr 3.0 4.2 3.0\n",
        )
        .unwrap();
        std::fs::write(
            topology_dir.join("trimer.toml"),
            r#"
            charges = [0.5, 0.0, -0.5]

            [[bonds]]
            i = 0
            j = 1
            k = 25.0
            r0 = 1.2

            [[bonds]]
            i = 1
            j = 2
            k = 25.0
            r0 = 1.2

            [[angles]]
            i = 0
            j = 1
            k = 2
            force_constant = 10.0
            theta0 = 1.5707963267948966
            "#,
        )
        .unwrap();
        let config_path = config_dir.join("topology-neighbor.toml");
        std::fs::write(
            &config_path,
            r#"
            [simulation]
            steps = 2
            dt = 0.001
            output_interval = 1
            temperature = 0.0

            [input]
            path = "../molecules/trimer.xyz"
            format = "xyz"
            frame = 0

            [topology]
            path = "../topology/trimer.toml"

            [system]
            mass = 1.0
            seed = 123

            [box]
            x = 6.0
            y = 6.0
            z = 6.0

            [force]
            type = "lennard-jones"
            epsilon = 0.02
            sigma = 1.0
            cutoff = 2.5

            [coulomb]
            enabled = true
            constant = 0.1
            cutoff = 3.0

            [neighbor]
            enabled = true
            skin = 0.3
            rebuild_interval = 10
            "#,
        )
        .unwrap();

        benchmark_neighbor_command(&config_path, 2).unwrap();

        std::fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn repeated_cli_runs_with_same_config_write_identical_energy_logs() {
        let root = unique_temp_dir("md-cli-determinism");
        let output_dir = root.join("runs");
        std::fs::create_dir_all(&root).unwrap();
        let config_path = root.join("deterministic.toml");
        std::fs::write(&config_path, deterministic_config_toml(&output_dir)).unwrap();

        run_command(&config_path).unwrap();
        run_command(&config_path).unwrap();

        let first_energy =
            std::fs::read_to_string(output_dir.join("deterministic-001/energy.csv")).unwrap();
        let second_energy =
            std::fs::read_to_string(output_dir.join("deterministic-002/energy.csv")).unwrap();

        assert_eq!(first_energy, second_energy);
        assert!(first_energy.contains("20,0.0200000000"));

        std::fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn report_command_includes_analysis_outputs_when_available() {
        let root = unique_temp_dir("md-cli-report");
        let output_dir = root.join("runs");
        std::fs::create_dir_all(&root).unwrap();
        let config_path = root.join("report.toml");
        std::fs::write(&config_path, deterministic_config_toml(&output_dir)).unwrap();

        run_command(&config_path).unwrap();
        let run_dir = output_dir.join("deterministic-001");
        analyze_command(&run_dir, 0, Some("0,1")).unwrap();
        report_command(&run_dir).unwrap();

        let report = std::fs::read_to_string(run_dir.join("run-report.md")).unwrap();
        let manifest = std::fs::read_to_string(run_dir.join("run-manifest.json")).unwrap();

        assert!(report.contains("## Analysis"));
        assert!(report.contains("- Pair 0-1 final distance:"));
        assert!(manifest.contains(r#""analysis_summary_file": "analysis-summary.json""#));
        assert!(manifest.contains(r#""run-report.md""#));

        std::fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn analyze_command_unwraps_periodic_frames_from_run_config() {
        let root = unique_temp_dir("md-cli-periodic-analysis");
        let run_dir = root.join("runs/periodic-analysis-001");
        std::fs::create_dir_all(&run_dir).unwrap();
        std::fs::write(
            run_dir.join("config.toml"),
            r#"
            [simulation]
            steps = 1
            dt = 0.001
            output_interval = 1
            temperature = 0.0

            [system]
            particles = 3
            mass = 1.0
            element = "Ar"
            seed = 1
            lattice_spacing = 1.0

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
            "#,
        )
        .unwrap();
        std::fs::write(
            run_dir.join("energy.csv"),
            "step,time,kinetic,potential,total,temperature\n0,0.0000000000,0.0000000000,-1.0000000000,-1.0000000000,0.0000000000\n1,0.0010000000,0.0000000000,-1.0000000000,-1.0000000000,0.0000000000\n",
        )
        .unwrap();
        std::fs::write(
            run_dir.join("trajectory.xyz"),
            "3\nstep=0\nAr 9.5 1.0 1.0\nAr 9.5 2.0 1.0\nAr 8.5 1.0 1.0\n3\nstep=1\nAr 0.5 1.0 1.0\nAr 0.5 2.0 1.0\nAr 9.5 1.0 1.0\n",
        )
        .unwrap();

        analyze_command(&run_dir, 0, Some("0,2")).unwrap();

        let summary = std::fs::read_to_string(run_dir.join("analysis-summary.json")).unwrap();
        let report = std::fs::read_to_string(run_dir.join("analysis-report.md")).unwrap();
        let rmsd = std::fs::read_to_string(run_dir.join("analysis-rmsd.csv")).unwrap();
        let distance = std::fs::read_to_string(run_dir.join("analysis-distance.csv")).unwrap();

        assert!(summary.contains(r#""periodic_unwrap": true"#));
        assert!(summary.contains(r#""rmsd_alignment": "centered-kabsch""#));
        assert!(report.contains("- Periodic unwrap: true"));
        assert!(report.contains("- RMSD alignment: centered-kabsch"));
        assert!(rmsd.ends_with("1,1,0.0000000000\n"));
        assert!(distance.ends_with("1,1,0,2,1.0000000000\n"));

        std::fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn resume_command_extends_run_from_checkpoint() {
        let root = unique_temp_dir("md-cli-resume");
        let output_dir = root.join("runs");
        std::fs::create_dir_all(&root).unwrap();
        let config_path = root.join("resume.toml");
        std::fs::write(&config_path, deterministic_config_toml(&output_dir)).unwrap();

        run_command(&config_path).unwrap();
        let run_dir = output_dir.join("deterministic-001");
        resume_command(&run_dir, Some(5)).unwrap();

        let summary = std::fs::read_to_string(run_dir.join("summary.json")).unwrap();
        let checkpoint = std::fs::read_to_string(run_dir.join("checkpoint.json")).unwrap();
        let manifest = std::fs::read_to_string(run_dir.join("run-manifest.json")).unwrap();
        let energy = std::fs::read_to_string(run_dir.join("energy.csv")).unwrap();

        assert!(summary.contains(r#""steps": 25"#));
        assert!(checkpoint.contains(r#""step": 25"#));
        assert!(manifest.contains(r#""checkpoint_file": "checkpoint.json""#));
        assert!(manifest.contains(r#""checkpoint_format": "json""#));
        assert!(energy.contains("25,0.0250000000"));

        std::fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn binary_checkpoint_run_and_resume_use_manifest_format() {
        let root = unique_temp_dir("md-cli-binary-checkpoint");
        let output_dir = root.join("runs");
        std::fs::create_dir_all(&root).unwrap();
        let config_path = root.join("binary-checkpoint.toml");
        std::fs::write(&config_path, binary_checkpoint_config_toml(&output_dir)).unwrap();

        run_command(&config_path).unwrap();
        let run_dir = output_dir.join("binary-checkpoint-001");
        assert!(run_dir.join("checkpoint.bin").exists());
        assert!(!run_dir.join("checkpoint.json").exists());
        let checkpoint = read_checkpoint_binary(run_dir.join("checkpoint.bin")).unwrap();
        assert_eq!(checkpoint.step, 10);

        resume_command(&run_dir, Some(5)).unwrap();

        let checkpoint = read_checkpoint_binary(run_dir.join("checkpoint.bin")).unwrap();
        let manifest = std::fs::read_to_string(run_dir.join("run-manifest.json")).unwrap();
        let energy = std::fs::read_to_string(run_dir.join("energy.csv")).unwrap();

        assert_eq!(checkpoint.step, 15);
        assert!(manifest.contains(r#""checkpoint_file": "checkpoint.bin""#));
        assert!(manifest.contains(r#""checkpoint_format": "binary""#));
        assert!(energy.contains("15,0.0150000000"));

        std::fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn small_nve_validation_final_energy_is_regression_checked() {
        let config: RunConfig =
            toml::from_str(include_str!("../../../examples/validation/small-nve.toml")).unwrap();

        let final_sample = simulate_in_memory(&config).unwrap();

        assert_eq!(final_sample.step, 100);
        assert!((final_sample.total - -3.9901918753040113).abs() < 1.0e-12);
    }

    #[test]
    fn long_nve_validation_energy_drift_stays_bounded() {
        let config: RunConfig = toml::from_str(include_str!(
            "../../../examples/validation/long-nve-drift.toml"
        ))
        .unwrap();

        let (initial_sample, final_sample) = simulate_in_memory_energy_span(&config).unwrap();
        let total_drift = (final_sample.total - initial_sample.total).abs();
        let drift_per_time = total_drift / final_sample.time;

        assert_eq!(final_sample.step, 10_000);
        assert!(total_drift < 1.0e-3);
        assert!(drift_per_time < 1.0e-4);
    }

    fn minimal_test_config(seed: u64, output_directory: PathBuf) -> RunConfig {
        RunConfig {
            project: ProjectSection::default(),
            simulation: SimulationSection {
                steps: 20,
                dt: 0.001,
                output_interval: 5,
                temperature: 0.2,
            },
            input: None,
            topology: None,
            system: SystemSection {
                particles: 8,
                mass: 1.0,
                element: "Ar".to_string(),
                seed,
                lattice_spacing: Some(1.4),
            },
            simulation_box: BoxSection {
                x: 6.0,
                y: 6.0,
                z: 6.0,
                periodic: false,
            },
            force: ForceSection {
                kind: "lennard-jones".to_string(),
                epsilon: 1.0,
                sigma: 1.0,
                cutoff: 2.5,
                shift_potential: false,
            },
            force_field: ForceFieldSection::default(),
            coulomb: CoulombSection::default(),
            neighbor: NeighborSection {
                enabled: false,
                skin: default_neighbor_skin(),
                rebuild_interval: default_neighbor_rebuild_interval(),
            },
            thermostat: ThermostatSection::default(),
            minimization: MinimizationSection::default(),
            checkpoint: CheckpointSection::default(),
            execution: ExecutionSection { parallel: false },
            output: OutputSection {
                directory: output_directory,
                name: Some("deterministic".to_string()),
            },
        }
    }

    fn deterministic_config_toml(output_dir: &Path) -> String {
        format!(
            r#"
            [simulation]
            steps = 20
            dt = 0.001
            output_interval = 5
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

            [force]
            type = "lennard-jones"
            epsilon = 1.0
            sigma = 1.0
            cutoff = 2.5

            [output]
            directory = "{}"
            name = "deterministic"
            "#,
            toml_string(output_dir)
        )
    }

    fn binary_checkpoint_config_toml(output_dir: &Path) -> String {
        format!(
            r#"
            [simulation]
            steps = 10
            dt = 0.001
            output_interval = 5
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

            [force]
            type = "lennard-jones"
            epsilon = 1.0
            sigma = 1.0
            cutoff = 2.5

            [checkpoint]
            format = "binary"

            [output]
            directory = "{}"
            name = "binary-checkpoint"
            "#,
            toml_string(output_dir)
        )
    }

    fn simulate_in_memory(config: &RunConfig) -> Result<EnergySample> {
        Ok(simulate_in_memory_energy_span(config)?.1)
    }

    fn simulate_in_memory_energy_span(config: &RunConfig) -> Result<(EnergySample, EnergySample)> {
        config.validate()?;
        let (mut state, topology) = prepare_initial_state_and_topology(config, None)?;
        apply_topology_to_state(&mut state, topology.as_ref())?;
        let force_options = config.force_options()?;
        let mixed_force_options = config.mixed_lennard_jones_options(&state, &force_options)?;
        let coulomb_options = config.coulomb_options()?;
        let initial_potential_energy = compute_configured_forces(
            &mut state,
            &force_options,
            mixed_force_options.as_ref(),
            None,
            topology.as_ref(),
            coulomb_options.as_ref(),
            config.execution.parallel,
        )?;
        let integrator = VelocityVerlet::new(config.simulation.dt)?;
        let initial_sample = sample_energy(0, integrator.dt(), &state, initial_potential_energy);
        let mut final_sample = initial_sample;

        for step in 1..=config.simulation.steps {
            let potential_energy =
                integrator.step_with_force(&mut state, force_options.boundary, |state| {
                    compute_configured_forces(
                        state,
                        &force_options,
                        mixed_force_options.as_ref(),
                        None,
                        topology.as_ref(),
                        coulomb_options.as_ref(),
                        config.execution.parallel,
                    )
                    .map_err(IntegratorError::from)
                })?;
            final_sample = sample_energy(step, integrator.dt(), &state, potential_energy);
            ensure_finite_energy(&final_sample)?;
        }

        Ok((initial_sample, final_sample))
    }

    fn toml_string(path: &Path) -> String {
        path.to_string_lossy()
            .replace('\\', "\\\\")
            .replace('"', "\\\"")
    }

    fn unique_temp_dir(prefix: &str) -> PathBuf {
        let now = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        std::env::temp_dir().join(format!("{prefix}-{}-{now}", std::process::id()))
    }
}
