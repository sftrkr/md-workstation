//! Input and output helpers for trajectories, energy logs, and summaries.

use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Read, Write};
use std::path::Path;

use md_core::{EnergySample, SystemState};
use serde::{Deserialize, Serialize};
use thiserror::Error;

#[derive(Debug, Error)]
pub enum IoError {
    #[error(transparent)]
    Io(#[from] std::io::Error),
    #[error(transparent)]
    Json(#[from] serde_json::Error),
    #[error("failed to parse {format} at line {line}: {message}")]
    Parse {
        format: &'static str,
        line: usize,
        message: String,
    },
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct XyzAtom {
    pub element: String,
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct XyzFrame {
    pub step: usize,
    pub atoms: Vec<XyzAtom>,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct PdbAtom {
    pub serial: Option<i32>,
    pub atom_name: String,
    pub residue_name: Option<String>,
    pub residue_id: Option<i32>,
    pub chain_id: Option<String>,
    pub element: String,
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct PdbFrame {
    pub step: usize,
    pub atoms: Vec<PdbAtom>,
}

/// Write one XYZ trajectory frame.
pub fn write_xyz_frame<W: Write>(
    writer: &mut W,
    state: &SystemState,
    step: usize,
) -> Result<(), IoError> {
    writeln!(writer, "{}", state.particle_count())?;
    writeln!(writer, "step={step}")?;
    for i in 0..state.particle_count() {
        writeln!(
            writer,
            "{} {:.10} {:.10} {:.10}",
            state.element[i], state.x[i], state.y[i], state.z[i]
        )?;
    }
    Ok(())
}

pub fn write_atom_metadata_csv<W: Write>(
    writer: &mut W,
    state: &SystemState,
) -> Result<(), IoError> {
    writeln!(
        writer,
        "atom_index,element,atom_name,residue_name,residue_id,chain_id"
    )?;
    for index in 0..state.particle_count() {
        write!(writer, "{},", index)?;
        write_csv_field(writer, Some(&state.element[index]))?;
        write!(writer, ",")?;
        write_csv_field(writer, optional_string(&state.atom_name, index))?;
        write!(writer, ",")?;
        write_csv_field(writer, optional_string(&state.residue_name, index))?;
        write!(writer, ",")?;
        if let Some(residue_id) = optional_i32(&state.residue_id, index) {
            write!(writer, "{residue_id}")?;
        }
        write!(writer, ",")?;
        write_csv_field(writer, optional_string(&state.chain_id, index))?;
        writeln!(writer)?;
    }
    Ok(())
}

pub fn read_xyz_trajectory(path: impl AsRef<Path>) -> Result<Vec<XyzFrame>, IoError> {
    let file = File::open(path)?;
    read_xyz_trajectory_from_reader(file)
}

pub fn read_xyz_trajectory_from_reader<R: Read>(reader: R) -> Result<Vec<XyzFrame>, IoError> {
    read_xyz_frames_from_reader(reader, XyzStepPolicy::RequireExplicit)
}

/// Read standard XYZ coordinate frames.
///
/// This accepts normal molecule-style XYZ comments. If a comment line is not
/// in the form `step=<number>`, the frame index is used as the step value.
pub fn read_xyz_coordinate_frames(path: impl AsRef<Path>) -> Result<Vec<XyzFrame>, IoError> {
    let file = File::open(path)?;
    read_xyz_coordinate_frames_from_reader(file)
}

pub fn read_xyz_coordinate_frames_from_reader<R: Read>(
    reader: R,
) -> Result<Vec<XyzFrame>, IoError> {
    read_xyz_frames_from_reader(reader, XyzStepPolicy::InferFromFrameIndex)
}

/// Read a deliberately small PDB coordinate subset.
///
/// Supported records are `ATOM`, `HETATM`, `MODEL`, and `ENDMDL`. Coordinates
/// are read from fixed-width PDB columns, and metadata is preserved for atom
/// name, residue name, residue sequence ID, chain ID, and element when present.
pub fn read_pdb_coordinate_frames(path: impl AsRef<Path>) -> Result<Vec<PdbFrame>, IoError> {
    let file = File::open(path)?;
    read_pdb_coordinate_frames_from_reader(file)
}

pub fn read_pdb_coordinate_frames_from_reader<R: Read>(
    reader: R,
) -> Result<Vec<PdbFrame>, IoError> {
    let mut frames = Vec::new();
    let mut atoms = Vec::new();

    for (line_index, line) in BufReader::new(reader).lines().enumerate() {
        let line = line?;
        let line_number = line_index + 1;
        match pdb_record(&line) {
            "MODEL" => {
                if !atoms.is_empty() {
                    push_pdb_frame(&mut frames, &mut atoms);
                }
            }
            "ATOM" | "HETATM" => atoms.push(parse_pdb_atom(&line, line_number)?),
            "ENDMDL" => {
                if !atoms.is_empty() {
                    push_pdb_frame(&mut frames, &mut atoms);
                }
            }
            "END" => break,
            _ => {}
        }
    }

    if !atoms.is_empty() {
        push_pdb_frame(&mut frames, &mut atoms);
    }
    if frames.is_empty() {
        return Err(IoError::Parse {
            format: "pdb",
            line: 1,
            message: "expected at least one ATOM or HETATM record".to_string(),
        });
    }

    Ok(frames)
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum XyzStepPolicy {
    RequireExplicit,
    InferFromFrameIndex,
}

fn read_xyz_frames_from_reader<R: Read>(
    reader: R,
    step_policy: XyzStepPolicy,
) -> Result<Vec<XyzFrame>, IoError> {
    let mut lines = BufReader::new(reader).lines().enumerate();
    let mut frames = Vec::new();

    loop {
        let Some((count_line_index, count_line)) = lines.next() else {
            break;
        };
        let count_line = count_line?;
        if count_line.trim().is_empty() {
            continue;
        }

        let atom_count: usize = count_line.trim().parse().map_err(|error| IoError::Parse {
            format: "xyz",
            line: count_line_index + 1,
            message: format!("invalid atom count: {error}"),
        })?;

        let Some((comment_line_index, comment_line)) = lines.next() else {
            return Err(IoError::Parse {
                format: "xyz",
                line: count_line_index + 1,
                message: "missing frame comment line".to_string(),
            });
        };
        let comment_line = comment_line?;
        let step = match step_policy {
            XyzStepPolicy::RequireExplicit => parse_xyz_step(&comment_line, comment_line_index + 1),
            XyzStepPolicy::InferFromFrameIndex => {
                parse_xyz_step_or_frame_index(&comment_line, comment_line_index + 1, frames.len())
            }
        }?;

        let mut atoms = Vec::with_capacity(atom_count);
        for _ in 0..atom_count {
            let Some((atom_line_index, atom_line)) = lines.next() else {
                return Err(IoError::Parse {
                    format: "xyz",
                    line: comment_line_index + 1,
                    message: "missing atom line".to_string(),
                });
            };
            atoms.push(parse_xyz_atom(&atom_line?, atom_line_index + 1)?);
        }

        frames.push(XyzFrame { step, atoms });
    }

    Ok(frames)
}

/// CSV writer for simulation energy samples.
pub struct EnergyCsvWriter<W: Write> {
    writer: W,
}

pub fn read_energy_csv(path: impl AsRef<Path>) -> Result<Vec<EnergySample>, IoError> {
    let file = File::open(path)?;
    read_energy_csv_from_reader(file)
}

pub fn read_energy_csv_from_reader<R: Read>(reader: R) -> Result<Vec<EnergySample>, IoError> {
    let mut samples = Vec::new();
    for (line_index, line) in BufReader::new(reader).lines().enumerate() {
        let line = line?;
        if line_index == 0 {
            validate_energy_header(&line)?;
            continue;
        }
        if line.trim().is_empty() {
            continue;
        }
        samples.push(parse_energy_sample(&line, line_index + 1)?);
    }
    Ok(samples)
}

impl<W: Write> EnergyCsvWriter<W> {
    pub fn new(mut writer: W) -> Result<Self, IoError> {
        writeln!(writer, "step,time,kinetic,potential,total,temperature")?;
        Ok(Self { writer })
    }

    pub fn append(writer: W) -> Self {
        Self { writer }
    }

    pub fn write_sample(&mut self, sample: &EnergySample) -> Result<(), IoError> {
        writeln!(
            self.writer,
            "{},{:.10},{:.10},{:.10},{:.10},{:.10}",
            sample.step,
            sample.time,
            sample.kinetic,
            sample.potential,
            sample.total,
            sample.temperature
        )?;
        Ok(())
    }

    pub fn flush(&mut self) -> Result<(), IoError> {
        self.writer.flush()?;
        Ok(())
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RunSummary {
    pub run_name: String,
    #[serde(default)]
    pub project_name: Option<String>,
    #[serde(default)]
    pub project_description: Option<String>,
    pub particle_count: usize,
    pub steps: usize,
    pub dt: f64,
    pub output_interval: usize,
    pub seed: u64,
    #[serde(default = "default_workflow")]
    pub workflow: String,
    #[serde(default = "default_ensemble")]
    pub ensemble: String,
    #[serde(default = "default_thermostat")]
    pub thermostat: String,
    #[serde(default)]
    pub thermostat_target_temperature: Option<f64>,
    #[serde(default)]
    pub thermostat_tau: Option<f64>,
    #[serde(default)]
    pub minimization_step_size: Option<f64>,
    #[serde(default)]
    pub minimization_force_tolerance: Option<f64>,
    #[serde(default)]
    pub minimization_final_max_force: Option<f64>,
    pub input_file: Option<String>,
    pub input_format: Option<String>,
    pub input_frame: Option<usize>,
    #[serde(default)]
    pub atom_metadata_file: Option<String>,
    pub topology_file: Option<String>,
    #[serde(default)]
    pub bond_count: usize,
    #[serde(default)]
    pub angle_count: usize,
    #[serde(default)]
    pub dihedral_count: usize,
    #[serde(default)]
    pub excluded_pair_count: usize,
    #[serde(default)]
    pub coulomb: bool,
    pub coulomb_cutoff: Option<f64>,
    pub force: String,
    pub boundary: String,
    pub shifted_potential: bool,
    pub neighbor_list: bool,
    pub neighbor_skin: Option<f64>,
    pub neighbor_rebuild_interval: Option<usize>,
    pub parallel: bool,
    pub rayon_threads: Option<usize>,
    pub units: String,
    pub final_energy: EnergySample,
    pub outputs: Vec<String>,
}

fn default_workflow() -> String {
    "dynamics".to_string()
}

fn default_ensemble() -> String {
    "NVE".to_string()
}

fn default_thermostat() -> String {
    "none".to_string()
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct RunManifest {
    pub schema_version: u32,
    pub application: String,
    pub run_name: String,
    pub project_name: Option<String>,
    pub status: String,
    pub created_at_unix_seconds: u64,
    pub completed_at_unix_seconds: Option<u64>,
    pub config_file: String,
    pub summary_file: String,
    pub trajectory_file: String,
    pub energy_file: String,
    pub input_file: Option<String>,
    pub topology_file: Option<String>,
    #[serde(default)]
    pub config_hash: Option<String>,
    #[serde(default)]
    pub input_hash: Option<String>,
    #[serde(default)]
    pub topology_hash: Option<String>,
    #[serde(default)]
    pub engine_version: Option<String>,
    #[serde(default)]
    pub engine_git_commit: Option<String>,
    #[serde(default)]
    pub rust_target: Option<String>,
    #[serde(default)]
    pub platform: Option<String>,
    #[serde(default)]
    pub rayon_threads: Option<usize>,
    #[serde(default)]
    pub command_line: Option<Vec<String>>,
    #[serde(default)]
    pub checkpoint_file: Option<String>,
    #[serde(default)]
    pub checkpoint_format: Option<String>,
    pub analysis_summary_file: Option<String>,
    pub report_file: Option<String>,
    pub outputs: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct RunCheckpoint {
    pub schema_version: u32,
    pub run_name: String,
    pub step: usize,
    pub time: f64,
    pub potential_energy: f64,
    pub state: SystemState,
}

pub fn write_summary_json(path: impl AsRef<Path>, summary: &RunSummary) -> Result<(), IoError> {
    let file = File::create(path)?;
    let writer = BufWriter::new(file);
    serde_json::to_writer_pretty(writer, summary)?;
    Ok(())
}

pub fn read_summary_json(path: impl AsRef<Path>) -> Result<RunSummary, IoError> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    Ok(serde_json::from_reader(reader)?)
}

pub fn write_manifest_json(path: impl AsRef<Path>, manifest: &RunManifest) -> Result<(), IoError> {
    let file = File::create(path)?;
    let writer = BufWriter::new(file);
    serde_json::to_writer_pretty(writer, manifest)?;
    Ok(())
}

pub fn read_manifest_json(path: impl AsRef<Path>) -> Result<RunManifest, IoError> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    Ok(serde_json::from_reader(reader)?)
}

pub fn write_checkpoint_json(
    path: impl AsRef<Path>,
    checkpoint: &RunCheckpoint,
) -> Result<(), IoError> {
    let file = File::create(path)?;
    let writer = BufWriter::new(file);
    serde_json::to_writer_pretty(writer, checkpoint)?;
    Ok(())
}

pub fn read_checkpoint_json(path: impl AsRef<Path>) -> Result<RunCheckpoint, IoError> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    Ok(serde_json::from_reader(reader)?)
}

const CHECKPOINT_BINARY_MAGIC: &[u8; 8] = b"MDCKPT01";

pub fn write_checkpoint_binary(
    path: impl AsRef<Path>,
    checkpoint: &RunCheckpoint,
) -> Result<(), IoError> {
    checkpoint
        .state
        .validate()
        .map_err(|error| IoError::Parse {
            format: "checkpoint-binary",
            line: 0,
            message: error.to_string(),
        })?;

    let file = File::create(path)?;
    let mut writer = BufWriter::new(file);
    writer.write_all(CHECKPOINT_BINARY_MAGIC)?;
    write_u32(&mut writer, checkpoint.schema_version)?;
    write_string(&mut writer, &checkpoint.run_name)?;
    write_u64(&mut writer, checkpoint.step as u64)?;
    write_f64(&mut writer, checkpoint.time)?;
    write_f64(&mut writer, checkpoint.potential_energy)?;
    write_system_state_binary(&mut writer, &checkpoint.state)?;
    writer.flush()?;
    Ok(())
}

pub fn read_checkpoint_binary(path: impl AsRef<Path>) -> Result<RunCheckpoint, IoError> {
    let file = File::open(path)?;
    let mut reader = BufReader::new(file);
    let mut magic = [0; 8];
    reader.read_exact(&mut magic)?;
    if &magic != CHECKPOINT_BINARY_MAGIC {
        return Err(IoError::Parse {
            format: "checkpoint-binary",
            line: 0,
            message: "invalid checkpoint magic".to_string(),
        });
    }

    let checkpoint = RunCheckpoint {
        schema_version: read_u32(&mut reader)?,
        run_name: read_string(&mut reader)?,
        step: read_usize(&mut reader, "step")?,
        time: read_f64(&mut reader)?,
        potential_energy: read_f64(&mut reader)?,
        state: read_system_state_binary(&mut reader)?,
    };
    checkpoint
        .state
        .validate()
        .map_err(|error| IoError::Parse {
            format: "checkpoint-binary",
            line: 0,
            message: error.to_string(),
        })?;
    Ok(checkpoint)
}

fn write_system_state_binary<W: Write>(writer: &mut W, state: &SystemState) -> Result<(), IoError> {
    let count = state.particle_count();
    write_u64(writer, count as u64)?;
    write_f64_vec(writer, &state.x)?;
    write_f64_vec(writer, &state.y)?;
    write_f64_vec(writer, &state.z)?;
    write_f64_vec(writer, &state.vx)?;
    write_f64_vec(writer, &state.vy)?;
    write_f64_vec(writer, &state.vz)?;
    write_f64_vec(writer, &state.fx)?;
    write_f64_vec(writer, &state.fy)?;
    write_f64_vec(writer, &state.fz)?;
    write_f64_vec(writer, &state.mass)?;
    write_f64_vec(writer, &state.charge)?;
    write_string_vec(writer, &state.element)?;
    write_optional_string_vec(writer, &state.atom_name, count)?;
    write_optional_string_vec(writer, &state.residue_name, count)?;
    write_optional_i32_vec(writer, &state.residue_id, count)?;
    write_optional_string_vec(writer, &state.chain_id, count)?;
    Ok(())
}

fn read_system_state_binary<R: Read>(reader: &mut R) -> Result<SystemState, IoError> {
    let count = read_usize(reader, "particle_count")?;
    Ok(SystemState {
        x: read_f64_vec(reader, count, "x")?,
        y: read_f64_vec(reader, count, "y")?,
        z: read_f64_vec(reader, count, "z")?,
        vx: read_f64_vec(reader, count, "vx")?,
        vy: read_f64_vec(reader, count, "vy")?,
        vz: read_f64_vec(reader, count, "vz")?,
        fx: read_f64_vec(reader, count, "fx")?,
        fy: read_f64_vec(reader, count, "fy")?,
        fz: read_f64_vec(reader, count, "fz")?,
        mass: read_f64_vec(reader, count, "mass")?,
        charge: read_f64_vec(reader, count, "charge")?,
        element: read_string_vec(reader, count, "element")?,
        atom_name: read_optional_string_vec(reader, count, "atom_name")?,
        residue_name: read_optional_string_vec(reader, count, "residue_name")?,
        residue_id: read_optional_i32_vec(reader, count, "residue_id")?,
        chain_id: read_optional_string_vec(reader, count, "chain_id")?,
    })
}

fn write_f64_vec<W: Write>(writer: &mut W, values: &[f64]) -> Result<(), IoError> {
    for value in values {
        write_f64(writer, *value)?;
    }
    Ok(())
}

fn read_f64_vec<R: Read>(
    reader: &mut R,
    count: usize,
    field: &'static str,
) -> Result<Vec<f64>, IoError> {
    let mut values = Vec::with_capacity(count);
    for _ in 0..count {
        values.push(read_f64(reader).map_err(|error| IoError::Parse {
            format: "checkpoint-binary",
            line: 0,
            message: format!("failed to read {field}: {error}"),
        })?);
    }
    Ok(values)
}

fn write_string_vec<W: Write>(writer: &mut W, values: &[String]) -> Result<(), IoError> {
    for value in values {
        write_string(writer, value)?;
    }
    Ok(())
}

fn read_string_vec<R: Read>(
    reader: &mut R,
    count: usize,
    field: &'static str,
) -> Result<Vec<String>, IoError> {
    let mut values = Vec::with_capacity(count);
    for _ in 0..count {
        values.push(read_string(reader).map_err(|error| IoError::Parse {
            format: "checkpoint-binary",
            line: 0,
            message: format!("failed to read {field}: {error}"),
        })?);
    }
    Ok(values)
}

fn write_optional_string_vec<W: Write>(
    writer: &mut W,
    values: &[Option<String>],
    count: usize,
) -> Result<(), IoError> {
    for index in 0..count {
        match values.get(index).and_then(|value| value.as_deref()) {
            Some(value) => {
                write_bool(writer, true)?;
                write_string(writer, value)?;
            }
            None => write_bool(writer, false)?,
        }
    }
    Ok(())
}

fn read_optional_string_vec<R: Read>(
    reader: &mut R,
    count: usize,
    field: &'static str,
) -> Result<Vec<Option<String>>, IoError> {
    let mut values = Vec::with_capacity(count);
    for _ in 0..count {
        values.push(if read_bool(reader)? {
            Some(read_string(reader).map_err(|error| IoError::Parse {
                format: "checkpoint-binary",
                line: 0,
                message: format!("failed to read {field}: {error}"),
            })?)
        } else {
            None
        });
    }
    Ok(values)
}

fn write_optional_i32_vec<W: Write>(
    writer: &mut W,
    values: &[Option<i32>],
    count: usize,
) -> Result<(), IoError> {
    for index in 0..count {
        match values.get(index).copied().flatten() {
            Some(value) => {
                write_bool(writer, true)?;
                writer.write_all(&value.to_le_bytes())?;
            }
            None => write_bool(writer, false)?,
        }
    }
    Ok(())
}

fn read_optional_i32_vec<R: Read>(
    reader: &mut R,
    count: usize,
    field: &'static str,
) -> Result<Vec<Option<i32>>, IoError> {
    let mut values = Vec::with_capacity(count);
    for _ in 0..count {
        values.push(if read_bool(reader)? {
            Some(read_i32(reader).map_err(|error| IoError::Parse {
                format: "checkpoint-binary",
                line: 0,
                message: format!("failed to read {field}: {error}"),
            })?)
        } else {
            None
        });
    }
    Ok(values)
}

fn write_string<W: Write>(writer: &mut W, value: &str) -> Result<(), IoError> {
    write_u64(writer, value.len() as u64)?;
    writer.write_all(value.as_bytes())?;
    Ok(())
}

fn read_string<R: Read>(reader: &mut R) -> Result<String, IoError> {
    let len = read_usize(reader, "string length")?;
    let mut bytes = vec![0; len];
    reader.read_exact(&mut bytes)?;
    String::from_utf8(bytes).map_err(|error| IoError::Parse {
        format: "checkpoint-binary",
        line: 0,
        message: format!("invalid UTF-8 string: {error}"),
    })
}

fn write_bool<W: Write>(writer: &mut W, value: bool) -> Result<(), IoError> {
    writer.write_all(&[u8::from(value)])?;
    Ok(())
}

fn read_bool<R: Read>(reader: &mut R) -> Result<bool, IoError> {
    let mut value = [0; 1];
    reader.read_exact(&mut value)?;
    match value[0] {
        0 => Ok(false),
        1 => Ok(true),
        other => Err(IoError::Parse {
            format: "checkpoint-binary",
            line: 0,
            message: format!("invalid boolean tag {other}"),
        }),
    }
}

fn write_u32<W: Write>(writer: &mut W, value: u32) -> Result<(), IoError> {
    writer.write_all(&value.to_le_bytes())?;
    Ok(())
}

fn read_u32<R: Read>(reader: &mut R) -> Result<u32, IoError> {
    let mut bytes = [0; 4];
    reader.read_exact(&mut bytes)?;
    Ok(u32::from_le_bytes(bytes))
}

fn write_u64<W: Write>(writer: &mut W, value: u64) -> Result<(), IoError> {
    writer.write_all(&value.to_le_bytes())?;
    Ok(())
}

fn read_u64<R: Read>(reader: &mut R) -> Result<u64, IoError> {
    let mut bytes = [0; 8];
    reader.read_exact(&mut bytes)?;
    Ok(u64::from_le_bytes(bytes))
}

fn read_usize<R: Read>(reader: &mut R, field: &'static str) -> Result<usize, IoError> {
    let value = read_u64(reader)?;
    usize::try_from(value).map_err(|_| IoError::Parse {
        format: "checkpoint-binary",
        line: 0,
        message: format!("{field} is too large for this platform"),
    })
}

fn write_f64<W: Write>(writer: &mut W, value: f64) -> Result<(), IoError> {
    writer.write_all(&value.to_le_bytes())?;
    Ok(())
}

fn read_f64<R: Read>(reader: &mut R) -> Result<f64, IoError> {
    let mut bytes = [0; 8];
    reader.read_exact(&mut bytes)?;
    Ok(f64::from_le_bytes(bytes))
}

fn read_i32<R: Read>(reader: &mut R) -> Result<i32, IoError> {
    let mut bytes = [0; 4];
    reader.read_exact(&mut bytes)?;
    Ok(i32::from_le_bytes(bytes))
}

fn parse_xyz_step(comment: &str, line: usize) -> Result<usize, IoError> {
    let trimmed = comment.trim();
    let step_value = trimmed
        .strip_prefix("step=")
        .ok_or_else(|| IoError::Parse {
            format: "xyz",
            line,
            message: "expected comment line in the form step=<number>".to_string(),
        })?;
    step_value.parse().map_err(|error| IoError::Parse {
        format: "xyz",
        line,
        message: format!("invalid step value: {error}"),
    })
}

fn parse_xyz_step_or_frame_index(
    comment: &str,
    line: usize,
    frame_index: usize,
) -> Result<usize, IoError> {
    if comment.trim().starts_with("step=") {
        parse_xyz_step(comment, line)
    } else {
        Ok(frame_index)
    }
}

fn parse_xyz_atom(line: &str, line_number: usize) -> Result<XyzAtom, IoError> {
    let fields: Vec<&str> = line.split_whitespace().collect();
    if fields.len() != 4 {
        return Err(IoError::Parse {
            format: "xyz",
            line: line_number,
            message: "expected element x y z".to_string(),
        });
    }
    Ok(XyzAtom {
        element: fields[0].to_string(),
        x: parse_f64(fields[1], "xyz", line_number, "x")?,
        y: parse_f64(fields[2], "xyz", line_number, "y")?,
        z: parse_f64(fields[3], "xyz", line_number, "z")?,
    })
}

fn parse_pdb_atom(line: &str, line_number: usize) -> Result<PdbAtom, IoError> {
    let atom_name = pdb_field(line, 12, 16).to_string();
    if atom_name.is_empty() {
        return Err(IoError::Parse {
            format: "pdb",
            line: line_number,
            message: "atom name is empty".to_string(),
        });
    }
    let element =
        pdb_element(pdb_field(line, 76, 78), &atom_name).ok_or_else(|| IoError::Parse {
            format: "pdb",
            line: line_number,
            message: "element is missing and could not be inferred from atom name".to_string(),
        })?;

    Ok(PdbAtom {
        serial: parse_optional_i32(pdb_field(line, 6, 11), "pdb", line_number, "serial")?,
        atom_name,
        residue_name: optional_pdb_string(pdb_field(line, 17, 20)),
        residue_id: parse_optional_i32(pdb_field(line, 22, 26), "pdb", line_number, "residue_id")?,
        chain_id: optional_pdb_string(pdb_field(line, 21, 22)),
        element,
        x: parse_f64(pdb_field(line, 30, 38), "pdb", line_number, "x")?,
        y: parse_f64(pdb_field(line, 38, 46), "pdb", line_number, "y")?,
        z: parse_f64(pdb_field(line, 46, 54), "pdb", line_number, "z")?,
    })
}

fn push_pdb_frame(frames: &mut Vec<PdbFrame>, atoms: &mut Vec<PdbAtom>) {
    frames.push(PdbFrame {
        step: frames.len(),
        atoms: std::mem::take(atoms),
    });
}

fn pdb_record(line: &str) -> &str {
    pdb_field(line, 0, 6)
}

fn pdb_field(line: &str, start: usize, end: usize) -> &str {
    if line.len() <= start {
        ""
    } else {
        let end = end.min(line.len());
        line[start..end].trim()
    }
}

fn optional_pdb_string(value: &str) -> Option<String> {
    (!value.trim().is_empty()).then(|| value.trim().to_string())
}

fn pdb_element(element_field: &str, atom_name: &str) -> Option<String> {
    let raw = if element_field.trim().is_empty() {
        atom_name
            .chars()
            .find(|character| character.is_ascii_alphabetic())?
            .to_string()
    } else {
        element_field.trim().to_string()
    };
    let mut chars = raw
        .chars()
        .filter(|character| character.is_ascii_alphabetic());
    let first = chars.next()?.to_ascii_uppercase();
    let second = chars.next().map(|character| character.to_ascii_lowercase());
    Some(match second {
        Some(second) => format!("{first}{second}"),
        None => first.to_string(),
    })
}

fn optional_string(values: &[Option<String>], index: usize) -> Option<&str> {
    values.get(index).and_then(|value| value.as_deref())
}

fn optional_i32(values: &[Option<i32>], index: usize) -> Option<i32> {
    values.get(index).copied().flatten()
}

fn write_csv_field<W: Write>(writer: &mut W, value: Option<&str>) -> Result<(), IoError> {
    let value = value.unwrap_or("");
    if value.contains(',') || value.contains('"') || value.contains('\n') {
        write!(writer, "\"")?;
        for character in value.chars() {
            if character == '"' {
                write!(writer, "\"\"")?;
            } else {
                write!(writer, "{character}")?;
            }
        }
        write!(writer, "\"")?;
    } else {
        write!(writer, "{value}")?;
    }
    Ok(())
}

fn validate_energy_header(line: &str) -> Result<(), IoError> {
    let expected = "step,time,kinetic,potential,total,temperature";
    if line.trim() == expected {
        Ok(())
    } else {
        Err(IoError::Parse {
            format: "energy-csv",
            line: 1,
            message: format!("expected header {expected:?}"),
        })
    }
}

fn parse_energy_sample(line: &str, line_number: usize) -> Result<EnergySample, IoError> {
    let fields: Vec<&str> = line.split(',').collect();
    if fields.len() != 6 {
        return Err(IoError::Parse {
            format: "energy-csv",
            line: line_number,
            message: "expected 6 comma-separated fields".to_string(),
        });
    }
    Ok(EnergySample {
        step: parse_usize(fields[0], "energy-csv", line_number, "step")?,
        time: parse_f64(fields[1], "energy-csv", line_number, "time")?,
        kinetic: parse_f64(fields[2], "energy-csv", line_number, "kinetic")?,
        potential: parse_f64(fields[3], "energy-csv", line_number, "potential")?,
        total: parse_f64(fields[4], "energy-csv", line_number, "total")?,
        temperature: parse_f64(fields[5], "energy-csv", line_number, "temperature")?,
    })
}

fn parse_usize(
    value: &str,
    format: &'static str,
    line: usize,
    field: &'static str,
) -> Result<usize, IoError> {
    value.trim().parse().map_err(|error| IoError::Parse {
        format,
        line,
        message: format!("invalid {field}: {error}"),
    })
}

fn parse_optional_i32(
    value: &str,
    format: &'static str,
    line: usize,
    field: &'static str,
) -> Result<Option<i32>, IoError> {
    let trimmed = value.trim();
    if trimmed.is_empty() {
        Ok(None)
    } else {
        trimmed.parse().map(Some).map_err(|error| IoError::Parse {
            format,
            line,
            message: format!("invalid {field}: {error}"),
        })
    }
}

fn parse_f64(
    value: &str,
    format: &'static str,
    line: usize,
    field: &'static str,
) -> Result<f64, IoError> {
    value.trim().parse().map_err(|error| IoError::Parse {
        format,
        line,
        message: format!("invalid {field}: {error}"),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn xyz_writer_emits_particle_count_and_atoms() {
        let mut state = SystemState::new(1);
        state.element[0] = "Ar".to_string();
        state.x[0] = 1.0;
        state.y[0] = 2.0;
        state.z[0] = 3.0;

        let mut output = Vec::new();
        write_xyz_frame(&mut output, &state, 7).unwrap();
        let output = String::from_utf8(output).unwrap();

        assert!(output.starts_with("1\nstep=7\n"));
        assert!(output.contains("Ar 1.0000000000 2.0000000000 3.0000000000"));
    }

    #[test]
    fn energy_writer_emits_header_and_sample() {
        let mut output = Vec::new();
        let mut writer = EnergyCsvWriter::new(&mut output).unwrap();
        writer
            .write_sample(&EnergySample::new(1, 0.1, 2.0, -1.0, 0.5))
            .unwrap();
        writer.flush().unwrap();
        drop(writer);

        let output = String::from_utf8(output).unwrap();
        assert!(output.contains("step,time,kinetic,potential,total,temperature"));
        assert!(
            output.contains("1,0.1000000000,2.0000000000,-1.0000000000,1.0000000000,0.5000000000")
        );
    }

    #[test]
    fn xyz_reader_parses_multiple_frames() {
        let input = b"1\nstep=0\nAr 1.0 2.0 3.0\n1\nstep=5\nAr 1.5 2.5 3.5\n";

        let frames = read_xyz_trajectory_from_reader(&input[..]).unwrap();

        assert_eq!(frames.len(), 2);
        assert_eq!(frames[1].step, 5);
        assert_eq!(frames[1].atoms[0].element, "Ar");
        assert!((frames[1].atoms[0].x - 1.5).abs() < 1.0e-12);
    }

    #[test]
    fn xyz_coordinate_reader_accepts_molecule_comments() {
        let input = b"2\nargon dimer\nAr 1.0 2.0 3.0\nAr 2.5 2.0 3.0\n2\nsecond frame\nAr 1.1 2.0 3.0\nAr 2.6 2.0 3.0\n";

        let frames = read_xyz_coordinate_frames_from_reader(&input[..]).unwrap();

        assert_eq!(frames.len(), 2);
        assert_eq!(frames[0].step, 0);
        assert_eq!(frames[1].step, 1);
        assert_eq!(frames[0].atoms[0].element, "Ar");
        assert!((frames[1].atoms[1].x - 2.6).abs() < 1.0e-12);
    }

    #[test]
    fn xyz_coordinate_reader_preserves_explicit_step_comments() {
        let input = b"1\nstep=12\nAr 1.0 2.0 3.0\n";

        let frames = read_xyz_coordinate_frames_from_reader(&input[..]).unwrap();

        assert_eq!(frames.len(), 1);
        assert_eq!(frames[0].step, 12);
    }

    #[test]
    fn pdb_coordinate_reader_preserves_subset_metadata() {
        let input = b"MODEL        1\nATOM      1  N   ALA A   7       1.000   2.000   3.000  1.00 10.00           N  \nATOM      2  CA  ALA A   7       2.000   2.100   3.200  1.00 10.00           C  \nENDMDL\nMODEL        2\nHETATM    3  O   HOH B  12       3.000   4.000   5.000  1.00 10.00           O  \nENDMDL\n";

        let frames = read_pdb_coordinate_frames_from_reader(&input[..]).unwrap();

        assert_eq!(frames.len(), 2);
        assert_eq!(frames[0].atoms.len(), 2);
        assert_eq!(frames[0].atoms[0].serial, Some(1));
        assert_eq!(frames[0].atoms[0].atom_name, "N");
        assert_eq!(frames[0].atoms[0].residue_name.as_deref(), Some("ALA"));
        assert_eq!(frames[0].atoms[0].residue_id, Some(7));
        assert_eq!(frames[0].atoms[0].chain_id.as_deref(), Some("A"));
        assert_eq!(frames[0].atoms[0].element, "N");
        assert!((frames[1].atoms[0].z - 5.0).abs() < 1.0e-12);
    }

    #[test]
    fn atom_metadata_writer_emits_optional_pdb_fields() {
        let mut state = SystemState::new(1);
        state.element[0] = "C".to_string();
        state.atom_name[0] = Some("CA".to_string());
        state.residue_name[0] = Some("ALA".to_string());
        state.residue_id[0] = Some(7);
        state.chain_id[0] = Some("A".to_string());

        let mut output = Vec::new();
        write_atom_metadata_csv(&mut output, &state).unwrap();
        let output = String::from_utf8(output).unwrap();

        assert!(output.contains("atom_index,element,atom_name,residue_name,residue_id,chain_id"));
        assert!(output.contains("0,C,CA,ALA,7,A"));
    }

    #[test]
    fn binary_checkpoint_round_trips_state_and_metadata() {
        let mut state = SystemState::new(2);
        state.x = vec![1.0, 2.0];
        state.y = vec![3.0, 4.0];
        state.z = vec![5.0, 6.0];
        state.vx = vec![0.1, 0.2];
        state.vy = vec![0.3, 0.4];
        state.vz = vec![0.5, 0.6];
        state.fx = vec![-1.0, 1.0];
        state.fy = vec![-2.0, 2.0];
        state.fz = vec![-3.0, 3.0];
        state.mass = vec![12.0, 16.0];
        state.charge = vec![0.1, -0.1];
        state.element = vec!["C".to_string(), "O".to_string()];
        state.atom_name = vec![Some("CA".to_string()), Some("O".to_string())];
        state.residue_name = vec![Some("ALA".to_string()), Some("ALA".to_string())];
        state.residue_id = vec![Some(7), Some(7)];
        state.chain_id = vec![Some("A".to_string()), Some("A".to_string())];
        let checkpoint = RunCheckpoint {
            schema_version: 1,
            run_name: "binary-test".to_string(),
            step: 12,
            time: 0.024,
            potential_energy: -1.5,
            state,
        };
        let path = unique_temp_path("md-io-checkpoint", "bin");

        write_checkpoint_binary(&path, &checkpoint).unwrap();
        let decoded = read_checkpoint_binary(&path).unwrap();

        assert_eq!(decoded, checkpoint);
        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn energy_reader_parses_samples() {
        let input = b"step,time,kinetic,potential,total,temperature\n0,0.0,2.0,-1.0,1.0,0.5\n";

        let samples = read_energy_csv_from_reader(&input[..]).unwrap();

        assert_eq!(samples.len(), 1);
        assert_eq!(samples[0].step, 0);
        assert!((samples[0].total - 1.0).abs() < 1.0e-12);
    }

    fn unique_temp_path(prefix: &str, extension: &str) -> std::path::PathBuf {
        let now = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        std::env::temp_dir().join(format!("{prefix}-{}-{now}.{extension}", std::process::id()))
    }
}
