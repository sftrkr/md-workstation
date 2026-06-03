//! Analysis tools for MD run outputs.

use md_core::{EnergySample, SimulationBox};
use md_io::{XyzAtom, XyzFrame};
use serde::{Deserialize, Serialize};
use thiserror::Error;

#[derive(Debug, Error)]
pub enum AnalysisError {
    #[error("analysis requires at least one energy sample")]
    EmptyEnergySeries,
    #[error("analysis requires at least one trajectory frame")]
    EmptyTrajectory,
    #[error("frame {frame} has no atoms")]
    EmptyAtomFrame { frame: usize },
    #[error("frame {frame} has {actual} atoms, expected {expected}")]
    AtomCountMismatch {
        frame: usize,
        expected: usize,
        actual: usize,
    },
    #[error("periodic unwrap box dimensions must be positive and finite")]
    InvalidPeriodicBox,
    #[error("reference frame index {index} is out of range for {frame_count} frames")]
    ReferenceFrameOutOfRange { index: usize, frame_count: usize },
    #[error("atom index {index} is out of range for {atom_count} atoms")]
    AtomIndexOutOfRange { index: usize, atom_count: usize },
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct EnergySummary {
    pub sample_count: usize,
    pub first_step: usize,
    pub last_step: usize,
    pub first_time: f64,
    pub last_time: f64,
    pub initial_total: f64,
    pub final_total: f64,
    pub total_drift: f64,
    pub total_drift_per_time: Option<f64>,
    pub min_total: f64,
    pub max_total: f64,
    pub mean_total: f64,
    pub min_temperature: f64,
    pub max_temperature: f64,
    pub mean_temperature: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct RmsdSample {
    pub frame_index: usize,
    pub step: usize,
    pub rmsd: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct RmsdSummary {
    pub reference_frame: usize,
    pub sample_count: usize,
    pub min: f64,
    pub max: f64,
    pub mean: f64,
    pub final_value: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct PairDistanceSample {
    pub frame_index: usize,
    pub step: usize,
    pub atom_i: usize,
    pub atom_j: usize,
    pub distance: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct PairDistanceSummary {
    pub atom_i: usize,
    pub atom_j: usize,
    pub sample_count: usize,
    pub min: f64,
    pub max: f64,
    pub mean: f64,
    pub final_value: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct AnalysisSummary {
    pub frame_count: usize,
    pub atom_count: usize,
    #[serde(default)]
    pub periodic_unwrap: bool,
    #[serde(default)]
    pub periodic_box: Option<[f64; 3]>,
    #[serde(default = "default_rmsd_alignment")]
    pub rmsd_alignment: String,
    pub energy: EnergySummary,
    pub rmsd: RmsdSummary,
    pub pair_distance: Option<PairDistanceSummary>,
}

pub fn summarize_energy(samples: &[EnergySample]) -> Result<EnergySummary, AnalysisError> {
    let first = samples.first().ok_or(AnalysisError::EmptyEnergySeries)?;
    let last = samples.last().ok_or(AnalysisError::EmptyEnergySeries)?;

    let mut min_total = f64::INFINITY;
    let mut max_total = f64::NEG_INFINITY;
    let mut sum_total = 0.0;
    let mut min_temperature = f64::INFINITY;
    let mut max_temperature = f64::NEG_INFINITY;
    let mut sum_temperature = 0.0;

    for sample in samples {
        min_total = min_total.min(sample.total);
        max_total = max_total.max(sample.total);
        sum_total += sample.total;
        min_temperature = min_temperature.min(sample.temperature);
        max_temperature = max_temperature.max(sample.temperature);
        sum_temperature += sample.temperature;
    }

    let duration = last.time - first.time;
    let total_drift = last.total - first.total;

    Ok(EnergySummary {
        sample_count: samples.len(),
        first_step: first.step,
        last_step: last.step,
        first_time: first.time,
        last_time: last.time,
        initial_total: first.total,
        final_total: last.total,
        total_drift,
        total_drift_per_time: (duration != 0.0).then_some(total_drift / duration),
        min_total,
        max_total,
        mean_total: sum_total / samples.len() as f64,
        min_temperature,
        max_temperature,
        mean_temperature: sum_temperature / samples.len() as f64,
    })
}

pub fn rmsd_series(
    frames: &[XyzFrame],
    reference_frame: usize,
) -> Result<Vec<RmsdSample>, AnalysisError> {
    validate_frames(frames)?;
    let reference = frames
        .get(reference_frame)
        .ok_or(AnalysisError::ReferenceFrameOutOfRange {
            index: reference_frame,
            frame_count: frames.len(),
        })?;

    let mut samples = Vec::with_capacity(frames.len());
    for (frame_index, frame) in frames.iter().enumerate() {
        let rmsd = if frame_index == reference_frame {
            0.0
        } else {
            aligned_rmsd(frame, reference)
        };
        samples.push(RmsdSample {
            frame_index,
            step: frame.step,
            rmsd,
        });
    }

    Ok(samples)
}

pub fn unwrap_periodic_frames(
    frames: &[XyzFrame],
    simulation_box: SimulationBox,
) -> Result<Vec<XyzFrame>, AnalysisError> {
    validate_frames(frames)?;
    simulation_box
        .validate()
        .map_err(|_| AnalysisError::InvalidPeriodicBox)?;

    let Some(first) = frames.first() else {
        return Err(AnalysisError::EmptyTrajectory);
    };
    let mut unwrapped = Vec::with_capacity(frames.len());
    unwrapped.push(first.clone());

    for frame_index in 1..frames.len() {
        let previous_wrapped = &frames[frame_index - 1];
        let previous_unwrapped = &unwrapped[frame_index - 1];
        let current_wrapped = &frames[frame_index];
        let mut atoms = Vec::with_capacity(current_wrapped.atoms.len());

        for atom_index in 0..current_wrapped.atoms.len() {
            let previous_raw = &previous_wrapped.atoms[atom_index];
            let previous_continuous = &previous_unwrapped.atoms[atom_index];
            let current = &current_wrapped.atoms[atom_index];
            let (dx, dy, dz) = simulation_box.minimum_image_delta(
                current.x - previous_raw.x,
                current.y - previous_raw.y,
                current.z - previous_raw.z,
            );
            atoms.push(XyzAtom {
                element: current.element.clone(),
                x: previous_continuous.x + dx,
                y: previous_continuous.y + dy,
                z: previous_continuous.z + dz,
            });
        }

        unwrapped.push(XyzFrame {
            step: current_wrapped.step,
            atoms,
        });
    }

    Ok(unwrapped)
}

pub fn summarize_rmsd(
    samples: &[RmsdSample],
    reference_frame: usize,
) -> Result<RmsdSummary, AnalysisError> {
    if samples.is_empty() {
        return Err(AnalysisError::EmptyTrajectory);
    }

    let values: Vec<f64> = samples.iter().map(|sample| sample.rmsd).collect();
    Ok(RmsdSummary {
        reference_frame,
        sample_count: samples.len(),
        min: min_value(&values),
        max: max_value(&values),
        mean: mean_value(&values),
        final_value: values[values.len() - 1],
    })
}

pub fn pair_distance_series(
    frames: &[XyzFrame],
    atom_i: usize,
    atom_j: usize,
) -> Result<Vec<PairDistanceSample>, AnalysisError> {
    validate_frames(frames)?;
    let atom_count = frames[0].atoms.len();
    if atom_i >= atom_count {
        return Err(AnalysisError::AtomIndexOutOfRange {
            index: atom_i,
            atom_count,
        });
    }
    if atom_j >= atom_count {
        return Err(AnalysisError::AtomIndexOutOfRange {
            index: atom_j,
            atom_count,
        });
    }

    let mut samples = Vec::with_capacity(frames.len());
    for (frame_index, frame) in frames.iter().enumerate() {
        let a = &frame.atoms[atom_i];
        let b = &frame.atoms[atom_j];
        let dx = a.x - b.x;
        let dy = a.y - b.y;
        let dz = a.z - b.z;
        samples.push(PairDistanceSample {
            frame_index,
            step: frame.step,
            atom_i,
            atom_j,
            distance: (dx * dx + dy * dy + dz * dz).sqrt(),
        });
    }

    Ok(samples)
}

pub fn summarize_pair_distances(
    samples: &[PairDistanceSample],
) -> Result<PairDistanceSummary, AnalysisError> {
    if samples.is_empty() {
        return Err(AnalysisError::EmptyTrajectory);
    }

    let values: Vec<f64> = samples.iter().map(|sample| sample.distance).collect();
    let first = &samples[0];
    Ok(PairDistanceSummary {
        atom_i: first.atom_i,
        atom_j: first.atom_j,
        sample_count: samples.len(),
        min: min_value(&values),
        max: max_value(&values),
        mean: mean_value(&values),
        final_value: values[values.len() - 1],
    })
}

pub fn summarize_analysis(
    energy_samples: &[EnergySample],
    frames: &[XyzFrame],
    reference_frame: usize,
    pair_distance_samples: Option<&[PairDistanceSample]>,
) -> Result<AnalysisSummary, AnalysisError> {
    summarize_analysis_with_metadata(
        energy_samples,
        frames,
        reference_frame,
        pair_distance_samples,
        false,
        None,
    )
}

pub fn summarize_analysis_with_metadata(
    energy_samples: &[EnergySample],
    frames: &[XyzFrame],
    reference_frame: usize,
    pair_distance_samples: Option<&[PairDistanceSample]>,
    periodic_unwrap: bool,
    periodic_box: Option<SimulationBox>,
) -> Result<AnalysisSummary, AnalysisError> {
    validate_frames(frames)?;
    let rmsd_samples = rmsd_series(frames, reference_frame)?;
    Ok(AnalysisSummary {
        frame_count: frames.len(),
        atom_count: frames[0].atoms.len(),
        periodic_unwrap,
        periodic_box: periodic_box
            .map(|simulation_box| [simulation_box.x, simulation_box.y, simulation_box.z]),
        rmsd_alignment: default_rmsd_alignment(),
        energy: summarize_energy(energy_samples)?,
        rmsd: summarize_rmsd(&rmsd_samples, reference_frame)?,
        pair_distance: pair_distance_samples
            .map(summarize_pair_distances)
            .transpose()?,
    })
}

fn validate_frames(frames: &[XyzFrame]) -> Result<(), AnalysisError> {
    let first = frames.first().ok_or(AnalysisError::EmptyTrajectory)?;
    let expected = first.atoms.len();
    if expected == 0 {
        return Err(AnalysisError::EmptyAtomFrame { frame: 0 });
    }
    for (frame_index, frame) in frames.iter().enumerate() {
        if frame.atoms.len() != expected {
            return Err(AnalysisError::AtomCountMismatch {
                frame: frame_index,
                expected,
                actual: frame.atoms.len(),
            });
        }
    }
    Ok(())
}

fn min_value(values: &[f64]) -> f64 {
    values.iter().copied().fold(f64::INFINITY, f64::min)
}

fn max_value(values: &[f64]) -> f64 {
    values.iter().copied().fold(f64::NEG_INFINITY, f64::max)
}

fn mean_value(values: &[f64]) -> f64 {
    values.iter().sum::<f64>() / values.len() as f64
}

fn default_rmsd_alignment() -> String {
    "centered-kabsch".to_string()
}

fn aligned_rmsd(frame: &XyzFrame, reference: &XyzFrame) -> f64 {
    let mobile_centroid = centroid(&frame.atoms);
    let reference_centroid = centroid(&reference.atoms);
    let mobile_centered = centered_points(&frame.atoms, mobile_centroid);
    let reference_centered = centered_points(&reference.atoms, reference_centroid);

    let quaternion = best_fit_quaternion(&mobile_centered, &reference_centered);
    let mut sum_sq = 0.0;
    for (mobile, reference) in mobile_centered.iter().zip(&reference_centered) {
        let rotated = rotate_vector(*mobile, quaternion);
        let dx = rotated[0] - reference[0];
        let dy = rotated[1] - reference[1];
        let dz = rotated[2] - reference[2];
        sum_sq += dx * dx + dy * dy + dz * dz;
    }

    (sum_sq / frame.atoms.len() as f64).sqrt()
}

fn centroid(atoms: &[XyzAtom]) -> [f64; 3] {
    let mut center = [0.0; 3];
    for atom in atoms {
        center[0] += atom.x;
        center[1] += atom.y;
        center[2] += atom.z;
    }
    let count = atoms.len() as f64;
    [center[0] / count, center[1] / count, center[2] / count]
}

fn centered_points(atoms: &[XyzAtom], centroid: [f64; 3]) -> Vec<[f64; 3]> {
    atoms
        .iter()
        .map(|atom| {
            [
                atom.x - centroid[0],
                atom.y - centroid[1],
                atom.z - centroid[2],
            ]
        })
        .collect()
}

fn best_fit_quaternion(mobile: &[[f64; 3]], reference: &[[f64; 3]]) -> [f64; 4] {
    let mut covariance = [[0.0; 3]; 3];
    for (mobile, reference) in mobile.iter().zip(reference) {
        for row in 0..3 {
            for col in 0..3 {
                covariance[row][col] += mobile[row] * reference[col];
            }
        }
    }

    let sxx = covariance[0][0];
    let sxy = covariance[0][1];
    let sxz = covariance[0][2];
    let syx = covariance[1][0];
    let syy = covariance[1][1];
    let syz = covariance[1][2];
    let szx = covariance[2][0];
    let szy = covariance[2][1];
    let szz = covariance[2][2];
    let trace = sxx + syy + szz;
    let horn = [
        [trace, syz - szy, szx - sxz, sxy - syx],
        [syz - szy, sxx - syy - szz, sxy + syx, szx + sxz],
        [szx - sxz, sxy + syx, -sxx + syy - szz, syz + szy],
        [sxy - syx, szx + sxz, syz + szy, -sxx - syy + szz],
    ];

    dominant_eigenvector_4x4(horn)
}

fn dominant_eigenvector_4x4(matrix: [[f64; 4]; 4]) -> [f64; 4] {
    let shift = matrix
        .iter()
        .flat_map(|row| row.iter())
        .map(|value| value.abs())
        .sum::<f64>();
    let mut shifted = matrix;
    for (index, row) in shifted.iter_mut().enumerate() {
        row[index] += shift;
    }

    let initial_norm = quaternion_norm([1.0, 0.5, 0.25, 0.125]);
    let mut vector = [
        1.0 / initial_norm,
        0.5 / initial_norm,
        0.25 / initial_norm,
        0.125 / initial_norm,
    ];
    for _ in 0..500 {
        let next = [
            shifted[0][0] * vector[0]
                + shifted[0][1] * vector[1]
                + shifted[0][2] * vector[2]
                + shifted[0][3] * vector[3],
            shifted[1][0] * vector[0]
                + shifted[1][1] * vector[1]
                + shifted[1][2] * vector[2]
                + shifted[1][3] * vector[3],
            shifted[2][0] * vector[0]
                + shifted[2][1] * vector[1]
                + shifted[2][2] * vector[2]
                + shifted[2][3] * vector[3],
            shifted[3][0] * vector[0]
                + shifted[3][1] * vector[1]
                + shifted[3][2] * vector[2]
                + shifted[3][3] * vector[3],
        ];
        let norm = quaternion_norm(next);
        if norm <= 1.0e-14 || !norm.is_finite() {
            return [1.0, 0.0, 0.0, 0.0];
        }
        vector = [
            next[0] / norm,
            next[1] / norm,
            next[2] / norm,
            next[3] / norm,
        ];
    }
    vector
}

fn quaternion_norm(quaternion: [f64; 4]) -> f64 {
    quaternion
        .iter()
        .map(|value| value * value)
        .sum::<f64>()
        .sqrt()
}

fn rotate_vector(vector: [f64; 3], quaternion: [f64; 4]) -> [f64; 3] {
    let [w, x, y, z] = quaternion;
    let uv = cross([x, y, z], vector);
    let uuv = cross([x, y, z], uv);
    [
        vector[0] + 2.0 * (w * uv[0] + uuv[0]),
        vector[1] + 2.0 * (w * uv[1] + uuv[1]),
        vector[2] + 2.0 * (w * uv[2] + uuv[2]),
    ]
}

fn cross(left: [f64; 3], right: [f64; 3]) -> [f64; 3] {
    [
        left[1] * right[2] - left[2] * right[1],
        left[2] * right[0] - left[0] * right[2],
        left[0] * right[1] - left[1] * right[0],
    ]
}

#[cfg(test)]
mod tests {
    use super::*;
    use md_io::{XyzAtom, XyzFrame};

    #[test]
    fn summarizes_energy_drift_and_temperature() {
        let samples = vec![
            EnergySample::new(0, 0.0, 1.0, -2.0, 0.5),
            EnergySample::new(10, 1.0, 1.5, -2.4, 0.7),
            EnergySample::new(20, 2.0, 1.4, -2.2, 0.6),
        ];

        let summary = summarize_energy(&samples).unwrap();

        assert_eq!(summary.sample_count, 3);
        assert!((summary.total_drift - 0.2).abs() < 1.0e-12);
        assert!((summary.mean_temperature - 0.6).abs() < 1.0e-12);
    }

    #[test]
    fn computes_aligned_rmsd_against_reference_frame() {
        let frames = vec![
            triangle_frame(0, 0.0, 0.0, false),
            triangle_frame(1, 5.0, -2.0, true),
        ];

        let samples = rmsd_series(&frames, 0).unwrap();
        let summary = summarize_rmsd(&samples, 0).unwrap();

        assert_eq!(samples[0].rmsd, 0.0);
        assert!(
            samples[1].rmsd < 1.0e-10,
            "aligned RMSD was {}",
            samples[1].rmsd
        );
        assert!(summary.final_value < 1.0e-10);
    }

    #[test]
    fn unwrap_periodic_frames_reconstructs_boundary_crossing() {
        let frames = vec![
            XyzFrame {
                step: 0,
                atoms: vec![XyzAtom {
                    element: "Ar".to_string(),
                    x: 9.8,
                    y: 1.0,
                    z: 1.0,
                }],
            },
            XyzFrame {
                step: 1,
                atoms: vec![XyzAtom {
                    element: "Ar".to_string(),
                    x: 0.2,
                    y: 1.1,
                    z: 1.0,
                }],
            },
        ];
        let simulation_box = SimulationBox::new(10.0, 10.0, 10.0).unwrap();

        let unwrapped = unwrap_periodic_frames(&frames, simulation_box).unwrap();

        assert!((unwrapped[1].atoms[0].x - 10.2).abs() < 1.0e-12);
        assert!((unwrapped[1].atoms[0].y - 1.1).abs() < 1.0e-12);
    }

    #[test]
    fn computes_pair_distance_series() {
        let frames = vec![frame(0, 0.0), frame(1, 1.0)];

        let samples = pair_distance_series(&frames, 0, 1).unwrap();
        let summary = summarize_pair_distances(&samples).unwrap();

        assert!((samples[0].distance - 1.0).abs() < 1.0e-12);
        assert!((samples[1].distance - 1.0).abs() < 1.0e-12);
        assert_eq!(summary.sample_count, 2);
    }

    fn frame(step: usize, offset: f64) -> XyzFrame {
        XyzFrame {
            step,
            atoms: vec![
                XyzAtom {
                    element: "Ar".to_string(),
                    x: offset,
                    y: 0.0,
                    z: 0.0,
                },
                XyzAtom {
                    element: "Ar".to_string(),
                    x: offset + 1.0,
                    y: 0.0,
                    z: 0.0,
                },
            ],
        }
    }

    fn triangle_frame(step: usize, tx: f64, ty: f64, rotate: bool) -> XyzFrame {
        let points = [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]];
        let atoms = points
            .into_iter()
            .map(|[x, y]| {
                let (x, y) = if rotate { (-y, x) } else { (x, y) };
                XyzAtom {
                    element: "Ar".to_string(),
                    x: x + tx,
                    y: y + ty,
                    z: 0.0,
                }
            })
            .collect();
        XyzFrame { step, atoms }
    }
}
