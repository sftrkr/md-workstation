//! Force calculations for MD Workstation.

use md_core::{CoreError, SimulationBox, SystemState};
use md_neighbor::{NeighborBoundary, NeighborError, NeighborList};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use thiserror::Error;

const PARALLEL_MIN_PAIR_CHUNK_SIZE: usize = 256;
const PARALLEL_TARGET_CHUNKS_PER_THREAD: usize = 4;
const MIN_ANGLE_SIN: f64 = 1.0e-8;
const MIN_DIHEDRAL_NORM: f64 = 1.0e-8;
const DIHEDRAL_FORCE_STEP: f64 = 1.0e-6;

/// Lennard-Jones non-bonded parameters in reduced units.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct LennardJonesParams {
    pub epsilon: f64,
    pub sigma: f64,
    pub cutoff: f64,
}

impl LennardJonesParams {
    pub fn validate(&self) -> Result<(), ForceError> {
        if !self.epsilon.is_finite() || self.epsilon <= 0.0 {
            return Err(ForceError::InvalidParameter {
                name: "epsilon",
                value: self.epsilon,
            });
        }
        if !self.sigma.is_finite() || self.sigma <= 0.0 {
            return Err(ForceError::InvalidParameter {
                name: "sigma",
                value: self.sigma,
            });
        }
        if !self.cutoff.is_finite() || self.cutoff <= 0.0 {
            return Err(ForceError::InvalidParameter {
                name: "cutoff",
                value: self.cutoff,
            });
        }
        Ok(())
    }
}

/// Per-particle Lennard-Jones parameters used by mixed-type calculations.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LennardJonesParticleParams {
    pub epsilon: f64,
    pub sigma: f64,
}

impl LennardJonesParticleParams {
    pub fn validate(&self) -> Result<(), ForceError> {
        if !self.epsilon.is_finite() || self.epsilon <= 0.0 {
            return Err(ForceError::InvalidParameter {
                name: "particle epsilon",
                value: self.epsilon,
            });
        }
        if !self.sigma.is_finite() || self.sigma <= 0.0 {
            return Err(ForceError::InvalidParameter {
                name: "particle sigma",
                value: self.sigma,
            });
        }
        Ok(())
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LennardJonesMixingRule {
    LorentzBerthelot,
}

/// Lennard-Jones options for systems with per-particle LJ parameters.
#[derive(Debug, Clone, PartialEq)]
pub struct MixedLennardJonesOptions {
    pub particle_params: Vec<LennardJonesParticleParams>,
    pub cutoff: f64,
    pub boundary: BoundaryCondition,
    pub shift_potential: bool,
    pub mixing_rule: LennardJonesMixingRule,
}

impl MixedLennardJonesOptions {
    pub fn validate(&self, particle_count: usize) -> Result<(), ForceError> {
        if self.particle_params.len() != particle_count {
            return Err(ForceError::InvalidTopology {
                reason: format!(
                    "mixed LJ parameter count {} does not match particle count {particle_count}",
                    self.particle_params.len()
                ),
            });
        }
        if !self.cutoff.is_finite() || self.cutoff <= 0.0 {
            return Err(ForceError::InvalidParameter {
                name: "cutoff",
                value: self.cutoff,
            });
        }
        for params in &self.particle_params {
            params.validate()?;
        }
        if let BoundaryCondition::Periodic(simulation_box) = self.boundary {
            simulation_box.validate()?;
            let max_cutoff = 0.5 * simulation_box.min_dimension();
            if self.cutoff > max_cutoff {
                return Err(ForceError::CutoffTooLargeForPeriodicBox {
                    cutoff: self.cutoff,
                    max_cutoff,
                });
            }
        }
        Ok(())
    }

    fn pair_params(&self, i: usize, j: usize) -> (f64, f64) {
        let left = self.particle_params[i];
        let right = self.particle_params[j];
        match self.mixing_rule {
            LennardJonesMixingRule::LorentzBerthelot => (
                0.5 * (left.sigma + right.sigma),
                (left.epsilon * right.epsilon).sqrt(),
            ),
        }
    }
}

/// Boundary handling for non-bonded force calculations.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum BoundaryCondition {
    Open,
    Periodic(SimulationBox),
}

/// Lennard-Jones calculation options.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LennardJonesOptions {
    pub params: LennardJonesParams,
    pub boundary: BoundaryCondition,
    pub shift_potential: bool,
}

impl LennardJonesOptions {
    pub fn open(params: LennardJonesParams) -> Self {
        Self {
            params,
            boundary: BoundaryCondition::Open,
            shift_potential: false,
        }
    }

    pub fn periodic(
        params: LennardJonesParams,
        simulation_box: SimulationBox,
        shift_potential: bool,
    ) -> Self {
        Self {
            params,
            boundary: BoundaryCondition::Periodic(simulation_box),
            shift_potential,
        }
    }

    pub fn validate(&self) -> Result<(), ForceError> {
        self.params.validate()?;

        if let BoundaryCondition::Periodic(simulation_box) = self.boundary {
            simulation_box.validate()?;
            let max_cutoff = 0.5 * simulation_box.min_dimension();
            if self.params.cutoff > max_cutoff {
                return Err(ForceError::CutoffTooLargeForPeriodicBox {
                    cutoff: self.params.cutoff,
                    max_cutoff,
                });
            }
        }

        Ok(())
    }
}

/// Harmonic bond term in reduced units.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct HarmonicBond {
    pub i: usize,
    pub j: usize,
    pub k: f64,
    pub r0: f64,
}

impl HarmonicBond {
    pub fn validate(&self, particle_count: usize) -> Result<(), ForceError> {
        if self.i == self.j {
            return Err(ForceError::InvalidTopology {
                reason: "bond atom indices must be different".to_string(),
            });
        }
        if self.i >= particle_count || self.j >= particle_count {
            return Err(ForceError::InvalidTopology {
                reason: format!(
                    "bond indices ({}, {}) are out of range for {particle_count} particles",
                    self.i, self.j
                ),
            });
        }
        if !self.k.is_finite() || self.k <= 0.0 {
            return Err(ForceError::InvalidParameter {
                name: "bond.k",
                value: self.k,
            });
        }
        if !self.r0.is_finite() || self.r0 <= 0.0 {
            return Err(ForceError::InvalidParameter {
                name: "bond.r0",
                value: self.r0,
            });
        }
        Ok(())
    }
}

/// Harmonic angle term in reduced units.
///
/// The angle is formed by atoms i-j-k with atom j at the vertex. theta0 is in
/// radians and force_constant has reduced energy/radian^2 units.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct HarmonicAngle {
    pub i: usize,
    pub j: usize,
    pub k: usize,
    pub force_constant: f64,
    pub theta0: f64,
}

impl HarmonicAngle {
    pub fn validate(&self, particle_count: usize) -> Result<(), ForceError> {
        if self.i == self.j || self.i == self.k || self.j == self.k {
            return Err(ForceError::InvalidTopology {
                reason: "angle atom indices must be different".to_string(),
            });
        }
        if self.i >= particle_count || self.j >= particle_count || self.k >= particle_count {
            return Err(ForceError::InvalidTopology {
                reason: format!(
                    "angle indices ({}, {}, {}) are out of range for {particle_count} particles",
                    self.i, self.j, self.k
                ),
            });
        }
        if !self.force_constant.is_finite() || self.force_constant <= 0.0 {
            return Err(ForceError::InvalidParameter {
                name: "angle.force_constant",
                value: self.force_constant,
            });
        }
        if !self.theta0.is_finite() || self.theta0 <= 0.0 || self.theta0 >= std::f64::consts::PI {
            return Err(ForceError::InvalidParameter {
                name: "angle.theta0",
                value: self.theta0,
            });
        }
        Ok(())
    }
}

/// Periodic dihedral term in reduced units.
///
/// The dihedral is formed by atoms i-j-k-l. phase is in radians and the
/// potential is force_constant * (1 + cos(multiplicity * phi - phase)).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct PeriodicDihedral {
    pub i: usize,
    pub j: usize,
    pub k: usize,
    pub l: usize,
    pub force_constant: f64,
    pub multiplicity: u32,
    pub phase: f64,
}

impl PeriodicDihedral {
    pub fn validate(&self, particle_count: usize) -> Result<(), ForceError> {
        if self.i == self.j
            || self.i == self.k
            || self.i == self.l
            || self.j == self.k
            || self.j == self.l
            || self.k == self.l
        {
            return Err(ForceError::InvalidTopology {
                reason: "dihedral atom indices must be different".to_string(),
            });
        }
        if self.i >= particle_count
            || self.j >= particle_count
            || self.k >= particle_count
            || self.l >= particle_count
        {
            return Err(ForceError::InvalidTopology {
                reason: format!(
                    "dihedral indices ({}, {}, {}, {}) are out of range for {particle_count} particles",
                    self.i, self.j, self.k, self.l
                ),
            });
        }
        if !self.force_constant.is_finite() || self.force_constant <= 0.0 {
            return Err(ForceError::InvalidParameter {
                name: "dihedral.force_constant",
                value: self.force_constant,
            });
        }
        if self.multiplicity == 0 {
            return Err(ForceError::InvalidTopology {
                reason: "dihedral multiplicity must be greater than zero".to_string(),
            });
        }
        if !self.phase.is_finite() {
            return Err(ForceError::InvalidParameter {
                name: "dihedral.phase",
                value: self.phase,
            });
        }
        Ok(())
    }
}

/// Non-bonded pair exclusion.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub struct ExcludedPair {
    pub i: usize,
    pub j: usize,
}

impl ExcludedPair {
    pub fn normalized(&self, particle_count: usize) -> Result<(usize, usize), ForceError> {
        normalize_exclusion_tuple((self.i, self.j), particle_count)
    }
}

/// Validate, sort, and deduplicate non-bonded exclusions.
pub fn normalize_excluded_pairs(
    exclusions: &[ExcludedPair],
    particle_count: usize,
) -> Result<Vec<(usize, usize)>, ForceError> {
    let pairs = exclusions
        .iter()
        .map(|exclusion| exclusion.normalized(particle_count))
        .collect::<Result<Vec<_>, _>>()?;
    normalize_exclusion_tuples(&pairs, particle_count)
}

/// Coulomb non-bonded calculation options in reduced units.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CoulombOptions {
    pub constant: f64,
    pub cutoff: f64,
    pub boundary: BoundaryCondition,
    pub shift_potential: bool,
}

impl CoulombOptions {
    pub fn open(constant: f64, cutoff: f64, shift_potential: bool) -> Self {
        Self {
            constant,
            cutoff,
            boundary: BoundaryCondition::Open,
            shift_potential,
        }
    }

    pub fn periodic(
        constant: f64,
        cutoff: f64,
        simulation_box: SimulationBox,
        shift_potential: bool,
    ) -> Self {
        Self {
            constant,
            cutoff,
            boundary: BoundaryCondition::Periodic(simulation_box),
            shift_potential,
        }
    }

    pub fn validate(&self) -> Result<(), ForceError> {
        if !self.constant.is_finite() {
            return Err(ForceError::InvalidParameter {
                name: "coulomb.constant",
                value: self.constant,
            });
        }
        if !self.cutoff.is_finite() || self.cutoff <= 0.0 {
            return Err(ForceError::InvalidParameter {
                name: "coulomb.cutoff",
                value: self.cutoff,
            });
        }
        if let BoundaryCondition::Periodic(simulation_box) = self.boundary {
            simulation_box.validate()?;
            let max_cutoff = 0.5 * simulation_box.min_dimension();
            if self.cutoff > max_cutoff {
                return Err(ForceError::CutoffTooLargeForPeriodicBox {
                    cutoff: self.cutoff,
                    max_cutoff,
                });
            }
        }
        Ok(())
    }
}

/// Summary produced by a force calculation.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ForceReport {
    pub potential_energy: f64,
    pub pair_count: usize,
}

/// Summary produced by a harmonic bond force calculation.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BondForceReport {
    pub potential_energy: f64,
    pub bond_count: usize,
}

/// Summary produced by a harmonic angle force calculation.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AngleForceReport {
    pub potential_energy: f64,
    pub angle_count: usize,
}

/// Summary produced by a periodic dihedral force calculation.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DihedralForceReport {
    pub potential_energy: f64,
    pub dihedral_count: usize,
}

#[derive(Debug, Error)]
pub enum ForceError {
    #[error(transparent)]
    Core(#[from] CoreError),
    #[error("{name} must be positive and finite: {value}")]
    InvalidParameter { name: &'static str, value: f64 },
    #[error("particles {i} and {j} overlap or are too close for a stable LJ calculation")]
    OverlappingParticles { i: usize, j: usize },
    #[error("force calculation produced a non-finite value")]
    NonFiniteForce,
    #[error("periodic cutoff {cutoff} exceeds half the shortest box length: {max_cutoff}")]
    CutoffTooLargeForPeriodicBox { cutoff: f64, max_cutoff: f64 },
    #[error("neighbor list is incompatible with force options: {reason}")]
    IncompatibleNeighborList { reason: &'static str },
    #[error("invalid topology: {reason}")]
    InvalidTopology { reason: String },
    #[error(transparent)]
    Neighbor(#[from] NeighborError),
}

/// Compute unshifted Lennard-Jones forces with a simple O(N^2) pair loop.
pub fn compute_lennard_jones_forces(
    system: &mut SystemState,
    params: &LennardJonesParams,
) -> Result<ForceReport, ForceError> {
    compute_lennard_jones_forces_with_options(system, &LennardJonesOptions::open(*params))
}

/// Compute Lennard-Jones forces with optional periodic boundaries and potential shifting.
pub fn compute_lennard_jones_forces_with_options(
    system: &mut SystemState,
    options: &LennardJonesOptions,
) -> Result<ForceReport, ForceError> {
    compute_lennard_jones_forces_with_options_and_exclusions(system, options, &[])
}

/// Compute Lennard-Jones forces with non-bonded pair exclusions.
pub fn compute_lennard_jones_forces_with_options_and_exclusions(
    system: &mut SystemState,
    options: &LennardJonesOptions,
    excluded_pairs: &[(usize, usize)],
) -> Result<ForceReport, ForceError> {
    options.validate()?;
    system.validate()?;
    let excluded_pairs = normalize_exclusion_tuples(excluded_pairs, system.particle_count())?;
    system.clear_forces();

    let params = options.params;
    let cutoff_sq = params.cutoff * params.cutoff;
    let sigma_sq = params.sigma * params.sigma;
    let potential_shift = if options.shift_potential {
        lennard_jones_potential_at_distance(params.cutoff, &params)
    } else {
        0.0
    };
    let mut potential_energy = 0.0;
    let mut pair_count = 0;

    for i in 0..system.particle_count() {
        for j in (i + 1)..system.particle_count() {
            if is_excluded_pair(&excluded_pairs, i, j) {
                continue;
            }
            accumulate_lennard_jones_pair(
                system,
                i,
                j,
                options.boundary,
                cutoff_sq,
                sigma_sq,
                params.epsilon,
                potential_shift,
                &mut potential_energy,
                &mut pair_count,
            )?;
        }
    }

    Ok(ForceReport {
        potential_energy,
        pair_count,
    })
}

/// Compute mixed-parameter Lennard-Jones forces with non-bonded pair exclusions.
pub fn compute_lennard_jones_forces_with_mixed_options_and_exclusions(
    system: &mut SystemState,
    options: &MixedLennardJonesOptions,
    excluded_pairs: &[(usize, usize)],
) -> Result<ForceReport, ForceError> {
    options.validate(system.particle_count())?;
    system.validate()?;
    let excluded_pairs = normalize_exclusion_tuples(excluded_pairs, system.particle_count())?;
    system.clear_forces();

    let cutoff_sq = options.cutoff * options.cutoff;
    let mut potential_energy = 0.0;
    let mut pair_count = 0;

    for i in 0..system.particle_count() {
        for j in (i + 1)..system.particle_count() {
            if is_excluded_pair(&excluded_pairs, i, j) {
                continue;
            }
            let (sigma, epsilon) = options.pair_params(i, j);
            let sigma_sq = sigma * sigma;
            let potential_shift = if options.shift_potential {
                lennard_jones_potential(options.cutoff, epsilon, sigma)
            } else {
                0.0
            };
            accumulate_lennard_jones_pair(
                system,
                i,
                j,
                options.boundary,
                cutoff_sq,
                sigma_sq,
                epsilon,
                potential_shift,
                &mut potential_energy,
                &mut pair_count,
            )?;
        }
    }

    Ok(ForceReport {
        potential_energy,
        pair_count,
    })
}

/// Compute Lennard-Jones forces using a prebuilt neighbor list.
pub fn compute_lennard_jones_forces_with_neighbor_list(
    system: &mut SystemState,
    options: &LennardJonesOptions,
    neighbor_list: &NeighborList,
) -> Result<ForceReport, ForceError> {
    compute_lennard_jones_forces_with_neighbor_list_and_exclusions(
        system,
        options,
        neighbor_list,
        &[],
    )
}

/// Compute Lennard-Jones forces using a prebuilt neighbor list and exclusions.
pub fn compute_lennard_jones_forces_with_neighbor_list_and_exclusions(
    system: &mut SystemState,
    options: &LennardJonesOptions,
    neighbor_list: &NeighborList,
    excluded_pairs: &[(usize, usize)],
) -> Result<ForceReport, ForceError> {
    options.validate()?;
    validate_neighbor_list_compatibility(system, options, neighbor_list)?;
    system.validate()?;
    let excluded_pairs = normalize_exclusion_tuples(excluded_pairs, system.particle_count())?;
    system.clear_forces();

    let params = options.params;
    let cutoff_sq = params.cutoff * params.cutoff;
    let sigma_sq = params.sigma * params.sigma;
    let potential_shift = if options.shift_potential {
        lennard_jones_potential_at_distance(params.cutoff, &params)
    } else {
        0.0
    };
    let mut potential_energy = 0.0;
    let mut pair_count = 0;

    for &(i, j) in neighbor_list.pairs() {
        if is_excluded_pair(&excluded_pairs, i, j) {
            continue;
        }
        accumulate_lennard_jones_pair(
            system,
            i,
            j,
            options.boundary,
            cutoff_sq,
            sigma_sq,
            params.epsilon,
            potential_shift,
            &mut potential_energy,
            &mut pair_count,
        )?;
    }

    Ok(ForceReport {
        potential_energy,
        pair_count,
    })
}

/// Compute mixed-parameter Lennard-Jones forces using a prebuilt neighbor list.
pub fn compute_lennard_jones_forces_with_mixed_neighbor_list_and_exclusions(
    system: &mut SystemState,
    options: &MixedLennardJonesOptions,
    neighbor_list: &NeighborList,
    excluded_pairs: &[(usize, usize)],
) -> Result<ForceReport, ForceError> {
    options.validate(system.particle_count())?;
    validate_mixed_neighbor_list_compatibility(system, options, neighbor_list)?;
    system.validate()?;
    let excluded_pairs = normalize_exclusion_tuples(excluded_pairs, system.particle_count())?;
    system.clear_forces();

    let cutoff_sq = options.cutoff * options.cutoff;
    let mut potential_energy = 0.0;
    let mut pair_count = 0;

    for &(i, j) in neighbor_list.pairs() {
        if is_excluded_pair(&excluded_pairs, i, j) {
            continue;
        }
        let (sigma, epsilon) = options.pair_params(i, j);
        let sigma_sq = sigma * sigma;
        let potential_shift = if options.shift_potential {
            lennard_jones_potential(options.cutoff, epsilon, sigma)
        } else {
            0.0
        };
        accumulate_lennard_jones_pair(
            system,
            i,
            j,
            options.boundary,
            cutoff_sq,
            sigma_sq,
            epsilon,
            potential_shift,
            &mut potential_energy,
            &mut pair_count,
        )?;
    }

    Ok(ForceReport {
        potential_energy,
        pair_count,
    })
}

/// Compute Lennard-Jones forces using Rayon over all unique pairs.
pub fn compute_lennard_jones_forces_parallel_with_options(
    system: &mut SystemState,
    options: &LennardJonesOptions,
) -> Result<ForceReport, ForceError> {
    compute_lennard_jones_forces_parallel_with_options_and_exclusions(system, options, &[])
}

/// Compute Lennard-Jones forces using Rayon over all unique pairs with exclusions.
pub fn compute_lennard_jones_forces_parallel_with_options_and_exclusions(
    system: &mut SystemState,
    options: &LennardJonesOptions,
    excluded_pairs: &[(usize, usize)],
) -> Result<ForceReport, ForceError> {
    options.validate()?;
    system.validate()?;
    let excluded_pairs = normalize_exclusion_tuples(excluded_pairs, system.particle_count())?;
    compute_lennard_jones_forces_parallel_for_all_pairs(system, options, &excluded_pairs)
}

/// Compute mixed-parameter Lennard-Jones forces using Rayon over all unique pairs.
pub fn compute_lennard_jones_forces_parallel_with_mixed_options_and_exclusions(
    system: &mut SystemState,
    options: &MixedLennardJonesOptions,
    excluded_pairs: &[(usize, usize)],
) -> Result<ForceReport, ForceError> {
    options.validate(system.particle_count())?;
    system.validate()?;
    let excluded_pairs = normalize_exclusion_tuples(excluded_pairs, system.particle_count())?;
    compute_lennard_jones_forces_parallel_for_all_mixed_pairs(system, options, &excluded_pairs)
}

/// Compute Lennard-Jones forces using Rayon over a prebuilt neighbor list.
pub fn compute_lennard_jones_forces_parallel_with_neighbor_list(
    system: &mut SystemState,
    options: &LennardJonesOptions,
    neighbor_list: &NeighborList,
) -> Result<ForceReport, ForceError> {
    compute_lennard_jones_forces_parallel_with_neighbor_list_and_exclusions(
        system,
        options,
        neighbor_list,
        &[],
    )
}

/// Compute Lennard-Jones forces using Rayon over a neighbor list with exclusions.
pub fn compute_lennard_jones_forces_parallel_with_neighbor_list_and_exclusions(
    system: &mut SystemState,
    options: &LennardJonesOptions,
    neighbor_list: &NeighborList,
    excluded_pairs: &[(usize, usize)],
) -> Result<ForceReport, ForceError> {
    options.validate()?;
    validate_neighbor_list_compatibility(system, options, neighbor_list)?;
    system.validate()?;
    let excluded_pairs = normalize_exclusion_tuples(excluded_pairs, system.particle_count())?;
    compute_lennard_jones_forces_parallel_for_pairs(
        system,
        options,
        neighbor_list.pairs(),
        &excluded_pairs,
    )
}

/// Compute mixed-parameter Lennard-Jones forces using Rayon over a neighbor list.
pub fn compute_lennard_jones_forces_parallel_with_mixed_neighbor_list_and_exclusions(
    system: &mut SystemState,
    options: &MixedLennardJonesOptions,
    neighbor_list: &NeighborList,
    excluded_pairs: &[(usize, usize)],
) -> Result<ForceReport, ForceError> {
    options.validate(system.particle_count())?;
    validate_mixed_neighbor_list_compatibility(system, options, neighbor_list)?;
    system.validate()?;
    let excluded_pairs = normalize_exclusion_tuples(excluded_pairs, system.particle_count())?;
    compute_lennard_jones_forces_parallel_for_mixed_pairs(
        system,
        options,
        neighbor_list.pairs(),
        &excluded_pairs,
    )
}

/// Compute harmonic bond forces after clearing existing forces.
pub fn compute_harmonic_bond_forces(
    system: &mut SystemState,
    bonds: &[HarmonicBond],
    boundary: BoundaryCondition,
) -> Result<BondForceReport, ForceError> {
    system.validate()?;
    system.clear_forces();
    add_harmonic_bond_forces(system, bonds, boundary)
}

/// Add harmonic bond forces to the existing force buffers.
pub fn add_harmonic_bond_forces(
    system: &mut SystemState,
    bonds: &[HarmonicBond],
    boundary: BoundaryCondition,
) -> Result<BondForceReport, ForceError> {
    system.validate()?;
    let mut potential_energy = 0.0;

    for bond in bonds {
        bond.validate(system.particle_count())?;
        let contribution = harmonic_bond_contribution(system, bond, boundary)?;
        system.fx[bond.i] += contribution.fx;
        system.fy[bond.i] += contribution.fy;
        system.fz[bond.i] += contribution.fz;
        system.fx[bond.j] -= contribution.fx;
        system.fy[bond.j] -= contribution.fy;
        system.fz[bond.j] -= contribution.fz;
        potential_energy += contribution.potential;
    }

    Ok(BondForceReport {
        potential_energy,
        bond_count: bonds.len(),
    })
}

/// Compute harmonic angle forces after clearing existing forces.
pub fn compute_harmonic_angle_forces(
    system: &mut SystemState,
    angles: &[HarmonicAngle],
    boundary: BoundaryCondition,
) -> Result<AngleForceReport, ForceError> {
    system.validate()?;
    system.clear_forces();
    add_harmonic_angle_forces(system, angles, boundary)
}

/// Add harmonic angle forces to the existing force buffers.
pub fn add_harmonic_angle_forces(
    system: &mut SystemState,
    angles: &[HarmonicAngle],
    boundary: BoundaryCondition,
) -> Result<AngleForceReport, ForceError> {
    system.validate()?;
    let mut potential_energy = 0.0;

    for angle in angles {
        angle.validate(system.particle_count())?;
        let contribution = harmonic_angle_contribution(system, angle, boundary)?;
        system.fx[angle.i] += contribution.fix;
        system.fy[angle.i] += contribution.fiy;
        system.fz[angle.i] += contribution.fiz;
        system.fx[angle.j] += contribution.fjx;
        system.fy[angle.j] += contribution.fjy;
        system.fz[angle.j] += contribution.fjz;
        system.fx[angle.k] += contribution.fkx;
        system.fy[angle.k] += contribution.fky;
        system.fz[angle.k] += contribution.fkz;
        potential_energy += contribution.potential;
    }

    Ok(AngleForceReport {
        potential_energy,
        angle_count: angles.len(),
    })
}

/// Compute periodic dihedral forces after clearing existing forces.
pub fn compute_periodic_dihedral_forces(
    system: &mut SystemState,
    dihedrals: &[PeriodicDihedral],
    boundary: BoundaryCondition,
) -> Result<DihedralForceReport, ForceError> {
    system.validate()?;
    system.clear_forces();
    add_periodic_dihedral_forces(system, dihedrals, boundary)
}

/// Add periodic dihedral forces to the existing force buffers.
pub fn add_periodic_dihedral_forces(
    system: &mut SystemState,
    dihedrals: &[PeriodicDihedral],
    boundary: BoundaryCondition,
) -> Result<DihedralForceReport, ForceError> {
    system.validate()?;
    let mut potential_energy = 0.0;

    for dihedral in dihedrals {
        dihedral.validate(system.particle_count())?;
        let contribution = periodic_dihedral_contribution(system, dihedral, boundary)?;
        for force in contribution.forces {
            system.fx[force.index] += force.fx;
            system.fy[force.index] += force.fy;
            system.fz[force.index] += force.fz;
        }
        potential_energy += contribution.potential;
    }

    Ok(DihedralForceReport {
        potential_energy,
        dihedral_count: dihedrals.len(),
    })
}

/// Compute Coulomb forces after clearing existing forces.
pub fn compute_coulomb_forces_with_options(
    system: &mut SystemState,
    options: &CoulombOptions,
) -> Result<ForceReport, ForceError> {
    compute_coulomb_forces_with_options_and_exclusions(system, options, &[])
}

/// Compute Coulomb forces with non-bonded pair exclusions after clearing forces.
pub fn compute_coulomb_forces_with_options_and_exclusions(
    system: &mut SystemState,
    options: &CoulombOptions,
    excluded_pairs: &[(usize, usize)],
) -> Result<ForceReport, ForceError> {
    system.validate()?;
    system.clear_forces();
    add_coulomb_forces_with_options_and_exclusions(system, options, excluded_pairs)
}

/// Add Coulomb forces to the existing force buffers.
pub fn add_coulomb_forces_with_options(
    system: &mut SystemState,
    options: &CoulombOptions,
) -> Result<ForceReport, ForceError> {
    add_coulomb_forces_with_options_and_exclusions(system, options, &[])
}

/// Add Coulomb forces to existing buffers while skipping excluded pairs.
pub fn add_coulomb_forces_with_options_and_exclusions(
    system: &mut SystemState,
    options: &CoulombOptions,
    excluded_pairs: &[(usize, usize)],
) -> Result<ForceReport, ForceError> {
    options.validate()?;
    system.validate()?;
    let excluded_pairs = normalize_exclusion_tuples(excluded_pairs, system.particle_count())?;

    let cutoff_sq = options.cutoff * options.cutoff;
    let mut potential_energy = 0.0;
    let mut pair_count = 0;

    for i in 0..system.particle_count() {
        for j in (i + 1)..system.particle_count() {
            if is_excluded_pair(&excluded_pairs, i, j) {
                continue;
            }
            let Some(contribution) = coulomb_pair_contribution(system, i, j, options, cutoff_sq)?
            else {
                continue;
            };
            system.fx[i] += contribution.fx;
            system.fy[i] += contribution.fy;
            system.fz[i] += contribution.fz;
            system.fx[j] -= contribution.fx;
            system.fy[j] -= contribution.fy;
            system.fz[j] -= contribution.fz;
            potential_energy += contribution.potential;
            pair_count += 1;
        }
    }

    Ok(ForceReport {
        potential_energy,
        pair_count,
    })
}

fn compute_lennard_jones_forces_parallel_for_pairs(
    system: &mut SystemState,
    options: &LennardJonesOptions,
    pairs: &[(usize, usize)],
    excluded_pairs: &[(usize, usize)],
) -> Result<ForceReport, ForceError> {
    let params = options.params;
    let cutoff_sq = params.cutoff * params.cutoff;
    let sigma_sq = params.sigma * params.sigma;
    let potential_shift = if options.shift_potential {
        lennard_jones_potential_at_distance(params.cutoff, &params)
    } else {
        0.0
    };
    let particle_count = system.particle_count();
    let chunk_size = parallel_work_chunk_size(pairs.len());

    let local_buffers: Result<Vec<_>, ForceError> = pairs
        .par_chunks(chunk_size)
        .map(|chunk| {
            let mut local = LocalForceBuffer::new(particle_count);
            for &(i, j) in chunk {
                if is_excluded_pair(excluded_pairs, i, j) {
                    continue;
                }
                if let Some(contribution) = lennard_jones_pair_contribution(
                    system,
                    i,
                    j,
                    options.boundary,
                    cutoff_sq,
                    sigma_sq,
                    params.epsilon,
                    potential_shift,
                )? {
                    local.add(contribution);
                }
            }
            Ok(local)
        })
        .collect();

    let total = merge_local_buffers(local_buffers?, particle_count);

    Ok(finish_local_force_buffer(system, total))
}

fn compute_lennard_jones_forces_parallel_for_all_pairs(
    system: &mut SystemState,
    options: &LennardJonesOptions,
    excluded_pairs: &[(usize, usize)],
) -> Result<ForceReport, ForceError> {
    let params = options.params;
    let cutoff_sq = params.cutoff * params.cutoff;
    let sigma_sq = params.sigma * params.sigma;
    let potential_shift = if options.shift_potential {
        lennard_jones_potential_at_distance(params.cutoff, &params)
    } else {
        0.0
    };
    let particle_count = system.particle_count();
    let chunks = parallel_index_chunks(unique_pair_count(particle_count));

    let local_buffers: Result<Vec<_>, ForceError> = chunks
        .par_iter()
        .map(|&(start, end)| {
            let mut local = LocalForceBuffer::new(particle_count);
            for_each_linear_pair(particle_count, start, end, |i, j| {
                if is_excluded_pair(excluded_pairs, i, j) {
                    return Ok(());
                }
                if let Some(contribution) = lennard_jones_pair_contribution(
                    system,
                    i,
                    j,
                    options.boundary,
                    cutoff_sq,
                    sigma_sq,
                    params.epsilon,
                    potential_shift,
                )? {
                    local.add(contribution);
                }
                Ok(())
            })?;
            Ok(local)
        })
        .collect();

    let total = merge_local_buffers(local_buffers?, particle_count);

    Ok(finish_local_force_buffer(system, total))
}

fn compute_lennard_jones_forces_parallel_for_mixed_pairs(
    system: &mut SystemState,
    options: &MixedLennardJonesOptions,
    pairs: &[(usize, usize)],
    excluded_pairs: &[(usize, usize)],
) -> Result<ForceReport, ForceError> {
    let cutoff_sq = options.cutoff * options.cutoff;
    let particle_count = system.particle_count();
    let chunk_size = parallel_work_chunk_size(pairs.len());

    let local_buffers: Result<Vec<_>, ForceError> = pairs
        .par_chunks(chunk_size)
        .map(|chunk| {
            let mut local = LocalForceBuffer::new(particle_count);
            for &(i, j) in chunk {
                if is_excluded_pair(excluded_pairs, i, j) {
                    continue;
                }
                let (sigma, epsilon) = options.pair_params(i, j);
                let sigma_sq = sigma * sigma;
                let potential_shift = if options.shift_potential {
                    lennard_jones_potential(options.cutoff, epsilon, sigma)
                } else {
                    0.0
                };
                if let Some(contribution) = lennard_jones_pair_contribution(
                    system,
                    i,
                    j,
                    options.boundary,
                    cutoff_sq,
                    sigma_sq,
                    epsilon,
                    potential_shift,
                )? {
                    local.add(contribution);
                }
            }
            Ok(local)
        })
        .collect();

    let total = merge_local_buffers(local_buffers?, particle_count);

    Ok(finish_local_force_buffer(system, total))
}

fn compute_lennard_jones_forces_parallel_for_all_mixed_pairs(
    system: &mut SystemState,
    options: &MixedLennardJonesOptions,
    excluded_pairs: &[(usize, usize)],
) -> Result<ForceReport, ForceError> {
    let cutoff_sq = options.cutoff * options.cutoff;
    let particle_count = system.particle_count();
    let chunks = parallel_index_chunks(unique_pair_count(particle_count));

    let local_buffers: Result<Vec<_>, ForceError> = chunks
        .par_iter()
        .map(|&(start, end)| {
            let mut local = LocalForceBuffer::new(particle_count);
            for_each_linear_pair(particle_count, start, end, |i, j| {
                if is_excluded_pair(excluded_pairs, i, j) {
                    return Ok(());
                }
                let (sigma, epsilon) = options.pair_params(i, j);
                let sigma_sq = sigma * sigma;
                let potential_shift = if options.shift_potential {
                    lennard_jones_potential(options.cutoff, epsilon, sigma)
                } else {
                    0.0
                };
                if let Some(contribution) = lennard_jones_pair_contribution(
                    system,
                    i,
                    j,
                    options.boundary,
                    cutoff_sq,
                    sigma_sq,
                    epsilon,
                    potential_shift,
                )? {
                    local.add(contribution);
                }
                Ok(())
            })?;
            Ok(local)
        })
        .collect();

    let total = merge_local_buffers(local_buffers?, particle_count);

    Ok(finish_local_force_buffer(system, total))
}

fn accumulate_lennard_jones_pair(
    system: &mut SystemState,
    i: usize,
    j: usize,
    boundary: BoundaryCondition,
    cutoff_sq: f64,
    sigma_sq: f64,
    epsilon: f64,
    potential_shift: f64,
    potential_energy: &mut f64,
    pair_count: &mut usize,
) -> Result<(), ForceError> {
    let Some(contribution) = lennard_jones_pair_contribution(
        system,
        i,
        j,
        boundary,
        cutoff_sq,
        sigma_sq,
        epsilon,
        potential_shift,
    )?
    else {
        return Ok(());
    };

    system.fx[i] += contribution.fx;
    system.fy[i] += contribution.fy;
    system.fz[i] += contribution.fz;
    system.fx[j] -= contribution.fx;
    system.fy[j] -= contribution.fy;
    system.fz[j] -= contribution.fz;
    *potential_energy += contribution.potential;
    *pair_count += 1;

    Ok(())
}

fn lennard_jones_pair_contribution(
    system: &SystemState,
    i: usize,
    j: usize,
    boundary: BoundaryCondition,
    cutoff_sq: f64,
    sigma_sq: f64,
    epsilon: f64,
    potential_shift: f64,
) -> Result<Option<PairContribution>, ForceError> {
    let (dx, dy, dz) = displacement(system, i, j, boundary);
    let r_sq = dx * dx + dy * dy + dz * dz;

    if r_sq < 1.0e-24 {
        return Err(ForceError::OverlappingParticles { i, j });
    }

    if r_sq > cutoff_sq {
        return Ok(None);
    }

    let inv_r_sq = 1.0 / r_sq;
    let sr2 = sigma_sq * inv_r_sq;
    let sr6 = sr2 * sr2 * sr2;
    let sr12 = sr6 * sr6;
    let potential = 4.0 * epsilon * (sr12 - sr6) - potential_shift;
    let force_scale = 24.0 * epsilon * inv_r_sq * (2.0 * sr12 - sr6);

    let fx = force_scale * dx;
    let fy = force_scale * dy;
    let fz = force_scale * dz;

    if !(potential.is_finite() && fx.is_finite() && fy.is_finite() && fz.is_finite()) {
        return Err(ForceError::NonFiniteForce);
    }

    Ok(Some(PairContribution {
        i,
        j,
        fx,
        fy,
        fz,
        potential,
    }))
}

fn harmonic_bond_contribution(
    system: &SystemState,
    bond: &HarmonicBond,
    boundary: BoundaryCondition,
) -> Result<BondContribution, ForceError> {
    let (dx, dy, dz) = displacement(system, bond.i, bond.j, boundary);
    let r_sq = dx * dx + dy * dy + dz * dz;
    if r_sq < 1.0e-24 {
        return Err(ForceError::OverlappingParticles {
            i: bond.i,
            j: bond.j,
        });
    }

    let r = r_sq.sqrt();
    let displacement_from_equilibrium = r - bond.r0;
    let potential = 0.5 * bond.k * displacement_from_equilibrium * displacement_from_equilibrium;
    let force_scale = -bond.k * displacement_from_equilibrium / r;
    let fx = force_scale * dx;
    let fy = force_scale * dy;
    let fz = force_scale * dz;

    if !(potential.is_finite() && fx.is_finite() && fy.is_finite() && fz.is_finite()) {
        return Err(ForceError::NonFiniteForce);
    }

    Ok(BondContribution {
        fx,
        fy,
        fz,
        potential,
    })
}

fn harmonic_angle_contribution(
    system: &SystemState,
    angle: &HarmonicAngle,
    boundary: BoundaryCondition,
) -> Result<AngleContribution, ForceError> {
    let (a_x, a_y, a_z) = displacement(system, angle.i, angle.j, boundary);
    let (b_x, b_y, b_z) = displacement(system, angle.k, angle.j, boundary);
    let a_sq = a_x * a_x + a_y * a_y + a_z * a_z;
    let b_sq = b_x * b_x + b_y * b_y + b_z * b_z;
    if a_sq < 1.0e-24 {
        return Err(ForceError::OverlappingParticles {
            i: angle.i,
            j: angle.j,
        });
    }
    if b_sq < 1.0e-24 {
        return Err(ForceError::OverlappingParticles {
            i: angle.k,
            j: angle.j,
        });
    }

    let a = a_sq.sqrt();
    let b = b_sq.sqrt();
    let a_hat_x = a_x / a;
    let a_hat_y = a_y / a;
    let a_hat_z = a_z / a;
    let b_hat_x = b_x / b;
    let b_hat_y = b_y / b;
    let b_hat_z = b_z / b;
    let cos_theta = (a_hat_x * b_hat_x + a_hat_y * b_hat_y + a_hat_z * b_hat_z).clamp(-1.0, 1.0);
    let sin_theta_sq = (1.0 - cos_theta * cos_theta).max(0.0);
    let sin_theta = sin_theta_sq.sqrt();
    if sin_theta < MIN_ANGLE_SIN {
        return Err(ForceError::InvalidTopology {
            reason: format!(
                "angle ({}, {}, {}) is too close to 0 or pi radians",
                angle.i, angle.j, angle.k
            ),
        });
    }

    let theta = cos_theta.acos();
    let displacement_from_equilibrium = theta - angle.theta0;
    let potential =
        0.5 * angle.force_constant * displacement_from_equilibrium * displacement_from_equilibrium;
    let d_u_dtheta = angle.force_constant * displacement_from_equilibrium;
    let scale_i = d_u_dtheta / (a * sin_theta);
    let scale_k = d_u_dtheta / (b * sin_theta);

    let fix = scale_i * (b_hat_x - cos_theta * a_hat_x);
    let fiy = scale_i * (b_hat_y - cos_theta * a_hat_y);
    let fiz = scale_i * (b_hat_z - cos_theta * a_hat_z);
    let fkx = scale_k * (a_hat_x - cos_theta * b_hat_x);
    let fky = scale_k * (a_hat_y - cos_theta * b_hat_y);
    let fkz = scale_k * (a_hat_z - cos_theta * b_hat_z);
    let fjx = -(fix + fkx);
    let fjy = -(fiy + fky);
    let fjz = -(fiz + fkz);

    if !(potential.is_finite()
        && fix.is_finite()
        && fiy.is_finite()
        && fiz.is_finite()
        && fjx.is_finite()
        && fjy.is_finite()
        && fjz.is_finite()
        && fkx.is_finite()
        && fky.is_finite()
        && fkz.is_finite())
    {
        return Err(ForceError::NonFiniteForce);
    }

    Ok(AngleContribution {
        fix,
        fiy,
        fiz,
        fjx,
        fjy,
        fjz,
        fkx,
        fky,
        fkz,
        potential,
    })
}

fn periodic_dihedral_contribution(
    system: &SystemState,
    dihedral: &PeriodicDihedral,
    boundary: BoundaryCondition,
) -> Result<DihedralContribution, ForceError> {
    let potential = periodic_dihedral_potential(system, dihedral, boundary)?;
    let atom_indices = [dihedral.i, dihedral.j, dihedral.k, dihedral.l];
    let mut forces = [
        ForceOnAtom::zero(dihedral.i),
        ForceOnAtom::zero(dihedral.j),
        ForceOnAtom::zero(dihedral.k),
        ForceOnAtom::zero(dihedral.l),
    ];

    for (slot, &atom_index) in atom_indices.iter().enumerate() {
        forces[slot].fx =
            periodic_dihedral_force_component(system, dihedral, boundary, atom_index, Axis::X)?;
        forces[slot].fy =
            periodic_dihedral_force_component(system, dihedral, boundary, atom_index, Axis::Y)?;
        forces[slot].fz =
            periodic_dihedral_force_component(system, dihedral, boundary, atom_index, Axis::Z)?;
    }

    if !(potential.is_finite()
        && forces
            .iter()
            .all(|force| force.fx.is_finite() && force.fy.is_finite() && force.fz.is_finite()))
    {
        return Err(ForceError::NonFiniteForce);
    }

    Ok(DihedralContribution { forces, potential })
}

fn periodic_dihedral_force_component(
    system: &SystemState,
    dihedral: &PeriodicDihedral,
    boundary: BoundaryCondition,
    atom_index: usize,
    axis: Axis,
) -> Result<f64, ForceError> {
    let mut plus = system.clone();
    shift_coordinate(&mut plus, atom_index, axis, DIHEDRAL_FORCE_STEP);
    let plus_potential = periodic_dihedral_potential(&plus, dihedral, boundary)?;

    let mut minus = system.clone();
    shift_coordinate(&mut minus, atom_index, axis, -DIHEDRAL_FORCE_STEP);
    let minus_potential = periodic_dihedral_potential(&minus, dihedral, boundary)?;

    Ok(-(plus_potential - minus_potential) / (2.0 * DIHEDRAL_FORCE_STEP))
}

fn shift_coordinate(system: &mut SystemState, atom_index: usize, axis: Axis, delta: f64) {
    match axis {
        Axis::X => system.x[atom_index] += delta,
        Axis::Y => system.y[atom_index] += delta,
        Axis::Z => system.z[atom_index] += delta,
    }
}

fn periodic_dihedral_potential(
    system: &SystemState,
    dihedral: &PeriodicDihedral,
    boundary: BoundaryCondition,
) -> Result<f64, ForceError> {
    let phi = periodic_dihedral_angle(system, dihedral, boundary)?;
    let argument = dihedral.multiplicity as f64 * phi - dihedral.phase;
    let potential = dihedral.force_constant * (1.0 + argument.cos());
    if !potential.is_finite() {
        return Err(ForceError::NonFiniteForce);
    }
    Ok(potential)
}

fn periodic_dihedral_angle(
    system: &SystemState,
    dihedral: &PeriodicDihedral,
    boundary: BoundaryCondition,
) -> Result<f64, ForceError> {
    let b0 = displacement(system, dihedral.i, dihedral.j, boundary);
    let b1 = displacement(system, dihedral.k, dihedral.j, boundary);
    let b2 = displacement(system, dihedral.l, dihedral.k, boundary);

    let b1_norm = vector_norm(b1);
    if b1_norm < 1.0e-24 {
        return Err(ForceError::OverlappingParticles {
            i: dihedral.j,
            j: dihedral.k,
        });
    }

    let b1_hat = scale_vector(b1, 1.0 / b1_norm);
    let v = subtract_vector(b0, scale_vector(b1_hat, dot_vector(b0, b1_hat)));
    let w = subtract_vector(b2, scale_vector(b1_hat, dot_vector(b2, b1_hat)));
    if vector_norm(v) < MIN_DIHEDRAL_NORM || vector_norm(w) < MIN_DIHEDRAL_NORM {
        return Err(ForceError::InvalidTopology {
            reason: format!(
                "dihedral ({}, {}, {}, {}) is degenerate or nearly linear",
                dihedral.i, dihedral.j, dihedral.k, dihedral.l
            ),
        });
    }

    let x = dot_vector(v, w);
    let y = dot_vector(cross_vector(b1_hat, v), w);
    let phi = y.atan2(x);
    if !phi.is_finite() {
        return Err(ForceError::NonFiniteForce);
    }
    Ok(phi)
}

fn coulomb_pair_contribution(
    system: &SystemState,
    i: usize,
    j: usize,
    options: &CoulombOptions,
    cutoff_sq: f64,
) -> Result<Option<PairContribution>, ForceError> {
    let charge_product = system.charge[i] * system.charge[j];
    if charge_product == 0.0 || options.constant == 0.0 {
        return Ok(None);
    }

    let (dx, dy, dz) = displacement(system, i, j, options.boundary);
    let r_sq = dx * dx + dy * dy + dz * dz;
    if r_sq < 1.0e-24 {
        return Err(ForceError::OverlappingParticles { i, j });
    }
    if r_sq > cutoff_sq {
        return Ok(None);
    }

    let r = r_sq.sqrt();
    let unshifted_potential = options.constant * charge_product / r;
    let potential_shift = if options.shift_potential {
        options.constant * charge_product / options.cutoff
    } else {
        0.0
    };
    let potential = unshifted_potential - potential_shift;
    let force_scale = options.constant * charge_product / (r_sq * r);
    let fx = force_scale * dx;
    let fy = force_scale * dy;
    let fz = force_scale * dz;

    if !(potential.is_finite() && fx.is_finite() && fy.is_finite() && fz.is_finite()) {
        return Err(ForceError::NonFiniteForce);
    }

    Ok(Some(PairContribution {
        i,
        j,
        fx,
        fy,
        fz,
        potential,
    }))
}

#[derive(Debug, Clone, Copy)]
struct PairContribution {
    i: usize,
    j: usize,
    fx: f64,
    fy: f64,
    fz: f64,
    potential: f64,
}

#[derive(Debug, Clone, Copy)]
struct BondContribution {
    fx: f64,
    fy: f64,
    fz: f64,
    potential: f64,
}

#[derive(Debug, Clone, Copy)]
struct AngleContribution {
    fix: f64,
    fiy: f64,
    fiz: f64,
    fjx: f64,
    fjy: f64,
    fjz: f64,
    fkx: f64,
    fky: f64,
    fkz: f64,
    potential: f64,
}

#[derive(Debug, Clone, Copy)]
struct DihedralContribution {
    forces: [ForceOnAtom; 4],
    potential: f64,
}

#[derive(Debug, Clone, Copy)]
struct ForceOnAtom {
    index: usize,
    fx: f64,
    fy: f64,
    fz: f64,
}

impl ForceOnAtom {
    fn zero(index: usize) -> Self {
        Self {
            index,
            fx: 0.0,
            fy: 0.0,
            fz: 0.0,
        }
    }
}

#[derive(Debug, Clone, Copy)]
enum Axis {
    X,
    Y,
    Z,
}

#[derive(Debug, Clone)]
struct LocalForceBuffer {
    fx: Vec<f64>,
    fy: Vec<f64>,
    fz: Vec<f64>,
    potential_energy: f64,
    pair_count: usize,
}

impl LocalForceBuffer {
    fn new(particle_count: usize) -> Self {
        Self {
            fx: vec![0.0; particle_count],
            fy: vec![0.0; particle_count],
            fz: vec![0.0; particle_count],
            potential_energy: 0.0,
            pair_count: 0,
        }
    }

    fn add(&mut self, contribution: PairContribution) {
        self.fx[contribution.i] += contribution.fx;
        self.fy[contribution.i] += contribution.fy;
        self.fz[contribution.i] += contribution.fz;
        self.fx[contribution.j] -= contribution.fx;
        self.fy[contribution.j] -= contribution.fy;
        self.fz[contribution.j] -= contribution.fz;
        self.potential_energy += contribution.potential;
        self.pair_count += 1;
    }

    fn merge(&mut self, other: Self) {
        for i in 0..self.fx.len() {
            self.fx[i] += other.fx[i];
            self.fy[i] += other.fy[i];
            self.fz[i] += other.fz[i];
        }
        self.potential_energy += other.potential_energy;
        self.pair_count += other.pair_count;
    }

    fn write_to_system(self, system: &mut SystemState) {
        system.fx = self.fx;
        system.fy = self.fy;
        system.fz = self.fz;
    }
}

fn merge_local_buffers(
    local_buffers: Vec<LocalForceBuffer>,
    particle_count: usize,
) -> LocalForceBuffer {
    let mut total = LocalForceBuffer::new(particle_count);
    for local in local_buffers {
        total.merge(local);
    }
    total
}

fn finish_local_force_buffer(system: &mut SystemState, total: LocalForceBuffer) -> ForceReport {
    let report = ForceReport {
        potential_energy: total.potential_energy,
        pair_count: total.pair_count,
    };
    total.write_to_system(system);
    report
}

fn parallel_work_chunk_size(item_count: usize) -> usize {
    if item_count == 0 {
        return 1;
    }
    let target_chunks = rayon::current_num_threads()
        .max(1)
        .saturating_mul(PARALLEL_TARGET_CHUNKS_PER_THREAD)
        .max(1);
    let balanced = item_count.saturating_add(target_chunks - 1) / target_chunks;
    balanced.max(PARALLEL_MIN_PAIR_CHUNK_SIZE)
}

fn parallel_index_chunks(item_count: usize) -> Vec<(usize, usize)> {
    index_chunks(item_count, parallel_work_chunk_size(item_count))
}

fn index_chunks(item_count: usize, chunk_size: usize) -> Vec<(usize, usize)> {
    let chunk_size = chunk_size.max(1);
    let mut chunks = Vec::new();
    let mut start = 0;
    while start < item_count {
        let end = start.saturating_add(chunk_size).min(item_count);
        chunks.push((start, end));
        start = end;
    }
    chunks
}

fn unique_pair_count(particle_count: usize) -> usize {
    let count =
        (particle_count as u128).saturating_mul(particle_count.saturating_sub(1) as u128) / 2;
    count.min(usize::MAX as u128) as usize
}

fn pair_index_offset(row: usize, particle_count: usize) -> usize {
    let row = row as u128;
    let particle_count = particle_count as u128;
    let offset = row.saturating_mul(
        particle_count
            .saturating_mul(2)
            .saturating_sub(row)
            .saturating_sub(1),
    ) / 2;
    offset.min(usize::MAX as u128) as usize
}

fn pair_indices_from_linear_index(index: usize, particle_count: usize) -> (usize, usize) {
    debug_assert!(particle_count >= 2);
    debug_assert!(index < unique_pair_count(particle_count));

    let mut low = 0;
    let mut high = particle_count.saturating_sub(2);
    while low < high {
        let mid = low + (high - low + 1) / 2;
        if pair_index_offset(mid, particle_count) <= index {
            low = mid;
        } else {
            high = mid - 1;
        }
    }

    let i = low;
    let row_start = pair_index_offset(i, particle_count);
    let j = i + 1 + (index - row_start);
    (i, j)
}

fn advance_pair_indices(i: &mut usize, j: &mut usize, particle_count: usize) {
    if *j + 1 < particle_count {
        *j += 1;
    } else {
        *i += 1;
        *j = *i + 1;
    }
}

fn for_each_linear_pair<F>(
    particle_count: usize,
    start: usize,
    end: usize,
    mut visit: F,
) -> Result<(), ForceError>
where
    F: FnMut(usize, usize) -> Result<(), ForceError>,
{
    if start >= end {
        return Ok(());
    }

    let (mut i, mut j) = pair_indices_from_linear_index(start, particle_count);
    for _ in start..end {
        visit(i, j)?;
        advance_pair_indices(&mut i, &mut j, particle_count);
    }
    Ok(())
}

fn dot_vector(a: (f64, f64, f64), b: (f64, f64, f64)) -> f64 {
    a.0 * b.0 + a.1 * b.1 + a.2 * b.2
}

fn cross_vector(a: (f64, f64, f64), b: (f64, f64, f64)) -> (f64, f64, f64) {
    (
        a.1 * b.2 - a.2 * b.1,
        a.2 * b.0 - a.0 * b.2,
        a.0 * b.1 - a.1 * b.0,
    )
}

fn subtract_vector(a: (f64, f64, f64), b: (f64, f64, f64)) -> (f64, f64, f64) {
    (a.0 - b.0, a.1 - b.1, a.2 - b.2)
}

fn scale_vector(vector: (f64, f64, f64), scale: f64) -> (f64, f64, f64) {
    (vector.0 * scale, vector.1 * scale, vector.2 * scale)
}

fn vector_norm(vector: (f64, f64, f64)) -> f64 {
    dot_vector(vector, vector).sqrt()
}

fn normalize_exclusion_tuples(
    exclusions: &[(usize, usize)],
    particle_count: usize,
) -> Result<Vec<(usize, usize)>, ForceError> {
    let mut pairs = exclusions
        .iter()
        .copied()
        .map(|pair| normalize_exclusion_tuple(pair, particle_count))
        .collect::<Result<Vec<_>, _>>()?;
    pairs.sort_unstable();
    pairs.dedup();
    Ok(pairs)
}

fn normalize_exclusion_tuple(
    pair: (usize, usize),
    particle_count: usize,
) -> Result<(usize, usize), ForceError> {
    let (i, j) = pair;
    if i == j {
        return Err(ForceError::InvalidTopology {
            reason: "excluded pair atom indices must be different".to_string(),
        });
    }
    if i >= particle_count || j >= particle_count {
        return Err(ForceError::InvalidTopology {
            reason: format!(
                "excluded pair indices ({i}, {j}) are out of range for {particle_count} particles"
            ),
        });
    }
    Ok(if i < j { (i, j) } else { (j, i) })
}

fn is_excluded_pair(exclusions: &[(usize, usize)], i: usize, j: usize) -> bool {
    let pair = if i < j { (i, j) } else { (j, i) };
    exclusions.binary_search(&pair).is_ok()
}

fn validate_neighbor_list_compatibility(
    system: &SystemState,
    options: &LennardJonesOptions,
    neighbor_list: &NeighborList,
) -> Result<(), ForceError> {
    if neighbor_list.particle_count() != system.particle_count() {
        return Err(ForceError::IncompatibleNeighborList {
            reason: "particle count changed",
        });
    }

    let neighbor_config = neighbor_list.config();
    if neighbor_config.cutoff < options.params.cutoff {
        return Err(ForceError::IncompatibleNeighborList {
            reason: "neighbor cutoff is smaller than force cutoff",
        });
    }

    if !boundaries_match(options.boundary, neighbor_config.boundary) {
        return Err(ForceError::IncompatibleNeighborList {
            reason: "boundary condition differs",
        });
    }

    Ok(())
}

fn validate_mixed_neighbor_list_compatibility(
    system: &SystemState,
    options: &MixedLennardJonesOptions,
    neighbor_list: &NeighborList,
) -> Result<(), ForceError> {
    if neighbor_list.particle_count() != system.particle_count() {
        return Err(ForceError::IncompatibleNeighborList {
            reason: "particle count changed",
        });
    }

    let neighbor_config = neighbor_list.config();
    if neighbor_config.cutoff < options.cutoff {
        return Err(ForceError::IncompatibleNeighborList {
            reason: "neighbor cutoff is smaller than force cutoff",
        });
    }

    if !boundaries_match(options.boundary, neighbor_config.boundary) {
        return Err(ForceError::IncompatibleNeighborList {
            reason: "boundary condition differs",
        });
    }

    Ok(())
}

fn boundaries_match(
    force_boundary: BoundaryCondition,
    neighbor_boundary: NeighborBoundary,
) -> bool {
    match (force_boundary, neighbor_boundary) {
        (BoundaryCondition::Open, NeighborBoundary::Open) => true,
        (BoundaryCondition::Periodic(a), NeighborBoundary::Periodic(b)) => a == b,
        _ => false,
    }
}

fn displacement(
    system: &SystemState,
    i: usize,
    j: usize,
    boundary: BoundaryCondition,
) -> (f64, f64, f64) {
    let dx = system.x[i] - system.x[j];
    let dy = system.y[i] - system.y[j];
    let dz = system.z[i] - system.z[j];

    match boundary {
        BoundaryCondition::Open => (dx, dy, dz),
        BoundaryCondition::Periodic(simulation_box) => {
            simulation_box.minimum_image_delta(dx, dy, dz)
        }
    }
}

fn lennard_jones_potential_at_distance(distance: f64, params: &LennardJonesParams) -> f64 {
    lennard_jones_potential(distance, params.epsilon, params.sigma)
}

fn lennard_jones_potential(distance: f64, epsilon: f64, sigma: f64) -> f64 {
    let sr = sigma / distance;
    let sr2 = sr * sr;
    let sr6 = sr2 * sr2 * sr2;
    let sr12 = sr6 * sr6;
    4.0 * epsilon * (sr12 - sr6)
}

#[cfg(test)]
mod tests {
    use super::*;
    use md_neighbor::{NeighborBoundary, NeighborList, NeighborListConfig};

    fn two_particle_system(distance: f64) -> SystemState {
        let mut state = SystemState::new(2);
        state.element[0] = "Ar".to_string();
        state.element[1] = "Ar".to_string();
        state.x[0] = 0.0;
        state.x[1] = distance;
        state
    }

    fn right_angle_system() -> SystemState {
        let mut state = SystemState::new(3);
        state.x[0] = 1.0;
        state.y[0] = 0.0;
        state.x[1] = 0.0;
        state.y[1] = 0.0;
        state.x[2] = 0.0;
        state.y[2] = 1.0;
        state
    }

    fn four_atom_dihedral_system() -> SystemState {
        let mut state = SystemState::new(4);
        state.x[0] = 1.8;
        state.y[0] = 3.0;
        state.z[0] = 3.0;
        state.x[1] = 3.0;
        state.y[1] = 3.0;
        state.z[1] = 3.0;
        state.x[2] = 3.0;
        state.y[2] = 4.2;
        state.z[2] = 3.0;
        state.x[3] = 4.0;
        state.y[3] = 4.6;
        state.z[3] = 3.8;
        state
    }

    #[test]
    fn two_particles_at_lj_equilibrium_have_near_zero_force() {
        let equilibrium_distance = 2.0_f64.powf(1.0 / 6.0);
        let mut state = two_particle_system(equilibrium_distance);
        let params = LennardJonesParams {
            epsilon: 1.0,
            sigma: 1.0,
            cutoff: 2.5,
        };

        let report = compute_lennard_jones_forces(&mut state, &params).unwrap();

        assert_eq!(report.pair_count, 1);
        assert!((report.potential_energy + 1.0).abs() < 1.0e-12);
        assert!(state.fx[0].abs() < 1.0e-12);
        assert!(state.fx[1].abs() < 1.0e-12);
    }

    #[test]
    fn pair_forces_are_equal_and_opposite() {
        let mut state = two_particle_system(1.5);
        let params = LennardJonesParams {
            epsilon: 1.0,
            sigma: 1.0,
            cutoff: 2.5,
        };

        compute_lennard_jones_forces(&mut state, &params).unwrap();

        assert!((state.fx[0] + state.fx[1]).abs() < 1.0e-12);
        assert!((state.fy[0] + state.fy[1]).abs() < 1.0e-12);
        assert!((state.fz[0] + state.fz[1]).abs() < 1.0e-12);
    }

    #[test]
    fn particles_outside_cutoff_are_ignored() {
        let mut state = two_particle_system(3.0);
        let params = LennardJonesParams {
            epsilon: 1.0,
            sigma: 1.0,
            cutoff: 2.5,
        };

        let report = compute_lennard_jones_forces(&mut state, &params).unwrap();

        assert_eq!(report.pair_count, 0);
        assert_eq!(report.potential_energy, 0.0);
        assert_eq!(state.fx[0], 0.0);
    }

    #[test]
    fn excluded_lj_pair_is_ignored() {
        let mut state = two_particle_system(1.0);
        let options = LennardJonesOptions::open(LennardJonesParams {
            epsilon: 1.0,
            sigma: 1.0,
            cutoff: 2.5,
        });

        let report = compute_lennard_jones_forces_with_options_and_exclusions(
            &mut state,
            &options,
            &[(0, 1)],
        )
        .unwrap();

        assert_eq!(report.pair_count, 0);
        assert_eq!(report.potential_energy, 0.0);
        assert_eq!(state.fx, vec![0.0, 0.0]);
    }

    #[test]
    fn potential_is_zero_at_sigma() {
        let mut state = two_particle_system(1.0);
        let params = LennardJonesParams {
            epsilon: 1.0,
            sigma: 1.0,
            cutoff: 2.5,
        };

        let report = compute_lennard_jones_forces(&mut state, &params).unwrap();

        assert_eq!(report.pair_count, 1);
        assert!(report.potential_energy.abs() < 1.0e-12);
        assert!((state.fx[0] + 24.0).abs() < 1.0e-12);
        assert!((state.fx[1] - 24.0).abs() < 1.0e-12);
    }

    #[test]
    fn mixed_lj_uses_lorentz_berthelot_pair_parameters() {
        let mixed_sigma = 1.5;
        let mixed_epsilon = 2.0;
        let equilibrium_distance = mixed_sigma * 2.0_f64.powf(1.0 / 6.0);
        let mut state = two_particle_system(equilibrium_distance);
        let options = MixedLennardJonesOptions {
            particle_params: vec![
                LennardJonesParticleParams {
                    epsilon: 1.0,
                    sigma: 1.0,
                },
                LennardJonesParticleParams {
                    epsilon: 4.0,
                    sigma: 2.0,
                },
            ],
            cutoff: 4.0,
            boundary: BoundaryCondition::Open,
            shift_potential: false,
            mixing_rule: LennardJonesMixingRule::LorentzBerthelot,
        };

        let report = compute_lennard_jones_forces_with_mixed_options_and_exclusions(
            &mut state,
            &options,
            &[],
        )
        .unwrap();

        assert_eq!(report.pair_count, 1);
        assert!((report.potential_energy + mixed_epsilon).abs() < 1.0e-12);
        assert!(state.fx[0].abs() < 1.0e-10);
        assert!(state.fx[1].abs() < 1.0e-10);
    }

    #[test]
    fn total_force_is_zero_for_three_particle_system() {
        let mut state = SystemState::new(3);
        state.x[0] = 0.0;
        state.y[0] = 0.0;
        state.x[1] = 1.4;
        state.y[1] = 0.0;
        state.x[2] = 0.7;
        state.y[2] = 1.212435565;

        let params = LennardJonesParams {
            epsilon: 1.0,
            sigma: 1.0,
            cutoff: 2.5,
        };

        compute_lennard_jones_forces(&mut state, &params).unwrap();

        let total_fx: f64 = state.fx.iter().sum();
        let total_fy: f64 = state.fy.iter().sum();
        let total_fz: f64 = state.fz.iter().sum();
        assert!(total_fx.abs() < 1.0e-12);
        assert!(total_fy.abs() < 1.0e-12);
        assert!(total_fz.abs() < 1.0e-12);
    }

    #[test]
    fn overlapping_particles_are_rejected() {
        let mut state = two_particle_system(0.0);
        let params = LennardJonesParams {
            epsilon: 1.0,
            sigma: 1.0,
            cutoff: 2.5,
        };

        let error = compute_lennard_jones_forces(&mut state, &params).unwrap_err();

        assert!(matches!(
            error,
            ForceError::OverlappingParticles { i: 0, j: 1 }
        ));
    }

    #[test]
    fn invalid_parameters_are_rejected() {
        let params = LennardJonesParams {
            epsilon: 0.0,
            sigma: 1.0,
            cutoff: 2.5,
        };

        let error = params.validate().unwrap_err();

        assert!(matches!(
            error,
            ForceError::InvalidParameter {
                name: "epsilon",
                ..
            }
        ));
    }

    #[test]
    fn harmonic_bond_at_equilibrium_has_zero_force_and_energy() {
        let mut state = two_particle_system(1.3);
        let bonds = [HarmonicBond {
            i: 0,
            j: 1,
            k: 25.0,
            r0: 1.3,
        }];

        let report =
            compute_harmonic_bond_forces(&mut state, &bonds, BoundaryCondition::Open).unwrap();

        assert_eq!(report.bond_count, 1);
        assert!(report.potential_energy.abs() < 1.0e-12);
        assert!(state.fx[0].abs() < 1.0e-12);
        assert!(state.fx[1].abs() < 1.0e-12);
    }

    #[test]
    fn stretched_harmonic_bond_pulls_atoms_together() {
        let mut state = two_particle_system(1.5);
        let bonds = [HarmonicBond {
            i: 0,
            j: 1,
            k: 10.0,
            r0: 1.0,
        }];

        let report =
            compute_harmonic_bond_forces(&mut state, &bonds, BoundaryCondition::Open).unwrap();

        assert!((report.potential_energy - 1.25).abs() < 1.0e-12);
        assert!((state.fx[0] - 5.0).abs() < 1.0e-12);
        assert!((state.fx[1] + 5.0).abs() < 1.0e-12);
    }

    #[test]
    fn invalid_harmonic_bond_indices_are_rejected() {
        let mut state = two_particle_system(1.0);
        let bonds = [HarmonicBond {
            i: 0,
            j: 2,
            k: 10.0,
            r0: 1.0,
        }];

        let error =
            compute_harmonic_bond_forces(&mut state, &bonds, BoundaryCondition::Open).unwrap_err();

        assert!(matches!(error, ForceError::InvalidTopology { .. }));
    }

    #[test]
    fn harmonic_angle_at_equilibrium_has_zero_force_and_energy() {
        let mut state = right_angle_system();
        let angles = [HarmonicAngle {
            i: 0,
            j: 1,
            k: 2,
            force_constant: 12.0,
            theta0: std::f64::consts::FRAC_PI_2,
        }];

        let report =
            compute_harmonic_angle_forces(&mut state, &angles, BoundaryCondition::Open).unwrap();

        assert_eq!(report.angle_count, 1);
        assert!(report.potential_energy.abs() < 1.0e-12);
        assert!(state.fx.iter().all(|force| force.abs() < 1.0e-12));
        assert!(state.fy.iter().all(|force| force.abs() < 1.0e-12));
        assert!(state.fz.iter().all(|force| force.abs() < 1.0e-12));
    }

    #[test]
    fn compressed_harmonic_angle_pushes_atoms_apart() {
        let mut state = right_angle_system();
        let angles = [HarmonicAngle {
            i: 0,
            j: 1,
            k: 2,
            force_constant: 10.0,
            theta0: 2.0 * std::f64::consts::PI / 3.0,
        }];

        let report =
            compute_harmonic_angle_forces(&mut state, &angles, BoundaryCondition::Open).unwrap();

        assert!(report.potential_energy > 0.0);
        assert!(state.fy[0] < 0.0);
        assert!(state.fx[2] < 0.0);
        assert!((state.fx[0] + state.fx[1] + state.fx[2]).abs() < 1.0e-12);
        assert!((state.fy[0] + state.fy[1] + state.fy[2]).abs() < 1.0e-12);
    }

    #[test]
    fn periodic_dihedral_potential_matches_periodic_formula() {
        let mut state = four_atom_dihedral_system();
        let dihedrals = [PeriodicDihedral {
            i: 0,
            j: 1,
            k: 2,
            l: 3,
            force_constant: 0.5,
            multiplicity: 3,
            phase: 0.2,
        }];
        let phi = periodic_dihedral_angle(&state, &dihedrals[0], BoundaryCondition::Open).unwrap();
        let expected = dihedrals[0].force_constant
            * (1.0 + (dihedrals[0].multiplicity as f64 * phi - dihedrals[0].phase).cos());

        let report =
            compute_periodic_dihedral_forces(&mut state, &dihedrals, BoundaryCondition::Open)
                .unwrap();

        assert_eq!(report.dihedral_count, 1);
        assert!((report.potential_energy - expected).abs() < 1.0e-12);
        assert!(state.fx.iter().all(|force| force.is_finite()));
        assert!(state.fy.iter().all(|force| force.is_finite()));
        assert!(state.fz.iter().all(|force| force.is_finite()));
        assert!(state.fx.iter().sum::<f64>().abs() < 1.0e-7);
        assert!(state.fy.iter().sum::<f64>().abs() < 1.0e-7);
        assert!(state.fz.iter().sum::<f64>().abs() < 1.0e-7);
    }

    #[test]
    fn periodic_dihedral_at_minimum_has_near_zero_force_and_energy() {
        let mut state = four_atom_dihedral_system();
        let template = PeriodicDihedral {
            i: 0,
            j: 1,
            k: 2,
            l: 3,
            force_constant: 0.5,
            multiplicity: 1,
            phase: 0.0,
        };
        let phi = periodic_dihedral_angle(&state, &template, BoundaryCondition::Open).unwrap();
        let dihedrals = [PeriodicDihedral {
            phase: phi - std::f64::consts::PI,
            ..template
        }];

        let report =
            compute_periodic_dihedral_forces(&mut state, &dihedrals, BoundaryCondition::Open)
                .unwrap();

        assert!(report.potential_energy.abs() < 1.0e-12);
        assert!(state.fx.iter().all(|force| force.abs() < 1.0e-6));
        assert!(state.fy.iter().all(|force| force.abs() < 1.0e-6));
        assert!(state.fz.iter().all(|force| force.abs() < 1.0e-6));
    }

    #[test]
    fn invalid_periodic_dihedral_indices_are_rejected() {
        let mut state = four_atom_dihedral_system();
        let dihedrals = [PeriodicDihedral {
            i: 0,
            j: 1,
            k: 2,
            l: 4,
            force_constant: 0.5,
            multiplicity: 3,
            phase: 0.0,
        }];

        let error =
            compute_periodic_dihedral_forces(&mut state, &dihedrals, BoundaryCondition::Open)
                .unwrap_err();

        assert!(matches!(error, ForceError::InvalidTopology { .. }));
    }

    #[test]
    fn coulomb_like_charges_repel() {
        let mut state = two_particle_system(2.0);
        state.charge[0] = 1.0;
        state.charge[1] = 1.0;
        let options = CoulombOptions::open(1.0, 3.0, false);

        let report = compute_coulomb_forces_with_options(&mut state, &options).unwrap();

        assert_eq!(report.pair_count, 1);
        assert!((report.potential_energy - 0.5).abs() < 1.0e-12);
        assert!((state.fx[0] + 0.25).abs() < 1.0e-12);
        assert!((state.fx[1] - 0.25).abs() < 1.0e-12);
    }

    #[test]
    fn coulomb_opposite_charges_attract() {
        let mut state = two_particle_system(2.0);
        state.charge[0] = 1.0;
        state.charge[1] = -1.0;
        let options = CoulombOptions::open(1.0, 3.0, false);

        let report = compute_coulomb_forces_with_options(&mut state, &options).unwrap();

        assert_eq!(report.pair_count, 1);
        assert!((report.potential_energy + 0.5).abs() < 1.0e-12);
        assert!((state.fx[0] - 0.25).abs() < 1.0e-12);
        assert!((state.fx[1] + 0.25).abs() < 1.0e-12);
    }

    #[test]
    fn coulomb_shifted_potential_is_zero_at_cutoff() {
        let mut state = two_particle_system(2.0);
        state.charge[0] = 1.0;
        state.charge[1] = -1.0;
        let options = CoulombOptions::open(1.0, 2.0, true);

        let report = compute_coulomb_forces_with_options(&mut state, &options).unwrap();

        assert_eq!(report.pair_count, 1);
        assert!(report.potential_energy.abs() < 1.0e-12);
    }

    #[test]
    fn excluded_coulomb_pair_is_ignored() {
        let mut state = two_particle_system(2.0);
        state.charge[0] = 1.0;
        state.charge[1] = -1.0;
        let options = CoulombOptions::open(1.0, 3.0, false);

        let report =
            add_coulomb_forces_with_options_and_exclusions(&mut state, &options, &[(1, 0)])
                .unwrap();

        assert_eq!(report.pair_count, 0);
        assert_eq!(report.potential_energy, 0.0);
        assert_eq!(state.fx, vec![0.0, 0.0]);
    }

    #[test]
    fn periodic_minimum_image_interacts_across_boundary() {
        let mut state = SystemState::new(2);
        state.x[0] = 0.75;
        state.x[1] = 9.75;

        let options = LennardJonesOptions::periodic(
            LennardJonesParams {
                epsilon: 1.0,
                sigma: 1.0,
                cutoff: 1.5,
            },
            SimulationBox::new(10.0, 10.0, 10.0).unwrap(),
            false,
        );

        let report = compute_lennard_jones_forces_with_options(&mut state, &options).unwrap();

        assert_eq!(report.pair_count, 1);
        assert!(report.potential_energy.abs() < 1.0e-12);
        assert!((state.fx[0] - 24.0).abs() < 1.0e-12);
        assert!((state.fx[1] + 24.0).abs() < 1.0e-12);
    }

    #[test]
    fn open_boundary_ignores_same_particles_without_minimum_image() {
        let mut state = SystemState::new(2);
        state.x[0] = 0.75;
        state.x[1] = 9.75;
        let params = LennardJonesParams {
            epsilon: 1.0,
            sigma: 1.0,
            cutoff: 1.5,
        };

        let report = compute_lennard_jones_forces(&mut state, &params).unwrap();

        assert_eq!(report.pair_count, 0);
        assert_eq!(report.potential_energy, 0.0);
    }

    #[test]
    fn shifted_potential_is_zero_at_cutoff() {
        let mut state = two_particle_system(2.5);
        let options = LennardJonesOptions {
            params: LennardJonesParams {
                epsilon: 1.0,
                sigma: 1.0,
                cutoff: 2.5,
            },
            boundary: BoundaryCondition::Open,
            shift_potential: true,
        };

        let report = compute_lennard_jones_forces_with_options(&mut state, &options).unwrap();

        assert_eq!(report.pair_count, 1);
        assert!(report.potential_energy.abs() < 1.0e-12);
    }

    #[test]
    fn pbc_rejects_cutoff_larger_than_half_shortest_box_length() {
        let options = LennardJonesOptions::periodic(
            LennardJonesParams {
                epsilon: 1.0,
                sigma: 1.0,
                cutoff: 3.1,
            },
            SimulationBox::new(10.0, 8.0, 6.0).unwrap(),
            false,
        );

        let error = options.validate().unwrap_err();

        assert!(matches!(
            error,
            ForceError::CutoffTooLargeForPeriodicBox { .. }
        ));
    }

    #[test]
    fn neighbor_list_matches_naive_for_open_boundary() {
        let mut naive = SystemState::new(5);
        naive.x = vec![0.0, 0.9, 1.8, 4.0, 4.9];
        naive.y = vec![0.0, 0.1, -0.1, 0.0, 0.2];
        let mut neighbor = naive.clone();

        let options = LennardJonesOptions::open(LennardJonesParams {
            epsilon: 1.0,
            sigma: 1.0,
            cutoff: 2.0,
        });
        let neighbor_list = NeighborList::build(
            &neighbor,
            NeighborListConfig {
                cutoff: 2.0,
                skin: 0.3,
                boundary: NeighborBoundary::Open,
            },
        )
        .unwrap();

        let naive_report = compute_lennard_jones_forces_with_options(&mut naive, &options).unwrap();
        let neighbor_report = compute_lennard_jones_forces_with_neighbor_list(
            &mut neighbor,
            &options,
            &neighbor_list,
        )
        .unwrap();

        assert_eq!(neighbor_report.pair_count, naive_report.pair_count);
        assert!((neighbor_report.potential_energy - naive_report.potential_energy).abs() < 1.0e-12);
        assert_forces_close(&neighbor, &naive);
    }

    #[test]
    fn neighbor_list_matches_naive_for_periodic_boundary() {
        let simulation_box = SimulationBox::new(8.0, 8.0, 8.0).unwrap();
        let mut naive = SystemState::new(5);
        naive.x = vec![0.2, 7.6, 3.0, 4.2, 4.8];
        naive.y = vec![0.2, 0.3, 3.0, 3.1, 4.8];
        naive.z = vec![0.2, 0.2, 3.0, 2.9, 4.8];
        let mut neighbor = naive.clone();

        let options = LennardJonesOptions::periodic(
            LennardJonesParams {
                epsilon: 1.0,
                sigma: 1.0,
                cutoff: 2.0,
            },
            simulation_box,
            true,
        );
        let neighbor_list = NeighborList::build(
            &neighbor,
            NeighborListConfig {
                cutoff: 2.0,
                skin: 0.3,
                boundary: NeighborBoundary::Periodic(simulation_box),
            },
        )
        .unwrap();

        let naive_report = compute_lennard_jones_forces_with_options(&mut naive, &options).unwrap();
        let neighbor_report = compute_lennard_jones_forces_with_neighbor_list(
            &mut neighbor,
            &options,
            &neighbor_list,
        )
        .unwrap();

        assert_eq!(neighbor_report.pair_count, naive_report.pair_count);
        assert!((neighbor_report.potential_energy - naive_report.potential_energy).abs() < 1.0e-12);
        assert_forces_close(&neighbor, &naive);
    }

    #[test]
    fn linear_pair_indices_cover_unique_pairs_in_order() {
        let pairs: Vec<_> = (0..unique_pair_count(5))
            .map(|index| pair_indices_from_linear_index(index, 5))
            .collect();

        assert_eq!(
            pairs,
            vec![
                (0, 1),
                (0, 2),
                (0, 3),
                (0, 4),
                (1, 2),
                (1, 3),
                (1, 4),
                (2, 3),
                (2, 4),
                (3, 4),
            ]
        );
        assert_eq!(
            index_chunks(unique_pair_count(5), 4),
            vec![(0, 4), (4, 8), (8, 10)]
        );
    }

    #[test]
    fn parallel_naive_matches_serial_naive_for_periodic_boundary() {
        let simulation_box = SimulationBox::new(9.0, 9.0, 9.0).unwrap();
        let mut serial = lattice_system(4, 1.4);
        let mut parallel = serial.clone();
        let options = LennardJonesOptions::periodic(
            LennardJonesParams {
                epsilon: 1.0,
                sigma: 1.0,
                cutoff: 2.5,
            },
            simulation_box,
            true,
        );

        let serial_report =
            compute_lennard_jones_forces_with_options(&mut serial, &options).unwrap();
        let parallel_report =
            compute_lennard_jones_forces_parallel_with_options(&mut parallel, &options).unwrap();

        assert_eq!(parallel_report.pair_count, serial_report.pair_count);
        assert!(
            (parallel_report.potential_energy - serial_report.potential_energy).abs() < 1.0e-10
        );
        assert_forces_close_with_tolerance(&parallel, &serial, 1.0e-10);
    }

    #[test]
    fn parallel_neighbor_list_matches_serial_neighbor_list() {
        let simulation_box = SimulationBox::new(9.0, 9.0, 9.0).unwrap();
        let mut serial = lattice_system(4, 1.4);
        let mut parallel = serial.clone();
        let options = LennardJonesOptions::periodic(
            LennardJonesParams {
                epsilon: 1.0,
                sigma: 1.0,
                cutoff: 2.5,
            },
            simulation_box,
            true,
        );
        let neighbor_list = NeighborList::build(
            &serial,
            NeighborListConfig {
                cutoff: 2.5,
                skin: 0.3,
                boundary: NeighborBoundary::Periodic(simulation_box),
            },
        )
        .unwrap();

        let serial_report =
            compute_lennard_jones_forces_with_neighbor_list(&mut serial, &options, &neighbor_list)
                .unwrap();
        let parallel_report = compute_lennard_jones_forces_parallel_with_neighbor_list(
            &mut parallel,
            &options,
            &neighbor_list,
        )
        .unwrap();

        assert_eq!(parallel_report.pair_count, serial_report.pair_count);
        assert!(
            (parallel_report.potential_energy - serial_report.potential_energy).abs() < 1.0e-10
        );
        assert_forces_close_with_tolerance(&parallel, &serial, 1.0e-10);
    }

    fn assert_forces_close(left: &SystemState, right: &SystemState) {
        assert_forces_close_with_tolerance(left, right, 1.0e-12);
    }

    fn assert_forces_close_with_tolerance(left: &SystemState, right: &SystemState, tolerance: f64) {
        for i in 0..left.particle_count() {
            assert!((left.fx[i] - right.fx[i]).abs() < tolerance);
            assert!((left.fy[i] - right.fy[i]).abs() < tolerance);
            assert!((left.fz[i] - right.fz[i]).abs() < tolerance);
        }
    }

    fn lattice_system(axis_count: usize, spacing: f64) -> SystemState {
        let particle_count = axis_count * axis_count * axis_count;
        let mut state = SystemState::new(particle_count);
        for i in 0..particle_count {
            let ix = i % axis_count;
            let iy = (i / axis_count) % axis_count;
            let iz = i / (axis_count * axis_count);
            state.x[i] = 1.0 + ix as f64 * spacing;
            state.y[i] = 1.0 + iy as f64 * spacing;
            state.z[i] = 1.0 + iz as f64 * spacing;
        }
        state
    }
}
