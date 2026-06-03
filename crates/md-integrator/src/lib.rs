//! Time integration algorithms.

use md_core::{CoreError, SystemState};
use md_force::{
    compute_lennard_jones_forces_parallel_with_neighbor_list,
    compute_lennard_jones_forces_parallel_with_options,
    compute_lennard_jones_forces_with_neighbor_list, compute_lennard_jones_forces_with_options,
    BoundaryCondition, ForceError, LennardJonesOptions, LennardJonesParams,
};
use md_neighbor::{NeighborError, NeighborList};
use thiserror::Error;

#[derive(Debug, Error)]
pub enum IntegratorError {
    #[error("time step must be positive and finite: {0}")]
    InvalidTimeStep(f64),
    #[error(transparent)]
    Core(#[from] CoreError),
    #[error(transparent)]
    Force(#[from] ForceError),
    #[error(transparent)]
    Neighbor(#[from] NeighborError),
}

/// Velocity Verlet integrator for fixed-timestep NVE simulations.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct VelocityVerlet {
    dt: f64,
}

impl VelocityVerlet {
    pub fn new(dt: f64) -> Result<Self, IntegratorError> {
        if !dt.is_finite() || dt <= 0.0 {
            return Err(IntegratorError::InvalidTimeStep(dt));
        }
        Ok(Self { dt })
    }

    pub fn dt(&self) -> f64 {
        self.dt
    }

    /// Advance the system one step using the current forces as the t forces.
    ///
    /// Call a force calculation once before the first integration step.
    pub fn step_lennard_jones(
        &self,
        system: &mut SystemState,
        params: &LennardJonesParams,
    ) -> Result<f64, IntegratorError> {
        self.step_lennard_jones_with_options(system, &LennardJonesOptions::open(*params))
    }

    pub fn step_lennard_jones_with_options(
        &self,
        system: &mut SystemState,
        options: &LennardJonesOptions,
    ) -> Result<f64, IntegratorError> {
        system.validate()?;
        self.first_half(system);
        if let BoundaryCondition::Periodic(simulation_box) = options.boundary {
            system.wrap_positions(&simulation_box);
        }
        let report = compute_lennard_jones_forces_with_options(system, options)?;
        self.second_half(system);
        system.validate()?;
        Ok(report.potential_energy)
    }

    pub fn step_with_force<F>(
        &self,
        system: &mut SystemState,
        boundary: BoundaryCondition,
        mut compute_forces: F,
    ) -> Result<f64, IntegratorError>
    where
        F: FnMut(&mut SystemState) -> Result<f64, IntegratorError>,
    {
        system.validate()?;
        self.first_half(system);
        if let BoundaryCondition::Periodic(simulation_box) = boundary {
            system.wrap_positions(&simulation_box);
        }
        let potential_energy = compute_forces(system)?;
        self.second_half(system);
        system.validate()?;
        Ok(potential_energy)
    }

    pub fn step_lennard_jones_parallel_with_options(
        &self,
        system: &mut SystemState,
        options: &LennardJonesOptions,
    ) -> Result<f64, IntegratorError> {
        system.validate()?;
        self.first_half(system);
        if let BoundaryCondition::Periodic(simulation_box) = options.boundary {
            system.wrap_positions(&simulation_box);
        }
        let report = compute_lennard_jones_forces_parallel_with_options(system, options)?;
        self.second_half(system);
        system.validate()?;
        Ok(report.potential_energy)
    }

    pub fn step_lennard_jones_with_neighbor_list(
        &self,
        system: &mut SystemState,
        options: &LennardJonesOptions,
        neighbor_list: &mut NeighborList,
        force_rebuild: bool,
    ) -> Result<f64, IntegratorError> {
        system.validate()?;
        self.first_half(system);
        if let BoundaryCondition::Periodic(simulation_box) = options.boundary {
            system.wrap_positions(&simulation_box);
        }
        if force_rebuild || neighbor_list.needs_rebuild(system)? {
            neighbor_list.rebuild(system)?;
        }
        let report =
            compute_lennard_jones_forces_with_neighbor_list(system, options, neighbor_list)?;
        self.second_half(system);
        system.validate()?;
        Ok(report.potential_energy)
    }

    pub fn step_lennard_jones_parallel_with_neighbor_list(
        &self,
        system: &mut SystemState,
        options: &LennardJonesOptions,
        neighbor_list: &mut NeighborList,
        force_rebuild: bool,
    ) -> Result<f64, IntegratorError> {
        system.validate()?;
        self.first_half(system);
        if let BoundaryCondition::Periodic(simulation_box) = options.boundary {
            system.wrap_positions(&simulation_box);
        }
        if force_rebuild || neighbor_list.needs_rebuild(system)? {
            neighbor_list.rebuild(system)?;
        }
        let report = compute_lennard_jones_forces_parallel_with_neighbor_list(
            system,
            options,
            neighbor_list,
        )?;
        self.second_half(system);
        system.validate()?;
        Ok(report.potential_energy)
    }

    fn first_half(&self, system: &mut SystemState) {
        let half_dt = 0.5 * self.dt;
        for i in 0..system.particle_count() {
            let inv_mass = 1.0 / system.mass[i];
            system.vx[i] += half_dt * system.fx[i] * inv_mass;
            system.vy[i] += half_dt * system.fy[i] * inv_mass;
            system.vz[i] += half_dt * system.fz[i] * inv_mass;

            system.x[i] += self.dt * system.vx[i];
            system.y[i] += self.dt * system.vy[i];
            system.z[i] += self.dt * system.vz[i];
        }
    }

    fn second_half(&self, system: &mut SystemState) {
        let half_dt = 0.5 * self.dt;
        for i in 0..system.particle_count() {
            let inv_mass = 1.0 / system.mass[i];
            system.vx[i] += half_dt * system.fx[i] * inv_mass;
            system.vy[i] += half_dt * system.fy[i] * inv_mass;
            system.vz[i] += half_dt * system.fz[i] * inv_mass;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use md_force::{compute_harmonic_bond_forces, compute_lennard_jones_forces, HarmonicBond};
    use md_neighbor::{NeighborBoundary, NeighborListConfig};

    #[test]
    fn rejects_invalid_time_step() {
        assert!(VelocityVerlet::new(0.0).is_err());
        assert!(VelocityVerlet::new(f64::NAN).is_err());
    }

    #[test]
    fn stable_equilibrium_pair_has_small_energy_drift() {
        let mut state = SystemState::new(2);
        state.x[0] = 0.0;
        state.x[1] = 2.0_f64.powf(1.0 / 6.0);

        let params = LennardJonesParams {
            epsilon: 1.0,
            sigma: 1.0,
            cutoff: 2.5,
        };
        let initial = compute_lennard_jones_forces(&mut state, &params)
            .unwrap()
            .potential_energy;
        let integrator = VelocityVerlet::new(0.001).unwrap();

        let mut potential = initial;
        for _ in 0..100 {
            potential = integrator.step_lennard_jones(&mut state, &params).unwrap();
        }

        assert!((potential - initial).abs() < 1.0e-10);
        assert!(state.kinetic_energy() < 1.0e-10);
    }

    #[test]
    fn vibrating_pair_conserves_energy_with_small_time_step() {
        let mut state = SystemState::new(2);
        state.x[0] = 0.0;
        state.x[1] = 1.25;
        state.vx[0] = 0.02;
        state.vx[1] = -0.02;

        let params = LennardJonesParams {
            epsilon: 1.0,
            sigma: 1.0,
            cutoff: 2.5,
        };
        let initial_potential = compute_lennard_jones_forces(&mut state, &params)
            .unwrap()
            .potential_energy;
        let initial_total = state.kinetic_energy() + initial_potential;
        let integrator = VelocityVerlet::new(0.0005).unwrap();

        let mut potential = initial_potential;
        for _ in 0..1000 {
            potential = integrator.step_lennard_jones(&mut state, &params).unwrap();
        }

        let final_total = state.kinetic_energy() + potential;
        assert!((final_total - initial_total).abs() < 1.0e-6);
    }

    #[test]
    fn single_particle_drifts_with_constant_velocity() {
        let mut state = SystemState::new(1);
        state.vx[0] = 0.5;
        state.vy[0] = -0.25;
        state.vz[0] = 0.125;

        let params = LennardJonesParams {
            epsilon: 1.0,
            sigma: 1.0,
            cutoff: 2.5,
        };
        compute_lennard_jones_forces(&mut state, &params).unwrap();
        let integrator = VelocityVerlet::new(0.01).unwrap();

        for _ in 0..10 {
            integrator.step_lennard_jones(&mut state, &params).unwrap();
        }

        assert!((state.x[0] - 0.05).abs() < 1.0e-12);
        assert!((state.y[0] + 0.025).abs() < 1.0e-12);
        assert!((state.z[0] - 0.0125).abs() < 1.0e-12);
        assert!((state.vx[0] - 0.5).abs() < 1.0e-12);
        assert!((state.vy[0] + 0.25).abs() < 1.0e-12);
        assert!((state.vz[0] - 0.125).abs() < 1.0e-12);
    }

    #[test]
    fn generic_force_step_supports_harmonic_bond_terms() {
        let mut state = SystemState::new(2);
        state.x[0] = 0.0;
        state.x[1] = 1.5;
        let bonds = [HarmonicBond {
            i: 0,
            j: 1,
            k: 10.0,
            r0: 1.0,
        }];
        compute_harmonic_bond_forces(&mut state, &bonds, BoundaryCondition::Open).unwrap();
        let integrator = VelocityVerlet::new(0.001).unwrap();

        let potential = integrator
            .step_with_force(&mut state, BoundaryCondition::Open, |state| {
                Ok(
                    compute_harmonic_bond_forces(state, &bonds, BoundaryCondition::Open)?
                        .potential_energy,
                )
            })
            .unwrap();

        assert!(state.x[0] > 0.0);
        assert!(state.x[1] < 1.5);
        assert!(potential < 1.25);
    }

    #[test]
    fn periodic_step_wraps_positions_after_drift() {
        let mut state = SystemState::new(1);
        state.x[0] = 0.99;
        state.vx[0] = 0.2;

        let options = LennardJonesOptions::periodic(
            LennardJonesParams {
                epsilon: 1.0,
                sigma: 0.1,
                cutoff: 0.4,
            },
            md_core::SimulationBox::new(1.0, 1.0, 1.0).unwrap(),
            false,
        );
        compute_lennard_jones_forces_with_options(&mut state, &options).unwrap();
        let integrator = VelocityVerlet::new(0.1).unwrap();

        integrator
            .step_lennard_jones_with_options(&mut state, &options)
            .unwrap();

        assert!((state.x[0] - 0.01).abs() < 1.0e-12);
        assert!((state.vx[0] - 0.2).abs() < 1.0e-12);
    }

    #[test]
    fn neighbor_list_step_rebuilds_when_skin_is_exceeded() {
        let mut state = SystemState::new(2);
        state.x[0] = 0.0;
        state.x[1] = 3.0;
        state.vx[1] = -9.0;

        let options = LennardJonesOptions::open(LennardJonesParams {
            epsilon: 1.0,
            sigma: 1.0,
            cutoff: 1.2,
        });
        let neighbor_config = NeighborListConfig {
            cutoff: 1.2,
            skin: 0.2,
            boundary: NeighborBoundary::Open,
        };
        let mut neighbor_list = NeighborList::build(&state, neighbor_config).unwrap();
        assert!(neighbor_list.pairs().is_empty());
        let integrator = VelocityVerlet::new(0.2).unwrap();

        let potential = integrator
            .step_lennard_jones_with_neighbor_list(&mut state, &options, &mut neighbor_list, false)
            .unwrap();

        assert_eq!(neighbor_list.rebuild_count(), 2);
        assert_eq!(neighbor_list.pairs(), &[(0, 1)]);
        assert!(potential < 0.0);
    }
}
