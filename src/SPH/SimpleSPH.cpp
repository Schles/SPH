#include "SimpleSPH.h"

SimpleSPH::SimpleSPH(ParticleManager* particleManager, Parameters* params) : SmoothedParticleHydrodynamics(particleManager, params){

  //m_particles_velocity_viscosity_double.resize(particleManager->particlesFluid);
}


void SimpleSPH::step(double deltaTime){

  nhSearch();
  /*
  // compute density at particle positions
  #pragma omp parallel for
  for (int i = 0; i < particleManager->particlesFluid; ++i)
    compute_density(i);

	
  // compute pressure force
  #pragma omp parallel for
  for (int i = 0; i < particleManager->particlesFluid; ++i)
    compute_acceleration(i);

    
  // compute viscosity
  #pragma omp parallel for
  for (int i = 0; i < particleManager->particlesFluid; ++i)
    compute_viscosity(0,i);
  */
  update_semi_implicit_euler(deltaTime);



}

// Computes the density and density factorial
// for the given particle id

void SimpleSPH::compute_density(int _particle_id)
{
  /*
	m_particles_density[_particle_id] = 0;

	CompactNSearch::PointSet const& ps = m_neighborhoodsearch.point_set(0);

	for (int j = 0; j < ps.n_neighbors(_particle_id); j++) {

	  CompactNSearch::PointID nid = ps.neighbor(_particle_id, j);

	  Eigen::Vector3d q = getPosition(0, _particle_id) - getPosition(nid.point_set_id, nid.point_id);

	  // Calculate particle density
	  m_particles_density[_particle_id] += m_kernel->w(q, _parameters->h) * _parameters->particleMass;
	}

	// density factorial
	double density_square = pow(m_particles_density[_particle_id], 2);
	m_particles_density_fac[_particle_id] = (std::max)(_parameters->stiffnessFluid * (m_particles_density[_particle_id] - _parameters->restDensity), 0.0) / density_square;
  */
}




/*
 * Computes the acceleraction of a particles
 * a = a_pressure + f_gravition/m + f_boundry/m
 */
 
void SimpleSPH::compute_acceleration(int _particle_id)
{
  /*
	CompactNSearch::PointSet const& ps = m_neighborhoodsearch.point_set(0);

	m_particles_acceleration[_particle_id] = { 0.0, 0.0, 0.0 };

	bool boundary_hit = false;
	Eigen::Vector3d boundary_force = { 0.0, 0.0, 0.0 };

	for (int j = 0; j < ps.n_neighbors(_particle_id); j++) {

	  CompactNSearch::PointID nid = ps.neighbor(_particle_id, j);
	  Eigen::Vector3d vec_ij = getPosition(0, _particle_id) - getPosition(nid.point_set_id, nid.point_id);
	  
	  if (nid.point_set_id == 0)
	    {
	      // pressure force


	      Eigen::Vector3d q_grad_weighted =  m_kernel->gradient_w(vec_ij, _parameters->h);

	      Eigen::Vector3d f_pressure_i = -1.0 * q_grad_weighted * (m_particles_density_fac[_particle_id] + m_particles_density_fac[nid.point_id]) * _parameters->particleMass;

	      m_particles_acceleration[_particle_id] += f_pressure_i * _parameters->particleMass;
	    }
	  else
	    {
	      // boundary force		  


	      double norm = vec_ij.norm();
		      
	      Eigen::Vector3d bbForce = vec_ij;
 
	      bbForce /= norm;
	      bbForce *= _parameters->particleMass/(_parameters->particleMass + _parameters->particleMassBB); // constant
	      bbForce *= _parameters->stiffnessBB * m_kernel->computeDensityCoeffient(vec_ij, _parameters->h);

	      boundary_force += bbForce;
	    }

	}

	Eigen::Vector3d f_g(0,0,0);

	// test: dont apply gravity
	if(_parameters->useGravity)
	  f_g = Eigen::Vector3d{0, -9.81, 0}; // * _parameters->particleMass;

	//m_color[_particle_id] = m_particles_acceleration[_particle_id] * 50;
	/*
	if(_particle_id == observeParticle)
      	  debugParticle(_particle_id, boundary_force);
	  */
  /*
	m_particles_acceleration[_particle_id] += f_g; // / m_particles_density[_particle_id];

	m_particles_acceleration[_particle_id] += boundary_force / _parameters->particleMass; //  / m_particles_density[_particle_id];
	*/
}

// computes viscosity for a particle id
void SimpleSPH::compute_viscosity(int fluidIndex, int pid){

  Fluid* fluidObject = getFluidObject(fluidIndex);
  
  fluidObject->velocity_viscosity[pid] = Eigen::Vector3d(0.0, 0.0, 0.0);

  CompactNSearch::PointSet const& ps = m_neighborhoodsearch.point_set(fluidIndex);
	
  for (int j = 0; j < ps.n_neighbors(pid); j++) {
	  
    CompactNSearch::PointID nid = ps.neighbor(pid, j);

    if(isFluid(nid.point_set_id)){

      double neighborDensity = getFluidObject(nid.point_set_id)->density[nid.point_id];
      
      if(neighborDensity > 0.0){
	
	Eigen::Vector3d vec_ij = getPosition(fluidIndex, pid) - getPosition(nid.point_set_id, nid.point_id);

	double q = m_kernel->w(vec_ij, _parameters->h);
	
	double a = _parameters->particleMass  * q * _parameters->viscosity / neighborDensity;

	Eigen::Vector3d delta = (particleManager->m_particleObjects[nid.point_set_id]->velocity[nid.point_id] - fluidObject->velocity[pid]);
	    
	// apply viscosity
	fluidObject->velocity[pid] += a * (delta);

      }
    } 
  }

}

