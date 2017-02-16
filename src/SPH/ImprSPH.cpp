

#include "ImprSPH.h"

ImprSPH::ImprSPH(ParticleManager* particleManager, Parameters* params) : SimpleSPH(particleManager, params) {
  calc_boundary_psi();
}


void ImprSPH::step(double deltaTime){

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
    SimpleSPH::compute_viscosity(0,i);
  */
  update_semi_implicit_euler(deltaTime);


}

void ImprSPH::compute_acceleration(int _particle_id){
  /*
	CompactNSearch::PointSet const& ps = m_neighborhoodsearch.point_set(0);

	m_particles_acceleration[_particle_id] = { 0.0, 0.0, 0.0 };

	bool boundary_hit = false;
	Eigen::Vector3d boundary_force = { 0.0, 0.0, 0.0 };
	
	for (int j = 0; j < ps.n_neighbors(_particle_id); j++) {

	  CompactNSearch::PointID nid = ps.neighbor(_particle_id, j);
		
	  if (nid.point_set_id == 0) {
		  // pressure force

	    Eigen::Vector3d vec_ij = getPosition(0, _particle_id) - getPosition(nid.point_set_id, nid.point_id);
	    Eigen::Vector3d f_pressure_i = -1.0* m_kernel->gradient_w(vec_ij, _parameters->h) * (m_particles_density_fac[_particle_id] + m_particles_density_fac[nid.point_id]) * _parameters->particleMass;

	    m_particles_acceleration[_particle_id] += f_pressure_i;
	  }
	  else {
	      // boundary force
	      //		  boundary_force += compute_boundary_force(_particle_id, neighbor_id); // / _parameters->particleMass;
	      boundary_hit = true;

	  }
	
	}
	
	Eigen::Vector3d f_g(0,0,0);
	

	m_particles_acceleration[_particle_id] += f_g;

	m_particles_acceleration[_particle_id] += boundary_force / _parameters->particleMass;
  */	
  
}

void ImprSPH::compute_density(int _particle_id){
  /*
  CompactNSearch::PointSet const& ps = m_neighborhoodsearch.point_set(0);

	m_particles_density[_particle_id] = 0.0;
	
	for (int j = 0; j < ps.n_neighbors(_particle_id); j++) {

	  CompactNSearch::PointID nid = ps.neighbor(_particle_id, j);
	  Eigen::Vector3d q = getPosition(0, _particle_id) - getPosition(nid.point_set_id, nid.point_id);
	
		// fluid
	  if (nid.point_set_id == 0){
	    m_particles_density[_particle_id] += _parameters->particleMass * m_kernel->w(q, _parameters->h);
	  } else { //boundry
	    m_particles_density[_particle_id] += _parameters->particleMassBB * m_kernel->w(q, _parameters->h);
	  }

	}
	double density_square = pow(m_particles_density[_particle_id], 2);
	double pressure = _parameters->stiffnessFluid * (pow(m_particles_density[_particle_id] / _parameters->restDensity, 7) - 1.0);
	m_particles_density_fac[_particle_id] = pressure / density_square;

	//	m_particles_density_fac[_particle_id] = (std::max)(_parameters->stiffnessFluid * (m_particles_density[_particle_id] - _parameters->restDensity), 0.0) / density_square;
	
	*/	
}

// computes the boundry force between two given particle ids

Eigen::Vector3d ImprSPH::compute_boundary_force(int _i, int _j)
{
  
  Eigen::Vector3d res;
  
  /*        unsigned int bbForceCase = 2;
	Eigen::Vector3d vec_ij = deltaVector(_i, _j);


	if(bbForceCase == 1){
	  res = m_kernel->gradient_w(vec_ij, _parameters->h);
	  res *= -1.0 * _parameters->particleMass * _parameters->particleMassBB * m_particles_density_fac[_i];
	} else if(bbForceCase == 2){
	  res = m_kernel->gradient_w(vec_ij, _parameters->h);
	  res *= -1.0 * _parameters->particleMass * m_boundary_psi[_j] * m_particles_density_fac[_i];
	}
  */
	return res;
}

double ImprSPH::compute_boundary_volume(int _particle_id){
  double volume = 0;
}


void ImprSPH::calc_boundary_psi_i(unsigned int pointSet){
    CompactNSearch::PointSet const& ps = m_neighborhoodsearch.point_set(pointSet);

    Boundary* boundaryObject = getBoundaryObject(pointSet);
    
    for(unsigned int i = 0; i < (*m_particleObjects)[pointSet]->position.size(); i++){

      double delta = 0;
      
      for (int k = 0; k < ps.n_neighbors(i); k++) {

	CompactNSearch::PointID nid = ps.neighbor(i, k);
	
	Eigen::Vector3d q = getPosition(pointSet, i) - getPosition(nid.point_set_id, nid.point_id);

	//		std::cout << q.norm() << std::endl;
		
	// boundary
	//      std::cout << "nID " << neighbor_id << "type " << (*m_particles_attribute)[neighbor_id] << std::endl;
	if(!isFluid(nid.point_set_id)){
	  delta += m_kernel->w(q, _parameters->h);
	  //	std::cout << "Hallo" << delta << std::endl;
	}

	
      }
      
      double volume = 0.0;

      if(delta > 0)
	volume = 1.0 / delta;
      
      boundaryObject->m_boundaryPsi[i] = _parameters->restDensity * volume;
    }
    
}

double ImprSPH::calc_boundary_psi(){

  // disable all 
  for(int pointSet = 0; pointSet < m_particleObjects->size(); pointSet++){
      m_neighborhoodsearch.point_set(pointSet).enable_neighborsearch(false);
  }

  // search only static
  for(int pointSet = 0; pointSet < m_particleObjects->size(); pointSet++){
    if(! ( (*m_particleObjects)[pointSet]->isDynamic || (*m_particleObjects)[pointSet]->isFluid))
      m_neighborhoodsearch.point_set(pointSet).enable_neighborsearch(true);
  }

  m_neighborhoodsearch.find_neighbors();

  for(int pointSet = 0; pointSet < m_particleObjects->size(); pointSet++){
    if(! ( (*m_particleObjects)[pointSet]->isDynamic || (*m_particleObjects)[pointSet]->isFluid))
      calc_boundary_psi_i(pointSet);
  }

  
  
  for(int pointSet = 0; pointSet < m_particleObjects->size(); pointSet++){
    for(int j = 0; j < m_particleObjects->size(); j++){
      m_neighborhoodsearch.point_set(j).enable_neighborsearch(false);
    }

    if((*m_particleObjects)[pointSet]->isDynamic && !(*m_particleObjects)[pointSet]->isFluid){
      m_neighborhoodsearch.point_set(pointSet).enable_neighborsearch(true);
      m_neighborhoodsearch.find_neighbors();
      calc_boundary_psi_i(pointSet);
    }

  }

  // Activate Fluids
  for(int pointSet = 0; pointSet < m_particleObjects->size(); pointSet++){
    if( (*m_particleObjects)[pointSet]->isFluid)
      m_neighborhoodsearch.point_set(pointSet).enable_neighborsearch(true);
    else
      m_neighborhoodsearch.point_set(pointSet).enable_neighborsearch(false);

  }


}
