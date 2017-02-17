#include "PBSPH.h"

PBSPH::PBSPH(ParticleManager* particleManager, Parameters* params) : SPH(particleManager, params) {
  boundaryPsi = new BoundaryPsi(&m_neighborhoodsearch, particleManager, params, m_kernel);
}

void PBSPH::compute_position_update(int fluidIndex, int _particle_id){
  Fluid* fluidObject = particleManager->getFluidObject(fluidIndex);
  
  fluidObject->position[_particle_id] += fluidObject->delta_x[_particle_id];
  fluidObject->delta_x[_particle_id] = Eigen::Vector3d(0.0, 0.0, 0.0);
}

void PBSPH::compute_prev_position_update(int fluidIndex, int _particle_id){
  Fluid* fluidObject = particleManager->getFluidObject(fluidIndex);
  
  fluidObject->last_pos[_particle_id] = fluidObject->old_pos[_particle_id];
  fluidObject->old_pos[_particle_id] = fluidObject->position[_particle_id];
}

void PBSPH::compute_velocity(int set, int i, double timeStep){
  Fluid* fluid = particleManager->getFluidObject(set);

  switch(_parameters->veloUpdateMethod){
    
  default:
    fluid->velocity[i] = (1.0 / timeStep) * (  fluid->position[i] - fluid->old_pos[i]);
    break;
  case 1:
    fluid->velocity[i] = (1.0 / timeStep) * ( 1.5 * fluid->position[i] - 2.0 * fluid->old_pos[i] + 0.5 * fluid->last_pos[i]);
    break;
  }

}

void PBSPH::compute_density(int fluidIndex, int _particle_id){
    
  CompactNSearch::PointSet const& ps = m_neighborhoodsearch.point_set(fluidIndex);
  double &density = particleManager->getFluidObject(fluidIndex)->density[_particle_id];

  density = 0.0;
    
  for (int j = 0; j < ps.n_neighbors(_particle_id); j++) {
	  
    CompactNSearch::PointID pid = ps.neighbor(_particle_id, j);

    Eigen::Vector3d q = particleManager->getPosition(fluidIndex, _particle_id) - particleManager->getPosition(pid.point_set_id, pid.point_id);

    // fluid
    if(particleManager->isFluid(ps.neighbor(_particle_id, j).point_set_id)){
      density += m_kernel->w(q, _parameters->h) * _parameters->particleMass;
    } else { //boundry
      density += m_kernel->w(q, _parameters->h) * particleManager->getBoundaryObject(pid.point_set_id)->m_boundaryPsi[pid.point_id];
    }
  }
}

void PBSPH::compute_external_forces(int fluidIndex, int _particle_id){
  
  Eigen::Vector3d &accel = particleManager->getFluidObject(fluidIndex)->acceleration[_particle_id];
  accel = Eigen::Vector3d(0.0, 0.0, 0.0 );
	
  if(_parameters->useGravity)
    accel += _parameters->graviationalForce;

  if(_parameters->add)
    accel += Eigen::Vector3d(10.0, 0.0, 0.0);

  // Friction
  if(_parameters->useFriction)
    compute_friction(fluidIndex, _particle_id);
  
}

void PBSPH::compute_friction(int fluidIndex, int _particle_id){

  CompactNSearch::PointSet const &ps = m_neighborhoodsearch.point_set(fluidIndex);

  Fluid* fluidObject = particleManager->getFluidObject(fluidIndex);
    
  Eigen::Vector3d frictionForce = Eigen::Vector3d(0.0, 0.0, 0.0);

  for (int j = 0; j < ps.n_neighbors(_particle_id); j++)
    {
      CompactNSearch::PointID pid = ps.neighbor(_particle_id, j);

      if (particleManager->isFluid(pid.point_set_id)) continue;

      Eigen::Vector3d x = particleManager->getPosition(fluidIndex, _particle_id) - particleManager->getPosition(pid.point_set_id, pid.point_id);
      Eigen::Vector3d v = particleManager->m_particleObjects[0]->velocity[_particle_id];

      double fac1 = (std::min)(v.dot(x), 0.0);
      double fac2 = x.dot(x) + 0.01 * _parameters->h * _parameters->h;

      double vis_coeff = _parameters->viscosityCoefficient;
      if(fluidObject->density[_particle_id] <= 0.0) continue;
	
      double Pi = -vis_coeff * _parameters->h / (2 * fluidObject->density[_particle_id]) * fac1 / fac2;

      double boundryPsi = particleManager->getBoundaryObject(pid.point_set_id)->m_boundaryPsi[pid.point_id];
      frictionForce += -_parameters->particleMass * boundryPsi * _parameters->restDensity * Pi * m_kernel->gradient_w(x, _parameters->h);
    }

    
  fluidObject->acceleration[_particle_id] += frictionForce;
}


void PBSPH::compute_lambda(int fluidIndex, int particleId){
  CompactNSearch::PointSet const &ps = m_neighborhoodsearch.point_set(fluidIndex);
	
  Fluid* fluidObject = particleManager->getFluidObject(fluidIndex);
	

  double constraint = std::max( (fluidObject->density[particleId] / _parameters->restDensity) - 1.0, 0.0);
	  
  if( constraint != 0.0){	
    double gradC2 = 0.0;
    Eigen::Vector3d gradC_i(0.0, 0.0, 0.0);
    Eigen::Vector3d gradC_j;

    for (int j = 0; j < ps.n_neighbors(particleId); j++) {
      CompactNSearch::PointID nid = ps.neighbor(particleId, j);
	  
      Eigen::Vector3d q = particleManager->getPosition(fluidIndex, particleId) - particleManager->getPosition(nid.point_set_id, nid.point_id);

	      
      if (particleManager->isFluid(nid.point_set_id)){
	gradC_j = -1.0 * _parameters->particleMass * m_kernel->gradient_w(q, _parameters->h) /  _parameters->restDensity;
      } else {
	double boundaryPsi = particleManager->getBoundaryObject(nid.point_set_id)->m_boundaryPsi[nid.point_id];
	gradC_j = -1.0 * boundaryPsi * m_kernel->gradient_w(q, _parameters->h) /  _parameters->restDensity;
      }
	  
      gradC2 += gradC_j.squaredNorm();
      gradC_i -= gradC_j;
    }

    gradC2 += gradC_i.squaredNorm();
	
    fluidObject->lambda[particleId] = -constraint / ( gradC2 + _parameters->pb_epsilon);
  } else
    fluidObject->lambda[particleId] = 0.0;
            
	
}

void PBSPH::compute_delta_x(int fluidIndex, int _particle_id){

  CompactNSearch::PointSet const& ps = m_neighborhoodsearch.point_set(fluidIndex);

  Eigen::Vector3d deltaX(0.0, 0.0, 0.0);
  Eigen::Vector3d gradC_j(0.0, 0.0, 0.0);
  
  //  #pragma omp parallel for
  for (int j = 0; j < ps.n_neighbors(_particle_id); j++) {
    CompactNSearch::PointID nid = ps.neighbor(_particle_id, j);


    unsigned int neighborSet = nid.point_set_id;
    unsigned int neighborId = nid.point_id;
    
    Eigen::Vector3d q = particleManager->getPosition(fluidIndex, _particle_id) - particleManager->getPosition(neighborSet, neighborId);   

    if (particleManager->isFluid(nid.point_set_id)){
      gradC_j = -1.0 * _parameters->particleMass * m_kernel->gradient_w(q, _parameters->h) /  _parameters->restDensity;
      deltaX -= ( particleManager->getFluidObject(fluidIndex)->lambda[_particle_id] + particleManager->getFluidObject(nid.point_set_id)->lambda[nid.point_id] ) * gradC_j;
    } else {
      double boundaryPsi = particleManager->getBoundaryObject(nid.point_set_id)->m_boundaryPsi[nid.point_id];
      double alpha =  ( 2.0 *  particleManager->getFluidObject(fluidIndex)->lambda[_particle_id] );
      gradC_j = -1.0 * boundaryPsi * m_kernel->gradient_w(q, _parameters->h) /  _parameters->restDensity;    
      deltaX -= alpha * gradC_j;

    }

  }

  particleManager->getFluidObject(fluidIndex)->delta_x[_particle_id] += deltaX;
}


void PBSPH::compute_viscosity(int fluidIndex, int pid){

  Fluid* fluidObject = particleManager->getFluidObject(fluidIndex);
  
  CompactNSearch::PointSet const& ps = m_neighborhoodsearch.point_set(fluidIndex);
	
  for (int j = 0; j < ps.n_neighbors(pid); j++) {
	  
    CompactNSearch::PointID nid = ps.neighbor(pid, j);

    if(particleManager->isFluid(nid.point_set_id)){

      double neighborDensity = particleManager->getFluidObject(nid.point_set_id)->density[nid.point_id];
      
      if(neighborDensity > 0.0){
	
	Eigen::Vector3d vec_ij = particleManager->getPosition(fluidIndex, pid) - particleManager->getPosition(nid.point_set_id, nid.point_id);

	double q = m_kernel->w(vec_ij, _parameters->h);
	
	double a = _parameters->particleMass  * q * _parameters->viscosity / neighborDensity;

	Eigen::Vector3d delta = (particleManager->m_particleObjects[nid.point_set_id]->velocity[nid.point_id] - fluidObject->velocity[pid]);
	    
	// apply viscosity
	fluidObject->velocity[pid] += a * (delta);

      }
    } 
  }
}
