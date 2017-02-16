#include "PBSPH.h"

PBSPH::PBSPH(ParticleManager* particleManager, Parameters* params) : ImprSPH(particleManager, params) {
  /*
  m_lambda.resize(particleManager->particlesFluid);
  m_delta_x.resize(particleManager->particlesFluid);
  m_old_pos.resize(particleManager->particlesFluid);
  m_last_pos.resize(particleManager->particlesFluid);
  */
  m_vector_output.resize(visualizeAmount);

  for(int i = 0; i < visualizeAmount; i++){
    m_vector_output[i].resize(2 * getObjectSize(0));
  }
  /*
  m_delta_x_fluid.resize(particleManager->particlesFluid);
  m_delta_x_bb.resize(particleManager->particlesFluid);
  */

  init();
}

void PBSPH::step(double deltaTime){

  long t1 = Time::getMilliseconds();
  
  // update by delta x, reset delta_x vectors
  //  #pragma omp parallel for


  
  for(int fI = 0; fI < particleManager->fluidIndicies.size(); fI++){
    unsigned int fluidIndex = particleManager->fluidIndicies[fI];

    Fluid* fluidObject = getFluidObject(fluidIndex);
    
    for (int i = 0; i < getObjectSize(fluidIndex); ++i){
      
      fluidObject->position[i] += fluidObject->delta_x[i];
      
      fluidObject->delta_x[i] = Eigen::Vector3d(0.0, 0.0, 0.0);
      /*
      m_delta_x_bb[i] = Eigen::Vector3d(0.0, 0.0, 0.0);
      m_delta_x_fluid[i] = Eigen::Vector3d(0.0, 0.0, 0.0);
      */
    }
  }
  
  // change velocity, based on delta x
  //  #pragma omp parallel for
  for(int fI = 0; fI < particleManager->fluidIndicies.size(); fI++){
    unsigned int fluidIndex = particleManager->fluidIndicies[fI];

    for (int i = 0; i < getObjectSize(fluidIndex); ++i){
      updateVelocity(fluidIndex, i, deltaTime);
    }
  }

  // compute viscosity
  //#pragma omp parallel for
  for(int fI = 0; fI < particleManager->fluidIndicies.size(); fI++){
    unsigned int fluidIndex = particleManager->fluidIndicies[fI];
    
    for (int i = 0; i < getObjectSize(fluidIndex); ++i){
      SimpleSPH::compute_viscosity(fluidIndex, i); 
    }
  }

  
  // store old positions for next velocity update
  //  #pragma omp parallel for
  for(int fI = 0; fI < particleManager->fluidIndicies.size(); fI++){
    unsigned int fluidIndex = particleManager->fluidIndicies[fI];

    Fluid* fluidObject = getFluidObject(fluidIndex);
    
    for (int i = 0; i < getObjectSize(fluidIndex); ++i){

      
      fluidObject->last_pos[i] = fluidObject->old_pos[i];
      //      m_last_pos[i] = m_old_pos[i];

      fluidObject->old_pos[i] = fluidObject->position[i];
      //      m_old_pos[i] = (*m_particleObjects)[0]->position[i];
    }
  }


  if(_parameters->enableAdaptiveTimeStep)
    updateTimeStepCFL();
  else
    _parameters->stepSize = _parameters->m_defaultStepSize;

  // Apply external forces to particles
  //  #pragma omp parallel for

  for(int fI = 0; fI < particleManager->fluidIndicies.size(); fI++){
    unsigned int fluidIndex = particleManager->fluidIndicies[fI];

    for (int i = 0; i < getObjectSize(fluidIndex); ++i)
      externalForces(fluidIndex, i);
  }
  
  
  //Advance timestep
  update_semi_implicit_euler(deltaTime);  
  
  // Neighborhoodsearch
  nhSearch();

  
  // Position based pressure solving
  solvePressure();
  
  // Compute color codings
  compute_color();

  // compute stats for TW Bar
  compute_stats();

  // compute displacement vectors to render
  updateDeltaXOutput();

  long t2 = Time::getMilliseconds();

  _parameters->timePerFrame = t2 - t1;
}



void PBSPH::updateVelocity(unsigned int set, unsigned int i, double timeStep){
  int strat = 0;

  Fluid* fluid = getFluidObject(set);
  
  if(strat == 0)
    fluid->velocity[i] = (1.0 / timeStep) * (  fluid->position[i] - fluid->old_pos[i]);
  else
    fluid->velocity[i] = (1.0 / timeStep) * ( 1.5 * fluid->position[i] - 2.0 * fluid->old_pos[i] + 0.5 * fluid->last_pos[i]);

}

void PBSPH::compute_density(int fluidIndex, int _particle_id){
    
    CompactNSearch::PointSet const& ps = m_neighborhoodsearch.point_set(fluidIndex);

    double tmp;

    double &density = getFluidObject(fluidIndex)->density[_particle_id];

    density = 0.0;
    
    for (int j = 0; j < ps.n_neighbors(_particle_id); j++) {
	  
      CompactNSearch::PointID pid = ps.neighbor(_particle_id, j);

      Eigen::Vector3d q = getPosition(fluidIndex, _particle_id) - getPosition(pid.point_set_id, pid.point_id);

      // fluid
      if(isFluid(ps.neighbor(_particle_id, j).point_set_id)){
	tmp = _parameters->particleMass * m_kernel->w(q, _parameters->h);
	density += tmp;
      } else { //boundry
	double boundryPsi = getBoundaryObject(pid.point_set_id)->m_boundaryPsi[pid.point_id];
	tmp = boundryPsi * m_kernel->w(q, _parameters->h);
	density += tmp;
      }
    }
}

void PBSPH::externalForces(int fluidIndex, int _particle_id){
  
  Eigen::Vector3d &accel = getFluidObject(fluidIndex)->acceleration[_particle_id];
  accel = Eigen::Vector3d(0.0, 0.0, 0.0 );
	
  if(_parameters->useGravity)
    accel += _parameters->graviationalForce;

  if(_parameters->useFriction)
    friction(fluidIndex, _particle_id);
    /*
  if(_parameters->add)
    m_particles_acceleration[_particle_id] += Eigen::Vector3d(10.0, 0.0, 0.0);
  */
}

void PBSPH::friction(int fluidIndex, int _particle_id)
{

    CompactNSearch::PointSet const &ps = m_neighborhoodsearch.point_set(fluidIndex);

    Fluid* fluidObject = getFluidObject(fluidIndex);
    
    Eigen::Vector3d frictionForce = Eigen::Vector3d(0.0, 0.0, 0.0);

    for (int j = 0; j < ps.n_neighbors(_particle_id); j++)
    {
        CompactNSearch::PointID pid = ps.neighbor(_particle_id, j);

        if (isFluid(pid.point_set_id)) continue;

        Eigen::Vector3d x = getPosition(fluidIndex, _particle_id) - getPosition(pid.point_set_id, pid.point_id);
        Eigen::Vector3d v = particleManager->m_particleObjects[0]->velocity[_particle_id];

        double fac1 = (std::min)(v.dot(x), 0.0);
        double fac2 = x.dot(x) + 0.01 * _parameters->h * _parameters->h;

        double vis_coeff = _parameters->viscosityCoefficient;
	if(fluidObject->density[_particle_id] <= 0.0) continue;
	
	double Pi = -vis_coeff * _parameters->h / (2 * fluidObject->density[_particle_id]) * fac1 / fac2;

        double boundryPsi = getBoundaryObject(pid.point_set_id)->m_boundaryPsi[pid.point_id];
        frictionForce += -_parameters->particleMass * boundryPsi * _parameters->restDensity * Pi * m_kernel->gradient_w(x, _parameters->h);
    }

    
    fluidObject->acceleration[_particle_id] += frictionForce;
}

void PBSPH::solvePressure(){
  unsigned int iter = 0;

  unsigned int maxIter = _parameters->pb_max_iter;



  while( iter < maxIter){


    //    #pragma omp parallel default(shared)
    {

      // compute density at particle positions
      //      #pragma omp for schedule(static)

      for(int fI = 0; fI < particleManager->fluidIndicies.size(); fI++){
	unsigned int fluidIndex = particleManager->fluidIndicies[fI];
	
	CompactNSearch::PointSet const &ps = m_neighborhoodsearch.point_set(fluidIndex);
	
	Fluid* fluidObject = getFluidObject(fluidIndex);
	
	for (int i = 0; i < getObjectSize(fluidIndex); ++i){
	  compute_density(fluidIndex, i);

	  double constraint = std::max( (fluidObject->density[i] / _parameters->restDensity) - 1.0, 0.0);
	  
	  if( constraint != 0.0){
	
	    double gradC2 = 0.0;
	    Eigen::Vector3d gradC_i(0.0, 0.0, 0.0);
	    Eigen::Vector3d gradC_j;

	    for (int j = 0; j < ps.n_neighbors(i); j++) {
	      CompactNSearch::PointID nid = ps.neighbor(i, j);
	  
	      Eigen::Vector3d q = getPosition(fluidIndex, i) - getPosition(nid.point_set_id, nid.point_id);

	      
	      if (isFluid(nid.point_set_id)){
		gradC_j = -1.0 * _parameters->particleMass * m_kernel->gradient_w(q, _parameters->h) /  _parameters->restDensity;
	      } else {
		double boundaryPsi = getBoundaryObject(nid.point_set_id)->m_boundaryPsi[nid.point_id];
		gradC_j = -1.0 * boundaryPsi * m_kernel->gradient_w(q, _parameters->h) /  _parameters->restDensity;
	      }
	  
	      gradC2 += gradC_j.squaredNorm();
	      gradC_i -= gradC_j;
	    }

	    gradC2 += gradC_i.squaredNorm();
	
	    fluidObject->lambda[i] = -constraint / ( gradC2 + _parameters->pb_epsilon);
	  } else
	    fluidObject->lambda[i] = 0.0;
            
	}
      }
      
	// Compute delta X
      //  #pragma omp for schedule(static)
      for(int fI = 0; fI < particleManager->fluidIndicies.size(); fI++){
	unsigned int fluidIndex = particleManager->fluidIndicies[fI];
	
	for (int i = 0; i < getObjectSize(fluidIndex); ++i){
	  getFluidObject(fluidIndex)->delta_x[i] += compute_delta_x(fluidIndex, i);
	}
      }
	iter++;
      
    }
  }
   
  
}

Eigen::Vector3d PBSPH::compute_delta_x(int fluidIndex, int _particle_id){

  CompactNSearch::PointSet const& ps = m_neighborhoodsearch.point_set(fluidIndex);

  Eigen::Vector3d deltaX(0.0, 0.0, 0.0);
  Eigen::Vector3d gradC_j(0.0, 0.0, 0.0);
  
  //  #pragma omp parallel for
  for (int j = 0; j < ps.n_neighbors(_particle_id); j++) {
    CompactNSearch::PointID nid = ps.neighbor(_particle_id, j);


    unsigned int neighborSet = nid.point_set_id;
    unsigned int neighborId = nid.point_id;
    
    Eigen::Vector3d q = getPosition(fluidIndex, _particle_id) - getPosition(neighborSet, neighborId);   

    if (isFluid(nid.point_set_id)){
      gradC_j = -1.0 * _parameters->particleMass * m_kernel->gradient_w(q, _parameters->h) /  _parameters->restDensity;
      deltaX -= ( getFluidObject(fluidIndex)->lambda[_particle_id] + getFluidObject(nid.point_set_id)->lambda[nid.point_id] ) * gradC_j;
    } else {
      double boundaryPsi = getBoundaryObject(nid.point_set_id)->m_boundaryPsi[nid.point_id];
      double alpha =  ( 2.0 *  getFluidObject(fluidIndex)->lambda[_particle_id] );
      gradC_j = -1.0 * boundaryPsi * m_kernel->gradient_w(q, _parameters->h) /  _parameters->restDensity;    
      deltaX -= alpha * gradC_j;

    }

  }

  return deltaX;
}


void PBSPH::init(){
  std::cout << "Reset particles" << std::endl;


  initParticles();
  calc_boundary_psi();

  for(int fI = 0; fI < particleManager->fluidIndicies.size(); fI++){
    unsigned int fluidIndex = particleManager->fluidIndicies[fI];

    Fluid* fluid = getFluidObject(fluidIndex);
    
    for (int i = 0; i < getObjectSize(fluidIndex); ++i){
      fluid->last_pos[i] = getPosition(fluidIndex,i);
      fluid->old_pos[i] = getPosition(fluidIndex,i);
      fluid->delta_x[i] = Eigen::Vector3d(0.0, 0.0, 0.0);
      
      for(int j = 0; j < visualizeAmount; j++){
	m_vector_output[j][2*i] = Eigen::Vector3f(0.0, 0.0, 0.0);
	m_vector_output[j][2*i+1] = Eigen::Vector3f(0.0, 0.0, 0.0);
      }
      
    }
    
  }


  //  #pragma omp parallel for
  for(int fI = 0; fI < particleManager->fluidIndicies.size(); fI++){
    unsigned int fluidIndex = particleManager->fluidIndicies[fI];
    
    for (int i = 0; i < getObjectSize(fluidIndex); ++i){
      compute_density(fluidIndex, i);
    }
  }
  
  compute_stats();
}


void PBSPH::updateDeltaXOutput(){

  Fluid* fluid = getFluidObject(0);
  
  for(int j = 0; j < visualizeAmount; j++){
    for(int i = 0; i < getObjectSize(0); i++){
      m_vector_output[j][2*i][0] = (float) fluid->position[i][0];
      m_vector_output[j][2*i][1] = (float) fluid->position[i][1];
      m_vector_output[j][2*i][2] = (float) fluid->position[i][2];

      if(j == 0){
	m_vector_output[j][2*i + 1][0] = (float) (fluid->position[i][0] + fluid->delta_x[i][0]);
	m_vector_output[j][2*i + 1][1] = (float) (fluid->position[i][1] + fluid->delta_x[i][1]);
	m_vector_output[j][2*i + 1][2] = (float) (fluid->position[i][2] + fluid->delta_x[i][2]);
      }
      /*
      if(j == 1){
	m_vector_output[j][2*i + 1][0] = (float) (particleManager->m_particleObjects[0]->position[i][0] + m_delta_x_fluid[i][0]);
	m_vector_output[j][2*i + 1][1] = (float) (particleManager->m_particleObjects[0]->position[i][1] + m_delta_x_fluid[i][1]);
	m_vector_output[j][2*i + 1][2] = (float) (particleManager->m_particleObjects[0]->position[i][2] + m_delta_x_fluid[i][2]);
      }

      if(j == 2){
	m_vector_output[j][2*i + 1][0] = (float) (particleManager->m_particleObjects[0]->position[i][0] + m_delta_x_bb[i][0]);
	m_vector_output[j][2*i + 1][1] = (float) (particleManager->m_particleObjects[0]->position[i][1] + m_delta_x_bb[i][1]);
	m_vector_output[j][2*i + 1][2] = (float) (particleManager->m_particleObjects[0]->position[i][2] + m_delta_x_bb[i][2]);
      }
*/
    }
  }
  

}
