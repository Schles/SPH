#include "TimeStep.h"

TimeStepPBSPH::TimeStepPBSPH(ParticleManager* particleManager, Parameters *params) : PBSPH(particleManager, params){
  m_vector_output.resize(visualizeAmount);

  for(int i = 0; i < visualizeAmount; i++){
    m_vector_output[i].resize(2 * particleManager->getObjectSize(0));
  }

  init();
}


void TimeStepPBSPH::init(){
  std::cout << "Reset particles" << std::endl;

  initParticles();

  // Calculate Boundary Psi
  boundaryPsi->calc_boundary_psi();

  for(int fI = 0; fI < particleManager->fluidIndicies.size(); fI++){
    unsigned int fluidIndex = particleManager->fluidIndicies[fI];
    Fluid* fluid = particleManager->getFluidObject(fluidIndex);
    for (int i = 0; i < particleManager->getObjectSize(fluidIndex); ++i){
      
      fluid->last_pos[i] = particleManager->getPosition(fluidIndex,i);
      fluid->old_pos[i] = particleManager->getPosition(fluidIndex,i);
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
    for (int i = 0; i < particleManager->getObjectSize(fluidIndex); ++i){

      //Compute density
      compute_density(fluidIndex, i);
      
    }
  }
  
  compute_stats();
}


void TimeStepPBSPH::step(double deltaTime){

  long t1 = Time::getMilliseconds();

  helperStep(deltaTime);
  
  //  #pragma omp parallel for  
  for(int fI = 0; fI < particleManager->fluidIndicies.size(); fI++){
    unsigned int fluidIndex = particleManager->fluidIndicies[fI];
    for (int i = 0; i < particleManager->getObjectSize(fluidIndex); ++i){

      compute_prev_position_update(fluidIndex, i);
      compute_external_forces(fluidIndex, i);    
      compute_semi_implicit_euler(fluidIndex, i, deltaTime);

      
    }
  }
    
  // Neighborhoodsearch
  nhSearch();

  
  // Position based pressure solving
  solvePressure();

  //  #pragma omp parallel for
  for(int fI = 0; fI < particleManager->fluidIndicies.size(); fI++){
    unsigned int fluidIndex = particleManager->fluidIndicies[fI];
    for (int i = 0; i < particleManager->getObjectSize(fluidIndex); ++i){

      // change velocity, based on delta x
       compute_velocity(fluidIndex, i, deltaTime);
      
    }
  }
  

  //#pragma omp parallel for
  for(int fI = 0; fI < particleManager->fluidIndicies.size(); fI++){
    unsigned int fluidIndex = particleManager->fluidIndicies[fI];
    for (int i = 0; i < particleManager->getObjectSize(fluidIndex); ++i){

      // compute viscosity
      compute_viscosity(fluidIndex, i);
      
    }
  }
  
  
  for(int fI = 0; fI < particleManager->fluidIndicies.size(); fI++){
    unsigned int fluidIndex = particleManager->fluidIndicies[fI];
    for (int i = 0; i < particleManager->getObjectSize(fluidIndex); ++i){

      // Compute color codings
      compute_color(fluidIndex, i);

    }
  }

  // compute stats for TW Bar
  compute_stats();

  // compute displacement vectors to render
  updateDeltaXOutput();

  long t2 = Time::getMilliseconds();

  m_statistics->timePerFrame = t2 - t1;
}


void TimeStepPBSPH::solvePressure(){
  unsigned int iter = 0;
  
  while( iter < _parameters->pb_max_iter){
    
    //    #pragma omp parallel default(shared)
    {

      //      #pragma omp for schedule(static)

      for(int fI = 0; fI < particleManager->fluidIndicies.size(); fI++){
	unsigned int fluidIndex = particleManager->fluidIndicies[fI];	
	for (int i = 0; i < particleManager->getObjectSize(fluidIndex); ++i){

	  // compute density at particle positions
	  compute_density(fluidIndex, i);

	  // compute magnitude of position displacement
	  compute_lambda(fluidIndex, i);

	}
      }
      
      //  #pragma omp for schedule(static)
      for(int fI = 0; fI < particleManager->fluidIndicies.size(); fI++){
	unsigned int fluidIndex = particleManager->fluidIndicies[fI];	
	for (int i = 0; i < particleManager->getObjectSize(fluidIndex); ++i){

	  // Compute delta x, position displacement
	  compute_delta_x(fluidIndex, i);
	  
	}
      }

            //  #pragma omp for schedule(static)
      for(int fI = 0; fI < particleManager->fluidIndicies.size(); fI++){
	unsigned int fluidIndex = particleManager->fluidIndicies[fI];	
	for (int i = 0; i < particleManager->getObjectSize(fluidIndex); ++i){

	  // update by delta x, reset delta_x vectors
	  compute_position_update(fluidIndex, i);
	  
	}
      }

      iter++;
    }
  }
}

void TimeStepPBSPH::updateDeltaXOutput(){

  Fluid* fluid = particleManager->getFluidObject(0);
  
  for(int j = 0; j < visualizeAmount; j++){
    for(int i = 0; i < particleManager->getObjectSize(0); i++){
      m_vector_output[j][2*i][0] = (float) fluid->position[i][0];
      m_vector_output[j][2*i][1] = (float) fluid->position[i][1];
      m_vector_output[j][2*i][2] = (float) fluid->position[i][2];

      if(j == 0){
	m_vector_output[j][2*i + 1][0] = (float) (fluid->position[i][0] + fluid->delta_x[i][0]);
	m_vector_output[j][2*i + 1][1] = (float) (fluid->position[i][1] + fluid->delta_x[i][1]);
	m_vector_output[j][2*i + 1][2] = (float) (fluid->position[i][2] + fluid->delta_x[i][2]);
      }
    }
  }
  

}

void TimeStepPBSPH::helperStep(double timeStep){
  
  if(_parameters->enableAdaptiveTimeStep)
    updateTimeStepCFL();
  else
    _parameters->stepSize = _parameters->m_defaultStepSize;

  double diameter = 2.0 * _parameters->particleRadius;
  _parameters->particleMass = _parameters->particleMassScaling * pow(diameter, 3) * _parameters->restDensity;
}
