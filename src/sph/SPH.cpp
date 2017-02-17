
#include "SPH.h"

#include <iostream>
#include <memory>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>    // std::max

#include "../util/fluidtime.h"

SPH::SPH(ParticleManager* manager, Parameters *params)
	: m_neighborhoodsearch( 2.0 * params->h),
	  _parameters(params),
	  particleManager(manager)
{
  m_statistics = new Statistics();
  m_kernel = new Kernel(params->kernelFunctionId, params->h);
  m_particleObjects = &manager->m_particleObjects;

  // Add particles
  for(int i = 0; i < m_particleObjects->size(); i++){
    if(particleManager->isFluid(i)){
      m_neighborhoodsearch.add_point_set(&particleManager->getPosition(i,0)[0], particleManager->getObjectSize(i), true, true);
    } else {    
      m_neighborhoodsearch.add_point_set( &(particleManager->getBoundaryObject(i)->position[0][0]), particleManager->getObjectSize(i), (*m_particleObjects)[i]->isDynamic, false);
    }
  }
}

//compute neighborhood
void SPH::nhSearch(){
  m_neighborhoodsearch.find_neighbors();
}


// semi-implicit euler update
void SPH::compute_semi_implicit_euler(int fluidIndex, int i, double _timestep){	

  Fluid* fluid = particleManager->getFluidObject(fluidIndex);

  // Update velocity
  particleManager->m_particleObjects[fluidIndex]->velocity[i] += fluid->acceleration[i] * _timestep; 
      
  // Update postion based on delta t
  particleManager->m_particleObjects[fluidIndex]->position[i] +=  particleManager->m_particleObjects[fluidIndex]->velocity[i] * _timestep;
}


// Utility, computes the colors for each particle
void SPH::compute_color(int fluidIndex, int pid){

  Fluid* fluid = particleManager->getFluidObject(fluidIndex);
  particleManager->m_particleObjects[fluidIndex]->color[pid] = Eigen::Vector3d(0.0, 0.0, 1.0) * fluid->density[pid] / ( 2.0 * _parameters->restDensity);

}

void SPH::initParticles(){

  // calc mass of particles, based on particle radius
  double diameter = 2.0 * _parameters->particleRadius;
  _parameters->particleMass = _parameters->particleMassScaling * pow(diameter, 3) * _parameters->restDensity;

  
  // reposition particle at start
  for(int fI = 0; fI < particleManager->fluidIndicies.size(); fI++){
    unsigned int fluidIndex = particleManager->fluidIndicies[fI];
    for(int i = 0; i < particleManager->getObjectSize(fluidIndex); i++){
      
      (*m_particleObjects)[fluidIndex]->position[i] = (*m_particleObjects)[fluidIndex]->position0[i];
      (*m_particleObjects)[fluidIndex]->velocity[i] = Eigen::Vector3d(0.0, 0.0, 0.0);
      (*m_particleObjects)[fluidIndex]->color[i] = Eigen::Vector3d(0.0, 0.0, 0.0);
      
    }
  }
}


void SPH::updateTimeStepCFL(){
  double stepSize = _parameters->stepSize;

  double maxVelo = 0.1;

  double diameter = 2.0 * _parameters->particleRadius;

  for(int fI = 0; fI < particleManager->fluidIndicies.size(); fI++){
    unsigned int fluidIndex = particleManager->fluidIndicies[fI];
  
    for(int i = 0; i < particleManager->getObjectSize(fluidIndex); i++){
      const Eigen::Vector3d a = particleManager->getFluidObject(fluidIndex)->acceleration[i];
      const Eigen::Vector3d v = particleManager->getFluidObject(fluidIndex)->velocity[i];
    
      double velo = (v + a * stepSize).squaredNorm();
    
      if(velo > maxVelo)
	maxVelo = velo;
    }
  }


  stepSize = _parameters->m_cflFactor * 0.4 * (diameter / sqrt(maxVelo));
  
  stepSize = std::min(stepSize, _parameters->m_cflMaxTimeStep);
  stepSize = std::max(stepSize, _parameters->m_cflMinTimeStep);

  _parameters->stepSize = stepSize;
  

}

void SPH::compute_stats(){
      double avg_density = 0;
      double max_density = 0;

      double min_velo = 1000;
      double avg_velo = 0;
      double max_velo = 0;

      double avg_fluid_neighbors = 0.0;
      double avg_neighbors = 0.0;

      double fluidSize = 0.0000;
      
      for(int fI = 0; fI < particleManager->fluidIndicies.size(); fI++){
	unsigned int fluidIndex = particleManager->fluidIndicies[fI];

	Fluid* fluid = particleManager->getFluidObject(fluidIndex);

	fluidSize += particleManager->getObjectSize(fluidIndex);
	
	CompactNSearch::PointSet const &ps = m_neighborhoodsearch.point_set(fluidIndex);
	
	for (int i = 0; i < particleManager->getObjectSize(fluidIndex); ++i){

	  // density
	  avg_density += fluid->density[i];

	  if(fluid->density[i] > max_density)
	    max_density = fluid->density[i];
	
	  // velo
	
	  double velo = particleManager->m_particleObjects[fluidIndex]->velocity[i].norm();
	
	  avg_velo += velo;

	  if(velo < min_velo)
	    min_velo = velo;

	  if(velo > max_velo)
	    max_velo = velo;


	  for (int j = 0; j < ps.n_neighbors(i); j++) {
	    CompactNSearch::PointID nid = ps.neighbor(i, j);
	    if (nid.point_set_id == 0){
	      avg_fluid_neighbors++;
	    } 
	    avg_neighbors++;
	  }
	}
      }
      
      m_statistics->avg_density = avg_density / fluidSize;
      m_statistics->max_density = max_density;


      m_statistics->min_velo = min_velo;
      m_statistics->avg_velo = avg_velo / fluidSize;
      m_statistics->max_velo = max_velo;
      
      m_statistics->avg_neighbors = avg_neighbors / fluidSize;
      m_statistics->avg_fluid_neighbors = avg_fluid_neighbors / fluidSize;
      
}


bool SPH::debugP(int id){
  if( id == _parameters->observeParticle)
    return true;

  return false;
}
