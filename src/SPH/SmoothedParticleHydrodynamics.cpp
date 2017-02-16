
#include "SmoothedParticleHydrodynamics.h"

#include <iostream>
#include <memory>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>    // std::max

#include "../util/fluidtime.h"

SmoothedParticleHydrodynamics::SmoothedParticleHydrodynamics(ParticleManager* manager, Parameters *params)
	: m_neighborhoodsearch( 2.0 * params->h),
	  _parameters(params),
	  particleManager(manager)
{


  m_kernel = new Kernel(params->kernelFunctionId, params->h);
  m_particleObjects = &manager->m_particleObjects;
  /*	
  m_particles_pressure_force.resize(m_particles_size);
  
  m_particles_density.resize(manager->particlesFluid);
  
  m_particles_density_fac.resize(m_particles_size);
  
  m_particles_acceleration.resize(manager->particlesFluid);
  m_particles_velocity_viscosity_double.resize(m_particles_size);
*/

  //  initParticles();
  
  // Add fluids

  // add boundrys

  for(int i = 0; i < m_particleObjects->size(); i++){

    if(isFluid(i)){
      m_neighborhoodsearch.add_point_set(&getPosition(i,0)[0], getObjectSize(i), true, true);
    } else {    
      m_neighborhoodsearch.add_point_set( &((*m_particleObjects)[i]->position[0][0]), getObjectSize(i), (*m_particleObjects)[i]->isDynamic, false);
    }
  }
}

Eigen::Vector3d &SmoothedParticleHydrodynamics::getPosition(unsigned int point_set, unsigned int i){
  return (*m_particleObjects)[point_set]->position[i];
}

Boundary* SmoothedParticleHydrodynamics::getBoundaryObject(unsigned int i){
  return static_cast<Boundary*>((*m_particleObjects)[i]);
}

Fluid* SmoothedParticleHydrodynamics::getFluidObject(unsigned int i){
  return static_cast<Fluid*>((*m_particleObjects)[i]);
}

unsigned int SmoothedParticleHydrodynamics::getObjectSize(unsigned int i){
  return (*m_particleObjects)[i]->position.size();
}


bool SmoothedParticleHydrodynamics::isFluid(unsigned int point_set_id){
  for(int i = 0; i < particleManager->fluidIndicies.size(); i++){
    if(point_set_id == particleManager->fluidIndicies[i])
      return true;
  }

  return false;
}


//compute neighborhood
void SmoothedParticleHydrodynamics::nhSearch(){
  m_neighborhoodsearch.find_neighbors();
}


void SmoothedParticleHydrodynamics::update_semi_implicit_euler(double _timestep)
{	
  	// semi-implicit euler update
        for(int fI = 0; fI < particleManager->fluidIndicies.size(); fI++){
	  unsigned int fluidIndex = particleManager->fluidIndicies[fI];

	  Fluid* fluid = getFluidObject(fluidIndex);
	  
	  CompactNSearch::PointSet const& ps = m_neighborhoodsearch.point_set(fluidIndex);
	  for (int i = 0; i < getObjectSize(fluidIndex); ++i)
	    {

	      
	      
	      // Update velocity
	      particleManager->m_particleObjects[fluidIndex]->velocity[i] += fluid->acceleration[i] * _timestep; 
	  
	      // Update postion based on delta t
	      particleManager->m_particleObjects[fluidIndex]->position[i] +=  particleManager->m_particleObjects[fluidIndex]->velocity[i] * _timestep;
	      

	    }
	}
}


// Utility, computes the colors for each particle
void SmoothedParticleHydrodynamics::compute_color(){

  for(int bb = 1; bb < particleManager->m_particleObjects.size(); bb++){
    
    for(int i = 0; i < particleManager->m_particleObjects[bb]->color.size(); i++){
      particleManager->m_particleObjects[bb]->color[i] = Eigen::Vector3d(0.0, 0.0, 0.0);
    }
	  
  }

  for(int fI = 0; fI < particleManager->fluidIndicies.size(); fI++){
	unsigned int fluidIndex = particleManager->fluidIndicies[fI];

	Fluid* fluid = getFluidObject(fluidIndex);
	
	for (int pid = 0; pid < getObjectSize(fluidIndex); ++pid){
	  particleManager->m_particleObjects[fluidIndex]->color[pid] = Eigen::Vector3d(0.0, 0.0, 1.0) * fluid->density[pid] / ( 2.0 * _parameters->restDensity);
	}
  }


  
  
  if( _parameters->observeParticle >= 0 ){
    int pid = _parameters->observeParticle;
    
    particleManager->m_particleObjects[0]->color[pid] = Eigen::Vector3d(1.0, 0.0, 0.0);

    CompactNSearch::PointSet const& ps = m_neighborhoodsearch.point_set(0);

    for (int j = 0; j < ps.n_neighbors(pid); j++) {
	  
      CompactNSearch::PointID nid = ps.neighbor(pid, j);

      Eigen::Vector3d q = getPosition(0, pid) - getPosition(nid.point_set_id, nid.point_id);


      if(q.norm() <= _parameters->h)
	particleManager->m_particleObjects[nid.point_set_id]->color[nid.point_id] = Eigen::Vector3d(1.0, 0.5, 0.0);
      
    }

    
    return;
  }
    
  


  


  
  /*
  CompactNSearch::PointSet const& ps = m_neighborhoodsearch.point_set(0);

  double color = 0;
 
  for (int j = 0; j < ps.n_neighbors(pid); j++) {
    int neighbor_id = ps.neighbor(pid, j).point_id;

    if ((*m_particles_attribute)[neighbor_id] == 0)
      continue;

    Eigen::Vector3d vec_ij = deltaVector(pid, j);
    double q = m_kernel->w(vec_ij, _parameters->h);

    if(m_particles_density[j] != 0)
      color += q * _parameters->particleMass / m_particles_density[j];

  }

  m_color[pid] = Eigen::Vector3d(0,0,255) * color;
  */
}

void SmoothedParticleHydrodynamics::initParticles(){

  // calc mass of particles
  // based on particle radius

  double diameter = 2.0 * _parameters->particleRadius;
  _parameters->particleMass = _parameters->particleMassScaling * pow(diameter, 3) * _parameters->restDensity;

  //  _parameters->particleMass = 1;
  
  // reposition particle at start

  for(int fI = 0; fI < particleManager->fluidIndicies.size(); fI++){
    unsigned int fluidIndex = particleManager->fluidIndicies[fI];
  
    for(int i = 0; i < getObjectSize(fluidIndex); i++){
      (*m_particleObjects)[fluidIndex]->position[i] = (*m_particleObjects)[fluidIndex]->position0[i];
      (*m_particleObjects)[fluidIndex]->velocity[i] = Eigen::Vector3d(0.0, 0.0, 0.0);
      (*m_particleObjects)[fluidIndex]->color[i] = Eigen::Vector3d(0.0, 0.0, 0.0);	  
    }
  }
  /*
  for(int j = 1; j < m_particleObjects->size(); j++){
    if((*m_particleObjects)[j]->isDynamic){
      std::cout << "sii" << (*m_particleObjects)[j]->position0.size() << std::endl;
      for(int i = 0; i < (*m_particleObjects)[j]->position0.size(); i++){
	(*m_particleObjects)[j]->position[i] = (*m_particleObjects)[j]->position0[i];
      }
    }
  }
  */
}


void SmoothedParticleHydrodynamics::updateTimeStepCFL(){
  double stepSize = _parameters->stepSize;

  double maxVelo = 0.1;

  double diameter = 2.0 * _parameters->particleRadius;

  for(int fI = 0; fI < particleManager->fluidIndicies.size(); fI++){
    unsigned int fluidIndex = particleManager->fluidIndicies[fI];
  
    for(int i = 0; i < getObjectSize(fluidIndex); i++){
      const Eigen::Vector3d a = getFluidObject(fluidIndex)->acceleration[i];
      const Eigen::Vector3d v = getFluidObject(fluidIndex)->velocity[i];
    
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

void SmoothedParticleHydrodynamics::compute_stats(){
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

	Fluid* fluid = getFluidObject(fluidIndex);

	fluidSize += getObjectSize(fluidIndex);
	
	CompactNSearch::PointSet const &ps = m_neighborhoodsearch.point_set(fluidIndex);
	
	for (int i = 0; i < getObjectSize(fluidIndex); ++i){

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
      _parameters->avg_density = avg_density / fluidSize;
      _parameters->max_density = max_density;


      _parameters->min_velo = min_velo;
      _parameters->avg_velo = avg_velo / fluidSize;
      _parameters->max_velo = max_velo;

      _parameters->avg_neighbors = avg_neighbors / fluidSize;
      _parameters->avg_fluid_neighbors = avg_fluid_neighbors / fluidSize;
      
}


bool SmoothedParticleHydrodynamics::debugP(int id){
  if( id == _parameters->observeParticle)
    return true;

  return false;
}
