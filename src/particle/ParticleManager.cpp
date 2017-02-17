#include "ParticleManager.h"
#include <math.h>

ParticleManager::ParticleManager(Parameters *params){
    m_Params = params;   
}


void ParticleManager::updateDynamicBoundary(double time){
  // dirty for boundries checking

  double amp = 0.05;

  Eigen::Vector3d deltaX(-1.0, 0.0, 0.0);

  double radian = 3.14 / 180.0;

  double input = time * 50;
  
  double y = 1 - cos(input);

  deltaX *= amp *  y;
  
  for (int i = 2; i < m_particleObjects.size(); ++i){
      for(int j = 0; j < m_particleObjects[i]->position.size(); j++){
	m_particleObjects[i]->position[j] = m_particleObjects[i]->position0[j] + deltaX;
      }
  }
}

void ParticleManager::castFluidPositions(unsigned int point_set){
  
  for (int i = 0; i < m_particleObjects[point_set]->position.size(); ++i)
	{
	  // Set output postion
	  m_fluid_position_output[point_set][i][0] = (float) m_particleObjects[point_set]->position[i][0];
	  m_fluid_position_output[point_set][i][1] = (float) m_particleObjects[point_set]->position[i][1];
	  m_fluid_position_output[point_set][i][2] = (float) m_particleObjects[point_set]->position[i][2];
	  
	  // Set Output colors
	  m_fluid_color_output[point_set][i][0] = (float) m_particleObjects[point_set]->color[i][0];
	  m_fluid_color_output[point_set][i][1] = (float) m_particleObjects[point_set]->color[i][1];
	  m_fluid_color_output[point_set][i][2] = (float) m_particleObjects[point_set]->color[i][2];

	}

}



void ParticleManager::resizeOutputBuffer(unsigned int bufferId, unsigned int size){
  unsigned int bsize = bufferId; //m_fluid_position_output.size();

  m_fluid_position_output.resize(bsize + 1);
  m_fluid_position_output[bsize].resize(size);

  m_fluid_color_output.resize(bsize + 1);
  m_fluid_color_output[bsize].resize(size);
}




unsigned int ParticleManager::addParticles(Particles* particles){

  unsigned int index = m_particleObjects.size();

  m_particleObjects.push_back(particles);
  
  unsigned int size = particles->position.size();

  resizeOutputBuffer(index, size);

  return index;
}


Eigen::Vector3d &ParticleManager::getPosition(unsigned int point_set, unsigned int i){
  return (m_particleObjects)[point_set]->position[i];
}

Particles* ParticleManager::getObject(unsigned int i){
  return m_particleObjects[i];
}

Boundary* ParticleManager::getBoundaryObject(unsigned int i){
  return static_cast<Boundary*>((m_particleObjects)[i]);
}

Fluid* ParticleManager::getFluidObject(unsigned int i){
  return static_cast<Fluid*>((m_particleObjects)[i]);
}

unsigned int ParticleManager::getObjectSize(unsigned int i){
  return (m_particleObjects)[i]->position.size();
}

unsigned int ParticleManager::getParticleGroupSize(){
  return m_particleObjects.size();
}


bool ParticleManager::isFluid(unsigned int point_set_id){
  for(int i = 0; i < fluidIndicies.size(); i++){
    if(point_set_id == fluidIndicies[i])
      return true;
  }

  return false;
}
