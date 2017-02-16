#include "ParticleManager.h"
#include <math.h>

ParticleManager::ParticleManager(Parameters *params){
  //    particlesBB = 0;
    //    particlesFluid = 0;
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

void ParticleManager::updateFluidPositions(unsigned int point_set){



  
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


void ParticleManager::initParticleBuffer(unsigned int size){
  std::cout << "!!!!!!!!!!!!!!!!!!!! " << std::endl;
  Particles* particles = new Particles();
  
  particles->position.resize(size);
  particles->position0.resize(size);
  particles->velocity.resize(size);
  particles->color.resize(size);
  
   
  m_particleObjects.push_back(particles);


  // Init output buffer
  resizeOutputBuffer(0, size);
}

unsigned int ParticleManager::addParticles(Particles* particles){

  unsigned int index = m_particleObjects.size();

  m_particleObjects.push_back(particles);
  
  unsigned int size = particles->position.size();

  resizeOutputBuffer(index, size);


  return index;
}

void ParticleManager::simpleTest(){

  Boundary* particles = new Boundary();
  
  int memAllocs = 4;
  
  //    particlesBB += memAllocs;

    particles->position.resize(memAllocs);
    particles->color.resize(memAllocs);
    particles->m_boundaryPsi.resize(memAllocs);

    resizeOutputBuffer(1, memAllocs);

    int memOffset = 0;

    Eigen::Vector3d boundryColor(0.0,0.0,0.0);

    // Z Planes

    particles->position[0] = Eigen::Vector3d( 0.02, 0.0, 0.02 );
    particles->color[0] = boundryColor;

    particles->position[1] = Eigen::Vector3d( 0.02, 0.0, -0.02 );
    particles->color[1] = boundryColor;

    particles->position[2] = Eigen::Vector3d( -0.02, 0.0, 0.02 );
    particles->color[2] = boundryColor;

    particles->position[3] = Eigen::Vector3d( -0.02, 0.0, -0.02 );
    particles->color[3] = boundryColor;

	    
    m_particleObjects.push_back(particles);

}

