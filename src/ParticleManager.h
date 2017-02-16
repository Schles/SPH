#ifndef PARTICLE_MANAGER_H
#define PARTICLE_MANAGER_H
#include <iostream>
#include <vector>

#include <Eigen/Core>

#include "Particles.h"
#include "Parameters.h"

class ParticleManager {
public:
        ParticleManager(Parameters* params);
	/*
        std::vector<Eigen::Vector3d>			   m_vertices_fluid_particles;
        std::vector<int>                                   m_vertices_attribute;
        std::vector<Eigen::Vector3d>			   m_materials;
        std::vector<Eigen::Vector3d>			   m_velocities;
	*/

	std::vector<Particles*> m_particleObjects;
		
	//        unsigned int particlesFluid;
        //unsigned int particlesBB;

	std::vector<int> fluidIndicies;
	
	void updateFluidPositions(unsigned int point_set);

	void updateDynamicBoundary(double time);

	void resizeOutputBuffer(unsigned int bufferId, unsigned int size);
	
	void initParticleBuffer(unsigned int size);
	
	std::vector<std::vector<Eigen::Vector3f>> m_fluid_position_output;
	std::vector<std::vector<Eigen::Vector3f>> m_fluid_color_output;


	unsigned int addParticles(Particles* particles);
	
	Parameters *m_Params;
	
 private:
	void simpleTest();

};


#endif
