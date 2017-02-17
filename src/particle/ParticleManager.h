#ifndef PARTICLE_MANAGER_H
#define PARTICLE_MANAGER_H
#include <iostream>
#include <vector>

#include <Eigen/Core>

#include "Particles.h"
#include "../Parameters.h"

class ParticleManager {
public:
        ParticleManager(Parameters* params);
	
	void updateFluidPositions(unsigned int point_set);

	void updateDynamicBoundary(double time);

	unsigned int addParticles(Particles* particles);

	void resizeOutputBuffer(unsigned int bufferId, unsigned int size);
	


	Parameters *m_Params;
	
	std::vector<Particles*> m_particleObjects;
	std::vector<int> fluidIndicies;
	
	std::vector<std::vector<Eigen::Vector3f>> m_fluid_position_output;
	std::vector<std::vector<Eigen::Vector3f>> m_fluid_color_output;

	bool isFluid(unsigned int point_set_id);

	Particles* getObject(unsigned int i);
	Boundary* getBoundaryObject(unsigned int i);
	Fluid* getFluidObject(unsigned int i);

	unsigned int getParticleGroupSize();
	unsigned int getObjectSize(unsigned int i);
	Eigen::Vector3d &getPosition(unsigned int point_set, unsigned int i);
	



	


};


#endif
