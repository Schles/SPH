#ifndef PARTICLE_FACTORY_H
#define PARTICLE_FACTORY_H
#include <iostream>
#include <vector>

#include <Eigen/Core>

#include "Particles.h"
#include "../Parameters.h"

#include "ParticleManager.h"

class ParticleFactory {
public:
        ParticleFactory(Parameters* params);
	
        Fluid* createFluidParticles(Eigen::Vector3d start, Eigen::Vector3d end, Eigen::Vector3d density);
        Boundary* createBoundingBox(Eigen::Vector3d start, Eigen::Vector3d end, double h, bool isDynamic);
        Boundary* createBoundingWall(Eigen::Vector3d p, Eigen::Vector3d u, Eigen::Vector3d v, double h);	
        Boundary* createBoundingBall(Eigen::Vector3d start,  double h);

	Fluid* createFluidTranslate(Eigen::Vector3d vec, double size, double density);

	Fluid* createFluidObject(unsigned int size);
	Boundary* createBoundingObject(unsigned int size);


	Parameters *m_Params;
};


#endif
