#ifndef SCENES_H
#define SCENES_H

#include <iostream>
#include <vector>

#include <Eigen/Core>

#include "Particles.h"
#include "Parameters.h"

#include "ParticleFactory.h"

class Scenes : public ParticleFactory {
public:
        Scenes(Parameters* params);
		
        ParticleManager* createSetup(unsigned int i, double h);

	void createBoundingBox(ParticleManager* manager, Eigen::Vector3d start, Eigen::Vector3d end, double h);
};


#endif
