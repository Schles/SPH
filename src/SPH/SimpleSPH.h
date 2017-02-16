#ifndef SIMPLE_SPH_H
#define SIMPLE_SPH_H

#include "SmoothedParticleHydrodynamics.h"

class SimpleSPH : protected SmoothedParticleHydrodynamics {

public:
  SimpleSPH(ParticleManager* particleManager, Parameters *params);
  void step(double deltaTime);

protected:
  void compute_viscosity(int fluidIndex, int _particle_id);  
  
private:
  void compute_density(int _particle_id);
  void compute_acceleration(int _particle_id);
  Eigen::Vector3d compute_boundary_force(int _i, int _j);
  
};

#endif
