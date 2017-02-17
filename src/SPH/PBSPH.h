#ifndef PD_SPH_H
#define PD_SPH_H


#include "SmoothedParticleHydrodynamics.h"

#include "BoundaryPsi.h"

#include "../util/fluidtime.h"

class PBSPH : protected SmoothedParticleHydrodynamics {

 public:
  PBSPH(ParticleManager* particleManager, Parameters *params);
  
 protected:
  BoundaryPsi* boundaryPsi;
  
  void compute_position_update(int fluidIndex, int particle);
  void compute_prev_position_update(int fluidIndex, int particle);
  void compute_velocity(int fluidIndex, int pId, double timeStep);
  
  void compute_friction(int fluidIndex, int particle);
  void compute_external_forces(int fluidIndex, int particle);
  
  void compute_density(int fluidIndex, int particle);
  void compute_delta_x(int fluidIndex, int particleId);
  void compute_lambda(int fluidIndex, int particleId);
  void compute_viscosity(int fluidIndex, int pid); 
};


#endif
