#ifndef PD_SPH_H
#define PD_SPH_H

#include "ImprSPH.h"
#include "SmoothedParticleHydrodynamics.h"

#include "../util/fluidtime.h"

class PBSPH : protected ImprSPH {

 public:
  PBSPH(ParticleManager* particleManager, Parameters *params);

  void init();
  
  void step(double deltaTime);

  std::vector<std::vector<Eigen::Vector3f>> m_vector_output;
  
 private:
  void updateVelocity(unsigned int set, unsigned int pId, double timeStep);
  void solvePressure();
  void friction(int fluidIndex, int particle);
  void externalForces(int fluidIndex, int particle);
  void compute_density(int fluidIndex, int particle);

  Eigen::Vector3d compute_delta_x(int fluidIndex, int particleId);

  std::vector<Eigen::Vector3d> m_delta_x_fluid;
  std::vector<Eigen::Vector3d> m_delta_x_bb;

  void updateDeltaXOutput();

  int visualizeAmount = 1;
  
};


#endif
