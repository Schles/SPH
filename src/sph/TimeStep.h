#ifndef TIMESTEP_PB_SPH_H
#define TIMESTEP_PB_SPH_H

#include "PBSPH.h"

class TimeStepPBSPH : protected PBSPH {

 public:
  TimeStepPBSPH(ParticleManager* particleManager, Parameters *params);

  void init();
  
  void step(double deltaTime);

  std::vector<std::vector<Eigen::Vector3f>> m_vector_output;
  
 private:
  void solvePressure();
  void updateDeltaXOutput();  
  int visualizeAmount = 1;
  void helperStep(double deltaTime);
  
};


#endif
