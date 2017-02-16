#ifndef IMPR_SPH_H
#define IMPR_SPH_H

#include "SimpleSPH.h"
#include "SmoothedParticleHydrodynamics.h"

class ImprSPH : protected SimpleSPH {

public:
  ImprSPH(ParticleManager* particleManager, Parameters *params);
  void step(double deltaTime);

 protected:
   void compute_density(int _particle_id);
   Eigen::Vector3d compute_boundary_force(int _i, int _j);
   double calc_boundary_psi();
private:
 
  void compute_acceleration(int _particle_id);
  void compute_viscosity(int _particle_id);

  void calc_boundary_psi_i(unsigned int pointSet);
  double compute_boundary_volume(int _particle_id);


  
};

#endif
