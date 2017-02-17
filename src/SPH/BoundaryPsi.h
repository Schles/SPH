#ifndef BOUNDARY_PSI_H
#define BOUNDARY_PSI_H


#include "SmoothedParticleHydrodynamics.h"

#include "../util/fluidtime.h"

class BoundaryPsi {

 public:
  BoundaryPsi(CompactNSearch::NeighborhoodSearch *nhSearch, ParticleManager* particleManager, Parameters *params, Kernel *m_kernel);

  void calc_boundary_psi();
  
 private:

  ParticleManager *particleManager;
  Parameters* _parameters;
  CompactNSearch::NeighborhoodSearch *m_neighborhoodsearch;
  Kernel *m_kernel;
  
  void calc_boundary_psi_i(unsigned int pointSet);
  
  
};


#endif
