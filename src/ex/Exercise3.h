#ifndef _EXERCISE3_H_
#define _EXERCISE3_H_

#include "../ParticleManager.h"
#include "../Parameters.h"

class Exercise3 {
 public:
  Exercise3(ParticleManager* manager, Parameters *params);
  double getNaiveDensity(int particle_id);
  double getNaiveBBForce(int particle_id);
  
 private:
  ParticleManager* m_manager;
  Parameters* m_params;
};


#endif /* _EXERCISE3_H_ */
