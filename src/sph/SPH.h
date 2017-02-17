#ifndef SMOOTHEDPARTICLEHYDRODYNAMICS_H
#define SMOOTHEDPARTICLEHYDRODYNAMICS_H
#include <Eigen/Core>

#include <CompactNSearch>
#include <omp.h>
#include <array>

#include "../Statistics.h"
#include "../kernel/Kernel.h"
#include "../kernel/KernelDensity.h"
#include "../Parameters.h"
#include "../particle/ParticleManager.h"


class SPH
{
public:
	// initialize new SPH system
	SPH(ParticleManager *manager, Parameters *params);

	Statistics* m_statistics;
	
	// update particle positions using semi-implicit euler method
	void compute_semi_implicit_euler(int fluidIndex, int particle_id, double _timestep);	
	void nhSearch();
	void initParticles();

	void addParticleSet(int fluidIndex);
	
	void updateTimeStepCFL();
protected:
	Kernel *m_kernel;
	Parameters* _parameters;
	ParticleManager* particleManager;
	
	CompactNSearch::NeighborhoodSearch m_neighborhoodsearch;
	
	void compute_color(int fluidIndex, int particleId);
	void compute_stats();

	bool debugP(int id);			  
};




#endif
