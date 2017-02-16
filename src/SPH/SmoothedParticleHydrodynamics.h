#ifndef SMOOTHEDPARTICLEHYDRODYNAMICS_H
#define SMOOTHEDPARTICLEHYDRODYNAMICS_H
#include <Eigen/Core>

#include <CompactNSearch>
#include <omp.h>
#include <array>
#include "../kernel/Kernel.h"
#include "../kernel/KernelDensity.h"
#include "../Parameters.h"
#include "../ParticleManager.h"


class SmoothedParticleHydrodynamics
{
public:
	// initialize new SPH system
	SmoothedParticleHydrodynamics(ParticleManager *manager, Parameters *params);
	
	std::vector<Particles*> *m_particleObjects;
	
	// update particle positions using semi-implicit euler method
	void update_semi_implicit_euler(double _timestep);	
	void nhSearch();
	void initParticles();
protected:
	Kernel *m_kernel;
	Parameters* _parameters;
	ParticleManager* particleManager;

	void compute_color();
	void compute_stats();

	/*
	std::vector<Eigen::Vector3d> m_particles_velocity_viscosity_double;

	std::vector<double> m_particles_density;
	std::vector<double> m_particles_density_fac;
	std::vector<double> m_particles_pressure_force;
	
	std::vector<Eigen::Vector3d> m_particles_acceleration;
	 */
	
	Eigen::Vector3d &getPosition(unsigned int point_set, unsigned int i);

	std::vector<int>* m_particles_attribute;
	
	bool isFluid(unsigned int point_set_id);
	
	std::vector<double> m_boundary_psi;
	
	CompactNSearch::NeighborhoodSearch m_neighborhoodsearch;

	int m_particles_size;

	bool debugP(int id);
       
	Boundary* getBoundaryObject(unsigned int i);
	Fluid* getFluidObject(unsigned int i);

	unsigned int getObjectSize(unsigned int i);
	
	void updateTimeStepCFL();
	
private:



	bool use_visity = true;



	  
};




#endif
