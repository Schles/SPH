#ifndef KERNEL_DENSITY_H
#define KERNEL_DENSITY_H

#include <Eigen/Core>
#include "Kernel.h"
#include <CompactNSearch>
#include <omp.h>
#include <array>
/*
    Kernelfunction numbers:
    1: Poly6
    2: Spiky
    3: cubic spline
 */

class KernelDensity {
    public:
        KernelDensity();
        KernelDensity(double _particleMass, double _compactConstant);
        void initParticleField(Eigen::Vector3d start, Eigen::Vector3d end, double particleCount);
        void addSampleLine(unsigned int samplePoints);
        void neighborhood_search(double neighborhoodRadius);
        double skalarQuantity(unsigned int i, double h);
        double density(unsigned int i, double h);

        double computeDensityCoeffient(Eigen::Vector3d x, double h);

        std::vector<std::array<double, 3>> positions;
        std::vector<unsigned int> sampleLine; // holds index of sample line points



    private:
        Kernel *kernel;
        CompactNSearch::NeighborhoodSearch *compactNSearch;

        unsigned int point_set_id;
        double compactConstant = 1.0;
        double particleMass = 1.0;
	
	
        unsigned int w_func = 3;



};


#endif
