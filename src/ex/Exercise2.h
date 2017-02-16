#ifndef EXERCISE_2_H
#define EXERCISE_2_H

#include <Eigen/Core>

#include "Kernel.h"
#include "KernelDensity.h"


#include "../util/fluidtime.h"

class Exercise2{
    public:
        Exercise2(double h, unsigned int particleCount, double particleMass, double compactConstant);
        void assignment1();
    private:
        double h;
        unsigned int particleCount;
        double particleMass;
        double compactConstant;
};


#endif
