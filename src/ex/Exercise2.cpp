#include "Exercise2.h"
#include <iostream>
#include <fstream>
#include "../util/fluidtime.h"


Exercise2::Exercise2(double _h, unsigned int _particleCount, double _particleMass, double _compactConstant){
    h = _h;
    particleCount = _particleCount;
    particleMass = _particleMass;
    compactConstant = _compactConstant;

    std::cout << "h: " << h << std::endl;
    std::cout << "mass: " << particleMass << std::endl;
    std::cout << "kappa: " << compactConstant << std::endl;
    std::cout << "particlePerAxis: " << particleCount << std::endl;
    std::cout << "#particles: " << particleCount * particleCount * particleCount << std::endl;

}

void Exercise2::assignment1(){

        std::ofstream myfile;

        double samplePoints = 3 * particleCount;

        long t1 = Time::getMilliseconds();

        std::cout << "creating particlefield..." << std::endl;

        KernelDensity kernelDensity(particleMass, compactConstant);

        // create particle field
        kernelDensity.initParticleField(
                Eigen::Vector3d(-1.0, -1.0, -1.0),
                Eigen::Vector3d(1.0, 1.0, 1.0),
                particleCount);

        long t2 = Time::getMilliseconds();
        std::cout << t2 - t1 << "ms" << std::endl;

        // 1c)
        kernelDensity.addSampleLine(samplePoints);

        std::cout << "neighborhood search..." << std::endl;

        // calculate neighbors
        kernelDensity.neighborhood_search( 2.0 * h);

        long t3 = Time::getMilliseconds();

        std::cout << t3 - t2 << "ms" << std::endl;

        // calc density and write to file
        myfile.open ("../plot/density.dat");

        Eigen::Vector3d ref = Eigen::Vector3d(0.0,0.0,0.0);


        for(int i = 0; i < kernelDensity.sampleLine.size(); i++){
            unsigned int pid = kernelDensity.sampleLine[i];

            //double res = kernelDensity.skalarQuantity(pid, h);
            double res = kernelDensity.density(pid, h);

            myfile << kernelDensity.positions[pid][0] << " " << res << std::endl;
        }


        myfile.close();
}
