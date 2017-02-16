#include "KernelDensity.h"
#include <iostream>
#include <math.h>

KernelDensity::KernelDensity(){
    particleMass = 1;

}

KernelDensity::KernelDensity(double _particleMass, double _compactConstant){
    particleMass = _particleMass;
    compactConstant = _compactConstant;
    kernel = new Kernel(3, 0.1);

}


void KernelDensity::initParticleField(Eigen::Vector3d start, Eigen::Vector3d end, double particleCount){
    double x,y,z;
    double dir_x, dir_y, dir_z;

    dir_x = (end[0] - start[0]) / (particleCount - 1);
    dir_y = (end[1] - start[1]) / (particleCount - 1);
    dir_z = (end[2] - start[2]) / (particleCount - 1);

    for(int i_x = 0; i_x < particleCount; i_x++){
        x = start[0] + dir_x * i_x;

        for(int i_y = 0; i_y < particleCount; i_y++){
            y = start[1] + dir_y * i_y;

            for(int i_z = 0; i_z < particleCount; i_z++){
                z = start[2] + dir_z * i_z;

                positions.push_back({x,y,z});

            }
        }
    }
}

void KernelDensity::addSampleLine(unsigned int samplePoints){
    double dir_sl = 4.0 / (samplePoints - 1);

    for(int i = 0; i < samplePoints; i++){
        double x = -2 + dir_sl * i;

        positions.push_back({x,0,0});
        sampleLine.push_back(positions.size() - 1);
    }
}

void KernelDensity::neighborhood_search(double neighborhoodRadius){
    compactNSearch = new CompactNSearch::NeighborhoodSearch(neighborhoodRadius);
    point_set_id = compactNSearch->add_point_set(positions.front().data(), positions.size());
    compactNSearch->find_neighbors();
}


double KernelDensity::skalarQuantity(unsigned int i, double h){

    double res = 0;

    Eigen::Vector3d deltaVector;
    Eigen::Vector3d r = Eigen::Vector3d(positions[i][0], positions[i][1], positions[i][2]);

    CompactNSearch::PointSet const& ps = compactNSearch->point_set(point_set_id);

    for (int j = 0; j < ps.n_neighbors(i); j++) {

        CompactNSearch::PointID const& pid = ps.neighbor(i, j);

        Eigen::Vector3d r_j = Eigen::Vector3d(positions[pid.point_id][0], positions[pid.point_id][1], positions[pid.point_id][2]);

        deltaVector = r - r_j;

        double w = kernel->w(deltaVector, h);

        //std::cout << "dX: " << r_j[0] << ", w: " << w <<  " n: " << ps.n_neighbors(pid.point_id) << std::endl;

        double vLen = sqrt( deltaVector.dot(deltaVector) );

        // Compact Condition / Support

        if(vLen <= compactConstant * h){
            res += (particleMass / density(j, h)) * ps.n_neighbors(pid.point_id) * w;
        }

    }

    return res;
}

double KernelDensity::computeDensityCoeffient(Eigen::Vector3d x, double h){

    double norm = x.norm();

	double q = norm / h;

    std::cout << "q " << q << std::endl;
    std::cout << "norm " << norm << std::endl;

	if (0 < q && q < 2.0 / 3.0)
	{
		q = 2.0 / 3.0;
	}
	else if (q < 1.0)
	{
		q = (2.0 * q - 3.0 / 2.0 * pow(q, 2));
	}
	else if (q < 2.0)
	{
		q = 1.0 / 2.0 * pow(2.0 - q, 2);
	}
	else
	{
		q = 0;
	}

    std::cout << q << std::endl;

    return q / norm;

}

double KernelDensity::density(unsigned int i, double h){
    double res = 0;

    Eigen::Vector3d deltaVector;
    Eigen::Vector3d r = Eigen::Vector3d(positions[i][0], positions[i][1], positions[i][2]);

    CompactNSearch::PointSet const& ps = compactNSearch->point_set(point_set_id);

    for (int j = 0; j < ps.n_neighbors(i); j++) {

        CompactNSearch::PointID const& pid = ps.neighbor(i, j);

        Eigen::Vector3d r_j = Eigen::Vector3d(positions[pid.point_id][0], positions[pid.point_id][1], positions[pid.point_id][2]);

        deltaVector = r - r_j;

        //double vLen = sqrt( deltaVector.dot(deltaVector) );
        double vLen = deltaVector.norm();

        // Compact Condition / Support
        if(vLen <= compactConstant * h){
        }
        res += particleMass * kernel->w(deltaVector, h);

    }

    return res;
}
