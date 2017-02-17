#pragma once

#include <string>
#include <vector>
#include <Eigen/Core>

class ParticleExporter {

public:
    
    static void exportParticles(std::string _path, std::vector<Eigen::Vector3f> _particles, int _count);
};