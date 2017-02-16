#include "ParticleExporter.h"

#include <fstream>

void ParticleExporter::exportParticles(std::string _path, std::vector<Eigen::Vector3f> _particles, int _count)
{
    float* data = (float*)malloc(_count * 3 * sizeof(float));

    for (int i = 0; i < _count; i++)
    {
        data[3 * i + 0] = _particles[i][0];
        data[3 * i + 1] = _particles[i][1];
        data[3 * i + 2] = _particles[i][2];
    }

    std::ofstream out(_path, std::ios_base::binary);

    if (out.good())
    {
        out.write((char*)data, _count * 3 * sizeof(float));
        out.close();
    }
    else {
        printf("Could not open output file.");
    }

    free(data);
}
