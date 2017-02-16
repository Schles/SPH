// This file is part of GLviz.
//
// Copyright(c) 2014, 2015 Sebastian Lipponer
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files(the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and / or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions :
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
// IN THE SOFTWARE.

#include <GLviz>

#include "config.hpp"

#include <Eigen/Core>

#include <iostream>
#include <memory>
#include <fstream>
#include <sstream>
#include <vector>
#include <array>
#include <exception>
#include <cmath>
#include <ctime>

#include "../CompactNSearch/CompactNSearch.h"

#ifdef WIN32
#include "Windows.h"
#endif

using namespace Eigen;

namespace
{

    GLviz::Camera camera;

    float g_time(0.0f);
    bool g_stop_simulation(true);

    bool g_enable_mesh3(true);
    bool g_enable_wireframe(false);

    bool g_enable_points(false);
    float g_point_radius(0.0014f);
    float g_projection_radius(0.0f);

    float g_wireframe[4] = {
        0.0f, 0.0f, 0.0f, 1.0f
    };

    float g_mesh_material[4] = {
        0.0f, 0.25f, 1.0f, 8.0f
    };

    float g_points_material[4] = {
        1.0f, 1.0f, 1.0f, 8.0f
    };

    int g_shading_method(0);

    // returns the time passed since t0 (in seconds)
    double createTimeStamp(double t0)
    {
#ifdef WIN32
        LONGLONG currentCount;
        QueryPerformanceCounter((LARGE_INTEGER*)&currentCount);

        LONGLONG frequency;
        QueryPerformanceFrequency((LARGE_INTEGER*)&frequency);

        return (double)currentCount / frequency - t0;
#endif

    }

    class BruteForceSearch
    {
    public:
        BruteForceSearch(Real r)
        {
            neighborhood_size = r;
        }

        void add_discretization(Eigen::Vector3d* x_, std::size_t n_)
        {
            points_count = n_;
            points_neighbor_count.resize(n_);
            points = x_;
        }

        void neighborhood_search()
        {
            for (int i = 0; i < points_count; i++)
            {
                points_neighbor_count[i] = 0;
                for (int j = 0; j < points_count; j++)
                {
                    if (points_in_neighborhood(i, j))
                        points_neighbor_count[i]++;
                }
            }
        }

        unsigned int n_neighbors(unsigned int i) const
        {
            if (i < points_count)
            {
                return points_neighbor_count[i];;
            }
            else
            {
                return 0;
            }
        }


    private:
        bool points_in_neighborhood(int id0_, int id1_)
        {
            if (id0_ == id1_)
                return false;

            Eigen::Vector3d distance_vector = points[id0_] - points[id1_];

            double distance = sqrt(pow(distance_vector[0], 2) + pow(distance_vector[1], 2) + pow(distance_vector[2], 2));

            if (distance < neighborhood_size)
            {
                return true;
            }
            else
            {
                return false;
            }
        }

        std::size_t points_count = 0;
        Eigen::Vector3d* points;
        std::vector<std::size_t> points_neighbor_count;
        Real neighborhood_size;
    };

    struct MyViz
    {
        MyViz(unsigned int particleCount_, unsigned int boundaryBoxSize_)
        {
            particle_count = particleCount_;

            srand(std::time(NULL));

            // Create randomly distributed point cloud
            particle_positions = (Vector3d*)malloc(sizeof(Vector3d) * particle_count);
            for (unsigned int i = 0; i < particle_count; ++i)
            {
                particle_positions[i][0] = rand() % boundaryBoxSize_;
                particle_positions[i][1] = rand() % boundaryBoxSize_;
                particle_positions[i][2] = rand() % boundaryBoxSize_;
            }

            // Compute neighborhood information
            CompactNSearch search = CompactNSearch(neighborhood_radius);
            BruteForceSearch brute = BruteForceSearch(neighborhood_radius);
            search.add_discretization(particle_positions, particle_count, false, true);
            brute.add_discretization(particle_positions, particle_count);

            printf("\nBenchmark for n=%i, boundary=%i\n", particle_count, boundaryBoxSize_);

            double t0 = createTimeStamp(0);
            search.neighborhood_search();
            double t1 = createTimeStamp(t0);
            printf("CompactNSearch: %f s\n", t1);

            double t2 = createTimeStamp(0);
            brute.neighborhood_search();
            double t3 = createTimeStamp(t0);
            printf("BruteForceSearch: %f s\n", t3);

            // Check results
            std::vector<SPHDiscretization> neighborhood;
            neighborhood = search.discretizations();

            for (int i = 0; i < particle_count; i++)
            {
                int countCompact = neighborhood[0].n_neighbors(i);
                int countBruteForce = brute.n_neighbors(i);
                if (countCompact != countBruteForce)
                {
                    printf("Error: %i != %i\n", countCompact, countBruteForce);
                }
            }
        }

        void FreeViz()
        {
            free(particle_positions);
        }

        unsigned int particle_count = 5000;
        Real neighborhood_radius = 20.0;
        Vector3d* particle_positions;

    };

    std::unique_ptr<MyViz> viz;

    void
        displayFunc()
    {
        glEnable(GL_MULTISAMPLE);
        glEnable(GL_DEPTH_TEST);

        glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
        glClearDepth(1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);


        if (g_enable_mesh3)
        {

        }

        if (g_enable_points)
        {

        }
    }

    void
        reshapeFunc(int width, int height)
    {
        const float aspect = static_cast<float>(width) /
            static_cast<float>(height);

        glViewport(0, 0, width, height);
        camera.set_perspective(60.0f, aspect, 0.005f, 5.0f);
    }

    void
        closeFunc()
    {
        viz->FreeViz();
        viz = nullptr;
    }

    void
        timerFunc(int delta_t_msec)
    {
        float delta_t_sec = static_cast<float>(delta_t_msec) / 1000.0f;

        if (!g_stop_simulation)
        {
            g_time += delta_t_sec;

            const float k = 50.0f;
            const float a = 0.03f;
            const float v = 10.0f;


        }
    }

    void TW_CALL reset_simulation(void*)
    {
        g_time = 0.0f;
    }

    void TW_CALL
        export_frame(void*)
    {
        GLviz::export_frame("export");
    }

}

int
main(int argc, char* argv[])
{
    GLviz::init(argc, argv);

    viz = std::unique_ptr<MyViz>(new MyViz(1000, 100));
    viz->FreeViz();

    viz = std::unique_ptr<MyViz>(new MyViz(2000, 100));
    viz->FreeViz();

    viz = std::unique_ptr<MyViz>(new MyViz(3000, 100));
    viz->FreeViz();

    viz = std::unique_ptr<MyViz>(new MyViz(4000, 100));
    viz->FreeViz();

    viz = std::unique_ptr<MyViz>(new MyViz(5000, 100));
    viz->FreeViz();

    viz = std::unique_ptr<MyViz>(new MyViz(1000, 100));
    viz->FreeViz();

    viz = std::unique_ptr<MyViz>(new MyViz(1000, 200));
    viz->FreeViz();

    viz = std::unique_ptr<MyViz>(new MyViz(1000, 300));
    viz->FreeViz();

    viz = std::unique_ptr<MyViz>(new MyViz(1000, 400));
    viz->FreeViz();

    viz = std::unique_ptr<MyViz>(new MyViz(1000, 500));
    viz->FreeViz();

    // Setup AntTweakBar.
    {
        TwBar* bar = GLviz::twbar();

        TwAddVarRO(bar, "time", TW_TYPE_FLOAT,
            &g_time, " precision=3 label='t in sec' group='Simulation' ");

        TwAddButton(bar, "Reset",
            reset_simulation, NULL,
            " key=r help='Reset simulation' group='Simulation' ");

        TwAddVarRW(bar, "Stop", TW_TYPE_BOOLCPP,
            &g_stop_simulation,
            " key=SPACE help='Stop simulation' group='Simulation' ");

        TwAddVarRW(bar, "Draw Triangle Mesh", TW_TYPE_BOOLCPP,
            &g_enable_mesh3,
            " key=1 help='Draw Triangle Mesh' group='Triangle Mesh' ");

        TwType shading_type = TwDefineEnumFromString("shading_type", "Flat,Phong");
        TwAddVarRW(bar, "Shading", shading_type, &g_shading_method, " key=5 group='Triangle Mesh' ");

        TwAddSeparator(bar, NULL, " group='Triangle Mesh' ");

        TwAddVarRW(bar, "Wireframe", TW_TYPE_BOOLCPP,
            &g_enable_wireframe,
            " key=w help='Draw Wireframe' group='Triangle Mesh' ");

        TwAddVarRW(bar, "Wireframe Color", TW_TYPE_COLOR3F,
            g_wireframe,
            " help='Wireframe Color' group='Triangle Mesh' ");

        TwAddSeparator(bar, NULL, " group='Triangle Mesh' ");

        TwAddVarRW(bar, "Mesh Material", TW_TYPE_COLOR3F,
            g_mesh_material,
            " help='Mesh Ambient' group='Triangle Mesh' ");

        TwAddVarRW(bar, "Mesh Shininess", TW_TYPE_FLOAT,
            &g_mesh_material[3],
            " min=1e-12 max=1000 help='Mesh Shininess' group='Triangle Mesh' ");

        TwAddVarRW(bar, "Draw Points", TW_TYPE_BOOLCPP,
            &g_enable_points,
            " key=2 help='Draw Points' group='Points' ");

        TwAddVarRW(bar, "Radius", TW_TYPE_FLOAT,
            &g_point_radius,
            " min=0 max=0.1 step=0.0001 key=2 help='Radius' group='Points' ");

        TwAddVarRW(bar, "Points Material", TW_TYPE_COLOR3F,
            g_points_material,
            " help='Points Ambient' group='Points' ");

        TwAddVarRW(bar, "Points Shininess", TW_TYPE_FLOAT,
            &g_points_material[3],
            " min=1e-12 max=1000 help='Points Shininess' group='Points' ");

        TwAddButton(bar, "Export Frame", export_frame, nullptr, nullptr);
    }

    GLviz::display_callback(displayFunc);
    GLviz::reshape_callback(reshapeFunc);
    GLviz::close_callback(closeFunc);
    GLviz::timer_callback(timerFunc, 15);

    return GLviz::exec(camera);
}
