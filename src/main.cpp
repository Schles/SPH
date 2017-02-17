
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
#include <string>
#include <iomanip>
#include <cstdlib>

#include "Exercise2.h"

#include "Scenes.h"
#include "ParticleManager.h"
#include "ParticleExporter.h"
#include "Parameters.h"

#include "SPH/TimeStep.h"


using namespace Eigen;

namespace
{

GLviz::Camera camera;

float g_time(0.0f);
bool g_stop_simulation(true);
float g_point_radius(0.0025f);
float g_bb_radius(0.0009f);
float g_projection_radius(0.01f);
int g_shading_method(0);

int g_test_setup(2);

float g_points_material[4] = {
    1.0f, 1.0f, 1.0f, 8.0f
};


  bool g_enable_wireframe(false);

  float g_wireframe[4] = {
    0.0f, 0.0f, 0.0f, 1.0f
};

float g_mesh_material[4] = {
    0.0f, 0.25f, 1.0f, 8.0f
};

float g_displacement_color[3][4] = {
  {1.0f, 0.00f, 0.0f, 1.0f}, // Red, delta_x
  {0.0f, 1.00f, 0.0f, 1.0f}, // Green, delta_x_fluid
  {0.0f, 0.00f, 1.0f, 1.0f}  // Blue, delta_x_bb
};


struct MyViz
{
    MyViz()
    {
        // Setup vertex array v.
        vertex_array_v.bind();

        vertex_array_buffer.bind();
        glEnableVertexAttribArray(0);
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(GLfloat), reinterpret_cast<const GLvoid*>(0));

        vertex_color_buffer.bind();
        glEnableVertexAttribArray(1);
        glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(GLfloat), reinterpret_cast<const GLvoid*>(0));
	
        vertex_array_v.unbind();


       // Setup vertex array vf.
        vertex_array_v2.bind();

        vertex_displacement_buffer.bind();
        glEnableVertexAttribArray(0);
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(GLfloat), reinterpret_cast<const GLvoid*>(0));

        vertex_displacement_buffer.unbind();

	vertex_array_v2.unbind();

	
        // Bind uniforms to their binding points.
        uniform_camera.bind_buffer_base(0);
        uniform_material.bind_buffer_base(1);
        uniform_wireframe.bind_buffer_base(2);
        uniform_sphere.bind_buffer_base(3);

        camera.translate(Eigen::Vector3f(0.0f, 0.0f, -2.0f));
    }


  void draw_vector(GLfloat color[], GLsizei nf)
    {
      
        GLuint id = program_line.m_program_obj;
		
	GLint loc = glGetUniformLocation(id, "uColor");

	glProgramUniform4fv(id, loc, 1, color); //1.0f, 0.0f, 0.5f, 1.0f);
	
	program_line.use();


	
	vertex_array_v2.bind();
	
	
        glDrawArrays(GL_LINES, 0, nf);

        vertex_array_v2.unbind();        

        program_line.unuse();


    }

  
    void draw_spheres(GLsizei nv)
    {
        glEnable(GL_PROGRAM_POINT_SIZE);
        glPointParameterf(GL_POINT_SPRITE_COORD_ORIGIN, GL_LOWER_LEFT);

	program_sphere.use();

        vertex_array_v.bind();
        vertex_color_buffer.bind();

        glDrawArrays(GL_POINTS, 0, nv);

	vertex_color_buffer.unbind();
        vertex_array_v.unbind();

        program_sphere.unuse();


    }

  //    GLviz::glVertexArray  vertex_array_v;
      GLviz::glVertexArray  vertex_array_v, vertex_array_v2;
    GLviz::glArrayBuffer  vertex_array_buffer, normal_array_buffer, vertex_color_buffer, vertex_displacement_buffer;
    GLviz::glElementArrayBuffer  index_array_buffer;

    GLviz::UniformBufferCamera      uniform_camera;
    GLviz::UniformBufferMaterial    uniform_material;
    GLviz::UniformBufferWireframe   uniform_wireframe;
    GLviz::UniformBufferSphere      uniform_sphere;

  GLviz::ProgramLine   program_line;
    GLviz::ProgramSphere  program_sphere;
};

std::unique_ptr<MyViz> viz;

std::vector<Eigen::Vector3f>               m_ref_vertices;
std::vector<Eigen::Vector3f>               m_ref_normals;

std::vector<Eigen::Vector3f>               m_vertices;
std::vector<Eigen::Vector3f>               m_normals;
std::vector<std::array<unsigned int, 3> >  m_faces;

Parameters                                 m_Parameters;
ParticleManager*                           manager;
  Scenes *sceneManager;
  
TimeStepPBSPH* m_SPH; //Pressure Solver



  
void
displayFluidParticles()
{
	// display particles
	glEnable(GL_MULTISAMPLE);
	glEnable(GL_DEPTH_TEST);

	glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
	glClearDepth(1.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
   
	viz->uniform_camera.set_buffer_data(camera);

	GLviz::Frustum view_frustum = camera.get_frustum();

	float g_projection_radius = view_frustum.near_()*10 * (GLviz::screen_height()*10 / (view_frustum.top()*10 - view_frustum.bottom()*10));
	
	// boundary box
	for(int bb = 0; bb < manager->m_particleObjects.size(); bb++){


	  manager->updateFluidPositions(bb);

	  Particles* pp = manager->m_particleObjects[bb];

	  double particleSize = manager->m_fluid_position_output[bb].size();
	  
	  if(pp->render){
	  
	    viz->vertex_array_buffer.set_buffer_data(3 * sizeof(GLfloat) *  manager->m_fluid_position_output[bb].size(), manager->m_fluid_position_output[bb].front().data());

	    double radius = g_bb_radius;

	    if(pp->isFluid)
	      radius = m_Parameters.particleRadius / 8;
	    
	    viz->uniform_sphere.set_buffer_data(radius, g_projection_radius);

	    viz->vertex_color_buffer.set_buffer_data(3 * sizeof(GLfloat) * particleSize, manager->m_fluid_color_output[bb].front().data());
	    
	    viz->draw_spheres(static_cast<GLsizei>(manager->m_fluid_position_output[bb].size()));
	  }
	}

    if(m_Parameters.showDisplacement){

      for(int i = 0; i < m_SPH->m_vector_output.size(); i++){
	viz->vertex_displacement_buffer.set_buffer_data(3 * sizeof(GLfloat) * m_SPH->m_vector_output[i].size(), m_SPH->m_vector_output[i].front().data());

	viz->draw_vector(g_displacement_color[i], static_cast<GLsizei>(m_SPH->m_vector_output[i].size() * 2));
      } 
    }


}

void TW_CALL
addParticles(void*){
  std::cout << "add" << std::endl;
}

void
reshapeFunc(int width, int height)
{
    const float aspect = static_cast<float>(width) / static_cast<float>(height);

    glViewport(0, 0, width, height);
    camera.set_perspective(60.0f, aspect, 0.005f, 5.0f);
}

void
closeFunc()
{
    viz = nullptr;
}


unsigned int fileIter = 0;

void
timerFunc(int delta_t_msec)
{
    float delta_t_sec = static_cast<float>(delta_t_msec) / 1000.0f;

    if (!g_stop_simulation)
    {
        g_time += m_Parameters.stepSize * 1.0;

        // update particle positions
        m_SPH->step(m_Parameters.stepSize);

	//	manager->updateDynamicBoundary(g_time);
	
	double fps = 1/30.0;

	if(m_Parameters.exportFrames){

	  // Only write fps frames each second
	  if(fileIter * fps < g_time){

	    //Export to Video
	    std::stringstream ss;
	    ss << "frame" << std::setw(4) << std::setfill('0') << fileIter;
	    std::string s = ss.str();	
	    GLviz::export_frame("export/" + s);


        // Export fluid particle positions
        ParticleExporter::exportParticles("export/particles" + s, manager->m_fluid_position_output[0], manager->m_fluid_position_output[0].size());

	    
	    fileIter++;
	  }
	}

    }

}


void TW_CALL
export_frame(void*)
{
  GLviz::export_frame("export");
}

void
create_fluid_particles(unsigned int setup, unsigned int observeParticle)
{

  sceneManager = new Scenes(&m_Parameters);
  manager = sceneManager->createSetup(setup, m_Parameters.h);

  //  Scenes scene(&m_Parameters);

  std::cout << "open " << manager->m_particleObjects.size() << std::endl;
  m_SPH = new TimeStepPBSPH(manager, &m_Parameters);

  std::cout << "ende " << std::endl;
  m_Parameters.observeParticle = observeParticle;

}



void TW_CALL reset_simulation(void*)
{
    g_time = 0.0f;


    m_SPH->init();

}

void TW_CALL step_simulation(void*)
{
    g_time += m_Parameters.stepSize;
    m_SPH->step(m_Parameters.stepSize);
    //    manager->updateDynamicBoundary(g_time);
}
}


int
main(int argc, char* argv[])
{
  unsigned int observeParticle = -1;
  g_test_setup = 3;

  if(argc > 2){
    //  observeParticle = atoi(argv[2]);
    m_Parameters.particlesPerAxis = atoi(argv[2]);
  }
  if (argc > 1){
    g_test_setup = atoi(argv[1]);
  } 
	// Setup fluid particles
  create_fluid_particles(g_test_setup, observeParticle);

    GLviz::init(argc, argv);

    viz = std::unique_ptr<MyViz>(new MyViz());

    
    // Setup AntTweakBar.
    {
        TwBar* bar = GLviz::twbar();

        TwAddVarRO(bar, "FrameTime", TW_TYPE_DOUBLE,
            &m_Parameters.timePerFrame, " precision=0 label='ms per Frame' group='Simulation' ");
	
        TwAddVarRW(bar, "Stop", TW_TYPE_BOOLCPP,
            &g_stop_simulation,
            " key=p help='Stop simulation' group='Simulation' ");

        TwAddVarRO(bar, "time", TW_TYPE_FLOAT,
            &g_time, " precision=4 label='t in sec' group='Simulation' ");

        TwAddButton(bar, "Reset",
            reset_simulation, NULL,
            " key=r help='Reset simulation' group='Simulation' ");
	/*
	TwAddButton(bar, "Add",
            addParticles, NULL,
            " key=r help='Reset simulation' group='Simulation' ");
	*/

	TwAddVarRW(bar, "Add", TW_TYPE_BOOLCPP,
		   &m_Parameters.add,
            " key=p help='Stop simulation' group='Simulation' ");
		
        TwAddButton(bar, "Step forward",
            step_simulation, NULL,
            " key=SPACE help='Step forward' group='Simulation' ");

        TwAddVarRW(bar, "Gravity", TW_TYPE_BOOLCPP,
            &m_Parameters.useGravity,
            " key=g help='Use Gravitation' group='Simulation' ");

	TwAddVarRW(bar, "Displacement", TW_TYPE_BOOLCPP,
            &m_Parameters.showDisplacement,
            " key=d help='Show displacement' group='Simulation' ");

	TwAddVarRW(bar, "Variable Stepsize", TW_TYPE_BOOLCPP,
            &m_Parameters.enableAdaptiveTimeStep,
            " key=t help='Show displacement' group='Simulation' ");

	TwAddVarRO(bar, "Stepsize", TW_TYPE_DOUBLE,
            &m_Parameters.stepSize, " precision=6 label='Stepsize' group='Simulation' ");
		
        TwAddVarRW(bar, "Mass", TW_TYPE_DOUBLE,
            &m_Parameters.particleMassScaling,
            " min=0 max=10 step=0.01 key=4 help='B' group='Stiffness Coefficient' ");
	
        TwAddVarRW(bar, "Rest density", TW_TYPE_DOUBLE,
            &m_Parameters.restDensity,
            " min=0 max=100000 step=1 key=4 help='B' group='Stiffness Coefficient' ");
	
	TwAddVarRW(bar, "Viscosity epsilon", TW_TYPE_DOUBLE,
            &m_Parameters.viscosity,
            " min=0 max=1 step=0.01 key=3 help='B' group='Stiffness Coefficient' ");
	
	
        TwAddVarRO(bar, "Avg. Density", TW_TYPE_DOUBLE,
            &m_Parameters.avg_density, " precision=3 label='Avg. Density' group='Info' ");

	TwAddVarRO(bar, "Max. Density", TW_TYPE_DOUBLE,
            &m_Parameters.max_density, " precision=3 label='Max. Density' group='Info' ");

	TwAddVarRO(bar, "Min. Velo", TW_TYPE_DOUBLE,
            &m_Parameters.min_velo, " precision=3 label='Min. Velo' group='Info' ");
	
	TwAddVarRO(bar, "Avg. Velo", TW_TYPE_DOUBLE,
            &m_Parameters.avg_velo, " precision=3 label='Avg. Velo' group='Info' ");

	TwAddVarRO(bar, "Max. Velo", TW_TYPE_DOUBLE,
            &m_Parameters.max_velo, " precision=3 label='Max. Velo' group='Info' ");

	TwAddVarRO(bar, "Avg. Fluid Neighbors", TW_TYPE_DOUBLE,
		   &m_Parameters.avg_fluid_neighbors, " precision=3 label='Avg. Fluid Neighbors' group='Info' ");
	
	TwAddVarRO(bar, "Avg. Neighbors", TW_TYPE_DOUBLE,
		   &m_Parameters.avg_neighbors, " precision=3 label='Avg. Neighbors' group='Info' ");
	

	/*
        TwAddVarRW(bar, "Fluid", TW_TYPE_DOUBLE,
            &m_Parameters.stiffnessFluid,
            " min=0 max=100000 step=0.001 key=1 help='B' group='Stiffness Coefficient' ");


        TwAddVarRW(bar, "Bounding Box", TW_TYPE_DOUBLE,
            &m_Parameters.stiffnessBB,
            " min=0 max=10000 step=0.001 key=2 help='B' group='Stiffness Coefficient' ");
	*/



	TwAddVarRW(bar, "Export Frames", TW_TYPE_BOOLCPP,
            &m_Parameters.exportFrames,
            " key=e help='Export Frames' group='Export' ");
    }

    GLviz::display_callback(displayFluidParticles);
    GLviz::reshape_callback(reshapeFunc);
    GLviz::close_callback(closeFunc);
    //    GLviz::timer_callback(timerFunc, m_Parameters.stepSize * 1000);
    GLviz::timer_callback(timerFunc, 5);

    return GLviz::exec(camera);


}
