#include "ParticleFactory.h"

#include <math.h>
#define PI 3.14159265

ParticleFactory::ParticleFactory(Parameters* params) {
  m_Params = params;
}


/*
 *  start: Start Vector of the Particle Field
 *  end: End Vector
 *  density: Vector with amount of particle in the axis. density[1] = 10, 10 particles on the y axis
 */

Fluid* ParticleFactory::createFluidParticles(Eigen::Vector3d start, Eigen::Vector3d end, Eigen::Vector3d density){



  unsigned int memAllocs = density[0] * density[1] * density[2];


  Fluid* particles = createFluidObject(memAllocs);
  particles->isFluid = true;

//  unsigned int offset = particlesFluid;
  unsigned int offset = 0;

  
//  particlesFluid += memAllocs;

  int memAddr = 0;

  // create water particles
  for(int i = 0; i < density[0]; ++i)
    for (int j = 0; j < density[1]; ++j)
      for(int k = 0; k < density[2]; ++k)
	{

	  double x;
	  if(density[0] == 1)
	    x = start[0] + (end[0] - start[0]) / 2.0;
	  else
	    x = start[0] + i * (end[0] - start[0]) / (density[0] - 1);


	  double y;
	  if(density[1] == 1)
	    y = start[1] + (end[1] - start[1]) / 2.0;
	  else
	    y = start[1] + j * (end[1] - start[1]) / (density[1] - 1);

	  double z;
	  if(density[2] <= 1)
	    z = start[2] + (end[2] - start[2]) / 2.0;
	  else
	    z = start[2] + k * (end[2] - start[2]) / (density[2] - 1);


	  memAddr = offset + i * density[1] * density[2] + j * density[2] + k;


	  particles->position[memAddr] = Eigen::Vector3d(x, y, z);
	  particles->position0[memAddr] = Eigen::Vector3d(x, y, z);

	  // particles->velocity[memAddr] = Eigen::Vector3d(0.0, 0.0, 0.0);
	  particles->color[memAddr] = Eigen::Vector3d(0.0, 0.0, 0.0);;

	}

  //  m_particleObjects.push_back(particles);
    
std::cout << "Created " << memAllocs << " particles." << std::endl;
return particles;
}

Boundary* ParticleFactory::createBoundingBox(Eigen::Vector3d start, Eigen::Vector3d end, double h, bool isDynamic){


  int density[3] = {
    (end[0] - start[0]) / h,
    (end[1] - start[1]) / h,
    (end[2] - start[2]) / h 
  };
  unsigned int memAllocs = 2 * density[0] * (density[1] + 1) + 2 * (density[1] + 1) * density[2] + 2 * (density[0] - 1) * (density[2] - 1);

 
  if(isDynamic)
    memAllocs -= (density[1] + 1) * density[2] - 2 * density[2] - density[1] + 1;
  
  //  particlesBB += memAllocs;

  Boundary* particles = createBoundingObject(memAllocs);
  particles->isDynamic = false;

//  resizeOutputBuffer(m_particleObjects.size(), memAllocs);
    
  int memOffset = 0;

  Eigen::Vector3d boundryColor(0.0,0.0,0.0);

  // Z Planes
  for(int i = 0; i < density[1] + 1; ++i)
    for (int j = 0; j < density[0]; ++j)
      {
	double x = start[0] + j * (end[0] - start[0]) / density[0];
	double y = start[1] + i * (end[1] - start[1]) / density[1];
	double z = start[2];
	    
	int memAddr = memOffset + i * density[0] + j;
	//std::cout << memAddr << ": " << i << ", " << j << std::endl;
	    
	particles->position[memAddr] = { x, y, z };
	particles->color[memAddr] = boundryColor;
	    

	int a = density[0]*(1 + density[1]);
	int b = i * density[0];


	memAddr = memOffset + ((density[0])*(density[1] + 1)) + (i * density[0]) + j;

	x += (end[0] - start[0]) / density[0];
	z = end[2];
	    
	particles->position[memAddr] = { x, y, z };
	particles->color[memAddr] = boundryColor;

      }

  unsigned int deltaOffset = 0;

  unsigned int debug = 0;
  
  // X Planes
  memOffset += 2 * (density[0])*(density[1]+1);

  unsigned int x_memAddr = memOffset;
  
  for(int i = 0; i < density[1] + 1; ++i)
    for (int j = 0; j < density[2]; ++j)
      {
	double x = end[0];
	double y = start[1] + i * (end[1] - start[1]) / density[1];
	double z = start[2] + j * (end[2] - start[2]) / density[2];

	bool skip = isDynamic && !(i == 0 || i == density[1] || j == 0);

	if(skip){
	  deltaOffset++;
	} else {	    
	  particles->position[x_memAddr] = { x, y, z };
	  particles->color[x_memAddr] = boundryColor;
	  x_memAddr++;
	}



	//memAddr = memOffset + (density[1] + 1)*(density[2]) + i * density[2] + j - deltaOffset;
	
	x = start[0];
	z += (end[2] - start[2]) / density[2];

	particles->position[x_memAddr] = { x, y, z };
	particles->color[x_memAddr] = boundryColor;
	x_memAddr++;

      }

  // Y Planes
  memOffset += 2 * (density[1] + 1) * density[2] - deltaOffset;

  for (int i = 0; i < density[0] - 1; ++i)
    for (int j = 0; j < density[2] - 1; ++j)
      {
	double x = start[0] + ( i + 1 ) * (end[0] - start[0]) / density[0];;
	double y = start[1];
	double z = start[2] + ( j + 1 ) * (end[2] - start[2]) / density[2];

	int memAddr = memOffset + i * (density[2] - 1) + j;

	particles->position[memAddr] = { x, y, z };
	particles->color[memAddr] = boundryColor;

	memAddr = memOffset + (density[0] - 1 )*(density[2] - 1) + i * (density[2] - 1) + j;
	y = end[1];

	particles->position[memAddr] = { x, y, z };
	particles->color[memAddr] = boundryColor;

      }

    
  //  m_particleObjects.push_back(particles);
  std::cout << "Created " << memAllocs << " bounding particles." << std::endl;
  return particles;
}

Boundary* ParticleFactory::createBoundingWall(Eigen::Vector3d p, Eigen::Vector3d u, Eigen::Vector3d v, double h){
  unsigned int xSize = u.norm() / h;
  unsigned int ySize = v.norm() / h;
  
  unsigned int memAllocs = xSize * ySize;

  std::cout << "u " << u.norm() << std::endl;
  std::cout << " h" << h << std::endl;
  
  Boundary* particles = createBoundingObject(memAllocs);

  int memAddr = 0;

  Eigen::Vector3d res;
  
  Eigen::Vector3d boundryColor(0.0,0.0,0.0);
  
  for(int i = 0; i < xSize; i ++)
    for (int j = 0; j < ySize; ++j)
      {
       	res = p + (1.0 * i / xSize) * u + (1.0 * j/ ySize) * v;
	    
	particles->position[memAddr] = res;
	particles->position0[memAddr] = res;
	particles->color[memAddr] = boundryColor;
	memAddr++;
      }


  return particles;
}

Fluid* ParticleFactory::createFluidObject(unsigned int memAllocs){
  Fluid* particles = new Fluid();
  particles->isFluid = true;
  particles->render = true;
  
  particles->position.resize(memAllocs);
  particles->position0.resize(memAllocs);
  particles->velocity.resize(memAllocs);
  particles->color.resize(memAllocs);


  particles->density.resize(memAllocs);


  particles->acceleration.resize(memAllocs);

  particles->lambda.resize(memAllocs);

  particles->delta_x.resize(memAllocs);
  particles->old_pos.resize(memAllocs);
  particles->last_pos.resize(memAllocs);

  return particles;
}

Boundary* ParticleFactory::createBoundingObject(unsigned int memAllocs){
  Boundary* particles = new Boundary();
  particles->isFluid = false;
  particles->render = true;
  
  particles->position.resize(memAllocs);
  particles->position0.resize(memAllocs);
  particles->velocity.resize(memAllocs);
  particles->color.resize(memAllocs);

  particles->m_boundaryPsi.resize(memAllocs);

  return particles;
}


Boundary* ParticleFactory::createBoundingBall(Eigen::Vector3d start, double h){

  Boundary* particles = new Boundary();
  particles->isDynamic = true;
 
  int memAddr = 0;
  
  std::cout << "create " << std::endl;

   double param, result;


   unsigned int density = 360 / 16;

  unsigned int memAllocs = density * density;

   
  //  particlesBB += memAllocs;

  particles->position.resize(memAllocs);
  particles->position0.resize(memAllocs);
  particles->color.resize(memAllocs);
  particles->m_boundaryPsi.resize(memAllocs);

  //  resizeOutputBuffer(m_particleObjects.size(), memAllocs);
    
  Eigen::Vector3d boundryColor(1.0,0.5,0.5);

   double r = 0.1;


   
  for(double alpha = 0; alpha <= 360; alpha += density){
    for( double beta = 0; beta <= 360; beta += density){

      double x = cos(alpha * PI / 180.0) * sin(beta * PI / 180.0) * r;
      double y = sin(alpha * PI / 180.0) * sin(beta * PI / 180.0) * r;

      double z = cos(beta * PI / 180.0) * r;
      
      std::cout << x << " , " << y << "," << z << std::endl;
      particles->position0[memAddr] = { x, y, z };
      particles->color[memAddr] = boundryColor;
      memAddr++;
      
    }
  }
  

  //m_particleObjects.push_back(particles);
  return particles;
  
  //  std::cout << "Created " << memAllocs << " moving bounding particles." << std::endl; 
  
}

Fluid* ParticleFactory::createFluidTranslate(Eigen::Vector3d vec, double fluidSize, double d){
  Eigen::Vector3d start = Eigen::Vector3d(-fluidSize, -fluidSize, -fluidSize);
  Eigen::Vector3d end = Eigen::Vector3d(fluidSize, fluidSize, fluidSize);

  Eigen::Vector3d density = Eigen::Vector3d(d, d, d);

  return createFluidParticles(vec + start, vec + end, density);
}
