#include "Scenes.h"

Scenes::Scenes(Parameters* params) : ParticleFactory(params) {

}


ParticleManager* Scenes::createSetup(unsigned int i, double h){

  ParticleManager* particleManager = new ParticleManager(m_Params);
  
    Eigen::Vector3d start;
    Eigen::Vector3d end;
    Eigen::Vector3d density;

    Eigen::Vector3d u;
    Eigen::Vector3d v;

    double s = 0.2;

    unsigned int size;
    
    double fluidSize = 0.3;
    double bbSize = 0.5;

    double particleRadius = 0.1;

    int fluidDensity = 10;

    
    Fluid* particles;
    Boundary* boundary;

    unsigned int particleIndex;
    
    switch(i){
        case 1: // Presentation example
	  fluidSize = 0.25;
	  fluidDensity = m_Params->particlesPerAxis;

	  particles = createFluidTranslate(Eigen::Vector3d(0.25, 0.25, 0.25), fluidSize, fluidDensity);

	  particleIndex = particleManager->addParticles(particles);

	  particleManager->fluidIndicies.push_back(particleIndex);
	  
	  particleRadius = fluidSize / fluidDensity;

	  m_Params->particleRadius = particleRadius;
	  
	  bbSize = 0.6;
	  
          start = Eigen::Vector3d(-0.05, -0.1, -0.05);;
	  
          end = Eigen::Vector3d(1.25, 0.8, 0.55);


	  createBoundingBox(particleManager, start, end, h/2);
	  
	  


	  
            break;

    case 2:
      	  fluidSize = 0.25;
	  fluidDensity = m_Params->particlesPerAxis;

	  particles = createFluidTranslate(Eigen::Vector3d(-0.5, 0.0, 0.0), fluidSize, fluidDensity);
	  particleIndex = particleManager->addParticles(particles);
	  particleManager->fluidIndicies.push_back(particleIndex);

	  particles = createFluidTranslate(Eigen::Vector3d(0.5, 0.0, 0.0), fluidSize, fluidDensity);
	  particleIndex = particleManager->addParticles(particles);
	  particleManager->fluidIndicies.push_back(particleIndex);
	  

	  
	  /*
	  particles = createFluidTranslate(Eigen::Vector3d(0.0, -0.1, 0.0), fluidSize, fluidDensity);

	  particleIndex = particleManager->addParticles(particles);

	  particleManager->fluidIndicies.push_back(particleIndex);

*/	  

	  start = Eigen::Vector3d(-0.8, -0.6, -0.3);;
	  
          end = Eigen::Vector3d(0.85, 1.0, 0.3);


	  createBoundingBox(particleManager, start, end, h/2);
	  
      break;
	    /*
        case 2: // Breaking Dam
	  fluidSize = 0.25;
	  fluidDensity = m_Params->particlesPerAxis;	  

	  
	  
	  initParticleBuffer(pow(fluidDensity, 3));
	  
	  Particles* particles = createFluidTranslate(Eigen::Vector3d(0.25, 0.25, 0.25), fluidSize, fluidDensity);


	  
	  particleRadius = fluidSize / fluidDensity;

	  m_Params->particleRadius = particleRadius;
	  
	  bbSize = 0.6;
	  
          start = Eigen::Vector3d(-0.05, -0.1, -0.05);;
          end = Eigen::Vector3d(1.25, 0.8, 0.55);


	  createBoundingBox(start, end, h/2, false);

	  break;


        case 3: // Double Breaking Dam

	  fluidSize = 0.25;
	  fluidDensity = m_Params->particlesPerAxis;	  

	  initParticleBuffer(2 * pow(fluidDensity, 3));
	  
	  createFluidTranslate(Eigen::Vector3d(-0.5, 0.5, 0.0), fluidSize, fluidDensity);
	  createFluidTranslate(Eigen::Vector3d(0.5, 0.1, 0.0), fluidSize, fluidDensity);

	  
	  particleRadius = fluidSize / fluidDensity;


	  m_Params->particleRadius = particleRadius;
	  
	  bbSize = 0.6;

          start = Eigen::Vector3d(-0.8, -0.3, -0.3);
          end = Eigen::Vector3d(0.8, 1.6, 0.3);

	  createBoundingBox(start, end, h/2, false);

          break;    

        case 4: // Double Breaking Dam

	  fluidSize = 0.25;
	  fluidDensity = m_Params->particlesPerAxis;	  

	  
	  
	  initParticleBuffer(pow(fluidDensity, 3));
	  
	  createFluidTranslate(Eigen::Vector3d(0.25, 0.25, 0.25), fluidSize, fluidDensity);


	  
	  particleRadius = fluidSize / fluidDensity;

	  m_Params->particleRadius = particleRadius;
	  
	  bbSize = 0.6;
	  
          start = Eigen::Vector3d(-0.05, -0.1, -0.05);;
          end = Eigen::Vector3d(1.25, 0.8, 0.55);



	  createBoundingBox(start, end, h/2, true);
	  //	  createMovingBoundingBox(start, end, h/2);
	  
          break;    


    case 5: // Double Breaking Dam

	  fluidSize = 0.25;
	  fluidDensity = m_Params->particlesPerAxis;	  

	  initParticleBuffer(2 * pow(fluidDensity, 3));
	  
	  createFluidTranslate(Eigen::Vector3d(-0.5, 0.5, 0.0), fluidSize, fluidDensity);
	  createFluidTranslate(Eigen::Vector3d(0.5, 0.1, 0.0), fluidSize, fluidDensity);

	  
	  particleRadius = fluidSize / fluidDensity;


	  m_Params->particleRadius = particleRadius;
	  
	  bbSize = 0.6;

          start = Eigen::Vector3d(-0.8, -0.3, -0.3);
          end = Eigen::Vector3d(0.8, 1.6, 0.3);

	  createBoundingBox(start, end, h/2, false);

	 
	  
	  createBoundingBall(Eigen::Vector3d(0.0, 0.0, 0.0), h/2);

          break;    
	    */
        default:
            break;
    }

    return particleManager;
}

void Scenes::createBoundingBox(ParticleManager* manager, Eigen::Vector3d start, Eigen::Vector3d end, double h){

  Boundary* boundary;

  
  // Y Achse unten
  boundary= createBoundingWall(start, Eigen::Vector3d(end[0] - start[0],0,0), Eigen::Vector3d(0,0,end[2] - start[2]), h);
  manager->addParticles(boundary);

  // Y Achse oben
  boundary= createBoundingWall(Eigen::Vector3d(start[0],end[1],start[2]), Eigen::Vector3d(end[0] - start[0],0,0), Eigen::Vector3d(0,0,end[2] - start[2]), h);
  manager->addParticles(boundary);

  // X Achse links

  boundary= createBoundingWall(start, Eigen::Vector3d(0,end[1] - start[1],0), Eigen::Vector3d(0,0,end[2] - start[2]), h);
  manager->addParticles(boundary);
  
  // X Achse rechts

  boundary= createBoundingWall(Eigen::Vector3d(end[0], start[1],start[2]), Eigen::Vector3d(0,end[1] - start[1],0), Eigen::Vector3d(0,0,end[2] - start[2]), h);
  manager->addParticles(boundary);

  // Z Achse hinten

  boundary= createBoundingWall(start, Eigen::Vector3d(0,end[1] - start[1],0), Eigen::Vector3d(end[0] - start[0],0,0), h);
  manager->addParticles(boundary);
  
  // Z Achse vorne

  boundary= createBoundingWall(Eigen::Vector3d(start[0], start[1], end[2]), Eigen::Vector3d(0,end[1] - start[1],0), Eigen::Vector3d(end[0] - start[0],0,0), h);
  boundary->render = false;
  manager->addParticles(boundary);
}
