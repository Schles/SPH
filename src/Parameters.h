#pragma once

#include <Eigen/Core>

struct Parameters {
  const double h = 0.1;                 // Kernel smoothing size
  
  double particleRadius = 0.025; // Will be set by code
  
  double particleMass = 1.00;    // particle mass, will be set by code

  double particleMassScaling = 2.0;
  
  double restDensity = 150;       // Rest Density
  double viscosity = 0.05;        // viscosity epsilon

  unsigned int particlesPerAxis = 10;
  
  // Position based dynamics
  
  const double pb_epsilon = 1.0e-6;
  const double pb_max_iter = 1;

  // Friction

  double viscosityCoefficient = 0.1;
  
  // WCSPH
  
  double stiffnessFluid = 1.0;   // B
  double stiffnessBB = 0.7;       // B_r
  double particleMassBB = 1.0;    // boundary box particle mass  

  // Settings
  
  bool useGravity = true;         // enable gravitational force
  bool useFriction = true;
  bool showDisplacement = true;
  bool exportFrames = false;
  bool enableAdaptiveTimeStep = true;

  bool add = false;

  // Constants
  
  const double cs = 2.0;                // Speed of sound in medium
  const int veloUpdateMethod = 0;// velocity update method
  const unsigned int kernelFunctionId = 4;  // use precached cubic spline kernel fkt
  
  const Eigen::Vector3d graviationalForce = Eigen::Vector3d(0.0, -98.1, 0.0); // scale with 10 to increase sim speed

  
  // step size settings
  
  double stepSize = 0.0001;      // Iteration Step size, scale by 0.1 if more particles then 10

  const double m_defaultStepSize = 0.0001;
  const double m_cflFactor = 0.5;
  const double m_cflMaxTimeStep = 0.0005;   // minmal step size
  const double m_cflMinTimeStep = 0.0001;  // maximal step size

  // Debug

  int observeParticle = -1;
  
};
