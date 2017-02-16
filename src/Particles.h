#ifndef PARTICLES_H
#define PARTICLES_H

struct Particles {
  std::vector<Eigen::Vector3d> position;
  std::vector<Eigen::Vector3d> position0;
  std::vector<Eigen::Vector3d> velocity;
  std::vector<Eigen::Vector3d> color;
  bool isDynamic;
  bool isFluid;
  bool render;
};

struct Boundary : public Particles {

  std::vector<double> m_boundaryPsi;

};

struct Fluid : public Particles {
  std::vector<Eigen::Vector3d> velocity_viscosity;

  std::vector<double> density;
  std::vector<double> density_fac;
  std::vector<double> pressure_force;
	
  std::vector<Eigen::Vector3d> acceleration;

  std::vector<double> lambda;
  std::vector<Eigen::Vector3d> delta_x;
  std::vector<Eigen::Vector3d> old_pos;
  std::vector<Eigen::Vector3d> last_pos;
};

#endif
