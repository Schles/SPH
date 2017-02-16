#ifndef KERNEL_H
#define KERNEL_H
#include <Eigen/Core>
#include <vector>

/*
    Kernelfunction numbers:
    1: Poly6
    2: Spiky
    3: cubic spline
 */

class Kernel {
    public:
  Kernel(int kernelFunction, double h);
        double w(Eigen::Vector3d r, double h);

        Eigen::Vector3d gradient_w(Eigen::Vector3d rVec, double h);
        Eigen::Vector3d centralDifferences(Eigen::Vector3d rVec, double h);



	double w_cubicSpline(double r, double h);
	double gradient_w_cubicSpline(double r, double h);
	
        double computeDensityCoeffient(Eigen::Vector3d x, double h);

    private:
        double w_poly6(Eigen::Vector3d r, double h);
        double w_spiky(Eigen::Vector3d r, double h);
	double w_precached_cubicSpline(double r, double h);

        double gradient_w_poly6(double r, double h);
        double gradient_w_spiky(double r, double h);
	double gradient_precached_cubicSpline(double r, double h);

	double m_h;
	unsigned int m_kernelFunction = 3;

	double pi = 3.14;

	void preCacheW_CubicSpline();
	
	double m_precached_resolution;

	std::vector<double> m_precached_W_cubicSpline;
	std::vector<double> m_precached_grad_W_cubicSpline;
};


#endif
