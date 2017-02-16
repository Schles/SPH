#include "Kernel.h"
#include <iostream>

#define _USE_MATH_DEFINES

#include <math.h>
#include <cmath>


Kernel::Kernel(int kernelFunction, double h){
  m_h = h;
  m_kernelFunction = kernelFunction;
  pi = static_cast<double>(M_PI);

  m_precached_resolution = 0.0001;
  
  preCacheW_CubicSpline();
 
}

double Kernel::w_cubicSpline(double r, double h){
    double q = r / h;

    double alpha = 8.0 / ( pi * pow(h,3));
    
    if(q <= 1){
      if(q <= 0.5)
        return alpha * ( 6.0 * pow(q, 3) - 6.0 * pow(q, 2) + 1.0 );
      else
        return alpha * ( 2.0 * pow( 1.0 - q , 3));
    }

    return 0;
}

double Kernel::gradient_w_cubicSpline(double r, double h){
 
  //  double r = rVec.norm();
  
  double q = r / h;
  double alpha = ( 1.0 / (r * h)) * (48.0 / ( pi * pow(h, 3) ) ) ;


  
  if(q <= 1.0){
    if(r > 1.0e-6){
      if(q <= 0.5){
	return alpha * q * (3.0 * q - 2.0);
      } else {
	return alpha * ( -1.0 * pow(1.0 - q, 2) );
      }
    }
  }


  
  return 0.0;

}


double Kernel::w_poly6(Eigen::Vector3d rVec, double h){

    double r = rVec.squaredNorm();
    if( 0 <= r && r <= h * h){
        double factor = 315.0 / ( 64 * M_PI * pow(h, 9));
        double delta = pow(h, 2) - r ;
        double term = pow(delta, 3);
		double result = factor * term;
		return result;
    }

    return 0;
}


double Kernel::gradient_w_poly6(double r, double h){

    if( 0 <= r && r <= h){
        double factor = 315 / ( 64 * 3.14 * pow(h, 9));
        double delta = pow(h, 2) - pow(r, 2);
        double term = -6.0 * pow(delta, 2);

        double result = factor * term;

        return result;
    }

    return 0;

}


double Kernel::w_spiky(Eigen::Vector3d rVec, double h){

    double r = rVec.norm();
    if( 0 <= r && r <= h){
        double factor = 15 / ( M_PI * pow(h, 6));
        double term = pow(h - r, 3);

        double result = factor * term;

        return result;
    }

    return 0;

}

double Kernel::gradient_w_spiky(double r, double h){

    if( 0 <= r && r <= h){
        double factor = 15 / ( M_PI * pow(h, 6));
        double term = -1.0 * pow(h - r, 2);

        double result = factor * term;

        return result;
    }

    return 0;

}



double Kernel::w(Eigen::Vector3d rVec, double h){
    double res;

    if(m_kernelFunction == 1) { // poly 6
        res = w_poly6(rVec, h);
    } else if(m_kernelFunction == 2) { // Spiky
        res = w_spiky(rVec, h);
    } else if(m_kernelFunction == 3) { // Cubic Spline
      res = w_cubicSpline(rVec.norm(), h);
    } else if(m_kernelFunction == 4) {
      res = w_precached_cubicSpline(rVec.norm(), h);
    }

    return res;
}


Eigen::Vector3d Kernel::gradient_w(Eigen::Vector3d rVec, double h){
    Eigen::Vector3d resultVec;


    if(m_kernelFunction == 1) { // poly 6
        resultVec[0] = gradient_w_poly6(fabs(rVec[0]), h);
        resultVec[1] = gradient_w_poly6(fabs(rVec[1]), h);
        resultVec[2] = gradient_w_poly6(fabs(rVec[2]), h);
    } else if(m_kernelFunction == 2) { // Spiky
        resultVec[0] = gradient_w_spiky(fabs(rVec[0]), h);
        resultVec[1] = gradient_w_spiky(fabs(rVec[1]), h);
        resultVec[2] = gradient_w_spiky(fabs(rVec[2]), h);
    } else if(m_kernelFunction == 3) { // Cubic Spline
      return gradient_w_cubicSpline(rVec.norm(), h) * rVec;
    } else if(m_kernelFunction == 4) { // Precached Cubic Spline
      return gradient_precached_cubicSpline(rVec.norm(), h) * rVec;
    }

    return resultVec;
}

Eigen::Vector3d Kernel::centralDifferences(Eigen::Vector3d rVec, double h){

    double epsilon = pow(10, -6);
    Eigen::Vector3d e_x = Eigen::Vector3d(1.0, 0.0, 0.0);
    Eigen::Vector3d e_y = Eigen::Vector3d(0.0, 1.0, 0.0);
    Eigen::Vector3d e_z = Eigen::Vector3d(0.0, 0.0, 1.0);


    Eigen::Vector3d resultVec;
    /*
    if(m_kernelFunction == 1) { // poly 6
        resultVec[0] = w_poly6( rVec + epsilon * e_x , h) - w_poly6( rVec - epsilon * e_x , h);
        resultVec[1] = w_poly6( rVec + epsilon * e_y , h) - w_poly6( rVec - epsilon * e_y , h);
        resultVec[2] = w_poly6( rVec + epsilon * e_z , h) - w_poly6( rVec - epsilon * e_z , h);
    } else if(m_kernelFunction == 2) { // Spiky
        resultVec[0] = w_spiky( rVec + epsilon * e_x , h) - w_spiky( rVec - epsilon * e_x , h);
        resultVec[1] = w_spiky( rVec + epsilon * e_y , h) - w_spiky( rVec - epsilon * e_y , h);
        resultVec[2] = w_spiky( rVec + epsilon * e_z , h) - w_spiky( rVec - epsilon * e_z , h);
    } else if(m_kernelFunction == 3) { // Cubic Spline
        resultVec[0] = w( rVec + epsilon * e_x , h) - w_cubicSpline( rVec - epsilon * e_x , h);
        resultVec[1] = w( rVec + epsilon * e_y , h) - w_cubicSpline( rVec - epsilon * e_y , h);
        resultVec[2] = w( rVec + epsilon * e_z , h) - w_cubicSpline( rVec - epsilon * e_z , h);
    }
*/
    resultVec *= 1.0 / (2.0 * epsilon);

    return resultVec;
}


double Kernel::computeDensityCoeffient(Eigen::Vector3d x, double h){

    double norm = x.norm();

	double q = norm / h;

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

    return q / norm;

}

double Kernel::gradient_precached_cubicSpline(double r, double h){
  unsigned int i = floor(r / m_precached_resolution);
  
  unsigned int size = ceil(m_h / m_precached_resolution);
  
  if( i < size - 1){   
    return m_precached_grad_W_cubicSpline[i] - (m_precached_grad_W_cubicSpline[i] - m_precached_grad_W_cubicSpline[i + 1]) / 2.0;
  }

  return 0.0;
}


double Kernel::w_precached_cubicSpline(double r, double h){
   
  unsigned int i = floor(r / m_precached_resolution);
  
  unsigned int size = ceil(m_h / m_precached_resolution);
  
  if( i < size - 1){
    return m_precached_W_cubicSpline[i] - (m_precached_W_cubicSpline[i] - m_precached_W_cubicSpline[i + 1]) / 2.0;
  }

  return 0.0;
}

void Kernel::preCacheW_CubicSpline(){
  unsigned int size = ceil(m_h / m_precached_resolution);
  
  m_precached_W_cubicSpline.resize(size);
  m_precached_grad_W_cubicSpline.resize(size);
  
  for(int i = 0; i < size; i++){
    m_precached_W_cubicSpline[i] = w_cubicSpline(i * m_precached_resolution, m_h);
    m_precached_grad_W_cubicSpline[i] = gradient_w_cubicSpline(i * m_precached_resolution, m_h);
  }
}
