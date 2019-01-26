#include "tools.h"

using Eigen::VectorXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

// Convert polar to cartesian coordinates

void Tools::convertPolarToCart(double &px, double &py)
{
  float rho;
  float phi;

  rho = px;
  phi = py;

  px = rho*sin(phi);
  py = rho*cos(phi);
}

// calculation of RMSE

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {

  // declare additional variables that will be used

  Eigen::VectorXd xSub(4); // subtracted vectors
  Eigen::VectorXd squared(4); // squared vectors
  VectorXd rmse(4);
  rmse << 0,0,0,0;

  ///Error Checks ///

  // vector size must not be zeros

  if (estimations.size() == 0 || ground_truth.size() == 0 )
  {
    std::cerr << "The size of either the estimation or gound truth vector is 0" << std::endl;
    return rmse;
  }

  // estimation vector size should be equal to ground truth size

  if (estimations.size() != ground_truth.size())
  {
    std::cerr << "The size of the estimation matrix does not equal the ground_truth matrix" << std::endl;
    return rmse;
  }

  // accumulate squared residuals
  for (unsigned int i=0; i < estimations.size(); ++i) {

    VectorXd residual = estimations[i] - ground_truth[i];

    // coefficient-wise multiplication
    residual = residual.array()*residual.array();
    rmse += residual;
  }

  // calculate the mean
  rmse = rmse/estimations.size();

  // calculate the squared root
  rmse = rmse.array().sqrt();

  // return the result
  return rmse;

}
