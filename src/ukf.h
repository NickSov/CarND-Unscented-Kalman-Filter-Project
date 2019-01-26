#ifndef UKF_H
#define UKF_H

#include "Eigen/Dense"
#include "measurement_package.h"

class UKF {
 public:
  /**
   * Constructor
   */
  UKF();

  /**
   * Destructor
   */
  virtual ~UKF();

  /**
   * ProcessMeasurement
   * @param meas_package The latest measurement data of either radar or laser
   */
  void ProcessMeasurement(MeasurementPackage meas_package);

  /**
   * Prediction Predicts sigma points, the state, and the state covariance
   * matrix
   * @param delta_t Time between k and k+1 in s
   */
  void Prediction(double delta_t);

  /**
   * Updates the state and the state covariance matrix using a laser measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateLidar(MeasurementPackage meas_package);

  /**
   * Updates the state and the state covariance matrix using a radar measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateRadar(MeasurementPackage meas_package);


  // initially set to false, set to true in first call of ProcessMeasurement
  bool is_initialized;

  // previous timestamp

  double previous_timestamp;

  // if this is false, laser measurements will be ignored (except for init)
  bool use_laser;

  // if this is false, radar measurements will be ignored (except for init)
  bool use_radar;

  // state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
  Eigen::VectorXd x;

  // create augmented mean vector
  Eigen::VectorXd x_aug;

  // state covariance matrix
  Eigen::MatrixXd P;

  // create augmented state covariance
  Eigen::MatrixXd P_aug;

  // augmented square root matrix
  Eigen::MatrixXd A;

  // sigma points matrix
  Eigen::MatrixXd Xsig;

  // predicted sigma points matrix
  Eigen::MatrixXd Xsig_pred;

  // covariance matrix for prediction

  Eigen::MatrixXd P_pred;

  // vector for predicted states

  Eigen::VectorXd x_pred;

  // process covariance matrix

  Eigen::MatrixXd Q;

  // time when the state is true, in us
  long long time_us;

  // Process noise standard deviation longitudinal acceleration in m/s^2
  double std_a;

  // Process noise standard deviation yaw acceleration in rad/s^2
  double std_yawdd;

  // Laser measurement noise standard deviation position1 in m
  double stdLaspx;

  // Laser measurement noise standard deviation position2 in m
  double stdLaspy;

  // Radar measurement noise standard deviation radius in m
  double stdRadr;

  // Radar measurement noise standard deviation angle in rad
  double stdRadphi;

  // Radar measurement noise standard deviation radius change in m/s
  double stdRadrd;

  // Weights of sigma points
  Eigen::VectorXd weights;

  // State dimension
  int n_x;

  // Augmented state dimension
  int n_aug;

  // Sigma point spreading parameter
  double lambda;
};

#endif  // UKF_H
