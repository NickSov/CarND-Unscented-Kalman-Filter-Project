#include "ukf.h"
#include "Eigen/Dense"
#include "measurement_package.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // set start init state
  is_initialized = false;

  // previous timestamp
  previous_timestamp = 0;

  // if this is false, laser measurements will be ignored (except during init)
  use_laser = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar = true;

  // initial state vector
  x = VectorXd(5);

  // augmented state vector
  VectorXd x_aug = VectorXd(7);

  // initial covariance matrix
  P = MatrixXd(5, 5);

  // create augmented state covariance
  P_aug = MatrixXd(7, 7);

  // process covariance matrix
  Q = MatrixXd(2,2);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a = 0.2; // original was 30

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd = 0.2; // original was 30

  // State dimension
  n_x = 5;

  // Augmented state dimension
  n_aug = 7;

  // Sigma point spreading parameter
  lambda = 3 - n_x;

  // Weights of sigma points
  weights = VectorXd(2*n_aug+1);

  // augmented sigma points matrix
  Xsig = MatrixXd(n_aug, 2 * n_aug + 1);

  // vector for predicted states
  x_pred = VectorXd(n_x);

  // augmented sigma points prediction matrix
  Xsig_pred = MatrixXd(n_x, 2 * n_aug + 1);

  // covariance matrix for prediction

  P_pred = MatrixXd(n_x, n_x);

  // create augmented mean state vector
  x_aug.head(n_x) = x;

  /**
   * DO NOT MODIFY measurement noise values below.
   * These are provided by the sensor manufacturer.
   */

  // Laser measurement noise standard deviation position1 in m
  stdLaspx = 0.15;

  // Laser measurement noise standard deviation position2 in m
  stdLaspy = 0.15;

  // Radar measurement noise standard deviation radius in m
  stdRadr = 0.3;

  // Radar measurement noise standard deviation angle in rad
  stdRadphi = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  stdRadrd = 0.3;

  /**
   * End DO NOT MODIFY section for measurement noise values
   */



  /**
   * TODO: Complete the initialization. See ukf.h for other member properties.
   * Hint: one or more values initialized above might be wildly off...
   */
}

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage meas_package) {

  //// Initialization ///////////////////////////////

  if (!is_initialized) {

    // first measurement

    x << 1, 1, 1, 1, 1; //initial state matrix, filler

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      // initialize state vector via radar measurement

      double px = measurement_pack.raw_measurements_[0];
      double py = measurement_pack.raw_measurements_[1];

      tools.convertPolarToCart(px,py); // convert to cartesian coordinate

      // initialize state based on first radar measurement

      x(0) = px;
      x(1) = py;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      // initialize state vector via lidar measurements

      x(0) = measurement_pack.raw_measurements_[0];
      x(1) = measurement_pack.raw_measurements_[1];

    }

    // initial time

    previous_timestamp = measurement_pack.timestamp_;

    // done initializing, no need to predict or update
    is_initialized = true;

  /**
   * TODO: Complete this function! Make sure you switch between lidar and radar
   * measurements.
   */
 }

 // calculate the time delta

 double dt = (measurement_pack.timestamp_ - previous_timestamp) / 1000000.0;	//dt - expressed in seconds
 previous_timestamp = measurement_pack.timestamp_;

 // predict

  Prediction(dt);

 // update

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) // radar update
  {
    UpdateRadar(measurement_pack.raw_measurements_);
  } else {  // lidar update
    UpdateLidar(measurement_pack.raw_measurements_);
  }

}

void UKF::Prediction(double delta_t) {
  /**
   * TODO: Complete this function! Estimate the object's location.
   * Modify the state vector, x_. Predict sigma points, the state,
   * and the state covariance matrix.
   */

   // construct process covariance matrix

   Q << pow(std_a,2), 0,
        0, pow(std_yawdd,2);

   // augmented square root matrix

   A = P_aug.llt().matrixL();

   // create augmented covariance matrix

   P_aug.topLeftCorner(n_x,n_x) = P;
   P_aug.bottomRightCorner(q_dim,q_dim) = Q;

   //// generate augmented sigma points ///////////////////////////////

   // create temporary matrices for sigma point calculation

   MatrixXd tempX(n_aug,1);
   MatrixXd tempNeg(n_aug, n_aug);
   MatrixXd tempPos(n_aug, n_aug);

     for(int i = 0; i < n_aug; i++)
   {
     for(int j = 0; j < n_aug; j++)
     {
         tempPos(i,j) = 0;
         tempNeg(i,j) = 0;
     }
   }

   for (int i =0; i<n_aug; i++)
   {
       for(int j=0; j < n_aug; j++)
       {
         tempPos(i,j) = x_aug(i)+ (sqrt(lambda+n_x) * A(i,j));
       }
       for(int k=0; k < n_aug; k++)
       {
         tempNeg(i,k) = x_aug(i)- (sqrt(lambda+n_x) * A(i,k));
       }
   }

   // concatenate matrices to form Xsig

   Xsig << x_aug, tempPos, tempNeg;

   //// predict sigma points ///////////////////////////////

   // declare vectors to be used in prediction

   VectorXd xVector(n_x);
   VectorXd pNoiseVector(n_x);
   VectorXd integralVector(n_x);
   VectorXd xkPlusOneVector(n_x);

   // loop through the augmented sigma matrix

   for(int i = 0; i < Xsig_pred.cols(); i++)
   {
     if(Xsig(4,i)!= 0)
     {
       // declare variables to be used in prediction

       double vk = Xsig(2,i);
       double psi = Xsig(3,i);
       double psiDot = Xsig(4,i);
       double nu_a = Xsig(5,i);
       double nuDoubleDot = Xsig(6,i);

       xVector << Xsig(0,i),
                 Xsig(1,i),
                 vk,
                 psi,
                 psiDot;

       integralVector << (vk/psiDot)*(sin(psi+psiDot*dt) - sin(psi)),
                      (vk/psiDot)*(-cos(psi+psiDot*dt) + cos(psi)),
                      0,
                      psiDot*dt,
                      0;

       pNoiseVector << 0.5*pow(dt,2)*cos(psi)*nu_a,
                        0.5*pow(dt,2)*sin(psi)*nu_a,
                        dt*nu_a,
                        0.5*pow(dt,2)*nuDoubleDot,
                        dt*nuDoubleDot;

       xkPlusOneVector = xVector + pNoiseVector + integralVector;

       Xsig_pred.col(i) = xkPlusOneVector;

     } else {

       // declare variables to be used in prediction

       double vk = Xsig(2,i);
       double psi = Xsig(3,i);
       double psiDot = Xsig(4,i);
       double nu_a = Xsig(5,i);
       double nuDoubleDot = Xsig(6,i);

       xVector << Xsig(0,i),
                 Xsig(1,i),
                 vk,
                 psi,
                 psiDot;

       integralVector << vk*cos(psi)*dt,
                        vk*sin(psi)*dt,
                        0,
                        psiDot*dt,
                        0;

       pNoiseVector << 0.5*pow(dt,2)*cos(psi)*nu_a,
                      0.5*pow(dt,2)*sin(psi)*nu_a,
                      dt*nu_a,
                      0.5*pow(dt,2)*nuDoubleDot,
                      dt*nuDoubleDot;

       xkPlusOneVector = xVector + pNoiseVector + integralVector;

       Xsig_pred.col(i) = xkPlusOneVector;

     }
   }

   //// predict mean and convariance ///////////////////////////////

   double wFirst = lambda_aug/(lambda_aug+n_aug);
   double wSecond = 1/(2*(lambda_aug+n_aug));

   for(int i = 0; i < (2 * n_aug + 1); i++)
   {
     if(i == 0)
     {
       weights(i) = wFirst;
     } else {
       weights(i) = wSecond;
     }
   }

   // predict state mean

   for(int i = 0; i < (2 * n_aug + 1); i++) // iterate over columns
   {
     x_pred = x_pred + weights(i)*Xsig_pred.col(i);
   }

   // predict state convariance

   for(int i = 0; i < (2 * n_aug + 1); i++) // iterate over columns
   {
     VectorXd diff = Xsig_pred.col(i)-x_pred;
     while (diff(3)>M_PI)
     {
       diff(3)-=2.0*M_PI;
     }
     while (diff(3)<0-M_PI)
     {
       diff(3)+=2.0*M_PI;
     }
     P_pred = P_pred + weights(i)*(diff)*(diff).transpose();
   }
}

void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use lidar data to update the belief
   * about the object's position. Modify the state vector, x_, and
   * covariance, P_.
   * You can also calculate the lidar NIS, if desired.
   */

   //// Predict lidar measurement ///////////////////////////////

    // matrices and variables for lidar update

     int n_z = 2; // measurement dimension - px, py
     MatrixXd zSigLid = MatrixXd(n_z, 2*n_aug+1);   // matrix for sigma points in the lidar measurement space
     VectorXd zPredLid = VectorXd(n_z); // mean predicted measurement matrix for lidar measurements
     MatrixXd S = MatrixXd(n_z,n_z); // lidar measurement convariance matrix S
     MatrixXd crosCorMat = MatrixXd(n_x, n_z); // lidar cross correlation matrix

    // transform sigma points into the measurement space

    for(int i=0; i<(2*n_aug+1); i++)
    {
      double px = Xsig_pred(0,i);
      double py = Xsig_pred(1,i);

      zSigLid.col(i) << rhoMat,phiMat,phiDotMat;
     }

      // calculate the mean predicted measurement

     for(int i = 0; i < (2 * n_aug + 1); i++) // iterate over columns
     {
       zPredLid = zPredLid + weights(i)*zSigLid.col(i);
     }

     MatrixXd R = MatrixXd(2,2);

     R << stdLaspx*stdLaspx, 0,
          0, stdLaspy*stdLaspy;

     for(int i = 0; i < (2 * n_aug + 1); i++) // iterate over columns
     {
       VectorXd diffMeasPredict = zSigLid.col(i)-zPredLid;
       while (diffMeasPredict(1)>M_PI)
       {
         diffMeasPredict(1)-=2.0*M_PI;
       }
       while (diffMeasPredict(1)<0-M_PI)
       {
         diffMeasPredict(1)+=2.0*M_PI;
       }
       S = S + weights(i)*(diffMeasPredict)*(diffMeasPredict).transpose();
     }

      // Add measurement noise covariance

     S = S + R;

   //// Update radar state ///////////////////////////////

     /// calculate the cross correlation MatrixXd

     for(int i = 0; i < (2 * n_aug + 1); i++) // iterate over columns
     {
       VectorXd diff = Xsig_pred.col(i)-x_pred;
       VectorXd diffMeasPredict = zSigLid.col(i)-zPredLid;

       while (diff(3)>M_PI)
       {
         diff(3)-=2.0*M_PI;
       }
       while (diff(3)<0-M_PI)
       {
         diff(3)+=2.0*M_PI;
       }

       while (diffMeasPredict(1)>M_PI)
       {
         diffMeasPredict(1)-=2.0*M_PI;
       }
       while (diffMeasPredict(1)<0-M_PI)
       {
         diffMeasPredict(1)+=2.0*M_PI;
       }

       // cross correlation matrix

       crosCorMat = crosCorMat + weights(i)*diff*diffMeasPredict.transpose();

     }

     /// calculate the Kalman gain

     //MatrixXd K = MatrixXd(crosCorMat.rows(),S.cols());
     MatrixXd K = crosCorMat*S.inverse();

     // update the state mean matrix

     x = x + K*(z-zPredLid);

     // update the covariance matrix

     P = P - K*S*K.transpose();

}

void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use radar data to update the belief
   * about the object's position. Modify the state vector, x_, and
   * covariance, P_.
   * You can also calculate the radar NIS, if desired.
  */

  //// Predict radar measurement ///////////////////////////////

   // matrices and variables

   int n_z = 3; // measurement dimension - r, phi, r_dot
   MatrixXd zSigRad = MatrixXd(n_z, 2*n_aug+1);   // matrix for sigma points in the radar measurement space
   VectorXd zPredRad = VectorXd(n_z); // mean predicted measurement matrix for radar measurements
   MatrixXd S = MatrixXd(n_z,n_z); // radar measurement convariance matrix S
   MatrixXd crosCorMat = MatrixXd(n_x, n_z); // radar cross correlation matrix

   // transform sigma points into the measurement space

   for(int i=0; i<(2*n_aug+1); i++)
   {
     double px = Xsig_pred(0,i);
     double py = Xsig_pred(1,i);
     double phi = Xsig_pred(3,i);
     double nu = Xsig_pred(2,i);

     double rhoMat = sqrt(pow(px,2) + pow(py,2));
     double phiMat = atan(py/px);
     double phiDotMat = (px*cos(phi)*nu + py*sin(phi)*nu)/rhoMat;

     // division by zero? account for rho somehow...

     while (phi>M_PI)
     {
       phi-=2.0*M_PI;
     }
     while (phi<0-M_PI)
     {
       phi+=2.0*M_PI;
     }

     zSigRad.col(i) << rhoMat,phiMat,phiDotMat;

    }

     // calculate the mean predicted measurement

    for(int i = 0; i < (2 * n_aug + 1); i++) // iterate over columns
    {
      zPredRad = zPredRad + weights(i)*zSigRad.col(i);
    }

    MatrixXd R = MatrixXd(3,3);

    R << stdRadr*stdRadr, 0, 0,
         0, stdRadphi*stdRadphi, 0,
         0, 0, stdRadrd*stdRadrd;

    for(int i = 0; i < (2 * n_aug + 1); i++) // iterate over columns
    {
      VectorXd diffMeasPredict = zSigRad.col(i)-zPredRad;
      while (diffMeasPredict(1)>M_PI)
      {
        diffMeasPredict(1)-=2.0*M_PI;
      }
      while (diffMeasPredict(1)<0-M_PI)
      {
        diffMeasPredict(1)+=2.0*M_PI;
      }
      S = S + weights(i)*(diffMeasPredict)*(diffMeasPredict).transpose();
    }

     // Add measurement noise covariance

    S = S + R;

  //// Update radar state ///////////////////////////////

    /// calculate the cross correlation MatrixXd

    for(int i = 0; i < (2 * n_aug + 1); i++) // iterate over columns
    {
      VectorXd diff = Xsig_pred.col(i)-x_pred;
      VectorXd diffMeasPredict = zSigRad.col(i)-zPredRad;

      while (diff(3)>M_PI)
      {
        diff(3)-=2.0*M_PI;
      }
      while (diff(3)<0-M_PI)
      {
        diff(3)+=2.0*M_PI;
      }

      while (diffMeasPredict(1)>M_PI)
      {
        diffMeasPredict(1)-=2.0*M_PI;
      }
      while (diffMeasPredict(1)<0-M_PI)
      {
        diffMeasPredict(1)+=2.0*M_PI;
      }

      // cross correlation matrix

      crosCorMat = crosCorMat + weights(i)*diff*diffMeasPredict.transpose();

    }

    /// calculate the Kalman gain

    //MatrixXd K = MatrixXd(crosCorMat.rows(),S.cols());
    MatrixXd K = crosCorMat*S.inverse();

    // update the state mean matrix

    x = x + K*(z-zPredRad);

    // update the covariance matrix

    P = P - K*S*K.transpose();
}
