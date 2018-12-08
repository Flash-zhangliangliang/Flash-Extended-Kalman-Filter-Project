#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

#define EPS 0.0001 // A very small number


using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  /*
  add
  */
  H_laser_ << 1, 0, 0, 0,			  
  	          0, 1, 0, 0;
  /**/
  
  Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
        0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
        0, 0.0009, 0,
        0, 0, 0.09;

  /**
  TODO:
    * Finish initializing the FusionEKF.
    * Set the process and measurement noises
  */
  // Initial state transition matrix  
  ekf_.F_ = MatrixXd(4, 4);  
  ekf_.F_ << 1, 0, 1, 0,			
  	         0, 1, 0, 1,			
  	         0, 0, 1, 0,			
  	         0, 0, 0, 1;  
  
  // Initial covariance matrix    
  ekf_.P_ = MatrixXd(4, 4);    
  ekf_.P_ << 1, 0, 0, 0,			   
		     0, 1, 0, 0,			   
		     0, 0, 1000, 0,			   
		     0, 0, 0, 1000;

}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    /**
    TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      float rho = measurement_pack.raw_measurements_[0]; // range	  
      float phi = measurement_pack.raw_measurements_[1]; // bearing	  
      float rho_dot = measurement_pack.raw_measurements_[2]; // velocity of rho	  
      // Coordinates convertion from polar to cartesian	  
      float px = rho * cos(phi); 	  
	  float py = rho * sin(phi);	  
	  float vx = rho_dot * cos(phi);	  
	  float vy = rho_dot * sin(phi);	  
	  ekf_.x_ << px, py, vx , vy;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      // We don't know velocities from the first measurement of the LIDAR, so, we use zeros      
      ekf_.x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;
    }

    //---
	// Deal with the special case initialisation problems    
	if (fabs(ekf_.x_(0)) < EPS and fabs(ekf_.x_(1)) < EPS){		
		ekf_.x_(0) = EPS;		
		ekf_.x_(1) = EPS;	
	}	
	//------
	
	// Print the initialization results    
	cout << "EKF init: " << ekf_.x_ << endl;    
	// Save the initiall timestamp for dt calculation    
	previous_timestamp_ = measurement_pack.timestamp_;

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
   TODO:
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */
  // Calculate the timestep between measurements in seconds  
  float dt = (measurement_pack.timestamp_ - previous_timestamp_);  
  dt /= 1000000.0; // convert micros to s  
  previous_timestamp_ = measurement_pack.timestamp_;    

  // State transition matrix update  
  ekf_.F_(0, 2) = dt;
  ekf_.F_(1, 3) = dt;  

  // Noise covariance matrix computation	 
  // Noise values from the task  
  float noise_ax = 9.0;  
  float noise_ay = 9.0; 
 
  // Precompute some usefull values to speed up calculations of Q  
  float dt_2 = pow(dt, 2.0); //dt^2  
  float dt_3 = pow(dt, 3.0); //dt^3  
  float dt_4 = pow(dt, 4.0); //dt^4  
  float dt_4_4 = dt_4 / 4; //dt^4/4  
  float dt_3_2 = dt_3 / 2; //dt^3/2  
  ekf_.Q_ = MatrixXd(4, 4);  
  ekf_.Q_ << dt_4_4 * noise_ax, 0, dt_3_2 * noise_ax, 0,	         
  	         0, dt_4_4 * noise_ay, 0, dt_3_2 * noise_ay,	         
  	         dt_3_2 * noise_ax, 0, dt_2 * noise_ax, 0, 	         
  	         0, dt_3_2 * noise_ay, 0, dt_2 * noise_ay;
  //----------
  

  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    // Use Jacobian instead of H	
    ekf_.H_ = tools.CalculateJacobian(ekf_.x_);	
	ekf_.R_ = R_radar_;	
	ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  } else {
    // Laser updates
    ekf_.H_ = H_laser_;	
	ekf_.R_ = R_laser_;	
	ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
