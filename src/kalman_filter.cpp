
#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
  
}

void KalmanFilter::Predict() {
   
   
   x_ = F_* x_ ;
	 P_ = F_ * P_ * F_.transpose() + Q_;
 
}

void KalmanFilter::Update(const VectorXd &z) {
  
		// KF Measurement update step
		VectorXd y = z - H_ * x_;
    MatrixXd H_T = H_.transpose();
		MatrixXd s = H_ * P_ * H_T +R_;
		MatrixXd K = P_ * H_T * s.inverse();
    
            
		// new state
		x_ = x_ + (K * y);
    long x_size = x_.size();
	  I = MatrixXd::Identity(x_size, x_size);
		P_ = (I - K * H_)*P_;

}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  	// KF Measurement update step
    Tools tl = Tools();
    MatrixXd Hj_ = tl.CalculateJacobian(x_);

    float c1 = sqrt(x_[0]*x_[0] + x_[1]*x_[1]);
    float c2 = x_[0]*x_[2] + x_[1]*x_[3];
	  float rho = sqrt(c1);

    
    VectorXd  h(3);
    
    h << rho,atan2(x_[1],x_[0]),c2/c1;


		VectorXd y = z - h;
    MatrixXd Hj_T = Hj_.transpose();
		MatrixXd s = Hj_  * P_ * Hj_T +R_;
    MatrixXd S_I = s.inverse();
		MatrixXd K = P_ * Hj_T * S_I;
		// new state
		x_ = x_ + (K * y);
    long x_size = x_.size();
	  I = MatrixXd::Identity(x_size, x_size);
		P_ = (I - K * Hj_)*P_;
}
