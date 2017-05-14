
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
	 P_ = F_ * P_ * F_.transpose();
}

void KalmanFilter::Update(const VectorXd &z) {
  	
		// KF Measurement update step
		VectorXd y = z - H_ * x_;
		MatrixXd s = H_ * P_ * H_.transpose() +R_;
		MatrixXd K = P_ * H_.transpose() * s.inverse();
	  
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

    float c1 = x_[0]*x_[0] + x_[1]*x_[1];
    float c2 = x_[0]*x_[3] + x_[1]*x_[4];
	  float rho = sqrt(c1);
    MatrixXd  h = MatrixXd(1,3);
    h << rho,atan2(x_[0],x_[1]),c1/c2;

		VectorXd y = z - h;
		MatrixXd s = Hj_  * P_ * Hj_.transpose() +R_;
		MatrixXd K = P_ * Hj_.transpose() * s.inverse();
	  
		// new state
		x_ = x_ + (K * y);
    long x_size = x_.size();
	  I = MatrixXd::Identity(x_size, x_size);
		P_ = (I - K * H_)*P_;
}
