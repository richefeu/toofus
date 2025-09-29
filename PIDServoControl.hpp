// Copyright (C) <vincent.richefeu@3sr-grenoble.fr>
//
// This file is part of TOOFUS (TOols OFten USued)
//
// It can not be copied and/or distributed without the express
// permission of the authors.
// It is coded for academic purposes.
//
// Note
// Without a license, the code is copyrighted by default.
// People can read the code, but they have no legal right to use it.
// To use the code, you must contact the author directly and ask permission.

#ifndef PID_SERVO_CONTROL_HPP
#define PID_SERVO_CONTROL_HPP

class PIDServoControl {
public:
  // Kp -  proportional gain
  // Ki -  Integral gain
  // Kd -  derivative gain
  PIDServoControl(double t_Kp, double t_Ki, double t_Kd)
      : m_Kp(t_Kp), m_Ki(t_Ki), m_Kd(t_Kd), m_previous_error(0.0), m_integral(0.0) {}

  // dt -  loop interval time
  // min - minimum correction of manipulated variable
  // max - maximum correction of manipulated variable
  double correction(double target, double currentValue, double dt, double min, double max) {
    // Calculate error
    double error = target - currentValue;

    // Proportional term
    double Pcorr = m_Kp * error;

    // Integral term
    m_integral += error * dt;
    double Icorr = m_Ki * m_integral;

    // Derivative term (finite difference)
    double derivative = (error - m_previous_error) / dt;
    double Dcorr      = m_Kd * derivative;

    // Calculate total correction
    double total_corr = Pcorr + Icorr + Dcorr;

    // Restrict the correction to max/min
    if (total_corr > max) total_corr = max;
    else if (total_corr < min) total_corr = min;

    // Save error to previous error
    m_previous_error = error;

    return total_corr;
  }

private:
  double m_Kp;
  double m_Ki;
  double m_Kd;
  double m_previous_error;
  double m_integral;
};

#endif /* end of include guard: PID_SERVO_CONTROL_HPP */

#if 0
// === TEST ===
#include <cmath>
#include <fstream>
#include <iostream>

int main() {

  double dt = 0.05;
  double mini = -1.;
  double maxi = 1.;
  double Kp = 0.7;
  double Ki = 0.1;
  double Kd = 0.0;
  PIDServoControl pid = PIDServoControl(Kp, Ki, Kd);

  std::ofstream file("testPIDControl.txt");
  double val = 2.5;
  double target = 0.0;
  double t = 0.0;
  for (int i = 0; i < 200; i++) {
    t = i * dt;
    target = cos(t * t * 0.5) + t;

    val += pid.correction(target, val, dt, mini, maxi);

    file << t << ' ' << val << '\n';
  }

  return 0;
}

#endif
