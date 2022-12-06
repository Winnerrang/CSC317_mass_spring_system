#include "fast_mass_springs_step_dense.h"
#include <igl/matlab_format.h>

void fast_mass_springs_step_dense(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & E,
  const double k,
  const Eigen::VectorXi & b,
  const double delta_t,
  const Eigen::MatrixXd & fext,
  const Eigen::VectorXd & r,
  const Eigen::MatrixXd & M,
  const Eigen::MatrixXd & A,
  const Eigen::MatrixXd & C,
  const Eigen::LLT<Eigen::MatrixXd> & prefactorization,
  const Eigen::MatrixXd & Uprev,
  const Eigen::MatrixXd & Ucur,
  Eigen::MatrixXd & Unext)
{
  //////////////////////////////////////////////////////////////////////////////
  // Replace with your code 

    double w = 1e10;
    Eigen::MatrixXd cur_p(Ucur), prev_p(Uprev), next_p;

    for (int iter = 0; iter < 50; iter++) {
        Eigen::MatrixXd d = Eigen::MatrixXd::Zero(A.rows(), Ucur.cols());
        for (int spring_idx = 0; spring_idx < E.rows(); spring_idx++) {
            auto spring = Ucur.row(cur_p(spring_idx, 0)) - cur_p.row(E(spring_idx, 1));
            d.row(spring_idx) = spring.normalized() * r(spring_idx);
        }

        //std::cout << d << std::endl;

        const Eigen::MatrixXd y = k * A.transpose() * d +
            1 / (delta_t * delta_t) * M * (2 * cur_p - prev_p) +
            fext +
            w * C.transpose() * C * V;
        next_p = prefactorization.solve(y);
        prev_p = cur_p;
        cur_p = next_p;
    }
    Unext = next_p;
    
  //////////////////////////////////////////////////////////////////////////////
}
