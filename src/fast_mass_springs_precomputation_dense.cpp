#include "fast_mass_springs_precomputation_dense.h"
#include "signed_incidence_matrix_dense.h"
#include <Eigen/Dense>

bool fast_mass_springs_precomputation_dense(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & E,
  const double k,
  const Eigen::VectorXd & m,
  const Eigen::VectorXi & b,
  const double delta_t,
  Eigen::VectorXd & r,
  Eigen::MatrixXd & M,
  Eigen::MatrixXd & A,
  Eigen::MatrixXd & C,
  Eigen::LLT<Eigen::MatrixXd> & prefactorization)
{
  /////////////////////////////////////////////////////////////////////////////
  double w = 1e10;
  
  signed_incidence_matrix_dense((int) V.rows(), E, A);

  r.resize(E.rows());
  for (int spring_idx = 0; spring_idx < E.rows(); spring_idx++) {
      int spring_start(E(spring_idx, 0)), spring_end(E(spring_idx, 1));
      r(spring_idx) = (V.row(spring_start) - V.row(spring_end)).norm();
  }


  M = Eigen::MatrixXd::Zero(m.size(), m.size());
  for (int vertex_idx = 0; vertex_idx < m.size(); vertex_idx++) {
      M(vertex_idx, vertex_idx) = m(vertex_idx);
  }


  C = Eigen::MatrixXd::Zero(b.size(), V.rows());
  for (int idx = 0; idx < b.size(); idx++) {
      C(idx, b(idx)) = 1;
  }

  //std::cout << C <<std::endl;
  Eigen::MatrixXd Q = k * A.transpose() * A +
          M /(delta_t * delta_t)  +
          w * C.transpose() * C;

  /////////////////////////////////////////////////////////////////////////////
  prefactorization.compute(Q);
  return prefactorization.info() != Eigen::NumericalIssue;
}
