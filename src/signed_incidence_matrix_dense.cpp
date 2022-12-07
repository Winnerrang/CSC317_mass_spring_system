#include "signed_incidence_matrix_dense.h"

void signed_incidence_matrix_dense(
  const int n,
  const Eigen::MatrixXi & E,
  Eigen::MatrixXd & A)
{
  //////////////////////////////////////////////////////////////////////////////

  A = Eigen::MatrixXd::Zero(E.rows(),n);
  for (int spring_idx = 0; spring_idx < E.rows(); spring_idx++) {
      A(spring_idx, E(spring_idx, 0)) = 1;
      A(spring_idx, E(spring_idx, 1)) = -1;
  }

  //////////////////////////////////////////////////////////////////////////////
}
