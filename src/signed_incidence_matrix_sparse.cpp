#include "signed_incidence_matrix_sparse.h"
#include <vector>

void signed_incidence_matrix_sparse(
  const int n,
  const Eigen::MatrixXi & E,
  Eigen::SparseMatrix<double>  & A)
{
  //////////////////////////////////////////////////////////////////////////////
  // Replace with your code
      std::vector<Eigen::Triplet<double> > ijv;
      for (int spring_idx = 0; spring_idx < E.rows(); spring_idx++) {
          ijv.emplace_back(spring_idx, E(spring_idx, 0), 1);
          ijv.emplace_back(spring_idx, E(spring_idx, 1), -1);
      }

      A.resize(E.rows(),n);
      A.setFromTriplets(ijv.begin(),ijv.end());
  //////////////////////////////////////////////////////////////////////////////
}
