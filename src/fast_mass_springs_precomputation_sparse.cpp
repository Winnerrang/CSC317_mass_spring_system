#include "fast_mass_springs_precomputation_sparse.h"
#include "signed_incidence_matrix_sparse.h"
#include <vector>

bool fast_mass_springs_precomputation_sparse(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & E,
  const double k,
  const Eigen::VectorXd & m,
  const Eigen::VectorXi & b,
  const double delta_t,
  Eigen::VectorXd & r,
  Eigen::SparseMatrix<double>  & M,
  Eigen::SparseMatrix<double>  & A,
  Eigen::SparseMatrix<double>  & C,
  Eigen::SimplicialLLT<Eigen::SparseMatrix<double> > & prefactorization)
{
  /////////////////////////////////////////////////////////////////////////////
  // Replace with your code
  std::vector<Eigen::Triplet<double> > ijv;
  const int n = V.rows();
  double w = 1e10;

  signed_incidence_matrix_sparse(n, E, A);

  // construct r
  r.resize(E.rows());
  for (int spring_idx = 0; spring_idx < E.rows(); spring_idx++) {
      int spring_start(E(spring_idx, 0)), spring_end(E(spring_idx, 1));
      r(spring_idx) = (V.row(spring_start) - V.row(spring_end)).norm();
  }

  // Construct M
  ijv.clear();
  M.resize(n, n);
  for (int vertex_idx = 0; vertex_idx < m.size(); vertex_idx++) {
      ijv.emplace_back(vertex_idx, vertex_idx, m(vertex_idx));
  }
  M.setFromTriplets(ijv.begin(), ijv.end());


  // Construct C
  ijv.clear();
  C.resize(b.size(), V.rows());
  for (int idx = 0; idx < b.size(); idx++) {
      ijv.emplace_back(idx, b(idx), 1);
  }
  C.setFromTriplets(ijv.begin(), ijv.end());

  Eigen::SparseMatrix<double> Q(n, n);
  Q = k * A.transpose() * A +
      M / (delta_t * delta_t) +
      w * C.transpose() * C;

  /////////////////////////////////////////////////////////////////////////////
  prefactorization.compute(Q);
  return prefactorization.info() != Eigen::NumericalIssue;
}
