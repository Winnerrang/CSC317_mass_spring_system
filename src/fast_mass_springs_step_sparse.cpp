#include "fast_mass_springs_step_sparse.h"
#include <igl/matlab_format.h>

void fast_mass_springs_step_sparse(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & E,
  const double k,
  const Eigen::VectorXi & b,
  const double delta_t,
  const Eigen::MatrixXd & fext,
  const Eigen::VectorXd & r,
  const Eigen::SparseMatrix<double>  & M,
  const Eigen::SparseMatrix<double>  & A,
  const Eigen::SparseMatrix<double>  & C,
  const Eigen::SimplicialLLT<Eigen::SparseMatrix<double> > & prefactorization,
  const Eigen::MatrixXd & Uprev,
  const Eigen::MatrixXd & Ucur,
  Eigen::MatrixXd & Unext)
{
  //////////////////////////////////////////////////////////////////////////////
  // Replace with your code
    double w = 1e10;
    Eigen::MatrixXd next_p(Ucur);
    Eigen::MatrixXd y;
    Eigen::RowVector3d spring;
    Eigen::MatrixXd d(A.rows(), Ucur.cols());
    const Eigen::MatrixXd const_y = M / (delta_t * delta_t) * (2 * Ucur - Uprev) +
        fext +
        w * C.transpose() * C * V;

    for (int iter = 0; iter < 50; iter++) {

        // local step that optimize spring direction
        for (int spring_idx = 0; spring_idx < E.rows(); spring_idx++) {
            auto spring = next_p.row(E(spring_idx, 0)) - next_p.row(E(spring_idx, 1));
            d.row(spring_idx) = spring.normalized() * r(spring_idx);
        }

        // global step that optimize future position
        y = k * A.transpose() * d + const_y;
        next_p = prefactorization.solve(y);
    }

    // can not put Unext inside the for loop
    // because it will create some weird behavior
    // reason unknown
    Unext = next_p;

  //////////////////////////////////////////////////////////////////////////////
}
