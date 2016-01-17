/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/
#include <element_promotion/QuadratureRule.h>

#include <Teuchos_RCP.hpp>
#include <Teuchos_LAPACK.hpp>
#include <Teuchos_SerialDenseMatrix.hpp>
#include <Teuchos_SerialDenseSolver.hpp>
#include <Teuchos_SerialDenseVector.hpp>

#include <cmath>
#include <vector>

namespace sierra{
namespace nalu{

std::pair<Teuchos::SerialDenseVector<int, double>, Teuchos::SerialDenseVector<int, double>>
jacobi_recursion_coefficients(
    const double alpha,
    const double beta,
    const int order)
{
  int N = order;
  Teuchos::SerialDenseVector<int, double> a(N);
  Teuchos::SerialDenseVector<int, double> b(N);

  double nu = (beta - alpha) / (alpha + beta + 2.0);
  double mu = std::pow(2.0, alpha + beta + 1.0) * std::tgamma(alpha + 1.0)
  * std::tgamma(beta + 1.0) / std::tgamma(alpha + beta + 2.0);
  double nab;
  double sqdif = beta * beta - alpha * alpha;

  a[0] = nu;
  b[0] = mu;

  if (N > 1) {
    for (int n = 1; n < N; ++n) {
      nab = 2 * n + alpha + beta;
      a[n] = sqdif / (nab * (nab + 2));
    }

    b[1] = 4.0 * (alpha + 1.0) * (beta + 1.0)
                    / (std::pow(alpha + beta + 2.0, 2) * (alpha + beta + 3.0));

    if (N > 2) {
      for (int n = 2; n < N; ++n) {
        nab = 2 * n + alpha + beta;
        b[n] = 4.0 * (n + alpha) * (n + beta) * n * (n + alpha + beta)
                        / (nab * nab * (nab + 1.0) * (nab - 1.0));
      }
    }
  }
  return std::make_pair(a, b);
}
//--------------------------------------------------------------------
std::pair<std::vector<double>, std::vector<double>>
gauss_legendre_rule(int order)
{
  int INFO;
  const int N = order;
  const char COMPZ = 'I';

  Teuchos::SerialDenseVector<int, double> D(N);
  Teuchos::SerialDenseVector<int, double> b(N);
  Teuchos::SerialDenseVector<int, double> E(N);
  Teuchos::SerialDenseVector<int, double> work(4 * N);
  Teuchos::SerialDenseMatrix<int, double> Z(N, N);

  std::tie(D, b) = jacobi_recursion_coefficients(0.0, 0.0, order);

  for (int i = 0; i < N - 1; ++i) {
    E[i] = std::sqrt(b[i + 1]);
  }

  Teuchos::LAPACK<int, double>().STEQR(
    COMPZ, N, D.values(), E.values(),
    Z.values(), N, work.values(), &INFO
  );

  std::vector<double> x(N);
  std::vector<double> w(N);
  for (int i = 0; i < N; ++i) {
    x[i] = D(i);
    w[i] = b[0] * Z(0, i) * Z(0, i);
  }
  return std::make_pair(x, w);
}
//--------------------------------------------------------------------
std::pair<Teuchos::SerialDenseVector<int, double>, Teuchos::SerialDenseVector<int, double>>
coefficients_for_lobatto(
  int order,
  double xl1,
  double xl2)
{
  const int N = order - 1;

  Teuchos::SerialDenseVector<int, double> en(N);
  Teuchos::SerialDenseVector<int, double> g(N);
  Teuchos::SerialDenseMatrix<int, double> M(N, N);
  Teuchos::SerialDenseSolver<int, double> solver;
  Teuchos::SerialDenseVector<int, double> a;
  Teuchos::SerialDenseVector<int, double> b;
  std::tie(a, b) = jacobi_recursion_coefficients(0.0, 0.0, N);

  // Nth canonical vector
  en[N - 1] = 1;

  for (int i = 0; i < N; ++i) {
    M(i, i) = a[i] - xl1;
  }

  for (int i = 0; i < N - 1; ++i) {
    double offdiag = std::sqrt(b[i + 1]);
    M(i + 1, i) = offdiag;
    M(i, i + 1) = offdiag;
  }

  solver.setMatrix(Teuchos::rcp(&M, false));
  solver.setVectors(Teuchos::rcp(&g, false), Teuchos::rcp(&en, false));
  solver.solve();
  double g1 = g[N - 1];

  M.putScalar(0.0);
  en.putScalar(0.0);
  g.putScalar(0.0);
  en[N - 1] = 1;
  for (int i = 0; i < N; ++i) {
    M(i, i) = a[i] - xl2;
  }

  for (int i = 0; i < N - 1; ++i) {
    double offdiag = std::sqrt(b[i + 1]);
    M(i + 1, i) = offdiag;
    M(i, i + 1) = offdiag;
  }

  solver.setMatrix(Teuchos::rcp(&M, false));
  solver.setVectors(Teuchos::rcp(&g, false), Teuchos::rcp(&en, false));
  solver.solve();
  double g2 = g[N - 1];

  Teuchos::SerialDenseVector<int, double> amod(N + 1);
  Teuchos::SerialDenseVector<int, double> bmod(N + 1);

  for (int i = 0; i < N; ++i) {
    amod[i] = a[i];
    bmod[i] = b[i];
  }
  amod[N] = (g1 * xl2 - g2 * xl1) / (g1 - g2);
  bmod[N] = (xl2 - xl1) / (g1 - g2);

  return std::make_pair(amod, bmod);
}
//--------------------------------------------------------------------
std::pair<std::vector<double>, std::vector<double>>
gauss_lobatto_legendre_rule(
  int order,
  double xleft,
  double xright)
{
  int INFO;
  const int N = order;
  const char COMPZ = 'I';

  Teuchos::SerialDenseVector<int, double> D(N);
  Teuchos::SerialDenseVector<int, double> b(N);
  Teuchos::SerialDenseVector<int, double> E(N);
  Teuchos::SerialDenseVector<int, double> work(4 * N);
  Teuchos::SerialDenseMatrix<int, double> Z(N, N);

  std::tie(D, b) = coefficients_for_lobatto(order, xleft, xright);

  for (int i = 0; i < N - 1; ++i) {
    E[i] = std::sqrt(b[i + 1]);
  }

  Teuchos::LAPACK<int, double>().STEQR(
    COMPZ, N, D.values(), E.values(),
    Z.values(), N, work.values(), &INFO
  );

  std::vector<double> x(N);
  std::vector<double> w(N);
  for (int i = 0; i < N; ++i) {
    x[i] = D(i);
    w[i] = b[0] * Z(0, i) * Z(0, i);
  }

  //remove machine eps from end points
  x[0] = xleft;
  x[N - 1] = xright;

  return std::make_pair(x, w);
}

}  // namespace nalu
}  // namespace sierra
