/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/
#ifndef TensorProductBasis_h
#define TensorProductBasis_h

#include <vector>

namespace sierra{
namespace nalu{

class LagrangeBasis
{
public:
  LagrangeBasis(
    std::vector<std::vector<unsigned>>&  indicesMap,
    const std::vector<double>& nodeLocs
  );

  virtual ~LagrangeBasis() {};

  void set_lagrange_weights();

  std::vector<double> eval_basis_weights(
    const std::vector<double>& intgLoc) const;

  std::vector<double> eval_deriv_weights(
    const std::vector<double>& intgLoc) const;

  double tensor_lagrange_derivative(
    unsigned dimension,
    const double* x,
    const unsigned* nodes,
    unsigned derivativeDirection
  ) const;
  double tensor_lagrange_interpolant(unsigned dimension, const double* x, const unsigned* nodes) const;
  double lagrange_1D(double x, unsigned nodeNumber) const;
  double lagrange_deriv_1D(double x, unsigned nodeNumber) const;

  std::vector<std::vector<unsigned>> indicesMap_;
  unsigned numNodes1D_;
  std::vector<double> nodeLocs_;
  const unsigned dimension_;
  std::vector<double> lagrangeWeights_;
};


} // namespace nalu
} // namespace Sierra

#endif
