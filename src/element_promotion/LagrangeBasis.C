/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/
#include <element_promotion/LagrangeBasis.h>

#include <stk_util/environment/ReportHandler.hpp>

#include <algorithm>
#include <cmath>
#include <memory>
#include <string>

namespace sierra{
namespace nalu{

LagrangeBasis::LagrangeBasis(
  std::vector<std::vector<unsigned>>&  indicesMap,
  const std::vector<double>& nodeLocs)
  :  indicesMap_(indicesMap),
     numNodes1D_(nodeLocs.size()),
     nodeLocs_(nodeLocs),
     dimension_(indicesMap[0].size())
{
  for (auto& indices : indicesMap) {
    ThrowRequire(indices.size() == dimension_);
  }
  set_lagrange_weights();
}
//--------------------------------------------------------------------------
void
LagrangeBasis::set_lagrange_weights()
{
  lagrangeWeights_.assign(numNodes1D_,1.0);
  for (unsigned i = 0; i < numNodes1D_; ++i) {
    for (unsigned j = 0; j < numNodes1D_; ++j) {
      if ( i != j ) {
        lagrangeWeights_[i] *= (nodeLocs_[i]-nodeLocs_[j]);
      }
    }
    lagrangeWeights_[i] = 1.0 / lagrangeWeights_[i];
  }
}
//--------------------------------------------------------------------------
std::vector<double>
LagrangeBasis::eval_basis_weights(
  const std::vector<double>& intgLoc) const
{
  auto numIps = intgLoc.size() / dimension_;
  ThrowAssert(numIps * dimension_ == intgLoc.size());

  auto numNodes = std::pow(numNodes1D_, dimension_);
  std::vector<double> interpolationWeights(numIps*numNodes);

  for (unsigned ip = 0; ip < numIps; ++ip) {
    unsigned scalar_ip_offset = ip*numNodes;
    for (unsigned nodeNumber = 0; nodeNumber < numNodes; ++nodeNumber) {
      unsigned scalar_offset = scalar_ip_offset+nodeNumber;
      unsigned vector_offset = ip * dimension_;
      interpolationWeights[scalar_offset]
          = tensor_lagrange_interpolant( dimension_,
                                        &intgLoc[vector_offset],
                                         indicesMap_[nodeNumber].data() );
    }
  }
  return interpolationWeights;
}
//--------------------------------------------------------------------------
std::vector<double>
LagrangeBasis::eval_deriv_weights(const std::vector<double>& intgLoc) const
{
  auto numIps = intgLoc.size()/dimension_;
  auto numNodes = std::pow(numNodes1D_,dimension_);
  std::vector<double> derivWeights(numIps * numNodes * dimension_);

  unsigned derivIndex = 0;
  for (unsigned ip = 0; ip < numIps; ++ip) {
    for (unsigned nodeNumber = 0; nodeNumber < numNodes; ++nodeNumber) {
      unsigned vector_offset = ip*dimension_;
      for (unsigned derivDirection = 0; derivDirection < dimension_; ++derivDirection) {
        derivWeights[derivIndex]
            = tensor_lagrange_derivative( dimension_,
                                         &intgLoc[vector_offset],
                                          indicesMap_[nodeNumber].data(),
                                          derivDirection );
        ++derivIndex;
      }
    }
  }
  return derivWeights;
}
//--------------------------------------------------------------------------
double
LagrangeBasis::tensor_lagrange_interpolant(
  unsigned dimension,
  const double* x,
  const unsigned* nodes) const
{
  double interpolant_weight = 1.0;
  for (unsigned j = 0; j < dimension; ++j) {
    interpolant_weight *= lagrange_1D(x[j], nodes[j]);
  }
  return interpolant_weight;
}
//--------------------------------------------------------------------------
double
LagrangeBasis::tensor_lagrange_derivative(
  unsigned dimension,
  const double* x,
  const unsigned* nodes,
  unsigned derivativeDirection) const
{
  double derivativeWeight = 1.0;
  for (unsigned j = 0; j < dimension; ++j) {
    if (j == derivativeDirection) {
      derivativeWeight *= lagrange_deriv_1D(x[j], nodes[j]);
    }
    else {
      derivativeWeight *= lagrange_1D(x[j], nodes[j]);
    }
  }
  return derivativeWeight;
}
//--------------------------------------------------------------------------
double
LagrangeBasis::lagrange_1D(double x, unsigned nodeNumber) const
{
  double numerator = 1.0;
  for (unsigned j = 0; j < numNodes1D_; ++j) {
    if (j != nodeNumber) {
      numerator *= (x - nodeLocs_[j]);
    }
  }
  return (numerator * lagrangeWeights_[nodeNumber]);
}
//--------------------------------------------------------------------------
double
LagrangeBasis::lagrange_deriv_1D(double x, unsigned nodeNumber) const
{
  double outer = 0.0;
  for (unsigned j = 0; j < numNodes1D_; ++j) {
    if (j != nodeNumber) {
      double inner = 1.0;
      for (unsigned i = 0; i < numNodes1D_; ++i) {
        if (i != j && i != nodeNumber) {
          inner *= (x - nodeLocs_[i]);
        }
      }
      outer += inner;
    }
  }
  return (outer * lagrangeWeights_[nodeNumber]);
}

}  // namespace nalu
} // namespace sierra
