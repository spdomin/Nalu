/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/
#include <element_promotion/TensorProductQuadratureRule.h>
#include <element_promotion/QuadratureRule.h>

#include <cmath>
#include <stdexcept>

namespace sierra{
namespace nalu{

TensorProductQuadratureRule::TensorProductQuadratureRule(
  std::string  /*type*/,
  int numQuad,
  std::vector<double>& scsLocs)
{
  scsEndLoc_.resize(scsLocs.size()+2);
  scsEndLoc_[0] = -1.0;

  for (unsigned j = 0; j < scsLocs.size();++j) {
    scsEndLoc_[j+1] = scsLocs[j];
  }
  scsEndLoc_[scsLocs.size()+1] = +1.0;

  std::tie(abscissae_, weights_) = gauss_legendre_rule(numQuad);
  double isoparametricFactor = 0.5;
  for (auto& weight : weights_) {
    weight *= isoparametricFactor;
  }
}
//--------------------------------------------------------------------------
double
TensorProductQuadratureRule::isoparametric_mapping(
  const double b,
  const double a,
  const double xi) const
{
  return (0.5*(xi*(b-a) + (a+b)));
}
//--------------------------------------------------------------------------
double
TensorProductQuadratureRule::gauss_point_location(
  int nodeOrdinal,
  int gaussPointOrdinal) const
{

  double location1D =
      isoparametric_mapping( scsEndLoc_[nodeOrdinal+1],
                             scsEndLoc_[nodeOrdinal],
                             abscissae_[gaussPointOrdinal] );
   return location1D;
}
//--------------------------------------------------------------------------
double
TensorProductQuadratureRule::tensor_product_weight(
  int s1Node, int s2Node, int s3Node,
  int s1Ip, int s2Ip, int s3Ip) const
{
  const double Ls1 = scsEndLoc_[s1Node+1]-scsEndLoc_[s1Node];
  const double Ls2 = scsEndLoc_[s2Node+1]-scsEndLoc_[s2Node];
  const double Ls3 = scsEndLoc_[s3Node+1]-scsEndLoc_[s3Node];
  const double isoparametricArea = Ls1 * Ls2 * Ls3;

  const double weight =
      isoparametricArea * weights_[s1Ip] * weights_[s2Ip] * weights_[s3Ip];

  return weight;

}
//--------------------------------------------------------------------------
double
TensorProductQuadratureRule::tensor_product_weight(
  int s1Node, int s2Node,
  int s1Ip, int s2Ip) const
{
  //surface integration
  const double Ls1 = scsEndLoc_[s1Node+1]-scsEndLoc_[s1Node];
  const double Ls2 = scsEndLoc_[s2Node+1]-scsEndLoc_[s2Node];
  const double isoparametricArea = Ls1 * Ls2;
  const double weight = isoparametricArea * weights_[s1Ip] * weights_[s2Ip];

  return weight;
}
//--------------------------------------------------------------------------
double
TensorProductQuadratureRule::tensor_product_weight(int s1Node, int s1Ip) const
{
  const double isoparametricLength = scsEndLoc_[s1Node+1]-scsEndLoc_[s1Node];
  const double weight = isoparametricLength * weights_[s1Ip];

  return weight;
}

}  // namespace nalu
} // namespace sierra
