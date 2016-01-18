/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/
#ifndef TensorProductQuadratureRule_h
#define TensorProductQuadratureRule_h

#include <string>
#include <vector>

namespace sierra{
namespace nalu{

class TensorProductQuadratureRule
{
public:
  TensorProductQuadratureRule(
    std::string type,
    int numQuad,
    std::vector<double>& scsLocs
  );
  ~TensorProductQuadratureRule() {};

  std::vector<double>& abscissae() { return abscissae_; };
  std::vector<double>& weights() { return weights_; };
  std::vector<double>& scsEndLoc() { return scsEndLoc_; };

  double abscissa(unsigned j) const { return abscissae_[j]; };
  double weight(unsigned j) const { return weights_[j]; };

  double gauss_point_location(
    int nodeOrdinal,
    int gaussPointOrdinal) const;

  double tensor_product_weight(
    int s1Node, int s2Node, int s3Node,
    int s1Ip, int s2Ip, int s3Ip) const;

  double tensor_product_weight(
    int s1Node, int s2Node,
    int s1Ip, int s2Ip) const;

  double tensor_product_weight(int s1Node, int s1Ip) const;

  double isoparametric_mapping(
    const double b,
    const double a,
    const double xi) const;

private:
  std::vector<double> abscissae_;
  std::vector<double> weights_;
  std::vector<double> scsEndLoc_;
};

} // namespace nalu
} // namespace Sierra

#endif
