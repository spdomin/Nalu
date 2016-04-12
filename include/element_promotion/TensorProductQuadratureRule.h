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

  const std::vector<double>& abscissae() const { return abscissae_; };
  const std::vector<double>& weights() const { return weights_; };
  const std::vector<double>& scsEndLoc() const { return scsEndLoc_; };

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
  int numQuad_;
  std::vector<double> abscissae_;
  std::vector<double> weights_;
  std::vector<double> scsEndLoc_;
  bool useSGL_;
};

} // namespace nalu
} // namespace Sierra

#endif
