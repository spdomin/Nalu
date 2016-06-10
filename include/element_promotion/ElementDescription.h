/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level nalu      */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#ifndef ElementDescription_h
#define ElementDescription_h

#include <element_promotion/LagrangeBasis.h>
#include <element_promotion/TensorProductQuadratureRule.h>

#include <stddef.h>
#include <map>
#include <memory>
#include <vector>

namespace sierra {
namespace nalu {

typedef std::map<std::vector<size_t>, std::vector<size_t>> AddedConnectivityOrdinalMap;
typedef std::map<std::vector<size_t>, std::vector<std::vector<double>>> AddedNodeLocationsMap;
typedef std::vector<std::vector<size_t>> SubElementConnectivity;

struct ElementDescription
{
public:
  static std::unique_ptr<ElementDescription> create(
    int dimension,
    int order,
    std::string quadType = "GaussLegendre",
    bool useReducedGeometricBasis = false
  );
  virtual ~ElementDescription() = default;

  int tensor_product_node_map(int i, int j, int k) const
  {
    return nodeMap[i+nodes1D*(j+nodes1D*k)];
  }

  int tensor_product_node_map(int i, int j) const
  {
    return nodeMap[i+nodes1D*j];
  }

  int tensor_product_node_map_bc(int i, int j) const
  {
    return nodeMapBC[i+nodes1D*j];
  }

  int tensor_product_node_map_bc(int i) const
  {
    return nodeMapBC[i];
  }

  int tensor_index(int nodeNumber) const
  {
    const auto& ords = inverseNodeMap[nodeNumber];
    if (dimension == 2) {
      return (ords[0] + nodes1D*ords[1]);
    }
    return (ords[0] + nodes1D*(ords[1]+nodes1D*ords[2]));
  }

  const int* side_ordinals_for_face(int face_ordinal) const {
    return sideOrdinalMap[face_ordinal].data();
  }

  int side_ordinal(int face_ordinal, int node_ordinal_index) const {
    return sideOrdinalMap[face_ordinal][node_ordinal_index];
  }

  double gauss_point_location(
    int nodeOrdinal,
    int gaussPointOrdinal) const
  {
    return quadrature->gauss_point_location(nodeOrdinal,gaussPointOrdinal);
  }

  double tensor_product_weight(
    int s1Node, int s2Node, int s3Node,
    int s1Ip, int s2Ip, int s3Ip) const
  {
    return quadrature->tensor_product_weight(s1Node,s2Node,s3Node,s1Ip,s2Ip,s3Ip);
  }

  double tensor_product_weight(
    int s1Node, int s2Node,
    int s1Ip, int s2Ip) const
  {
    return quadrature->tensor_product_weight(s1Node,s2Node, s1Ip,s2Ip);
  }

  double tensor_product_weight(int s1Node, int s1Ip) const
  {
    return quadrature->tensor_product_weight(s1Node, s1Ip);
  }

  std::vector<double>
  eval_basis_weights(const std::vector<double>& intgLoc) const
  {
     return basis->eval_basis_weights(intgLoc);
  }

  std::vector<double>
  eval_deriv_weights(const std::vector<double>& intgLoc) const
  {
    return basis->eval_deriv_weights(intgLoc);
  }

  std::vector<double> scsLoc;
  std::vector<double> nodeLocs;

  size_t dimension;
  size_t nodes1D;
  size_t nodesPerElement;
  size_t nodesPerFace;
  size_t nodesInBaseElement;
  size_t nodesPerSubElement;
  AddedConnectivityOrdinalMap addedConnectivities;
  AddedConnectivityOrdinalMap edgeNodeConnectivities;
  AddedConnectivityOrdinalMap faceNodeConnectivities;
  AddedConnectivityOrdinalMap volumeNodeConnectivities;
  AddedNodeLocationsMap locationsForNewNodes;
  SubElementConnectivity subElementConnectivity;

  std::string quadType;
  bool useReducedGeometricBasis;
  unsigned polyOrder;
  unsigned numQuad;
  std::unique_ptr<TensorProductQuadratureRule> quadrature;
  std::unique_ptr<LagrangeBasis> basis;
  std::unique_ptr<LagrangeBasis> basisBoundary;
  std::unique_ptr<LagrangeBasis> linearBasis;
  std::vector<unsigned> nodeMap;
  std::vector<unsigned> nodeMapBC;
  std::vector<std::vector<unsigned>> inverseNodeMap;
  std::vector<std::vector<unsigned>> inverseNodeMapBC;
  std::vector<std::vector<size_t>> faceNodeMap;
  std::vector<std::vector<int>> sideOrdinalMap;
protected:
  ElementDescription() = default;
};

struct QuadMElementDescription: public ElementDescription
{
public:
  QuadMElementDescription(
    std::vector<double> nodeLocs,
    std::vector<double> scsLoc,
    std::string quadType,
    bool useReducedGeometricBasis
  );
private:
  void set_node_connectivity();
  void set_subelement_connectivity();
};

struct HexMElementDescription: public ElementDescription
{
public:
  HexMElementDescription(
    std::vector<double> nodeLocs,
    std::vector<double> scsLoc,
    std::string quadType,
    bool useReducedGeometricBasis
  );
private:
  void set_node_connectivity();
  void set_subelement_connectivity();
};

} // namespace nalu
} // namespace Sierra

#endif
