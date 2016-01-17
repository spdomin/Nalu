#include <element_promotion/ElementDescription.h>

#include <element_promotion/FaceOperations.h>
#include <element_promotion/LagrangeBasis.h>
#include <element_promotion/QuadratureRule.h>
#include <element_promotion/TensorProductQuadratureRule.h>
#include <NaluEnv.h>
#include <nalu_make_unique.h>

#include <stk_util/environment/ReportHandler.hpp>

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <tuple>
#include <unordered_map>
#include <utility>

namespace sierra {
namespace nalu {

  //TODO(rcknaus): ElementDescription has become pretty bulky
  // separate out the CVFEM-specific stuff

std::unique_ptr<ElementDescription>
ElementDescription::create(int dimension,int order)
{
  bool symmWeights = false;

  if (dimension == 2 && order == 2) {
    std::vector<double> in_nodeLocs = { -1.0, 0.0, +1.0 };
    std::vector<double> in_scsLoc = { -std::sqrt(3.0)/3.0, std::sqrt(3.0)/3.0 };
    return make_unique<QuadMElementDescription>(in_nodeLocs,in_scsLoc);
  }
  if (dimension == 2 && order == 3 && symmWeights) {
    // symmetric mass matrix
    // I can't find a symmetric mass matrix for P > 3
    double xgll    = 0.4487053820572093009546164613323186035;
    double scsDist = 0.8347278713337825805263131558586123084;
    std::vector<double> in_nodeLocs = { -1.0, -xgll, +xgll, +1.0 };
    std::vector<double> in_scsLoc = { -scsDist, 0.0, scsDist };

    return make_unique<QuadMElementDescription>(in_nodeLocs,in_scsLoc);
  }

  if (dimension == 2 && order >= 3) {
     std::vector<double> lobattoNodes;
     std::vector<double> legendreSCSLocations;
     std::tie(lobattoNodes,std::ignore) = gauss_lobatto_legendre_rule(order+1);
     std::tie(legendreSCSLocations,std::ignore) = gauss_legendre_rule(order);

     return make_unique<QuadMElementDescription>(lobattoNodes,legendreSCSLocations);
  }

  if (dimension == 3 && order == 2) {
    double scsDist = std::sqrt(3.0)/3.0;
    std::vector<double> in_nodeLocs = {-1.0, 0.0, +1.0};
    std::vector<double> in_scsLoc = { -scsDist, +scsDist };
    return make_unique<HexMElementDescription>(in_nodeLocs, in_scsLoc);
  }

  if (dimension == 3 && order == 3 && symmWeights) {
    double xgll    = 0.4487053820572093009546164613323186035;
    double scsDist = 0.8347278713337825805263131558586123084;
    std::vector<double> in_nodeLocs = { -1.0, -xgll, +xgll, +1.0 };
    std::vector<double> in_scsLoc = { -scsDist, 0.0, scsDist };
    return make_unique<HexMElementDescription>(in_nodeLocs, in_scsLoc);
  }

  if (dimension == 3 && order >= 3) {
     std::vector<double> lobattoNodes;
     std::vector<double> legendreSCSLocations;
     std::tie(lobattoNodes,std::ignore) = gauss_lobatto_legendre_rule(order+1);
     std::tie(legendreSCSLocations,std::ignore) = gauss_legendre_rule(order);

     return make_unique<HexMElementDescription>(lobattoNodes,legendreSCSLocations);
  }

  throw std::runtime_error("Element type not implemented");
  return nullptr;
}

QuadMElementDescription::QuadMElementDescription(
  std::vector<double> in_nodeLocs, std::vector<double> in_scsLoc)
  : ElementDescription()
{
  nodeLocs = in_nodeLocs;
  scsLoc = in_scsLoc;
  ThrowRequire(nodeLocs.size()-1 == scsLoc.size());

  polyOrder = nodeLocs.size()-1;
  nodes1D = nodeLocs.size();
  nodesPerElement = nodes1D*nodes1D;
  dimension = 2;
  numQuad = (polyOrder % 2 == 0) ? polyOrder/2 + 1 : (polyOrder+1)/2;
  useGLLGLL = false;

  set_node_connectivity();
  set_subelement_connectivity();

  quadrature = make_unique<TensorProductQuadratureRule>("GaussLegendre", numQuad, scsLoc);
  basis = make_unique<LagrangeBasis>(inverseNodeMap, nodeLocs);
  basisBoundary = make_unique<LagrangeBasis>(inverseNodeMapBC,nodeLocs);
}
//--------------------------------------------------------------------------
void
QuadMElementDescription::set_node_connectivity()
{
  unsigned baseNodesPerElement = 4;
  unsigned nodeNumber = baseNodesPerElement;

  nodeMap.resize(nodes1D*nodes1D);
  auto nmap = [&] (unsigned i, unsigned j) -> unsigned&
  {
    return (nodeMap[i+nodes1D*j]);
  };

  struct EdgeInfo
  {
    int direction;
    int xloc;
    int yloc;
  };

  std::vector<std::pair<std::vector<size_t>, EdgeInfo>> baseEdgeInfo = {
      {{0,1}, {+1, 0,0}},
      {{1,2}, {+2, 1,0}},
      {{2,3}, {-1, 1,1}},
      {{3,0}, {-2, 0,1}}
  };

  int faceMap[4] = {0,1,2,3};

  unsigned jmax = nodes1D-1;
  std::vector<size_t> baseFaceNodes = {0, 1, 2, 3};

  nmap(0,0)       = 0;
  nmap(jmax,0)    = 1;
  nmap(jmax,jmax) = 2;
  nmap(0,jmax)    = 3;

  faceNodeMap.resize(4);
  unsigned faceOrdinal = 0;
  for (auto& baseEdge : baseEdgeInfo) {
    auto direction = baseEdge.second.direction;
    std::vector<size_t> nodesToAdd(nodes1D-2);
    for (unsigned j =0; j < nodes1D-2; ++j) {
      nodesToAdd[j] = nodeNumber;
      ++nodeNumber;
    }
    edgeNodeConnectivities.insert({nodesToAdd,baseEdge.first});

    auto reorderedNodes = nodesToAdd;
    if (direction < 0) {
      std::reverse(reorderedNodes.begin(), reorderedNodes.end());
    }

    unsigned il = (baseEdge.second.xloc == 1) ? jmax : 0;
    unsigned jl = (baseEdge.second.yloc == 1) ? jmax : 0;

    if (std::abs(direction) == 1) {
      for (unsigned j = 1; j < polyOrder; ++j) {
        nmap(j,jl) = reorderedNodes.at(j-1);
      }
    }
    else {
      for (unsigned j = 1; j < polyOrder; ++j) {
        nmap(il,j) = reorderedNodes.at(j-1);
      }
    }

    std::vector<std::vector<double>> locs(polyOrder-1);
    for (unsigned i = 0; i < polyOrder-1; ++i) {
      locs[i].push_back(nodeLocs[1+i]);
    }
    locationsForNewNodes.insert({nodesToAdd,locs});

    std::vector<size_t> faceNodes(nodes1D);
     if (std::abs(direction) == 1) {
       for (unsigned j = 0; j < nodes1D; ++j) {
         faceNodes[j] = nmap(j,jl);
       }
     }

     if (std::abs(direction) == 2) {
       for (unsigned j = 0; j < nodes1D; ++j) {
         faceNodes[j] = nmap(il,j);
       }
     }

     faceNodeMap[faceMap[faceOrdinal]] = faceNodes;
     ++faceOrdinal;
  }

  // reverse edges oriented in the negative direction in isoparametric space
  std::reverse(faceNodeMap[faceMap[2]].begin(), faceNodeMap[faceMap[2]].end());
  std::reverse(faceNodeMap[faceMap[3]].begin(), faceNodeMap[faceMap[3]].end());

  unsigned faceNodeNumber = nodeNumber;
  unsigned nodesLeft = (nodes1D*nodes1D) - faceNodeNumber;
  std::vector<size_t> faceNodesToAdd(nodesLeft);
  for (unsigned j = 0; j < nodesLeft;++j) {
    faceNodesToAdd[j] = faceNodeNumber;
    ++faceNodeNumber;
  }
  faceNodeConnectivities.insert({faceNodesToAdd,baseFaceNodes});

  for (unsigned j = 1; j < polyOrder; ++j) {
    for (unsigned i = 1; i < polyOrder; ++i) {
      nmap(i,j) = faceNodesToAdd.at((i-1)+(polyOrder-1)*(j-1));
    }
  }

  std::vector<std::vector<double>> locs((polyOrder-1)*(polyOrder-1));
  for (unsigned j = 0; j < polyOrder-1; ++j) {
    for (unsigned i = 0; i < polyOrder-1; ++i) {
      locs[i+(polyOrder-1)*j] = {nodeLocs[i+1],nodeLocs[j+1]};
    }
  }
  locationsForNewNodes.insert({faceNodesToAdd, locs});

  for (const auto& edgeNode : edgeNodeConnectivities) {
    addedConnectivities.insert(edgeNode);
  }

  for (const auto& faceNode : faceNodeConnectivities) {
    addedConnectivities.insert(faceNode);
  }

  nodeMapBC.resize(nodes1D);
  nodeMapBC[0] = 0;
  nodeMapBC[nodes1D-1] = 1;

  nodeNumber = 2;
  for (unsigned j = 1; j < polyOrder; ++j) {
    nodeMapBC[j] = nodeNumber;
    ++nodeNumber;
  }

  //inverse maps
  inverseNodeMap.resize(nodes1D*nodes1D);
  for (unsigned i = 0; i < nodes1D; ++i) {
    for (unsigned j = 0; j < nodes1D; ++j) {
      inverseNodeMap[tensor_product_node_map(i,j)] = {i, j};
    }
  }

  inverseNodeMapBC.resize(nodes1D);
  for (unsigned j = 0; j < nodes1D; ++j) {
    inverseNodeMapBC[tensor_product_node_map(j)] = { j };
  }

  sideOrdinalMap.resize(4);
  for (unsigned face_ordinal = 0; face_ordinal < 4; ++face_ordinal) {
    sideOrdinalMap[face_ordinal].resize(nodes1D);
    for (unsigned j = 0; j < nodes1D; ++j) {
      auto& ord = inverseNodeMapBC[j];
      sideOrdinalMap[face_ordinal][j] = faceNodeMap[face_ordinal][ord[0]];
    }
  }
}
//--------------------------------------------------------------------------
void
QuadMElementDescription::set_subelement_connectivity()
{
  subElementConnectivity.resize((nodes1D-1)*(nodes1D-1));
  for (unsigned j = 0; j < nodes1D-1; ++j) {
    for (unsigned i = 0; i < nodes1D-1; ++i) {
      subElementConnectivity[i+(nodes1D-1)*j] =
      {
          static_cast<size_t>(tensor_product_node_map(i,j)),
          static_cast<size_t>(tensor_product_node_map(i+1,j)),
          static_cast<size_t>(tensor_product_node_map(i+1,j+1)),
          static_cast<size_t>(tensor_product_node_map(i,j+1))
      };
    }
  }
}
//--------------------------------------------------------------------------
HexMElementDescription::HexMElementDescription(
  std::vector<double> in_nodeLocs, std::vector<double> in_scsLoc)
:  ElementDescription()
{
  scsLoc = in_scsLoc;
  nodeLocs = in_nodeLocs;
  nodes1D = nodeLocs.size();
  polyOrder = nodes1D-1;
  nodesPerElement = nodes1D*nodes1D*nodes1D;
  dimension = 3;
  numQuad = (polyOrder % 2 == 0) ? polyOrder/2 + 1 : (polyOrder+1)/2;
  useGLLGLL = false;

  set_node_connectivity();
  set_subelement_connectivity();

  quadrature = make_unique<TensorProductQuadratureRule>("GaussLegendre", numQuad, scsLoc);
  basis = make_unique<LagrangeBasis>(inverseNodeMap, nodeLocs);
  basisBoundary = make_unique<LagrangeBasis>(inverseNodeMapBC, nodeLocs);
}
//--------------------------------------------------------------------------
void
HexMElementDescription::set_subelement_connectivity()
{
  subElementConnectivity.resize((nodes1D-1)*(nodes1D-1)*(nodes1D-1));
  for (unsigned k = 0; k < nodes1D-1; ++k) {
    for (unsigned j = 0; j < nodes1D-1; ++j) {
      for (unsigned i = 0; i < nodes1D-1; ++i) {
        subElementConnectivity[i+(nodes1D-1)*(j+(nodes1D-1)*k)] =
        {
            static_cast<size_t>(tensor_product_node_map(i+0,j+0,k+0)),
            static_cast<size_t>(tensor_product_node_map(i+1,j+0,k+0)),
            static_cast<size_t>(tensor_product_node_map(i+1,j+0,k+1)),
            static_cast<size_t>(tensor_product_node_map(i+0,j+0,k+1)),
            static_cast<size_t>(tensor_product_node_map(i+0,j+1,k+0)),
            static_cast<size_t>(tensor_product_node_map(i+1,j+1,k+0)),
            static_cast<size_t>(tensor_product_node_map(i+1,j+1,k+1)),
            static_cast<size_t>(tensor_product_node_map(i+0,j+1,k+1))
        };
      }
    }
  }
}
//--------------------------------------------------------------------------
void
HexMElementDescription::set_node_connectivity()
{
  unsigned baseNodesPerElement = 8;
  unsigned nodeNumber = baseNodesPerElement;

  nodeMap.assign(nodes1D*nodes1D*nodes1D,0);
  auto nmap = [&] (unsigned i, unsigned j, unsigned k) -> unsigned&
  {
    return (nodeMap[i+nodes1D*(j+nodes1D*k)]);
  };

  struct EdgeInfo
  {
    int direction;
    std::vector<int> baseLoc;
  };

  std::vector<std::pair<std::vector<size_t>, EdgeInfo>>
  baseEdgeInfo = {
      {{0,1}, {+1, {0,0,0} }},
      {{1,2}, {+2, {1,0,0} }},
      {{2,3}, {-1, {1,1,0} }},
      {{3,0}, {-2, {0,1,0} }},
      {{0,4}, {+3, {0,0,0} }},
      {{1,5}, {+3, {1,0,0} }},
      {{2,6}, {+3, {1,1,0} }},
      {{3,7}, {+3, {0,1,0} }},
      {{4,5}, {+1, {0,0,1} }},
      {{5,6}, {+2, {1,0,1} }},
      {{6,7}, {-1, {1,1,1} }},
      {{7,4}, {-2, {0,1,1} }},
  };

  //add the base nodes to the map
  unsigned jmax = nodes1D-1;
  nmap(0,0,0)          = 0;
  nmap(jmax,0,0)       = 1;
  nmap(jmax,jmax,0)    = 2;
  nmap(0,jmax,0)       = 3;
  nmap(0,0,jmax)       = 4;
  nmap(jmax,0,jmax)    = 5;
  nmap(jmax,jmax,jmax) = 6;
  nmap(0,jmax,jmax)    = 7;

  std::vector<std::vector<size_t>> baseVolumeNodes = {
      {0,1,2,3,4,5,6,7}
  };

  unsigned nodes1DAdded = nodes1D-2;
  for (auto& baseEdge : baseEdgeInfo) {
    std::vector<size_t> nodesToAdd(nodes1DAdded);
    for (unsigned j =0; j < nodes1DAdded; ++j) {
      nodesToAdd[j] = nodeNumber;
      ++nodeNumber;
    }
    edgeNodeConnectivities.insert({nodesToAdd,baseEdge.first});

    const auto& edgeInfo = baseEdge.second;
    const auto direction = edgeInfo.direction;
    const auto& baseLoc = edgeInfo.baseLoc;

    auto reorderedNodes = nodesToAdd;
    if (direction < 0) {
      std::reverse(reorderedNodes.begin(), reorderedNodes.end());
    }

    unsigned il = (baseLoc[0] == 1) ? jmax : 0;
    unsigned jl = (baseLoc[1] == 1) ? jmax : 0;
    unsigned kl = (baseLoc[2] == 1) ? jmax : 0;

    switch(std::abs(direction))
    {
      case 1:
      {
        for (unsigned j = 1; j < polyOrder; ++j) {
          nmap(j,jl,kl) = reorderedNodes.at(j-1);
        }
        break;
      }
      case 2:
      {
        for (unsigned j = 1; j < polyOrder; ++j) {
          nmap(il,j,kl) = reorderedNodes.at(j-1);
        }
        break;
      }
      case 3:
      {
        for (unsigned j = 1; j < polyOrder; ++j) {
          nmap(il,jl,j) = reorderedNodes.at(j-1);
        }
        break;
      }
      default:
      {
        throw std::runtime_error("Invalid direction");
        break;
      }
    }
    std::vector<std::vector<double>> locs(polyOrder-1);
    for (unsigned i = 0; i < polyOrder-1; ++i) {
      locs[i].push_back(nodeLocs[1+i]);
    }
    locationsForNewNodes.insert({nodesToAdd,locs});
  }

  // volume nodes are inserted second to be consistent with exodus format for P=2
  // (likely for consistency with Hex20 elements), but awkward here
  unsigned volumeNodeNumber = nodeNumber;
  for (auto& baseVolume : baseVolumeNodes) {
    std::vector<size_t> volumeNodesToAdd(nodes1DAdded*nodes1DAdded*nodes1DAdded);
    for (auto& volumeNodeOrdinal : volumeNodesToAdd) {
       volumeNodeOrdinal = volumeNodeNumber;
      ++volumeNodeNumber;
    }
    volumeNodeConnectivities.insert({volumeNodesToAdd,baseVolume});

    for (unsigned k = 1; k < polyOrder; ++k) {
      for (unsigned j = 1; j < polyOrder; ++j) {
        for (unsigned i = 1; i < polyOrder; ++i) {
          nmap(i,j,k) =
            volumeNodesToAdd.at((i-1)+(polyOrder-1)*((j-1)+(polyOrder-1)*(k-1)));
        }
      }
    }

    std::vector<std::vector<double>> locs((polyOrder-1)*(polyOrder-1)*(polyOrder-1));
    for (unsigned k = 0; k < polyOrder-1; ++k) {
      for (unsigned j = 0; j < polyOrder-1; ++j) {
        for (unsigned i = 0; i < polyOrder-1; ++i) {
          locs[i+(polyOrder-1)*(j+(polyOrder-1)*k)] =
            {nodeLocs[i+1],nodeLocs[k+1], nodeLocs[j+1]};
        }
      }
    }
    locationsForNewNodes.insert({volumeNodesToAdd, locs});
  }

  struct FaceInfo
  {
    bool yreflected; // about y
    bool rotated;
    int xnormal;
    int ynormal;
    int znormal;
  };

  std::vector<std::pair<std::vector<size_t>, FaceInfo>> baseFaceInfo =
  {
      {{0, 3, 2, 1}, {false,true, 0,0,-1}},
      {{4, 5, 6, 7}, {false,false, 0,0,+1}},
      {{0, 4, 7, 3}, {false,true, -1,0,0}},
      {{1, 2, 6, 5}, {false,false, +1,0,0}},
      {{0, 1, 5, 4}, {false,false,  0,-1,0}},
      {{2, 3, 7, 6}, {true,false, 0,+1,0}}
  };

  int faceMap[6] = { 4, 5, 3, 1, 0, 2 };

  unsigned faceNodeNumber = volumeNodeNumber;
  for (auto& baseFace : baseFaceInfo) {
    std::vector<size_t> faceNodesToAdd(nodes1DAdded*nodes1DAdded);
    for (auto& faceNodeOrdinal : faceNodesToAdd) {
      faceNodeOrdinal = faceNodeNumber;
      ++faceNodeNumber;
    }
    faceNodeConnectivities.insert({faceNodesToAdd,baseFace.first});

    const auto& faceInfo = baseFace.second;
    const auto xnormal = faceInfo.xnormal;
    const auto ynormal = faceInfo.ynormal;
    const auto znormal = faceInfo.znormal;
    ThrowAssert(std::abs(xnormal)+std::abs(ynormal) + std::abs(znormal) == 1);

    std::vector<size_t> reorderedFaceNodes = faceNodesToAdd;
    bool isReflected = faceInfo.yreflected;
    if (isReflected) {
      reorderedFaceNodes = flip_x<size_t>(faceNodesToAdd, polyOrder-1);
    }

    bool isRotated = faceInfo.rotated;
    if (isRotated) {
      reorderedFaceNodes = transpose_ordinals<size_t>(faceNodesToAdd, polyOrder-1);
    }

    if (xnormal != 0) {
      const int il = (xnormal > 0) ? jmax : 0;
      for (unsigned j = 1; j < polyOrder; ++j) {
        for (unsigned i = 1; i < polyOrder; ++i) {
          nmap(il,i,j) = reorderedFaceNodes.at((i-1)+(polyOrder-1)*(j-1));
        }
      }
    }

    if (ynormal != 0) {
      const int jl = (ynormal > 0) ? jmax : 0;
      for (unsigned j = 1; j < polyOrder; ++j) {
        for (unsigned i = 1; i < polyOrder; ++i) {
          nmap(i,jl,j) = reorderedFaceNodes.at((i-1)+(polyOrder-1)*(j-1));
        }
      }
    }

    if (znormal != 0) {
      const int kl = (znormal > 0) ? jmax : 0;
      for (unsigned j = 1; j < polyOrder; ++j) {
        for (unsigned i = 1; i < polyOrder; ++i) {
          nmap(i,j,kl) = reorderedFaceNodes.at((i-1)+(polyOrder-1)*(j-1));
        }
      }
    }

    std::vector<std::vector<double>> locs((polyOrder-1)*(polyOrder-1));
    for (unsigned j = 0; j < polyOrder-1; ++j) {
      for (unsigned i = 0; i < polyOrder-1; ++i) {
        locs[i+(polyOrder-1)*j] = {nodeLocs[i+1],nodeLocs[j+1]};
      }
    }
    locationsForNewNodes.insert({faceNodesToAdd, locs});
  }

  faceNodeMap.resize(6);
  unsigned faceOrdinal = 0;
  for (auto& baseFace : baseFaceInfo) {
    const auto& faceInfo = baseFace.second;
    const auto xnormal = faceInfo.xnormal;
    const auto ynormal = faceInfo.ynormal;
    const auto znormal = faceInfo.znormal;
    ThrowAssert(std::abs(xnormal)+std::abs(ynormal) + std::abs(znormal) == 1);

    std::vector<size_t> faceNodes(nodes1D*nodes1D);
    if (xnormal != 0) {
      const int il = (xnormal > 0) ? jmax : 0;
      for (unsigned j = 0; j < nodes1D; ++j) {
        for (unsigned i = 0; i < nodes1D; ++i) {
          faceNodes[i+nodes1D*j] = nmap(il,i,j);
        }
      }
    }

    if (ynormal != 0) {
      const int jl = (ynormal > 0) ? jmax : 0;
      for (unsigned j = 0; j < nodes1D; ++j) {
        for (unsigned i = 0; i < nodes1D; ++i) {
          faceNodes[i+nodes1D*j] = nmap(i,jl,j);
        }
      }
    }

    if (znormal != 0) {
      const int kl = (znormal > 0) ? jmax : 0;
      for (unsigned j = 0; j < nodes1D; ++j) {
        for (unsigned i = 0; i < nodes1D; ++i) {
          faceNodes[i+nodes1D*j] = nmap(i,j,kl);
        }
      }
    }

    faceNodeMap[faceMap[faceOrdinal]] = faceNodes;
    ++faceOrdinal;
  }

  // transform face ordinals depending on how they're arranged
  // in isoparametric coordinates
  faceNodeMap[2] = flip_x<size_t>(faceNodeMap[2], nodes1D);
  faceNodeMap[3] = transpose_ordinals<size_t>(faceNodeMap[3], nodes1D);
  faceNodeMap[4] = transpose_ordinals<size_t>(faceNodeMap[4], nodes1D);


  for (const auto& edgeNode : edgeNodeConnectivities) {
    addedConnectivities.insert(edgeNode);
  }

  for (const auto& volumeNode : volumeNodeConnectivities) {
    addedConnectivities.insert(volumeNode);
  }

  for (const auto& faceNode : faceNodeConnectivities) {
    addedConnectivities.insert(faceNode);
  }

  nodeMapBC = QuadMElementDescription(nodeLocs,scsLoc).nodeMap;

  //inverse maps
  inverseNodeMap.resize(nodes1D*nodes1D*nodes1D);
  for (unsigned i = 0; i < nodes1D; ++i) {
    for (unsigned j = 0; j < nodes1D; ++j) {
      for (unsigned k = 0; k < nodes1D; ++k) {
        inverseNodeMap[tensor_product_node_map(i,j,k)] = {i, j, k};
      }
    }
  }

  inverseNodeMapBC.resize(nodes1D*nodes1D);
  for (unsigned i = 0; i < nodes1D; ++i) {
    for (unsigned j = 0; j < nodes1D; ++j) {
      inverseNodeMapBC[tensor_product_node_map_bc(i,j)] = { i,j };
    }
  }

  sideOrdinalMap.resize(6);
  for (unsigned face_ordinal = 0; face_ordinal < 6; ++face_ordinal) {
    sideOrdinalMap[face_ordinal].resize(nodes1D*nodes1D);
    for (unsigned j = 0; j < nodes1D*nodes1D; ++j) {
      auto& ords = inverseNodeMapBC[j];
      sideOrdinalMap[face_ordinal][j] = faceNodeMap[face_ordinal][ords[0]+nodes1D*ords[1]];
    }
  }
}

} // namespace nalu
}  // namespace sierra
