/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporatlion.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/
#include <element_promotion/MasterElementHO.h>
#include <element_promotion/ElementDescription.h>
#include <element_promotion/LagrangeBasis.h>

#include <master_element/MasterElement.h>
#include <FORTRAN_Proto.h>

#include <array>
#include <cmath>
#include <memory>
#include <stdexcept>

namespace sierra{
namespace nalu{

HigherOrderHexSCV::HigherOrderHexSCV(const ElementDescription& elem)
  : MasterElement(),
    elem_(elem)
{
  // TODO(rcknaus): need to pick a modal basis (Hex8 or Hex27, really)
  // and use it only for geometric calculations (i.e. volume & area)
  // otherwise, computing the mass matrix scales far too viciously with p, O(p^9)
  nDim_ = elem_.dimension;
  nodesPerElement_ = elem_.nodesPerElement;

  // set up integration rule and relevant maps for scvs
  set_interior_info();

  // compute and save shape functions and derivatives at ips
  shapeFunctions_ = elem_.eval_basis_weights(intgLoc_);
  shapeDerivs_ = elem_.eval_deriv_weights(intgLoc_);
}

//--------------------------------------------------------------------------
//-------- set_interior_info -----------------------------------------------
//--------------------------------------------------------------------------
void
HigherOrderHexSCV::set_interior_info()
{
  //1D integration rule per sub-control volume
  numIntPoints_ = (elem_.nodes1D * elem_.nodes1D  * elem_.nodes1D)
                * ( elem_.numQuad * elem_.numQuad * elem_.numQuad);

  // define ip node mappings
  ipNodeMap_.resize(numIntPoints_);
  intgLoc_.resize(numIntPoints_*nDim_);
  intgLocShift_.resize(numIntPoints_*nDim_);
  ipWeight_.resize(numIntPoints_);

  // tensor product nodes x tensor product quadrature
  int vector_index = 0; int scalar_index = 0;
  for (unsigned n = 0; n < elem_.nodes1D; ++n) {
    for (unsigned m = 0; m < elem_.nodes1D; ++m) {
      for (unsigned l = 0; l < elem_.nodes1D; ++l) {

        // current node number
        const int nodeNumber = elem_.tensor_product_node_map(l,m,n);

        //tensor-product quadrature for a particular sub-cv
        for (unsigned k = 0; k < elem_.numQuad; ++k) {
          for (unsigned j = 0; j < elem_.numQuad; ++j) {
            for (unsigned i = 0; i < elem_.numQuad; ++i) {
              //integration point location
              intgLoc_[vector_index]     = elem_.gauss_point_location(l,i);
              intgLoc_[vector_index + 1] = elem_.gauss_point_location(m,j);
              intgLoc_[vector_index + 2] = elem_.gauss_point_location(n,k);

              //weight
              ipWeight_[scalar_index] = elem_.tensor_product_weight(l,m,n,i,j,k);

              //sub-control volume association
              ipNodeMap_[scalar_index] = nodeNumber;

              // increment indices
              ++scalar_index;
              vector_index += nDim_;
            }
          }
        }
      }
    }
  }
}
//--------------------------------------------------------------------------
void
HigherOrderHexSCV::shape_fcn(double *shpfc)
{
  int numShape = shapeFunctions_.size();
  for (int j = 0; j < numShape; ++j) {
    shpfc[j] = shapeFunctions_[j];
  }
}
//--------------------------------------------------------------------------
const int *
HigherOrderHexSCV::ipNodeMap(
  int /*ordinal*/)
{
  // define scv->node mappings
  return &ipNodeMap_[0];
}
//--------------------------------------------------------------------------
void HigherOrderHexSCV::determinant(
  const int nelem,
  const double *coords,
  double *volume,
  double *error)
{
  *error = 0.0;
  for (int k = 0; k < nelem; ++k) {
    const int scalar_elem_offset = numIntPoints_ * k;
    const int coord_elem_offset = nDim_ * nodesPerElement_ * k;
    for (int ip = 0; ip < numIntPoints_; ++ip) {
      const int grad_offset = nDim_ * nodesPerElement_ * ip;

      //weighted jacobian determinant
      const double det_j = jacobian_determinant(
        &coords[coord_elem_offset],
        &shapeDerivs_[grad_offset]
      );

      //apply weight and store to volume
      volume[scalar_elem_offset + ip] = ipWeight_[ip] * det_j;

      //flag error
      if (det_j < 0.0) {
        *error = 1.0;
      }
    }
  }
}
//--------------------------------------------------------------------------
double
HigherOrderHexSCV::jacobian_determinant(
  const double *elemNodalCoords,
  const double *shapeDerivs) const
{
  double dx_ds1 = 0.0;  double dx_ds2 = 0.0; double dx_ds3 = 0.0;
  double dy_ds1 = 0.0;  double dy_ds2 = 0.0; double dy_ds3 = 0.0;
  double dz_ds1 = 0.0;  double dz_ds2 = 0.0; double dz_ds3 = 0.0;
  for (int node = 0; node < nodesPerElement_; ++node) {
    const int vector_offset = nDim_ * node;

    const double xCoord = elemNodalCoords[vector_offset+0];
    const double yCoord = elemNodalCoords[vector_offset+1];
    const double zCoord = elemNodalCoords[vector_offset+2];

    const double dn_ds1 = shapeDerivs[vector_offset+0];
    const double dn_ds2 = shapeDerivs[vector_offset+1];
    const double dn_ds3 = shapeDerivs[vector_offset+2];

    dx_ds1 += dn_ds1 * xCoord;
    dx_ds2 += dn_ds2 * xCoord;
    dx_ds3 += dn_ds3 * xCoord;

    dy_ds1 += dn_ds1 * yCoord;
    dy_ds2 += dn_ds2 * yCoord;
    dy_ds3 += dn_ds3 * yCoord;

    dz_ds1 += dn_ds1 * zCoord;
    dz_ds2 += dn_ds2 * zCoord;
    dz_ds3 += dn_ds3 * zCoord;
  }

  const double det_j = dx_ds1 * ( dy_ds2 * dz_ds3 - dz_ds2 * dy_ds3 )
                     + dy_ds1 * ( dz_ds2 * dx_ds3 - dx_ds2 * dz_ds3 )
                     + dz_ds1 * ( dx_ds2 * dy_ds3 - dy_ds2 * dx_ds3 );

  return det_j;
}
//--------------------------------------------------------------------------
HigherOrderHexSCS::HigherOrderHexSCS(const ElementDescription& elem)
: MasterElement(),
  elem_(elem)
{
  nDim_ = elem_.dimension;
  nodesPerElement_ = elem_.nodesPerElement;

  // set up integration rule and relevant maps on scs
  set_interior_info();

  // set up integration rule and relevant maps on faces
  set_boundary_info();

  shapeFunctions_ = elem_.eval_basis_weights(intgLoc_);
  shapeDerivs_ = elem_.eval_deriv_weights(intgLoc_);
  expFaceShapeDerivs_ = elem_.eval_deriv_weights(intgExpFace_);
}
//--------------------------------------------------------------------------
void
HigherOrderHexSCS::set_interior_info()
{
  const int surfacesPerDirection = elem_.nodes1D - 1;
  const int ipsPerSurface = (elem_.numQuad*elem_.numQuad)*(elem_.nodes1D*elem_.nodes1D);
  const int numSurfaces = surfacesPerDirection * nDim_;

  numIntPoints_ = numSurfaces*ipsPerSurface;
  const int numVectorPoints = numIntPoints_*nDim_;

  // define L/R mappings
  lrscv_.resize(2*numIntPoints_);

  // standard integration location
  intgLoc_.resize(numVectorPoints);

  // shifted
  intgLocShift_.resize(numVectorPoints);

  // Save quadrature weight and directionality information
  ipInfo_.resize(numIntPoints_);

  // specify integration point locations in a dimension-by-dimension manner
  //u direction: bottom-top (0-1)
  int vector_index = 0; int lrscv_index = 0; int scalar_index = 0;
  for (int m = 0; m < surfacesPerDirection; ++m) {
    for (unsigned l = 0; l < elem_.nodes1D; ++l) {
      for (unsigned k = 0; k < elem_.nodes1D; ++k) {

        int leftNode; int rightNode;
        if (m % 2 == 0) {
          leftNode = elem_.tensor_product_node_map(k,l,m);
          rightNode = elem_.tensor_product_node_map(k,l,m+1);
        }
        else {
          leftNode = elem_.tensor_product_node_map(k,l,m+1);
          rightNode = elem_.tensor_product_node_map(k,l,m);
        }

        for (unsigned j = 0; j < elem_.numQuad; ++j) {
          for (unsigned i = 0; i < elem_.numQuad; ++i) {
            lrscv_[lrscv_index]     = leftNode;
            lrscv_[lrscv_index + 1] = rightNode;

            intgLoc_[vector_index]     = elem_.gauss_point_location(k,i);
            intgLoc_[vector_index + 1] = elem_.gauss_point_location(l,j);
            intgLoc_[vector_index + 2] = elem_.scsLoc.at(m);

            //compute the quadrature weight
            ipInfo_[scalar_index].weight = std::pow(-1.0,m+1) * elem_.tensor_product_weight(k,l,i,j);

            //direction
            ipInfo_[scalar_index].direction = Jacobian::U_DIRECTION;

            ++scalar_index;
            lrscv_index += 2;
            vector_index += nDim_;
          }
        }
      }
    }
  }

  //t direction: front-back (2-3)
  for (int m = 0; m < surfacesPerDirection; ++m) {
    for (unsigned l = 0; l < elem_.nodes1D; ++l) {
      for (unsigned k = 0; k < elem_.nodes1D; ++k) {

        int leftNode; int rightNode;
        if (m % 2 == 0) {
          leftNode = elem_.tensor_product_node_map(k,m,l);
          rightNode = elem_.tensor_product_node_map(k,m+1,l);
        }
        else {
          leftNode = elem_.tensor_product_node_map(k,m+1,l);
          rightNode = elem_.tensor_product_node_map(k,m,l);
        }

        for (unsigned j = 0; j < elem_.numQuad; ++j) {
          for (unsigned i = 0; i < elem_.numQuad; ++i) {
            lrscv_[lrscv_index]     = leftNode;
            lrscv_[lrscv_index + 1] = rightNode;

            intgLoc_[vector_index]     = elem_.gauss_point_location(k,i);
            intgLoc_[vector_index + 1] = elem_.scsLoc[m];
            intgLoc_[vector_index + 2] = elem_.gauss_point_location(l,j);

            //compute the quadrature weight
            ipInfo_[scalar_index].weight = std::pow(-1.0,m+1) * elem_.tensor_product_weight(k,l,i,j);

            //direction
            ipInfo_[scalar_index].direction = Jacobian::T_DIRECTION;

            ++scalar_index;
            lrscv_index += 2;
            vector_index += nDim_;
          }
        }
      }
    }
  }

  //s direction: left-right (4-5)
  for (int m = 0; m < surfacesPerDirection; ++m) {
    for (unsigned l = 0; l < elem_.nodes1D; ++l) {
      for (unsigned k = 0; k < elem_.nodes1D; ++k) {

        int leftNode; int rightNode;
        if (m % 2 == 0) {
          leftNode = elem_.tensor_product_node_map(m,k,l);
          rightNode = elem_.tensor_product_node_map(m+1,k,l);
        }
        else {
          leftNode = elem_.tensor_product_node_map(m+1,k,l);
          rightNode = elem_.tensor_product_node_map(m,k,l);
        }

        for (unsigned j = 0; j < elem_.numQuad; ++j) {
          for (unsigned i = 0; i < elem_.numQuad; ++i) {
            lrscv_[lrscv_index]     = leftNode;
            lrscv_[lrscv_index + 1] = rightNode;

            intgLoc_[vector_index]     = elem_.scsLoc[m];
            intgLoc_[vector_index + 1] = elem_.gauss_point_location(k,i);
            intgLoc_[vector_index + 2] = elem_.gauss_point_location(l,j);

            //compute the quadrature weight
            ipInfo_[scalar_index].weight = std::pow(-1.0, m) * elem_.tensor_product_weight(k,l,i,j);

            //direction
            ipInfo_[scalar_index].direction = Jacobian::S_DIRECTION;

            ++scalar_index;
            lrscv_index += 2;
            vector_index += nDim_;
          }
        }
      }
    }
  }
}

//--------------------------------------------------------------------------
//-------- set_boundary_info -----------------------------------------------
//--------------------------------------------------------------------------
void
HigherOrderHexSCS::set_boundary_info()
{
  const int numFaces = 2 * nDim_;
  int numQuad = static_cast<int>(elem_.numQuad);
  int nodes1D = static_cast<int>(elem_.nodes1D);
  const int nodesPerFace = nodes1D * nodes1D;
  ipsPerFace_ = nodesPerFace * (numQuad * numQuad);
  int surfacesPerDirection = nodes1D - 1;
  const int numFaceIps = numFaces * ipsPerFace_;

  oppFace_.resize(numFaceIps);
  ipNodeMap_.resize(numFaceIps);
  oppNode_.resize(numFaceIps);
  intgExpFace_.resize(numFaceIps*nDim_);

  // tensor-product style access to the map
  auto face_node_number = [&] (int i, int j, int faceOrdinal)
  {
    return elem_.faceNodeMap[faceOrdinal][i + nodes1D * j];
  };

  // map face ip ordinal to nearest sub-control surface ip ordinal
  // sub-control surface renumbering
  const std::vector<int> faceToSurface = {
      surfacesPerDirection,     // nearest scs face to t=-1.0
      3*surfacesPerDirection-1, // nearest scs face to s=+1.0, the last face
      2*surfacesPerDirection-1, // nearest scs face to t=+1.0
      2*surfacesPerDirection,   // nearest scs face to s=-1.0
      0,                        // nearest scs face to u=-1.0, the first face
      surfacesPerDirection-1    // nearest scs face to u=+1.0, the first face
  };
  auto opp_face_map = [&] ( int k, int l, int i, int j, int face_index)
  {
    int face_offset = faceToSurface[face_index] * ipsPerFace_;

    int node_index = k + nodes1D * l;
    int node_offset = node_index * (numQuad * numQuad);

    int ip_index = face_offset+node_offset+i+numQuad*j;

    return ip_index;
   };

  // location of the faces in the correct order
  const std::vector<double> faceLoc = {-1.0, +1.0, +1.0, -1.0, -1.0, +1.0};

  // Set points face-by-face
  int vector_index = 0; int scalar_index = 0; int faceOrdinal = 0;

  // front face: t = -1.0: counter-clockwise
  faceOrdinal = 0;
  for (int l = 0; l < nodes1D; ++l) {
    for (int k = 0; k < nodes1D; ++k) {
      const int nearNode = face_node_number(k,l,faceOrdinal);
      int oppNode = elem_.tensor_product_node_map(k,1,l);

      //tensor-product quadrature for a particular sub-cv
      for (int j = 0; j < numQuad; ++j) {
        for (int i = 0; i < numQuad; ++i) {
          // set maps
          ipNodeMap_[scalar_index] = nearNode;
          oppNode_[scalar_index] = oppNode;
          oppFace_[scalar_index] = opp_face_map(k,l,i,j,faceOrdinal);

          //integration point location
          intgExpFace_[vector_index]     = intgLoc_[oppFace_[scalar_index]*nDim_+0];
          intgExpFace_[vector_index + 1] = faceLoc[faceOrdinal];
          intgExpFace_[vector_index + 2] = intgLoc_[oppFace_[scalar_index]*nDim_+2];

          // increment indices
          ++scalar_index;
          vector_index += nDim_;
        }
      }
    }
  }

  // right face: s = +1.0: counter-clockwise
  faceOrdinal = 1;
  for (int l = 0; l < nodes1D; ++l) {
    for (int k = 0; k < nodes1D; ++k) {
      const int nearNode = face_node_number(k,l,faceOrdinal);
      int oppNode = elem_.tensor_product_node_map(nodes1D-2,k,l);

      //tensor-product quadrature for a particular sub-cv
      for (int j = 0; j < numQuad; ++j) {
        for (int i = 0; i < numQuad; ++i) {
          // set maps
          ipNodeMap_[scalar_index] = nearNode;
          oppNode_[scalar_index] = oppNode;
          oppFace_[scalar_index] = opp_face_map(k,l,i,j,faceOrdinal);

          //integration point location
          intgExpFace_[vector_index]     = faceLoc[faceOrdinal];
          intgExpFace_[vector_index + 1] = intgLoc_[oppFace_[scalar_index]*nDim_+1];
          intgExpFace_[vector_index + 2] = intgLoc_[oppFace_[scalar_index]*nDim_+2];

          // increment indices
          ++scalar_index;
          vector_index += nDim_;
        }
      }
    }
  }

  // back face: t = +1.0: s-direction reversed
  faceOrdinal = 2;
  for (int l = 0; l < nodes1D; ++l) {
    for (int k = nodes1D-1; k >= 0; --k) {
      const int nearNode = face_node_number(k,l,faceOrdinal);
      int oppNode = elem_.tensor_product_node_map(k,nodes1D-2,l);

      //tensor-product quadrature for a particular sub-cv
      for (int j = 0; j < numQuad; ++j) {
        for (int i = numQuad-1; i >= 0; --i) {
          // set maps
          ipNodeMap_[scalar_index] = nearNode;
          oppNode_[scalar_index] = oppNode;
          oppFace_[scalar_index] = opp_face_map(k,l,i,j,faceOrdinal);

          //integration point location
          intgExpFace_[vector_index]     = intgLoc_[oppFace_[scalar_index]*nDim_+0];
          intgExpFace_[vector_index + 1] = faceLoc[faceOrdinal];
          intgExpFace_[vector_index + 2] = intgLoc_[oppFace_[scalar_index]*nDim_+2];

          // increment indices
          ++scalar_index;
          vector_index += nDim_;
        }
      }
    }
  }

  //left face: x = -1.0 swapped t and u
  faceOrdinal = 3;
  for (int l = 0; l < nodes1D; ++l) {
    for (int k = 0; k < nodes1D; ++k) {
      const int nearNode = face_node_number(l,k,faceOrdinal);
      int oppNode = elem_.tensor_product_node_map(1,l,k);

      //tensor-product quadrature for a particular sub-cv
      for (int j = 0; j < numQuad; ++j) {
        for (int i = 0; i < numQuad; ++i) {
          // set maps
          ipNodeMap_[scalar_index] = nearNode;
          oppNode_[scalar_index]   = oppNode;
          oppFace_[scalar_index]   = opp_face_map(l,k,j,i,faceOrdinal);

          //integration point location
          intgExpFace_[vector_index]     = faceLoc[faceOrdinal];
          intgExpFace_[vector_index + 1] = intgLoc_[oppFace_[scalar_index]*nDim_+1];
          intgExpFace_[vector_index + 2] = intgLoc_[oppFace_[scalar_index]*nDim_+2];

          // increment indices
          ++scalar_index;
          vector_index += nDim_;
        }
      }
    }
  }

  //bottom face: u = -1.0: swapped s and t
  faceOrdinal = 4;
  for (int l = 0; l < nodes1D; ++l) {
    for (int k = 0; k < nodes1D; ++k) {
      const int nearNode = face_node_number(l,k,faceOrdinal);
      int oppNode = elem_.tensor_product_node_map(l,k,1);

      //tensor-product quadrature for a particular sub-cv
      for (int j = 0; j < numQuad; ++j) {
        for (int i = 0; i < numQuad; ++i) {
          // set maps
          ipNodeMap_[scalar_index] = nearNode;
          oppNode_[scalar_index] = oppNode;
          oppFace_[scalar_index] = opp_face_map(l,k,j,i,faceOrdinal);

          //integration point location
          intgExpFace_[vector_index]     = intgLoc_[oppFace_[scalar_index]*nDim_+0];
          intgExpFace_[vector_index + 1] = intgLoc_[oppFace_[scalar_index]*nDim_+1];
          intgExpFace_[vector_index + 2] = faceLoc[faceOrdinal];

          // increment indices
          ++scalar_index;
          vector_index += nDim_;
        }
      }
    }
  }

  //top face: u = +1.0: counter-clockwise
  faceOrdinal = 5;
  for (int l = 0; l < nodes1D; ++l) {
    for (int k = 0; k < nodes1D; ++k) {
      const int nearNode = face_node_number(k,l,faceOrdinal);
      int oppNode = elem_.tensor_product_node_map(k,l,nodes1D-2);

      //tensor-product quadrature for a particular sub-cv
      for (int j = 0; j < numQuad; ++j) {
        for (int i = 0; i < numQuad; ++i) {
          // set maps
          ipNodeMap_[scalar_index] = nearNode;
          oppNode_[scalar_index] = oppNode;
          oppFace_[scalar_index] = opp_face_map(k,l,i,j,faceOrdinal);

          //integration point location
          intgExpFace_[vector_index]     = intgLoc_[oppFace_[scalar_index]*nDim_+0];
          intgExpFace_[vector_index + 1] = intgLoc_[oppFace_[scalar_index]*nDim_+1];
          intgExpFace_[vector_index + 2] = faceLoc[faceOrdinal];

          // increment indices
          ++scalar_index;
          vector_index += nDim_;
        }
      }
    }
  }
  //throw std::runtime_error("check");
}
//--------------------------------------------------------------------------
void
HigherOrderHexSCS::shape_fcn(double* shpfc)
{
  int numShape = shapeFunctions_.size();
  for (int j = 0; j < numShape; ++j) {
    shpfc[j] = shapeFunctions_[j];
  }
}
//--------------------------------------------------------------------------
const int *
HigherOrderHexSCS::adjacentNodes()
{
  // define L/R mappings
  return &lrscv_[0];
}
//--------------------------------------------------------------------------
const int *
HigherOrderHexSCS::ipNodeMap(
  int ordinal)
{
  // define ip->node mappings for each face (ordinal);
  return &ipNodeMap_[ordinal*ipsPerFace_];
}
//--------------------------------------------------------------------------
int
HigherOrderHexSCS::opposingNodes(
  const int ordinal,
  const int node)
{
  return oppNode_[ordinal*ipsPerFace_+node];
}
//--------------------------------------------------------------------------
int
HigherOrderHexSCS::opposingFace(
  const int ordinal,
  const int node)
{
  return oppFace_[ordinal*ipsPerFace_+node];
}
//--------------------------------------------------------------------------
void
HigherOrderHexSCS::determinant(
  const int nelem,
  const double *coords,
  double *areav,
  double *error)
{
  //returns the normal vector x_t x x_u for constant s curves
  //returns the normal vector x_u x x_s for constant t curves
  //returns the normal vector x_s x x_t for constant u curves
  std::array<double,3> areaVector;

  for (int k = 0; k < nelem; ++k) {
    const int coord_elem_offset = nDim_ * nodesPerElement_ * k;
    const int vector_elem_offset = nDim_ * numIntPoints_ * k;

    for (int ip = 0; ip < numIntPoints_; ++ip) {
      const int grad_offset = nDim_ * nodesPerElement_ * ip;
      const int offset = nDim_ * ip + vector_elem_offset;

      //compute area vector for this ip
      area_vector( ipInfo_[ip].direction,
        &coords[coord_elem_offset],
        &shapeDerivs_[grad_offset],
        areaVector );

      // apply quadrature weight and orientation (combined as weight)
      for (int j = 0; j < nDim_; ++j) {
        areav[offset+j]  = ipInfo_[ip].weight * areaVector[j];
      }
    }
  }

  *error = 0;
}
//--------------------------------------------------------------------------
void
HigherOrderHexSCS::area_vector(
  const Jacobian::Direction direction,
  const double *elemNodalCoords,
  double *shapeDeriv,
  std::array<double,3>& areaVector) const
{

  int s1Component; int s2Component;
  switch (direction) {
    case Jacobian::S_DIRECTION:
      s1Component = static_cast<int>(Jacobian::T_DIRECTION);
      s2Component = static_cast<int>(Jacobian::U_DIRECTION);
      break;
    case Jacobian::T_DIRECTION:
      s1Component = static_cast<int>(Jacobian::S_DIRECTION);
      s2Component = static_cast<int>(Jacobian::U_DIRECTION);
      break;
    case Jacobian::U_DIRECTION:
      s1Component = static_cast<int>(Jacobian::T_DIRECTION);
      s2Component = static_cast<int>(Jacobian::S_DIRECTION);
      break;
    default:
      throw std::runtime_error("Not a valid direction for this element!");
  }

  // return the normal area vector given shape derivatives dnds OR dndt
  double dx_ds1 = 0.0; double dy_ds1 = 0.0; double dz_ds1 = 0.0;
  double dx_ds2 = 0.0; double dy_ds2 = 0.0; double dz_ds2 = 0.0;

  for (int node = 0; node < nodesPerElement_; ++node) {
    const int vector_offset = nDim_ * node;
    const double xCoord = elemNodalCoords[vector_offset+0];
    const double yCoord = elemNodalCoords[vector_offset+1];
    const double zCoord = elemNodalCoords[vector_offset+2];

    const double dn_ds1 = shapeDeriv[vector_offset+s1Component];
    const double dn_ds2 = shapeDeriv[vector_offset+s2Component];

    dx_ds1 += dn_ds1 * xCoord;
    dx_ds2 += dn_ds2 * xCoord;

    dy_ds1 += dn_ds1 * yCoord;
    dy_ds2 += dn_ds2 * yCoord;

    dz_ds1 += dn_ds1 * zCoord;
    dz_ds2 += dn_ds2 * zCoord;
  }

  //cross product
  areaVector[0] = dy_ds1*dz_ds2 - dz_ds1*dy_ds2;
  areaVector[1] = dz_ds1*dx_ds2 - dx_ds1*dz_ds2;
  areaVector[2] = dx_ds1*dy_ds2 - dy_ds1*dx_ds2;
}
//--------------------------------------------------------------------------
void HigherOrderHexSCS::grad_op(
  const int nelem,
  const double *coords,
  double *gradop,
  double *deriv,
  double *det_j,
  double *error)
{
  for (int k = 0; k < nelem; ++k) {
    const int coord_elem_offset = nDim_ * nodesPerElement_ * k;
    const int scalar_elem_offset = numIntPoints_ * k;
    const int grad_elem_offset = numIntPoints_ * nDim_ * nodesPerElement_ * k;

    for (int ip = 0; ip < numIntPoints_; ++ip) {
      const int grad_offset = nDim_ * nodesPerElement_ * ip;
      const int offset = grad_offset + grad_elem_offset;

      for (int j = 0; j < nodesPerElement_ * nDim_; ++j) {
        deriv[offset + j] = shapeDerivs_[grad_offset +j];
      }

      gradient( &coords[coord_elem_offset],
        &shapeDerivs_[grad_offset],
        &gradop[offset],
        &det_j[scalar_elem_offset+ip] );

      if (det_j[ip] <= 0.0) {
        *error = 1.0;
      }
    }
  }
}

//--------------------------------------------------------------------------
//-------- face_grad_op ----------------------------------------------------
//--------------------------------------------------------------------------
void HigherOrderHexSCS::face_grad_op(
  const int nelem,
  const int face_ordinal,
  const double *coords,
  double *gradop,
  double *det_j,
  double *error)
{
  const int face_offset =  nDim_ * ipsPerFace_ * nodesPerElement_ * face_ordinal;
  for (int k = 0; k < nelem; ++k) {
    const int coord_elem_offset = nDim_ * nodesPerElement_ * k;
    const int scalar_elem_offset = ipsPerFace_ * k;
    const int grad_elem_offset = ipsPerFace_ * nDim_ * nodesPerElement_ * k;

    for (int ip = 0; ip < ipsPerFace_; ++ip) {
      const int grad_offset = nDim_ * nodesPerElement_ * ip;
      const int offset = grad_offset + grad_elem_offset;

      gradient( &coords[coord_elem_offset],
        &expFaceShapeDerivs_[face_offset+grad_offset],
        &gradop[offset],
        &det_j[scalar_elem_offset+ip] );

      if (det_j[ip] <= 0.0) {
        *error = 1.0;
      }
    }
  }
}

//--------------------------------------------------------------------------
//-------- gradient --------------------------------------------------------
//--------------------------------------------------------------------------
void
HigherOrderHexSCS::gradient(
  const double* elemNodalCoords,
  const double* shapeDeriv,
  double* grad,
  double* det_j) const
{
  double dx_ds1 = 0.0;  double dx_ds2 = 0.0; double dx_ds3 = 0.0;
  double dy_ds1 = 0.0;  double dy_ds2 = 0.0; double dy_ds3 = 0.0;
  double dz_ds1 = 0.0;  double dz_ds2 = 0.0; double dz_ds3 = 0.0;

  //compute Jacobian
  for (int node = 0; node < nodesPerElement_; ++node) {
    const int vector_offset = nDim_ * node;

    const double xCoord = elemNodalCoords[vector_offset + 0];
    const double yCoord = elemNodalCoords[vector_offset + 1];
    const double zCoord = elemNodalCoords[vector_offset + 2];

    const double dn_ds1 = shapeDeriv[vector_offset + 0];
    const double dn_ds2 = shapeDeriv[vector_offset + 1];
    const double dn_ds3 = shapeDeriv[vector_offset + 2];

    dx_ds1 += dn_ds1 * xCoord;
    dx_ds2 += dn_ds2 * xCoord;
    dx_ds3 += dn_ds3 * xCoord;

    dy_ds1 += dn_ds1 * yCoord;
    dy_ds2 += dn_ds2 * yCoord;
    dy_ds3 += dn_ds3 * yCoord;

    dz_ds1 += dn_ds1 * zCoord;
    dz_ds2 += dn_ds2 * zCoord;
    dz_ds3 += dn_ds3 * zCoord;
  }

  *det_j = dx_ds1 * ( dy_ds2 * dz_ds3 - dz_ds2 * dy_ds3 )
               + dy_ds1 * ( dz_ds2 * dx_ds3 - dx_ds2 * dz_ds3 )
               + dz_ds1 * ( dx_ds2 * dy_ds3 - dy_ds2 * dx_ds3 );

  const double inv_det_j = (*det_j > 0.0) ? 1.0 / (*det_j) : 0.0;

  const double ds1_dx = inv_det_j*(dy_ds2 * dz_ds3 - dz_ds2 * dy_ds3);
  const double ds2_dx = inv_det_j*(dz_ds1 * dy_ds3 - dy_ds1 * dz_ds3);
  const double ds3_dx = inv_det_j*(dy_ds1 * dz_ds2 - dz_ds1 * dy_ds2);

  const double ds1_dy = inv_det_j*(dz_ds2 * dx_ds3 - dx_ds2 * dz_ds3);
  const double ds2_dy = inv_det_j*(dx_ds1 * dz_ds3 - dz_ds1 * dx_ds3);
  const double ds3_dy = inv_det_j*(dz_ds1 * dx_ds2 - dx_ds1 * dz_ds2);

  const double ds1_dz = inv_det_j*(dx_ds2 * dy_ds3 - dy_ds2 * dx_ds3);
  const double ds2_dz = inv_det_j*(dy_ds1 * dx_ds3 - dx_ds1 * dy_ds3);
  const double ds3_dz = inv_det_j*(dx_ds1 * dy_ds2 - dy_ds1 * dx_ds2);

  // metrics
  for (int node = 0; node < nodesPerElement_; ++node) {
    const int vector_offset = nDim_ * node;

    const double dn_ds1 = shapeDeriv[vector_offset + 0];
    const double dn_ds2 = shapeDeriv[vector_offset + 1];
    const double dn_ds3 = shapeDeriv[vector_offset + 2];

    grad[vector_offset + 0] = dn_ds1 * ds1_dx + dn_ds2 * ds2_dx + dn_ds3 * ds3_dx;
    grad[vector_offset + 1] = dn_ds1 * ds1_dy + dn_ds2 * ds2_dy + dn_ds3 * ds3_dy;
    grad[vector_offset + 2] = dn_ds1 * ds1_dz + dn_ds2 * ds2_dz + dn_ds3 * ds3_dz;
  }
}
//--------------------------------------------------------------------------
//-------- gij -------------------------------------------------------------
//--------------------------------------------------------------------------
void HigherOrderHexSCS::gij(
  const double *coords,
  double *gupperij,
  double *glowerij,
  double *deriv)
{
  SIERRA_FORTRAN(threed_gij)
    ( &nodesPerElement_,
      &numIntPoints_,
      deriv,
      coords, gupperij, glowerij);
}
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
HigherOrderQuad3DSCS::HigherOrderQuad3DSCS(const ElementDescription& elem)
: MasterElement(),
  elem_(elem)
{
  surfaceDimension_ = 2;
  nDim_ = 3;
  nodesPerElement_ = elem_.nodes1D * elem_.nodes1D;

  // set up integration rule and relevant maps on scs
  set_interior_info();

  // compute and save shape functions and derivatives at ips
  shapeFunctions_ = elem_.basisBoundary->eval_basis_weights(intgLoc_);
  shapeDerivs_ = elem_.basisBoundary->eval_deriv_weights(intgLoc_);
}
//--------------------------------------------------------------------------
void
HigherOrderQuad3DSCS::set_interior_info()
{
  nodesPerElement_ = elem_.nodes1D * elem_.nodes1D;

  //1D integration rule per sub-control volume
  numIntPoints_ = (elem_.nodes1D * elem_.nodes1D) * ( elem_.numQuad * elem_.numQuad );

  // define ip node mappings
  ipNodeMap_.resize(numIntPoints_);
  intgLoc_.resize(numIntPoints_*surfaceDimension_);
  intgLocShift_.resize(numIntPoints_*surfaceDimension_);
  ipWeight_.resize(numIntPoints_);

  // tensor product nodes x tensor product quadrature
  int vector_index_2D = 0; int scalar_index = 0;
  for (unsigned l = 0; l < elem_.nodes1D; ++l) {
    for (unsigned k = 0; k < elem_.nodes1D; ++k) {
      for (unsigned j = 0; j < elem_.numQuad; ++j) {
        for (unsigned i = 0; i < elem_.numQuad; ++i) {
          //integration point location
          intgLoc_[vector_index_2D]     = elem_.gauss_point_location(k,i);
          intgLoc_[vector_index_2D + 1] = elem_.gauss_point_location(l,j);

          //weight
          ipWeight_[scalar_index] = elem_.tensor_product_weight(k,l,i,j);

          //sub-control volume association
          ipNodeMap_[scalar_index] = elem_.tensor_product_node_map_bc(k,l);

          // increment indices
          ++scalar_index;
          vector_index_2D += surfaceDimension_;
        }
      }
    }
  }
}
//--------------------------------------------------------------------------
void
HigherOrderQuad3DSCS::shape_fcn(double* shpfc)
{
  int numShape = shapeFunctions_.size();
  for (int j = 0; j < numShape; ++j) {
    shpfc[j] = shapeFunctions_[j];
  }
}
//--------------------------------------------------------------------------
const int *
HigherOrderQuad3DSCS::ipNodeMap(
  int /*ordinal*/)
{
  // define ip->node mappings for each face (single ordinal);
  return &ipNodeMap_[0];
}
//--------------------------------------------------------------------------
void
HigherOrderQuad3DSCS::determinant(
  const int nelem,
  const double *coords,
  double *areav,
  double * /*error*/)
{
  std::array<double,3> areaVector;

  for (int k = 0; k < nelem; ++k) {
    const int coord_elem_offset = nDim_ * nodesPerElement_ * k;
    const int vector_elem_offset = nDim_ * numIntPoints_ * k;

    for (int ip = 0; ip < numIntPoints_; ++ip) {
      const int grad_offset = surfaceDimension_ * nodesPerElement_ * ip;
      const int offset = nDim_ * ip + vector_elem_offset;

      //compute area vector for this ip
      area_vector( &coords[coord_elem_offset],
        &shapeDerivs_[grad_offset],
        areaVector );

      // apply quadrature weight and orientation (combined as weight)
      for (int j = 0; j < nDim_; ++j) {
        areav[offset+j]  = ipWeight_[ip] * areaVector[j];
      }
    }
  }
}
//--------------------------------------------------------------------------
void
HigherOrderQuad3DSCS::area_vector(
  const double *elemNodalCoords,
  const double *shapeDeriv,
  std::array<double,3>& areaVector) const
{
  // return the normal area vector given shape derivatives dnds OR dndt
  double dx_ds1 = 0.0; double dy_ds1 = 0.0; double dz_ds1 = 0.0;
  double dx_ds2 = 0.0; double dy_ds2 = 0.0; double dz_ds2 = 0.0;

  for (int node = 0; node < nodesPerElement_; ++node) {
    const int vector_offset = nDim_ * node;
    const int surface_vector_offset = surfaceDimension_ * node;

    const double xCoord = elemNodalCoords[vector_offset+0];
    const double yCoord = elemNodalCoords[vector_offset+1];
    const double zCoord = elemNodalCoords[vector_offset+2];

    const double dn_ds1 = shapeDeriv[surface_vector_offset+0];
    const double dn_ds2 = shapeDeriv[surface_vector_offset+1];

    dx_ds1 += dn_ds1 * xCoord;
    dx_ds2 += dn_ds2 * xCoord;

    dy_ds1 += dn_ds1 * yCoord;
    dy_ds2 += dn_ds2 * yCoord;

    dz_ds1 += dn_ds1 * zCoord;
    dz_ds2 += dn_ds2 * zCoord;
  }

  //cross product
  areaVector[0] = dy_ds1 * dz_ds2 - dz_ds1 * dy_ds2;
  areaVector[1] = dz_ds1 * dx_ds2 - dx_ds1 * dz_ds2;
  areaVector[2] = dx_ds1 * dy_ds2 - dy_ds1 * dx_ds2;
}
//--------------------------------------------------------------------------
HigherOrderQuad2DSCV::HigherOrderQuad2DSCV(const ElementDescription& elem)
: MasterElement(),
  elem_(elem)
{
  nDim_ = elem_.dimension;
  nodesPerElement_ = elem_.nodesPerElement;

  // set up integration rule and relevant maps for scvs
  set_interior_info();

  // compute and save shape functions and derivatives at ips
  shapeFunctions_ = elem_.eval_basis_weights(intgLoc_);
  shapeDerivs_ = elem_.eval_deriv_weights(intgLoc_);
}
//--------------------------------------------------------------------------
void
HigherOrderQuad2DSCV::set_interior_info()
{
  //1D integration rule per sub-control volume
  numIntPoints_ = (elem_.nodes1D * elem_.nodes1D) * ( elem_.numQuad * elem_.numQuad );

  // define ip node mappings
  ipNodeMap_.resize(numIntPoints_);
  intgLoc_.resize(numIntPoints_*nDim_);
  intgLocShift_.resize(numIntPoints_*nDim_);
  ipWeight_.resize(numIntPoints_);

  // tensor product nodes x tensor product quadrature
  int vector_index = 0; int scalar_index = 0;
  for (unsigned l = 0; l < elem_.nodes1D; ++l) {
    for (unsigned k = 0; k < elem_.nodes1D; ++k) {
      for (unsigned j = 0; j < elem_.numQuad; ++j) {
        for (unsigned i = 0; i < elem_.numQuad; ++i) {
          intgLoc_[vector_index]     = elem_.gauss_point_location(k,i);
          intgLoc_[vector_index + 1] = elem_.gauss_point_location(l,j);
          ipWeight_[scalar_index] = elem_.tensor_product_weight(k,l,i,j);
          ipNodeMap_[scalar_index] = elem_.tensor_product_node_map(k,l);

          // increment indices
          ++scalar_index;
          vector_index += nDim_;
        }
      }
    }
  }
}
//--------------------------------------------------------------------------
void
HigherOrderQuad2DSCV::shape_fcn(double *shpfc)
{
  int numShape = shapeFunctions_.size();
  for (int j = 0; j < numShape; ++j) {
    shpfc[j] = shapeFunctions_[j];
  }
}
//--------------------------------------------------------------------------
const int *
HigherOrderQuad2DSCV::ipNodeMap(int /*ordinal*/)
{
  return &ipNodeMap_[0];
}
//--------------------------------------------------------------------------
void
HigherOrderQuad2DSCV::determinant(
  const int nelem,
  const double *coords,
  double *volume,
  double *error)
{
  *error = 0.0;
  for (int k = 0; k < nelem; ++k) {
    const int scalar_elem_offset = numIntPoints_ * k;
    const int coord_elem_offset = nDim_ * nodesPerElement_ * k;
    for (int ip = 0; ip < numIntPoints_; ++ip) {
      const int grad_offset = nDim_ * nodesPerElement_ * ip;

      //weighted jacobian determinant
      const double det_j = jacobian_determinant( &coords[coord_elem_offset],
        &shapeDerivs_[grad_offset] );

      //apply weight and store to volume
      volume[scalar_elem_offset + ip] = ipWeight_[ip] * det_j;

      //flag error
      if (det_j <= 0.0) {
        *error = 1.0;
      }
    }
  }
}
//--------------------------------------------------------------------------
double
HigherOrderQuad2DSCV::jacobian_determinant(
  const double *elemNodalCoords,
  const double *shapeDerivs) const
{
  double dx_ds1 = 0.0;  double dx_ds2 = 0.0;
  double dy_ds1 = 0.0;  double dy_ds2 = 0.0;

  for (int node = 0; node < nodesPerElement_; ++node) {
    const int vector_offset = node * nDim_;

    const double xCoord = elemNodalCoords[vector_offset + 0];
    const double yCoord = elemNodalCoords[vector_offset + 1];

    const double dn_ds1  = shapeDerivs[vector_offset + 0];
    const double dn_ds2  = shapeDerivs[vector_offset + 1];

    dx_ds1 += dn_ds1 * xCoord;
    dx_ds2 += dn_ds2 * xCoord;

    dy_ds1 += dn_ds1 * yCoord;
    dy_ds2 += dn_ds2 * yCoord;
  }

  const double det_j = dx_ds1 * dy_ds2 - dy_ds1 * dx_ds2;
  return det_j;
}
//--------------------------------------------------------------------------
HigherOrderQuad2DSCS::HigherOrderQuad2DSCS(const ElementDescription& elem)
: MasterElement(),
  elem_(elem)
{
  nDim_ = 2;
  nodesPerElement_ = elem_.nodesPerElement;
  // set up integration rule and relevant maps for scs
  set_interior_info();

  // set up integration rule and relevant maps for faces
  set_boundary_info();

  // compute and save shape functions and derivatives at ips
  shapeFunctions_ = elem_.eval_basis_weights(intgLoc_);
  shapeDerivs_ = elem_.eval_deriv_weights(intgLoc_);
  expFaceShapeDerivs_ = elem_.eval_deriv_weights(intgExpFace_);
}
//--------------------------------------------------------------------------
void
HigherOrderQuad2DSCS::set_interior_info()
{
  const int linesPerDirection = elem_.nodes1D - 1;
  const int ipsPerLine = elem_.numQuad * elem_.nodes1D;
  const int numLines = linesPerDirection * nDim_;

  numIntPoints_ = numLines * ipsPerLine;

  // define L/R mappings
  lrscv_.resize(2*numIntPoints_);

  // standard integration location
  intgLoc_.resize(numIntPoints_*nDim_);

  // shifted
  intgLocShift_.resize(numIntPoints_*nDim_);

  ipInfo_.resize(numIntPoints_);

  // specify integration point locations in a dimension-by-dimension manner

  //u-direction
  int vector_index = 0;
  int lrscv_index = 0;
  int scalar_index = 0;
  for (int m = 0; m < linesPerDirection; ++m) {
    for (unsigned l = 0; l < elem_.nodes1D; ++l) {

      int leftNode;
      int rightNode;
      int orientation;
      if (m % 2 == 0) {
        leftNode  = elem_.tensor_product_node_map(l,m);
        rightNode = elem_.tensor_product_node_map(l,m + 1);
        orientation = -1.0;
      }
      else {
        leftNode  = elem_.tensor_product_node_map(l,m + 1);
        rightNode = elem_.tensor_product_node_map(l,m);
        orientation = +1.0;
      }

      for (unsigned j = 0; j < elem_.numQuad; ++j) {

        lrscv_[lrscv_index] = leftNode;
        lrscv_[lrscv_index + 1] = rightNode;

        intgLoc_[vector_index] = elem_.gauss_point_location(l,j);
        intgLoc_[vector_index + 1] = elem_.scsLoc[m];

        //compute the quadrature weight
        ipInfo_[scalar_index].weight = orientation*elem_.tensor_product_weight(l,j);

        //direction
        ipInfo_[scalar_index].direction = Jacobian::T_DIRECTION;

        ++scalar_index;
        lrscv_index += 2;
        vector_index += nDim_;
      }
    }
  }

  //t-direction
  for (int m = 0; m < linesPerDirection; ++m) {
    for (unsigned l = 0; l < elem_.nodes1D; ++l) {

      int leftNode;
      int rightNode;
      int orientation;
      if (m % 2 == 0) {
        leftNode  = elem_.tensor_product_node_map(m,l);
        rightNode = elem_.tensor_product_node_map(m+1,l);
        orientation = +1.0;
      }
      else {
        leftNode  = elem_.tensor_product_node_map(m+1,l);
        rightNode = elem_.tensor_product_node_map(m,l);
        orientation = -1.0;
      }

      for (unsigned j = 0; j < elem_.numQuad; ++j) {

        lrscv_[lrscv_index]   = leftNode;
        lrscv_[lrscv_index+1] = rightNode;

        intgLoc_[vector_index] = elem_.scsLoc[m];
        intgLoc_[vector_index+1] = elem_.gauss_point_location(l,j);

        //compute the quadrature weight
        ipInfo_[scalar_index].weight = orientation*elem_.tensor_product_weight(l,j);

        //direction
        ipInfo_[scalar_index].direction = Jacobian::S_DIRECTION;

        ++scalar_index;
        lrscv_index += 2;
        vector_index += nDim_;
      }
    }
  }
}
//--------------------------------------------------------------------------
void
HigherOrderQuad2DSCS::set_boundary_info()
{
  const int numFaces = 2*nDim_;
  const int nodesPerFace = elem_.nodes1D;
  const int linesPerDirection = elem_.nodes1D-1;
  ipsPerFace_ = nodesPerFace*elem_.numQuad;

  const int numFaceIps = numFaces*ipsPerFace_;

  oppFace_.resize(numFaceIps);
  ipNodeMap_.resize(numFaceIps);
  oppNode_.resize(numFaceIps);
  intgExpFace_.resize(numFaceIps*nDim_);

  auto faceMap = elem_.faceNodeMap;
  auto nodeMap = elem_.nodeMap;
  auto nodeMap1D = elem_.nodeMapBC;

  auto face_node_number = [&] (int number,int faceOrdinal)
  {
    return elem_.faceNodeMap[faceOrdinal][number];
  };

  const std::array<int, 4> faceToLine = {{
      0,
      2*linesPerDirection-1,
      linesPerDirection-1,
      linesPerDirection
  }};

  const std::array<double, 4> faceLoc = {{-1.0, +1.0, +1.0, -1.0}};

  int scalar_index = 0; int vector_index = 0;
  int faceOrdinal = 0; //bottom face
  int oppFaceIndex = 0;
  for (unsigned k = 0; k < elem_.nodes1D; ++k) {
    const int nearNode = face_node_number(k,faceOrdinal);
    int oppNode = elem_.tensor_product_node_map(k, 1);

    for (unsigned j = 0; j < elem_.numQuad; ++j) {
      ipNodeMap_[scalar_index] = nearNode;
      oppNode_[scalar_index] = oppNode;
      oppFace_[scalar_index] = oppFaceIndex + faceToLine[faceOrdinal]*ipsPerFace_;

      intgExpFace_[vector_index]   = intgLoc_[oppFace_[scalar_index]*nDim_+0];
      intgExpFace_[vector_index+1] = faceLoc[faceOrdinal];

      ++scalar_index;
      vector_index += nDim_;
      ++oppFaceIndex;
    }
  }

  faceOrdinal = 1; //right face
  oppFaceIndex = 0;
  for (unsigned k = 0; k < elem_.nodes1D; ++k) {
    const int nearNode = face_node_number(k,faceOrdinal);
    int oppNode = elem_.tensor_product_node_map(elem_.nodes1D-2,k);

    for (unsigned j = 0; j < elem_.numQuad; ++j) {
      ipNodeMap_[scalar_index] = nearNode;
      oppNode_[scalar_index] = oppNode;
      oppFace_[scalar_index] = oppFaceIndex + faceToLine[faceOrdinal]*ipsPerFace_;

      intgExpFace_[vector_index]   = faceLoc[faceOrdinal];
      intgExpFace_[vector_index+1] = intgLoc_[oppFace_[scalar_index]*nDim_+1];

      ++scalar_index;
      vector_index += nDim_;
      ++oppFaceIndex;
    }
  }


  faceOrdinal = 2; //top face
  oppFaceIndex = 0;
  //NOTE: this face is reversed
  int elemNodeM1 = static_cast<int>(elem_.nodes1D-1);
  for (int k = elemNodeM1; k >= 0; --k) {
    const int nearNode = face_node_number(elem_.nodes1D-k-1,faceOrdinal);
    int oppNode = elem_.tensor_product_node_map(k, elem_.nodes1D-2);
    for (unsigned j = 0; j < elem_.numQuad; ++j) {
      ipNodeMap_[scalar_index] = nearNode;
      oppNode_[scalar_index] = oppNode;
      oppFace_[scalar_index] = (ipsPerFace_-1) - oppFaceIndex + faceToLine[faceOrdinal]*ipsPerFace_;

      intgExpFace_[vector_index] = intgLoc_[oppFace_[scalar_index]*nDim_+0];
      intgExpFace_[vector_index+1] = faceLoc[faceOrdinal];

      ++scalar_index;
      vector_index += nDim_;
      ++oppFaceIndex;
    }
  }

  faceOrdinal = 3; //left face
  oppFaceIndex = 0;
  //NOTE: this face is reversed
  for (int k = elemNodeM1; k >= 0; --k) {
    const int nearNode = face_node_number(elem_.nodes1D-k-1,faceOrdinal);
    int oppNode = elem_.tensor_product_node_map(1,k);
    for (unsigned j = 0; j < elem_.numQuad; ++j) {
      ipNodeMap_[scalar_index] = nearNode;
      oppNode_[scalar_index] = oppNode;
      oppFace_[scalar_index] = (ipsPerFace_-1) - oppFaceIndex + faceToLine[faceOrdinal]*ipsPerFace_;

      intgExpFace_[vector_index]   = faceLoc[faceOrdinal];
      intgExpFace_[vector_index+1] = intgLoc_[oppFace_[scalar_index]*nDim_+1];

      ++scalar_index;
      vector_index += nDim_;
      ++oppFaceIndex;
    }
  }
}
//--------------------------------------------------------------------------
void
HigherOrderQuad2DSCS::shape_fcn(double *shpfc)
{
  int numShape = shapeFunctions_.size();
  for (int j = 0; j < numShape; ++j) {
    shpfc[j] = shapeFunctions_[j];
  }
}
//--------------------------------------------------------------------------
const int *
HigherOrderQuad2DSCS::ipNodeMap(int ordinal)
{
  // define ip->node mappings for each face (ordinal);
  return &ipNodeMap_[ordinal*ipsPerFace_];
}
//--------------------------------------------------------------------------
void
HigherOrderQuad2DSCS::determinant(
  const int nelem,
  const double *coords,
  double *areav,
  double *error)
{
  //returns the normal vector (dyds,-dxds) for constant t curves
  //returns the normal vector (dydt,-dxdt) for constant s curves

  std::array<double, 2> areaVector;

  for (int k = 0; k < nelem; ++k) {
    const int coord_elem_offset = nDim_ * nodesPerElement_ * k;
    const int vector_elem_offset = nDim_*numIntPoints_*k;

    for (int ip = 0; ip < numIntPoints_; ++ip) {
      const int grad_offset = nDim_ * nodesPerElement_ * ip;
      const int offset = nDim_ * ip + vector_elem_offset;

      //compute area vector for this ip
      area_vector( ipInfo_[ip].direction,
                   &coords[coord_elem_offset],
                   &shapeDerivs_[grad_offset],
                   areaVector );

      // apply quadrature weight and orientation (combined as weight)
      for (int j = 0; j < nDim_; ++j) {
        areav[offset+j]  = ipInfo_[ip].weight * areaVector[j];
      }
    }
  }

  *error = 0;
}
//--------------------------------------------------------------------------
void HigherOrderQuad2DSCS::grad_op(
  const int nelem,
  const double *coords,
  double *gradop,
  double *deriv,
  double *det_j,
  double *error)
{
  *error = 0.0;

  for (int k = 0; k < nelem; ++k) {
    const int coord_elem_offset = nDim_ * nodesPerElement_ * k;
    const int scalar_elem_offset = numIntPoints_ * k;
    const int grad_elem_offset = numIntPoints_ * nDim_ * nodesPerElement_ * k;

    for (int ip = 0; ip < numIntPoints_; ++ip) {
      const int grad_offset = nDim_ * nodesPerElement_ * ip;
      const int offset = grad_offset + grad_elem_offset;

      for (int j = 0; j < nodesPerElement_ * nDim_; ++j) {
        deriv[offset + j] = shapeDerivs_[grad_offset +j];
      }

      gradient( &coords[coord_elem_offset],
                &shapeDerivs_[grad_offset],
                &gradop[offset],
                &det_j[scalar_elem_offset+ip] );

      if (det_j[ip] <= 0.0) {
        *error = 1.0;
      }
    }
  }
}
//--------------------------------------------------------------------------
void
HigherOrderQuad2DSCS::face_grad_op(
  const int nelem,
  const int face_ordinal,
  const double *coords,
  double *gradop,
  double *det_j,
  double *error)
{
  *error = 0.0;

  const int face_offset =  nDim_ * ipsPerFace_ * nodesPerElement_ * face_ordinal;
  for (int k = 0; k < nelem; ++k) {
    const int coord_elem_offset = nDim_ * nodesPerElement_ * k;
    const int scalar_elem_offset = ipsPerFace_ * k;
    const int grad_elem_offset = ipsPerFace_ * nDim_ * nodesPerElement_ * k;

    for (int ip = 0; ip < ipsPerFace_; ++ip) {
      const int grad_offset = nDim_ * nodesPerElement_ * ip;
      const int offset = grad_offset + grad_elem_offset;

      gradient( &coords[coord_elem_offset],
                &expFaceShapeDerivs_[face_offset+grad_offset],
                &gradop[offset],
                &det_j[scalar_elem_offset+ip] );

      if (det_j[ip] <= 0.0) {
        *error = 1.0;
      }
    }
  }
}
//--------------------------------------------------------------------------
void
HigherOrderQuad2DSCS::gradient(
  const double* elemNodalCoords,
  const double* shapeDeriv,
  double* grad,
  double* det_j) const
{
  double dx_ds1 = 0.0;  double dx_ds2 = 0.0;
  double dy_ds1 = 0.0;  double dy_ds2 = 0.0;

  //compute Jacobian
  for (int node = 0; node < nodesPerElement_; ++node) {
    const int vector_offset = nDim_ * node;

    const double xCoord = elemNodalCoords[vector_offset + 0];
    const double yCoord = elemNodalCoords[vector_offset + 1];
    const double dn_ds1 = shapeDeriv[vector_offset + 0];
    const double dn_ds2 = shapeDeriv[vector_offset + 1];

    dx_ds1 += dn_ds1 * xCoord;
    dx_ds2 += dn_ds2 * xCoord;

    dy_ds1 += dn_ds1 * yCoord;
    dy_ds2 += dn_ds2 * yCoord;
  }

  *det_j = dx_ds1*dy_ds2 - dy_ds1*dx_ds2;

  const double inv_det_j = (*det_j > 0.0) ? 1.0 / (*det_j) : 0.0;

  const double ds1_dx =  inv_det_j*dy_ds2;
  const double ds2_dx = -inv_det_j*dy_ds1;

  const double ds1_dy = -inv_det_j*dx_ds2;
  const double ds2_dy =  inv_det_j*dx_ds1;

  for (int node = 0; node < nodesPerElement_; ++node) {
    const int vector_offset = nDim_ * node;

    const double dn_ds1 = shapeDeriv[vector_offset + 0];
    const double dn_ds2 = shapeDeriv[vector_offset + 1];

    grad[vector_offset + 0] = dn_ds1 * ds1_dx + dn_ds2 * ds2_dx;
    grad[vector_offset + 1] = dn_ds1 * ds1_dy + dn_ds2 * ds2_dy;
  }
}
//--------------------------------------------------------------------------
const int *
HigherOrderQuad2DSCS::adjacentNodes()
{
  // define L/R mappings
  return &lrscv_[0];
}
//--------------------------------------------------------------------------
int
HigherOrderQuad2DSCS::opposingNodes(
  const int ordinal,
  const int node)
{
  return oppNode_[ordinal*ipsPerFace_+node];
}
//--------------------------------------------------------------------------
int
HigherOrderQuad2DSCS::opposingFace(
  const int ordinal,
  const int node)
{
  return oppFace_[ordinal*ipsPerFace_+node];
}
//--------------------------------------------------------------------------
void
HigherOrderQuad2DSCS::area_vector(
  const Jacobian::Direction direction,
  const double *elemNodalCoords,
  const double *shapeDeriv,
  std::array<double,2>& areaVector ) const
{
  int s1Component;
  switch (direction) {
    case Jacobian::S_DIRECTION:
      s1Component = static_cast<int>(Jacobian::T_DIRECTION);
      break;
    case Jacobian::T_DIRECTION:
      s1Component = static_cast<int>(Jacobian::S_DIRECTION);
      break;
    default:
      throw std::runtime_error("Not a valid direction for this element!");
  }

  double dxdr = 0.0;  double dydr = 0.0;
  for (int node = 0; node < nodesPerElement_; ++node) {
    const int vector_offset = nDim_ * node;
    const double xCoord = elemNodalCoords[vector_offset+0];
    const double yCoord = elemNodalCoords[vector_offset+1];

    dxdr += shapeDeriv[vector_offset+s1Component] * xCoord;
    dydr += shapeDeriv[vector_offset+s1Component] * yCoord;
  }
  areaVector[0] =  dydr;
  areaVector[1] = -dxdr;
}
//--------------------------------------------------------------------------
void HigherOrderQuad2DSCS::gij(
  const double *coords,
  double *gupperij,
  double *glowerij,
  double *deriv)
{
  SIERRA_FORTRAN(twod_gij)
    ( &nodesPerElement_,
      &numIntPoints_,
      deriv,
      coords, gupperij, glowerij);
}
//--------------------------------------------------------------------------
HigherOrderEdge2DSCS::HigherOrderEdge2DSCS(const ElementDescription& elem)
: MasterElement(),
  elem_(elem)
{
  nDim_ = 2;
  nodesPerElement_ = elem_.nodes1D;
  numIntPoints_ = elem_.numQuad * elem_.nodes1D;

  ipNodeMap_.resize(numIntPoints_);
  intgLoc_.resize(numIntPoints_);

  ipWeight_.resize(numIntPoints_);

  int scalar_index = 0;
  for (unsigned k = 0; k < elem_.nodes1D; ++k) {
    for (unsigned i = 0; i < elem_.numQuad; ++i) {
      intgLoc_[scalar_index]  = elem_.gauss_point_location(k,i);
      ipWeight_[scalar_index] = elem_.tensor_product_weight(k,i);
      ipNodeMap_[scalar_index] = elem_.tensor_product_node_map(k);
      ++scalar_index;
    }
  }

  shapeFunctions_ = elem_.basisBoundary->eval_basis_weights(intgLoc_);
  shapeDerivs_ = elem_.basisBoundary->eval_deriv_weights(intgLoc_);
}
//--------------------------------------------------------------------------
const int *
HigherOrderEdge2DSCS::ipNodeMap(int /*ordinal*/)
{
  return &ipNodeMap_[0];
}
//--------------------------------------------------------------------------
void
HigherOrderEdge2DSCS::determinant(
  const int nelem,
  const double *coords,
  double *areav,
  double *error)
{
  std::array<double,2> areaVector;

  for (int k = 0; k < nelem; ++k) {
    const int coord_elem_offset = nDim_ * nodesPerElement_ * k;

    for (int ip = 0; ip < numIntPoints_; ++ip) {
      const int offset = nDim_ * ip + coord_elem_offset;
      const int grad_offset = ip * nodesPerElement_; // times edgeDim = 1

      // calculate the area vector
      area_vector( &coords[coord_elem_offset],
                   &shapeDerivs_[grad_offset],
                   areaVector );

      // weight the area vector with the Gauss-quadrature weight for this IP
      areav[offset + 0] = ipWeight_[ip] * areaVector[0];
      areav[offset + 1] = ipWeight_[ip] * areaVector[1];
    }
  }

  // check
  *error = 0.0;
}
//--------------------------------------------------------------------------
void
HigherOrderEdge2DSCS::shape_fcn(double *shpfc)
{
  int numShape = shapeFunctions_.size();
   for (int j = 0; j < numShape; ++j) {
     shpfc[j] = shapeFunctions_[j];
   }
}
//--------------------------------------------------------------------------
void
HigherOrderEdge2DSCS::area_vector(
  const double* elemNodalCoords,
  const double* shapeDeriv,
  std::array<double,2>& areaVector) const
{
  double dxdr = 0.0;  double dydr = 0.0;
  for (int node = 0; node < nodesPerElement_; ++node) {
    const int vector_offset = nDim_ * node;
    const double xCoord = elemNodalCoords[vector_offset+0];
    const double yCoord = elemNodalCoords[vector_offset+1];

    dxdr += shapeDeriv[node] * xCoord;
    dydr += shapeDeriv[node] * yCoord;
  }
  areaVector[0] =  dydr;
  areaVector[1] = -dxdr;
}

}  // namespace nalu
} // namespace sierra
