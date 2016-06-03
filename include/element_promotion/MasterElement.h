/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef MasterElement_h
#define MasterElement_h

#include <stdexcept>
#include <vector>

namespace sierra{
namespace naluUnit{

namespace Jacobian{
enum Direction
  {
    S_DIRECTION = 0,
    T_DIRECTION = 1,
    U_DIRECTION = 2
  };
  }

class MasterElement
{
public:

  MasterElement();
  virtual ~MasterElement();

  virtual void determinant(
    const int nelem,
    const double *coords,
    double *volume,
    double * error ) {
    throw std::runtime_error("determinant not implemented");}

  virtual void grad_op(
    const int nelem,
    const double *coords,
    double *gradop,
    double *deriv,
    double *det_j,
    double * error ) {
    throw std::runtime_error("grad_op not implemented");}

  virtual void shifted_grad_op(
    const int nelem,
    const double *coords,
    double *gradop,
    double *deriv,
    double *det_j,
    double * error ) {
    throw std::runtime_error("grad_op not implemented");}

  virtual void gij(
    const double *coords,
    double *gupperij,
    double *glowerij,
    double *deriv) {
    throw std::runtime_error("gij not implemented");}

  virtual void nodal_grad_op(
    const int nelem,
    double *deriv,
    double * error ) {
    throw std::runtime_error("nodal_grad_op not implemented");}

  virtual void face_grad_op(
    const int nelem,
    const int face_ordinal,
    const double *coords,
    double *gradop,
    double *det_j,
    double * error ) {
    throw std::runtime_error("face_grad_op not implemented; avoid this element type at open bcs, walls and symms");}

  virtual const int * adjacentNodes() {
    throw std::runtime_error("adjacentNodes not implementedunknown bc");
    return NULL;}

  virtual const int * ipNodeMap(int ordinal = 0) {
    throw std::runtime_error("ipNodeMap not implemented");
    return NULL;}

  virtual void shape_fcn(
    double *shpfc) {
    throw std::runtime_error("shape_fcn not implemented"); }

  virtual void shifted_shape_fcn(
    double *shpfc) {
    throw std::runtime_error("shifted_shape_fcn not implemented"); }

  virtual int opposingNodes(
    const int ordinal, const int node) {
    throw std::runtime_error("adjacentNodes not implemented"); }

  virtual int opposingFace(
    const int ordinal, const int node) {
    throw std::runtime_error("opposingFace not implemented");
    return 0; }

  virtual double isInElement(
    const double *elemNodalCoord,
    const double *pointCoord,
    double *isoParCoord) {
    throw std::runtime_error("isInElement not implemented");
    return 1.0e6; }

  virtual void interpolatePoint(
    const int &nComp,
    const double *isoParCoord,
    const double *field,
    double *result) {
    throw std::runtime_error("interpolatePoint not implemented"); }

  virtual void general_shape_fcn(
    const int numIp,
    const double *isoParCoord,
    double *shpfc) {
    throw std::runtime_error("general_shape_fcn not implement"); }

  virtual void general_face_grad_op(
    const int face_ordinal,
    const double *isoParCoord,
    const double *coords,
    double *gradop,
    double *det_j,
    double * error ) {
    throw std::runtime_error("general_face_grad_op not implemented");}

  virtual void sidePcoords_to_elemPcoords(
    const int & side_ordinal,
    const int & npoints,
    const double *side_pcoords,
    double *elem_pcoords) {
    throw std::runtime_error("sidePcoords_to_elemPcoords");}

  virtual const int * faceNodeOnExtrudedElem() {
    throw std::runtime_error("faceNodeOnExtrudedElem not implement"); }

  virtual const int * opposingNodeOnExtrudedElem() {
    throw std::runtime_error("opposingNodeOnExtrudedElem not implement"); }

  virtual const int * faceScsIpOnExtrudedElem() {
    throw std::runtime_error("faceScsIpOnExtrudedElem not implement"); }

  virtual const int * faceScsIpOnFaceEdges() {
    throw std::runtime_error("faceScsIpOnFaceEdges not implement"); }

  virtual const double * edgeAlignedArea() {
    throw std::runtime_error("edgeAlignedArea not implement"); }

  double isoparametric_mapping(const double b, const double a, const double xi) const;

  int nDim_;
  int nodesPerElement_;
  int numIntPoints_;
  double scaleToStandardIsoFac_;

  std::vector<int> lrscv_;
  std::vector<int> ipNodeMap_;
  std::vector<int> oppNode_;
  std::vector<int> oppFace_;
  std::vector<double> intgLoc_;
  std::vector<double> intgLocShift_;
  std::vector<double> intgExpFace_;
  std::vector<double> nodeLoc_;
  // extrusion-based scheme
  std::vector<int> faceNodeOnExtrudedElem_;
  std::vector<int> opposingNodeOnExtrudedElem_;
  std::vector<int> faceScsIpOnExtrudedElem_;
  std::vector<int> faceScsIpOnFaceEdges_;
  std::vector<double> edgeAlignedArea_;

};

} // namespace nalu
} // namespace Sierra

#endif
