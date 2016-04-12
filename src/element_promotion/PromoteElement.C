/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/
#include <element_promotion/PromoteElement.h>

#include <element_promotion/ElementDescription.h>
#include <element_promotion/FaceOperations.h>
#include <element_promotion/PromotedPartHelper.h>
#include <NaluEnv.h>

#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <stk_util/parallel/CommSparse.hpp>
#include <stk_mesh/base/Bucket.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/FieldBase.hpp>
#include <stk_topology/topology.hpp>
#include <stk_util/environment/ReportHandler.hpp>
#include <stk_util/parallel/ParallelComm.hpp>
#include <stk_mesh/base/FEMHelpers.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <iterator>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <numeric>

namespace sierra{
namespace nalu{

//==========================================================================
// Class Definition
//==========================================================================
// PromoteElement - Promotes a mesh based on a description of the new element
// connectivities and node locations
// TODO(rcknaus): allow some parts not to be promoted
// TODO(rcknaus): Get rid of "ordinal reversing" methods
//===============================c===========================================
PromoteElement::PromoteElement(ElementDescription& elemDescription)
: elemDescription_(elemDescription),
  nodesPerElement_(elemDescription.nodesPerElement),
  dimension_(elemDescription.dimension)
{
 //do nothing
}
//--------------------------------------------------------------------------
void
PromoteElement::promote_elements(
  const stk::mesh::PartVector& baseParts,
  VectorFieldType& coordinates,
  stk::mesh::BulkData& mesh)
{
  ThrowRequireMsg(mesh.in_modifiable_state(),
    "Mesh must be in a modifiable state for element promotion");

  // holds the original elements and original nodes
  baseParts_ = baseParts;

  // holds only the new nodes.  Naming convention for new parts has to follow
  // the "promotion_suffix()" convention
  promotedParts_ = promote_part_vector(baseParts);
  ThrowAssertMsg(part_vector_is_valid(promotedParts_),
    "One or more promoted parts were not declared");

  auto basePartSelector = stk::mesh::selectUnion(baseParts);
  auto nodeRequests = create_child_node_requests(elemDescription_, mesh, basePartSelector);
  determine_child_ordinals(elemDescription_,mesh, nodeRequests);
  batch_create_child_nodes(elemDescription_, mesh, nodeRequests, baseParts);

  ElemRelationsMap elemNodeMap;
  populate_elem_node_relations(elemDescription_, mesh, basePartSelector, nodeRequests, elemNodeMap);
  create_elements(mesh, baseParts, elemNodeMap);

  auto baseAndSuperElemParts = append_super_elems_to_part_vector(baseParts);
  populate_boundary_elem_node_relations(mesh, baseAndSuperElemParts);

  if (dimension_ == 2) {
    set_new_node_coords<2>(coordinates, elemDescription_, mesh, nodeRequests, elemNodeMap);
  }
  else {
    set_new_node_coords<3>(coordinates, elemDescription_, mesh, nodeRequests, elemNodeMap);
  }
}
//--------------------------------------------------------------------------
PromoteElement::NodeRequests
PromoteElement::create_child_node_requests(
  const ElementDescription& elemDescription,
  stk::mesh::BulkData& mesh,
  const stk::mesh::Selector& selector) const
{
  // Creates a set of parentids with the number of children attached to them
  // Result is passed to the batch_create_child_nodes method where the
  // nodes are actually created
  const auto& connectivities = elemDescription.addedConnectivities;
  NodeRequests requestSet;
  const stk::mesh::BucketVector& elem_buckets = mesh.get_buckets(
    stk::topology::ELEM_RANK, selector);

  for (const auto* ib : elem_buckets) {
    const stk::mesh::Bucket& b = *ib;
    const stk::mesh::Bucket::size_type length = b.size();
    for (stk::mesh::Bucket::size_type k = 0; k < length; ++k) {
      const auto& elem = b[k];
      const stk::mesh::Entity* nodes = b.begin_nodes(k);
      for (const auto& relation : connectivities) {
        const auto& parentOrdinals = relation.second;
        size_t numParents = parentOrdinals.size();
        std::vector<stk::mesh::EntityId> parentIds(numParents);

        // convert from nodes to entity ids
        for (size_t j = 0; j < numParents; ++j) {
          parentIds[j] = mesh.identifier(nodes[parentOrdinals[j]]);
        }

        // Attempt to insert a new request.  Requests with the same
        // parentIds as another in the set will not be added to the set
        auto result = requestSet.insert(ChildNodeRequest{parentIds});

        // add a shared elem regardless of whether the request is new
        result.first->add_shared_elem(elem);

        // if it's the first time a set of parents was added, also set
        // the number of children / save off a copy of the unsorted ids
        if (result.second) {
          result.first->set_num_children(relation.first.size());
          result.first->unsortedParentIds_ = std::move(parentIds);
        }
      }
    }
  }
  return requestSet;
}
//--------------------------------------------------------------------------
void
PromoteElement::determine_child_ordinals(
  const ElementDescription& elemDescription,
  const stk::mesh::BulkData& mesh,
  NodeRequests& requests) const
{
  /*  FIXME(rcknaus): clean-up / remove the ordinal reversing logic
   *
   *  For P > 2, we have to worry about the orientation of faces/edges
   *
   *  The coordinates in the addedlocations map are kept constant
   *  and the ordinals themselves are reordered such that the coordinates
   *  refer to the correct nodes
   *
   */
  unsigned nodes1D = elemDescription.nodes1D;
  unsigned numAddedNodes1D = nodes1D-2;
  unsigned numParents1D = nodes1D-numAddedNodes1D; //2
  for (auto& request : requests) {
    unsigned numShared = request.sharedElems_.size();
    request.childOrdinalsForElem_.resize(numShared);
    request.reorderedChildOrdinalsForElem_.resize(numShared);
    for (unsigned elemNumber = 0; elemNumber < numShared; ++elemNumber) {
      const auto unsortedOrdinals =
          request.determine_child_node_ordinals(mesh, elemDescription, elemNumber);
      const auto& ordinals = request.childOrdinalsForElem_[elemNumber];
      const auto& canonicalOrdinals = elemDescription.addedConnectivities.at(ordinals);

      request.reorderedChildOrdinalsForElem_[elemNumber] = reorder_ordinals<size_t>(
        ordinals,
        unsortedOrdinals,
        canonicalOrdinals,
        numParents1D,
        numAddedNodes1D
      );
    }
  }
}
//--------------------------------------------------------------------------
void
PromoteElement::batch_create_child_nodes(
  const ElementDescription& elemDescription,
  stk::mesh::BulkData& mesh,
  NodeRequests& requests,
  const stk::mesh::PartVector& node_parts) const
{
  size_t num_nodes_requested = count_requested_nodes(requests);
  std::vector<stk::mesh::EntityId> available_node_ids(num_nodes_requested);
  mesh.generate_new_ids(stk::topology::NODE_RANK, num_nodes_requested,
    available_node_ids);

  size_t it_req = 0;
  for (auto& request : requests) {
    for (unsigned j = 0; j < request.num_children(); ++j, ++it_req) {
      request.add_proc_id_pair(
        mesh.parallel_rank(), available_node_ids[it_req], j
      );
    }
    request.determine_sharing_procs(mesh);
  }

  if (mesh.parallel_size() > 1) {
    parallel_communicate_ids(elemDescription, mesh, requests);
  }

  for (auto& request : requests) {
    auto meshPartsForChildNodes = request.mesh_parts_for_child_nodes(mesh, node_parts);
    request.set_node_entity_for_request(mesh, meshPartsForChildNodes);
  }
}
//--------------------------------------------------------------------------
void
PromoteElement::parallel_communicate_ids(
  const ElementDescription& elemDescription,
  const stk::mesh::BulkData& mesh,
  NodeRequests& requests) const
{
  stk::CommSparse comm_spec(mesh.parallel());

  // If the parent nodes were on a parallel boundary,
  // send some information that will allow us to
  // find the request on the other processor
  // and decide which global_ids the new nodes should have
  for (int phase = 0; phase < 2; ++phase) {
    for (const auto& request : requests) {
      for (auto other_proc : request.sharingProcs_) {
        if (other_proc != mesh.parallel_rank()) {
          const auto& request_parents = request.unsortedParentIds_;
          const auto numChildren = request.num_children();
          comm_spec.send_buffer(other_proc).pack(request.num_parents());
          comm_spec.send_buffer(other_proc).pack(numChildren);
          for (auto parentId : request_parents) {
            comm_spec.send_buffer(other_proc).pack(parentId);
          }

          for (unsigned j = 0; j < numChildren; ++j) {
            comm_spec.send_buffer(other_proc).pack(
              request.suggested_node_id(j));
          }
        }
      }
    }

    if (phase == 0) {
      comm_spec.allocate_buffers();
    }
    else {
      comm_spec.communicate();
    }
  }

  unsigned numAddedNodes1D = elemDescription.nodes1D-2;
  for (int i = 0; i < mesh.parallel_size(); ++i) {
    if (i != mesh.parallel_rank()) {
      while (comm_spec.recv_buffer(i).remaining() != 0) {

        size_t num_parents;
        size_t num_children;
        stk::mesh::EntityId suggested_node_id;
        comm_spec.recv_buffer(i).unpack(num_parents);
        comm_spec.recv_buffer(i).unpack(num_children);
        std::vector<stk::mesh::EntityId> parentIds(num_parents);

        for (auto& parentId : parentIds) {
          comm_spec.recv_buffer(i).unpack(parentId);
        }

        // Check that this proc has a request to create nodes on the
        // edge/face sent from another proc
        auto iter = requests.find(ChildNodeRequest{parentIds});
        bool hasParents = iter != requests.end();

        std::vector<size_t> indices;
        if (hasParents) {
          indices.resize(num_children);
          std::iota(indices.begin(), indices.end(), 0);

          // nodes are compared against a single set of parent ordinals
          // for edges, the sorted parent ordinals still form a chain and can be used,
          // making the ordinals parallel consistent by construction

          // For faces/volumes, I use the fact that the ordinals are not randomly ordered
          // so I can't just use the sorted parentIds atm and have to enforce parallel consistency
          // by sending over the reference parentIds and then match indices with the ordinals to ensure
          // consistent global node ids in parallel

          if (num_children == numAddedNodes1D*numAddedNodes1D && dimension_ == 3) {
            auto request = *iter;
            unsigned elemNumber = 0u;
            unsigned numParents1D = 2u;

            // lower rank processes lead
            if (i < mesh.parallel_rank()) {
              request.unsortedParentIds_ = std::move(parentIds);
            }

            const auto childOrdinals =
                request.determine_child_node_ordinals(mesh, elemDescription, elemNumber);

            const auto& canonicalOrdinals =
                elemDescription.addedConnectivities.at(request.childOrdinalsForElem_[elemNumber]);

            indices = reorder_ordinals<size_t>(
              indices,
              childOrdinals,
              canonicalOrdinals,
              numParents1D,
              numAddedNodes1D
            );
          }
        }

        for (unsigned j = 0; j < num_children; ++j) {
          //always unpack to keep the correct place in buffer
          comm_spec.recv_buffer(i).unpack(suggested_node_id);

          // Add a proc_id pair between coincident shared nodes
          if (hasParents) {
            iter->add_proc_id_pair(i, suggested_node_id, indices[j]);
          }
        }
      }
    }
  }
}
//--------------------------------------------------------------------------
void
PromoteElement::populate_elem_node_relations(
  const ElementDescription& elemDescription,
  stk::mesh::BulkData& mesh,
  const stk::mesh::Selector& selector,
  const NodeRequests& requests,
  ElemRelationsMap& elemNodeMap)
{

  const stk::mesh::BucketVector& elem_buckets = mesh.get_buckets(
    stk::topology::ELEM_RANK, selector);

  elemNodeMap.reserve(count_entities(elem_buckets));
  nodeElemMap_.reserve(count_requested_nodes(requests));

  populate_original_elem_node_relations(mesh, selector, requests, elemNodeMap);
  populate_new_elem_node_relations(requests,elemNodeMap);
}
//--------------------------------------------------------------------------
void
PromoteElement::populate_original_elem_node_relations(
  const stk::mesh::BulkData& mesh,
  const stk::mesh::Selector& selector,
  const NodeRequests& requests,
  ElemRelationsMap& elemNodeMap)
{
  const stk::mesh::BucketVector& elem_buckets = mesh.get_buckets(
    stk::topology::ELEM_RANK, selector);

  // initialize base downward relationships
  for (const auto* ib : elem_buckets) {
    const stk::mesh::Bucket& b = *ib;
    const stk::mesh::Bucket::size_type length = b.size();
    for (stk::mesh::Bucket::size_type k = 0; k < length; ++k) {
      const stk::mesh::Entity elem = b[k];
      const stk::mesh::Entity* nodes = b.begin_nodes(k);
      elemNodeMap.insert({elem, std::vector<stk::mesh::Entity>(nodesPerElement_) });
      for (size_t j = 0; j < b.num_nodes(k); ++j) {
        elemNodeMap[elem][j] = nodes[j];
      }
    }
  }

  const stk::mesh::BucketVector& node_buckets = mesh.get_buckets(
    stk::topology::NODE_RANK, selector);

  // initialize base upward relationships
  for (const auto* ib : node_buckets) {
    const stk::mesh::Bucket& b = *ib;
    const stk::mesh::Bucket::size_type length = b.size();
    for (stk::mesh::Bucket::size_type k = 0; k < length; ++k) {
      const stk::mesh::Entity node = b[k];
      const stk::mesh::Entity* elem_rels = b.begin_elements(k);
      const size_t num_elems = b.num_elements(k);
      nodeElemMap_.insert({ node, std::vector<stk::mesh::Entity>(num_elems) });
      for (size_t j = 0; j < b.num_elements(k); ++j) {
        nodeElemMap_.at(node)[j] = elem_rels[j];
      }
    }
  }
}
//--------------------------------------------------------------------------
void
PromoteElement::populate_new_elem_node_relations(
  const NodeRequests& requests,
  ElemRelationsMap& elemNodeMap)
{
  for (const auto& request : requests) {
    unsigned numShared = request.sharedElems_.size();
    for (unsigned elemNumber = 0; elemNumber < numShared; ++elemNumber) {
      auto sharedElem = request.sharedElems_[elemNumber];
      auto& ordinals = request.reorderedChildOrdinalsForElem_[elemNumber];

      // Place the newly created nodes in the connectivity map depending on
      // the assigned ordinal for each element shared by the face
      for (unsigned j = 0; j < request.num_children(); ++j) {
        elemNodeMap.at(sharedElem)[ordinals[j]] = request.children_[j];
      }
    }

    // Save the upward connectivity for the new nodes as well
    for (const auto child : request.children_) {
      nodeElemMap_.insert({ child, request.sharedElems_ });
    }
  }
}
//--------------------------------------------------------------------------
void
PromoteElement::create_elements(
  stk::mesh::BulkData& mesh,
  const stk::mesh::PartVector& baseParts,
  ElemRelationsMap& elemNodeMap) const
{
  ThrowAssert(check_elem_node_relations(mesh, elemNodeMap));

  auto baseElemParts = base_elem_parts(baseParts);
  // Generate all new ids up front
  const auto numNewElem = count_entities(mesh.get_buckets(
    stk::topology::ELEM_RANK,
    stk::mesh::selectUnion(baseElemParts))
  );

  std::vector<stk::mesh::EntityId> availableElemIds(numNewElem);
  mesh.generate_new_ids(stk::topology::ELEM_RANK, numNewElem, availableElemIds);

  // declare super element copies for each base element
  for (const auto* ibasePart : baseElemParts) {
    const stk::mesh::Part& baseElemPart = *ibasePart;
    ThrowAssertMsg(baseElemPart.topology().rank() == stk::topology::ELEM_RANK,
      "Tried to create elements on a super-element part from a part that was not of element rank");
    const auto& elem_buckets = mesh.get_buckets(stk::topology::ELEM_RANK, baseElemPart);

    ThrowAssertMsg(super_elem_part(baseElemPart) != nullptr, "Super element part not declared");
    stk::mesh::Part& superElemPart = *super_elem_part(baseElemPart);

    std::vector<stk::mesh::EntityId> connectedNodeIds(nodesPerElement_);
    size_t elemIdIndex = 0;
    for (const auto* ib : elem_buckets) {
      const stk::mesh::Bucket& b = *ib;
      for (size_t k = 0; k < b.size(); ++k) {
        const std::vector<stk::mesh::Entity>& nodes = elemNodeMap.at(b[k]);
        std::transform(nodes.begin(), nodes.end(), connectedNodeIds.begin(),
          [&mesh] (stk::mesh::Entity e) {
            return mesh.identifier(e);
          }
        );

        stk::mesh::declare_element(
          mesh,
          superElemPart,
          availableElemIds[elemIdIndex],
          connectedNodeIds
        );
        ++elemIdIndex;
      }
    }
  }
}
//--------------------------------------------------------------------------
bool
PromoteElement::check_elem_node_relations(
  const stk::mesh::BulkData& mesh,
  ElemRelationsMap& elemNodeMap) const
{
  if (elemNodeMap.empty()) {
    return false;
  }

  if (nodeElemMap_.empty()) {
    return false;
  }

  for (const auto& elemNodePair : elemNodeMap) {
    for (const auto& node : elemNodePair.second) {
      if (!(mesh.is_valid(node))) {
        return false;
      }
    }
  }

  for (const auto& nodeElemPair : nodeElemMap_) {
    if (nodeElemPair.second.size() < 1) {
      return false;
    }
    for (const auto& elem : nodeElemPair.second) {
      if (!(mesh.is_valid(elem))) {
        return false;
      }
    }
  }
  return true;
}
//--------------------------------------------------------------------------
template<unsigned embedding_dimension> void
PromoteElement::set_new_node_coords(
  VectorFieldType& coordinates,
  const ElementDescription& elemDescription,
  const stk::mesh::BulkData& mesh,
  const NodeRequests& requests,
  const ElemRelationsMap& elemNodeMap) const
{
  //  hex/quad specific method for interpolating coordinates
  static_assert(embedding_dimension == 2 || embedding_dimension == 3,"");

  for (auto& request : requests) {
    unsigned elemIndex = 0u; // based on the first face entered in map

    auto numParents = request.parentIds_.size();
    const auto& ordinals = request.childOrdinalsForElem_[elemIndex];
    const auto* node_rels = elemNodeMap.at(request.sharedElems_[elemIndex]).data();
    const auto& childLocations = elemDescription.locationsForNewNodes.at(ordinals);
    const auto& unsortedParentIds = request.unsortedParentIds_;

    switch (numParents)
    {
      case 2:
      {
        set_coords_for_child<embedding_dimension, 1>(
          mesh, coordinates, node_rels,
          ordinals, unsortedParentIds,
          childLocations);
        break;
      }
      case 4:
      {
        set_coords_for_child<embedding_dimension, 2>(
          mesh, coordinates, node_rels,
          ordinals, unsortedParentIds,
          childLocations);
        break;
      }
      case 8:
      {
        set_coords_for_child<3, 3>(
          mesh, coordinates, node_rels,
          ordinals, unsortedParentIds,
          childLocations);
        break;
      }
      default:
      {
        throw std::runtime_error("invalid parent number");
      }
    }
  }
}
//--------------------------------------------------------------------------
template<unsigned embedding_dimension, unsigned dimension> void
PromoteElement::set_coords_for_child(
  const stk::mesh::BulkData& mesh,
  VectorFieldType& coordinates,
  const stk::mesh::Entity* node_rels,
  const std::vector<size_t>& childOrdinals,
  const std::vector<stk::mesh::EntityId>& parentNodeIds,
  const std::vector<std::vector<double>>& isoParCoords) const
{
  // Gathers the information needed for interpolation, then calls the interpolation method

  constexpr unsigned numParents = ipow(2,dimension);
  ThrowAssert(parentNodeIds.size() == numParents);

  std::array<stk::mesh::Entity,numParents> parentNodes;
  for (unsigned j = 0; j < numParents; ++j) {
    parentNodes[j] = mesh.get_entity(stk::topology::NODE_RANK, parentNodeIds[j]);
  }

  std::array<double*, numParents> parentCoordPtrs;
  for (unsigned m = 0; m < numParents; ++m) {
    parentCoordPtrs[m] = static_cast<double*>(stk::mesh::field_data(
      coordinates, parentNodes[m]
    ));
  }

  std::array<double, embedding_dimension * numParents> parentCoords;
  for (size_t m = 0; m < numParents; ++m) {
    for (size_t j = 0; j < embedding_dimension; ++j) {
      parentCoords[j + m * embedding_dimension] = parentCoordPtrs[m][j];
    }
  }

  for (unsigned j = 0; j < childOrdinals.size(); ++j) {
    auto* coords =
        static_cast<double*>(
            stk::mesh::field_data(coordinates, node_rels[childOrdinals[j]]
        )
    );

    interpolate_coords<embedding_dimension, dimension>(
      isoParCoords[j],
      parentCoords,
      coords
    );
  }
}
//--------------------------------------------------------------------------
template<unsigned embedding_dimension, unsigned dimension> void
PromoteElement::interpolate_coords(
  const std::vector<double>& isoParCoord,
  const std::array<double, embedding_dimension*ipow(2,dimension)>& parentCoords,
  double* interpolatedCoords) const
{
  // Interpolates edge, face, and volume data specifically for quad/hex elements

  static_assert ( embedding_dimension == 2 || embedding_dimension == 3, "");
  static_assert ( dimension <= embedding_dimension, "");

  constexpr unsigned num_shape = ipow(2,dimension);
  std::array<double, num_shape> shape_function;

  auto shape1D = [](double x, double xi) { return 0.5*(1.0+xi*x); };
  switch (dimension) {
    case 1:
    {
      shape_function[0] = shape1D(isoParCoord[0],-1.0);
      shape_function[1] = shape1D(isoParCoord[0],+1.0);
      break;
    }
    case 2:
    {
      const double s1 = isoParCoord[0];
      const double s2 = isoParCoord[1];
      shape_function[0] = shape1D(s1,-1.0)*shape1D(s2,-1.0);
      shape_function[1] = shape1D(s1,+1.0)*shape1D(s2,-1.0);
      shape_function[2] = shape1D(s1,+1.0)*shape1D(s2,+1.0);
      shape_function[3] = shape1D(s1,-1.0)*shape1D(s2,+1.0);
      break;
    }
    case 3:
    {
      const double s1 = isoParCoord[0];
      const double s2 = isoParCoord[1];
      const double s3 = isoParCoord[2];
      shape_function[0] = shape1D(s1,-1.0)*shape1D(s2,-1.0)*shape1D(s3,-1.0);
      shape_function[1] = shape1D(s1,+1.0)*shape1D(s2,-1.0)*shape1D(s3,-1.0);
      shape_function[2] = shape1D(s1,+1.0)*shape1D(s2,-1.0)*shape1D(s3,+1.0);
      shape_function[3] = shape1D(s1,-1.0)*shape1D(s2,-1.0)*shape1D(s3,+1.0);
      shape_function[4] = shape1D(s1,-1.0)*shape1D(s2,+1.0)*shape1D(s3,-1.0);
      shape_function[5] = shape1D(s1,+1.0)*shape1D(s2,+1.0)*shape1D(s3,-1.0);
      shape_function[6] = shape1D(s1,+1.0)*shape1D(s2,+1.0)*shape1D(s3,+1.0);
      shape_function[7] = shape1D(s1,-1.0)*shape1D(s2,+1.0)*shape1D(s3,+1.0);
      break;
    }
  }

  for (unsigned j = 0; j < embedding_dimension; ++j) {
    interpolatedCoords[j] = 0.0;
  }

  for (unsigned m = 0; m < num_shape; ++m) {
    for (unsigned j = 0; j < embedding_dimension; ++j) {
      interpolatedCoords[j] += shape_function[m] *
          parentCoords[j + m * embedding_dimension];
    }
  }
}
//--------------------------------------------------------------------------
template<typename T> std::vector<T>
PromoteElement::reorder_ordinals(
  const std::vector<T>& ordinals,
  const std::vector<T>& unsortedOrdinals,
  const std::vector<T>& canonicalOrdinals,
  unsigned numParents1D,
  unsigned numAddedNodes1D) const
{
  // unnecessary if P < 3
  if (numAddedNodes1D < 2) {
    return ordinals;
  }

  // Changes the ordinals so that the coordinate interpolation is consistent
  // Idea is to list orientations of edges/faces and then apply a transformation
  // to the ordinals that is consistent with the orientation of the edge / face
  // e.g., a "reversed edge" has its ordinals reversed.

  std::vector<T> reorderedOrdinals;
  if (parents_are_reversed<T>(unsortedOrdinals, canonicalOrdinals)) {
    reorderedOrdinals = ordinals;
    std::reverse(reorderedOrdinals.begin(), reorderedOrdinals.end());
    if(unsortedOrdinals.size() == numParents1D*numParents1D) {
      reorderedOrdinals = flip_x<T>(reorderedOrdinals, numAddedNodes1D);
    }
  }
  else if (parents_are_flipped_x<T>(unsortedOrdinals, canonicalOrdinals, numParents1D)) {
    reorderedOrdinals = flip_x<T>(ordinals,numAddedNodes1D);
  }
  else if (parents_are_flipped_y<T>(unsortedOrdinals, canonicalOrdinals, numParents1D)) {
    reorderedOrdinals = flip_y<T>(ordinals,numAddedNodes1D);
  }
  else if (should_transpose<T>(unsortedOrdinals, canonicalOrdinals)) {
    reorderedOrdinals = transpose_ordinals<T>(ordinals,numAddedNodes1D);
  }
  else if (should_invert<T>(unsortedOrdinals, canonicalOrdinals)) {
    reorderedOrdinals = invert_ordinals_yx<T>(ordinals,numAddedNodes1D);
  }
  else {
    // If all of the other checks fail, then the parent ordinals should be in
    // the canonical order.  If not, then some possible orientation was missed
    ThrowRequireMsg(unsortedOrdinals == canonicalOrdinals,
      "Element promotion: unexpected permutation of parent ordinals");

    reorderedOrdinals = ordinals;
  }
  return reorderedOrdinals;
}
//--------------------------------------------------------------------------
size_t
PromoteElement::count_requested_nodes(const NodeRequests& requests) const
{
  size_t numNodes = 0;
  for (const auto& request : requests) {
    numNodes += request.num_children();
  }
  return numNodes;
}
//--------------------------------------------------------------------------
void
PromoteElement::populate_boundary_elem_node_relations(
  const stk::mesh::BulkData& mesh,
  const stk::mesh::PartVector& mesh_parts)
{
  // Generates connectivity at exposed faces for the super elements (which lack that information)

  const auto superElemParts = only_super_elem_parts(mesh_parts);
  ThrowRequireMsg(part_vector_is_valid(superElemParts), "No super element part in part vector!");

  populate_exposed_face_to_super_elem_map(elemDescription_, mesh, mesh_parts, superElemParts);

  auto rank = (dimension_ == 2) ? stk::topology::EDGE_RANK : stk::topology::FACE_RANK;

  const stk::mesh::BucketVector& boundary_buckets = mesh.get_buckets(
    rank, stk::mesh::selectUnion(mesh_parts)
  );

  auto nodes1D = elemDescription_.nodes1D;
  std::vector<stk::mesh::Entity> faceNodes(std::pow(nodes1D,dimension_-1));
  for (const auto* ib : boundary_buckets) {
    const stk::mesh::Bucket& b = *ib;
    const stk::mesh::Bucket::size_type length = b.size();

    for (stk::mesh::Bucket::size_type k = 0; k < length; ++k) {
      const auto face = b[k];
      stk::mesh::Entity superElem = exposedFaceToSuperElemMap_.at(face)[0];

      const auto* face_elem_ords = mesh.begin_element_ordinals(face);
      const unsigned face_ordinal = face_elem_ords[0];
      const auto& elem_node_rels = mesh.begin_nodes(superElem);
      const auto& nodeOrdinalsForFace = elemDescription_.faceNodeMap[face_ordinal];

      if (dimension_ == 2) {
        for (unsigned j = 0; j < nodes1D; ++j) {
          int ordinal = elemDescription_.tensor_product_node_map_bc(j);
          faceNodes[ordinal] = elem_node_rels[nodeOrdinalsForFace[j]];
        }
      }
      else {
        for (unsigned j = 0; j < nodes1D; ++j) {
          for (unsigned i = 0; i < nodes1D; ++i) {
            int ordinal = elemDescription_.tensor_product_node_map_bc(i,j);
            faceNodes[ordinal] = elem_node_rels[nodeOrdinalsForFace[i+nodes1D*j]];
          }
        }
      }
      elemNodeMapBC_.insert({face, faceNodes});

      for (auto faceNode : faceNodes) {
        nodeElemMapBC_.insert({faceNode, {face}});
      }
    }
  }
}
//--------------------------------------------------------------------------
PromoteElement::NodesElemMap
PromoteElement::make_base_nodes_to_elem_map_at_boundary(
  const ElementDescription& elemDesc,
  const stk::mesh::BulkData& mesh,
  const stk::mesh::PartVector& mesh_parts) const
{
  /* For elements connected to a face, the method
   * generates a map between a (sorted) vector of the element's
   * node ids to the element itself
   */

  auto rank = (dimension_ == 2) ?
      stk::topology::EDGE_RANK : stk::topology::FACE_RANK;
  const auto& baseElemSideBuckets = mesh.get_buckets(
    rank,
    stk::mesh::selectUnion(mesh_parts)
  );

  const auto baseNumNodes = elemDesc.nodesInBaseElement;
  NodesElemMap nodesToElemMap;
  std::vector<stk::mesh::EntityId> parents(baseNumNodes);
  for (const auto* ib : baseElemSideBuckets) {
    const stk::mesh::Bucket& b = *ib;
    for (size_t k = 0; k < b.size(); ++k) {
      const auto face = b[k];
      ThrowAssert(mesh.num_elements(face) == 1);
      const stk::mesh::Entity parent_elem = mesh.begin_elements(face)[0];

      const auto* node_rels = mesh.begin_nodes(parent_elem);
      ThrowAssert(mesh.num_nodes(parent_elem) == parents.size());

      for (unsigned j = 0; j < baseNumNodes; ++j) {
        parents[j] = mesh.identifier(node_rels[j]);
      }
      std::sort(parents.begin(), parents.end());
      nodesToElemMap.insert({parents,parent_elem});
    }
  }
  return nodesToElemMap;
}
//--------------------------------------------------------------------------
void
PromoteElement::populate_exposed_face_to_super_elem_map(
  const ElementDescription& elemDesc,
  const stk::mesh::BulkData& mesh,
  const stk::mesh::PartVector& mesh_parts,
  const stk::mesh::PartVector& superElemParts)
{
  /* Generates a map between each exposed face and the super-element
   * notionally attached to that exposed face.
   */

  const auto& superElemBuckets = mesh.get_buckets(
    stk::topology::ELEM_RANK,
    stk::mesh::selectUnion(superElemParts)
  );
  auto nodesToElemMap =
      make_base_nodes_to_elem_map_at_boundary(elemDescription_, mesh, mesh_parts);

  std::unordered_map<stk::mesh::Entity, stk::mesh::Entity> elemToSuperElemMap;

  const auto baseNumNodes = elemDesc.nodesInBaseElement;
  std::vector<stk::mesh::EntityId> parents(baseNumNodes);
  for (const auto* ib : superElemBuckets) {
    const stk::mesh::Bucket& b = *ib;
    parents.resize(elemDesc.nodesInBaseElement);
    for (size_t k = 0; k < b.size(); ++k) {
      const auto* node_rels = b.begin_nodes(k);
      ThrowAssert(b.num_nodes(k) > baseNumNodes);

      // Requires the convention that the base nodes are stored
      // first in the elem node relations still holds
      for (unsigned j = 0; j < baseNumNodes ; ++j) {
        parents[j] = mesh.identifier(node_rels[j]);
      }
      std::sort(parents.begin(), parents.end());

      auto it = nodesToElemMap.find(parents);
      if (it != nodesToElemMap.end()) {
        const stk::mesh::Entity face = it->second;
        const stk::mesh::Entity superElem = b[k];
        elemToSuperElemMap.insert({face,superElem});
      }
    }
  }
  nodesToElemMap.clear();

  auto rank = (dimension_ == 2) ? stk::topology::EDGE_RANK : stk::topology::FACE_RANK;

  const stk::mesh::BucketVector& boundary_buckets = mesh.get_buckets(
    rank, stk::mesh::selectUnion(mesh_parts)
  );

  auto nodes1D = elemDescription_.nodes1D;
  std::vector<stk::mesh::Entity> faceNodes(std::pow(nodes1D,dimension_-1));
  for (const auto* ib : boundary_buckets) {
    const stk::mesh::Bucket& b = *ib;
    const stk::mesh::Bucket::size_type length = b.size();

    for (stk::mesh::Bucket::size_type k = 0; k < length; ++k) {
      const auto face = b[k];
      ThrowAssert(mesh.num_elements(face) == 1);
      const stk::mesh::Entity baseElem = mesh.begin_elements(face)[0];
      const stk::mesh::Entity superElem = elemToSuperElemMap.at(baseElem);
      exposedFaceToSuperElemMap_.insert({face,{superElem}});
    }
  }
}
//==========================================================================
// Class Definition
//==========================================================================
// ChildNodeRequest - Provides some utilities to help promote elements
//==========================================================================
PromoteElement::ChildNodeRequest::ChildNodeRequest(
  std::vector<stk::mesh::EntityId>  in_parentIds)
: parentIds_(std::move(in_parentIds))
{
  // The equality / hash operation for the unordered set
  // are based on the sorted parentIds.
  std::sort(parentIds_.begin(), parentIds_.end());
}
//--------------------------------------------------------------------------
void
PromoteElement::ChildNodeRequest::set_node_entity_for_request(
  stk::mesh::BulkData& mesh,
  const stk::mesh::PartVector& meshPartsForChildNodes) const
{
  // Creates the actual stk nodes on the parts indicated by
  // "meshPartsForChildNodes" and also sets their node sharing relations

  for (unsigned j = 0; j < children_.size(); ++j) {
    auto& idProcPairs = procGIdPairsFromAllProcs_[j];
    std::sort(idProcPairs.begin(), idProcPairs.end());

    children_[j] = mesh.declare_entity(
      stk::topology::NODE_RANK, get_id_for_child(j), meshPartsForChildNodes
    );

    for (auto& idProcPair : idProcPairs) {
      if (idProcPair.first != mesh.parallel_rank()) {
        mesh.add_node_sharing(children_[j], idProcPair.first);
      }
    }
  }
}
//--------------------------------------------------------------------------
stk::mesh::PartVector
PromoteElement::ChildNodeRequest::mesh_parts_for_child_nodes(
  const stk::mesh::BulkData& mesh,
  stk::mesh::PartVector parts) const
{
  // Determines the parts to place the child node on.
  // Child nodes are placed on whichever parts from the "parts" input are shared between the
  // parent nodes.

  //FIXME(rcknaus): nodes can get assigned to the wrong part in the event of a "keyhole" element
  // block
  ThrowAssert(!parts.empty());
  std::sort(parts.begin(), parts.end());

  // Iteratively remove parts that are not shared by the parent nodes
  for (unsigned i = 0; i < parentIds_.size(); ++i) {
    auto parentNode = mesh.get_entity(stk::topology::NODE_RANK, parentIds_[i]);
    stk::mesh::PartVector otherParts = mesh.bucket(parentNode).supersets();
    std::sort(otherParts.begin(), otherParts.end());

    stk::mesh::PartVector temp;
    std::set_intersection(
      parts.begin(), parts.end(),
      otherParts.begin(), otherParts.end(),
      std::back_inserter(temp)
    );
    parts = std::move(temp);
  }

  // get the promoted parts associated with the shared base parts, e.g.
  // if part_1 and part 2 is shared between all parent nodes, return
  // part_1_promoted and part_2_promoted
  transform_to_promoted_part_vector(parts);
  return parts;
}
//--------------------------------------------------------------------------
void
PromoteElement::ChildNodeRequest::determine_sharing_procs(
  const stk::mesh::BulkData& mesh) const
{
  // Sets the sharing procs for the request to be the sharing procs that
  // all parents have in common
  ThrowAssert(!parentIds_.empty());

  mesh.comm_shared_procs(
    { stk::topology::NODE_RANK, parentIds_[0] }, sharingProcs_
  );
  ThrowAssert(std::is_sorted(sharingProcs_.begin(), sharingProcs_.end()));

  std::vector<int> parentSharingProcs;
  for (unsigned i = 1; i < parentIds_.size(); ++i) {
    mesh.comm_shared_procs(
      { stk::topology::NODE_RANK, parentIds_[i] }, parentSharingProcs
    );
    ThrowAssert(std::is_sorted(parentSharingProcs.begin(), parentSharingProcs.end()));

    std::vector<int> temp;
    std::set_intersection(
      sharingProcs_.begin(), sharingProcs_.end(),
      parentSharingProcs.begin(), parentSharingProcs.end(),
      std::back_inserter(temp)
    );
    sharingProcs_ = std::move(temp);
  }
}
//--------------------------------------------------------------------------
void
PromoteElement::ChildNodeRequest::add_proc_id_pair(
  int proc_id,
  stk::mesh::EntityId id,
  int childNumber) const
{
   procGIdPairsFromAllProcs_[childNumber].emplace_back(proc_id, id);
}
//--------------------------------------------------------------------------
void
PromoteElement::ChildNodeRequest::add_shared_elem(
  const stk::mesh::Entity& elem) const
{
  sharedElems_.push_back(elem);
}
//--------------------------------------------------------------------------
stk::mesh::EntityId
PromoteElement::ChildNodeRequest::get_id_for_child(int childNumber) const
{
  ThrowAssert(std::is_sorted(
    procGIdPairsFromAllProcs_[childNumber].begin(),
    procGIdPairsFromAllProcs_[childNumber].end()
  ));
  return procGIdPairsFromAllProcs_[childNumber][0].second;
}
//--------------------------------------------------------------------------
stk::mesh::EntityId
PromoteElement::ChildNodeRequest::suggested_node_id(int childNumber) const
{
  return procGIdPairsFromAllProcs_[childNumber][0].second;
}
//--------------------------------------------------------------------------
std::vector<size_t>
PromoteElement::ChildNodeRequest::determine_child_node_ordinals(
  const stk::mesh::BulkData& mesh,
  const ElementDescription& elemDesc,
  unsigned elemNumber) const
{
  const auto& elem = sharedElems_[elemNumber];
  stk::mesh::Entity const* node_rels = mesh.begin_nodes(elem);
  const size_t numNodes = mesh.num_nodes(elem);
  unsigned numParents = parentIds_.size();
  std::vector<size_t> unsortedParentOrdinals(numParents);

  // nodes are compared against a single set of parent ordinals
  // for edges, the sorted parent ordinals still form a chain and can be used
  // making the ordinals parallel consistent by construction

  // For faces/volumes, I use the fact that the ordinals are not randomly ordered
  // so I can't just use the sorted parentIds atm and have to enforce parallel consistency
  // by sending over the reference parentIds
  const auto& referenceIds = (numParents > 2) ? unsortedParentIds_ : parentIds_;

  for (unsigned i = 0; i < numParents; ++i) {
    for (unsigned j = 0; j < numNodes; ++j) {
      if (mesh.identifier(node_rels[j]) == referenceIds[i]) {
        unsortedParentOrdinals[i] = j;
      }
    }
  }

  for (const auto& relation : elemDesc.addedConnectivities) {
    if (relation.second.size() == numParents) {
      bool isPermutation = std::is_permutation(
        relation.second.begin(),
        relation.second.end(),
        unsortedParentOrdinals.begin()
      );

      if (isPermutation) {
        childOrdinalsForElem_[elemNumber] = relation.first;
        break;
      }
    }
  }
  return unsortedParentOrdinals;
}

} // namespace nalu
}  // namespace sierra
