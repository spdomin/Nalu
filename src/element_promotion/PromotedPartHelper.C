/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/
#include <element_promotion/PromotedPartHelper.h>

#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Part.hpp>
#include <stk_util/environment/ReportHandler.hpp>

#include <algorithm>
#include <vector>
#include <string>

namespace sierra{
namespace nalu{

  // A set of functions to deal with part pairs having specific ending tags

  bool part_vector_is_valid(const stk::mesh::PartVector& parts) {
    return (std::find(parts.begin(), parts.end(), nullptr) == parts.end() && !parts.empty() );
  }
  //--------------------------------------------------------------------------
  std::string promotion_suffix()
  {
    return "_promoted";
  }
  //--------------------------------------------------------------------------
  std::string super_element_suffix()
  {
    return "_se";
  }
  //--------------------------------------------------------------------------
  std::string promote_part_name(const std::string& base_name)
  {
    ThrowAssertMsg(!base_name.empty(), "Empty base name for promoted part");
    return (base_name+promotion_suffix());
  }
  //--------------------------------------------------------------------------
  std::string super_element_part_name(const std::string& base_name)
  {
    ThrowAssertMsg(!base_name.empty(), "Empty base name for super elem part");
    return (base_name+super_element_suffix());
  }
  //--------------------------------------------------------------------------
  std::string base_element_part_name(const std::string& super_name)
  {
    ThrowRequireMsg(super_name.find(super_element_suffix()) != std::string::npos,
      "Not a super-element part name!");
    return (super_name.substr(0,super_name.length()-super_element_suffix().length()));
  }
  //--------------------------------------------------------------------------
  stk::mesh::Part* super_elem_part(const stk::mesh::Part& part)
  {
    return (part.mesh_meta_data().get_part(super_element_part_name(part.name())));
  }
  //--------------------------------------------------------------------------
  stk::mesh::Part* base_elem_part_from_super_elem_part(const stk::mesh::Part& super_elem_part)
  {
    return (super_elem_part.mesh_meta_data().get_part(base_element_part_name(super_elem_part.name())));
  }
  //--------------------------------------------------------------------------
  stk::mesh::Part* super_elem_part(const stk::mesh::Part* part)
  {
    ThrowAssert(part != nullptr);
    return (part->mesh_meta_data().get_part(super_element_part_name(part->name())));
  }
  //--------------------------------------------------------------------------
  stk::mesh::Part* promoted_part(const stk::mesh::Part& part)
  {
    return (part.mesh_meta_data().get_part(promote_part_name(part.name())));
  }
  //--------------------------------------------------------------------------
  stk::mesh::Part* promoted_part(const stk::mesh::Part* part)
  {
    ThrowAssert(part != nullptr);
    return (part->mesh_meta_data().get_part(promote_part_name(part->name())));
  }
  //--------------------------------------------------------------------------
  void transform_to_promoted_part_vector(stk::mesh::PartVector& parts)
  {
    ThrowAssert(part_vector_is_valid(parts));
    std::transform(parts.begin(), parts.end(), parts.begin(), [](stk::mesh::Part* part) {
      return promoted_part(part);
    });
  }
  //--------------------------------------------------------------------------
  void transform_to_super_elem_part_vector(stk::mesh::PartVector& parts)
  {
    ThrowAssert(part_vector_is_valid(parts));
    std::transform(parts.begin(), parts.end(), parts.begin(), [](stk::mesh::Part* part) {
      return super_elem_part(part);
    });
  }
  //--------------------------------------------------------------------------
  stk::mesh::PartVector promote_part_vector(stk::mesh::PartVector parts)
  {
    transform_to_promoted_part_vector(parts);
    return parts;
  }
  //--------------------------------------------------------------------------
  stk::mesh::PartVector base_elem_parts(const stk::mesh::PartVector& parts)
  {
    ThrowAssert(part_vector_is_valid(parts));

    stk::mesh::PartVector elemParts;
    std::copy_if(parts.begin(), parts.end(), std::back_inserter(elemParts), [](stk::mesh::Part* p) {
      return (p->topology().rank() == stk::topology::ELEM_RANK
           && p->topology() < stk::topology::SUPERELEMENT_START);
    });
    return elemParts;
  }
  //--------------------------------------------------------------------------
  stk::mesh::PartVector only_super_elem_parts(const stk::mesh::PartVector& parts)
  {
    ThrowAssert(part_vector_is_valid(parts));

    stk::mesh::PartVector elemParts;
    std::copy_if(parts.begin(), parts.end(), std::back_inserter(elemParts), [](stk::mesh::Part* p) {
      return (p->topology() >= stk::topology::SUPERELEMENT_START);
    });
    return elemParts;
  }
  //--------------------------------------------------------------------------
  stk::mesh::PartVector super_elem_part_vector(const stk::mesh::PartVector& parts)
  {
    auto baseElemParts = base_elem_parts(parts);
    transform_to_super_elem_part_vector(baseElemParts);
    return baseElemParts;
  }
  //--------------------------------------------------------------------------
  stk::mesh::PartVector append_super_elems_to_part_vector(stk::mesh::PartVector parts)
  {
    auto superElemParts = super_elem_part_vector(parts);
    parts.insert(parts.end(), superElemParts.begin(), superElemParts.end());
    return parts;
  }
  //--------------------------------------------------------------------------
  size_t
  num_sub_elements(
    const stk::mesh::MetaData& metaData,
    const stk::mesh::BucketVector& buckets,
    unsigned polyOrder)
  {
    unsigned numEntities = 0;
    for (const auto* ib : buckets) {
      unsigned subElemsPerElem =
          (ib->topology().rank() == metaData.side_rank()) ?
              std::pow(polyOrder, metaData.spatial_dimension() - 1) :
              std::pow(polyOrder, metaData.spatial_dimension());

      numEntities += ib->size()*subElemsPerElem;
    }
    return (numEntities);
  }
  //--------------------------------------------------------------------------
  size_t
  count_entities(const stk::mesh::BucketVector& buckets)
  {
    unsigned numEntities = 0;
    for (const auto* ib : buckets) {
      numEntities += ib->size();
    }
    return numEntities;
  }

} // namespace nalu
}  // namespace sierra
