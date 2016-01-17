#ifndef PromotedElementIO_h
#define PromotedElementIO_h

#include <stk_mesh/base/CoordinateSystems.hpp>
#include <stk_mesh/base/FieldBase.hpp>
#include <Ioss_Region.h>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/Types.hpp>

#include <stddef.h>
#include <iosfwd>
#include <map>
#include <memory>
#include <string>
#include <unordered_set>
#include <vector>

namespace Ioss {
class DatabaseIO;
class ElementBlock;
class NodeBlock;
class SideBlock;
}  // namespace Ioss
namespace sierra {
namespace nalu {
class PromoteElement;
struct ElementDescription;
}  // namespace nalu
}  // namespace sierra

// field types
typedef stk::mesh::Field<double>  ScalarFieldType;
typedef stk::mesh::Field<double, stk::mesh::SimpleArrayTag>  GenericFieldType;
typedef stk::mesh::Field<double, stk::mesh::Cartesian>  VectorFieldType;

namespace stk {
  namespace mesh {
    class BulkData;
    class MetaData;
    class Part;

    typedef std::vector<Part*> PartVector;
  }
}

namespace sierra {
namespace nalu {

class PromotedElementIO
{

public:
  // constructor/destructor
  PromotedElementIO(
    const PromoteElement& promoteElement,
    const stk::mesh::MetaData& metaData,
    const stk::mesh::BulkData& bulkData,
    const VectorFieldType& coordinates,
    const std::string& fileName
  );

  virtual ~PromotedElementIO() = default;

  void output_results(const std::vector<const stk::mesh::FieldBase*> fields) const;

  void write_database_data(double currentTime);

  void write_element_connectivity(
    const stk::mesh::PartVector& baseParts,
    const std::vector<stk::mesh::EntityId>& entityIds);

  void write_sideset_connectivity(
      const stk::mesh::PartVector& baseParts);

  size_t sub_element_global_id() const;
  void write_node_block_definitions(
    const stk::mesh::PartVector& baseParts,
    const stk::mesh::PartVector& promotedParts);
  void write_elem_block_definitions(const stk::mesh::PartVector& baseParts);
  void write_sideset_definitions(const stk::mesh::PartVector& baseParts);
  void write_coordinate_list(
    const stk::mesh::PartVector& baseParts,
    const stk::mesh::PartVector& promotedParts);

  void add_fields(const std::vector<stk::mesh::FieldBase*>& fields);
  bool check_topology(const stk::mesh::PartVector& baseParts) const;
  int maximum_field_length(const stk::mesh::FieldBase& field) const;

  template<typename T> void
  put_data_on_node_block(
    Ioss::NodeBlock& nodeBlock,
    const stk::mesh::FieldBase& field,
    const stk::mesh::BucketVector& buckets,
    size_t numNodes) const;

  std::string storage_name(const stk::mesh::FieldBase& field) const;

  // meta, bulk and io
  const PromoteElement& promoteElement_;
  const stk::mesh::MetaData& metaData_;
  const stk::mesh::BulkData& bulkData_;
  const VectorFieldType& coordinates_;
  const std::string& fileName_;
  const ElementDescription& elem_;
  const unsigned nDim_;

  struct FieldNameHash {
    std::size_t operator()(const stk::mesh::FieldBase* const  field) const  {
      return std::hash<std::string>()(field->name());
    }
  };

  struct FieldNameEqual {
    bool operator()(const stk::mesh::FieldBase* fieldA, const stk::mesh::FieldBase* fieldB) const {
      return (fieldA->name() == fieldB->name());
    }
  };
  std::unordered_set<const stk::mesh::FieldBase*, FieldNameHash, FieldNameEqual> fields_;

  std::map<const stk::mesh::Part*, Ioss::ElementBlock*> elementBlockPointers_;
  std::map<const stk::mesh::Part*, Ioss::SideBlock*> sideBlockPointers_;
  Ioss::NodeBlock* nodeBlock_;

  std::unique_ptr<Ioss::Region> output_;
  Ioss::DatabaseIO* databaseIO;

  // fields


};

} // namespace nalu
} // namespace Sierra

#endif
