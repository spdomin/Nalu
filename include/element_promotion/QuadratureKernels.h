#ifndef QuadratureKernels_h
#define QuadratureKernels_h

#include <Teuchos_BLAS.hpp>

namespace sierra {
namespace nalu {

  struct ElementDescription;

  class SGLQuadratureOps
  {
  public:
    SGLQuadratureOps(const ElementDescription& elem);

    void volume_2D(
      const double*  nodalValuesTensor,
      double* result);

    void volume_3D(const double*  nodalValue, double* result);

    void surface_2D(
      const double*  integrand,
      double* result,
      int line_offset);

    void surface_3D(
      const double* __restrict__ integrand,
      double* __restrict__ result,
      int face_offset);

  private:
    const Teuchos::BLAS<int,double> blas_;
    int nodes1D_;
    int nodesPerElement_;
    std::vector<double> work2D_;
    std::vector<double> weightTensor_;
    std::vector<double> weightMatrix_;
    double* p_weightTensor_;
    double* p_weightMatrix_;
    double* p_work2D_;

  };

} // namespace nalu
} // namespace Sierra

#endif
