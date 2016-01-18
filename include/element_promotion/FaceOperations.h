#ifndef FaceOperations_h
#define FaceOperations_h

#include <stddef.h>
#include <vector>
#include <stk_util/environment/ReportHandler.hpp>

namespace sierra {
namespace nalu {

  template<class T>
  constexpr T ipow(const T base, unsigned const exponent)
  {
    return (exponent == 0) ? 1 : (base * ipow(base, exponent-1));
  }

  template<typename T> bool
  parents_are_reversed(
    const std::vector<T>& test,
    const std::vector<T>& gold)
  {
    const unsigned numParents = gold.size();
    ThrowAssert(test.size() == numParents);

    bool parentAreFlipped = true;
    for (unsigned j = 0; j < numParents; ++j) {
      if (gold[j] != test[numParents-1-j]) {
        return false;
      }
    }
    return parentAreFlipped;
  }

  template<typename T> std::vector<T>
  flip_x(
    const std::vector<T>& childOrdinals,
    unsigned size1D)
  {
    ThrowAssert(childOrdinals.size() == size1D*size1D);

    std::vector<T> reorderedOrdinals(childOrdinals.size());
    for (unsigned j = 0; j < size1D; ++j) {
      for (unsigned i = 0; i < size1D; ++i) {
        int ir = size1D-i-1;
        reorderedOrdinals[i+size1D*j] = childOrdinals[ir+size1D*j];
      }
    }
    return reorderedOrdinals;
  }

  template<typename T> bool
  parents_are_flipped_x(
    const std::vector<T>& test,
    const std::vector<T>& gold,
    unsigned size1D)
  {
    ThrowAssert(test.size() == gold.size());
    if (test.size() != size1D*size1D) {
      return false;
    }

    bool parentAreFlipped = true;
    for (unsigned j = 0; j < size1D; ++j) {
      for (unsigned i = 0; i < size1D; ++i) {
        int ir = size1D-i-1;
        if (gold[i+size1D*j] != test[ir+size1D*j]) {
          return false;
        }
      }
    }
    return parentAreFlipped;
  }

  template<typename T> std::vector<T>
  flip_y(
    const std::vector<T>& childOrdinals,
    unsigned size1D)
  {
    ThrowAssert(childOrdinals.size() == size1D*size1D);

    std::vector<T> reorderedOrdinals(childOrdinals.size());
    for (unsigned j = 0; j < size1D; ++j) {
      int jr = size1D-j-1;
      for (unsigned i = 0; i < size1D; ++i) {
        reorderedOrdinals[i+size1D*j] = childOrdinals[i+size1D*jr];
      }
    }
    return reorderedOrdinals;
  }

  template<typename T> bool
  parents_are_flipped_y(
    const std::vector<T>& test,
    const std::vector<T>& gold,
    unsigned size1D)
  {
    ThrowAssert(test.size() == gold.size());
    if (test.size() != size1D*size1D) {
      return false;
    }
    bool parentAreFlipped = true;
    for (unsigned j = 0; j < size1D; ++j) {
      int jr = size1D-j-1;
      for (unsigned i = 0; i < size1D; ++i) {
        if (gold[i+size1D*j] != test[i+size1D*jr]) {
          return false;
        }
      }
    }
    return parentAreFlipped;
  }

  template<typename T> bool
  should_transpose(
    const std::vector<T>& test,
    const std::vector<T>& gold)
  {
    ThrowAssert(test.size() == gold.size());
    if (test.size() != 4) {
      return false;
    }

    //check if off-diagonal nodes ((1,0)->1 and (0,1)->3) are reversed
    // if so, transpose the ordinals
    return (test[1] == gold[3] && test[3] == gold[1]);
  }


  template<typename T> std::vector<T>
  transpose_ordinals(
    const std::vector<T>& childOrdinals,
    unsigned size1D)
  {
    ThrowAssert(childOrdinals.size() == size1D*size1D);

    std::vector<T> reorderedOrdinals(childOrdinals.size());
    for (unsigned j = 0; j < size1D; ++j) {
      for (unsigned i = 0; i < size1D; ++i) {
        reorderedOrdinals[i+size1D*j] = childOrdinals[j+size1D*i];
      }
    }
    return reorderedOrdinals;
  }

  template<typename T> bool
  should_invert(
    const std::vector<T>& test,
    const std::vector<T>& gold)
  {
    ThrowAssert(test.size() == gold.size());
    if (test.size() != 4) {
      return false;
    }

    //check if diagonal nodes ((0,0)->0 and (1,1)->2) are reversed
    // if so, transpose-reverse the ordinals
    return (test[0] == gold[2]
         && test[1] == gold[1]
         && test[2] == gold[0]
         && test[3] == gold[3] );
  }

  template<typename T> std::vector<T>
  invert_ordinals_yx(
    const std::vector<T>& childOrdinals,
    unsigned size1D)
  {
    ThrowAssert(childOrdinals.size() == size1D*size1D);

    std::vector<T> reorderedOrdinals(childOrdinals.size());
    for (unsigned j = 0; j < size1D; ++j) {
      int jr = size1D-j-1;
      for (unsigned i = 0; i < size1D; ++i) {
        int ir = size1D-i-1;
        reorderedOrdinals[i+size1D*j] = childOrdinals[jr+size1D*ir];
      }
    }
    return reorderedOrdinals;
  }

} // namespace nalu
} // namespace Sierra

#endif
