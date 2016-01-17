/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#ifndef SteadyThermalContact3DAuxFunction_h
#define SteadyThermalContact3DAuxFunction_h

#include <AuxFunction.h>

#include <vector>

namespace sierra{
namespace nalu{

class SteadyThermalContact3DAuxFunction : public AuxFunction
{
public:

  SteadyThermalContact3DAuxFunction();

  virtual ~SteadyThermalContact3DAuxFunction() {}
  
  virtual void do_evaluate(
    const double * coords,
    const double time,
    const unsigned spatialDimension,
    const unsigned numPoints,
    double * fieldPtr,
    const unsigned fieldSize,
    const unsigned beginPos,
    const unsigned endPos) const;
  
private:
  double a_;
  double k_;
  double pi_;

};

} // namespace nalu
} // namespace Sierra

#endif
