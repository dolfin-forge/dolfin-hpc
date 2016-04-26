//
//
// Modified by Aurelien Larcher, 2014-2016.
//
// This class was imported from UNICORN and did not possess any copyright.
// It was provided with a second argument 'p' to avoid the same mistakes as in
// the UNICORN implementation where p could have the value 2 or 4 (or even more)
// at the same time in different parts of the solver.
//

#ifndef MESH_QUALITY_FUNCTION_H
#define MESH_QUALITY_FUNCTION_H

#include <dolfin/function/Function.h>

#include <dolfin/mesh/MeshQuality.h>

namespace dolfin
{

/**
 *  DOCUMENTATION:
 *
 *  @class
 *
 *  @brief
 *
 *  Cell quality: q_K
 *
 *  The mesh quality function is defined as 1./q_K^p
 *
 *  with p = 2 by default (the bigger, the sharper).
 *
 *
 */


class MeshQualityFunction : public Function
{

public:

  // The default value for p ought to be 2
  MeshQualityFunction(Mesh& mesh, uint p);

  //
  ~MeshQualityFunction();

  //
  uint rank() const;

  //
  uint dim(uint i) const;

  //
  void evaluate(real* values, const real* x, const ufc::cell& cell) const;

private:

  uint const p_;
  MeshQuality mqual_;
};

}

#endif /* MESH_QUALITY_FUNCTION_H */
