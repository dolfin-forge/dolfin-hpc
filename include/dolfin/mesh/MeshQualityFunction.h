//
//
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

  //
  MeshQualityFunction(Mesh& mesh, uint p = 2);

  //
  ~MeshQualityFunction();

  //
  uint rank() const;

  //
  uint dim(uint i) const;

  //
  void eval(real* values, const real* x) const;

private:

  uint const p_;
  MeshQuality mqual_;
};

}

#endif /* MESH_QUALITY_FUNCTION_H */
