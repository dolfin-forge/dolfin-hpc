
#ifndef __DOLFIN_MESH_EUCLIDEAN_BASIS_H
#define __DOLFIN_MESH_EUCLIDEAN_BASIS_H

#include <dolfin/common/types.h>
#include <dolfin/common/Array.h>
#include <dolfin/mesh/Point.h>

namespace dolfin
{

class EuclideanBasis
{

public:

  ///
  static uint compute(uint gdim, Point B[], Array<real> const& N,
                      Array<real> const& W, real cosalpha_max, bool weighted);


};

} /* namespace dolfin */

#endif /* __DOLFIN_MESH_EUCLIDEAN_BASIS_H */
