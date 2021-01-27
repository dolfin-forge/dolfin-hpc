
#ifndef __DOLFIN_MESH_EUCLIDEAN_BASIS_H
#define __DOLFIN_MESH_EUCLIDEAN_BASIS_H

#include <dolfin/common/types.h>
#include <dolfin/mesh/Point.h>

namespace dolfin
{

class EuclideanBasis
{

public:
  ///
  static auto compute( uint                        gdim,
                       Point                       B[],
                       std::vector< real > const & N,
                       std::vector< real > const & W,
                       real                        cosalpha_max,
                       bool                        weighted ) -> uint;
};

} /* namespace dolfin */

#endif /* __DOLFIN_MESH_EUCLIDEAN_BASIS_H */
