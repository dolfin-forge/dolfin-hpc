// Copyright (C) 2007 Garth N. Wells.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_SPARSITY_PATTERN_BUILDER_H
#define __DOLFIN_SPARSITY_PATTERN_BUILDER_H

namespace dolfin
{

class DofMapSet;
class Mesh;
class GenericSparsityPattern;
class UFC;

/**
 *  @class  SparsityPatternBuilder
 *
 *  @brief   This class provides functions to compute the sparsity pattern.
 */

class SparsityPatternBuilder
{

public:

  /// Build sparsity pattern
  static void build(GenericSparsityPattern& sparsity_pattern, Mesh& mesh,
                    UFC& ufc, DofMapSet const& dof_map_set);

};

} /* namespace dolfin */

#endif /* __DOLFIN_SPARSITY_PATTERN_BUILDER_H */
