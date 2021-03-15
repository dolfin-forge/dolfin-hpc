// Copyright (C) 2010 Niclas Jansson.
// Licensed under the GNU LGPL Version 2.1.

#pragma once

#include <dolfin/common/types.h>
#include <dolfin/mesh/MeshValues.h>

namespace dolfin
{

class Function;

namespace AdaptiveRefinement
{

using FunctionMapping = std::vector< std::pair< std::string, Function * > >;

/// Refine mesh using "simple" of "rivara" strategy
void refine( Mesh & mesh, MeshValues< bool, Cell > & cell_marker );

///
void refine_and_project( Mesh &                     mesh,
                         FunctionMapping const &    functions,
                         MeshValues< bool, Cell > & cell_marker );

} // namespace AdaptiveRefinement

} // namespace dolfin
