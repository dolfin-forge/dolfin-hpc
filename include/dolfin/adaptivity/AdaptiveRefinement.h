// Copyright (C) 2010 Niclas Jansson.
// Licensed under the GNU LGPL Version 2.1.

#pragma once

#include <dolfin/common/types.h>
#include <dolfin/mesh/MeshValues.h>

namespace dolfin
{

template<class T> class Array;
class Function;

namespace AdaptiveRefinement
{

typedef Array< std::pair< std::string, Function * > > FunctionMapping;

/// Refine mesh using "simple" of "rivara" strategy
void refine( Mesh & mesh, MeshValues< bool, Cell > & cell_marker );

///
void refine_and_project( Mesh & mesh,
                         FunctionMapping const & functions,
                         MeshValues< bool, Cell > &  cell_marker );

} // namespace AdaptiveRefinement

} // namespace dolfin
