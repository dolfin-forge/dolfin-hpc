// Copyright (C) 2007 Garth N. Wells.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_SPARSITY_PATTERN_BUILDER_H
#define __DOLFIN_SPARSITY_PATTERN_BUILDER_H

namespace dolfin
{

class Form;
class GenericSparsityPattern;

/**
 *  @brief   provides functions to compute the sparsity pattern.
 */

namespace SparsityPatternBuilder
{

/// Build sparsity pattern
void build( GenericSparsityPattern & sparsity_pattern,
            Form &                   form );

} // end namespace SparsityPatternBuilder

} /* namespace dolfin */

#endif /* __DOLFIN_SPARSITY_PATTERN_BUILDER_H */
