// Copyright (C) 2007-2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_ASSEMBLER_H
#define __DOLFIN_ASSEMBLER_H

#include <dolfin/common/types.h>
#include <dolfin/mesh/MeshValues.h>

namespace dolfin
{

class Form;
class GenericTensor;
class SubDomain;

/**
 *  @namespace Assembler
 *
 *  @brief This namesapce provides assembly of linear systems, or more
 *         generally, assembly of a sparse tensor from a given variational form.
 */

namespace Assembler
{

/// Assemble tensor from given variational form
void assemble( GenericTensor & A, Form & form, bool reset_tensor );

/// Assemble tensor from given variational form over a sub domain
void assemble( GenericTensor & A, Form & form,
               SubDomain const & sub_domain, bool reset_tensor );

/// Assemble tensor from given variational form over sub domains
void assemble( GenericTensor & A, Form & form,
               MeshValues< size_t, Cell > const &  cell_domains,
               MeshValues< size_t, Facet > const & exterior_facet_domains,
               MeshValues< size_t, Facet > const & interior_facet_domains,
               bool reset_tensor );

} // end namespace Assembler

} /* namespace dolfin */

#endif /* __DOLFIN_ASSEMBLER_H */
