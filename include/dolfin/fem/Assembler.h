// Copyright (C) 2007-2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Garth N. Wells, 2007
// Modified by Niclas Jansson, 2010
//
// First added:  2007-01-17
// Last changed: 2010-03-18

#ifndef __DOLFIN_ASSEMBLER_H
#define __DOLFIN_ASSEMBLER_H

#include <ufc.h>

#include <dolfin/common/types.h>

namespace dolfin
{

class Coefficient;
class DofMapSet;
class Form;
class GenericTensor;
class Mesh;
class SubDomain;
class UFC;
template<class T> class Array;
template<class T> class MeshFunction;

/**
 *  @class Assembler
 *
 *  @brief  This class provides assembly of linear systems, or more generally,
 *          assembly of a sparse tensor from a given variational form.
 */

class Assembler
{

public:

  /// Constructor
  Assembler();

  /// Destructor
  ~Assembler();

  /// Assemble tensor from given variational form
  void assemble(GenericTensor& A, Form& form, bool reset_tensor);

  /// Assemble tensor from given variational form over a sub domain
  void assemble(GenericTensor& A, Form& form, SubDomain const& sub_domain,
                bool reset_tensor);

  /// Assemble tensor from given variational form over sub domains
  void assemble(GenericTensor& A, Form& form,
                MeshFunction<uint> const& cell_domains,
                MeshFunction<uint> const& exterior_facet_domains,
                MeshFunction<uint> const& interior_facet_domains,
                bool reset_tensor);

private:

  /// Assemble tensor from given (UFC) form, coefficients and sub domains.
  /// This is the main assembly function in DOLFIN. All other assembly functions
  /// end up calling this function.
  ///
  /// The MeshFunction arguments can be used to specify assembly over subdomains
  /// of the mesh cells, exterior facets and interior facets.
  /// Either a null pointer or an empty MeshFunction may be used to specify that
  /// the tensor should be assembled over the entire set of cells or facets.
  void assemble(GenericTensor& A, const Form& form,
                Array<Coefficient*> const& coefficients,
                DofMapSet const& dofmaps,
                MeshFunction<uint> const* cell_domains,
                MeshFunction<uint> const* exterior_facet_domains,
                MeshFunction<uint> const* interior_facet_domains,
                bool reset_tensor = true);

  // Assemble over cells
  void assembleCells(GenericTensor& A, Array<Coefficient*> const& coefficients,
                     DofMapSet const& dofmaps, UFC& data,
                     MeshFunction<uint> const* domains) const;

  // Assemble over exterior facets
  void assembleExteriorFacets(GenericTensor& A,
                              Array<Coefficient*> const& coefficients,
                              DofMapSet const& dofmaps, UFC& data,
                              MeshFunction<uint> const* domains) const;

  // Assemble over interior facets
  void assembleInteriorFacets(GenericTensor& A,
                              Array<Coefficient*> const& coefficients,
                              DofMapSet const& dofmaps, UFC& data,
                              MeshFunction<uint> const* domains) const;

  // Bogus-assemble periodic contributions
  void initializePeriodicDofs(GenericTensor& A,
                              Array<Coefficient*> const& coefficients,
                              DofMapSet const& dofmaps, UFC& data,
                              MeshFunction<uint> const* domains) const;

  // Initialize global tensor
  void initGlobalTensor(GenericTensor& A, DofMapSet const& dofmaps, UFC& ufc,
                        bool reset_tensor) const;

};

} /* namespace dolfin */

#endif /* __DOLFIN_ASSEMBLER_H */
