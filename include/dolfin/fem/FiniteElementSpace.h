// Copyright (C) 2013 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2014-02-03
// Last changed: 2014-02-03

#ifndef __FINITE_ELEMENT_SPACE_H_
#define __FINITE_ELEMENT_SPACE_H_

#include <dolfin/common/types.h>
#include <dolfin/fem/FiniteElement.h>
#include <dolfin/ufl/UFLFiniteElement.h>

#include <ufc.h>

#include <string>

namespace dolfin
{

class DofMap;
class Form;
class Mesh;
class SubFunction;

class FiniteElementSpace
{

  friend class DiscreteFunction;

public:

  ///
  FiniteElementSpace(Mesh& mesh, std::string const& signature);

  ///
  FiniteElementSpace(Mesh& mesh, std::string const& finite_element_signature,
                     std::string const& dof_map_signature);

  ///
  FiniteElementSpace(Mesh& mesh, Form& form, uint const i);

  ///
  FiniteElementSpace(Mesh& mesh, ufc::finite_element& finite_element,
                     bool const finite_element_local);

  ///
  FiniteElementSpace(FiniteElementSpace const& space, uint const& i);

  ///
  ~FiniteElementSpace();

  ///
  Mesh& mesh() const;

  ///
  FiniteElement const& element() const;

  ///
  DofMap const& dofmap() const;

  //--- UFL INTERFACE ---------------------------------------------------------

#if ENABLE_UFL

  /// Returns the family of the finite element
  /// UFL + FIAT
  ufl::Family::Type const family() const;

  /// Returns the degree of the finite element
  /// UFL + FIAT
  uint const degree() const;

  ///
  operator ufl::FiniteElement const&() const
  {
    return ufl_class_;
  }

#endif

private:

  Mesh& mesh_;
  FiniteElement const finite_element_;
  DofMap& dof_map_;

  // Scratch space
  class Scratch
  {

  public:

    // Constructor
    Scratch(FiniteElement const& finite_element);

    // Destructor
    ~Scratch();

    // Value size (number of entries in tensor value)
    uint size;

    // Reference finite element space dimension
    uint dimension;

    // Local array for mapping of dofs
    uint* dofs;

    // Local array for expansion coefficients
    real* coefficients;

    // Local array for values
    real* values;

    // Local array for coordinates
    real** coordinates;

  };

  // Scratch space
  Scratch scratch;

#if ENABLE_UFL

  ufl::FiniteElement const ufl_class_;

#endif

};

//-----------------------------------------------------------------------------
inline Mesh& FiniteElementSpace::mesh() const
{
  return mesh_;
}

#if ENABLE_UFL

//-----------------------------------------------------------------------------
inline ufl::Family::Type const FiniteElementSpace::family() const
{
  return ufl_class_.family().type();
}

//-----------------------------------------------------------------------------
inline uint const FiniteElementSpace::degree() const
{
  return ufl_class_.degree();
}

#endif

}
/* namespace icorne */
#endif /* __FINITE_ELEMENT_SPACE_H_ */
