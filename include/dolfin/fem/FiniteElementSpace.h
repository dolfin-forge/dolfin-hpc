// Copyright (C) 2013 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2014-02-03
// Last changed: 2014-02-03

#ifndef __FINITE_ELEMENT_SPACE_H_
#define __FINITE_ELEMENT_SPACE_H_

#include <dolfin/common/types.h>
#include <dolfin/fem/FiniteElement.h>
#include <dolfin/fem/ScratchSpace.h>
#include <dolfin/mesh/Cell.h>
#include <dolfin/ufl/UFLFiniteElementBase.h>

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

#if ENABLE_UFL
  ///
  explicit FiniteElementSpace(Mesh& mesh,
                              ufl::FiniteElementBase const& finite_element);
#endif

  ///
  FiniteElementSpace(FiniteElementSpace const& space, uint const& i);

  ///
  ~FiniteElementSpace();

  ///
  Mesh& mesh() const;

  ///
  Cell& cell() const; //FIXME: Cannot const this due to Mesh implementation

  ///
  FiniteElement const& element() const;

  ///
  DofMap const& dofmap() const;

  /// Display basic information
  void disp() const;

  //---------------------------------------------------------------------------

  ///
  bool is_cellwise_defined() const;

  ///
  bool is_cellwise_constant() const;

  //--- UFL INTERFACE ---------------------------------------------------------

#if ENABLE_UFL

  /// Returns the family of the finite element
  /// UFL + FIAT
  ufl::Family::Type const family() const;

  /// Returns the degree of the finite element
  /// UFL + FIAT
  uint const degree() const;

  ///
  operator ufl::FiniteElementBase const&() const
  {
    return *ufl_;
  }

#endif

private:

  Mesh& mesh_;
  mutable Cell cell_;
  FiniteElement const finite_element_;
  DofMap& dof_map_; // The dof map is owned by the DofMapCache instance.

  // Scratch space
  ScratchSpace scratch;

#if ENABLE_UFL
  ufl::FiniteElementBase const * const ufl_;
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
  return ufl_->family().type();
}

//-----------------------------------------------------------------------------
inline uint const FiniteElementSpace::degree() const
{
  return ufl_->degree();
}

#endif

}
/* namespace icorne */
#endif /* __FINITE_ELEMENT_SPACE_H_ */
