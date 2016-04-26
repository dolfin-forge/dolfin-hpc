// Copyright (C) 2005-2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2005-11-28
// Last changed: 2008-03-17

#ifndef __DOLFIN_GENERIC_FUNCTION_H
#define __DOLFIN_GENERIC_FUNCTION_H

#include <dolfin/fem/Coefficient.h>
#include <dolfin/common/Variable.h>

#include <ufc.h>

namespace dolfin
{

class Mesh;
class Cell;

/**
 *  @class  GenericFunction
 *
 *  @brief  This class serves as a base class/interface for implementations of
 *          specific function representations.
 */

class GenericFunction : public Coefficient, public Variable
{

public:

  /// Constructor
  GenericFunction() {}

  /// Destructor
  virtual ~GenericFunction() {}

  //--- UFC INTERFACE ---------------------------------------------------------

  /// Evaluate function at given point in cell
  virtual void evaluate(real* values, const real* coordinates,
                        const ufc::cell& cell) const = 0;

  //--- INTERFACE -------------------------------------------------------------

  /// Evaluate function at given point
  virtual void eval(real* values, const real* x) const = 0;

  /// Return the rank of the value space
  virtual uint rank() const = 0;

  /// Return the dimension of the value space for axis i
  virtual uint dim(uint i) const = 0;

  // Return the value size
  virtual uint value_size() const = 0;

  /// Interpolate function to vertices of mesh
  virtual void interpolate_vertex_values(real* values) const = 0;

  /// Interpolate function to finite element space on cell
  virtual void interpolate(real* coefficients, const ufc::cell& cell,
                           const ufc::finite_element& finite_element,
                           const Cell& dolfin_cell) const = 0;

  /// Interpolate function to finite element space on facet
  virtual void interpolate(real* coefficients, const ufc::cell& cell,
                           const ufc::finite_element& finite_element,
                           const Cell& dolfin_cell, uint facet) const = 0;

  /// Synchronize
  virtual void sync() = 0;

  /// Display basic information
  virtual void disp() const = 0;

  //---------------------------------------------------------------------------

  /// Return the mesh
  virtual Mesh& mesh() const = 0;

};

} /* namespace dolfin */

#endif /* __DOLFIN_GENERIC_FUNCTION_H */
