// Copyright (C) 2005-2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2005-11-28
// Last changed: 2008-03-17

#ifndef __DOLFIN_GENERIC_FUNCTION_H
#define __DOLFIN_GENERIC_FUNCTION_H

#include <dolfin/common/types.h>
#include <dolfin/common/Array.h>
#include <dolfin/fem/Coefficient.h>

#include <ufc.h>

namespace dolfin
{

class Mesh;
class Cell;
class UFCCell;

/// This class serves as a base class/interface for implementations
/// of specific function representations.

class GenericFunction : public Coefficient
{
public:

  /// Constructor
  GenericFunction()
  {
  }

  /// Destructor
  virtual ~GenericFunction()
  {
  }

  //--- UFC INTERFACE ---------------------------------------------------------

  /// Evaluate function at given point in cell
  virtual void evaluate(real* values, const real* coordinates,
                        const ufc::cell& cell) const = 0;

  //--- INTERFACE -------------------------------------------------------------

  /// Return the mesh
  virtual Mesh& mesh() const = 0;

  /// Return the rank of the value space
  virtual uint rank() const = 0;

  /// Return the dimension of the value space for axis i
  virtual uint dim(uint i) const = 0;

  /// Interpolate function to vertices of mesh
  virtual void interpolate_vertex_values(real* values) const = 0;

  /// Interpolate function to finite element space on cell
  virtual void interpolate(real* coefficients, const UFCCell& cell,
                           const ufc::finite_element& finite_element,
                           const Cell& dolfin_cell) const = 0;

  /// Interpolate function to finite element space on facet
  virtual void interpolate(real* coefficients, const UFCCell& cell,
                           const ufc::finite_element& finite_element,
                           const Cell& dolfin_cell, uint facet) const
  {
    interpolate(coefficients, cell, finite_element, dolfin_cell);
  }

  /// Evaluate function at given point
  virtual void eval(real* values, const real* x) const = 0;

  /// Display basic information
  virtual void disp() const = 0;

  /// Synchronize
  inline void sync() { this->sync_ghosts(); }
  virtual void sync_ghosts() = 0;

};

//-----------------------------------------------------------------------------
inline void GenericFunction::disp() const
{
  cout << "GenericFunction" << endl;
  cout << "---------------" << endl;

  // Begin indentation
  begin("");
  cout << "Value rank            : " << this->rank() << endl;
  cout << "Value dimension       : " << this->dim(0) << endl;
  // End indentation
  end();
  skip();
}

}

#endif
