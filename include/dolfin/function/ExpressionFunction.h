// Copyright (C) 2013 Aurélien Larcher
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2013-05-25 (merged from branch larcher)
// Last changed: 2013-05-25

#ifndef __DOLFIN_EXPRESSION_FUNCTION_H
#define __DOLFIN_EXPRESSION_FUNCTION_H

#include "GenericFunction.h"
#include "Function.h"
#include "Expression.h"

namespace dolfin
{

  ///

  class ExpressionFunction : public GenericFunction
  {

  public:

    /// Create user-defined function
    ExpressionFunction(Mesh& mesh, Expression const& expr);

    /// Destructor
    ~ExpressionFunction();

    //--- INTERFACE -----------------------------------------------------------

    /// Return the mesh
    Mesh& mesh() const;

    /// Return the rank of the value space
    uint rank() const;

    /// Return the dimension of the value space for axis i
    uint dim(uint i) const;

    /// Interpolate function to vertices of mesh
    void interpolate_vertex_values(real* values) const;

    /// Interpolate function to finite element space on cell
    void interpolate(real coefficients[],
                     const ufc::cell& cell,
                     const ufc::finite_element& finite_element,
                     const Cell& dolfin_cell) const;

    /// Evaluate function at given point
    void eval(real* values, const real* x) const;

    /// Evaluate function at given point in cell (UFC function interface)
    void evaluate(real* values,
                  const real* coordinates,
                  const ufc::cell& cell) const;

    /// Display basic information
    void disp() const;

    /// Synchronize ghosted entries
    void sync_ghosts() { return; }

  private:

    //
    Mesh& mesh_;

    // Expression
    Expression const& e;

  };

}

#endif
