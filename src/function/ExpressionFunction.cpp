// Copyright (C) 2013 Aurélien Larcher
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2013-05-25 (merged from branch larcher)
// Last changed: 2013-05-25

#include <dolfin/function/ExpressionFunction.h>

#include <dolfin/function/Function.h>
#include <dolfin/mesh/Cell.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/Vertex.h>

namespace dolfin {

//-----------------------------------------------------------------------------
ExpressionFunction::ExpressionFunction(Mesh& mesh, Expression const& expr) :
    GenericFunction(mesh),
    ufc::function(),
    e(expr)
{
  // Do nothing
}
//-----------------------------------------------------------------------------
ExpressionFunction::~ExpressionFunction()
{
  // Do nothing
}
//-----------------------------------------------------------------------------
dolfin::uint ExpressionFunction::rank() const
{
  return e.rank();
}
//-----------------------------------------------------------------------------
dolfin::uint ExpressionFunction::dim(uint i) const
{
  return e.dim(i);
}
//-----------------------------------------------------------------------------
void ExpressionFunction::interpolate_vertex_values(real* values) const
{
  dolfin_assert(values);

  // Compute size of value (number of entries in tensor value)
  uint size = 1;
  for (uint i = 0; i < e.rank(); i++)
  {
    size *= e.dim(i);
  }

  // Call overloaded eval function at each vertex
  real * local_values = new real[size];

  for (VertexIterator vertex(mesh); !vertex.end(); ++vertex)
  {
    // Evaluate at function at vertex
    e.eval(local_values, vertex->x());

    // Copy values to array of vertex values
    for (uint i = 0; i < size; i++)
    {
      values[i * mesh.numVertices() + vertex->index()] = local_values[i];
    }
  }
  delete[] local_values;
}
//-----------------------------------------------------------------------------
void ExpressionFunction::interpolate(real* coefficients, const ufc::cell& cell,
                const ufc::finite_element& finite_element,
                const Cell& dolfin_cell) const
{
  dolfin_assert(coefficients);

  // Evaluate each dof to get coefficients for nodal basis expansion
  for (uint i = 0; i < finite_element.space_dimension(); i++)
  {
    coefficients[i] = finite_element.evaluate_dof(i, *this, cell);
  }
}
//-----------------------------------------------------------------------------
void ExpressionFunction::eval(real* values, const real* x) const
{
  e.eval(values, x);
}
//-----------------------------------------------------------------------------
void ExpressionFunction::evaluate(real* values, const real* coordinates,
              const ufc::cell& cell) const
{
  dolfin_assert(values);
  dolfin_assert(coordinates);

  // Compute size of value (number of entries in tensor value)
  uint size = 1;
  for (uint i = 0; i < e.rank(); i++)
  {
    size *= e.dim(i);
  }

  e.eval(values, coordinates);

}
//-----------------------------------------------------------------------------

}
