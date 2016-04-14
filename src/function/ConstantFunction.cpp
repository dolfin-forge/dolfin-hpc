// Copyright (C) 2006-2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Martin Sandve Alnes, 2008.
// Modified by Aurélien Larcher, 2013.
//
// First added:  2006-02-09
// Last changed: 2008-07-08

#include <dolfin/function/ConstantFunction.h>

#include <dolfin/log/dolfin_log.h>
#include <dolfin/fem/UFCCell.h>
#include <dolfin/mesh/Mesh.h>

#include <cstring>

namespace dolfin
{

//-----------------------------------------------------------------------------
ConstantFunction::ConstantFunction(const ConstantFunction& f) :
    GenericFunction(),
    mesh_(f.mesh()),
    values(0),
    value_rank(f.value_rank),
    shape(0),
    size(f.size)
{
  values = new real[size];
  shape = new uint[value_rank];
  for (uint i = 0; i < value_rank; i++)
  {
    shape[i] = f.shape[i];
  }
  for (uint i = 0; i < size; i++)
  {
    values[i] = f.values[i];
  }
}

//-----------------------------------------------------------------------------
ConstantFunction::ConstantFunction(Mesh& mesh, real value) :
    GenericFunction(),
    mesh_(mesh),
    values(0),
    value_rank(0),
    shape(0),
    size(1)
{
  shape = new uint[1];
  shape[0] = 1;
  values = new real[1];
  values[0] = value;
}

//-----------------------------------------------------------------------------
ConstantFunction::ConstantFunction(Mesh& mesh, uint size, real value) :
    GenericFunction(),
    mesh_(mesh),
    values(0),
    value_rank((size > 1 ? 1 : 0)),
    shape(0),
    size(size)
{
  shape = new uint[1];
  shape[0] = size;
  values = new real[size];
  for (uint i = 0; i < size; i++)
  {
    values[i] = value;
  }
}

//-----------------------------------------------------------------------------
ConstantFunction::ConstantFunction(Mesh& mesh, const Array<real>& some_values) :
    GenericFunction(),
    mesh_(mesh),
    values(0),
    value_rank((some_values.size() > 1 ? 1 : 0)),
    shape(0),
    size(some_values.size())
{
  if (some_values.empty())
  {
    error("Constant function instantiation from an empty array of values.");
  }
  shape = new uint[1];
  shape[0] = size;
  values = new real[size];
  for (uint i = 0; i < size; i++)
  {
    values[i] = some_values[i];
  }
}

//-----------------------------------------------------------------------------
ConstantFunction::ConstantFunction(Mesh& mesh, const Array<uint>& some_shape,
                                   const Array<real>& some_values) :
    GenericFunction(),
    mesh_(mesh),
    values(0),
    value_rank(0),
    shape(0),
    size(0)
{
  if (some_values.empty())
  {
    error("Constant function instantiation from an empty array of values.");
  }
  if (some_shape.empty())
  {
    error("Constant function instantiation from an empty value_shape.");
  }
  value_rank = some_shape.size();
  shape = new uint[value_rank];
  size = 1;
  for (uint i = 0; i < value_rank; i++)
  {
    shape[i] = some_shape[i];
    size *= shape[i];
  }
  if (size != some_values.size())
  {
    error("Size of values does not match value shape of constant function.");
  }
  values = new real[size];
  for (uint i = 0; i < size; i++)
  {
    values[i] = some_values[i];
  }
}

//-----------------------------------------------------------------------------
ConstantFunction::~ConstantFunction()
{
  delete[] shape;
  delete[] values;
}

//-----------------------------------------------------------------------------
Mesh& ConstantFunction::mesh() const
{
  return mesh_;
}

//-----------------------------------------------------------------------------
dolfin::uint ConstantFunction::rank() const
{
  return value_rank;
}

//-----------------------------------------------------------------------------
dolfin::uint ConstantFunction::dim(uint i) const
{
  if (i > value_rank)
  {
    error("Too large dimension in dim requested to constant function.");
  }
  return shape[i];
}

//-----------------------------------------------------------------------------
void ConstantFunction::interpolate_vertex_values(real* v) const
{
  dolfin_assert(v);

  // Set all vertex values to the constant tensor value
  uint const num_verts = mesh_.size(0);
  for (uint i = 0; i < num_verts; ++i)
  {
    for (uint j = 0; j < size; ++j)
    {
      uint k = i * size + j;
      v[k] = values[j];
    }
  }
}

//-----------------------------------------------------------------------------
void ConstantFunction::interpolate(real* coefficients, const UFCCell& cell,
                                   const ufc::finite_element& finite_element,
                                   const Cell& dolfin_cell) const
{
  dolfin_assert(coefficients);

  // Assert same value shape
  // TODO: Slow to do this for every element, should probably remove later
  dolfin_assert(value_rank == finite_element.value_rank());
  for (uint i = 0; i < value_rank; i++)
  {
    dolfin_assert(shape[i] == finite_element.value_dimension(i));
  }

  // UFC 1.1 version:
  /// Evaluate linear functionals for all dofs on the function f
  finite_element.evaluate_dofs(coefficients, *this, cell);
}

//-----------------------------------------------------------------------------
void ConstantFunction::eval(real* v, const real* x) const
{
  dolfin_assert(v);

  // Set all values to the constant tensor value
  std::memcpy(&v[0], &values[0], sizeof(real) * size);
}

//-----------------------------------------------------------------------------
void ConstantFunction::evaluate(real* v, const real* coordinates,
                                const ufc::cell& cell) const
{
  // Call eval(), cell ignored
  eval(v, coordinates);
}

//-----------------------------------------------------------------------------
void ConstantFunction::disp() const
{
  cout << "ConstantFunction" << endl;
  cout << "----------------" << endl;

  // Begin indentation
  begin("");
  GenericFunction::disp();
  // End indentation
  end();
  skip();
}
//-----------------------------------------------------------------------------

}
