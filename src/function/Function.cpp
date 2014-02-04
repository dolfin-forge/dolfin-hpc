// Copyright (C) 2007-2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Garth N. Wells 2005-2007.
// Modified by Martin Sandve Alnes 2008.
// Modified by Aurélien Larcher 2013. (extension and partial rewrite)
//
// First added:  2003-11-28
// Last changed: 2013-09-13
//
// The class Function serves as the envelope class and holds a pointer
// to a letter class that is a subclass of GenericFunction. All the
// functionality is handled by the specific implementation (subclass).

#include <dolfin/function/Function.h>

#include <dolfin/elements/ElementLibrary.h>
#include <dolfin/io/File.h>
#include <dolfin/fem/DofMap.h>
#include <dolfin/fem/FiniteElement.h>
#include <dolfin/fem/UFCCell.h>
#include <dolfin/function/ConstantFunction.h>
#include <dolfin/function/DiscreteFunction.h>
#include <dolfin/function/ExpressionFunction.h>
#include <dolfin/function/UFCFunction.h>
#include <dolfin/function/UserFunction.h>

using namespace dolfin;

//-----------------------------------------------------------------------------
Function::Function() :
    Variable("*no name*", "empty function"),
    f(NULL),
    _type(empty),
    _cell(0),
    _facet(-1)
{
  // Do nothing
}
//-----------------------------------------------------------------------------
Function::Function(Mesh& mesh) :
    Variable("*no name*", "user-defined function"),
    f(NULL),
    _type(user),
    _cell(0),
    _facet(-1)
{
  f = new UserFunction(mesh, this);
}
//-----------------------------------------------------------------------------
Function::Function(Mesh& mesh, Expression const& expr) :
    Variable("*no name*", "expression function"),
    f(NULL),
    _type(expression),
    _cell(0),
    _facet(-1)
{
  f = new ExpressionFunction(mesh, expr);
}
//-----------------------------------------------------------------------------
Function::Function(Mesh& mesh, real value) :
    Variable("*no name*", "constant function"),
    f(NULL),
    _type(constant),
    _cell(0),
    _facet(-1)
{
  f = new ConstantFunction(mesh, value);
}
//-----------------------------------------------------------------------------
Function::Function(Mesh& mesh, uint size, real value) :
    Variable("*no name*", "constant function"),
    f(NULL),
    _type(constant),
    _cell(0),
    _facet(-1)
{
  f = new ConstantFunction(mesh, size, value);
}
//-----------------------------------------------------------------------------
Function::Function(Mesh& mesh, const Array<real>& values) :
    Variable("*no name*", "constant function"),
    f(NULL),
    _type(constant),
    _cell(0),
    _facet(-1)
{
  f = new ConstantFunction(mesh, values);
}
//-----------------------------------------------------------------------------
Function::Function(Mesh& mesh, const Array<uint>& shape,
                   const Array<real>& values) :
    Variable("*no name*", "constant function"),
    f(NULL),
    _type(constant),
    _cell(0),
    _facet(-1)
{
  f = new ConstantFunction(mesh, shape, values);
}
//-----------------------------------------------------------------------------
Function::Function(Mesh& mesh, GenericVector& x, Form& form, uint i) :
    Variable("*no name*", "discrete function"),
    f(NULL),
    _type(discrete),
    _cell(0),
    _facet(-1)
{
  f = new DiscreteFunction(mesh, x, form, i);
}
//-----------------------------------------------------------------------------
Function::Function(Mesh& mesh, Form& form, uint i) :
    Variable("*no name*", "discrete function"),
    f(NULL),
    _type(discrete),
    _cell(0),
    _facet(-1)
{
  f = new DiscreteFunction(mesh, form, i);
}
//-----------------------------------------------------------------------------
Function::Function(Mesh& mesh, GenericVector& x,
                   std::string const& finite_element_signature) :
    Variable("*no name*", "discrete function"),
    f(NULL),
    _type(discrete),
    _cell(0),
    _facet(-1)
{
  f = new DiscreteFunction(mesh, x, finite_element_signature,
      DofMap::dofmap_signature(finite_element_signature));
}
//-----------------------------------------------------------------------------
Function::Function(Mesh& mesh, std::string const& finite_element_signature) :
    Variable("*no name*", "discrete function"),
    f(NULL),
    _type(discrete),
    _cell(0),
    _facet(-1)
{
  f = new DiscreteFunction(mesh, finite_element_signature,
      DofMap::dofmap_signature(finite_element_signature));
}
//-----------------------------------------------------------------------------
Function::Function(Mesh& mesh, GenericVector& x,
                   std::string const& finite_element_signature,
                   std::string const& dof_map_signature) :
    Variable("*no name*", "discrete function"),
    f(NULL),
    _type(discrete),
    _cell(0),
    _facet(-1)
{
  f = new DiscreteFunction(mesh, x, finite_element_signature,
      dof_map_signature);
}
//-----------------------------------------------------------------------------
Function::Function(Mesh& mesh, std::string const& finite_element_signature,
                   std::string const& dof_map_signature) :
    Variable("*no name*", "discrete function"),
    f(NULL),
    _type(discrete),
    _cell(0),
    _facet(-1)
{
  f = new DiscreteFunction(mesh, finite_element_signature, dof_map_signature);
}
//-----------------------------------------------------------------------------
Function::Function(const std::string filename) :
    Variable("*no name*", "discrete function from data file"),
    f(NULL),
    _type(empty),
    _cell(0),
    _facet(-1)
{
  File file(filename);
  file >> *this;
}
////-----------------------------------------------------------------------------
Function::Function(SubFunction sub_function) :
    Variable("*no name*", "discrete function"),
    f(NULL),
    _type(discrete),
    _cell(0),
    _facet(-1)
{
  this->f = new DiscreteFunction(sub_function);
}
//-----------------------------------------------------------------------------
Function::Function(Function const& f) :
    f(NULL),
    _type(f.type()),
    _cell(0),
    _facet(-1)
{
  if (f.type() == discrete)
  {
    this->f = new DiscreteFunction(*static_cast<DiscreteFunction*>(f.f));
    rename(f.name(), "discrete function");
  }
  else if (f.type() == constant)
  {
    this->f = new ConstantFunction(*static_cast<ConstantFunction*>(f.f));
    rename(f.name(), "constant function");
  }
  else if (f.type() == empty)
  {
    rename(f.name(), "empty function");
  }
  else
  {
    error("Copy constructor works for discrete,"
	  "constant and empty functions only (so far).");
  }
}
//-----------------------------------------------------------------------------
Function::~Function()
{
  if (f)
  {
    delete f;
  }
}
//-----------------------------------------------------------------------------
void Function::init(Mesh& mesh, real value)
{
  if (f)
  {
    delete f;
  }

  f = new ConstantFunction(mesh, value);
  _type = constant;
}
//-----------------------------------------------------------------------------
void Function::init(Mesh& mesh, Expression const& expr)
{
  if (f)
  {
    delete f;
  }

  f = new ExpressionFunction(mesh, expr);
  _type = expression;
}
//-----------------------------------------------------------------------------
void Function::init(Mesh& mesh, GenericVector& x, Form& form, uint i)
{
  if (f)
  {
    delete f;
  }

  f = new DiscreteFunction(mesh, x, form, i);
  _type = discrete;
}
//-----------------------------------------------------------------------------
void Function::init(Mesh& mesh, Form& form, uint i)
{
  if (f)
  {
    delete f;
  }

  f = new DiscreteFunction(mesh, form, i);
  _type = discrete;
}
//-----------------------------------------------------------------------------
void Function::init(Mesh& mesh, GenericVector& x,
                    std::string const& finite_element_signature)
{
  if (f)
  {
    delete f;
  }

  f = new DiscreteFunction(mesh, x, finite_element_signature,
      DofMap::dofmap_signature(finite_element_signature));
  _type = discrete;
}

//-----------------------------------------------------------------------------
void Function::init(Mesh& mesh, std::string const& finite_element_signature)
{
  if (f)
  {
    delete f;
  }

  f = new DiscreteFunction(mesh, finite_element_signature,
      DofMap::dofmap_signature(finite_element_signature));
  _type = discrete;
}
//-----------------------------------------------------------------------------
void Function::init(Mesh& mesh, GenericVector& x,
                    std::string const& finite_element_signature,
                    std::string const& dof_map_signature)
{
  if (f)
  {
    delete f;
  }

  f = new DiscreteFunction(mesh, x, finite_element_signature,
      dof_map_signature);
  _type = discrete;
}
//-----------------------------------------------------------------------------
void Function::init(Mesh& mesh, std::string const& finite_element_signature,
                    std::string const& dof_map_signature)
{
  if (f)
  {
    delete f;
  }

  f = new DiscreteFunction(mesh, finite_element_signature, dof_map_signature);
  _type = discrete;
}
//-----------------------------------------------------------------------------
Function::Type Function::type() const
{
  return _type;
}
//-----------------------------------------------------------------------------
dolfin::uint Function::rank() const
{
  if (!f)
  {
    error("Function contains no data.");
  }

  return f->rank();
}
//-----------------------------------------------------------------------------
dolfin::uint Function::dim(unsigned int i) const
{
  if (!f)
  {
    error("Function contains no data.");
  }

  return f->dim(i);
}
//-----------------------------------------------------------------------------
Mesh& Function::mesh() const
{
  if (!f)
  {
    error("Function contains no data.");
  }

  return f->mesh;
}
//-----------------------------------------------------------------------------
std::string Function::signature() const
{
  if (!f)
  {
    error("Function contains no data.");
  }

  if (_type != discrete)
  {
    error("A signature can only be returned by discrete functions.");
  }

  return (static_cast<DiscreteFunction*>(f))->signature();
}
//-----------------------------------------------------------------------------
GenericVector& Function::vector() const
{
  if (!f)
    error("Function contains no data.");

  if (_type != discrete)
  {
    error("A vector can only be extracted from discrete functions.");
  }

  return (static_cast<DiscreteFunction*>(f))->vector();
}
//-----------------------------------------------------------------------------
DofMap const& Function::dofmap() const
{
  if (!f)
    error("Function contains no data.");

  if (_type != discrete)
  {
    error("The dofmap can only be extracted from discrete functions.");
  }

  return (static_cast<DiscreteFunction*>(f))->dofmap();
}
//-----------------------------------------------------------------------------
FiniteElement const& Function::finite_element() const
{
  if (!f)
    error("Function contains no data.");

  if (_type != discrete)
  {
    error(
        "The finite element space can only be extracted from discrete functions.");
  }

  return (static_cast<DiscreteFunction*>(f))->finite_element();
}
//-----------------------------------------------------------------------------
void Function::get(real *& values)
{
  if (!f)
    error("Function contains no data.");

  if (_type != discrete)
  {
    error("Values can be retrieved only from discrete functions.");
  }

  return (static_cast<DiscreteFunction*>(f))->get(values);
}
//-----------------------------------------------------------------------------
void Function::set(real *& values)
{
  if (!f)
    error("Function contains no data.");

  if (_type != discrete)
  {
    error("Values can be set only to discrete functions.");
  }

  return (static_cast<DiscreteFunction*>(f))->set(values);
}
//-----------------------------------------------------------------------------
dolfin::uint Function::numSubFunctions() const
{
  if (_type != discrete)
    error("Only discrete functions have sub functions.");

  return static_cast<DiscreteFunction*>(f)->numSubFunctions();
}
//-----------------------------------------------------------------------------
SubFunction Function::operator[](uint i)
{
  if (_type != discrete)
    error("Sub functions can only be extracted from discrete functions.");

  SubFunction sub_function(*static_cast<DiscreteFunction*>(f), i);
  return sub_function;
}
//-----------------------------------------------------------------------------
Function const& Function::operator=(Function& f)
{

  // FIXME: Handle other assignments
  if (f._type != discrete)
    error("Can only handle assignment from discrete functions (for now).");

  // Either create or copy discrete function
  if (_type == discrete)
  {
    *static_cast<DiscreteFunction*>(this->f) =
        *static_cast<DiscreteFunction*>(f.f);
  }
  else
  {
    delete this->f;
    this->f = new DiscreteFunction(*static_cast<DiscreteFunction*>(f.f));
    _type = discrete;
    rename(f.name(), "discrete function");
  }
  return *this;
}
//-----------------------------------------------------------------------------
Function const& Function::operator=(SubFunction sub_function)
{
  if (f)
  {
    delete f;
  }

  f = new DiscreteFunction(sub_function);

  rename("*no name*", "discrete function");
  _type = discrete;

  return *this;
}
//-----------------------------------------------------------------------------
void Function::interpolate(Function const& other_func)
{
  if (f && this->type() == Function::discrete)
  {
    static_cast<DiscreteFunction *>(f)->interpolate(other_func);
  }
  else
  {
    dolfin::error(
        "Function::interpolate(Function const&) can only be called on discrete Function");
  }

}
//-----------------------------------------------------------------------------
void Function::interpolate_vertex_values(real* values)
{
  if (!f)
    error("Function contains no data.");

  f->interpolate_vertex_values(values);
}
//-----------------------------------------------------------------------------
void Function::interpolate(real* coefficients, const ufc::cell& ufc_cell,
                           const ufc::finite_element& finite_element,
                           Cell& cell, int facet)
{
  if (!f)
    error("Function contains no data.");

  // Make current cell and facet are available to user-defined function
  _cell = &cell;
  _facet = facet;

  // Interpolate function
  f->interpolate(coefficients, ufc_cell, finite_element, cell);

  // Make cell and facet unavailable
  _cell = 0;
  _facet = -1;
}
//-----------------------------------------------------------------------------
void Function::sync_ghosts()
{
  if (f)
    f->sync_ghosts();
}
//-----------------------------------------------------------------------------
void Function::eval(real* values, const real* x) const
{
  if (!f)
    error("Function contains no data.");

  // Try scalar version for user-defined function if not overloaded.
  // Otherwise, call eval() function in implementation. Note that we
  // must check if we have a user-defined function or we will go into
  // a loop between Function and UserFunction...
  if (_type == user)
    values[0] = eval(x);
  else
    f->eval(values, x);
}
//-----------------------------------------------------------------------------
dolfin::real Function::eval(const real* x) const
{
  // Try vector-version for non-user-defined function if not
  // overloaded. Otherwise, raise an exception. Note that we must
  // check that we *don't* have a user-defined function or we will go
  // into a loop between Function and UserFunction...

  if (_type != user)
  {
    real values[1] = {0.0};
    eval(values, x);
    return values[0];
  }

  error("Missing eval() for user-defined function (must be overloaded).");
  return 0.0;
}
//-----------------------------------------------------------------------------
const Cell& Function::cell() const
{
  if (!_cell)
    error("Current cell is unknown (only available during assembly).");
  return *_cell;
}
//-----------------------------------------------------------------------------
Point Function::normal() const
{
  return cell().normal(_facet);
}
//-----------------------------------------------------------------------------
int Function::facet() const
{
  return _facet;
}
//-----------------------------------------------------------------------------


