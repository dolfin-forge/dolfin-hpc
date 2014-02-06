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
#include <dolfin/function/UserFunction.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
std::string Function::type2string(Function::Type type)
{
  switch (type)
    {
    case empty:
      return "empty";
    case constant:
      return "constant";
    case discrete:
      return "discrete";
    case expression:
      return "expression";
    case user:
      return "user";
    default:
      error("Unknown function type: %d.", type);
      break;
    }

  return "";
}
//-----------------------------------------------------------------------------
Function::Function() :
    Variable("*no name*", "empty function"),
    f_(NULL),
    type_(empty),
    cell_(0),
    facet_(-1)
{
  // Do nothing
}
//-----------------------------------------------------------------------------
Function::Function(Mesh& mesh, real value) :
    Variable("*no name*", "constant function"),
    f_(NULL),
    type_(constant),
    cell_(0),
    facet_(-1)
{
  f_ = new ConstantFunction(mesh, value);
}
//-----------------------------------------------------------------------------
Function::Function(Mesh& mesh, uint size, real value) :
    Variable("*no name*", "constant function"),
    f_(NULL),
    type_(constant),
    cell_(0),
    facet_(-1)
{
  f_ = new ConstantFunction(mesh, size, value);
}
//-----------------------------------------------------------------------------
Function::Function(Mesh& mesh, const Array<real>& values) :
    Variable("*no name*", "constant function"),
    f_(NULL),
    type_(constant),
    cell_(0),
    facet_(-1)
{
  f_ = new ConstantFunction(mesh, values);
}
//-----------------------------------------------------------------------------
Function::Function(Mesh& mesh, const Array<uint>& shape,
                   const Array<real>& values) :
    Variable("*no name*", "constant function"),
    f_(NULL),
    type_(constant),
    cell_(0),
    facet_(-1)
{
  f_ = new ConstantFunction(mesh, shape, values);
}
//-----------------------------------------------------------------------------
Function::Function(Mesh& mesh, GenericVector& x, Form& form, uint i) :
    Variable("*no name*", "discrete function"),
    f_(NULL),
    type_(discrete),
    cell_(0),
    facet_(-1)
{
  f_ = new DiscreteFunction(mesh, x, form, i);
}
//-----------------------------------------------------------------------------
Function::Function(Mesh& mesh, Form& form, uint i) :
    Variable("*no name*", "discrete function"),
    f_(NULL),
    type_(discrete),
    cell_(0),
    facet_(-1)
{
  f_ = new DiscreteFunction(mesh, form, i);
}
//-----------------------------------------------------------------------------
Function::Function(Mesh& mesh, GenericVector& x,
                   std::string const& finite_element_signature) :
    Variable("*no name*", "discrete function"),
    f_(NULL),
    type_(discrete),
    cell_(0),
    facet_(-1)
{
  f_ = new DiscreteFunction(mesh, x, finite_element_signature);
}
//-----------------------------------------------------------------------------
Function::Function(Mesh& mesh, GenericVector& x,
                   std::string const& finite_element_signature,
                   std::string const& dof_map_signature) :
    Variable("*no name*", "discrete function"),
    f_(NULL),
    type_(discrete),
    cell_(0),
    facet_(-1)
{
  f_ = new DiscreteFunction(mesh, x, finite_element_signature,
                            dof_map_signature);
}
//-----------------------------------------------------------------------------
Function::Function(Mesh& mesh, std::string const& finite_element_signature,
                   std::string const& dof_map_signature) :
    Variable("*no name*", "discrete function"),
    f_(NULL),
    type_(discrete),
    cell_(0),
    facet_(-1)
{
  f_ = new DiscreteFunction(mesh, finite_element_signature, dof_map_signature);
}
//-----------------------------------------------------------------------------
Function::Function(Mesh& mesh, std::string const& finite_element_signature) :
    Variable("*no name*", "discrete function"),
    f_(NULL),
    type_(discrete),
    cell_(0),
    facet_(-1)
{
  f_ = new DiscreteFunction(mesh, finite_element_signature);
}
#if ENABLE_UFL
//-----------------------------------------------------------------------------
Function::Function(Mesh& mesh, ufl::FiniteElementBase const& finite_element) :
    Variable("*no name*", "discrete function"),
    f_(NULL),
    type_(discrete),
    cell_(0),
    facet_(-1)
{
  f_ = new DiscreteFunction(mesh, finite_element);
}
#endif
////-----------------------------------------------------------------------------
Function::Function(SubFunction sub_function) :
    Variable("*no name*", "discrete function"),
    f_(NULL),
    type_(discrete),
    cell_(0),
    facet_(-1)
{
  this->f_ = new DiscreteFunction(sub_function);
}
//-----------------------------------------------------------------------------
Function::Function(Mesh& mesh, Expression const& expr) :
    Variable("*no name*", "expression function"),
    f_(NULL),
    type_(expression),
    cell_(0),
    facet_(-1)
{
  f_ = new ExpressionFunction(mesh, expr);
}
//-----------------------------------------------------------------------------
Function::Function(Mesh& mesh) :
    Variable("*no name*", "user-defined function"),
    f_(NULL),
    type_(user),
    cell_(0),
    facet_(-1)
{
  f_ = new UserFunction(mesh, this);
}
//-----------------------------------------------------------------------------
Function const& Function::operator=(Function& f)
{

  // FIXME: Handle other assignments
  if (f.type_ != discrete)
    error("Can only handle assignment from discrete functions (for now).");

  // Either create or copy discrete function
  if (type_ == discrete)
  {
    *static_cast<DiscreteFunction*>(this->f_) =
        *static_cast<DiscreteFunction*>(f.f_);
  }
  else
  {
    delete this->f_;
    this->f_ = new DiscreteFunction(*static_cast<DiscreteFunction*>(f.f_));
    type_ = discrete;
    rename(f.name(), "discrete function");
  }
  return *this;
}
//-----------------------------------------------------------------------------
Function::Function(const std::string filename) :
    Variable("*no name*", "discrete function from data file"),
    f_(NULL),
    type_(empty),
    cell_(0),
    facet_(-1)
{
  File file(filename);
  file >> *this;
}
//-----------------------------------------------------------------------------
Function::Function(Function const& f) :
    f_(NULL),
    type_(f.type()),
    cell_(0),
    facet_(-1)
{
  if (f.type() == discrete)
  {
    this->f_ = new DiscreteFunction(*static_cast<DiscreteFunction*>(f.f_));
    rename(f.name(), "discrete function");
  }
  else if (f.type() == constant)
  {
    this->f_ = new ConstantFunction(*static_cast<ConstantFunction*>(f.f_));
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
  if (f_)
  {
    delete f_;
  }
}
//-----------------------------------------------------------------------------
void Function::init(Mesh& mesh, real value)
{
  if (f_)
  {
    delete f_;
  }

  f_ = new ConstantFunction(mesh, value);
  type_ = constant;
}
//-----------------------------------------------------------------------------
void Function::init(Mesh& mesh, GenericVector& x, Form& form, uint i)
{
  if (f_)
  {
    delete f_;
  }

  f_ = new DiscreteFunction(mesh, x, form, i);
  type_ = discrete;
}
//-----------------------------------------------------------------------------
void Function::init(Mesh& mesh, Form& form, uint i)
{
  if (f_)
  {
    delete f_;
  }

  f_ = new DiscreteFunction(mesh, form, i);
  type_ = discrete;
}
//-----------------------------------------------------------------------------
void Function::init(Mesh& mesh, GenericVector& x,
                    std::string const& finite_element_signature,
                    std::string const& dof_map_signature)
{
  if (f_)
  {
    delete f_;
  }

  f_ = new DiscreteFunction(mesh, x, finite_element_signature,
                            dof_map_signature);
  type_ = discrete;
}
//-----------------------------------------------------------------------------
void Function::init(Mesh& mesh, std::string const& finite_element_signature,
                    std::string const& dof_map_signature)
{
  if (f_)
  {
    delete f_;
  }

  f_ = new DiscreteFunction(mesh, finite_element_signature, dof_map_signature);
  type_ = discrete;
}
//-----------------------------------------------------------------------------
void Function::init(Mesh& mesh, GenericVector& x,
                    std::string const& finite_element_signature)
{
  if (f_)
  {
    delete f_;
  }

  f_ = new DiscreteFunction(mesh, x, finite_element_signature);
  type_ = discrete;
}
//-----------------------------------------------------------------------------
void Function::init(Mesh& mesh, std::string const& finite_element_signature)
{
  if (f_)
  {
    delete f_;
  }

  f_ = new DiscreteFunction(mesh, finite_element_signature);
  type_ = discrete;
}
#if ENABLE_UFL
//-----------------------------------------------------------------------------
void Function::init(Mesh& mesh, ufl::FiniteElementBase const& finite_element)
{
  if (f_)
  {
    delete f_;
  }

  f_ = new DiscreteFunction(mesh, finite_element);
  type_ = discrete;
}
#endif
//-----------------------------------------------------------------------------
void Function::init(Mesh& mesh, Expression const& expr)
{
  if (f_)
  {
    delete f_;
  }

  f_ = new ExpressionFunction(mesh, expr);
  type_ = expression;
}
//--- UFC INTERFACE -----------------------------------------------------------
//-----------------------------------------------------------------------------
void Function::evaluate(real* values, const real* coordinates,
                        const ufc::cell& cell) const
{
  if (!f_)
  {
    error("Function contains no data.");
  }
  f_->evaluate(values, coordinates, cell);
}
//--- COMPOSITION GenericFunction ---------------------------------------------
//-----------------------------------------------------------------------------
uint Function::rank() const
{
  if (!f_)
  {
    error("Function contains no data.");
  }

  return f_->rank();
}
//-----------------------------------------------------------------------------
uint Function::dim(unsigned int i) const
{
  if (!f_)
  {
    error("Function contains no data.");
  }

  return f_->dim(i);
}
//-----------------------------------------------------------------------------
void Function::interpolate_vertex_values(real* values)
{
  if (!f_)
    error("Function contains no data.");

  f_->interpolate_vertex_values(values);
}
//-----------------------------------------------------------------------------
void Function::interpolate(real* coefficients, const ufc::cell& ufc_cell,
                           const ufc::finite_element& finite_element,
                           Cell& cell, int facet)
{
  if (!f_)
    error("Function contains no data.");

  // Make current cell and facet are available to user-defined function
  cell_ = &cell;
  facet_ = facet;

  // Interpolate function
  f_->interpolate(coefficients, ufc_cell, finite_element, cell);

  // Make cell and facet unavailable
  cell_ = 0;
  facet_ = -1;
}
//-----------------------------------------------------------------------------
void Function::eval(real* values, const real* x) const
{
  if (!f_)
    error("Function contains no data.");

  f_->eval(values, x);
}
//-----------------------------------------------------------------------------
void Function::sync_ghosts()
{
  if (f_)
    f_->sync_ghosts();
}
//-----------------------------------------------------------------------------
void Function::disp() const
{
  cout << "Function" << endl;
  cout << "------- " << endl;

  // Begin indentation
  begin("");
  cout << "Type: " << this->type() << " ("
      << Function::type2string(this->type()) << ")" << endl;
  if (f_ != NULL)
  {
    f_->disp();
  }
  // End indentation
  end();
}
//-----------------------------------------------------------------------------
Mesh& Function::mesh() const
{
  if (!f_)
  {
    error("Function contains no data.");
  }

  return f_->mesh;
}
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
Function::Type Function::type() const
{
  return type_;
}
//--- Wrapper Facade for DiscreteFunction -------------------------------------
//-----------------------------------------------------------------------------
std::string Function::signature() const
{
  if (!f_)
  {
    error("Function contains no data.");
  }

  if (type_ != discrete)
  {
    error("A signature can only be returned by discrete functions.");
  }

  return (static_cast<DiscreteFunction*>(f_))->signature();
}
//-----------------------------------------------------------------------------
GenericVector& Function::vector() const
{
  if (!f_)
    error("Function contains no data.");

  if (type_ != discrete)
  {
    error("A vector can only be extracted from discrete functions.");
  }

  return (static_cast<DiscreteFunction*>(f_))->vector();
}
//-----------------------------------------------------------------------------
FiniteElementSpace const& Function::space() const
{
  if (!f_)
    error("Function contains no data.");

  if (type_ != discrete)
  {
    error("The dofmap can only be extracted from discrete functions.");
  }

  return (static_cast<DiscreteFunction*>(f_))->space();
}
//-----------------------------------------------------------------------------
DofMap const& Function::dofmap() const
{
  if (!f_)
    error("Function contains no data.");

  if (type_ != discrete)
  {
    error("The dofmap can only be extracted from discrete functions.");
  }

  return (static_cast<DiscreteFunction*>(f_))->dofmap();
}
//-----------------------------------------------------------------------------
FiniteElement const& Function::finite_element() const
{
  if (!f_)
    error("Function contains no data.");

  if (type_ != discrete)
  {
    error("The finite element space can only be extracted from discrete "
          "functions.");
  }

  return (static_cast<DiscreteFunction*>(f_))->finite_element();
}
//-----------------------------------------------------------------------------
uint Function::num_sub_functions() const
{
  if (type_ != discrete)
    error("Only discrete functions have sub functions.");

  return static_cast<DiscreteFunction*>(f_)->numSubFunctions();
}
//-----------------------------------------------------------------------------
void Function::interpolate(Function const& other_func)
{
  if (f_ && this->type() == Function::discrete)
  {
    static_cast<DiscreteFunction *>(f_)->interpolate(other_func);
  }
  else
  {
    dolfin::error("Function::interpolate(Function const&) can only be called "
                  "on discrete Function");
  }

}
//-----------------------------------------------------------------------------
void Function::get(real *& values)
{
  if (!f_)
    error("Function contains no data.");

  if (type_ != discrete)
  {
    error("Values can be retrieved only from discrete functions.");
  }

  return (static_cast<DiscreteFunction*>(f_))->get(values);
}
//-----------------------------------------------------------------------------
void Function::set(real *& values)
{
  if (!f_)
    error("Function contains no data.");

  if (type_ != discrete)
  {
    error("Values can be set only to discrete functions.");
  }

  return (static_cast<DiscreteFunction*>(f_))->set(values);
}
//-----------------------------------------------------------------------------
SubFunction Function::operator[](uint i)
{
  if (type_ != discrete)
    error("Sub functions can only be extracted from discrete functions.");

  SubFunction sub_function(*static_cast<DiscreteFunction*>(f_), i);
  return sub_function;
}
//--- PROTECTED ---------------------------------------------------------------
//-----------------------------------------------------------------------------
Cell const& Function::cell() const
{
  if (!cell_)
    error("Current cell is unknown (only available during assembly).");
  return *cell_;
}
//-----------------------------------------------------------------------------
Point Function::normal() const
{
  return cell().normal(facet_);
}
//-----------------------------------------------------------------------------
int Function::facet() const
{
  return facet_;
}
//-----------------------------------------------------------------------------

}
