// Copyright (C) 2007-2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Garth N. Wells 2005-2007.
// Modified by Martin Sandve Alnes 2008.
// Modified by Aurélien Larcher 2013-2014. (extension and partial rewrite)
//
// First added:  2003-11-28
// Last changed: 2014-02-06
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
#include <dolfin/mesh/Vertex.h>

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
Function::Function(Mesh& mesh) :
    Variable("*no name*", "empty function"),
    f_(NULL),
    type_(user),
    cell_(0),
    facet_(-1)
{
  f_ = new UserFunction(mesh, this);
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
Function::Function(GenericVector& x, Form& form, uint i) :
    Variable("*no name*", "discrete function"),
    f_(NULL),
    type_(discrete),
    cell_(0),
    facet_(-1)
{
  f_ = new DiscreteFunction(x, form, i);
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
Function::Function(Form& form, uint i) :
    Variable("*no name*", "discrete function"),
    f_(NULL),
    type_(discrete),
    cell_(0),
    facet_(-1)
{
  f_ = new DiscreteFunction(form, i);
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
Function::Function(GenericVector& x, FiniteElementSpace const& space) :
    Variable("*no name*", "discrete function"),
    f_(NULL),
    type_(discrete),
    cell_(0),
    facet_(-1)
{
  f_ = new DiscreteFunction(x, space);
}
//-----------------------------------------------------------------------------
Function::Function(FiniteElementSpace const& space) :
    Variable("*no name*", "discrete function"),
    f_(NULL),
    type_(discrete),
    cell_(0),
    facet_(-1)
{
  f_ = new DiscreteFunction(space);
}
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
//-----------------------------------------------------------------------------
Function::Function(Mesh& mesh, std::string const& element,
                   std::string const& dofmap) :
    Variable("*no name*", "discrete function"),
    f_(NULL),
    type_(discrete),
    cell_(0),
    facet_(-1)
{
  ufc::finite_element * f = ElementLibrary::create_finite_element(element);
  ufc::dofmap * d = ElementLibrary::create_dof_map(dofmap);
  FiniteElementSpace space(mesh, *f, *d, true);
  f_ = new DiscreteFunction(space);
}
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
Function const& Function::operator=(SubFunction sub_function)
{
  if (f_)
  {
    delete f_;
  }

  f_ = new DiscreteFunction(sub_function);

  rename("*no name*", "discrete function");
  type_ = discrete;

  return *this;
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
Function const& Function::operator=(Function& f)
{

  // FIXME: Handle other assignments
  if (f.type_ != discrete)
  {
    error("Can only handle assignment from discrete functions (for now).");
  }

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
  else if (f.type() == user)
  {
    this->f_ = new UserFunction(*static_cast<UserFunction*>(f.f_));
    rename(f.name(), "user function");
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
  delete f_;
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
void Function::init(Mesh& mesh, uint i, real value)
{
  if (f_)
  {
    delete f_;
  }

  f_ = new ConstantFunction(mesh, i, value);
  type_ = constant;
}
//-----------------------------------------------------------------------------
void Function::init(GenericVector& x, Form& form, uint i)
{
  if (f_)
  {
    delete f_;
  }

  //FIXME: Assumes one mesh per form
  f_ = new DiscreteFunction(form.mesh(), x, form, i);
  type_ = discrete;
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
void Function::init(Form& form, uint i)
{
  if (f_)
  {
    delete f_;
  }

  //FIXME: Assumes one mesh per form
  f_ = new DiscreteFunction(form.mesh(), form, i);
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
void Function::init(GenericVector& x, FiniteElementSpace const& space)
{
  if (f_)
  {
    delete f_;
  }
  //TODO: Check mesh consistency
  f_ = new DiscreteFunction(x, space);
  type_ = discrete;
}
//-----------------------------------------------------------------------------
void Function::init(FiniteElementSpace const& space)
{
  if (f_)
  {
    delete f_;
  }
  //TODO: Check mesh consistency
  f_ = new DiscreteFunction(space);
  type_ = discrete;
}
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
//-----------------------------------------------------------------------------
void Function::init(Mesh& mesh, std::string const& element,
                    std::string const& dofmap)
{
  if (f_)
  {
    delete f_;
  }

  ufc::finite_element * f = ElementLibrary::create_finite_element(element);
  ufc::dofmap * d = ElementLibrary::create_dof_map(dofmap);
  FiniteElementSpace space(mesh, *f, *d, true);
  f_ = new DiscreteFunction(space);
  type_ = discrete;
}
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
//-----------------------------------------------------------------------------
Function::Type Function::type() const
{
  return type_;
}
//--- Wrapper Facade for DiscreteFunction -------------------------------------
//-----------------------------------------------------------------------------
GenericVector& Function::vector() const
{
  if (type_ != discrete)
  {
    error("A vector can only be extracted from discrete functions.");
  }

  return (static_cast<DiscreteFunction*>(f_))->vector();
}
//-----------------------------------------------------------------------------
FiniteElementSpace const& Function::space() const
{
  if (type_ != discrete)
  {
    error("The dofmap can only be extracted from discrete functions.");
  }

  return (static_cast<DiscreteFunction*>(f_))->space();
}
//-----------------------------------------------------------------------------
std::string const Function::signature() const
{
  if (type_ != discrete)
  {
    error("A signature can only be returned by discrete functions.");
  }

  return (static_cast<DiscreteFunction*>(f_))->signature();
}
//-----------------------------------------------------------------------------
uint const Function::num_sub_functions() const
{
  if (type_ != discrete) error("Only discrete functions have sub functions.");

  return static_cast<DiscreteFunction*>(f_)->num_sub_functions();
}
//-----------------------------------------------------------------------------
void Function::interpolate(Function const& other_func)
{
  if (type_ == discrete)
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
real * Function::create_block() const
{
  if (type_ != discrete)
  {
    error("Block array can be created only from discrete functions.");
  }

  return (static_cast<DiscreteFunction*>(f_))->create_block();
}
//-----------------------------------------------------------------------------
void Function::get_block(real *& values) const
{
  if (type_ != discrete)
  {
    error("Values can be retrieved only from discrete functions.");
  }

  return (static_cast<DiscreteFunction*>(f_))->get_block(values);
}
//-----------------------------------------------------------------------------
void Function::set_block(real *& values)
{
  if (type_ != discrete)
  {
    error("Values can be set only to discrete functions.");
  }

  return (static_cast<DiscreteFunction*>(f_))->set_block(values);
}
//-----------------------------------------------------------------------------
void Function::add_block(real *& values)
{
  if (type_ != discrete)
  {
    error("Values can be set only to discrete functions.");
  }

  return (static_cast<DiscreteFunction*>(f_))->add_block(values);
}
//-----------------------------------------------------------------------------
SubFunction Function::operator[](uint i)
{
  if (type_ != discrete)
  {
    error("Sub functions can only be extracted from discrete functions.");
  }

  SubFunction sub_function(*static_cast<DiscreteFunction*>(f_), i);
  return sub_function;
}

//-----------------------------------------------------------------------------
Array<Function *> Function::decompose()
{
  if (type_ != discrete)
  {
    error("Only discrete functions can be decomposed with decompose().");
  }

  Array<Function *> leaf_functions;

  //TODO: This implementation should belong to DiscreteFunction: let us consider
  //      that things are acceptable as long as high-level functions are OK.
  FiniteElementSpace const& space = this->space();
  DofMap const& dm = space.dofmap();

  Array<FiniteElementSpace *> spaces = space.flatten();
  for (uint s = 0; s < spaces.size(); ++s)
  {
    leaf_functions.push_back(new Function(*(spaces[s])));
  }

  if (space.is_vertex_based())
  {
    //NOTE: This implementation is based on the assumption that the dofmap for
    //      a CG1 function is indexed by the global indices of vertices.
    Mesh& mesh = this->mesh();
    MeshFunction<bool> marked(mesh, 0);
    real dof_value;
    uint * indices = new uint[dm.local_dimension()];
    CellIterator c(mesh);
    UFCCell ufc_cell(*c);
    for (; !c.end(); ++c)
    {
      ufc_cell.update(*c);

      for (VertexIterator v(*c); !v.end(); ++v)
      {

        uint * cvi = c->entities(0);
        uint ci = 0;
        for (ci = 0; ci < c->numEntities(0); ++ci)
        {
          if (cvi[ci] == v->index())
          {
            break;
          }
        }
        if (!v->is_ghost() && !marked.get(*v))
        {
          uint new_index = mesh.distdata().get_vertex_global(v->index());
          for (uint i = 0; i < leaf_functions.size(); ++i)
          {
            this->vector().get(&dof_value, 1, &indices[ci]);
            leaf_functions[i]->vector().set(&dof_value, 1, &new_index);
          }

          marked.set(*v, true);
          continue;

        }
      }
    }
  }
  else
  {
    error("Function decomposition is only implemented for vertex-based spaces");
  }

  for (uint i = 0; i < leaf_functions.size(); ++i)
  {
    leaf_functions[i]->sync_ghosts();
  }
  return leaf_functions;
}

//--- PROTECTED ---------------------------------------------------------------
//-----------------------------------------------------------------------------
Cell const& Function::cell() const
{
  if (!cell_)
  {
    error("Current cell is unknown (only available during assembly).");
  }
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
