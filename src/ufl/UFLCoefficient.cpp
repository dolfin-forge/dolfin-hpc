// Copyright (C) 2014 Bärbel Janssen.
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/ufl/UFLCoefficient.h>
#include <dolfin/ufl/UFLVectorElement.h>
#include <dolfin/ufl/UFLTensorElement.h>

namespace ufl
{

//-----------------------------------------------------------------------------
CoefficientBase::CoefficientBase(std::string const& name,
                                 FiniteElementSpace const& fe,
                                 dolfin::uint const& c) :
    Expression(name),
    finite_element_(&fe),
    count_(c)
{
}

//-----------------------------------------------------------------------------
CoefficientBase::CoefficientBase(std::string const& name, Cell const& cell,
                                 dolfin::uint const& c) :
    Expression(name),
    finite_element_(new FiniteElement(Family::R, cell, 0)),
    count_(c)
{
}

//-----------------------------------------------------------------------------
CoefficientBase::CoefficientBase(std::string const& name, Cell const& cell,
                                 dolfin::uint const& dim, dolfin::uint const& c) :
    Expression(name),
    finite_element_(new VectorElement(Family::R, cell, 0, dim)),
    count_(c)
{
}

//-----------------------------------------------------------------------------
CoefficientBase::CoefficientBase(
    std::string const& name, Cell const& cell, ValueArray const&,
    dolfin::_ordered_map<dolfin::uint, dolfin::uint> const&, dolfin::uint const& c) :
    Expression(name),
    finite_element_(
        new TensorElement(Family::R, cell, 0, cell.geometric_dimension())),
    count_(c)
{
}

//-----------------------------------------------------------------------------
CoefficientBase::CoefficientBase(std::string const& name, repr_t const& repr) :
    Expression(name, repr),
    finite_element_(
        name == "Coefficient" ? FiniteElementSpace::create(arg(0)) :
        name == "Constant" ? new FiniteElement(Family::R, Cell(arg(0)), 0) :
        name == "VectorConstant" ?
            new VectorElement(Family::R, Cell(arg(0)), 0,
                              type<dolfin::uint>(arg(1))) :
        name == "TensorConstant" ?
            new TensorElement(Family::R, Cell(arg(0)), 0,
                              Cell(arg(0)).geometric_dimension()) :
            FiniteElementSpace::create(arg(0))),
    count_(args().back())
{
//    finite_element_ = (name == "Coefficient" ? FiniteElementSpace::create(arg(0)) :
//       name == "Constant" ? new FiniteElement(Family::R, Cell(arg(0)), 0) :
//       name == "VectorConstant" ? new VectorElement(Family::R, Cell(arg(0)), 0, type<dolfin::uint>(arg(1))) :
//       name == "TensorConstant" ? new TensorElement(Family::R, Cell(arg(0)), 0, Cell(arg(0)).geometric_dimension()) :
//       FiniteElementSpace::create(arg(0)));
}

//-----------------------------------------------------------------------------
CoefficientBase::~CoefficientBase()
{
}

//-----------------------------------------------------------------------------
CoefficientBase const* CoefficientBase::create(Object::repr_t const& repr)
{
  std::string name = Class::make_name(repr);
  if (name == "Coefficient")
  {
    return new Coefficient(repr);
  }
  else if (name == "Constant")
  {
    return new Constant(repr);
  }
  else if (name == "VectorConstant")
  {
    return new VectorConstant(repr);
  }
  else if (name == "TensorConstant")
  {
    return new TensorConstant(repr);
  }
  else
  {
    error("Unknown type of ufl::CoefficientBase: '" + name + "'");
  }

  return nullptr;
}

//-----------------------------------------------------------------------------
std::vector<Expression const *> const CoefficientBase::operands() const
{
  error("Not yet implemented.");
  std::vector<Expression const*> expressions_;
  return expressions_;
}

//-----------------------------------------------------------------------------
FiniteElementSpace const& CoefficientBase::element() const
{
  return *finite_element_;
}

//-----------------------------------------------------------------------------
ValueArray const CoefficientBase::shape() const
{
  return finite_element_->value_shape();
}

//-----------------------------------------------------------------------------
tuple<Index> const CoefficientBase::free_indices() const
{
  return tuple<Index>();
}

//-----------------------------------------------------------------------------
dict<IndexBase, type<dolfin::uint> > const CoefficientBase::index_dimensions() const
{
  return dict<IndexBase, type<dolfin::uint> >();
}

//-----------------------------------------------------------------------------
std::vector<std::vector<std::vector<dolfin::real> > > const CoefficientBase::evaluate(
    dolfin::uint,
    std::vector<std::vector<std::vector<dolfin::real> > > const&,
    ufc::cell const& , std::vector<dolfin::real*> const&,
    const double * const *) const
{
  std::vector<std::vector<std::vector<dolfin::real> > > const new_vals0;
  return new_vals0;
}

//-----------------------------------------------------------------------------
Cell const CoefficientBase::cell() const
{
  return finite_element_->cell();
}

//-----------------------------------------------------------------------------
bool CoefficientBase::is_cellwise_constant() const
{
  return false;
}

//-----------------------------------------------------------------------------
Coefficient::Coefficient(FiniteElementSpace const& fe, dolfin::uint const& c) :
    CoefficientBase("Coefficient", fe, c),
    repr_(*this, *finite_element_, count_),
    str_((count_ < 10 ? "w_" + count_.str() : "w_{" + count_.str() + "}"))
{
}

//-----------------------------------------------------------------------------
Coefficient::Coefficient(repr_t const& repr) :
    CoefficientBase("Coefficient", repr),
    repr_(*this, *finite_element_, count_),
    str_((count_ < 10 ? "w_" + count_.str() : "w_{" + count_.str() + "}"))
{
}

//-----------------------------------------------------------------------------
Coefficient::~Coefficient()
{
}

//-----------------------------------------------------------------------------
std::vector<Class const*> const Coefficient::operands(
    std::string const& name) const
{
  std::vector<Class const*> expr0;
  if (name == this->name()) expr0.push_back(this);
  return expr0;
}

//-----------------------------------------------------------------------------
std::vector<std::vector<Class const*> > const Coefficient::level_operands(
    std::vector<std::vector<Class const*> > const& operands) const
{
  std::vector<std::vector<Class const*> > new_operands0;
  std::vector<Class const*> obj0;
  obj0.push_back(this);
  new_operands0.push_back(obj0);

  for (dolfin::uint i = 0; i < operands.size(); ++i)
    new_operands0.push_back(operands[i]);

  return new_operands0;
}

//-----------------------------------------------------------------------------
Coefficient const* Coefficient::create(Object::repr_t const& repr)
{
  return new Coefficient(repr);
}

//-----------------------------------------------------------------------------
Object::repr_t const& Coefficient::repr() const
{
  return repr_;
}

//-----------------------------------------------------------------------------
std::string const& Coefficient::str() const
{
  return str_;
}

//-----------------------------------------------------------------------------
void Coefficient::display() const
{
}

//-----------------------------------------------------------------------------
Constant::Constant(Cell const& cell, dolfin::uint const& c) :
    CoefficientBase("Constant", cell, c),
    repr_(*this, cell, count_),
    str_((count_ < 10 ? "c_" + count_.str() : "c_{" + count_.str() + "}"))
{
}

//-----------------------------------------------------------------------------
Constant::Constant(repr_t const& repr) :
    CoefficientBase("Constant", repr),
    repr_(*this, finite_element_->cell(), count_),
    str_((count_ < 10 ? "c_" + count_.str() : "c_{" + count_.str() + "}"))
{
}

//-----------------------------------------------------------------------------
Constant::~Constant()
{
}

//-----------------------------------------------------------------------------
std::vector<Class const*> const Constant::operands(
    std::string const& name) const
{
  std::vector<Class const*> expr0;
  if (name == this->name()) expr0.push_back(this);
  return expr0;
}

//-----------------------------------------------------------------------------
std::vector<std::vector<Class const*> > const Constant::level_operands(
    std::vector<std::vector<Class const*> > const& operands) const
{
  std::vector<std::vector<Class const*> > new_operands0;
  std::vector<Class const*> obj0;
  obj0.push_back(this);
  new_operands0.push_back(obj0);

  for (dolfin::uint i = 0; i < operands.size(); ++i)
    new_operands0.push_back(operands[i]);

  return new_operands0;
}

//-----------------------------------------------------------------------------
Constant const* Constant::create(Object::repr_t const& repr)
{
  return new Constant(repr);
}

//-----------------------------------------------------------------------------
Object::repr_t const& Constant::repr() const
{
  return repr_;
}

//-----------------------------------------------------------------------------
std::string const& Constant::str() const
{
  return str_;
}

//-----------------------------------------------------------------------------
void Constant::display() const
{
}

//-----------------------------------------------------------------------------
VectorConstant::VectorConstant(Cell const& cell, dolfin::uint const& dim,
                               dolfin::uint const& c) :
    CoefficientBase("VectorConstant", cell, dim, c),
//FIXME fill ValueArray    repr_(*this, cell, type<dolfin::uint>(finite_element_.value_shape()[0]), count_),
    repr_(*this, cell, type<dolfin::uint>(0), count_),
    str_((count_ < 10 ? "C_" + count_.str() : "C_{" + count_.str() + "}"))
{
}

//-----------------------------------------------------------------------------
VectorConstant::VectorConstant(repr_t const& repr) :
    CoefficientBase("VectorConstant", repr),
//FIXME fill ValueArray    repr_(*this, finite_element_.cell(), type<dolfin::uint>(finite_element_.value_shape()[0]), count_),
    repr_(*this, finite_element_->cell(), type<dolfin::uint>(0), count_),
    str_((count_ < 10 ? "C_" + count_.str() : "C_{" + count_.str() + "}"))
{
}

//-----------------------------------------------------------------------------
VectorConstant::~VectorConstant()
{
}

//-----------------------------------------------------------------------------
std::vector<Class const*> const VectorConstant::operands(
    std::string const& name) const
{
  std::vector<Class const*> expr0;
  if (name == this->name()) expr0.push_back(this);
  return expr0;
}

//-----------------------------------------------------------------------------
std::vector<std::vector<Class const*> > const VectorConstant::level_operands(
    std::vector<std::vector<Class const*> > const& operands) const
{
  std::vector<std::vector<Class const*> > new_operands0;
  std::vector<Class const*> obj0;
  obj0.push_back(this);
  new_operands0.push_back(obj0);

  for (dolfin::uint i = 0; i < operands.size(); ++i)
    new_operands0.push_back(operands[i]);

  return new_operands0;
}

//-----------------------------------------------------------------------------
VectorConstant const* VectorConstant::create(Object::repr_t const& repr)
{
  return new VectorConstant(repr);
}

//-----------------------------------------------------------------------------
Object::repr_t const& VectorConstant::repr() const
{
  return repr_;
}

//-----------------------------------------------------------------------------
std::string const& VectorConstant::str() const
{
  return str_;
}

//-----------------------------------------------------------------------------
void VectorConstant::display() const
{
}

//-----------------------------------------------------------------------------
TensorConstant::TensorConstant(
    Cell const& cell, ValueArray const& shape,
    dolfin::_ordered_map<dolfin::uint, dolfin::uint> const& symmetry, dolfin::uint const& c) :
    CoefficientBase("TensorConstant", cell, shape, symmetry, c),
    repr_(*this, finite_element_->cell(), /*finite_element_.value_shape(),
     finite_element_.symmetry(),*/
          count_),
    str_((count_ < 10 ? "C_" + count_.str() : "C_{" + count_.str() + "}"))
{
}

//-----------------------------------------------------------------------------
TensorConstant::TensorConstant(repr_t const& repr) :
    CoefficientBase("TensorConstant", repr),
    repr_(*this, finite_element_->cell(), /*finite_element_.value_shape(),
     finite_element_.symmetry(),*/
          count_),
    str_((count_ < 10 ? "C_" + count_.str() : "C_{" + count_.str() + "}"))
{
}

//-----------------------------------------------------------------------------
TensorConstant::~TensorConstant()
{
}

//-----------------------------------------------------------------------------
std::vector<Class const*> const TensorConstant::operands(
    std::string const& name) const
{
  std::vector<Class const*> expr0;
  if (name == this->name()) expr0.push_back(this);
  return expr0;
}

//-----------------------------------------------------------------------------
std::vector<std::vector<Class const*> > const TensorConstant::level_operands(
    std::vector<std::vector<Class const*> > const& operands) const
{
  std::vector<std::vector<Class const*> > new_operands0;
  std::vector<Class const*> obj0;
  obj0.push_back(this);
  new_operands0.push_back(obj0);

  for (dolfin::uint i = 0; i < operands.size(); ++i)
    new_operands0.push_back(operands[i]);

  return new_operands0;
}

//-----------------------------------------------------------------------------
TensorConstant const* TensorConstant::create(Object::repr_t const& repr)
{
  return new TensorConstant(repr);
}

//-----------------------------------------------------------------------------
Object::repr_t const& TensorConstant::repr() const
{
  return repr_;
}

//-----------------------------------------------------------------------------
std::string const& TensorConstant::str() const
{
  return str_;
}

//-----------------------------------------------------------------------------
void TensorConstant::display() const
{
}

}
