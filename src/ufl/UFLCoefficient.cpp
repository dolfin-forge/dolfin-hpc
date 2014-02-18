// Copyright (C) 2014 Bärbel Janssen.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  
// Last changed: 

#include <dolfin/ufl/UFLCoefficient.h>
#include <dolfin/ufl/UFLVectorElement.h>
#include <dolfin/ufl/UFLTensorElement.h>

namespace ufl
{

//-----------------------------------------------------------------------------
  CoefficientBase::CoefficientBase(std::string const& name,
      FiniteElementBase const& fe, dolfin::uint const& c) :
    Class(name),
    finite_element_(fe),
    count_(c)
  {
  }

//-----------------------------------------------------------------------------
  CoefficientBase::CoefficientBase(std::string const& name, repr_t const& repr) :
    Class(name),
    finite_element_(*FiniteElementBase::create(arg(0))),
    count_(arg(1))
  {
  }

//-----------------------------------------------------------------------------
  CoefficientBase::~CoefficientBase()
  {
  }
  
//-----------------------------------------------------------------------------
  FiniteElementBase const& CoefficientBase::element() const
  {
    return finite_element_;  
  }

//-----------------------------------------------------------------------------
  ValueArray const& CoefficientBase::shape() const
  {
    return finite_element_.value_shape();  
  }

//-----------------------------------------------------------------------------
  Cell const& CoefficientBase::cell() const
  {
    return finite_element_.cell();  
  }

//-----------------------------------------------------------------------------
  bool const CoefficientBase::is_cellwise_constant() const
  {
    return false;
  }
 
//-----------------------------------------------------------------------------
  Object::repr_t const CoefficientBase::repr() const
  {
    return repr_;
  }

//-----------------------------------------------------------------------------
  std::string const CoefficientBase::str() const
  {
    return str_;
  }
 
//-----------------------------------------------------------------------------
  void CoefficientBase::display() const
  {
  }

//-----------------------------------------------------------------------------
  Coefficient::Coefficient(FiniteElementBase const& fe, dolfin::uint const& c) :
    CoefficientBase("Coefficient", fe, c),
    repr_(*this, finite_element_, count_),
    str_((count_ < 10 ? "w_" + count_.str() : "w_{" + count_.str() + "}"))
  {
  }

//-----------------------------------------------------------------------------
  Coefficient::Coefficient(repr_t const& repr) :
    CoefficientBase("Coefficient", repr),//*FiniteElementBase::create(arg(0)), arg(1)),
    repr_(*this, finite_element_, count_),
    str_((count_ < 10 ? "w_" + count_.str() : "w_{" + count_.str() + "}"))
  {
  }

//-----------------------------------------------------------------------------
  Coefficient::~Coefficient()
  {
  }
  
//-----------------------------------------------------------------------------
  Object::repr_t const Coefficient::repr() const
  {
    return repr_;
  }

//-----------------------------------------------------------------------------
  std::string const Coefficient::str() const
  {
    return str_;
  }
 
//-----------------------------------------------------------------------------
  void Coefficient::display() const
  {
  }

//-----------------------------------------------------------------------------
  Constant::Constant(Cell const& cell, dolfin::uint const& c) :
    CoefficientBase("Constant", FiniteElement(Family::R, cell, 0), c),
    repr_(*this, finite_element_.cell(), count_),
    str_((count_ < 10 ? "c_" + count_.str() : "c_{" + count_.str() + "}"))
  {
  }

//-----------------------------------------------------------------------------
  Constant::Constant(repr_t const& repr) :
    CoefficientBase("Constant", FiniteElement(Family::R, Cell(arg(0)), 0), type<dolfin::uint>(arg(1))),
    repr_(*this, finite_element_.cell(), count_),
    str_((count_ < 10 ? "c_" + count_.str() : "c_{" + count_.str() + "}"))
  {
  }

//-----------------------------------------------------------------------------
  Constant::~Constant()
  {
  }
  
//-----------------------------------------------------------------------------
  Object::repr_t const Constant::repr() const
  {
    return repr_;
  }

//-----------------------------------------------------------------------------
  std::string const Constant::str() const
  {
    return str_;
  }
 
//-----------------------------------------------------------------------------
  void Constant::display() const
  {
  }

//-----------------------------------------------------------------------------
  VectorConstant::VectorConstant(Cell const& cell, 
      dolfin::uint const& dim, dolfin::uint const& c) :
    CoefficientBase("VectorConstant", VectorElement(Family::R, cell, 0, dim), c),
    repr_(*this, finite_element_.cell(), type<dolfin::uint>(finite_element_.value_shape()[0]), count_),
    str_((count_ < 10 ? "C_" + count_.str() : "C_{" + count_.str() + "}"))
  {
  }

//-----------------------------------------------------------------------------
  VectorConstant::VectorConstant(repr_t const& repr) :
    CoefficientBase("VectorConstant", VectorElement(Family::R, Cell(arg(0)),
          0, type<dolfin::uint>(arg(1))), type<dolfin::uint>(arg(2))),
    repr_(*this, finite_element_.cell(), type<dolfin::uint>(finite_element_.value_shape()[0]), count_),
    str_((count_ < 10 ? "C_" + count_.str() : "C_{" + count_.str() + "}"))
  {
  }

//-----------------------------------------------------------------------------
  VectorConstant::~VectorConstant()
  {
  }
  
//-----------------------------------------------------------------------------
  Object::repr_t const VectorConstant::repr() const
  {
    return repr_;
  }

//-----------------------------------------------------------------------------
  std::string const VectorConstant::str() const
  {
    return str_;
  }
 
//-----------------------------------------------------------------------------
  void VectorConstant::display() const
  {
  }

//-----------------------------------------------------------------------------
  TensorConstant::TensorConstant(Cell const& cell, ValueArray const& shape,
      std::map<dolfin::uint, dolfin::uint> const& symmetry, dolfin::uint const& c) :
    CoefficientBase("TensorConstant", TensorElement(Family::R, cell, 0, cell.geometric_dimension()), c),
    repr_(*this, finite_element_.cell(), /*finite_element_.value_shape(), 
        finite_element_.symmetry(),*/ count_),
    str_((count_ < 10 ? "C_" + count_.str() : "C_{" + count_.str() + "}"))
  {
  }

//-----------------------------------------------------------------------------
  TensorConstant::TensorConstant(repr_t const& repr) :
    CoefficientBase("TensorConstant", 
        TensorElement(Family::R, Cell(arg(0)), 0, Cell(arg(0)).geometric_dimension()),
        type<dolfin::uint>(arg(3))),
    repr_(*this, finite_element_.cell(), /*finite_element_.value_shape(), 
        finite_element_.symmetry(),*/ count_),
    str_((count_ < 10 ? "C_" + count_.str() : "C_{" + count_.str() + "}"))
  {
  }

//-----------------------------------------------------------------------------
  TensorConstant::~TensorConstant()
  {
  }
  
//-----------------------------------------------------------------------------
  Object::repr_t const TensorConstant::repr() const
  {
    return repr_;
  }

//-----------------------------------------------------------------------------
  std::string const TensorConstant::str() const
  {
    return str_;
  }
 
//-----------------------------------------------------------------------------
  void TensorConstant::display() const
  {
  }
}
