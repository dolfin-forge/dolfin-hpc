// Copyright (C) 2014 Bärbel Janssen.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  
// Last changed: 

#ifndef __UFL_COEFFICIENT_H
#define __UFL_COEFFICIENT_H

#include <dolfin/ufl/UFLClass.h>
#include <dolfin/ufl/UFLFiniteElement.h>

namespace ufl
{

/**
 *  DOCUMENTATION:
 *
 *  @class  CoefficientBase
 *
 *  @brief  Base class for interface complying with UFL Coefficient.
 */


  class CoefficientBase : public Class
  {
    public:
  
      //--- INTERFACE -------------------------------------------------------------

      /// Return a reference to the FiniteElementBase of this Coefficient
      FiniteElementBase const& element() const;

      /// Return a reference to the value shape of the FiniteElementBase of this Coefficient
      ValueArray const& shape() const;

      /// Return a reference to the cell of the FiniteElementBase of this Coefficient
      Cell const& cell() const;

      /// Return whether the basis functions of this element is spatially constant
      /// over each cell
      bool const is_cellwise_constant() const;

      /// __repr__
      repr_t const repr() const;

      /// __str__
      std::string const str() const;

      ///
      void display() const;

    protected:

      ///
      CoefficientBase(std::string const& name,
          FiniteElementBase const& finite_element, dolfin::uint const& count);

      ///
      CoefficientBase(std::string const& name, repr_t const& repr);

      ///
      ~CoefficientBase();

      FiniteElementBase const& finite_element_;
      type<dolfin::uint> const count_;

      repr_t const repr_;
      std::string const str_;

  };

/**
 *  DOCUMENTATION:
 *
 *  @class  Coefficient
 *
 *  @brief  Provides an interface complying with UFL Coefficient.
 */

  class Coefficient : public CoefficientBase
  {
    public:

      ///
      Coefficient(FiniteElementBase const& finite_element,
          dolfin::uint const& count);

      ///
      Coefficient(repr_t const& repr);

      ///
      ~Coefficient();
  
      //--- INTERFACE -------------------------------------------------------------

      /// __repr__
      repr_t const repr() const;

      /// __str__
      std::string const str() const;

      ///
      void display() const;

    protected:

      repr_t const repr_;
      std::string const str_;

  };


/**
 *  DOCUMENTATION:
 *
 *  @class  Constant
 *
 *  @brief  Provides an interface complying with UFL Constant.
 */

  class Constant : public CoefficientBase
  {
    public:

      ///
      Constant(Cell const& cell,
          dolfin::uint const& count);

      ///
      Constant(repr_t const& repr);

      ///
      ~Constant();
  
      /// __repr__
      repr_t const repr() const;

      /// __str__
      std::string const str() const;

      ///
      void display() const;

    protected:

      repr_t const repr_;
      std::string const str_;

  };


/**
 *  DOCUMENTATION:
 *
 *  @class  VectorConstant
 *
 *  @brief  Provides an interface complying with UFL VectorConstant.
 */

  class VectorConstant : public CoefficientBase
  {
    public:

      ///
      VectorConstant(Cell const& cell, dolfin::uint const& dim,
          dolfin::uint const& count);

      ///
      VectorConstant(repr_t const& repr);

      ///
      ~VectorConstant();
  
      /// __repr__
      repr_t const repr() const;

      /// __str__
      std::string const str() const;

      ///
      void display() const;

    protected:

      repr_t const repr_;
      std::string const str_;

  };


/**
 *  DOCUMENTATION:
 *
 *  @class  TensorConstant
 *
 *  @brief  Provides an interface complying with UFL TensorConstant.
 */

  class TensorConstant : public CoefficientBase
  {
    public:

      ///
      TensorConstant(Cell const& cell, ValueArray const& shape, 
          std::map<dolfin::uint, dolfin::uint> const& symmetry, dolfin::uint const& count);

      ///
      TensorConstant(repr_t const& repr);

      ///
      ~TensorConstant();
  
      /// __repr__
      repr_t const repr() const;

      /// __str__
      std::string const str() const;

      ///
      void display() const;

    protected:

      repr_t const repr_;
      std::string const str_;

  };

} /* namespace ufl */
#endif /* __UFL_COEFFICIENT_H */
