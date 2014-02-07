// Copyright (C) 2014 Bärbel Janssen.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  
// Last changed: 

#ifndef __UFL_ARGUMENT_H
#define __UFL_ARGUMENT_H

#include <dolfin/ufl/UFLClass.h>
#include <dolfin/ufl/UFLFiniteElement.h>

namespace ufl
{

/**
 *  DOCUMENTATION:
 *
 *  @class  Argument
 *
 *  @brief  Provides an interface complying with UFL Argument.
 */

  class Argument : public Class
  {
    public:

      ///
      Argument(FiniteElementBase const& finite_element,
          dolfin::uint const& count);

      ///
      ~Argument();
  
      //--- INTERFACE -------------------------------------------------------------

      /// Return a reference to the FiniteElementBase of this Argument
      FiniteElementBase const& element() const;

      /// Return a reference to the value shape of the FiniteElementBase of this Argument
      ValueArray const& shape() const;

      /// Return whether the basis functions of this element is spatially constant
      /// over each cell
      bool const is_cellwise_constant() const;

      /// __repr__
      repr_t const repr() const;

      /// __str__
      std::string const str() const;

      ///
      void display() const;

      ///
      Argument const* create(repr_t const & repr) const;

    protected:

      FiniteElementBase const& finite_element_;

      mutable repr_t repr_;
      mutable std::string str_;

      dolfin::uint const count_;
  };

} /* namespace ufl */
#endif /* __UFL_ARGUMENT_H */
