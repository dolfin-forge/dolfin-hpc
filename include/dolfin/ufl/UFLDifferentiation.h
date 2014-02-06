// Copyright (C) 2014 Bärbel Janssen.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  
// Last changed: 

#ifndef __UFL_DIFFERENTIATION_H_
#define __UFL_DIFFERENTIATION_H_

//#include <string>
//#include <vector>

#include <dolfin/ufl/UFLClass.h>
#include <dolfin/ufl/UFLExpression.h>
#include <dolfin/ufl/UFLIndex.h>

//#include <dolfin/common/types.h>

namespace ufl
{

  /**
   *  DOCUMENTATION:
   *
   *  @class  Class
   *
   *  @brief  Provides an interface complying with UFL Differentiation.
   */

  class SpatialDerivative : public Class
  {

    public:

      ///
      SpatialDerivative(Expression const& expression, IndexBase const& index);

      ///
      ~SpatialDerivative();

      //--- INTERFACE -------------------------------------------------------------

      ///
      std::pair<Expression const, IndexBase const> const& operands() const;

//      ///
//      free_indices() const;

//      ///
//      index_dimensions() const;

//      ///
//      shape() const;

      /// UFL: Return whether this expression is spatially constant over each cell
      bool const is_cellwise_constant() const;

      //--- INTERFACE inherited from UFLClass -------------------------------------
      
      /// __repr__
      repr_t const repr() const;

      /// __str__
      std::string const str() const;

      ///
      void display() const;

      ///
      SpatialDerivative const* create(repr_t const & repr) const;

    private:

      Expression const expression_;
      IndexBase const& index_;

      mutable repr_t repr_;
      mutable std::string str_;

  };
} /* namespace ufl */
#endif /* __UFL_DIFFERENTIATION_H_ */
