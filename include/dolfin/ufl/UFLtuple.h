// Copyright (C) 2014 Bärbel Janssen.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:
// Last changed:

#ifndef __UFL_TUPLE_H
#define __UFL_TUPLE_H

#include <dolfin/ufl/UFLExpression.h>

namespace ufl
{

/**
 *  DOCUMENTATION:
 *
 *  @class  Tuple
 *
 *  @brief  Provides an interface complying with Tuple.
 */

  class Tuple : public Class
  {
    public:

      ///
      Tuple(std::vector<Expression const *>& args);

      ///
      Tuple (repr_t const & repr);

      ///
      ~Tuple();

      //--- INTERFACE inherited from UFLClass -------------------------------------

      /// __repr__
      repr_t const repr() const;

      /// __str__
      std::string const str() const;

      ///
      void display() const;

    private:

      std::vector<Expression const *> expressions_;

      repr_t const repr_;
      std::string const str_;
  };

} /* namespace ufl */
#endif /* __UFL_LIST_H */
