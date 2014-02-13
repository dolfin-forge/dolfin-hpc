// Copyright (C) 2014 Bärbel Janssen.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:
// Last changed:

#ifndef __UFL_LIST_H
#define __UFL_LIST_H

#include <dolfin/ufl/UFLIntegral.h>

namespace ufl
{

/**
 *  DOCUMENTATION:
 *
 *  @class  List
 *
 *  @brief  Provides an interface complying with Lists.
 */

  class List : public Class
  {
    public:

      ///
      List(std::vector<Integral> const & integrals);

      ///
      List (repr_t const & repr);

      ///
      ~List();

      //--- INTERFACE inherited from UFLClass -------------------------------------

      /// __repr__
      repr_t const repr() const;

      /// __str__
      std::string const str() const;

      ///
      void display() const;

    private:

      std::vector<Integral const *> integrals_;

      repr_t const repr_;
      std::string const str_;
  };

} /* namespace ufl */
#endif /* __UFL_LIST_H */
