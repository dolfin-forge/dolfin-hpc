// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2014-01-21
// Last changed: 2014-01-21

#ifndef __UFL_CELL_H_
#define __UFL_CELL_H_

#include <dolfin/ufl/UFLClass.h>
#include <dolfin/ufl/UFLDomain.h>
#include <dolfin/ufl/UFLSpace.h>

namespace dolfin
{

/**
 *  DOCUMENTATION:
 *
 *  @class  UFLCell
 *
 *  @brief  Provides an interface complying with UFL Cell.
 */

class UFLCell : public UFLClass
{

  ///
  UFLCell(UFLDomain::Type const& domain);

  ///
  UFLCell(UFLDomain::Type const& domain, UFLSpace const& space );

  ///
  ~UFLCell();

  /// __repr__
  std::string const repr() const;

  /// __str__
  std::string const str() const;

private:

  UFLDomain::Type const domain_;
  UFLSpace const space_;
  bool const invalid_;
  uint const geometric_dimension_;
  uint const topological_dimension_;

  std::string const repr_;
  std::string const str_;

};

} /* namespace dolfin */
#endif /* __UFL_CELL_H */
