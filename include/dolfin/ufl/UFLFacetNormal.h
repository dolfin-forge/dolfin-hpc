// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2014-01-21
// Last changed: 2014-01-21

#ifndef __UFL_FACET_NORMAL_H_
#define __UFL_FACET_NORMAL_H_

#include <dolfin/ufl/UFLGeometricQuantity.h>

namespace dolfin
{

/**
 *  DOCUMENTATION:
 *
 *  @class
 *
 *  @brief
 */

class UFLCell;

class UFLFacetNormal : public UFLGeometricQuantity
{

public:

  ///
  UFLFacetNormal(UFLCell const& cell);

  ///
  ~UFLFacetNormal();

  ///
  std::vector<uint> const& shape() const;

  /// __repr__
  std::string const repr() const;

  /// __str__
  std::string const str() const;

private:

  std::vector<uint> const shape_;

  std::string const repr_;
  std::string const str_;

};

} /* namespace dolfin */
#endif /* __UFL_FACET_NORMAL_H_ */
