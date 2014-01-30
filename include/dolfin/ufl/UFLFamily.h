// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2014-01-28
// Last changed: 2014-01-28

#ifndef __UFL_FAMILY_H_
#define __UFL_FAMILY_H_

#include <dolfin/ufl/UFLtype.h>

#include <map>
#include <set>
#include <string>

namespace ufl
{

/**
 *  DOCUMENTATION:
 *
 *  @class  UFLFamily
 *
 *  @brief  Provides a C++ equivalent to family types.
 */

class Family : public type<std::string>
{

public:

  enum Type
  {
    ARG,
    AW,
    BDFM,
    BDM,
    CR,
    DG,
    HER,
    CG,
    MTW,
    MOR,
    N1curl,
    N2curl,
    RT,
    BQ,
    B,
    Q,
    R,
    U,
    Mixed,
    Vector,
    Tensor,
    Enriched,
    Restricted
  };

  ///
  Family(Family::Type const& t);

  ///
  ~Family();

  Family::Type const type() const;

private:

  Type const type_;

};

} /* namespace ufl */
#endif /* __UFL_FAMILY_H_ */
