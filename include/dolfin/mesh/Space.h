// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2014-11-05
// Last changed: 2014-11-05

#ifndef __DOLFIN_MESH_SPACE_H
#define __DOLFIN_MESH_SPACE_H

#include <dolfin/common/Clonable.h>

#include <dolfin/log/log.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
struct Space
{
  static uint const MAX_DIMENSION = 3;

  /// Constructor
  Space() {}

  /// Destructor
  virtual ~Space() {}

  /// Space dimension
  virtual uint dim() const = 0;

  /// Clone pattern
  virtual Space* clone() const = 0;

  /// Display info
  virtual void disp() const = 0;

};

//-----------------------------------------------------------------------------

class EuclideanSpace : public Space, public Clonable<EuclideanSpace>
{

public:

  EuclideanSpace(uint dim) :
    dim_(dim)
  {
  }

  ~EuclideanSpace()
  {
  }

  inline uint dim() const
  {
    return dim_;
  }

  inline EuclideanSpace * clone() const
  {
    return Clonable<EuclideanSpace>::clone();
  }

  void disp() const
  {
    section("EuclideanSpace");
    message("dimension : %u", dim_);
    endblock();
  }

private:

  uint dim_;

};

//-----------------------------------------------------------------------------

} /* namespace dolfin */

#endif /* __DOLFIN_MESH_SPACE_H */
