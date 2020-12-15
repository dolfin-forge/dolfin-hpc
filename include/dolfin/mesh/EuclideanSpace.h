// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_MESH_EUCLIDEAN_SPACE_H
#define __DOLFIN_MESH_EUCLIDEAN_SPACE_H

#include <dolfin/common/Clonable.h>
#include <dolfin/log/log.h>
#include <dolfin/mesh/Space.h>

namespace dolfin
{

class EuclideanSpace : public Space, public Clonable<EuclideanSpace>
{

public:

  EuclideanSpace(uint dim) :
    dim_(dim)
  {
  }

  ~EuclideanSpace() override = default;

  /// Equality
  bool operator==(EuclideanSpace const& other) const
  { return dim_ == other.dim_; }

  /// Non-equality
  bool operator!=(EuclideanSpace const& other) const
  { return dim_ != other.dim_; }

  /// Return dimension
  inline uint dim() const override
  {
    return dim_;
  }

  /// Return clone
  inline EuclideanSpace * clone() const override
  {
    return Clonable<EuclideanSpace>::clone();
  }

  void disp() const override
  {
    section("EuclideanSpace");
    message("dimension : %u", dim_);
    end();
  }

private:

  uint dim_;

};

//-----------------------------------------------------------------------------

} /* namespace dolfin */

#endif /* __DOLFIN_MESH_EUCLIDEAN_SPACE_H */
