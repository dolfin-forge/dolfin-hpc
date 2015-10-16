// Copyright (C) 2015 Aurelien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//

#ifndef __DOLFIN_BOUNDING_BOX_H
#define __DOLFIN_BOUNDING_BOX_H

#include <dolfin/common/types.h>
#include <dolfin/mesh/Point.h>

#include <cmath>

namespace dolfin
{

/**
 *  @class  BoundingBox
 *
 *  @brief  Implements a bounding box.
 *
 */

class BoundingBox
{

public:

  /// Create unit bounding box
  BoundingBox(uint d) :
    D_(d),
    BOX_(new Point[d+1])
  {
    for (uint i = 1; i <= D_; ++i)
    {
      BOX_[i][i - 1] = 1.0;
    }
  };

  /// Destructor
  ~BoundingBox()
  {
    delete [] BOX_;
  };

  /// Access i-the point
  Point& operator[](uint i)
  {
    dolfin_assert(i <= D_);
    return BOX_[i];
  }

  /// Access i-the point
  Point const& operator[](uint i) const
  {
    dolfin_assert(i <= D_);
    return BOX_[i];
  }

  /// Assignment
  BoundingBox& operator=(BoundingBox const& other)
  {
    if(this != &other)
    {
      for (uint i = 0; i <= D_; ++i)
      {
        BOX_[i] = other.BOX_[i];
      }
    }
    return *this;
  }

  /// Translation
  BoundingBox& operator+=(Point const& p)
  {
    for (uint i = 0; i <= D_; ++i)
    {
      BOX_[i] += p;
    }
    return *this;
  }

  /// Translation
  BoundingBox& operator-=(Point const& p)
  {
    for (uint i = 0; i <= D_; ++i)
    {
      BOX_[i] -= p;
    }
    return *this;
  }

  /// Homothety
  BoundingBox& operator*=(real const a)
  {
    Point c = this->centroid();
    for (uint i = 0; i <= D_; ++i)
    {
      BOX_[i] = c + a * (BOX_[i] - c);
    }
    return *this;
  }

  /// Dilatation
  BoundingBox& operator*=(Point const& p)
  {
    Point c = this->centroid();
    for (uint i = 0; i <= D_; ++i)
    {
      for (uint d = 0; d < Point::MAX_SIZE; ++d)
      {
        BOX_[i][d] = c[d] + p[d] * (BOX_[i][d] - c[d]);
      }
    }
    return *this;
  }

  /// Access i-the point
  Point centroid() const
  {
    Point c(BOX_[0]);
    for (uint i = 1; i <= D_; ++i)
    {
      c += 0.5 * (BOX_[i] - BOX_[0]);
    }
    return c;
  }

  // Display information
  void disp() const
  {
    message("BoundingBox");
    begin(  "-----------");
    message("Dimension : %d", D_);
    begin(  "Points    :");
    for (uint i = 0; i <= D_; ++i)
    {
      cout << BOX_[i] << endl;
    }
    end();
    begin(  "Centroid  :");
    cout << this->centroid() << endl;
    end();
    end();
    skip();
  }

private:

  uint const D_;
  Point * const BOX_;

};

} /* namespace dolfin */

#endif /* __DOLFIN_BOUNDING_BOX_H */
