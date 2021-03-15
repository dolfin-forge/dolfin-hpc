// Copyright (C) 2015 Aurelien Larcher.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_BOUNDING_BOX_H
#define __DOLFIN_BOUNDING_BOX_H

#include <dolfin/common/types.h>
#include <dolfin/mesh/Point.h>

#include <cmath>

namespace dolfin
{

//-----------------------------------------------------------------------------

/**
 *  @class  BoundingBox
 *
 *  @brief  Implements a bounding box.
 *
 */

class BoundingBox
{

public:
  /// Create unit bounding box in R^d
  BoundingBox( size_t d = dolfin::Space::MAX_DIMENSION )
    : D_( d )
    , BOX_( new Point[D_ + 1] )
  {
    for ( size_t i = 1; i <= D_; ++i )
    {
      BOX_[i][i - 1] = 1.0;
    }
  };

  /// Copy constructor
  BoundingBox( BoundingBox const & other )
    : D_( other.D_ )
    , BOX_( new Point[D_ + 1] )
  {
    for ( size_t i = 1; i <= D_; ++i )
    {
      BOX_[i][i - 1] = 1.0;
    }
  };

  /// Destructor
  ~BoundingBox()
  {
    delete[] BOX_;
  };

  /// Access i-the point
  auto operator[]( size_t i ) -> Point &;

  /// Access i-the point
  auto operator[]( size_t i ) const -> Point const &;

  /// Assignment
  auto operator=( BoundingBox const & other ) -> BoundingBox &;

  /// Translation
  auto operator+=( Point const & p ) -> BoundingBox &;

  /// Translation
  auto operator-=( Point const & p ) -> BoundingBox &;

  /// Homothety
  auto operator*=( real const a ) -> BoundingBox &;

  /// Scaling along axis
  auto scale( size_t axis, real const a ) -> BoundingBox &;

  /// Dilatation
  auto operator*=( Point const & p ) -> BoundingBox &;

  /// Return bounding box centroid
  auto centroid() const -> Point;

  // Display information
  void disp() const;

private:
  size_t const    D_;
  Point * const BOX_;
};

//-----------------------------------------------------------------------------

inline auto BoundingBox::operator[]( size_t i ) -> Point &
{
  dolfin_assert( i <= D_ );
  return BOX_[i];
}

//-----------------------------------------------------------------------------

inline auto BoundingBox::operator[]( size_t i ) const -> Point const &
{
  dolfin_assert( i <= D_ );
  return BOX_[i];
}

//-----------------------------------------------------------------------------

inline auto BoundingBox::operator=( BoundingBox const & other ) -> BoundingBox &
{
  if ( this != &other )
  {
    if ( this->D_ != other.D_ )
    {
      error( "BoundingBox : trying to assign with different dimensions" );
    }
    for ( size_t i = 0; i <= D_; ++i )
    {
      BOX_[i] = other.BOX_[i];
    }
  }
  return *this;
}

//-----------------------------------------------------------------------------

inline auto BoundingBox::operator+=( Point const & p ) -> BoundingBox &
{
  for ( size_t i = 0; i <= D_; ++i )
  {
    BOX_[i] += p;
  }
  return *this;
}

//-----------------------------------------------------------------------------

inline auto BoundingBox::operator-=( Point const & p ) -> BoundingBox &
{
  for ( size_t i = 0; i <= D_; ++i )
  {
    BOX_[i] -= p;
  }
  return *this;
}

//-----------------------------------------------------------------------------

inline auto BoundingBox::operator*=( real const a ) -> BoundingBox &
{
  Point c = this->centroid();
  for ( size_t i = 0; i <= D_; ++i )
  {
    BOX_[i] = c + a * ( BOX_[i] - c );
  }
  return *this;
}

//-----------------------------------------------------------------------------

inline auto BoundingBox::scale( size_t axis, real const a ) -> BoundingBox &
{
  for ( size_t i = 0; i <= D_; ++i )
  {
    BOX_[i][axis] *= a;
  }
  return *this;
}

//-----------------------------------------------------------------------------

inline auto BoundingBox::operator*=( Point const & p ) -> BoundingBox &
{
  Point c = this->centroid();
  for ( size_t i = 0; i <= D_; ++i )
  {
    for ( size_t d = 0; d < Point::MAX_SIZE; ++d )
    {
      BOX_[i][d] = c[d] + p[d] * ( BOX_[i][d] - c[d] );
    }
  }
  return *this;
}

//-----------------------------------------------------------------------------

inline auto BoundingBox::centroid() const -> Point
{
  Point c( BOX_[0] );
  for ( size_t i = 1; i <= D_; ++i )
  {
    c += 0.5 * ( BOX_[i] - BOX_[0] );
  }
  return c;
}

//-----------------------------------------------------------------------------

inline void BoundingBox::disp() const
{
  section( "BoundingBox" );
  message( "Dimension : %d", D_ );
  section( "Points" );
  for ( size_t i = 0; i <= D_; ++i )
  {
    cout << BOX_[i] << "\n";
  }
  end();
  section( "Centroid" );
  cout << this->centroid() << "\n";
  end();
  end();
}

//-----------------------------------------------------------------------------

} /* namespace dolfin */

#endif /* __DOLFIN_BOUNDING_BOX_H */
