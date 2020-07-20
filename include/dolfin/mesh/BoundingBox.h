// Copyright (C) 2015 Aurelien Larcher.
// Licensed under the GNU LGPL Version 2.1.

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
  /// Create unit bounding box in R^d
  BoundingBox( uint d = dolfin::Space::MAX_DIMENSION )
    : D_( d )
    , BOX_( new Point[D_ + 1] )
  {
    for ( uint i = 1; i <= D_; ++i )
    {
      BOX_[i][i - 1] = 1.0;
    }
  };

  /// Copy constructor
  BoundingBox( BoundingBox const & other )
    : D_( other.D_ )
    , BOX_( new Point[D_ + 1] )
  {
    for ( uint i = 1; i <= D_; ++i )
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
  Point & operator[]( uint i );

  /// Access i-the point
  Point const & operator[]( uint i ) const;

  /// Assignment
  BoundingBox & operator=( BoundingBox const & other );

  /// Translation
  BoundingBox & operator+=( Point const & p );

  /// Translation
  BoundingBox & operator-=( Point const & p );

  /// Homothety
  BoundingBox & operator*=( real const a );

  /// Scaling along axis
  BoundingBox & scale( uint axis, real const a );

  /// Dilatation
  BoundingBox & operator*=( Point const & p );

  /// Return bounding box centroid
  Point centroid() const;

  // Display information
  void disp() const;

private:
  uint const    D_;
  Point * const BOX_;
};

//-----------------------------------------------------------------------------
inline Point & BoundingBox::operator[]( uint i )
{
  dolfin_assert( i <= D_ );
  return BOX_[i];
}

//-----------------------------------------------------------------------------
inline Point const & BoundingBox::operator[]( uint i ) const
{
  dolfin_assert( i <= D_ );
  return BOX_[i];
}

//-----------------------------------------------------------------------------
inline BoundingBox & BoundingBox::operator=( BoundingBox const & other )
{
  if ( this != &other )
  {
    if ( this->D_ != other.D_ )
    {
      error( "BoundingBox : trying to assign with different dimensions" );
    }
    for ( uint i = 0; i <= D_; ++i )
    {
      BOX_[i] = other.BOX_[i];
    }
  }
  return *this;
}

//-----------------------------------------------------------------------------
inline BoundingBox & BoundingBox::operator+=( Point const & p )
{
  for ( uint i = 0; i <= D_; ++i )
  {
    BOX_[i] += p;
  }
  return *this;
}

//-----------------------------------------------------------------------------
inline BoundingBox & BoundingBox::operator-=( Point const & p )
{
  for ( uint i = 0; i <= D_; ++i )
  {
    BOX_[i] -= p;
  }
  return *this;
}

//-----------------------------------------------------------------------------
inline BoundingBox & BoundingBox::operator*=( real const a )
{
  Point c = this->centroid();
  for ( uint i = 0; i <= D_; ++i )
  {
    BOX_[i] = c + a * ( BOX_[i] - c );
  }
  return *this;
}

//-----------------------------------------------------------------------------
inline BoundingBox & BoundingBox::scale( uint axis, real const a )
{
  for ( uint i = 0; i <= D_; ++i )
  {
    BOX_[i][axis] *= a;
  }
  return *this;
}

//-----------------------------------------------------------------------------
inline BoundingBox & BoundingBox::operator*=( Point const & p )
{
  Point c = this->centroid();
  for ( uint i = 0; i <= D_; ++i )
  {
    for ( uint d = 0; d < Point::MAX_SIZE; ++d )
    {
      BOX_[i][d] = c[d] + p[d] * ( BOX_[i][d] - c[d] );
    }
  }
  return *this;
}

//-----------------------------------------------------------------------------
inline Point BoundingBox::centroid() const
{
  Point c( BOX_[0] );
  for ( uint i = 1; i <= D_; ++i )
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
  for ( uint i = 0; i <= D_; ++i )
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
