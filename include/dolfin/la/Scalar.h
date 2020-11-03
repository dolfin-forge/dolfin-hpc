// Copyright (C) 2007-2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_SCALAR_H
#define __DOLFIN_SCALAR_H

#include "GenericTensor.h"
#include <dolfin/common/maybe_unused.h>
#include <dolfin/config/dolfin_config.h>
#include <dolfin/main/MPI.h>
#include <dolfin/parameter/parameters.h>

#ifdef HAVE_PETSC
#include "PETScFactory.h"
#endif
#ifdef HAVE_JANPACK
#include "JANPACKFactory.h"
#endif

namespace dolfin
{

class GenericSparsityPattern;

/// This class represents a real-valued scalar quantity and
/// implements the GenericTensor interface for scalars.

class Scalar : public GenericTensor
{
public:
  /// Create zero scalar
  Scalar()
    : GenericTensor()
  {
  }

  /// Destructor
  ~Scalar() override = default;

  //--- Implementation of the GenericTensor interface ---

  /// Initialize zero tensor using sparsity pattern
  void init( const GenericSparsityPattern & sparsity_pattern ) override;

  /// Return copy of tensor
  Scalar * copy() const override;

  /// Return tensor rank (number of dimensions)
  uint rank() const override;

  /// Return size of given dimension
  uint size( uint ) const override;

  /// Get block of values
  void get( real * block, const uint *, const uint * const * ) const override;

  /// Set block of values
  void set( const real * block, const uint *, const uint * const * ) override;

  /// Add block of values
  void add( const real * block, const uint *, const uint * const * ) override;

  /// Set all entries to zero and keep any sparse structure
  void zero() override;

  /// Finalize assembly of tensor
  void apply( FinalizeType ) override;

  /// Display tensor
  void disp( uint ) const override;

  //--- Scalar interface ---

  /// Cast to real
  operator real() const;

  /// Assignment from real
  const Scalar & operator=( real value );

  //--- Special functions

  /// Return a factory for the default linear algebra backend
  inline LinearAlgebraFactory & factory() const override;

  /// Get value
  real getval() const;

private:
  // Value of scalar
  real value{ 0.0 };
};

//-----------------------------------------------------------------------------
inline void Scalar::init( const GenericSparsityPattern & sparsity_pattern )
{
  MAYBE_UNUSED( sparsity_pattern );
  value = 0.0;
}

//-----------------------------------------------------------------------------
inline Scalar * Scalar::copy() const
{
  Scalar * s = new Scalar();
  s->value   = value;
  return s;
}

//-----------------------------------------------------------------------------
inline uint Scalar::rank() const
{
  return 0;
}

//-----------------------------------------------------------------------------
inline uint Scalar::size( uint ) const
{
  error( "The size() function is not available for scalars." );
  return 0;
}

//-----------------------------------------------------------------------------
inline void
  Scalar::get( real * block, const uint *, const uint * const * ) const
{
  block[0] = value;
}

//-----------------------------------------------------------------------------
inline void
  Scalar::set( const real * block, const uint *, const uint * const * )
{
  value = block[0];
}

//-----------------------------------------------------------------------------
inline void
  Scalar::add( const real * block, const uint *, const uint * const * )
{
  value += block[0];
}

//-----------------------------------------------------------------------------
inline void Scalar::zero()
{
  value = 0.0;
}

//-----------------------------------------------------------------------------
inline void Scalar::apply( FinalizeType )
{
  MPI::all_reduce_in_place< MPI::sum >( value );
}

//-----------------------------------------------------------------------------
inline void Scalar::disp( uint ) const
{
  prm( "Scalar value", value );
}

//-----------------------------------------------------------------------------
inline Scalar::operator real() const
{
  return value;
}

//-----------------------------------------------------------------------------
inline const Scalar & Scalar::operator=( real value )
{
  this->value = value;
  return *this;
}

//-----------------------------------------------------------------------------
inline LinearAlgebraFactory & Scalar::factory() const
{
  // Get backend from parameter system
  std::string backend = dolfin_get< std::string >( "linear algebra backend" );

#if ( HAVE_PETSC && HAVE_JANPACK )
  if ( backend == "PETSc" )
  {
    return PETScFactory::instance();
  }
  else if ( backend == "JANPACK" )
  {
    return JANPACKFactory::instance();
  }
#elif HAVE_PETSC
  return PETScFactory::instance();
#elif HAVE_JANPACK
  return JANPACKFactory::instance();
#endif
  error( "Linear algebra backend \"" + backend + "\" not available." );
}

//-----------------------------------------------------------------------------
inline real Scalar::getval() const
{
  return value;
}

//-----------------------------------------------------------------------------

}

#endif
