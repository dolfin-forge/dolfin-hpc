// Copyright (C) 2004-2008 Johan Hoffman, Johan Jansson and Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_PETSC_VECTOR_H
#define __DOLFIN_PETSC_VECTOR_H

#include <dolfin/common/Variable.h>
#include <dolfin/config/dolfin_config.h>
#include <dolfin/la/GenericVector.h>
#include <dolfin/log/LogStream.h>

#ifdef HAVE_PETSC

#include <dolfin/la/PETScObject.h>

#include <petscvec.h>

namespace dolfin
{

/// This class provides a simple vector class based on PETSc.
/// It is a simple wrapper for a PETSc vector pointer (Vec)
/// implementing the GenericVector interface.
///
/// The interface is intentionally simple. For advanced usage,
/// access the PETSc Vec pointer using the function vec() and
/// use the standard PETSc interface.

class PETScVector : public GenericVector, public PETScObject, public Variable
{

public:

  /// Create empty vector
  PETScVector();

  /// Create vector of local size N
  explicit PETScVector(uint N, bool distributed = true);

  /// Copy constructor
  explicit PETScVector(const PETScVector& x);

  /// Create vector from given PETSc Vec pointer
  explicit PETScVector(Vec x);

  /// Destructor
  ~PETScVector() override;

  //--- Implementation of the GenericTensor interface ---

  /// Return copy of tensor
  PETScVector* copy() const override;

  /// Set all entries to zero and keep any sparse structure
  void zero() override;

  /// Finalize assembly of tensor
  void apply(FinalizeType finaltype = FINALIZE) override;

  /// Display tensor
  void disp(uint precision = 0) const override;

  //--- Implementation of the GenericVector interface ---

  /// Initialize vector of local size N, distributed by default
  void init(uint N) override;

  /// Initialize vector of local size N, distributed if specified
  void init(uint N, bool distributed) override;

  /// Initialize ghost entries
  void init_ghosted(uint n, _ordered_set<uint>& indices,
                    _ordered_map<uint, uint>& map) override;

  /// Return size of vector
  uint size() const override;

  /// Return local size of vector
  uint local_size() const override;

  /// Return rank's offset into vector
  uint offset() const override;

  /// Get block of values
  void get(real* block, uint m, const uint* rows) const override;

  /// Set block of values
  void set(const real* block, uint m, const uint* rows) override;

  /// Add block of values
  void add(const real* block, uint m, const uint* rows) override;

  /// Get all values
  void get(real* values) const override;

  /// Set all values
  void set(real* values) override;

  /// Add values to each entry
  void add(real* values) override;

  /// Add multiple of given vector (y=a*x+y)
  void axpy( real a, const GenericVector & x ) override;

  /// Add multiple of given vector (y=a*x+b*y)
  void axpby( real a, const GenericVector & x,
              real b ) override;

  /// Add multiple of given vector (w=a*x+y)
  void waxpy( real a, const GenericVector & x,
                      const GenericVector & y ) override;

  /// Add multiple of given vector (z=a*x+b*y+c*z)
  void axpbypcz( real a, const GenericVector & x,
                 real b, const GenericVector & y,
                 real c ) override;

  /// Return inner product with given vector
  real inner(const GenericVector& v) const override;

  /// Return norm of vector
  real norm(VectorNormType type = l2) const override;

  /// Return minimum value of vector
  real min() const override;

  /// Return maximum value of vector
  real max() const override;

  /// Return pointwise operator op of vector and given vector x
  void pointwise(const GenericVector& x, VectorPointwiseOp op=pw_min) const override;

  /// Multiply vector by given number
  PETScVector& operator*=(real a) override;

  /// Divide vector by given number
  PETScVector& operator/=(real a) override;

  /// Multiply vector by given vector component-wise
  PETScVector& operator*=(const GenericVector& x) override;

  /// Add given vector
  PETScVector& operator+=(const GenericVector& x) override;

  /// Subtract given vector
  PETScVector& operator-=(const GenericVector& x) override;

  /// Assignment operator
  PETScVector& operator=(const GenericVector& x) override;

  /// Assignment operator
  PETScVector& operator=(real a) override;

  //--- Special functions ---

  /// Return linear algebra backend factory
  LinearAlgebraFactory& factory() const override;

  //--- Special PETSc functions ---

  /// Assignment operator
  PETScVector& operator=(PETScVector const& x);

  /// Return PETSc Vec pointer
  Vec vec() const;

  inline bool ghosted()
  {
    return is_ghosted_;
  }

private:


  //
  void clear();

  // PETSc Vec pointer
  Vec x_{nullptr};

  // True if the vector is distributed
  bool is_distributed_{false};

  // True if the vector has ghost points
  bool is_ghosted_{false};

  using GhostMapping = _map<int, int>;
  GhostMapping mapping_;

};

//-----------------------------------------------------------------------------
inline void PETScVector::clear()
{
  if (x_)
  {
#if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR > 1
    VecDestroy(&x_);
#else
    VecDestroy(x_);
#endif
    is_ghosted_ = false;
    is_distributed_ = false;
  }
}
//-----------------------------------------------------------------------------
inline void PETScVector::init(uint N)
{
  init(N, true);
}
//-----------------------------------------------------------------------------
inline PETScVector* PETScVector::copy() const
{
  return new PETScVector(*this);
}
//-----------------------------------------------------------------------------
inline void PETScVector::set(const real* block, uint m, const uint* rows)
{
  dolfin_assert(x_);
  VecSetValues(x_, static_cast<int>(m),
               reinterpret_cast<int*>(const_cast<uint*>(rows)), block,
               INSERT_VALUES);
}
//-----------------------------------------------------------------------------
inline void PETScVector::add(const real* block, uint m, const uint* rows)
{
  dolfin_assert(x_);
  VecSetValues(x_, static_cast<int>(m),
               reinterpret_cast<int*>(const_cast<uint*>(rows)), block,
               ADD_VALUES);
}
//-----------------------------------------------------------------------------
inline void PETScVector::apply(FinalizeType)
{

  VecAssemblyBegin(x_);
  VecAssemblyEnd(x_);

  if (is_ghosted_)
  {
    VecGhostUpdateBegin(x_, INSERT_VALUES, SCATTER_FORWARD);
    VecGhostUpdateEnd(x_, INSERT_VALUES, SCATTER_FORWARD);
  }
}
//-----------------------------------------------------------------------------
inline void PETScVector::zero()
{
  dolfin_assert(x_);
  real a = 0.0;
  VecSet(x_, a);
}
//-----------------------------------------------------------------------------
inline uint PETScVector::size() const
{
  if(x_ == nullptr) return 0;
  int n = 0;
  VecGetSize(x_, &n);
  return static_cast<uint>(n);
}
//-----------------------------------------------------------------------------
inline uint PETScVector::local_size() const
{
  dolfin_assert(x_);
  int n = 0;
  VecGetLocalSize(x_, &n);
  return static_cast<uint>(n);
}
//-----------------------------------------------------------------------------
inline uint PETScVector::offset() const
{
  dolfin_assert(x_);
  int low, high;
  VecGetOwnershipRange(x_, &low, &high);
  return static_cast<uint>(low);
}
//-----------------------------------------------------------------------------
inline PETScVector& PETScVector::operator=(const GenericVector& v)
{
  *this = v.down_cast<PETScVector>();
  return *this;
}
//-----------------------------------------------------------------------------
inline PETScVector& PETScVector::operator=(PETScVector const& v)
{
  if (&v != this)
  {
    dolfin_assert(v.x_);
    init(v.local_size(), v.is_distributed_);
    VecCopy(v.x_, x_);
  }
  return *this;
}
//-----------------------------------------------------------------------------
inline PETScVector& PETScVector::operator=(real a)
{
  dolfin_assert(x_);
  VecSet(x_, a);
  return *this;
}
//-----------------------------------------------------------------------------
inline PETScVector& PETScVector::operator*=(const GenericVector& y)
{
  dolfin_assert(x_);
  PETScVector const& v = y.down_cast<PETScVector>();
  dolfin_assert(v.x_);

  if (size() != v.size())
  {
    error("Vectors must have the same size for componentwise multiplication.");
  }

  VecPointwiseMult(x_, x_, v.x_);

  return *this;
}
//-----------------------------------------------------------------------------
inline PETScVector& PETScVector::operator+=(const GenericVector& x)
{
  this->axpy(1.0, x);
  return *this;
}
//-----------------------------------------------------------------------------
inline PETScVector& PETScVector::operator-=(const GenericVector& x)
{
  this->axpy(-1.0, x);
  return *this;
}
//-----------------------------------------------------------------------------
inline PETScVector& PETScVector::operator*=(const real a)
{
  dolfin_assert(x_);
  VecScale(x_, a);

  return *this;
}
//-----------------------------------------------------------------------------
inline PETScVector& PETScVector::operator/=(const real a)
{
  dolfin_assert(x_);
  dolfin_assert(a != 0.0);

  const real b = 1.0 / a;
  VecScale(x_, b);

  return *this;
}
//-----------------------------------------------------------------------------
inline real PETScVector::inner(const GenericVector& y) const
{
  dolfin_assert(x_);

  PETScVector const& v = y.down_cast<PETScVector>();
  dolfin_assert(v.x_);

  real a;
  VecDot(v.x_, x_, &a);

  return a;
}
//-----------------------------------------------------------------------------
inline void PETScVector::axpy( real a, const GenericVector & x )
{
  dolfin_assert( x_ );

  PETScVector const & v = x.down_cast< PETScVector >();
  dolfin_assert( v.x_ );

  if ( size() != v.size() )
  {
    error( "The vectors must be of the same size to apply AXPY." );
  }

  VecAXPY( x_, a, v.x_ );
}
//-----------------------------------------------------------------------------
inline void PETScVector::axpby( real a, const GenericVector & x, real b )
{
  dolfin_assert( x_ );

  PETScVector const & v = x.down_cast< PETScVector >();
  dolfin_assert( v.x_ );

  if ( size() != v.size() )
  {
    error( "The vectors must be of the same size to apply AXPBY." );
  }

  VecAXPBY( x_, a, b, v.x_ );
}
//-----------------------------------------------------------------------------
inline void PETScVector::waxpy( real a, const GenericVector & x,
                                        const GenericVector & y )
{
  dolfin_assert( x_ );

  PETScVector const & v = x.down_cast< PETScVector >();
  dolfin_assert( v.x_ );

  PETScVector const & u = y.down_cast< PETScVector >();
  dolfin_assert( u.x_ );

  if ( size() != v.size() or size() != u.size() )
  {
    error( "The vectors must be of the same size to apply WAXPY." );
  }

  VecWAXPY( x_, a, v.x_, u.x_ );
}
//-----------------------------------------------------------------------------
inline void PETScVector::axpbypcz( real a, const GenericVector & x,
                                   real b, const GenericVector & y,
                                   real c )
{
  dolfin_assert( x_ );

  PETScVector const & v = x.down_cast< PETScVector >();
  dolfin_assert( v.x_ );

  PETScVector const & u = y.down_cast< PETScVector >();
  dolfin_assert( u.x_ );

  if ( size() != v.size() or size() != u.size() )
  {
    error( "The vectors must be of the same size to apply VecAXPBYPCZ." );
  }

  VecAXPBYPCZ( x_, a, b, c, v.x_, u.x_ );
}
//-----------------------------------------------------------------------------
inline real PETScVector::min() const
{
  real value = 0.0;

  VecMin(x_, PETSC_NULL, &value);

  return value;
}
//-----------------------------------------------------------------------------
inline real PETScVector::max() const
{
  real value = 0.0;

  VecMax(x_, PETSC_NULL, &value);

  return value;
}

//-----------------------------------------------------------------------------
inline Vec PETScVector::vec() const
{
  return x_;
}

//-----------------------------------------------------------------------------

} /* namespace dolfin */

#endif /* HAVE_PETSC */

#endif /* __DOLFIN_PETSC_VECTOR_H */
