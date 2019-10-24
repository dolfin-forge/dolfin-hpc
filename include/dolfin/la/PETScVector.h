// Copyright (C) 2004-2008 Johan Hoffman, Johan Jansson and Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_PETSC_VECTOR_H
#define __DOLFIN_PETSC_VECTOR_H

#include <dolfin/config/dolfin_config.h>

#ifdef HAVE_PETSC

#include <dolfin/la/PETScObject.h>
#include <dolfin/la/GenericVector.h>
#include <dolfin/common/Variable.h>

#include <dolfin/log/LogStream.h>

#include <petscvec.h>

#include <set>
#include <map>

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
  ~PETScVector();

  //--- Implementation of the GenericTensor interface ---

  /// Return copy of tensor
  PETScVector* copy() const;

  /// Set all entries to zero and keep any sparse structure
  void zero();

  /// Finalize assembly of tensor
  void apply(FinalizeType finaltype = FINALIZE);

  /// Display tensor
  void disp(uint precision = 0) const;

  //--- Implementation of the GenericVector interface ---

  /// Initialize vector of local size N, distributed by default
  void init(uint N);

  /// Initialize vector of local size N, distributed if specified
  void init(uint N, bool distributed);

  /// Initialize ghost entries
  void init_ghosted(uint n, std::set<uint>& indices, std::map<uint, uint>& map);

  /// Return size of vector
  uint size() const;

  /// Return local size of vector
  uint local_size() const;

  /// Return rank's offset into vector
  uint offset() const;

  /// Get block of values
  void get(real* block, uint m, const uint* rows) const;

  /// Set block of values
  void set(const real* block, uint m, const uint* rows);

  /// Add block of values
  void add(const real* block, uint m, const uint* rows);

  /// Get all values
  void get(real* values) const;

  /// Set all values
  void set(real* values);

  /// Add values to each entry
  void add(real* values);

  /// Add multiple of given vector (AXPY operation)
  void axpy(real a, const GenericVector& x);

  /// Return inner product with given vector
  real inner(const GenericVector& v) const;

  /// Return norm of vector
  real norm(VectorNormType type = l2) const;

  /// Return minimum value of vector
  real min() const;

  /// Return maximum value of vector
  real max() const;

  /// Multiply vector by given number
  PETScVector& operator*=(real a);

  /// Divide vector by given number
  PETScVector& operator/=(real a);

  /// Multiply vector by given vector component-wise
  PETScVector& operator*=(const GenericVector& x);

  /// Add given vector
  PETScVector& operator+=(const GenericVector& x);

  /// Subtract given vector
  PETScVector& operator-=(const GenericVector& x);

  /// Assignment operator
  PETScVector& operator=(const GenericVector& x);

  /// Assignment operator
  PETScVector& operator=(real a);

  //--- Special functions ---

  /// Return linear algebra backend factory
  LinearAlgebraFactory& factory() const;

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
  Vec x_;

  // True if the vector is distributed
  bool is_distributed_;

  // True if the vector has ghost points
  bool is_ghosted_;

#if (sun || __sun)
  typedef std::map<int, int> GhostMapping;
#else
  typedef std::map<const int, int> GhostMapping;
#endif
  GhostMapping mapping_;

};

} /* namespace dolfin */

#endif /* HAVE_PETSC */

#endif /* __DOLFIN_PETSC_VECTOR_H */
