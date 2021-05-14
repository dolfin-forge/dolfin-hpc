// Copyright (C) 2010 Niclas Jansson
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_JANPACK_VEC_H
#define __DOLFIN_JANPACK_VEC_H

#include <dolfin/config/dolfin_config.h>

#ifdef HAVE_JANPACK

#include <dolfin/la/GenericVector.h>
#include <dolfin/common/Variable.h>

#ifdef HAVE_JANPACK_MPI
#define jp_vec_type jp_vec_t
#else
#define jp_vec_type char
#include <janpack/hybrid.h>
#endif

#include <janpack/vec.h>

namespace dolfin
{

//-----------------------------------------------------------------------------

  class JANPACKVec : public GenericVector, public Variable
  {
  public:

    /// Create empty vector
    JANPACKVec();

    /// Create vector of size N
    explicit JANPACKVec(size_t N, bool distributed);

    /// Copy constructor
    explicit JANPACKVec(const JANPACKVec& x);

    /// Destructor
    ~JANPACKVec();

    //--- Implementation of the GenericTensor interface ---

    /// Return copy of tensor
    JANPACKVec* copy() const;

    /// Set all entries to zero and keep any sparse structure
    void zero();

    /// Finalize assembly of tensor
    void apply(FinalizeType finaltype=FINALIZE);

    /// Display tensor
    void disp(size_t precision=2) const;

    //--- Implementation of the GenericVector interface ---

    /// Initialize vector of size N
    void init(size_t N);

    /// Initialize vector of size N and distribute if specified
    void init(size_t N, bool distributed);

    ///
    void init_ghosted( _ordered_set<size_t>& indices );

    /// Return size of vector
    size_t size() const;

    /// Return local size of vector
    size_t local_size() const;

    /// Return rank's offset into vector
    size_t offset() const;

    /// Get block of values
    void get(real* block, size_t m, const size_t* rows) const;

    /// Set block of values
    void set(const real* block, size_t m, const size_t* rows);

    /// Add block of values
    void add(const real* block, size_t m, const size_t* rows);

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
    real norm(VectorNormType type=l2) const;

    /// Return minimum value of vector
    real min() const;

    /// Return maximum value of vector
    real max() const;

    /// Return pointwise operator op of vector and given vector x
    void pointwise(const GenericVector& x, VectorPointwiseOp op=pw_min) const;

    /// Multiply vector by given number
    JANPACKVec& operator*= (real a);

    /// Divide vector by given number
    JANPACKVec& operator/= (real a);

    /// Multiply vector by given vector component-wise
    JANPACKVec& operator*= (const GenericVector& x);

    /// Add given vector
    JANPACKVec& operator+= (const GenericVector& x);

    /// Subtract given vector
    JANPACKVec& operator-= (const GenericVector& x);

    /// Assignment operator
    GenericVector& operator= (const GenericVector& x);

    /// Assignment operator
    JANPACKVec& operator= (real a);

    //--- Special functions ---

    /// Return linear algebra backend factory
    LinearAlgebraFactory& factory() const;

    //--- Special JANPACK functions ---

    /// Return JANPACK jp_vec_t pointer
    jp_vec_type *vec() const;

    /// Assignment operator
    JANPACKVec& operator= (const JANPACKVec& x);


    inline bool ghosted() { return is_ghosted;}

  private:

    // JANPACK vector pointer
#ifdef HAVE_JANPACK_MPI
    jp_vec_t _x;
    jp_vec_t *x_;
#else
    char x_[JP_VEC_SIZE_T];
#endif

    // True if the vector has ghost points
    bool is_ghosted;

    // True if the vector has been initialized (note we can't test
    // against x due to the opaque interface
    bool is_init;

    //    std::vector<int> ghost_indices;
#if (sun || __sun)
    _ordered_map<int, int> mapping;
#else
    _ordered_map<const int, int> mapping;
#endif


  };

//-----------------------------------------------------------------------------

} // namespace dolfin

#endif

#endif
