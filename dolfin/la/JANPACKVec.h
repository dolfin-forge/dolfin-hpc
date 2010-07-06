// Copyright (C) 2010 Niclas Jansson
// Licensed under the GNU LGPL Version 2.1.
//

#ifndef __JANPACK_VEC_H
#define __JANPACK_VEC_H

#ifdef HAS_JANPACK

#include <vec.h>

#include <dolfin/log/LogStream.h>
#include <dolfin/common/Variable.h>
#include "GenericVector.h"

#include <dolfin/common/Array.h>
#include <set>
#include <map>

namespace dolfin
{

  class JANPACKVec : public GenericVector, public Variable
  {
  public:

    /// Create empty vector
    JANPACKVec();

    /// Create vector of size N
    explicit JANPACKVec(uint N);

    /// Copy constructor
    explicit JANPACKVec(const JANPACKVec& x);

    /// Destructor
    virtual ~JANPACKVec();

    //--- Implementation of the GenericTensor interface ---

    /// Return copy of tensor
    virtual JANPACKVec* copy() const;

    /// Set all entries to zero and keep any sparse structure
    virtual void zero();

    /// Finalize assembly of tensor
    virtual void apply(FinalizeType finaltype=FINALIZE);

    /// Display tensor
    virtual void disp(uint precision=2) const;

    //--- Implementation of the GenericVector interface ---

    /// Initialize vector of size N
    virtual void init(uint N);

    virtual void init_ghosted(uint n, std::set<uint>& indices,
			      std::map<uint, uint>& map);

    /// Return size of vector
    virtual uint size() const;

    /// Return local size of vector
    virtual uint local_size() const;

    /// Get block of values
    virtual void get(real* block, uint m, const uint* rows) const;

    /// Set block of values
    virtual void set(const real* block, uint m, const uint* rows);

    /// Add block of values
    virtual void add(const real* block, uint m, const uint* rows);

    /// Get all values
    virtual void get(real* values) const;

    /// Set all values
    virtual void set(real* values);

    /// Add values to each entry
    virtual void add(real* values);

    /// Add multiple of given vector (AXPY operation)
    virtual void axpy(real a, const GenericVector& x); 

    /// Return inner product with given vector
    virtual real inner(const GenericVector& v) const;

    /// Return norm of vector
    virtual real norm(VectorNormType type=l2) const;

    /// Return minimum value of vector
    virtual real min() const;

    /// Return maximum value of vector
    virtual real max() const;

    /// Multiply vector by given number
    virtual const JANPACKVec& operator*= (real a);

    /// Divide vector by given number
    virtual const JANPACKVec& operator/= (real a);

    /// Add given vector
    virtual const JANPACKVec& operator+= (const GenericVector& x);

    /// Subtract given vector
    virtual const JANPACKVec& operator-= (const GenericVector& x);

    /// Assignment operator
    virtual const GenericVector& operator= (const GenericVector& x);

    /// Assignment operator
    virtual const JANPACKVec& operator= (real a);

    //--- Special functions ---

    /// Return linear algebra backend factory
    virtual LinearAlgebraFactory& factory() const;

    //--- Special JANPACK functions ---

    /// Return JANPACK Vec_ pointer
    Vec_ *vec() const;

    /// Assignment operator
    const JANPACKVec& operator= (const JANPACKVec& x);


    inline bool ghosted() { return is_ghosted;}

  private:

    // JANPACK vector pointer
    Vec_ *x;
    Vec_ _x;
       
    // True if we don't own the vector x points to
    bool is_view;

    // True if the vector has ghost points
    bool is_ghosted;
    //    Array<int> ghost_indices;
#if (sun || __sun)    
    std::map<int, int> mapping;
#else
    std::map<const int, int> mapping;
#endif


  };
  
}

#endif

#endif
