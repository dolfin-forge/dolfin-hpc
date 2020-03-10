// Copyright (C) 2007-2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_SCALAR_H
#define __DOLFIN_SCALAR_H


#include <dolfin/config/dolfin_config.h>
#include <dolfin/parameter/parameters.h>
#include "GenericTensor.h"
#include <dolfin/main/MPI.h>
#include <dolfin/common/maybe_unused.h>

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
    Scalar() : GenericTensor(), value(0.0)
    {}

    /// Destructor
    virtual ~Scalar()
    {}

    //--- Implementation of the GenericTensor interface ---

    /// Initialize zero tensor using sparsity pattern
    void init(const GenericSparsityPattern& sparsity_pattern)
    { MAYBE_UNUSED(sparsity_pattern); value = 0.0; }

    /// Return copy of tensor
    virtual Scalar* copy() const
    { Scalar* s = new Scalar(); s->value = value; return s; }

    /// Return tensor rank (number of dimensions)
    uint rank() const
    { return 0; }

    /// Return size of given dimension
    uint size(uint) const
    { error("The size() function is not available for scalars."); return 0; }

    /// Get block of values
    void get(real* block, const uint*, const uint * const *) const
    { block[0] = value; }

    /// Set block of values
    void set(const real* block, const uint*, const uint * const *)
    { value = block[0]; }

    /// Add block of values
    void add(const real* block, const uint*, const uint * const *)
    { value += block[0]; }

    /// Set all entries to zero and keep any sparse structure
    void zero()
    { value = 0.0; }

    /// Finalize assembly of tensor
    void apply(FinalizeType)
    {
      real tmp = value;
      MPI::all_reduce<MPI::sum>(tmp, value);
    }


    /// Display tensor
    void disp(uint) const
    { prm("Scalar value", value); }

    //--- Scalar interface ---

    /// Cast to real
    operator real() const
    { return value; }

    /// Assignment from real
    const Scalar& operator= (real value)
    { this->value = value; return *this; }

    //--- Special functions

    /// Return a factory for the default linear algebra backend
    inline LinearAlgebraFactory& factory() const
    {

      // Get backend from parameter system
      std::string backend = dolfin_get<std::string>("linear algebra backend");

#if (HAVE_PETSC && HAVE_JANPACK)
      if (backend == "PETSc")
      {
	return PETScFactory::instance();
      }
      else if (backend == "JANPACK")
      {
	return JANPACKFactory::instance();
      }
#elif HAVE_PETSC
      return PETScFactory::instance();
#elif HAVE_JANPACK
      return JANPACKFactory::instance();
#endif
      error("Linear algebra backend \"" + backend + "\" not available.");
    }

    /// Get value
    real getval() const
    { return value; }

  private:

    // Value of scalar
    real value;

  };

}

#endif
