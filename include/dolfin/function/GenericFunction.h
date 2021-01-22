// Copyright (C) 2005-2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_GENERIC_FUNCTION_H
#define __DOLFIN_GENERIC_FUNCTION_H

#include <dolfin/common/Variable.h>
#include <dolfin/fem/Coefficient.h>
#include <dolfin/ufc/ufc.h>

namespace dolfin
{

class Mesh;
class Cell;

/**
 *  @class  GenericFunction
 *
 *  @brief  This class serves as a base class/interface for implementations of
 *          specific function representations.
 */

class GenericFunction : public Coefficient, public Variable
{

public:
  /// Constructor
  GenericFunction() = default;

  /// Destructor
  ~GenericFunction() override = default;

  //--- UFC INTERFACE ---------------------------------------------------------

  /// Evaluate function at given point in cell
  void evaluate( real *            values,
                 const real *      coordinates,
                 const ufc::cell & cell ) const override = 0;

  //--- INTERFACE -------------------------------------------------------------

  /// Evaluate function at given point
  void eval( real * values, const real * x ) const override = 0;

  /// Return the rank of the value space
  auto rank() const -> size_t override = 0;

  /// Return the dimension of the value space for axis i
  auto dim( size_t i ) const -> size_t override = 0;

  // Return the value size
  auto value_size() const -> size_t override = 0;

  /// Interpolate function to vertices of mesh
  virtual void interpolate_vertex_values( real * values ) const = 0;

  /// Interpolate function to finite element space on cell
  void interpolate( real *                      coefficients,
                    const ufc::cell &           cell,
                    const ufc::finite_element & finite_element,
                    const Cell & dolfin_cell ) const override = 0;

  /// Interpolate function to finite element space on facet
  void interpolate( real *                      coefficients,
                    const ufc::cell &           cell,
                    const ufc::finite_element & finite_element,
                    const Cell &                dolfin_cell,
                    size_t                      facet ) const override = 0;

  /// Synchronize
  void sync() override = 0;

  /// Display basic information
  void disp() const override = 0;

  //---------------------------------------------------------------------------

  /// Return the mesh
  virtual auto mesh() const -> Mesh & = 0;

  //---------------------------------------------------------------------------

  /// Delegate time dependency
  auto operator()( Time const & t ) -> GenericFunction &
  {
    this->sync( t );
    return *this;
  }

  //---------------------------------------------------------------------------

private:
  /// Time dependency hook
  void sync( Time const & t ) override = 0;
};

} /* namespace dolfin */

#endif /* __DOLFIN_GENERIC_FUNCTION_H */
