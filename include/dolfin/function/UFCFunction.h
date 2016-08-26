// Copyright (C) 2005-2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2005-11-28
// Last changed: 2008-03-17

#ifndef __DOLFIN_UFC_FUNCTION_H
#define __DOLFIN_UFC_FUNCTION_H

#include <dolfin/function/GenericFunction.h>
#include <dolfin/evolution/TimeDependent.h>

#include <dolfin/fem/UFCCell.h>
#include <dolfin/log/log.h>

namespace dolfin
{

class Mesh;
class Cell;

/**
 *  @class  UFCFunction
 *
 *  @brief  This class serves as a base class/interface for implementations of
 *          user-defined functions using UFC cell and facet information.
 */

template<class T>
class UFCFunction : public GenericFunction, public TimeDependent
{

public:

  /// Constructor
  UFCFunction(Mesh& mesh) :
  GenericFunction(),
  mesh_(mesh),
  facet_(-1)
  {
  }

  /// Destructor
  ~UFCFunction()
  {

  }

  /// Evaluate function at given point coordinate in UFC cell
  inline void eval(real* values, const real* x, UFCCell const& cell) const
  {
    evaluant_.eval(values, x, cell);
  }

  //--- UFC INTERFACE ---------------------------------------------------------

  ///
  inline void evaluate(real* values, const real* coordinates, const ufc::cell& cell) const
  {
    evaluant_.eval(values, coordinates, static_cast<UFCCell const&>(cell));
  }

  //--- GenericFunction INTERFACE ---------------------------------------------

  /// Return the mesh
  Mesh& mesh() const
  {
    return mesh_;
  }

  ///
  inline void evaluate(uint n, real* values, const real* coordinates,
                       const ufc::cell& cell) const
  {
    for (uint i = 0; i < n; ++i)
    {
      evaluant_.eval(&values[i * value_size()],
                     &coordinates[i * cell.geometric_dimension],
                     static_cast<UFCCell const&>(cell));
    }
  }

  /// Evaluate function at given point coordinate, cell is searched for
  void eval(real* values, const real* x) const
  {
    error("UFCFunction : eval unimplemented");
  }

  /// Return the rank of the value space
  inline uint rank() const
  {
    return evaluant_.rank();
  }

  /// Return the dimension of the value space for axis i
  inline uint dim(uint i) const
  {
    return evaluant_.dim(i);
  }

  // Return the value size
  uint value_size() const { return Coefficient::value_size(); }

  /// Interpolate function to vertices of mesh
  void interpolate_vertex_values(real* values) const
  {
    error("UFCFunction : interpolate_vertex_values unimplemented");
  }

  /// Interpolate function to finite element space on cell
  inline void interpolate(real* coefficients, const ufc::cell& cell,
                          const ufc::finite_element& finite_element,
                          const Cell& dolfin_cell) const
  {
    dolfin_assert(coefficients);
    finite_element.evaluate_dofs(coefficients, *this, cell);
  }

  /// Interpolate function to finite element space on facet
  inline  void interpolate(real* coefficients, const ufc::cell& cell,
                           const ufc::finite_element& finite_element,
                           const Cell& dolfin_cell, uint facet) const
  {
    this->facet_ = facet;
    dolfin_assert(coefficients);
    finite_element.evaluate_dofs(coefficients, *this, cell);
    this->facet_ = -1;
  }

  //---------------------------------------------------------------------------

  /// UFCFunction implements the time dependency
  inline UFCFunction<T> const& operator()(Time const& t) const
  {
    TimeDependent::operator ()(t);
    return *this;
  }

  //---------------------------------------------------------------------------

  /// Display basic information
  void disp() const
  {
    section("UFCFunction");
    end();
  }

  /// Synchronize
  virtual void sync() { /* No-op */ }

protected:

  inline int facet() const { return facet_; }

private:

  /// Time synchronization hook
  inline void sync(Time const& t) { TimeDependent::operator ()(t); }

  Mesh& mesh_;
  T evaluant_;

  mutable int facet_;

};

} /* namespace dolfin */

#endif /* __DOLFIN_UFC_FUNCTION_H */
