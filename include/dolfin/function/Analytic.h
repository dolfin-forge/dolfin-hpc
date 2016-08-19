#ifndef __DOLFIN_FUNCTION_ANALYTIC_H_
#define __DOLFIN_FUNCTION_ANALYTIC_H_

#include <dolfin/function/GenericFunction.h>

#include <dolfin/log/log.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/Vertex.h>

using dolfin::GenericFunction;
using dolfin::Cell;
using dolfin::Mesh;
using dolfin::VertexIterator;

namespace dolfin
{

template<class T>
class Analytic : public GenericFunction
{

public:

  /// Constructor
  Analytic<T>(Mesh& mesh) :
      mesh_(mesh),
      evaluant_()
  {
  }

  /// Constructor
  Analytic<T>(Mesh& mesh, T expr) :
      mesh_(mesh),
      evaluant_(expr)
  {
  }

  /// Destructor
  ~Analytic<T>()
  {
  }

  /// Copy constructor
  Analytic<T>(Analytic<T> const& other) :
      mesh_(const_cast<Mesh&>(other.mesh_)),
      evaluant_(other.evaluant_)
  {
  }

  //--- UFC INTERFACE ---------------------------------------------------------

  /// Evaluate function at given point in cell
  inline void evaluate(real* values, const real* coordinates,
                       const ufc::cell& cell) const
  {
    evaluant_.eval(values, coordinates);
  }

  //--- Expression INTERFACE --------------------------------------------------

  ///
  inline void evaluate(uint n, real* values, const real* coordinates,
                       const ufc::cell& cell) const
  {
    for (uint i = 0; i < n; ++i)
    {
      evaluant_.eval(&values[i * value_size()],
                     &coordinates[i * cell.geometric_dimension]);
    }
  }

  /// Evaluate function at given point
  inline void eval(real* values, const real* x) const
  {
    evaluant_.eval(values, x);
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

  /// Value size
  inline uint value_size() const
  {
    return evaluant_.value_size();
  }

  ///
  Analytic& operator()(Time const& t)
  {
    evaluant_(t);
    return *this;
  }

  //--- GenericFunction INTERFACE ---------------------------------------------

  /// Return the mesh
  inline Mesh& mesh() const
  {
    return mesh_;
  }

  /// Interpolate function to vertices of mesh
  inline void interpolate_vertex_values(real* values) const
  {
    dolfin_assert(values);

    // Get size of value (number of entries in tensor value)
    uint const size = this->value_size();
    real * local_values = new real[size];
    uint const num_verts = mesh_.size(0);
    for (VertexIterator vertex(mesh_); !vertex.end();
         ++vertex)
    {
      // Evaluate at function at vertex
      evaluant_.eval(local_values, vertex->x());

      // Copy values to array of vertex values
      for (uint i = 0; i < size; ++i)
      {
        values[i * num_verts + vertex->index()] = local_values[i];
      }
    }
    delete[] local_values;
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
  inline void interpolate(real* coefficients, const ufc::cell& cell,
                          const ufc::finite_element& finite_element,
                          const Cell& dolfin_cell, uint facet) const
  {
    dolfin_assert(coefficients);
    finite_element.evaluate_dofs(coefficients, *this, cell);
  }

  /// Display basic information
  inline void disp() const
  {
    section("Analytic");
    evaluant_.disp();
    end();
    skip();
  }

  /// Synchronize
  inline void sync()
  {
    // Do nothing
  }

public:

  Mesh& mesh_;
  T evaluant_;

};

} /* namespace licorne */

#endif /* __DOLFIN_FUNCTION_ANALYTIC_H_ */
