#ifndef __DOLFIN_FUNCTION_IMPL_OPS_H_
#define __DOLFIN_FUNCTION_IMPL_OPS_H_

#include <dolfin/function/ValueSpace.h>
#include <dolfin/fem/UFCCell.h>
#include <dolfin/mesh/Vertex.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
template <class T, class E, class V>
struct EntityOp : public V
{
  inline void evaluate(real* values, const real* x, const UFCCell& cell) const
  {
    static_cast<T const *>(this)->operator()(*cell.cell, values);
  }

  inline void operator()(E const& e, real *x) const
  {
    static_cast<T const *>(this)->operator()(e, x);
  }
};
//-----------------------------------------------------------------------------

//--- IMPLEMENTATION ----------------------------------------------------------

//-----------------------------------------------------------------------------
// Circumradius
//-----------------------------------------------------------------------------
template <class E>
struct Circumradius : public EntityOp<Circumradius<E>, E, ValueSpace<> >
{
  inline void operator()(E const& e, real *x) const { x[0] = e.circumradius(); }
};

//-----------------------------------------------------------------------------
// Coordinates
//-----------------------------------------------------------------------------
template<uint D>
struct Coordinate : public EntityOp<Coordinate<D>, Vertex, ValueSpace<> >
{
  inline void operator()(Vertex const& v, real *x) const { x[0] = v.x()[D]; }
};

template<uint D>
struct Coordinates : public EntityOp<Coordinates<D>, Vertex, ValueSpace<D> >
{
  inline void operator()(Vertex const& v, real *x) const
    { std::copy(v.x(), v.x() + D, x); }
};

//-----------------------------------------------------------------------------
// Diameter
//-----------------------------------------------------------------------------
template <class E>
struct Diameter: public EntityOp<Diameter<E>, E, ValueSpace<> >
{
  inline void operator()(E const& e, real *x) const { x[0] = e.diameter(); }
};

//-----------------------------------------------------------------------------
// Volume
//-----------------------------------------------------------------------------
template <class E>
struct Volume: public EntityOp<Volume<E>, E, ValueSpace<> >
{
  inline void operator()(E const& e, real *x) const { x[0] = e.volume(); }
};

//-----------------------------------------------------------------------------

} /* namespace dolfin */

#endif /* __DOLFIN_FUNCTION_IMPL_OPS_H_ */
