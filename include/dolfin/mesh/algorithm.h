
#ifndef DOLFIN_MESH_ALGORITHM
#define DOLFIN_MESH_ALGORITHM

#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/entities/iterators/CellIterator.h>
#include <dolfin/mesh/entities/iterators/FacetIterator.h>

namespace dolfin
{

/**
 *  A bunch of convenient functions for testing and development.
 *
 */

//-----------------------------------------------------------------------------
template<class Iterator, class Operator>
void foreach(Mesh& mesh, Operator transform)
{
  for (Iterator it(mesh); !it.end(); ++it)
  {
    transform(*it);
  }
}

//-----------------------------------------------------------------------------
// Functors
template<class Iterator, class Operator, class Transform, class Value>
void foreach(Mesh& mesh, Operator evaluator, Transform transform, Value& value)
{
  for (Iterator it(mesh); !it.end(); ++it)
  {
    value = transform(evaluator(*it), value);
  }
}

//-----------------------------------------------------------------------------
// Unary ops (size_t)
template<class Iterator, size_t const& (*T)(size_t const&), class Operator>
void foreach(Mesh& mesh, Operator evaluator, size_t& value)
{
  for (Iterator it(mesh); !it.end(); ++it)
  {
    value = T(evaluator(*it));
  }
}

//-----------------------------------------------------------------------------
// Unary ops (real)
template<class Iterator, real const& (*T)(real const&), class Operator>
void foreach(Mesh& mesh, Operator evaluator, real& value)
{
  for (Iterator it(mesh); !it.end(); ++it)
  {
    value = T(evaluator(*it));
  }
}

//-----------------------------------------------------------------------------
// Binary ops (size_t)
template<class Iterator, size_t const& (*T)(size_t const&, size_t const&), class Operator>
void foreach(Mesh& mesh, Operator evaluator, size_t& value)
{
  for (Iterator it(mesh); !it.end(); ++it)
  {
    value = T(evaluator(*it), value);
  }
}

//-----------------------------------------------------------------------------
// Binary ops (real)
template<class Iterator, real const& (*T)(real const&, real const&), class Operator>
void foreach(Mesh& mesh, Operator evaluator, real& value)
{
  for (Iterator it(mesh); !it.end(); ++it)
  {
    value = T(evaluator(*it), value);
  }
}

//-----------------------------------------------------------------------------
template<class Iterator, class Operator, class Value>
void accumulate(Mesh& mesh, Operator& evaluator, Value& value)
{
  for (Iterator it(mesh); !it.end(); ++it)
  {
    value += evaluator(*it);
  }
}

//-----------------------------------------------------------------------------

}

#endif /* DOLFIN_MESH_ALGORITHM */
