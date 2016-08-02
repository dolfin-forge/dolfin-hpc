//
//
//

#ifndef DOLFIN_MESH_ALGORITHM
#define DOLFIN_MESH_ALGORITHM

#include<dolfin/mesh/Mesh.h>

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
template<class Iterator, class Operator, class Transform, class Value>
void foreach(Mesh& mesh, Operator evaluator, Transform transform, Value& value)
{
  for (Iterator it(mesh); !it.end(); ++it)
  {
    value = transform(evaluator(*it), value);
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
