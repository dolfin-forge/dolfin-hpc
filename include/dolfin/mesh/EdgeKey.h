//
//
//

#ifndef __DOLFIN_MESH_EDGE_KEY
#define __DOLFIN_MESH_EDGE_KEY

#include <dolfin/mesh/Edge.h>

#include <map>

namespace dolfin
{

struct EdgeKey : public std::pair<uint, uint>
{

  /// An edge contains a pair of vertices
  EdgeKey() :
      std::pair<uint, uint>(),
      idx(0)
  {
  }

  /// An edge contains a pair of vertices
  EdgeKey(Edge const& e) :
      std::pair<uint, uint>(),
      idx(std::rand())
  {
    uint const * v = e.entities(0);
    this->first  = (v[0] < v[1] ? v[0] : v[1]);
    this->second = (v[0] < v[1] ? v[1] : v[0]);
  }

  /// An edge contains a pair of vertices
  EdgeKey(uint v0, uint v1) :
      std::pair<uint, uint>(v0 < v1 ? v0 : v1, v0 < v1 ? v1 : v0),
      idx(std::rand())
  {
  }

  /// Construct a key from edge vertices
  inline void set(Edge const& e)
  {
    uint const * v = e.entities(0);
    this->first  = (v[0] < v[1] ? v[0] : v[1]);
    this->second = (v[0] < v[1] ? v[1] : v[0]);
    idx = std::rand();
  }

  /// Construct a key from edge vertices
  inline void set(uint const * v)
  {
    this->first  = (v[0] < v[1] ? v[0] : v[1]);
    this->second = (v[0] < v[1] ? v[1] : v[0]);
    idx = std::rand();
  }

  /// Construct a key from edge vertices
  inline void set(uint v0, uint v1)
  {
    this->first  = (v0 < v1 ? v0 : v1);
    this->second = (v0 < v1 ? v1 : v0);
    idx = std::rand();
  }

  ///
  inline size_t hash() const
  {
    return (static_cast<std::size_t>(this->first)
        ^ static_cast<std::size_t>(this->second));
  }

  ///
  uint idx;

};

} /* namespace dolfin */

#if (HAVE_TR1_UNORDERED_MAP && HAVE_TR1_UNORDERED_SET) || \
    (__IBMCPP__ && __IBMCPP_TR1__)|| \
    (ENABLE_BOOST_TR1)

namespace std
{

namespace tr1
{

template<>
struct hash<dolfin::EdgeKey>
{
  inline std::size_t operator()(dolfin::EdgeKey const& e) const
  {
    return e.hash();
  }
};

} /* namespace tr1 */

} /* namespace std */

#endif

#endif /* __DOLFIN_MESH_EDGE_KEY */
