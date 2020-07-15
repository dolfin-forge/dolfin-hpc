
#ifndef __DOLFIN_MESH_EDGE_KEY
#define __DOLFIN_MESH_EDGE_KEY

#include <dolfin/mesh/Edge.h>

namespace dolfin
{

template <typename T=uint>
struct EdgeKey : public std::pair<T, T>
{

  /// An edge contains a pair of vertices
  EdgeKey() :
      std::pair<T, T>(),
      idx(0)
  {
  }

  /// An edge contains a pair of vertices
  EdgeKey(Edge const& e) :
      std::pair<T, T>(),
      idx(std::rand())
  {
    Array<uint> const & v = e.entities(0);
    this->first  = static_cast<T>( v[0] < v[1] ? v[0] : v[1] );
    this->second = static_cast<T>( v[0] < v[1] ? v[1] : v[0] );
  }

  /// An edge contains a pair of vertices
  EdgeKey(T v0, T v1) :
      std::pair<T, T>(v0 < v1 ? v0 : v1, v0 < v1 ? v1 : v0),
      idx(std::rand())
  {
  }

  /// Construct a key from edge vertices
  inline void set(Edge const& e)
  {
    Array<uint> const & v = e.entities(0);
    this->first  = static_cast<T>( v[0] < v[1] ? v[0] : v[1] );
    this->second = static_cast<T>( v[0] < v[1] ? v[1] : v[0] );
    idx = std::rand();
  }

  /// Construct a key from edge vertices
  inline void set(T const * v)
  {
    this->first  = (v[0] < v[1] ? v[0] : v[1]);
    this->second = (v[0] < v[1] ? v[1] : v[0]);
    idx = std::rand();
  }

  /// Construct a key from edge vertices
  inline void set(T v0, T v1)
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

namespace std
{

template<typename T>
struct hash<dolfin::EdgeKey<T> >
{
  inline std::size_t operator()(dolfin::EdgeKey<T> const& e) const
  {
    return e.hash();
  }
};

} /* namespace std */

#endif /* __DOLFIN_MESH_EDGE_KEY */
