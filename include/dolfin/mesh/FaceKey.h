//
//
//

#ifndef __DOLFIN_MESH_FACE_KEY
#define __DOLFIN_MESH_FACE_KEY

#include <dolfin/mesh/EdgeKey.h>
#include <dolfin/mesh/Face.h>

#include <set>

namespace dolfin
{

struct FaceKey : public std::set<EdgeKey>
{
  /// An face contains a bunch of vertices
  FaceKey() :
      std::set<EdgeKey>(),
      idx(0)
  {
  }

  /// An face contains a bunch of vertices
  FaceKey(Face const& f) :
    std::set<EdgeKey>(),
    idx(std::rand())
  {
    uint const * v = f.entities(0);
    uint const n = f.num_entities(0);
    for (uint i = 0; i < n; ++i)
    {
      this->insert(EdgeKey(v[i % n], v[(i + 1) % n]));
    }
  }

  /// An face contains a bunch of vertices
  FaceKey(uint n, uint const * v) :
    std::set<EdgeKey>(),
    idx(std::rand())
  {
    for (uint i = 0; i < n; ++i)
    {
      this->insert(EdgeKey(v[i % n], v[(i + 1) % n]));
    }
  }

  /// Construct a key from face vertices
  inline void set(Face const& f)
  {
    this->clear();
    uint const * v = f.entities(0);
    uint const n = f.num_entities(0);
    for (uint i = 0; i < n; ++i)
    {
      this->insert(EdgeKey(v[i % n], v[(i + 1) % n]));
    }
    idx = std::rand();
  }

  /// Construct a key from face vertices
  inline void set(uint n, uint const * v)
  {
    this->clear();
    for (uint i = 0; i < n; ++i)
    {
      this->insert(EdgeKey(v[i % n], v[(i + 1) % n]));
    }
    idx = std::rand();
  }

  ///
  inline size_t hash() const
  {
    size_t ret = 0;
    for (std::set<EdgeKey>::const_iterator it = this->begin();
         it != this->end(); ++it)
    {
      ret = ret ^ it->hash();
    }
    return ret;
  }

  ///
  uint idx;

};

} /* namespace dolfin */

#if ENABLE_UNORDERED_CXX
#if (HAVE_TR1_UNORDERED_MAP && HAVE_TR1_UNORDERED_SET) || \
    (__IBMCPP__ && __IBMCPP_TR1__)|| \
    (ENABLE_BOOST_TR1)

namespace std
{

namespace tr1
{

template<>
struct hash<dolfin::FaceKey>
{
  inline std::size_t operator()(dolfin::FaceKey const& f) const
  {
    return f.hash();
  }
};

} /* namespace tr1 */

} /* namespace std */

#elif (HAVE_UNORDERED_MAP && HAVE_UNORDERED_SET)

namespace std
{

template<>
struct hash<dolfin::FaceKey>
{
  inline std::size_t operator()(dolfin::FaceKey const& e) const
  {
    return e.hash();
  }
};

} /* namespace std */

#endif
#endif

#endif /* __DOLFIN_MESH_FACE_KEY */
