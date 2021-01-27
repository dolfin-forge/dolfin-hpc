
#ifndef __DOLFIN_MESH_FACE_KEY
#define __DOLFIN_MESH_FACE_KEY

#include <dolfin/mesh/EdgeKey.h>
#include <dolfin/mesh/Face.h>

namespace dolfin
{

//-----------------------------------------------------------------------------

struct FaceKey : public _ordered_set<EdgeKey>
{
  /// An face contains a bunch of vertices
  FaceKey() :
      _ordered_set<EdgeKey>(),
      idx(0)
  {
  }

  /// An face contains a bunch of vertices
  FaceKey(Face const& f) :
    _ordered_set<EdgeKey>(),
    idx(std::rand())
  {
    size_t const * v = f.entities(0);
    size_t const n = f.num_entities(0);
    for (size_t i = 0; i < n; ++i)
    {
      this->insert(EdgeKey(v[i % n], v[(i + 1) % n]));
    }
  }

  /// An face contains a bunch of vertices
  FaceKey(size_t n, size_t const * v) :
    _ordered_set<EdgeKey>(),
    idx(std::rand())
  {
    for (size_t i = 0; i < n; ++i)
    {
      this->insert(EdgeKey(v[i % n], v[(i + 1) % n]));
    }
  }

  /// Construct a key from face vertices
  inline void set(Face const& f)
  {
    this->clear();
    size_t const * v = f.entities(0);
    size_t const n = f.num_entities(0);
    for (size_t i = 0; i < n; ++i)
    {
      this->insert(EdgeKey(v[i % n], v[(i + 1) % n]));
    }
    idx = std::rand();
  }

  /// Construct a key from face vertices
  inline void set(size_t n, size_t const * v)
  {
    this->clear();
    for (size_t i = 0; i < n; ++i)
    {
      this->insert(EdgeKey(v[i % n], v[(i + 1) % n]));
    }
    idx = std::rand();
  }

  ///
  inline size_t hash() const
  {
    size_t ret = 0;
    for (_ordered_set<EdgeKey>::const_iterator it = this->begin();
         it != this->end(); ++it)
    {
      ret = ret ^ it->hash();
    }
    return ret;
  }

  ///
  size_t idx;

};

//-----------------------------------------------------------------------------

} /* namespace dolfin */

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

#endif /* __DOLFIN_MESH_FACE_KEY */
