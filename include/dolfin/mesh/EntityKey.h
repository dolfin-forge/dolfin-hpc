//
//
//

#ifndef __DOLFIN_MESH_ENTITY_KEY_H
#define __DOLFIN_MESH_ENTITY_KEY_H

#include <dolfin/common/types.h>
#include <dolfin/log/log.h>

#include <cstdlib>
#include <cstring>
#include <algorithm>

namespace dolfin
{

struct EntityKey
{

  ///
  EntityKey(uint D) :
      size(D),
      indices(NULL),
      idx(0)
  {
    if (D == 0)
    {
      error("EntityKey : trying to initialize a zero sized entity");
    }
    indices = new uint[size];
    std::fill_n(indices, size, 0);
  }

  ///
  EntityKey(uint D, uint const * v) :
      size(D),
      indices(NULL),
      idx(0)
  {
    if (D == 0)
    {
      error("EntityKey : trying to initialize a zero sized entity");
    }
    indices = new uint[size];
    set(v);
  }

  ///
  EntityKey(uint D, uint const * v, uint i) :
      size(D),
      indices(NULL),
      idx(i)
  {
    if (D == 0)
    {
      error("EntityKey : trying to initialize a zero sized entity");
    }
    indices = new uint[size];
    set(v);
  }

  ///
  EntityKey(EntityKey const& other) :
      size(other.size),
      indices(new uint[size]),
      idx(other.idx)
  {
    std::copy(other.indices, other.indices + size, indices);
  }

  ///
  ~EntityKey()
  {
    delete[] indices;
  }

  ///
  EntityKey& operator=(EntityKey const& other)
  {
    if(this != &other)
    {
      if(size != other.size)
      {
        error("EntityKey : purposely disabled assignment with different size");
      }
      std::copy(other.indices, other.indices + size, indices);
    }
    return *this;
  }

  ///
  inline bool operator<(EntityKey const& other) const
  {
    dolfin_assert(size == other.size);
    for (uint i = 0; i < size; ++i)
    {
      if (this->indices[i] >= other.indices[i])
      {
        return false;
      }
    }
    return true;
  }

  ///
  inline bool operator<=(EntityKey const& other) const
  {
    dolfin_assert(size == other.size);
    for (uint i = 0; i < size; ++i)
    {
      if (this->indices[i] > other.indices[i])
      {
        return false;
      }
    }
    return true;
  }

  ///
  inline bool operator==(EntityKey const& other) const
  {
    dolfin_assert(size == other.size);
    for (uint i = 0; i < size; ++i)
    {
      if (this->indices[i] != other.indices[i])
      {
        return false;
      }
    }
    return true;
  }

  ///
  inline bool operator!=(EntityKey const& other) const
  {
    dolfin_assert(size == other.size);
    for (uint i = 0; i < size; ++i)
    {
      if (this->indices[i] != other.indices[i])
      {
        return true;
      }
    }
    return false;
  }

  ///
  inline bool operator>=(EntityKey const& other) const
  {
    dolfin_assert(size == other.size);
    for (uint i = 0; i < size; ++i)
    {
      if (this->indices[i] < other.indices[i])
      {
        return false;
      }
    }
    return true;
  }

  ///
  inline bool operator>(EntityKey const& other) const
  {
    dolfin_assert(size == other.size);
    for (uint i = 0; i < size; ++i)
    {
      if (this->indices[i] <= other.indices[i])
      {
        return false;
      }
    }
    return true;
  }

  ///
  inline void set(uint const * v)
  {
    switch(size)
    {
      case 1:
        indices[0] = v[0];
        break;
      case 2:
        //Important: copy to avoid aliasing issues
        std::copy(v, v + 2, indices);
        if(indices[1] < indices[0])
        {
          std::swap(indices[0], indices[1]);
        }
        break;
      default:
        std::copy(v, v + size, indices);
        std::sort(indices, indices + size);
        break;
    }
  }

  ///
  inline void set(uint const * v, uint i)
  {
    set(v);
    idx = i;
  }

  ///
  inline size_t hash() const
  {
    size_t ret = static_cast<std::size_t>(indices[0]);
    for (uint i = 1; i < size; ++i)
    {
      ret = ret ^ static_cast<std::size_t>(indices[i]);
    }
    return ret;
  }

  ///
  inline void disp() const
  {
    section("EntityKey");
    message("idx  : %u", idx);
    message("size : %u", size);
    for (uint i = 0; i < size; ++i)
    {
      message("i%u : %u", i, indices[i]);
    }
    end();
  }

  uint const size;
  uint * indices;
  uint idx;

private:

  ///
  EntityKey() :
      size(0),
      indices(NULL),
      idx(0)
  {
    error("EntityKey : default constructor disabled");
  }

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
struct hash<dolfin::EntityKey>
{
  inline std::size_t operator()(dolfin::EntityKey const& e) const
  {
    return e.hash();
  }
};

} /* namespace tr1 */

} /* namespace std */

#endif

#endif /* __DOLFIN_MESH_ENTITY_KEY_H */
