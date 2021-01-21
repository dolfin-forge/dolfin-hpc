
#ifndef __DOLFIN_MESH_ENTITY_KEY_H
#define __DOLFIN_MESH_ENTITY_KEY_H

#include <dolfin/common/types.h>
#include <dolfin/log/log.h>

#include <cstdlib>
#include <string>
#include <algorithm>

namespace dolfin
{

/*
 * Use index array as key.
 * The tag idx is to carry information and not used for comparison.
 *
 */

struct EntityKey
{

  ///
  EntityKey(uint D) :
      size(D),
      indices(D ? new uint[size]() : nullptr),
      idx(0)
  {
  }

  ///
  EntityKey(uint D, uint const * v) :
      size(D),
      indices(D ? new uint[size] : nullptr),
      idx(0)
  {
    set(v);
  }

  ///
  EntityKey(uint D, uint const * v, uint i) :
      size(D),
      indices(D ? new uint[size] : nullptr),
      idx(i)
  {
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
  auto operator=(EntityKey const& other) -> EntityKey&
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
  inline auto operator<(EntityKey const& other) const -> bool
  {
    dolfin_assert(size == other.size);
    for (uint i = 0; i < size; ++i)
    {
      if (this->indices[i] != other.indices[i])
      {
        return this->indices[i] < other.indices[i];
      }
    }
    return false;
  }

  ///
  inline auto operator<=(EntityKey const& other) const -> bool
  {
    dolfin_assert(size == other.size);
    for (uint i = 0; i < size; ++i)
    {
      if (this->indices[i] != other.indices[i])
      {
        return this->indices[i] < other.indices[i];
      }
    }
    return true;
  }

  ///
  inline auto operator==(EntityKey const& other) const -> bool
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
  inline auto operator!=(EntityKey const& other) const -> bool
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
  inline auto operator>=(EntityKey const& other) const -> bool
  {
    dolfin_assert(size == other.size);
    for (uint i = 0; i < size; ++i)
    {
      if (this->indices[i] != other.indices[i])
      {
        return this->indices[i] > other.indices[i];
      }
    }
    return true;
  }

  ///
  inline auto operator>(EntityKey const& other) const -> bool
  {
    dolfin_assert(size == other.size);
    for (uint i = 0; i < size; ++i)
    {
      if (this->indices[i] != other.indices[i])
      {
        return this->indices[i] > other.indices[i];
      }
    }
    return false;
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
  inline auto hash() const -> size_t
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

  uint const size{0};
  uint * indices{nullptr};
  uint idx{0};

private:

  ///
  EntityKey() = delete;

};

} /* namespace dolfin */

namespace std
{

template<>
struct hash<dolfin::EntityKey>
{
  inline auto operator()(dolfin::EntityKey const& e) const -> std::size_t
  {
    return e.hash();
  }
};

} /* namespace std */

#endif /* __DOLFIN_MESH_ENTITY_KEY_H */
