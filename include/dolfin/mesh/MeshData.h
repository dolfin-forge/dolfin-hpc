// Copyright (C) 2017 Aurelien Larcher.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_MESH_DATA_H
#define __DOLFIN_MESH_DATA_H

#include <dolfin/mesh/MeshValues.h>

namespace dolfin
{

/**
 *  @class  MeshData
 *
 *  @brief  A container to hold mesh functions, for instance to define those
 *          to be transferred with load balancing.
 *
 */

//-----------------------------------------------------------------------------
class MeshData
{

  //--- Naive type erasure ----------------------------------------------------

  enum Values { Vbool, Vuint, Vreal };
  template<class V> static inline Values values();

  enum Entity { Evertex, Eedge, Eface, Efacet, Ecell };
  template<class E> static inline Entity entity();

  using key = std::pair<uint, uint>;

  template <class V, class E>
  struct type : public key { type() : key(values<V>(), entity<E>()) {} };

  struct data_array { virtual uint size() const = 0; };

  template <class V, class E>
  struct array : public data_array
  {
    using type = Array<MeshValues<V, E> *>; type data;
    uint size() const { return data.size(); };
  };

  //---------------------------------------------------------------------------

public:

  MeshData(Mesh& mesh) : M_(mesh) {}

  //
  Mesh& mesh() { return M_; }

  // Add mesh function to data, error if exists
  template <class V, class E> void add(MeshValues<V, E>& function)
  {
    if (&M_ != &function.mesh())
    {
      error("MeshData : adding mesh function from a different mesh");
    }
    typename array<V, E>::type& A = lookup<V, E>();
    if(std::find(A.begin(), A.end(), &function) != A.end())
    {
      error("MeshData : duplicate add of mesh function");
    }
    A.push_back(&function);
  }

  // Add mesh function to data if it does not exist
  template <class V, class E> void insert(MeshValues<V, E>& function)
  {
    if (&M_ != &function.mesh())
    {
      error("MeshData : adding mesh function from a different mesh");
    }
    typename array<V, E>::type& A = lookup<V, E>();
    if(std::find(A.begin(), A.end(), &function) == A.end())
    {
      A.push_back(&function);
    }
  }

  // Return number of function of given type in data
  template <class V, class E> uint size()
  {
    return lookup<V, E>().size();
  }

  // Return total number of functions
  uint size() const
  {
    uint s = 0;
    for (Store::const_iterator it = S_.begin(); it != S_.end(); ++it)
    {
      s += it->second->size();
    }
    return s;
  }

  // Count number of function instances in data
  template <class V, class E> uint count(MeshValues<V, E>& function)
  {
    if (&M_ != &function.mesh())
    {
      error("MeshData : looking for a mesh function from a different mesh");
    }
    typename array<V, E>::type& A = lookup<V, E>();
    return (std::find(A.begin(), A.end(), &function) != A.end());
  }

  //--- ITERATOR --------------------------------------------------------------

  template <class V, class E>
  struct iterator
  {

    iterator(MeshData& D) : b_(D.begin<V, E>()), e_(D.end<V, E>()), i_(b_) { };

    //
    bool valid() { return (i_ != e_); }

    //
    uint pos() { return (i_ - b_); }

    //
    iterator& operator++() { ++i_; return *this; }

    //
    MeshValues<V, E>* operator->() { return *i_; }

    //
    MeshValues<V, E>& operator*() { return **i_; }


  private:

    typename array<V, E>::type::iterator b_;
    typename array<V, E>::type::iterator e_;
    typename array<V, E>::type::iterator i_;

  };

  //---------------------------------------------------------------------------

private:

  template <class V, class E> typename array<V, E>::type& lookup()
  {
    type<V, E> K;
    std::pair<Store::iterator, bool> it(S_.find(K), false);
    if(it.first == S_.end()) { it = S_.insert(Item(K, new array<V, E>())); }
    return (static_cast<array<V, E> *>(it.first->second)->data);
  }

  template <class V, class E> typename array<V, E>::type::iterator begin()
  {
    return lookup<V, E>().begin();
  }

  template <class V, class E> typename array<V, E>::type::iterator end()
  {
    return lookup<V, E>().end();
  }

  using Store = _ordered_map<key, data_array *>;
  using Item  = std::pair<key, data_array *>;

  Mesh& M_;
  Store S_;

};

template<> inline MeshData::Values MeshData::values<bool>() { return MeshData::Vbool; }
template<> inline MeshData::Values MeshData::values<uint>() { return MeshData::Vuint; }
template<> inline MeshData::Values MeshData::values<real>() { return MeshData::Vreal; }

template<> inline MeshData::Entity MeshData::entity<Vertex>() { return MeshData::Evertex; }
template<> inline MeshData::Entity MeshData::entity<Edge>()   { return MeshData::Eedge;   }
template<> inline MeshData::Entity MeshData::entity<Face>()   { return MeshData::Eface;   }
template<> inline MeshData::Entity MeshData::entity<Facet>()  { return MeshData::Efacet;  }
template<> inline MeshData::Entity MeshData::entity<Cell>()   { return MeshData::Ecell;   }

//-----------------------------------------------------------------------------

} /* namespace dolfin */

#endif /* __DOLFIN_MESH_DATA_H */

