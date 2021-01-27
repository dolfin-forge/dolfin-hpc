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

  enum Values
  {
    Vbool,
    Vuint,
    Vsizet,
    Vreal
  };
  template < class V >
  static inline auto values() -> Values;

  enum Entity
  {
    Evertex,
    Eedge,
    Eface,
    Efacet,
    Ecell
  };
  template < class E >
  static inline auto entity() -> Entity;

  using key = std::pair< uint, uint >;

  template < class V, class E >
  struct type : public key
  {
    type()
      : key( values< V >(), entity< E >() )
    {
    }
  };

  struct data_array
  {
    virtual auto size() const -> uint = 0;
  };

  template < class V, class E >
  struct array : public data_array
  {
    using type = std::vector< MeshValues< V, E > * >;
    type data;
    auto size() const -> uint override
    {
      return data.size();
    };
  };

  //---------------------------------------------------------------------------

public:
  MeshData( Mesh & mesh )
    : M_( mesh )
  {
  }

  //
  auto mesh() -> Mesh &
  {
    return M_;
  }

  // Add mesh function to data, error if exists
  template < class V, class E >
  void add( MeshValues< V, E > & function )
  {
    if ( &M_ != &function.mesh() )
    {
      error( "MeshData : adding mesh function from a different mesh" );
    }
    typename array< V, E >::type & A = lookup< V, E >();
    if ( std::find( A.begin(), A.end(), &function ) != A.end() )
    {
      error( "MeshData : duplicate add of mesh function" );
    }
    A.push_back( &function );
  }

  // Add mesh function to data if it does not exist
  template < class V, class E >
  void insert( MeshValues< V, E > & function )
  {
    if ( &M_ != &function.mesh() )
    {
      error( "MeshData : adding mesh function from a different mesh" );
    }
    typename array< V, E >::type & A = lookup< V, E >();
    if ( std::find( A.begin(), A.end(), &function ) == A.end() )
    {
      A.push_back( &function );
    }
  }

  // Return number of function of given type in data
  template < class V, class E >
  auto size() -> uint
  {
    return lookup< V, E >().size();
  }

  // Return total number of functions
  auto size() const -> uint
  {
    uint s = 0;
    for ( Store::const_iterator it = S_.begin(); it != S_.end(); ++it )
    {
      s += it->second->size();
    }
    return s;
  }

  // Count number of function instances in data
  template < class V, class E >
  auto count( MeshValues< V, E > & function ) -> uint
  {
    if ( &M_ != &function.mesh() )
    {
      error( "MeshData : looking for a mesh function from a different mesh" );
    }
    typename array< V, E >::type & A = lookup< V, E >();
    return ( std::find( A.begin(), A.end(), &function ) != A.end() );
  }

  //--- ITERATOR --------------------------------------------------------------

  template < class V, class E >
  struct iterator
  {

    iterator( MeshData & D )
      : b_( D.begin< V, E >() )
      , e_( D.end< V, E >() )
      , i_( b_ ) {};

    //
    auto valid() -> bool
    {
      return ( i_ != e_ );
    }

    //
    auto pos() -> uint
    {
      return ( i_ - b_ );
    }

    //
    auto operator++() -> iterator &
    {
      ++i_;
      return *this;
    }

    //
    auto operator->() -> MeshValues< V, E > *
    {
      return *i_;
    }

    //
    auto operator*() -> MeshValues< V, E > &
    {
      return **i_;
    }

  private:
    typename array< V, E >::type::iterator b_;
    typename array< V, E >::type::iterator e_;
    typename array< V, E >::type::iterator i_;
  };

  //---------------------------------------------------------------------------

private:
  template < class V, class E >
  auto lookup() -> typename array< V, E >::type &
  {
    type< V, E >                       K;
    std::pair< Store::iterator, bool > it( S_.find( K ), false );
    if ( it.first == S_.end() )
    {
      it = S_.insert( Item( K, new array< V, E >() ) );
    }
    return ( static_cast< array< V, E > * >( it.first->second )->data );
  }

  template < class V, class E >
  auto begin() -> typename array< V, E >::type::iterator
  {
    return lookup< V, E >().begin();
  }

  template < class V, class E >
  auto end() -> typename array< V, E >::type::iterator
  {
    return lookup< V, E >().end();
  }

  using Store = _ordered_map< key, data_array * >;
  using Item  = std::pair< key, data_array * >;

  Mesh & M_;
  Store  S_;
};

template <>
inline auto MeshData::values< bool >() -> MeshData::Values
{
  return MeshData::Vbool;
}
template <>
inline auto MeshData::values< uint >() -> MeshData::Values
{
  return MeshData::Vuint;
}
template <>
inline auto MeshData::values< size_t >() -> MeshData::Values
{
  return MeshData::Vsizet;
}
template <>
inline auto MeshData::values< real >() -> MeshData::Values
{
  return MeshData::Vreal;
}

template <>
inline auto MeshData::entity< Vertex >() -> MeshData::Entity
{
  return MeshData::Evertex;
}
template <>
inline auto MeshData::entity< Edge >() -> MeshData::Entity
{
  return MeshData::Eedge;
}
template <>
inline auto MeshData::entity< Face >() -> MeshData::Entity
{
  return MeshData::Eface;
}
template <>
inline auto MeshData::entity< Facet >() -> MeshData::Entity
{
  return MeshData::Efacet;
}
template <>
inline auto MeshData::entity< Cell >() -> MeshData::Entity
{
  return MeshData::Ecell;
}

//-----------------------------------------------------------------------------

} /* namespace dolfin */

#endif /* __DOLFIN_MESH_DATA_H */
