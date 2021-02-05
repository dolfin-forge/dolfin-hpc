// Copyright (C) 2013-2014 Aurélien Larcher
// Licensed under the GNU LGPL Version 2.1.
//

#include <dolfin/fem/DofMap.h>
#include <dolfin/fem/DofMapCache.h>
#include <dolfin/fem/Form.h>
#include <dolfin/function/Function.h>

#include <iomanip>

namespace dolfin
{

//-----------------------------------------------------------------------------

DofMapCache::DofMapCache()
{
  cache_.clear();
  rlist_.clear();
}

//-----------------------------------------------------------------------------

DofMapCache::~DofMapCache()
{
  if ( cache_.size() != 0 )
  {
    warning( "DofMapCache is not empty: "
             "some dof maps have not been properly released" );
  }
  for ( container_t::value_type & token : cache_ )
  {
    delete token.second.item;
  }
  cache_.clear();
}

//-----------------------------------------------------------------------------

auto DofMapCache::acquire( Mesh & mesh, Form const & form, size_t const & i )
  -> DofMap &
{
  // Update dof maps of form:
  // This triggers creation of a dof map set if needed and acquistion of a token
  // for each coefficient of the form.
  // If the dof map set has been created the call just return without doing
  // anything.
  form.update_dofmaps();

  // Create UFC dof map for the i-th coefficient to get the signature
  ufc::dofmap * ufc_dofmap = form().create_dofmap( i );

  // Do not transfer ownership to the possibly created dofmap
  // The UFC dofmap is cloned as member attribute of the DOLFIN dofmap
  DofMap * ret = &acquire( mesh, *ufc_dofmap );

  // Delete UFC dof map as it is not longer used
  delete ufc_dofmap;

  return *ret;
}

//-----------------------------------------------------------------------------

auto DofMapCache::acquire( Mesh & mesh, ufc::dofmap const & dofmap )
  -> DofMap &
{
  DofMap *              ret = nullptr;
  std::string const     h   = DofMap::make_hash( mesh, dofmap );
  container_t::iterator it  = cache_.find( h );

  if ( it == cache_.end() )
  {
    // Create DOLFIN dofmap and insert in map
    ret = new DofMap( mesh, dofmap );
    cache_.insert( item_t( h, token_t( ret ) ) );

    if ( rlist_.find( ret ) == rlist_.end() )
    {
      rlist_.insert( ritem_t( ret, h ) );
    }
    else
    {
      error( "DofMap pointer has already been registered with another hash" );
    }
  }
  else
  {
    // The dofmap has already been created
    ret = it->second.item;
    dolfin_assert( ret );

    // Check hash consistency
    std::string const dm_h = ret->hash();
    if ( h != dm_h )
    {
      disp();
      std::stringstream ss;
      ss << std::endl
         << "DofMap@" << ret << std::endl
         << "Signature to be inserted  : " << h << std::endl
         << "Signature of DofMap       : " << dm_h << std::endl
         << "DofMap object refers to two different signatures";
      error( ss.str() );
    }
    it->second.count++;
  }
  return *ret;
}

//-----------------------------------------------------------------------------

void DofMapCache::release( DofMap & dofmap )
{
  std::string const h          = dofmap.hash();
  rlist_t::iterator it         = rlist_.find( &dofmap );
  std::string const expected_h = it->second;

  if ( it == rlist_.end() )
  {
    error( "Trying to release inexistent DofMap" );
  }
  else
  {
    if ( h != expected_h )
    {
      disp();
      std::stringstream ss;
      ss << std::endl
         << "DofMap@" << &dofmap << std::endl
         << "Signature to be released     : " << h << std::endl
         << "Signature already registered : " << expected_h << std::endl
         << "DofMap object refers to two different signatures";
      error( ss.str() );
    }
    container_t::iterator dm_it    = cache_.find( expected_h );
    token_t &             dm_token = dm_it->second;
    dm_token.count--;
    // message("Release dofmap: %s",h.c_str());
    if ( dm_token.count == 0 )
    {
      delete dm_token.item;
      rlist_.erase( it );
      cache_.erase( dm_it );
    }
  }
}

//-----------------------------------------------------------------------------

void DofMapCache::disp() const
{
  message( "Number of DofMaps in cache : %i", cache_.size() );
  message( "Cache :" );

  for ( container_t::value_type const & token : cache_ )
  {
    std::cout << std::setw( 128 ) << token.first << " : @" << token.second.item
              << " : " << token.second.count << std::endl;
  }
  message( "Reverse list :" );
  for ( rlist_t::value_type const & rlist : rlist_ )
  {
    std::cout << "@" << rlist.first << " : " << rlist.second << std::endl;
  }
}

//-----------------------------------------------------------------------------

} // namespace dolfin
