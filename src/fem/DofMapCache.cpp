// Copyright (C) 2013-2014 Aurélien Larcher
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2013-06-11
// Last changed: 2014-02-04

#include <dolfin/fem/DofMapCache.h>
#include <dolfin/fem/DofMap.h>
#include <dolfin/fem/Form.h>
#include <dolfin/function/Function.h>

#include <iomanip>

#include <ufc.h>

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
  if (cache_.size() != 0)
  {
    warning("DofMapCache is not empty: "
            "consider checking if the dof maps have been properly released");
  }
  for (dofmap_container_t::iterator it = cache_.begin(); it != cache_.end();
      ++it)
  {
    delete it->second.dofmap;
  }
  cache_.clear();
}

//-----------------------------------------------------------------------------
DofMap& DofMapCache::acquire_dofmap(Mesh& mesh, Form const& form, uint const& i)
{
  message(2, "Acquire DofMap for Form coefficient %i", i);

  // Update dof maps of form:
  // This triggers creation of a dof map set if needed and acquistion of a token
  // for each coefficient of the form.
  // If the dof map set has been created the call just return without doing
  // anything.
  form.update_dofmaps(mesh);

  DofMap * ret = NULL;
  // Create UFC dof map
  ufc::dof_map * ufc_dof_map = form.form().create_dof_map(i);
  dolfin_assert(ufc_dof_map);

  std::string const h = DofMap::make_hash(*ufc_dof_map, mesh);
  dofmap_container_t::iterator it = cache_.find(h);

  if (it == cache_.end())
  {

    // Create DOLFIN dof map
    ret = new DofMap(*ufc_dof_map, mesh, true);
    dolfin_assert(ret);

    cache_.insert(dofmap_item_t(h, dofmap_token_t(ret)));

    if (rlist_.find(ret) == rlist_.end())
    {
      rlist_.insert(dofmap_ritem_t(ret, h));
    }
    else
    {
      error("DofMap pointer has already been registered with another hash");
    }
  }
  else
  {
    ret = it->second.dofmap;
    it->second.count++;

    // Delete UFC dof map (not used)
    delete ufc_dof_map;
  }

  return *ret;
}

//-----------------------------------------------------------------------------
DofMap& DofMapCache::acquire_dofmap(Mesh& mesh,
                                    std::string const& dofmap_signature)
{
  message(2, "Acquire DofMap from signature: %s", dofmap_signature.c_str());
  DofMap * ret = NULL;

  std::string const h = DofMap::make_hash(dofmap_signature, mesh);
  dofmap_container_t::iterator it = cache_.find(h);

  if (it == cache_.end())
  {

    // Create DOLFIN dof map
    ret = new DofMap(dofmap_signature, mesh);
    dolfin_assert(ret);

    std::string const dm_h = ret->hash();
    if (h != dm_h)
    {
      disp();
      std::stringstream ss;
      ss << std::endl << "DofMap@" << ret << std::endl
          << "Signature to be inserted  : " << h << std::endl
          << "Signature of DofMap       : " << dm_h << std::endl
          << "DofMap object refers to two different signatures";
      error(ss.str());
    }

    cache_.insert(dofmap_item_t(h, dofmap_token_t(ret)));

    if (rlist_.find(ret) == rlist_.end())
    {
      rlist_.insert(dofmap_ritem_t(ret, h));
    }
    else
    {
      error("DofMap pointer has already been registered with another hash");
    }
  }
  else
  {
    ret = it->second.dofmap;
    it->second.count++;
  }
  return *ret;
}

//-----------------------------------------------------------------------------
DofMap& DofMapCache::acquire_dofmap(Function& f)
{
  if (f.type() != Function::discrete)
  {
    error("Can only acquire DofMap from a discrete function.");
  }
  return acquire_dofmap(f.mesh(), DofMap::dofmap_signature(f.signature()));
}

//-----------------------------------------------------------------------------
void DofMapCache::release_dofmap(DofMap& dof_map)
{
  std::string const h = dof_map.hash();
  dofmap_rlist_t::iterator it = rlist_.find(&dof_map);
  std::string const expected_h = it->second;

  if (it == rlist_.end())
  {
    error("Trying to release inexistent DofMap");
  }
  else
  {
    if (h != expected_h)
    {
      disp();
      std::stringstream ss;
      ss << std::endl << "DofMap@" << &dof_map << std::endl
          << "Signature to be released     : " << h << std::endl
          << "Signature already registered : " << expected_h << std::endl
          << "DofMap object refers to two different signatures";
      error(ss.str());
    }
    dofmap_container_t::iterator dm_it = cache_.find(expected_h);
    dofmap_token_t& dm_token = dm_it->second;
    dm_token.count--;
    if (dm_token.count == 0)
    {
      delete &dof_map;
      rlist_.erase(it);
      cache_.erase(dm_it);
    }
  }
}

//-----------------------------------------------------------------------------
void DofMapCache::disp() const
{
  message("Number of DofMaps in cache : %i", cache_.size());
  message("Cache :");
  for (dofmap_container_t::const_iterator it = cache_.begin();
      it != cache_.end(); ++it)
  {
    std::cout << std::setw(128) << it->first << " : @" << it->second.dofmap
        << " : " << it->second.count << std::endl;
  }
  message("Reverse list :");
  for (dofmap_rlist_t::const_iterator it = rlist_.begin(); it != rlist_.end();
      ++it)
  {
    std::cout << "@" << it->first << " : " << it->second << std::endl;
  }
}

}

