#include <dolfin/fem/DofMapCache.h>
#include <dolfin/fem/DofMap.h>

#include <iomanip>

#include <ufc.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
DofMapCache * const DofMapCache::instance_ = new DofMapCache();

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
		warning(
				"DofMapCache is not empty, consider checking if the dof maps have been properly released");
	}
	message("Cleanup DofMapCache");
	for (dofmap_container_t::iterator it = cache_.begin(); it != cache_.end(); ++it)
	{
		std::stringstream msg;
		msg << "Delete DofMap with hash : \'" << it->first
				<< "\'\n\thas still a count of " << it->second.count
				<< std::endl;
		warning(msg.str());
		delete it->second.dofmap;
	}
	cache_.clear();
}

//-----------------------------------------------------------------------------
DofMap * DofMapCache::acquire_dofmap(ufc::form const& form, uint const& i,
										Mesh& mesh)
{
	message("Acquire DofMap for Form coefficient %i", i);
	DofMap * ret = NULL;
	// Create UFC dof map
	ufc::dof_map * ufc_dof_map = form.create_dof_map(i);
	dolfin_assert(ufc_dof_map);

	std::string const h = DofMap::make_hash(*ufc_dof_map, mesh);
	dofmap_container_t::iterator it = cache_.find(h);
	message("Hash in cache is : %s", h.c_str());
	if (it == cache_.end())
	{
		message(0, "Creating dof map (not in cache): %s",
				ufc_dof_map->signature());


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
			error(
					"DofMap pointer has already been registered with another hash");
		}
	}
	else
	{
		message(0, "Reusing dof map (already in cache): %s",
				ufc_dof_map->signature());

		ret = it->second.dofmap;
		it->second.count++;


		// Delete UFC dof map (not used)
		delete ufc_dof_map;
	}
	return ret;
}

//-----------------------------------------------------------------------------
DofMap * DofMapCache::acquire_dofmap(std::string const& dofmap_signature,
										Mesh& mesh)
{
	message("Acquire DofMap from signature");
	DofMap * ret = NULL;

	std::string const h = DofMap::make_hash(dofmap_signature, mesh);
	message("Hash in cache is : %s", h.c_str());
	dofmap_container_t::iterator it = cache_.find(h);
	if (it == cache_.end())
	{
		message(0, "Creating dof map (not in cache): %s",
				dofmap_signature.c_str());


		// Create DOLFIN dof map
		ret = new DofMap(dofmap_signature, mesh);
		dolfin_assert(ret);

		std::string const dm_h = ret->hash();
		if (h != dm_h)
		{
			info();
			std::stringstream ss;
			ss << std::endl << "DofMap@" << ret << std::endl
					<< "Signature to be inserted  : " << h << std::endl
					<< "Signature of DofMap       : " << dm_h
					<< std::endl
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
			error(
					"DofMap pointer has already been registered with another hash");
		}
	}
	else
	{
		message(0, "Reusing dof map (already in cache): %s",
				dofmap_signature.c_str());

		ret = it->second.dofmap;
		it->second.count++;
	}
	return ret;
}

//-----------------------------------------------------------------------------
void DofMapCache::release_dofmap(DofMap& dof_map)
{
	message("Release DofMap");
	std::string const h = dof_map.hash();
	message("Hash to be released : %s", h.c_str());
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
			info();
			std::stringstream ss;
			ss << std::endl << "DofMap@" << &dof_map << std::endl
					<< "Signature to be released     : " << h << std::endl
					<< "Signature already registered : " << expected_h
					<< std::endl
					<< "DofMap object refers to two different signatures";
			error(ss.str());
		}
		dofmap_container_t::iterator dm_it = cache_.find(expected_h);
		dofmap_token_t& dm_token = dm_it->second;
		dm_token.count--;
		if (dm_token.count == 0)
		{
			message("Cleanup unused DofMap");
			delete &dof_map;
			rlist_.erase(it);
			cache_.erase(dm_it);
		}
	}
	message("State of cache after releasing DofMap");
	info();
}

//-----------------------------------------------------------------------------
void DofMapCache::info() const
{
	message("Number of DofMaps in cache : %i", cache_.size());
	for (dofmap_container_t::const_iterator it = cache_.begin();
			it != cache_.end(); ++it)
	{
		std::cout << std::setw(128) << it->first << " : @" << it->second.dofmap
				<< " : " << it->second.count << std::endl;
	}
	for (dofmap_rlist_t::const_iterator it = rlist_.begin(); it != rlist_.end();
			++it)
	{
		std::cout << "@" << it->first << " : " << it->second << std::endl;
	}
}

}

