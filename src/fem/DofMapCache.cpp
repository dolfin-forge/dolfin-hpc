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
	for (dofmap_container_t::iterator it = cache_.begin(); it != cache_.end();
			++it)
	{
		std::stringstream msg;
		msg << "DofMap with hash : \'" << it->first
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
	message("Acquire DofMap");
	std::cout << "Form @" << &form << " id = " << i << std::endl;
	DofMap * ret = NULL;
	// Create UFC dof map
	message("Create UFC dof map");
	std::cout << "Size of cache = " << cache_.size() << std::endl;
	ufc::dof_map * ufc_dof_map = form.create_dof_map(i);
	dolfin_assert(ufc_dof_map);

	std::string const h = DofMap::make_hash(*ufc_dof_map, mesh);
	info();
	dofmap_container_t::iterator it = cache_.find(h);
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
	message("Acquire DofMap");
	DofMap * ret = NULL;
	// Create UFC dof map
	message("Create UFC dof map");
	std::cout << "Size of cache = " << cache_.size() << std::endl;

	std::string const h = DofMap::make_hash(dofmap_signature, mesh);
	info();
	dofmap_container_t::iterator it = cache_.find(h);
	if (it == cache_.end())
	{
		message(0, "Creating dof map (not in cache): %s",
				dofmap_signature.c_str());

		// Create DOLFIN dof map
		ret = new DofMap(dofmap_signature, mesh);
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
	info();
	std::string h = dof_map.hash();
	dofmap_rlist_t::iterator it = rlist_.find(&dof_map);
	if (it == rlist_.end())
	{
		error("Trying to release inexistent DofMap");
	}
	else
	{
		if (h != it->second)
		{
			std::stringstream ss;
			ss << std::endl << "Signature to be released     : " << h
					<< std::endl << "Signature already registered : "
					<< it->second << std::endl
					<< "DofMap object refers to two different signatures";
			error(ss.str());
		}
		dofmap_container_t::iterator dm_it = cache_.find(h);
		dofmap_token_t& dm_token = dm_it->second;
		dm_token.count--;
		if (dm_token.count == 0)
		{
			message("Cleanup unused DofMap");
			rlist_.erase(rlist_.find(dm_token.dofmap));
			delete dm_token.dofmap;
			cache_.erase(dm_it);
		}
	}
	message("After releasing DofMap");
	info();
}

//-----------------------------------------------------------------------------
void DofMapCache::info() const
{
	std::cout << "Number of DofMaps in cache : " << cache_.size() << std::endl;
	for (dofmap_container_t::const_iterator it = cache_.begin();
			it != cache_.end(); ++it)
	{
		std::cout << std::setw(128) << it->first << " : " << it->second.count
				<< std::endl;
	}
}

}

