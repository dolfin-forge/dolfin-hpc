#ifndef __DOF_MAP_CACHE_H
#define __DOF_MAP_CACHE_H

#include <dolfin/common/types.h>
#include <dolfin/mesh/Mesh.h>

#include <ufc.h>

#include <map>
#include <string>

namespace dolfin
{

class DofMap;

class DofMapCache
{

	struct dofmap_token_t
	{
		uint count;
		DofMap * const dofmap;

		dofmap_token_t(DofMap * dm) :
				count(1), dofmap(dm)
		{
		}
	};

#if __SUNPRO_CC
	typedef std::map<std::string, dofmap_token_t > dofmap_container_t;
	typedef std::pair<std::string, dofmap_token_t> dofmap_item_t;
	typedef std::map<DofMap *, std::string > dofmap_rlist_t;
	typedef std::pair<DofMap *, std::string> dofmap_ritem_t;
#else
	typedef std::map<std::string const, dofmap_token_t> dofmap_container_t;
	typedef std::pair<std::string const, dofmap_token_t> dofmap_item_t;
	typedef std::map<DofMap *, std::string const> dofmap_rlist_t;
	typedef std::pair<DofMap *, std::string const> dofmap_ritem_t;
#endif

public:

	static DofMapCache& instance()
	{ return *instance_; }

	DofMap * acquire_dofmap(ufc::form const& form, uint const& i, Mesh& mesh);
	DofMap * acquire_dofmap(std::string const& dofmap_signature, Mesh& mesh);
	void release_dofmap(DofMap& dof_map);

	void info() const;

private:

	DofMapCache();
	~DofMapCache();

	void __request_dofmap(std::string const);

	dofmap_container_t cache_;
	dofmap_rlist_t rlist_;

	static DofMapCache * const instance_;

};

}

#endif /* __DOF_MAP_CACHE_H */


