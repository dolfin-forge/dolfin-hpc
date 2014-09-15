// Copyright (C) 2013-2014 Aurélien Larcher
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2013-06-13 (merged from branch larcher)
// Last changed: 2014-02-04

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
class Form;
class Function;

class DofMapCache
{

  friend class DofMap;

  struct token_t
  {
    uint count;
    DofMap * const item;

    token_t(DofMap * dm) :
        count(1),
        item(dm)
    {
    }
  };

#if __SUNPRO_CC
  typedef std::map<std::string, token_t > container_t;
  typedef std::pair<std::string, token_t> item_t;
  typedef std::map<DofMap *, std::string > rlist_t;
  typedef std::pair<DofMap *, std::string> ritem_t;
#else
  typedef std::map<std::string const, token_t> container_t;
  typedef std::pair<std::string const, token_t> item_t;
  typedef std::map<DofMap *, std::string const> rlist_t;
  typedef std::pair<DofMap *, std::string const> ritem_t;
#endif

public:

  ///
  void disp() const;

protected:

  /// Meyers singleton
  static DofMapCache& instance()
  {
    static DofMapCache instance_;
    return instance_;
  }

  /// Return a dofmap corresponding to the i-th coefficient space
  DofMap& acquire(Mesh& mesh, Form const& form, uint const& i);

  /// Return a dofmap for the UFC dofmap object
  DofMap& acquire(Mesh& mesh, ufc::dofmap& dofmap, bool owner);

  ///
  void release(DofMap& dofmap);

private:

  ///
  DofMapCache();

  ///
  ~DofMapCache();

  container_t cache_;
  rlist_t rlist_;

};

}

#endif /* __DOF_MAP_CACHE_H */

