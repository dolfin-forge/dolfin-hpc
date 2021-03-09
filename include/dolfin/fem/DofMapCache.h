// Copyright (C) 2013-2014 Aurélien Larcher
// Licensed under the GNU LGPL Version 2.1.
//

#ifndef __DOLFIN_DOF_MAP_CACHE_H
#define __DOLFIN_DOF_MAP_CACHE_H

#include <dolfin/common/types.h>
#include <dolfin/mesh/Mesh.h>

#include <string>

namespace ufc
{

class dofmap;

} // namespace ufc

namespace dolfin
{

class DofMap;
class Form;
class Function;

//-----------------------------------------------------------------------------

class DofMapCache
{

  friend class DofMap;

  struct token_t
  {
    size_t         count;
    DofMap * const item;

    token_t( DofMap * dm )
      : count( 1 )
      , item( dm )
    {
    }
  };

#if __SUNPRO_CC
  using container_t = _ordered_map< std::string, token_t >;
  using item_t      = std::pair< std::string, token_t >;
  using rlist_t     = _ordered_map< DofMap *, std::string >;
  using ritem_t     = std::pair< DofMap *, std::string >;
#else
  using container_t = _ordered_map< std::string const, token_t >;
  using item_t      = std::pair< std::string const, token_t >;
  using rlist_t     = _ordered_map< DofMap *, std::string const >;
  using ritem_t     = std::pair< DofMap *, std::string const >;
#endif

public:
  ///
  void disp() const;

protected:
  /// Meyers singleton
  static auto instance() -> DofMapCache &
  {
    static DofMapCache instance_;
    return instance_;
  }

  /// Return a dofmap corresponding to the i-th coefficient space
  auto acquire( Mesh & mesh, Form const & form, size_t const & i ) -> DofMap &;

  /// Return a dofmap for the UFC dofmap object
  auto acquire( Mesh & mesh, ufc::dofmap const & dofmap, bool owner ) -> DofMap &;

  ///
  void release( DofMap & dofmap );

private:
  ///
  DofMapCache();

  ///
  ~DofMapCache();

  container_t cache_;
  rlist_t     rlist_;
};

//-----------------------------------------------------------------------------

} // namespace dolfin

#endif /* __DOLFIN_DOF_MAP_CACHE_H */
