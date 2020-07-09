// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_FUNCTION_SPACE_MAP_H
#define __DOLFIN_FUNCTION_SPACE_MAP_H

#include <dolfin/common/types.h>
#include <dolfin/log/log.h>
#include <dolfin/ufl/UFLrepr.h>
#include <dolfin/ufl/UFLFiniteElementSpace.h>

#include <string>

namespace dolfin
{

/**
 *  @class  FunctionSpaceMap
 *
 *  @brief  This class serves as provider of function space signature to allow
 *          specifying the type of discrete space to instantiate from generated
 *          source code.
 */

class FunctionSpaceMap
{
  typedef _ordered_map< std::string, ufl::repr > Container;
  typedef std::pair< std::string, ufl::repr >    Item;

public:
  ///
  FunctionSpaceMap();

  ///
  ~FunctionSpaceMap();

  ///
  bool has( std::string const & label ) const;

  ///
  ufl::repr get( std::string const & label ) const;

  ///
  uint size() const;

  ///
  void add( std::string const & label, ufl::FiniteElementSpace const & space );

  ///
  void clear();

private:
  Container spaces_;
};

//-----------------------------------------------------------------------------
inline bool FunctionSpaceMap::has( std::string const & label ) const
{
  return ( spaces_.count( label ) > 0 );
}

//-----------------------------------------------------------------------------
inline ufl::repr FunctionSpaceMap::get( std::string const & label ) const
{
  Container::const_iterator it = spaces_.find( label );
  if ( spaces_.find( label ) == spaces_.end() )
  {
    error( "Space label '%s' not defined in function space map.",
           label.c_str() );
  }
  return it->second;
}

//-----------------------------------------------------------------------------
inline uint FunctionSpaceMap::size() const
{
  return spaces_.size();
}

//-----------------------------------------------------------------------------
inline void FunctionSpaceMap::add( std::string const &             label,
                                   ufl::FiniteElementSpace const & space )
{
  if ( this->has( label ) )
  {
    error( "Space label '%s' already defined in function space map.",
           label.c_str() );
  }
  spaces_.insert( Item( label, space.repr() ) );
}

//-----------------------------------------------------------------------------
inline void FunctionSpaceMap::clear()
{
  spaces_.clear();
}

}

#endif /* __DOLFIN_FUNCTION_SPACE_MAP_H */
