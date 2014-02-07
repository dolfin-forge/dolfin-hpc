// Copyright (C) 2014 Bärbel Janssen.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  
// Last changed: 

#include <dolfin/ufl/UFLIndex.h>

//#include <dolfin/common/types.h>
#include <dolfin/log/log.h>

namespace ufl
{

//-----------------------------------------------------------------------------
  IndexBase::IndexBase(std::string const& name,
      dolfin::uint const& count) :
    Class(name),
    count_(count)
  {
  }

//-----------------------------------------------------------------------------
  IndexBase::IndexBase(std::string const& name,
      IndexBase const& index) :
    Class(name),
    count_(index.count())
  {
  }

//-----------------------------------------------------------------------------
  IndexBase::IndexBase(std::string const& name,
      repr_t const & repr) : 
    Class(name, repr),
    count_(0)
  {
  }

//-----------------------------------------------------------------------------
  IndexBase::~IndexBase()
  {
  }
  
//-----------------------------------------------------------------------------
//  Object::repr_t const IndexBase::repr() const
//  {
//    return repr_;
//  }

//-----------------------------------------------------------------------------
//  std::string const IndexBase::str() const
//  {
//    return str_;
//  }

//-----------------------------------------------------------------------------
//  void IndexBase::display() const
//  {
//  }

//-----------------------------------------------------------------------------
  Index::Index(dolfin::uint const& count) :
    IndexBase("Index", count)
  {
  }
  
//-----------------------------------------------------------------------------
  Index::Index(repr_t const& repr) :
    IndexBase("Index", repr)
  {
  }

//-----------------------------------------------------------------------------
  Index::~Index()
  {
  }
  
//-----------------------------------------------------------------------------
  Object::repr_t const Index::repr() const
  {
    return repr_;
  }

//-----------------------------------------------------------------------------
  std::string const Index::str() const
  {
    return str_;
  }

//-----------------------------------------------------------------------------
  void Index::display() const
  {
  }

//-----------------------------------------------------------------------------
  FixedIndex::FixedIndex(dolfin::uint const& count) :
    IndexBase("FixedIndex", count)
  {
  }

//-----------------------------------------------------------------------------
  FixedIndex::FixedIndex(repr_t const& repr) :
    IndexBase("FixedIndex", repr)
  {
  }

//-----------------------------------------------------------------------------
  FixedIndex::~FixedIndex()
  {
  }
  
//-----------------------------------------------------------------------------
  Object::repr_t const FixedIndex::repr() const
  {
    return repr_;
  }

//-----------------------------------------------------------------------------
  std::string const FixedIndex::str() const
  {
    return str_;
  }

//-----------------------------------------------------------------------------
  void FixedIndex::display() const
  {
  }

//-----------------------------------------------------------------------------
  MultiIndex::MultiIndex(dolfin::uint const& count) :
    IndexBase("MultiIndex", count)    
  {
  }

//-----------------------------------------------------------------------------
  MultiIndex::MultiIndex(IndexBase const& index) :
    IndexBase("MultiIndex", index)
  {
  }

//-----------------------------------------------------------------------------
  MultiIndex::MultiIndex(repr_t const& repr) :
    IndexBase("MultiIndex", repr)
  {
  }

//-----------------------------------------------------------------------------
  MultiIndex::~MultiIndex()
  {
  }
  
//-----------------------------------------------------------------------------
  Object::repr_t const MultiIndex::repr() const
  {
    return repr_;
  }

//-----------------------------------------------------------------------------
  std::string const MultiIndex::str() const
  {
    return str_;
  }

//-----------------------------------------------------------------------------
  void MultiIndex::display() const
  {
  }
}
