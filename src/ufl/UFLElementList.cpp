// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/ufl/UFLElementList.h>

#include <dolfin/log/log.h>
#include <dolfin/ufl/UFLFiniteElementSpace.h>

namespace ufl
{

//-----------------------------------------------------------------------------
ElementList::ElementList() :
    it_(this->end())
{
}

//-----------------------------------------------------------------------------
ElementList::~ElementList()
{
  for (ElementList::iterator it = this->begin(); it != this->end(); ++it)
  {
    delete it->second;
  }
}

//-----------------------------------------------------------------------------
ElementList::ElementList(ElementList const &other) :
  std::map<Object::repr_t, FiniteElementSpace *>( other)
{
  this->clear();
  for (ElementList::const_iterator it = other.begin(); it != other.end(); ++it)
  {
    this->add(it->first);
  }
}

//-----------------------------------------------------------------------------
void ElementList::add(Object::repr_t const& signature)
{
  if (this->find(signature) == this->end())
  {
    this->insert(ElementItem(signature, FiniteElementSpace::create(signature)));
  }
}

//-----------------------------------------------------------------------------
bool ElementList::has(Object::repr_t const& signature) const
{
  return (this->find(signature) != this->end());
}

//-----------------------------------------------------------------------------
FiniteElementSpace const * ElementList::first() const
{
  it_ = this->begin();
  if (this->empty())
  {
    return NULL;
  }
  return it_->second;
}

//-----------------------------------------------------------------------------
FiniteElementSpace const * ElementList::next() const
{
  if (it_ != this->end())
  {
    ++it_;
  }
  return (it_ == this->end()) ? NULL : it_->second;
}

//-----------------------------------------------------------------------------
bool ElementList::valid() const
{
  return (it_ != this->end());
}

//-----------------------------------------------------------------------------
bool ElementList::check() const
{
  bool ret = true;
  for (ElementList::const_iterator it = this->begin(); it != this->end(); ++it)
  {
    ret &= (it->first == it->second->repr());
  }
  return ret;
}

//-----------------------------------------------------------------------------
void ElementList::disp() const
{
  for (ElementList::const_iterator it = this->begin(); it != this->end(); ++it)
  {
    dolfin::message("@%p : %s", it->second, it->first.c_str());
#if DEBUG
    it->second->display();
#endif
  }
}

//-----------------------------------------------------------------------------

}

