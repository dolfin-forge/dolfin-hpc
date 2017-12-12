// Copyright (C) 2003-2007 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2003-05-06
// Last changed: 2007-04-13

#include <dolfin/parameter/ParameterList.h>

#include <dolfin/log/log.h>

#include <string>

namespace dolfin
{

//-----------------------------------------------------------------------------
ParameterList::ParameterList()
{
  // Do nothing
}
//-----------------------------------------------------------------------------
ParameterList::~ParameterList()
{
  // Do nothing
}
//-----------------------------------------------------------------------------
void ParameterList::add(std::string key, Parameter value)
{
  if (defined(key))
  {
    error("Unable to add parameter \"%s\" (already defined).", key.c_str());
  }

  storage_.insert(Item(key, value));
}
//-----------------------------------------------------------------------------
void ParameterList::set(std::string key, Parameter value)
{
  iterator p = storage_.find(key);

  if (p == storage_.end())
  {
    error("Unknown parameter \"%s\".", key.c_str());
  }

  p->second = value;
}
//-----------------------------------------------------------------------------
Parameter ParameterList::get(std::string const& key) const
{
  const_iterator p = storage_.find(key);

  if (p == storage_.end())
  {
    error("Unknown parameter \"%s\".", key.c_str());
  }
  
  return p->second;
}
//-----------------------------------------------------------------------------
bool ParameterList::defined(std::string const& key) const
{
  return (storage_.find(key) != storage_.end());
}
//-----------------------------------------------------------------------------
bool ParameterList::operator==(ParameterList const& other) const
{
  return (this->storage_ == other.storage_);
}
//-----------------------------------------------------------------------------
bool ParameterList::operator!=(ParameterList const& other) const
{
  return (this->storage_ != other.storage_);
}
//-----------------------------------------------------------------------------
ParameterList& ParameterList::operator<<(ParameterList const& other)
{
  if (this != &other)
  {
    for (ParameterList::const_iterator it = other.begin(); it != other.end();
         ++it)
    {
      if (this->find(it->first) != this->end())
      {
        this->set(it->first, it->second);
      }
    }
  }
  return *this;
}
//-----------------------------------------------------------------------------
ParameterList const& ParameterList::operator>>(ParameterList& other) const
{
  if (this != &other)
  {
    for (ParameterList::const_iterator it = this->begin(); it != this->end();
         ++it)
    {
      if (other.find(it->first) != other.end())
      {
        other.set(it->first, it->second);
      }
    }
  }
  return *this;
}
//-----------------------------------------------------------------------------

} /* namespace dolfin */
