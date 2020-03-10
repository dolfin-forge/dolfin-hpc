// Copyright (C) 2005-2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/parameter/Parametrized.h>

#include <dolfin/log/dolfin_log.h>
#include <dolfin/parameter/parameters.h>
#include <dolfin/parameter/ParameterSystem.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
Parametrized::Parametrized() :
    parameters_(),
    parent_(NULL)
{
  // Do nothing
}
//-----------------------------------------------------------------------------
Parametrized::~Parametrized()
{
  // Do nothing
}
//-----------------------------------------------------------------------------
bool Parametrized::operator==(Parametrized const& other)
{
  return (parameters_ == other.parameters_);
}
//-----------------------------------------------------------------------------
bool Parametrized::operator!=(Parametrized const& other)
{
  return (parameters_ != other.parameters_);
}
//-----------------------------------------------------------------------------
void Parametrized::add(std::string key, Parameter value)
{
  parameters_.add(key, value);
}
//-----------------------------------------------------------------------------
void Parametrized::set(std::string key, Parameter value)
{
  if ( !has(key) )
  {
    parameters_.add(key, value);
  }
  else
  {
    parameters_.set(key, value);
  }
}
//-----------------------------------------------------------------------------
void Parametrized::set(std::string key, Parametrized const& parent)
{
  // Check that key is "parent"
  if ( !(key == "parent") )
  {
    error("Illegal value for parameter \"%s\".", key.c_str());
  }

  // Check if we already have a parent
  if ( this->parent_ )
  {
    error("Local parameter database can only have one parent.");
  }

  // Check that parent is not itself
  if ( this == &parent )
  {
    error("Local parameter database cannot be its own parent.");
  }

  // Set parent
  this->parent_ = &parent;
}
//-----------------------------------------------------------------------------
Parameter Parametrized::get(std::string key) const
{
  // First check local database
  if ( has(key) )
  {
    return parameters_.get(key);
  }

  // Check parent if any
  if ( parent_ )
  {
    return parent_->get(key);
  }

  // Fall back on global database
  return ParameterSystem::parameters.get( key );
}
//-----------------------------------------------------------------------------
bool Parametrized::has(std::string key) const
{
  return parameters_.defined(key);
}
//-----------------------------------------------------------------------------
Parametrized& Parametrized::operator<<(ParameterList const& p)
{
  parameters_ << p; return *this;
}
//-----------------------------------------------------------------------------
Parametrized const& Parametrized::operator>>(ParameterList& p) const
{
  if ( parent_ ) { parent_->parameters_ >> p; }
  parameters_ >> p; return *this;
}
//-----------------------------------------------------------------------------
Parametrized& Parametrized::operator<<(Parametrized const& p)
{
  if ( p.parent_ ) {  (*this) << *(p.parent_); }
  parameters_ << p.parameters_; return *this;
}
//-----------------------------------------------------------------------------
Parametrized const& Parametrized::operator>>(Parametrized& p) const
{
  if ( parent_ ) { *parent_ >> p; }
  parameters_ >> p.parameters_; return *this;
}
//-----------------------------------------------------------------------------
void Parametrized::disp() const
{
  parameters_.disp();
}
//-----------------------------------------------------------------------------

} /* namespace dolfin */
