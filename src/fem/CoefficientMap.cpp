// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2014-11-10
// Last changed: 2014-11-10

#include <dolfin/fem/CoefficientMap.h>

#include <dolfin/log/log.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
CoefficientMap::CoefficientMap()
{
}

//-----------------------------------------------------------------------------
CoefficientMap::~CoefficientMap()
{
  map_.clear();
}

//-----------------------------------------------------------------------------
CoefficientMap::CoefficientMap(CoefficientMap const& other)
{
  *this = other;
}

//-----------------------------------------------------------------------------
CoefficientMap& CoefficientMap::operator=(CoefficientMap const& other)
{
  if (this == &other)
  {
    return *this;
  }
  this->map_ = other.map_;
  return *this;
}

//-----------------------------------------------------------------------------
bool CoefficientMap::has(std::string const& label) const
{
  return (map_.count(label) > 0);
}

//-----------------------------------------------------------------------------
uint CoefficientMap::size() const
{
  return this->map_.size();
}

//-----------------------------------------------------------------------------
void CoefficientMap::clear()
{
  map_.clear();
}

//-----------------------------------------------------------------------------
void CoefficientMap::disp() const
{
  section("CoefficientMap");
  // Begin indentation
  message("Number of coefficients : %u", this->size());
  uint ii = 0;
  begin("");
  for (Container::const_iterator it = map_.begin(); it != map_.end(); ++it)
  {
    message("%8u : %s", ii, it->first.c_str());
    ++ii;
  }
  end();
  end();
}

}

