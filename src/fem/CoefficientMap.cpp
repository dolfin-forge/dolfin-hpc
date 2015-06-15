// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2014-11-10
// Last changed: 2014-11-10

#include <dolfin/fem/CoefficientMap.h>

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
Function * CoefficientMap::get(std::string const& label) const
{
  std::map<std::string, dolfin::Function *>::const_iterator it =
      map_.find(label);
  if (it == map_.end())
  {
    error("Coefficient label '%s' does not exist in map.", label.c_str());
  }
  return it->second;
}

//-----------------------------------------------------------------------------
uint CoefficientMap::size() const
{
  return this->map_.size();
}

//-----------------------------------------------------------------------------
void CoefficientMap::set(std::string const& label,
                         dolfin::Function& coefficient)
{
  std::map<std::string, dolfin::Function *>::iterator it = map_.find(
      label);
  if (it != map_.end())
  {
    error("Coefficient label '%s' is already assigned to a function.",
          label.c_str());
  }
  map_[label] = &coefficient;
}

//-----------------------------------------------------------------------------
void CoefficientMap::clear()
{
  map_.clear();
}

//-----------------------------------------------------------------------------
void CoefficientMap::disp() const
{
  cout << "CoefficientMap" << endl;
  cout << "--------------" << endl;

  // Begin indentation
  begin("");
  cout << "Number of coefficients      :" << this->size() << endl;
  uint ii = 0;
  for (std::map<std::string, dolfin::Function *>::const_iterator it =
      map_.begin(); it != map_.end(); ++it)
  {
    cout << ii << " : " << it->first << endl;
    ++ii;
  }
  end();
}

}

