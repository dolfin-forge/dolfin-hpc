// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/fem/CoefficientMap.h>

#include <dolfin/log/log.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
bool CoefficientMap::has(std::string const& label) const
{
  return (this->count(label) > 0);
}

//-----------------------------------------------------------------------------
void CoefficientMap::disp() const
{
  section("CoefficientMap");
  // Begin indentation
  message("Number of coefficients : %u", this->size());
  uint ii = 0;
  dolfin::begin("");
  for ( CoefficientMap::value_type const & coeff : *this )
  {
    message("%8u : %s", ii, coeff.first.c_str());
    ++ii;
  }
  end();
  end();
}

}

