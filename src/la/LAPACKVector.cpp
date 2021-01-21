// Copyright (C) 2009 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/la/LAPACKVector.h>

#include <cstring>
#include <sstream>

namespace dolfin
{

//-----------------------------------------------------------------------------
LAPACKVector::LAPACKVector(uint M) :
    M_(M),
    values_(new double[M])
{
  std::memset(values_, 0, sizeof(double)*M_);
}

//-----------------------------------------------------------------------------
LAPACKVector::~LAPACKVector()
{
  delete[] values_;
}

//-----------------------------------------------------------------------------
auto LAPACKVector::str(bool verbose) const -> std::string
{
  std::stringstream s;

  if (verbose)
  {
    s << str(false) << std::endl << std::endl;

    for (uint i = 0; i < M_; i++)
    {
      s << (*this)[i];
      s << std::endl;
    }
  }
  else
  {
    s << "<LAPACK vector of size " << M_ << ">";
  }

  return s.str();
}
//-----------------------------------------------------------------------------

}
