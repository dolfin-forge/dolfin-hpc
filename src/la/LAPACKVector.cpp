// Copyright (C) 2009 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2009-12-15
// Last changed: 2009-12-15

#include <dolfin/la/LAPACKVector.h>

#include <sstream>

namespace dolfin
{

//-----------------------------------------------------------------------------
LAPACKVector::LAPACKVector(uint M) :
    M_(M),
    values_(new double[M])
{
  memset(values_, 0, sizeof(double)*M_);
}

//-----------------------------------------------------------------------------
LAPACKVector::~LAPACKVector()
{
  delete[] values_;
}

//-----------------------------------------------------------------------------
std::string LAPACKVector::str(bool verbose) const
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
