// Copyright (C) 2009 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/la/LAPACKMatrix.h>

#include <sstream>

namespace dolfin
{

//-----------------------------------------------------------------------------
LAPACKMatrix::LAPACKMatrix(uint M, uint N) :
    M_(M),
    N_(N),
    values_(new double[N * M])
{
  memset(values_, 0, sizeof(double)*M_*N_);
}

//-----------------------------------------------------------------------------
LAPACKMatrix::~LAPACKMatrix()
{
  delete[] values_;
}

//-----------------------------------------------------------------------------
std::string LAPACKMatrix::str(bool verbose) const
{
  std::stringstream s;

  if (verbose)
  {
    s << str(false) << std::endl << std::endl;

    for (uint i = 0; i < M_; i++)
    {
      for (uint j = 0; j < N_; j++)
      {
        s << (*this)(i, j);
        if (j < N_ - 1)
          s << " ";
      }
      s << std::endl;
    }
  }
  else
  {
    s << "<LAPACK matrix of size " << M_ << " x " << N_ << ">";
  }

  return s.str();
}
//-----------------------------------------------------------------------------

}

