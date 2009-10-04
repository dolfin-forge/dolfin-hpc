// Copyright (C) 2009 Niclas Jansson.
// Licensed under the GNU LGPL Version 2.1.
//

#include <fstream>
#include <dolfin/main/MPI.h>
#include <dolfin/la/Vector.h>
#include "BinaryFile.h"

#ifdef HAS_MPI
#include <mpi.h>
#endif

using namespace dolfin;

//----------------------------------------------------------------------------
BinaryFile::BinaryFile(const std::string filename) : GenericFile(filename)
{
  type = "Binary";
}
//----------------------------------------------------------------------------
BinaryFile::~BinaryFile()
{
  // Do nothing
}
//----------------------------------------------------------------------------
void BinaryFile::operator>>(GenericVector& x)
{

  std::ifstream fp(filename.c_str(), std::ifstream::binary);

  uint local_size, size;
  fp.read((char *)&local_size, sizeof(uint));
  
#ifdef HAS_MPI
  MPI_Allreduce(&local_size, &size, 1, 
		MPI_UNSIGNED, MPI_SUM, MPI::DOLFIN_COMM);
#else
  local_size = size;
#endif
  
  real *values = new real[local_size];  

  fp.read((char *)values, local_size * sizeof(real));
  fp.close();

  x.init(size);
  x.set(values);
  delete[] values;  
}
//----------------------------------------------------------------------------
void BinaryFile::operator<<(GenericVector& x)
{
  
  std::ofstream fp(filename.c_str(), std::ofstream::binary);

  real *values = new real[x.local_size()];
  uint size = x.local_size();
  x.get(values);

  fp.write((char *)&size, sizeof(uint));
  fp.write((char *)values, x.local_size() * sizeof(real));
  fp.close();

  delete[] values;     

  message(1, "Saved vector  to file %s in binary format.", filename.c_str());  
}
//----------------------------------------------------------------------------

