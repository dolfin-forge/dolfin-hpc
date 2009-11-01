// Copyright (C) 2009 Niclas Jansson.
// Licensed under the GNU LGPL Version 2.1.
//

#include <fstream>
#include <dolfin/common/types.h>
#include <dolfin/la/Vector.h>
#include "BinaryFile.h"

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

  uint size;
  fp.read((char *)&size, sizeof(uint));
  
  real *values = new real[size];  

  fp.read((char *)values, size * sizeof(real));
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

