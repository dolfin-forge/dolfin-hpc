// Copyright (C) 2009-2011 Niclas Jansson.
// Licensed under the GNU LGPL Version 2.1.
//
// First  added: 2009
// Last changed: 2011-04-01


#ifndef __BINARY_FILE_H
#define __BINARY_FILE_H

#include <dolfin/common/types.h>
#include <dolfin/la/Vector.h>
#include "GenericFile.h"


namespace dolfin
{
  class BinaryFile : public GenericFile
  {
    
  public:
    BinaryFile(const std::string filename);
    ~BinaryFile();
    
    // Input
    void operator>> (GenericVector& x);
    void operator>> (Mesh& mesh);
    
    // Output
    void operator<< (GenericVector& x);
    void operator<< (Mesh& mesh);

  private:
    
#ifdef HAVE_MPI
    void partition_vertices(int dim, int local_vertices, real *vertices);
#endif
    
  };
}
#endif
