// Copyright (C) 2009-2011 Niclas Jansson.
// Licensed under the GNU LGPL Version 2.1.
//
// First  added: 2009
// Last changed: 2011-04-02


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
    
    typedef struct __apa__{
      uint v1; 
      uint v2;
      uint v3;
      uint v4;
    } atomic_cell;

    inline int vertex_owner(uint L, uint R, uint i ) 
    {
      return (int) fmax( floor( (double) i / (double) (L + 1) ),
			 floor( (double) ((double) i - (double) R) / (double)  L));
    };
  };
}
#endif
