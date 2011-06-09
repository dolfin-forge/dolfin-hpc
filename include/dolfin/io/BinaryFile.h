// Copyright (C) 2009-2011 Niclas Jansson.
// Licensed under the GNU LGPL Version 2.1.
//
// First  added: 2009
// Last changed: 2011-06-08


#ifndef __BINARY_FILE_H
#define __BINARY_FILE_H

#include <stdint.h>
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
    
#ifdef ENABLE_MPIIO
    typedef struct {
      uint32_t bendian; 
      uint32_t pe_size;
    } BinaryFileHeader;
#endif      


    typedef struct {
      uint v1; 
      uint v2;
      uint v3;
      uint v4;
    } atomic_cell;

    inline int vertex_owner(uint L, uint R, uint i ) 
    {
      return (int) std::max( floor( (double) i / (double) (L + 1) ),
			     floor( (double) ((double) i - (double) R) / 
				    (double)  L));
    };

    inline void hdr_check(BinaryFileHeader hdr, uint pe_size, bool check_pe_size)
    {     
#ifdef HAVE_BIG_ENDIAN
      if (!hdr.bendian || (check_pe_size && hdr.pe_size != pe_size))
	error("Shut her down, Scotty, she's sucking mud again!");
#else
      message("hdr.bendian: %d hdr.pe_size: %d sizeof(hdr): %d", hdr.bendian, hdr.pe_size, sizeof(BinaryFileHeader));
      if (hdr.bendian || (check_pe_size && hdr.pe_size != pe_size))
	error("Shut her down, Scotty, she's sucking mud again!");
#endif
    };
  };
}
#endif
