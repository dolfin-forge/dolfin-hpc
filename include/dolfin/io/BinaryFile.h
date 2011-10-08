// Copyright (C) 2009-2011 Niclas Jansson.
// Licensed under the GNU LGPL Version 2.1.
//
// First  added: 2009
// Last changed: 2011-06-16


#ifndef __BINARY_FILE_H
#define __BINARY_FILE_H

#include <stdint.h>
#include <dolfin/common/types.h>
#include <dolfin/la/Vector.h>
#include "GenericFile.h"


#define BINARY_MAGIC 0xBABE
#define FNAME_LENGTH 256

namespace dolfin
{
  class BinaryFile : public GenericFile
  {
    
  public:
    BinaryFile(const std::string filename);
    BinaryFile(const std::string filename, real& t);
    ~BinaryFile();
    
    // Input
    void operator>> (GenericVector& x);
    void operator>> (Mesh& mesh);
    
    // Output
    void operator<< (GenericVector& x);
    void operator<< (Mesh& mesh);
    void operator<< (Function& u);
    void operator<< (MeshFunction<int>& meshfunction);
    void operator<< (MeshFunction<unsigned int>& meshfunction);
    void operator<< (MeshFunction<double>& meshfunction);
    void operator<< (std::vector<std::pair<Function*, std::string> >& f);

  private:
    
    enum Binary_data_t { BINARY_MESH_DATA, 
			 BINARY_VECTOR_DATA,
			 BINARY_FUNCTION_DATA,
			 BINARY_MESH_FUNCTION_DATA};
    

    typedef struct {
      uint32_t magic;
      uint32_t bendian; 
      uint32_t pe_size;
      Binary_data_t type;
    } BinaryFileHeader;

#ifdef ENABLE_MPIIO
    typedef struct {
      uint32_t dim;
      uint32_t size;
      real t;
      char name[FNAME_LENGTH];      
    } BinaryFunctionHeader;
#endif      

    typedef struct {
      uint v1; 
      uint v2;
      uint v3;
      uint v4;
    } atomic_cell;

    template<class T>
    void write_meshfunction(T& meshfunction);    

    inline int vertex_owner(uint L, uint R, uint i ) 
    {
      return (int) std::max( floor( (double) i / (double) (L + 1) ),
			     floor( (double) ((double) i - (double) R) / 
				    (double)  L));
    };


    void nameUpdate(const int counter);

    void write_function(std::vector<std::pair<Function*, std::string> >& f);

    inline void hdr_check(BinaryFileHeader hdr, 
			  Binary_data_t type, uint pe_size)
    {     
      
      if (hdr.magic != BINARY_MAGIC)
	error("Corrupt header");

#ifdef HAVE_BIG_ENDIAN
      if (!hdr.bendian)
	error("File written in little endian");
#else
      if (hdr.bendian)
	error("File written in big endian");
#endif

      if (hdr.type != type)
	error("Invalid data type in file");
      
      if (hdr.type == BINARY_VECTOR_DATA && (hdr.pe_size != pe_size))
	error("File stored on %d PE's, currently running on %d PE's", 
	      hdr.pe_size, pe_size);
    };


    // function filename
    std::string bin_filename;

    // Current time
    real* _t;
  };
}
#endif
