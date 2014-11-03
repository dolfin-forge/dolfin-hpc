// Copyright (C) 2009-2012 Niclas Jansson.
// Licensed under the GNU LGPL Version 2.1.
//
// First  added: 2009
// Last changed: 2012-06-12

#ifndef __BINARY_FILE_H
#define __BINARY_FILE_H

#include <stdint.h>
#include <dolfin/common/types.h>
#include <dolfin/la/Vector.h>
#include <dolfin/mesh/CellType.h>
#include "GenericFile.h"

#define BINARY_MAGIC 0xBABE
#define FNAME_LENGTH 256

namespace dolfin
{

class BinaryFile : public GenericFile
{

public:

  ///
  BinaryFile(const std::string filename);

  ///
  BinaryFile(const std::string filename, real& t);

  ///
  ~BinaryFile();

  /// Input
  void operator>>(GenericVector& x);
  void operator>>(Mesh& mesh);
  void operator>>(Function& f);
  void operator>>(std::vector<std::pair<Function*, std::string> >& f);
  void operator>>(MeshFunction<bool>& meshfunction);
  void operator>>(MeshFunction<int>& meshfunction);
  void operator>>(MeshFunction<unsigned int>& meshfunction);
  void operator>>(MeshFunction<double>& meshfunction);

  /// Output
  void operator<<(GenericVector& x);
  void operator<<(Mesh& mesh);
  void operator<<(Function& u);
  void operator<<(std::vector<std::pair<Function*, std::string> >& f);
  void operator<<(MeshFunction<bool>& meshfunction);
  void operator<<(MeshFunction<int>& meshfunction);
  void operator<<(MeshFunction<unsigned int>& meshfunction);
  void operator<<(MeshFunction<double>& meshfunction);

  enum Binary_data_t
  {
    BINARY_MESH_DATA,
    BINARY_VECTOR_DATA,
    BINARY_FUNCTION_DATA,
    BINARY_MESH_FUNCTION_DATA
  };

  typedef struct
  {
    uint32_t magic;
    uint32_t bendian;
    uint32_t pe_size;
    Binary_data_t type;
  } BinaryFileHeader;

#ifdef ENABLE_MPIIO
  typedef struct
  {
    uint32_t dim;
    uint32_t size;
    real t;
    char name[FNAME_LENGTH];
  } BinaryFunctionHeader;
#endif

private:

  typedef struct
  {
    uint v1;
    uint v2;
    uint v3;
    uint v4;
  } atomic_cell;

  template<typename T>
    void write_meshfunction(MeshFunction<T>& meshfunction);

  template<class T>
    void read_meshfunction(MeshFunction<T>& meshfunction);

  int vertex_owner(uint L, uint R, uint i);

  void nameUpdate(const int counter);

  void write_function(std::vector<std::pair<Function*, std::string> >& f);

  void hdr_check(BinaryFileHeader hdr, Binary_data_t type, uint pe_size);

  ///
  int cell_type(CellType::Type const type);
  CellType::Type cell_type(uint const type);

  // Function filename
  std::string bin_filename_;

  // Current time
  real* t_;
};

//--- INLINES -----------------------------------------------------------------

inline int BinaryFile::vertex_owner(uint L, uint R, uint i)
{
  return (int) std::max(
      std::floor((double) i / (double) (L + 1)),
      std::floor((double) ((double) i - (double) R) / (double) L));
}

//-----------------------------------------------------------------------------
inline void BinaryFile::hdr_check(BinaryFileHeader hdr, Binary_data_t type,
                                  uint pe_size)
{
  if (hdr.magic != BINARY_MAGIC)
  {
    error("Corrupt header");
  }

#ifdef HAVE_BIG_ENDIAN
  if (!hdr.bendian)
  {
    error("File written in little endian");
  }
#else
  if (hdr.bendian)
  {
    error("File written in big endian");
  }
#endif

  if (hdr.type != type)
  {
    error("Invalid data type in file");
  }

  if ((hdr.type == BINARY_FUNCTION_DATA || hdr.type == BINARY_VECTOR_DATA)
      && (hdr.pe_size != pe_size))
  {
    error("File stored on %d PEs, currently running on %d PEs", hdr.pe_size,
          pe_size);
  }
}

//-----------------------------------------------------------------------------
inline int BinaryFile::cell_type(CellType::Type const type)
{
  switch (type)
    {
    case CellType::point:
      return -2;
      break;
    case CellType::interval:
      return -1;
      break;
    case CellType::triangle:
      return 0;
      break;
    case CellType::tetrahedron:
      return 1;
      break;
    default:
      error("Unsupported mesh cell type in BinaryFile.");
      break;
    }
  return 0;
}

//-----------------------------------------------------------------------------
inline CellType::Type BinaryFile::cell_type(uint const type)
{
  switch (type)
    {
    case -2:
      return CellType::point;
      break;
    case -1:
      return CellType::interval;
      break;
    case 0:
      return CellType::triangle;
      break;
    case 1:
      return CellType::tetrahedron;
      break;
    default:
      error("Unsupported binary cell type in BinaryFile.");
      break;
    }
  return CellType::point;
}

}
#endif
