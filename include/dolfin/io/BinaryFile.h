// Copyright (C) 2009-2015 Niclas Jansson.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Aurelien Larcher, 2014.
//
// First  added: 2009
// Last changed: 2015-07-25

#ifndef __DOLFIN_BINARY_FILE_H
#define __DOLFIN_BINARY_FILE_H

#include <stdint.h>
#include <dolfin/common/types.h>
#include <dolfin/common/byteswap.h>
#include <dolfin/la/Vector.h>
#include <dolfin/mesh/CellType.h>
#include "GenericFile.h"

#include <cstring>
#include <list>

#define BINARY_MAGIC_V1 0xBABE
#define BINARY_MAGIC_V2 0xB4B3
#define BINARY_MAGIC    BINARY_MAGIC_V2
#define FNAME_LENGTH    256

namespace dolfin
{

class BinaryFile : public GenericFile
{

public:

  ///
  BinaryFile(const std::string filename);

  ///
  BinaryFile(const std::string filename, real const& t);

  ///
  ~BinaryFile();

  /// Input
  void operator>>(GenericVector& x);
  void operator>>(Mesh& mesh);
  void operator>>(Function& f);
  void operator>>(std::vector<std::pair<Function*, std::string> >& f);
  void operator>>(MeshFunction<bool>& meshfunction);
  void operator>>(MeshFunction<int>& meshfunction);
  void operator>>(MeshFunction<uint>& meshfunction);
  void operator>>(MeshFunction<real>& meshfunction);

  /// Output
  void operator<<(GenericVector& x);
  void operator<<(Mesh& mesh);
  void operator<<(Function& u);
  void operator<<(std::vector<std::pair<Function*, std::string> >& f);
  void operator<<(MeshFunction<bool>& meshfunction);
  void operator<<(MeshFunction<int>& meshfunction);
  void operator<<(MeshFunction<uint>& meshfunction);
  void operator<<(MeshFunction<real>& meshfunction);

  /// Overload GenericFile
  void read();
  void write();

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

  typedef struct atomic_cell
  {
    uint const size;
    uint * const v;

    //-----------------------------------
    atomic_cell(uint d) :
        size(d),
        v(new uint[size])
    {
    }
    //-----------------------------------
    atomic_cell(atomic_cell const& other) :
        size(other.size),
        v(new uint[size])
    {
      std::copy(other.v, other.v + size, v);
    }
    //-----------------------------------
    atomic_cell& operator=(atomic_cell const& other)
    {
      if(&other == this)
      {
        return *this;
      }
      if(size != other.size)
      {
        error("Size of atomic_cells in assignment do not match");
      }
      std::copy(other.v, other.v + size, v);
      return *this;
    }
    //-----------------------------------
    ~atomic_cell()
    {
      delete[] v;
    }
  } atomic_cell;

  template<typename T>
    void write_meshfunction(MeshFunction<T>& meshfunction);

  template<class T>
    void read_meshfunction(MeshFunction<T>& meshfunction);

  uint vertex_owner(uint L, uint R, uint i);

  void nameUpdate(const int counter);

  void write_function(std::vector<std::pair<Function*, std::string> >& f);

  bool hdr_check(BinaryFileHeader& hdr, Binary_data_t type, uint pe_size);

#ifdef ENABLE_MPIIO
  void bswap_func_hdr(BinaryFunctionHeader& hdr);
#endif

  /// Returns binary file cell type identifier for given DOLFIN cell type
  uint cell_type(CellType::Type const type);

  /// Returns DOLFIN cell type for given binary file cell type identifier
  CellType::Type cell_type(uint const type);

  // Function filename
  std::string bin_filename_;

  // Current time
  real const * const t_;

  // Version number of the binary file
  uint version_;

};

//--- INLINES -----------------------------------------------------------------
inline uint BinaryFile::vertex_owner(uint L, uint R, uint i)
{
  return static_cast<uint>(std::max(
      std::floor((double) i / (double) (L + 1)),
      std::floor((double) ((double) i - (double) R) / (double) L)));
}

//-----------------------------------------------------------------------------
inline bool BinaryFile::hdr_check(BinaryFileHeader& hdr, Binary_data_t type,
                                  uint pe_size)
{

  bool byteswap = false;

  if (hdr.magic == BINARY_MAGIC_V2)
  {
    message(1, "Loading Binary File format version 2");
    version_ = 2;
  }
  else if (hdr.magic == BINARY_MAGIC_V1)
  {
    message(1, "Loading Binary File format version 1");
    version_ = 1;
  }
  else if (bswap(hdr.magic) == BINARY_MAGIC_V2)
  {
    message(1, "Loading Binary File format version 2 (endian conversion)");
    version_ = 2;
    byteswap = true;
  }
  else if (bswap(hdr.magic) == BINARY_MAGIC_V1)
  {
    message(1, "Loading Binary File format version 1 (endian conversion)");
    version_ = 1;
    byteswap = true;
  }
  else
  {
    error("Corrupt header: invalid magic number (%0x)", hdr.magic);
  }

  if (byteswap) 
  {
    hdr.magic = bswap(hdr.magic);
    hdr.bendian = bswap(hdr.bendian);
    hdr.pe_size = bswap(hdr.pe_size);
    hdr.type = static_cast<Binary_data_t>(bswap(hdr.type));
  }

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

  return byteswap;
}

//-----------------------------------------------------------------------------
#ifdef ENABLE_MPIIO
inline void BinaryFile::bswap_func_hdr(BinaryFunctionHeader& hdr)
{
  hdr.dim = bswap(hdr.dim);
  hdr.size = bswap(hdr.size);
  hdr.t = bswap(hdr.t);
}
#endif
//-----------------------------------------------------------------------------
inline uint BinaryFile::cell_type(CellType::Type const type)
{
  switch (version_)
    {
    case 2:
      switch (type)
        {
        case CellType::point:
          return 0;
          break;
        case CellType::interval:
          return 1;
          break;
        case CellType::triangle:
          return 2;
          break;
        case CellType::tetrahedron:
          return 3;
          break;
        case CellType::quadrilateral:
          return 4;
          break;
        case CellType::hexahedron:
          return 6;
          break;
        default:
          error("Unsupported mesh cell type in BinaryFile V2.");
          break;
        }
      break;
    case 1:
      switch (type)
        {
        case CellType::triangle:
          return 0;
          break;
        case CellType::tetrahedron:
          return 1;
          break;
        default:
          error("Unsupported mesh cell type in BinaryFile V1.");
          break;
        }
      break;
    default:
      error("Invalid version compatibility number for cell type detection.");
      break;
    }
  return 0;
}

//-----------------------------------------------------------------------------
inline CellType::Type BinaryFile::cell_type(uint const type)
{
  switch (version_)
    {
    case 2:
      switch (type)
        {
        case 0:
          return CellType::point;
          break;
        case 1:
          return CellType::interval;
          break;
        case 2:
          return CellType::triangle;
          break;
        case 3:
          return CellType::tetrahedron;
          break;
        case 4:
          return CellType::quadrilateral;
          break;
        case 6:
          return CellType::hexahedron;
          break;
        default:
          error("Unsupported binary cell type in BinaryFile V2.");
          break;
        }
      break;
    case 1:
      switch (type)
        {
        case 0:
          return CellType::triangle;
          break;
        case 1:
          return CellType::tetrahedron;
          break;
        default:
          error("Unsupported binary cell type in BinaryFile V1.");
          break;
        }
      break;
    default:
      error("Invalid version compatibility number for cell type detection.");
      break;
    }
  return CellType::point;
}

}
#endif
