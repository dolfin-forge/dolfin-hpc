// Copyright (C) 2009-2015 Niclas Jansson.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_BINARY_FILE_H
#define __DOLFIN_BINARY_FILE_H

#include <dolfin/io/GenericFile.h>

#include <dolfin/mesh/CellType.h>

#include <string>

#define BINARY_MAGIC_V1 0xBABE
#define BINARY_MAGIC_V2 0xB4B3
#define BINARY_MAGIC    BINARY_MAGIC_V2
#define FNAME_LENGTH    256
#define BINARY_VERSION  2

namespace dolfin
{

class CellType;
class Vector;

class BinaryFile : public GenericFile
{

public:

  ///
  BinaryFile(const std::string filename);

  ///
  BinaryFile(const std::string filename, real const& t);

  ///
  ~BinaryFile() = default;

  /// Input
  void operator>>(GenericVector& x);
  void operator>>(Mesh& mesh);
  void operator>>(Function& f);
  void operator>>(LabelList<Function>& f);
  void operator>>(MeshFunction<bool>& meshfunction);
  void operator>>(MeshFunction<int>& meshfunction);
  void operator>>(MeshFunction<uint>& meshfunction);
  void operator>>(MeshFunction<real>& meshfunction);

  /// Output
  void operator<<(GenericVector& x);
  void operator<<(Mesh& mesh);
  void operator<<(Function& u);
  void operator<<(LabelList<Function>& f);
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

  template<typename T>
    void write_meshfunction(MeshFunction<T>& meshfunction);

  template<class T>
    void read_meshfunction(MeshFunction<T>& meshfunction);

  void nameUpdate(const int counter);

  void write_function(LabelList<Function>& f);

  bool hdr_check(BinaryFileHeader& hdr, Binary_data_t type, uint pe_size);

#ifdef ENABLE_MPIIO
  void bswap_func_hdr(BinaryFunctionHeader& hdr);
#endif

  /// Returns binary file cell type identifier for given DOLFIN cell type
  uint cell_type(uint version, CellType::Type const type);

  /// Returns DOLFIN cell type for given binary file cell type identifier
  CellType::Type cell_type(uint version, uint const type);

  // Function filename
  std::string bin_filename_;

  // Current time
#if defined( HAVE_MPI )
  real const * const t_;
#endif

  // Version number of the binary file
  uint version_;

};

} /* namespace dolfin */

#endif /* __DOLFIN_BINARY_FILE_H */
