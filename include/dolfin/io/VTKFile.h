// Copyright (C) 2005-2007 Garth N. Wells.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Anders Logg 2006.
// Modified by Niclas Jansson 2008-2010.
// Modified by Aurélien Larcher 2014.
//
// First added:  2005-07-05
// Last changed: 2014-02-08

#ifndef __DOLFIN_PVTK_FILE_H
#define __DOLFIN_PVTK_FILE_H

#include <fstream>
#include <string>
#include <vector>
#include "GenericFile.h"

namespace dolfin
{

/**
 *  @class  VTKFile
 *
 *  @brief  Provides an output writer to Kitware's VTK file format for Meshes,
 *          MeshFunctions and Functions.
 */

class VTKFile : public GenericFile
{
public:

  ///
  VTKFile(const std::string filename);

  ///
  VTKFile(const std::string filename, real const& t);

  ///
  ~VTKFile();

  /// Write mesh
  void operator<<(Mesh& mesh);

  /// Write mesh function containing integers
  void operator<<(MeshFunction<int>& meshfunction);

  /// Write mesh function containing unsigned integers
  void operator<<(MeshFunction<uint>& meshfunction);

  /// Write mesh function containing doubles
  void operator<<(MeshFunction<real>& meshfunction);

  /// Write mesh function containing booleans
  void operator<<(MeshFunction<bool>& meshfunction);

  /// Write function
  void operator<<(Function& u);

  /// Write list of functions
  void operator<<(LabelList<Function>& f);

  /// Overload GenericFile
  void read();
  void write();

  /// Compute base64 encoded stream for VTK
  template<typename T>
    void encode_stream(std::stringstream& stream,
                       const std::vector<T>& data) const;

private:

  void write_dataset(LabelList<Function>& f);

  // Compute base64 encoded stream for VTK
  template<typename T>
    void encode_inline_base64(std::stringstream& stream,
                              const std::vector<T>& data) const;

  // Compute compressed base64 encoded stream for VTK
  template<typename T>
    void encode_inline_compressed_base64(std::stringstream& stream,
                                         const std::vector<T>& data) const;

  void MeshWrite(Mesh& mesh) const;
  void ResultsWrite(LabelList<Function> f) const;
  void pvdFileWrite(uint u);
  void pvtuFileWrite(bool mesh_function, uint const dim);
  void pvtuFileWriteFunction(LabelList<Function> f);
  void VTKHeaderOpen(Mesh& mesh) const;
  void VTKHeaderClose() const;
  void vtuNameUpdate(const int counter);
  void pvtuNameUpdate(const int counter);

  template<class T>
    void MeshFunctionWrite(MeshFunction<T>& meshfunction);

  // Most recent position in pvd file
  std::ios::pos_type mark;

  // vtu filename
  std::string vtu_filename;

  // pvtu filename
  std::string pvtu_filename;

  // Current time
  real const * const _t;

  uint _rank;

};

}

#endif
