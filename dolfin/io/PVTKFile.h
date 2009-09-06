// Copyright (C) 2005-2007 Garth N. Wells.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Anders Logg 2006.
// Modified by Niclas Jansson 2008-2009.
//
// First added:  2005-07-05
// Last changed: 2009-09-06

#ifndef __PVTK_FILE_H
#define __PVTK_FILE_H

#include <fstream>
#include <string>
#include <vector>
#include "GenericFile.h"

namespace dolfin
{

  class PVTKFile : public GenericFile
  {
  public:
    
    PVTKFile(const std::string filename);
    ~PVTKFile();
    
    void operator<< (Mesh& mesh);
    void operator<< (MeshFunction<int>& meshfunction);
    void operator<< (MeshFunction<unsigned int>& meshfunction);
    void operator<< (MeshFunction<double>& meshfunction);
    void operator<< (Function& u);
    
    void write();

    // Compute base64 encoded stream for VTK
    template<typename T>
    void encode_stream(std::stringstream& stream, const std::vector<T>& data) const;

  private:

    // Compute base64 encoded stream for VTK
    template<typename T>
    void encode_inline_base64(std::stringstream& stream, const std::vector<T>& data) const;

    // Compute compressed base64 encoded stream for VTK
    template<typename T>
    void encode_inline_compressed_base64(std::stringstream& stream, const std::vector<T>& data) const;

    void MeshWrite(Mesh& mesh) const;
    void ResultsWrite(Function& u) const;
    void pvdFileWrite(uint u);
    void pvtuFileWrite();
    void pvtuFileWrite_func(Function& u);
    void VTKHeaderOpen(Mesh& mesh) const;
    void VTKHeaderClose() const;
    void vtuNameUpdate(const int counter);
    void pvtuNameUpdate(const int counter);

    template<class T>
    void MeshFunctionWrite(T& meshfunction);    
    
    // Most recent position in pvd file
    std::ios::pos_type mark;
    
    // vtu filename
    std::string vtu_filename;

    // pvtu filename
    std::string pvtu_filename;
  };
  
}

#endif
