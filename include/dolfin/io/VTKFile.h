// Copyright (C) 2005-2007 Garth N. Wells.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Anders Logg 2006.
// Modified by Niclas Jansson 2008-2010.
//
// First added:  2005-07-05
// Last changed: 2010-06-09

#ifndef __PVTK_FILE_H
#define __PVTK_FILE_H

#include <fstream>
#include <string>
#include <vector>
#include "GenericFile.h"

namespace dolfin
{

  class VTKFile : public GenericFile
  {
  public:
    
    VTKFile(const std::string filename);
    VTKFile(const std::string filename, real& t);
    ~VTKFile();
    
    void operator<< (Mesh& mesh);
    void operator<< (MeshFunction<int>& meshfunction);
    void operator<< (MeshFunction<unsigned int>& meshfunction);
    void operator<< (MeshFunction<double>& meshfunction);
    void operator<< (Function& u);
    void operator<< (std::vector<std::pair<Function*, std::string> >& f);
    
    void write();

    // Compute base64 encoded stream for VTK
    template<typename T>
    void encode_stream(std::stringstream& stream, const std::vector<T>& data) const;

  private:

    void write_dataset(std::vector<std::pair<Function*, std::string> >& f);

    // Compute base64 encoded stream for VTK
    template<typename T>
    void encode_inline_base64(std::stringstream& stream, const std::vector<T>& data) const;

    // Compute compressed base64 encoded stream for VTK
    template<typename T>
    void encode_inline_compressed_base64(std::stringstream& stream, const std::vector<T>& data) const;

    void MeshWrite(Mesh& mesh) const;
    void ResultsWrite(std::vector<std::pair<Function*, std::string> > f) const;
    void pvdFileWrite(uint u);
    void pvtuFileWrite(bool mesh_function = false);
    void pvtuFileWrite_func(std::vector<std::pair<Function*, std::string> > f);
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

    // Current time
    real* _t;
    
  };
  
}

#endif
