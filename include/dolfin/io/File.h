// Copyright (C) 2002-2008 Johan Hoffman and Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Garth N. Wells, 2005, 2006.
// Modified by Magnus Vikstrom 2007
// Modified by Nuno Lopes 2008
// Modified by Niclas Jansson 2008-2012.
//
// First added:  2002-11-12
// Last changed: 2012-05-11

#ifndef __DOLFIN_FILE_H
#define __DOLFIN_FILE_H

#include <dolfin/common/Label.h>

#include <sstream>

namespace dolfin
{

  class Function;
  class Mesh;
  template <class T> class MeshFunction;
  class ParameterList;
  class GenericFile;
  class GenericMatrix;
  class GenericVector;

  /// A File represents a data file for reading and writing objects.
  /// Unless specified explicitly, the format is determined by the
  /// file name suffix.

  class File
  {
  public:

    /// File formats
    enum Type {xml, vtk, binary, off ,raw, stl, xyz};

    /// Create a file with given name
    File(const std::string& filename);

    /// Create a file with given name and type (format)
    File(const std::string& filename, Type type);

    /// Create a file with given name, at time t
    File(const std::string& filename, real const& t);

    /// Destructor
    ~File();

    //--- Input ---

    /// Read vector from file
    void operator>> (GenericVector& x);

    /// Read matrix from file
    void operator>> (GenericMatrix& A);

    /// Read mesh from file
    void operator>> (Mesh& mesh);

    /// Read mesh function from file
    void operator>> (MeshFunction<int>& meshfunction);
    void operator>> (MeshFunction<uint>& meshfunction);
    void operator>> (MeshFunction<real>& meshfunction);
    void operator>> (MeshFunction<bool>& meshfunction);

    /// Read function from file
    void operator>> (Function& u);

    /// Read parameter list from file
    void operator>> (ParameterList& parameters);

    /// Read a collection of funtion to file
    void operator>> (LabelList<Function>& list);

    //--- Output ---

    /// Write vector to file
    void operator<< (GenericVector& x);

    /// Write matrix to file
    void operator<< (GenericMatrix& A);

    /// Write mesh to file
    void operator<< (Mesh& mesh);

    /// Write mesh function to file
    void operator<< (MeshFunction<int>& meshfunction);
    void operator<< (MeshFunction<uint>& meshfunction);
    void operator<< (MeshFunction<real>& meshfunction);
    void operator<< (MeshFunction<bool>& meshfunction);

    /// Write any generic function to file
    void operator<< (Function& u);

    /// Write parameter list to file
    void operator<< (ParameterList& parameters);

    /// Write a collection of functions to file
    void operator<< (LabelList<Function>& f);

    void set_counter(uint new_value);

    //--- STATIC

    static std::string basename(std::string file);

  private:

    GenericFile * file_;

  };

} /* namespace dolfin */

#endif /* __DOLFIN_FILE_H */
