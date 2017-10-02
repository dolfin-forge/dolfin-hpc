// Copyright (C) 2003-2008 Johan Hoffman and Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Niclas Jansson, 2008-2012.
//
// First added:  2003-07-15
// Last changed: 2012-05-11

#ifndef __DOLFIN_GENERIC_FILE_H
#define __DOLFIN_GENERIC_FILE_H

#include <dolfin/common/Label.h>

#include <vector>

namespace dolfin
{
  
  class Function;
  class GenericMatrix;
  class GenericVector;
  class Mesh;
  template <class T> class MeshFunction;
  class ParameterList;
   
  class GenericFile
  {
  public:
    
    GenericFile(std::string const& type, std::string const& filename);
    virtual ~GenericFile();
    
    // Input

    virtual void operator>> (GenericVector& x);
    virtual void operator>> (GenericMatrix& A);
    virtual void operator>> (Mesh& mesh);
    virtual void operator>> (MeshFunction<int>& meshfunction);
    virtual void operator>> (MeshFunction<uint>& meshfunction);
    virtual void operator>> (MeshFunction<real>& meshfunction);
    virtual void operator>> (MeshFunction<bool>& meshfunction);
    virtual void operator>> (Function& mesh);
    virtual void operator>> (ParameterList& parameters);
    virtual void operator>> (LabelList<Function>& list);
    
    // Output
    
    virtual void operator<< (GenericVector& x);
    virtual void operator<< (GenericMatrix& A);
    virtual void operator<< (Mesh& mesh);
    virtual void operator<< (MeshFunction<int>& meshfunction);
    virtual void operator<< (MeshFunction<uint>& meshfunction);
    virtual void operator<< (MeshFunction<real>& meshfunction);
    virtual void operator<< (MeshFunction<bool>& meshfunction);
    virtual void operator<< (Function& u);
    virtual void operator<< (ParameterList& parameters);
    virtual void operator<< (LabelList<Function>& f);
    
    void set_counter(uint value);

    virtual void read();
    virtual void write();
    
  private:

    std::string const type_;

  protected:

    void read_not_impl(const std::string object);
    void write_not_impl(const std::string object);
    void parallel_read_not_impl(const std::string object);
    void parallel_write_not_impl(const std::string object);

    std::string filename;
    
    bool opened_read;
    bool opened_write;

    bool check_header; // True if we have written a header

    // Counters for the number of times various data has been written
    uint counter;

  };
  
}

#endif
