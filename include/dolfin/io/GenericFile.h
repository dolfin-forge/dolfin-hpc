// Copyright (C) 2003-2008 Johan Hoffman and Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Niclas Jansson, 2008-2012.
//
// First added:  2003-07-15
// Last changed: 2012-05-11

#ifndef __DOLFIN_GENERIC_FILE_H
#define __DOLFIN_GENERIC_FILE_H

#include <string>
#include <vector>

#include <dolfin/la/GenericVector.h>
#include <dolfin/la/GenericMatrix.h>

namespace dolfin
{
  
  class Mesh;
  class Graph;
  template <class T> class MeshFunction;
  class Function;
  class ParameterList;
   
  class GenericFile
  {
  public:
    
    GenericFile(const std::string filename);
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
    virtual void operator>> (Graph& graph);
    virtual void operator>> (std::vector<std::pair<Function*, std::string> >& f);
    
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
    virtual void operator<< (Graph& graph);
    virtual void operator<< (std::vector<std::pair<Function*, std::string> >& f);
    
    void set_counter(uint value);

    virtual void read();
    virtual void write();
    
  protected:
    
    void read_not_impl(const std::string object);
    void write_not_impl(const std::string object);
    void parallel_read_not_impl(const std::string object);
    void parallel_write_not_impl(const std::string object);

    std::string filename;
    std::string type;
    
    bool opened_read;
    bool opened_write;

    bool check_header; // True if we have written a header

    // Counters for the number of times various data has been written
    uint counter;
    uint counter1;
    uint counter2;

  };
  
}

#endif
