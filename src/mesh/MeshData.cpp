// Copyright (C) 2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2008-05-19
// Last changed: 2008-05-28

#include <dolfin/mesh/MeshData.h>

using namespace dolfin;

typedef std::map<std::string, MeshFunction<dolfin::uint>*>::iterator mf_iterator;
typedef std::map<std::string, MeshFunction<dolfin::uint>*>::const_iterator mf_const_iterator;

typedef std::map<std::string, Array<dolfin::uint>*>::iterator a_iterator;
typedef std::map<std::string, Array<dolfin::uint>*>::const_iterator a_const_iterator;
/*
typedef std::map<std::string, MeshFunction<int>*>::iterator mf_iterator_int;
typedef std::map<std::string, MeshFunction<int>*>::const_iterator mf_const_iterator_int;

typedef std::map<std::string, Array<int>*>::iterator a_iterator_int;
typedef std::map<std::string, Array<int>*>::const_iterator a_const_iterator_int;

typedef std::map<std::string, MeshFunction<dolfin::real>*>::iterator mf_iterator_real;
typedef std::map<std::string, MeshFunction<dolfin::real>*>::const_iterator mf_const_iterator_real;

typedef std::map<std::string, Array<dolfin::real>*>::iterator a_iterator_real;
typedef std::map<std::string, Array<dolfin::real>*>::const_iterator a_const_iterator_real;
*/
//-----------------------------------------------------------------------------
MeshData::MeshData(Mesh& mesh) : mesh(mesh)
{
  // Do nothing
}
//-----------------------------------------------------------------------------
MeshData::~MeshData()
{
  clear();
}
//-----------------------------------------------------------------------------
void MeshData::clear()
{ 
  for (mf_iterator it = meshfunctions.begin(); it != meshfunctions.end(); ++it)
    delete it->second;
  meshfunctions.clear();

  for (a_iterator it = arrays.begin(); it != arrays.end(); ++it)
    delete it->second;
  arrays.clear();
/*
  for (mf_iterator_int it = meshfunctions_int.begin(); it != meshfunctions_int.end(); ++it)
    delete it->second;
  meshfunctions_int.clear();

  for (a_iterator_int it = arrays_int.begin(); it != arrays_int.end(); ++it)
    delete it->second;
  arrays_int.clear();

  for (mf_iterator_real it = meshfunctions_real.begin(); it != meshfunctions_real.end(); ++it)
    delete it->second;
  meshfunctions_real.clear();

  for (a_iterator_real it = arrays_real.begin(); it != arrays_real.end(); ++it)
    delete it->second;
  arrays_real.clear();
*/
}
//-----------------------------------------------------------------------------
MeshFunction<dolfin::uint>* MeshData::createMeshFunction(std::string name)
{
  // Check if data already exists
  mf_iterator it = meshfunctions.find(name);
  if (it != meshfunctions.end())
  {
    warning("Mesh data named \"%s\" already exists.", name.c_str());
    return it->second;
  }

  // Create new data
  MeshFunction<uint>* f = new MeshFunction<uint>(mesh);
  dolfin_assert(f);

  // Add to map
  meshfunctions[name] = f;

  return f;
}
//-----------------------------------------------------------------------------
Array<dolfin::uint>* MeshData::createArray(std::string name, uint size)
{
  // Check if data already exists
  a_iterator it = arrays.find(name);
  if (it != arrays.end())
  {
    warning("Mesh data named \"%s\" already exists.", name.c_str());
    return it->second;
  }

  // Create new data
  Array<uint>* a = new Array<uint>(size);
  *a = 0;

  // Add to map
  arrays[name] = a;

  return a;
}
//-----------------------------------------------------------------------------
MeshFunction<dolfin::uint>* MeshData::meshFunction(std::string name)
{
  // Check if data exists
  mf_iterator it = meshfunctions.find(name);
  if (it == meshfunctions.end())
    return 0;
  
  return it->second;
}
//-----------------------------------------------------------------------------
Array<dolfin::uint>* MeshData::array(std::string name)
{
  // Check if data exists
  a_iterator it = arrays.find(name);
  if (it == arrays.end())
    return 0;
  
  return it->second;
}
//-----------------------------------------------------------------------------
/*MeshFunction<int>* MeshData::createMeshFunctionInt(std::string name)
{
  // Check if data already exists
  mf_iterator_int it = meshfunctions_int.find(name);
  if (it != meshfunctions_int.end())
  {
    warning("Mesh data named \"%s\" already exists.", name.c_str());
    return it->second;
  }

  // Create new data
  MeshFunction<int>* f = new MeshFunction<int>(mesh);
  dolfin_assert(f);

  // Add to map
  meshfunctions_int[name] = f;

  return f;
}
//-----------------------------------------------------------------------------
Array<int>* MeshData::createArrayInt(std::string name, uint size)
{
  // Check if data already exists
  a_iterator_int it = arrays_int.find(name);
  if (it != arrays_int.end())
  {
    warning("Mesh data named \"%s\" already exists.", name.c_str());
    return it->second;
  }

  // Create new data
  Array<int>* a = new Array<int>(size);
  *a = 0;

  // Add to map
  arrays_int[name] = a;

  return a;
}
//-----------------------------------------------------------------------------
MeshFunction<int>* MeshData::meshFunctionInt(std::string name)
{
  // Check if data exists
  mf_iterator_int it = meshfunctions_int.find(name);
  if (it == meshfunctions_int.end())
    return 0;
  
  return it->second;
}
//-----------------------------------------------------------------------------
Array<int>* MeshData::arrayInt(std::string name)
{
  // Check if data exists
  a_iterator_int it = arrays_int.find(name);
  if (it == arrays_int.end())
    return 0;
  
  return it->second;
}
//-----------------------------------------------------------------------------
MeshFunction<dolfin::real>* MeshData::createMeshFunctionReal(std::string name)
{
  // Check if data already exists
  mf_iterator_real it = meshfunctions_real.find(name);
  if (it != meshfunctions_real.end())
  {
    warning("Mesh data named \"%s\" already exists.", name.c_str());
    return it->second;
  }

  // Create new data
  MeshFunction<real>* f = new MeshFunction<real>(mesh);
  dolfin_assert(f);

  // Add to map
  meshfunctions_real[name] = f;

  return f;
}
//-----------------------------------------------------------------------------
Array<dolfin::real>* MeshData::createArrayReal(std::string name, uint size)
{
  // Check if data already exists
  a_iterator_real it = arrays_real.find(name);
  if (it != arrays_real.end())
  {
    warning("Mesh data named \"%s\" already exists.", name.c_str());
    return it->second;
  }

  // Create new data
  Array<real>* a = new Array<real>(size);
  *a = 0;

  // Add to map
  arrays_real[name] = a;

  return a;
}
//-----------------------------------------------------------------------------
MeshFunction<dolfin::real>* MeshData::meshFunctionReal(std::string name)
{
  // Check if data exists
  mf_iterator_real it = meshfunctions_real.find(name);
  if (it == meshfunctions_real.end())
    return 0;
  
  return it->second;
}
//-----------------------------------------------------------------------------
Array<dolfin::real>* MeshData::arrayReal(std::string name)
{
  // Check if data exists
  a_iterator_real it = arrays_real.find(name);
  if (it == arrays_real.end())
    return 0;
  
  return it->second;
}*/
//-----------------------------------------------------------------------------
void MeshData::disp() const
{
  // Begin indentation
  cout << "Auxiliary mesh data" << endl;
  begin("-------------------");
  cout << endl;

  for (mf_const_iterator it = meshfunctions.begin(); it != meshfunctions.end(); ++it)
  {
    cout << "MeshFunction<uint> of size "
         << it->second->size()
         << " on entities of topological dimension "
         << it->second->dim()
         << ": \"" << it->first << "\"" << endl;
  }

  for (a_const_iterator it = arrays.begin(); it != arrays.end(); ++it)
    cout << "Array<uint> of size " << static_cast<uint>(it->second->size())
         << ": \"" << it->first << "\"" << endl;
/*
	for (mf_const_iterator_int it = meshfunctions_int.begin(); it != meshfunctions_int.end(); ++it)
  {
    cout << "MeshFunction<int> of size "
         << it->second->size()
         << " on entities of topological dimension "
         << it->second->dim()
         << ": \"" << it->first << "\"" << endl;
  }

  for (a_const_iterator_int it = arrays_int.begin(); it != arrays_int.end(); ++it)
    cout << "Array<int> of size " << static_cast<uint>(it->second->size())
         << ": \"" << it->first << "\"" << endl;

	for (mf_const_iterator_real it = meshfunctions_real.begin(); it != meshfunctions_real.end(); ++it)
  {
    cout << "MeshFunction<real> of size "
         << it->second->size()
         << " on entities of topological dimension "
         << it->second->dim()
         << ": \"" << it->first << "\"" << endl;
  }

  for (a_const_iterator_real it = arrays_real.begin(); it != arrays_real.end(); ++it)
    cout << "Array<real> of size " << static_cast<uint>(it->second->size())
         << ": \"" << it->first << "\"" << endl;
*/
  // End indentation
  end();
}
//-----------------------------------------------------------------------------
