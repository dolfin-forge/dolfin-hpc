// Copyright (C) 2002-2006 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2002-12-06
// Last changed: 2006-05-23

#include <dolfin/common/types.h>
#include <dolfin/main/MPI.h>
#include <dolfin/log/dolfin_log.h>
#include <dolfin/la/Vector.h>
#include "XMLVector.h"

#ifdef HAS_MPI
#include <mpi.h>
#endif

using namespace dolfin;

//-----------------------------------------------------------------------------
XMLVector::XMLVector(GenericVector& vector)
  : XMLObject(), x(vector), values(0), size(0)
{
  state = OUTSIDE;
}
//-----------------------------------------------------------------------------
void XMLVector::startElement(const xmlChar *name, const xmlChar **attrs)
{
  switch ( state )
  {
  case OUTSIDE:
    
    if ( xmlStrcasecmp(name, (xmlChar *) "vector") == 0 )
    {
      startVector(name, attrs);
      state = INSIDE_VECTOR;
    }
    
    break;
    
  case INSIDE_VECTOR:
    
    if ( xmlStrcasecmp(name, (xmlChar *) "entry") == 0 )
      readEntry(name, attrs);
    
    break;
    
  default:
    ;
  }
}
//-----------------------------------------------------------------------------
void XMLVector::endElement(const xmlChar *name)
{
  switch ( state )
  {
  case INSIDE_VECTOR:
    
    if ( xmlStrcasecmp(name, (xmlChar *) "vector") == 0 )
    {
      endVector();
      state = DONE;
    }
    
    break;
    
  default:
    ;
  }
}
//-----------------------------------------------------------------------------
void XMLVector::startVector(const xmlChar *name, const xmlChar **attrs)
{
  // Parse size of vector
#ifdef HAS_MPI
  uint local_size = parseUnsignedInt(name, attrs, "size");
  MPI_Allreduce(&local_size, &size, 1, MPI_UNSIGNED, MPI_SUM, MPI::DOLFIN_COMM);
#else
  size = parseUnsignedInt(name, attrs, "size");
#endif
  
  // Initialize vector
  if (values)
    delete [] values;
#ifdef HAS_MPI
  values = new real[local_size];
#else
  values = new real[size];
#endif
}
//-----------------------------------------------------------------------------
void XMLVector::endVector()
{
  // Copy values to vector
  dolfin_assert(values);
  x.init(size);
  x.set(values);
  delete [] values;
  values = 0;
}
//-----------------------------------------------------------------------------
void XMLVector::readEntry(const xmlChar *name, const xmlChar **attrs)
{
  // Parse values
  uint row   = parseUnsignedInt(name, attrs, "row");
  real value = parseReal(name, attrs, "value");
  
  // Check values
  if (row >= size)
    error("Illegal XML data for Vector: row index %d out of range (0 - %d)",
          row, size - 1);
  
  // Set value
  dolfin_assert(values);
  values[row] = value;
}
//-----------------------------------------------------------------------------
