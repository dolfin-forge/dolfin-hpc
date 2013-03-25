// Copyright (C) 2013 Balthasar Reuter.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2013-03-25
// Last changed: 2013-03-25

#ifndef __MESH_FUNCTION_CONVERTER_H
#define __MESH_FUNCTION_CONVERTER_H

#include "MeshFunction.h"
#include <algorithm>

namespace dolfin
{
  class MeshFunctionConverter
  {
  public:
    /// Casts each entry of a MeshFunction to a different datatype and
    /// stores the result in a new MeshFunction. The method is specialized
    /// for casts from double to int or uint to ensure proper rounding insted
    /// of truncation
    ///
    /// *Arguments*
    ///
    ///   source_function (MeshFunction<T1>&)
    ///     The source function
    ///
    ///   target_function (MeshFunction<T2>&)
    ///     The target function
    template<typename T1, typename T2>
    static void cast(MeshFunction<T1>& source_function, 
                     MeshFunction<T2>& target_function)
    {
      // Initialize target function
      target_function.init( source_function.mesh(), source_function.dim() );

      // extract data arrays
      T1 *source_values = source_function.values();
      T2 *target_values = target_function.values();

      // number of values
      uint num_values = source_function.size();

      // copy data
      std::copy( source_values, source_values + num_values, target_values );
    }

  private:
    /// Helper function that performs symmetric rounding to closest integer
    template<typename T1, typename T2>
    static T2 castWithSymmetricRounding(T1 const& value)
    {
      return static_cast<T2>( 
                  (value > 0) ? floor(value + 0.5) : ceil(value - 0.5) 
                );
    }

    /// Helper function to cast floating point values with proper rounding to closest integer
    template<typename RealType, typename IntegerType>
    static void castRealToIntegerWithRounding(MeshFunction<RealType>& source_function, 
                                              MeshFunction<IntegerType>& target_function)
    {
      // Initialize target function
    target_function.init( source_function.mesh(), source_function.dim() );

    // extract data arrays
    RealType *source_values = source_function.values();
    IntegerType *target_values = target_function.values();

    // number of values
    uint num_values = source_function.size();

    // copy data
    std::transform( source_values, source_values + num_values, target_values, castWithSymmetricRounding<RealType,IntegerType> );
    }
  };

  template<>
  void MeshFunctionConverter::cast<double, uint>(MeshFunction<double>& source_function,
                                                 MeshFunction<uint>& target_function)
  {
    castRealToIntegerWithRounding(source_function, target_function);
  }

  template<>
  void MeshFunctionConverter::cast<float, uint>(MeshFunction<float>& source_function,
                                                 MeshFunction<uint>& target_function)
  {
    castRealToIntegerWithRounding(source_function, target_function);
  }

  template<>
  void MeshFunctionConverter::cast<double, int>(MeshFunction<double>& source_function,
                                                 MeshFunction<int>& target_function)
  {
    castRealToIntegerWithRounding(source_function, target_function);
  }

  template<>
  void MeshFunctionConverter::cast<float, int>(MeshFunction<float>& source_function,
                                                 MeshFunction<int>& target_function)
  {
    castRealToIntegerWithRounding(source_function, target_function);
  }
}

#endif