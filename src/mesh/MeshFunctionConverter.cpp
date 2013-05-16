// Copyright (C) 2013 Balthasar Reuter.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2013-03-25
// Last changed: 2013-03-25

#include <dolfin/mesh/MeshFunctionConverter.h>

namespace dolfin
{
//-----------------------------------------------------------------------------
template<>
void MeshFunctionConverter::cast<double, uint>(MeshFunction<double>& source_function,
                                               MeshFunction<uint>& target_function)
{
  castRealToIntegerWithRounding(source_function, target_function);
}
//-----------------------------------------------------------------------------
template<>
void MeshFunctionConverter::cast<float, uint>(MeshFunction<float>& source_function,
                                               MeshFunction<uint>& target_function)
{
  castRealToIntegerWithRounding(source_function, target_function);
}
//-----------------------------------------------------------------------------
template<>
void MeshFunctionConverter::cast<double, int>(MeshFunction<double>& source_function,
                                               MeshFunction<int>& target_function)
{
  castRealToIntegerWithRounding(source_function, target_function);
}
//-----------------------------------------------------------------------------
template<>
void MeshFunctionConverter::cast<float, int>(MeshFunction<float>& source_function,
                                               MeshFunction<int>& target_function)
{
  castRealToIntegerWithRounding(source_function, target_function);
}
//-----------------------------------------------------------------------------
template<>
void MeshFunctionConverter::cast<double, bool>(MeshFunction<double>& source_function,
                                               MeshFunction<bool>& target_function)
{
  castRealToIntegerWithRounding(source_function, target_function);
}
//-----------------------------------------------------------------------------
template<>
void MeshFunctionConverter::cast<float, bool>(MeshFunction<float>& source_function,
                                               MeshFunction<bool>& target_function)
{
  castRealToIntegerWithRounding(source_function, target_function);
}
//-----------------------------------------------------------------------------
}