// Automatically generated code mapping element and dof map signatures
// to the corresponding ufc::finite_element and ufc::dof_map classes

#include <dolfin/config/dolfin_config.h>
#include <cstring>
#include <dolfin/fem/DofMap.h>

#ifdef ENABLE_UFL		// UFC 2.1.1 headers

#include <dolfin/elements/ffc_Lagrange_1_1d.h>
#include <dolfin/elements/ffc_Lagrange_2_1d.h>
#include <dolfin/elements/ffc_Lagrange_1_2d.h>
#include <dolfin/elements/ffc_Lagrange_2_2d.h>
#include <dolfin/elements/ffc_Lagrange_1_3d.h>
#include <dolfin/elements/ffc_Lagrange_2_3d.h>
#include <dolfin/elements/ffc_Discontinuous_Lagrange_0_1d.h>
#include <dolfin/elements/ffc_Discontinuous_Lagrange_1_1d.h>
#include <dolfin/elements/ffc_Discontinuous_Lagrange_2_1d.h>
#include <dolfin/elements/ffc_Discontinuous_Lagrange_0_2d.h>
#include <dolfin/elements/ffc_Discontinuous_Lagrange_1_2d.h>
#include <dolfin/elements/ffc_Discontinuous_Lagrange_2_2d.h>
#include <dolfin/elements/ffc_Discontinuous_Lagrange_0_3d.h>
#include <dolfin/elements/ffc_Discontinuous_Lagrange_1_3d.h>
#include <dolfin/elements/ffc_Discontinuous_Lagrange_2_3d.h>
#include <dolfin/elements/ffc_Lagrange_1_2dVector.h>
#include <dolfin/elements/ffc_Lagrange_2_2dVector.h>
#include <dolfin/elements/ffc_Lagrange_1_3dVector.h>
#include <dolfin/elements/ffc_Lagrange_2_3dVector.h>
#include <dolfin/elements/ffc_Discontinuous_Lagrange_0_2dVector.h>
#include <dolfin/elements/ffc_Discontinuous_Lagrange_1_2dVector.h>
#include <dolfin/elements/ffc_Discontinuous_Lagrange_2_2dVector.h>
#include <dolfin/elements/ffc_Discontinuous_Lagrange_0_3dVector.h>
#include <dolfin/elements/ffc_Discontinuous_Lagrange_1_3dVector.h>
#include <dolfin/elements/ffc_Discontinuous_Lagrange_2_3dVector.h>
#include <dolfin/elements/ffc_Brezzi_Douglas_Marini_1_2d.h>

#else  // UFC 1.1 headers

#ifdef ENABLE_UFL		// UFC 2.1.1 headers

#include <dolfin/elements/ffc_Lagrange_1_1d.h>
#include <dolfin/elements/ffc_Lagrange_2_1d.h>
#include <dolfin/elements/ffc_Lagrange_1_2d.h>
#include <dolfin/elements/ffc_Lagrange_2_2d.h>
#include <dolfin/elements/ffc_Lagrange_1_3d.h>
#include <dolfin/elements/ffc_Lagrange_2_3d.h>
#include <dolfin/elements/ffc_Discontinuous_Lagrange_0_1d.h>
#include <dolfin/elements/ffc_Discontinuous_Lagrange_1_1d.h>
#include <dolfin/elements/ffc_Discontinuous_Lagrange_2_1d.h>
#include <dolfin/elements/ffc_Discontinuous_Lagrange_0_2d.h>
#include <dolfin/elements/ffc_Discontinuous_Lagrange_1_2d.h>
#include <dolfin/elements/ffc_Discontinuous_Lagrange_2_2d.h>
#include <dolfin/elements/ffc_Discontinuous_Lagrange_0_3d.h>
#include <dolfin/elements/ffc_Discontinuous_Lagrange_1_3d.h>
#include <dolfin/elements/ffc_Discontinuous_Lagrange_2_3d.h>
#include <dolfin/elements/ffc_Lagrange_1_2dVector.h>
#include <dolfin/elements/ffc_Lagrange_2_2dVector.h>
#include <dolfin/elements/ffc_Lagrange_1_3dVector.h>
#include <dolfin/elements/ffc_Lagrange_2_3dVector.h>
#include <dolfin/elements/ffc_Discontinuous_Lagrange_0_2dVector.h>
#include <dolfin/elements/ffc_Discontinuous_Lagrange_1_2dVector.h>
#include <dolfin/elements/ffc_Discontinuous_Lagrange_2_2dVector.h>
#include <dolfin/elements/ffc_Discontinuous_Lagrange_0_3dVector.h>
#include <dolfin/elements/ffc_Discontinuous_Lagrange_1_3dVector.h>
#include <dolfin/elements/ffc_Discontinuous_Lagrange_2_3dVector.h>
#include <dolfin/elements/ffc_Brezzi_Douglas_Marini_1_2d.h>

#else  // UFC 1.1 headers

#include <dolfin/elements/ffc_00.h>
#include <dolfin/elements/ffc_01.h>
#include <dolfin/elements/ffc_02.h>
#include <dolfin/elements/ffc_03.h>
#include <dolfin/elements/ffc_04.h>
#include <dolfin/elements/ffc_05.h>
#include <dolfin/elements/ffc_06.h>
#include <dolfin/elements/ffc_07.h>
#include <dolfin/elements/ffc_08.h>
#include <dolfin/elements/ffc_09.h>
#include <dolfin/elements/ffc_10.h>
#include <dolfin/elements/ffc_11.h>
#include <dolfin/elements/ffc_12.h>
#include <dolfin/elements/ffc_13.h>
#include <dolfin/elements/ffc_14.h>
#include <dolfin/elements/ffc_15.h>
#include <dolfin/elements/ffc_16.h>
#include <dolfin/elements/ffc_17.h>
#include <dolfin/elements/ffc_18.h>
#include <dolfin/elements/ffc_19.h>
#include <dolfin/elements/ffc_20.h>
#include <dolfin/elements/ffc_21.h>
#include <dolfin/elements/ffc_22.h>
#include <dolfin/elements/ffc_23.h>
#include <dolfin/elements/ffc_24.h>

#endif

#include <dolfin/elements/ElementLibrary.h>

#ifdef ENABLE_UFL

ufc::finite_element* dolfin::ElementLibrary::create_finite_element(const char* signature)
{
  if (strcmp(signature, "FiniteElement('Lagrange', Cell('interval', Space(1)), 1, None)") == 0)
    return new ffc_lagrange_1_1d_finite_element_0();
  if (strcmp(signature, "FiniteElement('Lagrange', Cell('interval', Space(1)), 2, None)") == 0)
    return new ffc_lagrange_2_1d_finite_element_0();
  if (strcmp(signature, "FiniteElement('Lagrange', Cell('triangle', Space(2)), 1, None)") == 0)
    return new ffc_lagrange_1_2d_finite_element_0();
  if (strcmp(signature, "FiniteElement('Lagrange', Cell('triangle', Space(2)), 2, None)") == 0)
    return new ffc_lagrange_2_2d_finite_element_0();
  if (strcmp(signature, "FiniteElement('Lagrange', Cell('tetrahedron', Space(3)), 1, None)") == 0)
    return new ffc_lagrange_1_3d_finite_element_0();
  if (strcmp(signature, "FiniteElement('Lagrange', Cell('tetrahedron', Space(3)), 2, None)") == 0)
    return new ffc_lagrange_2_3d_finite_element_0();
  if (strcmp(signature, "FiniteElement('Discontinuous Lagrange', Cell('interval', Space(1)), 0, None)") == 0)
    return new ffc_discontinuous_lagrange_0_1d_finite_element_0();
  if (strcmp(signature, "FiniteElement('Discontinuous Lagrange', Cell('interval', Space(1)), 1, None)") == 0)
    return new ffc_discontinuous_lagrange_1_1d_finite_element_0();
  if (strcmp(signature, "FiniteElement('Discontinuous Lagrange', Cell('interval', Space(1)), 2, None)") == 0)
    return new ffc_discontinuous_lagrange_2_1d_finite_element_0();
  if (strcmp(signature, "FiniteElement('Discontinuous Lagrange', Cell('triangle', Space(2)), 0, None)") == 0)
    return new ffc_discontinuous_lagrange_0_2d_finite_element_0();
  if (strcmp(signature, "FiniteElement('Discontinuous Lagrange', Cell('triangle', Space(2)), 1, None)") == 0)
    return new ffc_discontinuous_lagrange_1_2d_finite_element_0();
  if (strcmp(signature, "FiniteElement('Discontinuous Lagrange', Cell('triangle', Space(2)), 2, None)") == 0)
    return new ffc_discontinuous_lagrange_2_2d_finite_element_0();
  if (strcmp(signature, "FiniteElement('Discontinuous Lagrange', Cell('tetrahedron', Space(3)), 0, None)") == 0)
    return new ffc_discontinuous_lagrange_0_3d_finite_element_0();
  if (strcmp(signature, "FiniteElement('Discontinuous Lagrange', Cell('tetrahedron', Space(3)), 1, None)") == 0)
    return new ffc_discontinuous_lagrange_1_3d_finite_element_0();
  if (strcmp(signature, "FiniteElement('Discontinuous Lagrange', Cell('tetrahedron', Space(3)), 2, None)") == 0)
    return new ffc_discontinuous_lagrange_2_3d_finite_element_0();
  if (strcmp(signature, "VectorElement('Lagrange', Cell('triangle', Space(2)), 1, 2, None)") == 0)
    return new ffc_lagrange_1_2dvector_finite_element_1();
  if (strcmp(signature, "VectorElement('Lagrange', Cell('triangle', Space(2)), 2, 2, None)") == 0)
    return new ffc_lagrange_2_2dvector_finite_element_1();
  if (strcmp(signature, "VectorElement('Lagrange', Cell('tetrahedron', Space(3)), 1, 3, None)") == 0)
    return new ffc_lagrange_1_3dvector_finite_element_1();
  if (strcmp(signature, "VectorElement('Lagrange', Cell('tetrahedron', Space(3)), 2, 3, None)") == 0)
    return new ffc_lagrange_2_3dvector_finite_element_1();
  if (strcmp(signature, "VectorElement('Discontinuous Lagrange', Cell('triangle', Space(2)), 0, 2, None)") == 0)
    return new ffc_discontinuous_lagrange_0_2dvector_finite_element_1();
  if (strcmp(signature, "VectorElement('Discontinuous Lagrange', Cell('triangle', Space(2)), 1, 2, None)") == 0)
    return new ffc_discontinuous_lagrange_1_2dvector_finite_element_1();
  if (strcmp(signature, "VectorElement('Discontinuous Lagrange', Cell('triangle', Space(2)), 2, 2, None)") == 0)
    return new ffc_discontinuous_lagrange_2_2dvector_finite_element_1();
  if (strcmp(signature, "VectorElement('Discontinuous Lagrange', Cell('tetrahedron', Space(3)), 0, 3, None)") == 0)
    return new ffc_discontinuous_lagrange_0_3dvector_finite_element_1();
  if (strcmp(signature, "VectorElement('Discontinuous Lagrange', Cell('tetrahedron', Space(3)), 1, 3, None)") == 0)
    return new ffc_discontinuous_lagrange_1_3dvector_finite_element_1();
  if (strcmp(signature, "VectorElement('Discontinuous Lagrange', Cell('tetrahedron', Space(3)), 2, 3, None)") == 0)
    return new ffc_discontinuous_lagrange_2_3dvector_finite_element_1();
  if (strcmp(signature, "FiniteElement('Brezzi-Douglas-Marini', Cell('triangle', Space(2)), 1, None)") == 0)
    return new ffc_brezzi_douglas_marini_1_2d_finite_element_0();
  return 0;
}

ufc::dof_map* dolfin::ElementLibrary::create_dof_map(const char* signature)
{
  if (strcmp(signature, DofMap::dofmap_signature("FiniteElement('Lagrange', Cell('interval', Space(1)), 1, None)").c_str()) == 0)
    return new ffc_lagrange_1_1d_dofmap_0();
  if (strcmp(signature, DofMap::dofmap_signature("FiniteElement('Lagrange', Cell('interval', Space(1)), 2, None)").c_str()) == 0)
    return new ffc_lagrange_2_1d_dofmap_0();
  if (strcmp(signature, DofMap::dofmap_signature("FiniteElement('Lagrange', Cell('triangle', Space(2)), 1, None)").c_str()) == 0)
    return new ffc_lagrange_1_2d_dofmap_0();
  if (strcmp(signature, DofMap::dofmap_signature("FiniteElement('Lagrange', Cell('triangle', Space(2)), 2, None)").c_str()) == 0)
    return new ffc_lagrange_2_2d_dofmap_0();
  if (strcmp(signature, DofMap::dofmap_signature("FiniteElement('Lagrange', Cell('tetrahedron', Space(3)), 1, None)").c_str()) == 0)
    return new ffc_lagrange_1_3d_dofmap_0();
  if (strcmp(signature, DofMap::dofmap_signature("FiniteElement('Lagrange', Cell('tetrahedron', Space(3)), 2, None)").c_str()) == 0)
    return new ffc_lagrange_2_3d_dofmap_0();
  if (strcmp(signature, DofMap::dofmap_signature("FiniteElement('Discontinuous Lagrange', Cell('interval', Space(1)), 0, None)").c_str()) == 0)
    return new ffc_discontinuous_lagrange_0_1d_dofmap_0();
  if (strcmp(signature, DofMap::dofmap_signature("FiniteElement('Discontinuous Lagrange', Cell('interval', Space(1)), 1, None)").c_str()) == 0)
    return new ffc_discontinuous_lagrange_1_1d_dofmap_0();
  if (strcmp(signature, DofMap::dofmap_signature("FiniteElement('Discontinuous Lagrange', Cell('interval', Space(1)), 2, None)").c_str()) == 0)
    return new ffc_discontinuous_lagrange_2_1d_dofmap_0();
  if (strcmp(signature, DofMap::dofmap_signature("FiniteElement('Discontinuous Lagrange', Cell('triangle', Space(2)), 0, None)").c_str()) == 0)
    return new ffc_discontinuous_lagrange_0_2d_dofmap_0();
  if (strcmp(signature, DofMap::dofmap_signature("FiniteElement('Discontinuous Lagrange', Cell('triangle', Space(2)), 1, None)").c_str()) == 0)
    return new ffc_discontinuous_lagrange_1_2d_dofmap_0();
  if (strcmp(signature, DofMap::dofmap_signature("FiniteElement('Discontinuous Lagrange', Cell('triangle', Space(2)), 2, None)").c_str()) == 0)
    return new ffc_discontinuous_lagrange_2_2d_dofmap_0();
  if (strcmp(signature, DofMap::dofmap_signature("FiniteElement('Discontinuous Lagrange', Cell('tetrahedron', Space(3)), 0, None)").c_str()) == 0)
    return new ffc_discontinuous_lagrange_0_3d_dofmap_0();
  if (strcmp(signature, DofMap::dofmap_signature("FiniteElement('Discontinuous Lagrange', Cell('tetrahedron', Space(3)), 1, None)").c_str()) == 0)
    return new ffc_discontinuous_lagrange_1_3d_dofmap_0();
  if (strcmp(signature, DofMap::dofmap_signature("FiniteElement('Discontinuous Lagrange', Cell('tetrahedron', Space(3)), 2, None)").c_str()) == 0)
    return new ffc_discontinuous_lagrange_2_3d_dofmap_0();
  if (strcmp(signature, DofMap::dofmap_signature("VectorElement('Lagrange', Cell('triangle', Space(2)), 1, 2, None)").c_str()) == 0)
    return new ffc_lagrange_1_2dvector_dofmap_1();
  if (strcmp(signature, DofMap::dofmap_signature("VectorElement('Lagrange', Cell('triangle', Space(2)), 2, 2, None)").c_str()) == 0)
    return new ffc_lagrange_2_2dvector_dofmap_1();
  if (strcmp(signature, DofMap::dofmap_signature("VectorElement('Lagrange', Cell('tetrahedron', Space(3)), 1, 3, None)").c_str()) == 0)
    return new ffc_lagrange_1_3dvector_dofmap_1();
  if (strcmp(signature, DofMap::dofmap_signature("VectorElement('Lagrange', Cell('tetrahedron', Space(3)), 2, 3, None)").c_str()) == 0)
    return new ffc_lagrange_2_3dvector_dofmap_1();
  if (strcmp(signature, DofMap::dofmap_signature("VectorElement('Discontinuous Lagrange', Cell('triangle', Space(2)), 0, 2, None)").c_str()) == 0)
    return new ffc_discontinuous_lagrange_0_2dvector_dofmap_1();
  if (strcmp(signature, DofMap::dofmap_signature("VectorElement('Discontinuous Lagrange', Cell('triangle', Space(2)), 1, 2, None)").c_str()) == 0)
    return new ffc_discontinuous_lagrange_1_2dvector_dofmap_1();
  if (strcmp(signature, DofMap::dofmap_signature("VectorElement('Discontinuous Lagrange', Cell('triangle', Space(2)), 2, 2, None)").c_str()) == 0)
    return new ffc_discontinuous_lagrange_2_2dvector_dofmap_1();
  if (strcmp(signature, DofMap::dofmap_signature("VectorElement('Discontinuous Lagrange', Cell('tetrahedron', Space(3)), 0, 3, None)").c_str()) == 0)
    return new ffc_discontinuous_lagrange_0_3dvector_dofmap_1();
  if (strcmp(signature, DofMap::dofmap_signature("VectorElement('Discontinuous Lagrange', Cell('tetrahedron', Space(3)), 1, 3, None)").c_str()) == 0)
    return new ffc_discontinuous_lagrange_1_3dvector_dofmap_1();
  if (strcmp(signature, DofMap::dofmap_signature("VectorElement('Discontinuous Lagrange', Cell('tetrahedron', Space(3)), 2, 3, None)").c_str()) == 0)
    return new ffc_discontinuous_lagrange_2_3dvector_dofmap_1();
  if (strcmp(signature, DofMap::dofmap_signature("FiniteElement('Brezzi-Douglas-Marini', Cell('triangle', Space(2)), 1, None)").c_str()) == 0)
    return new ffc_brezzi_douglas_marini_1_2d_dofmap_0();
  return 0;
}

#else

ufc::finite_element* dolfin::ElementLibrary::create_finite_element(const char* signature)
{
  if (strcmp(signature, "Lagrange finite element of degree 1 on a interval") == 0)
    return new ffc_00_finite_element_0();
  if (strcmp(signature, "Lagrange finite element of degree 2 on a interval") == 0)
    return new ffc_01_finite_element_0();
  if (strcmp(signature, "Lagrange finite element of degree 1 on a triangle") == 0)
    return new ffc_02_finite_element_0();
  if (strcmp(signature, "Lagrange finite element of degree 2 on a triangle") == 0)
    return new ffc_03_finite_element_0();
  if (strcmp(signature, "Lagrange finite element of degree 1 on a tetrahedron") == 0)
    return new ffc_04_finite_element_0();
  if (strcmp(signature, "Lagrange finite element of degree 2 on a tetrahedron") == 0)
    return new ffc_05_finite_element_0();
  if (strcmp(signature, "Discontinuous Lagrange finite element of degree 0 on a interval") == 0)
    return new ffc_06_finite_element_0();
  if (strcmp(signature, "Discontinuous Lagrange finite element of degree 1 on a interval") == 0)
    return new ffc_07_finite_element_0();
  if (strcmp(signature, "Discontinuous Lagrange finite element of degree 2 on a interval") == 0)
    return new ffc_08_finite_element_0();
  if (strcmp(signature, "Discontinuous Lagrange finite element of degree 0 on a triangle") == 0)
    return new ffc_09_finite_element_0();
  if (strcmp(signature, "Discontinuous Lagrange finite element of degree 1 on a triangle") == 0)
    return new ffc_10_finite_element_0();
  if (strcmp(signature, "Discontinuous Lagrange finite element of degree 2 on a triangle") == 0)
    return new ffc_11_finite_element_0();
  if (strcmp(signature, "Discontinuous Lagrange finite element of degree 0 on a tetrahedron") == 0)
    return new ffc_12_finite_element_0();
  if (strcmp(signature, "Discontinuous Lagrange finite element of degree 1 on a tetrahedron") == 0)
    return new ffc_13_finite_element_0();
  if (strcmp(signature, "Discontinuous Lagrange finite element of degree 2 on a tetrahedron") == 0)
    return new ffc_14_finite_element_0();
  if (strcmp(signature, "Mixed finite element: [Lagrange finite element of degree 1 on a triangle, Lagrange finite element of degree 1 on a triangle]") == 0)
    return new ffc_15_finite_element_0();
  if (strcmp(signature, "Mixed finite element: [Lagrange finite element of degree 2 on a triangle, Lagrange finite element of degree 2 on a triangle]") == 0)
    return new ffc_16_finite_element_0();
  if (strcmp(signature, "Mixed finite element: [Lagrange finite element of degree 1 on a tetrahedron, Lagrange finite element of degree 1 on a tetrahedron, Lagrange finite element of degree 1 on a tetrahedron]") == 0)
    return new ffc_17_finite_element_0();
  if (strcmp(signature, "Mixed finite element: [Lagrange finite element of degree 2 on a tetrahedron, Lagrange finite element of degree 2 on a tetrahedron, Lagrange finite element of degree 2 on a tetrahedron]") == 0)
    return new ffc_18_finite_element_0();
  if (strcmp(signature, "Mixed finite element: [Discontinuous Lagrange finite element of degree 0 on a triangle, Discontinuous Lagrange finite element of degree 0 on a triangle]") == 0)
    return new ffc_19_finite_element_0();
  if (strcmp(signature, "Mixed finite element: [Discontinuous Lagrange finite element of degree 1 on a triangle, Discontinuous Lagrange finite element of degree 1 on a triangle]") == 0)
    return new ffc_20_finite_element_0();
  if (strcmp(signature, "Mixed finite element: [Discontinuous Lagrange finite element of degree 0 on a tetrahedron, Discontinuous Lagrange finite element of degree 0 on a tetrahedron, Discontinuous Lagrange finite element of degree 0 on a tetrahedron]") == 0)
    return new ffc_21_finite_element_0();
  if (strcmp(signature, "Mixed finite element: [Discontinuous Lagrange finite element of degree 1 on a tetrahedron, Discontinuous Lagrange finite element of degree 1 on a tetrahedron, Discontinuous Lagrange finite element of degree 1 on a tetrahedron]") == 0)
    return new ffc_22_finite_element_0();
  if (strcmp(signature, "Mixed finite element: [Discontinuous Lagrange finite element of degree 2 on a tetrahedron, Discontinuous Lagrange finite element of degree 2 on a tetrahedron, Discontinuous Lagrange finite element of degree 2 on a tetrahedron]") == 0)
    return new ffc_23_finite_element_0();
  if (strcmp(signature, "Brezzi-Douglas-Marini finite element of degree 1 on a triangle") == 0)
    return new ffc_24_finite_element_0();
  return 0;
}

ufc::dof_map* dolfin::ElementLibrary::create_dof_map(const char* signature)
{
  if (strcmp(signature, "FFC dof map for Lagrange finite element of degree 1 on a interval") == 0)
    return new ffc_00_dof_map_0();
  if (strcmp(signature, "FFC dof map for Lagrange finite element of degree 2 on a interval") == 0)
    return new ffc_01_dof_map_0();
  if (strcmp(signature, "FFC dof map for Lagrange finite element of degree 1 on a triangle") == 0)
    return new ffc_02_dof_map_0();
  if (strcmp(signature, "FFC dof map for Lagrange finite element of degree 2 on a triangle") == 0)
    return new ffc_03_dof_map_0();
  if (strcmp(signature, "FFC dof map for Lagrange finite element of degree 1 on a tetrahedron") == 0)
    return new ffc_04_dof_map_0();
  if (strcmp(signature, "FFC dof map for Lagrange finite element of degree 2 on a tetrahedron") == 0)
    return new ffc_05_dof_map_0();
  if (strcmp(signature, "FFC dof map for Discontinuous Lagrange finite element of degree 0 on a interval") == 0)
    return new ffc_06_dof_map_0();
  if (strcmp(signature, "FFC dof map for Discontinuous Lagrange finite element of degree 1 on a interval") == 0)
    return new ffc_07_dof_map_0();
  if (strcmp(signature, "FFC dof map for Discontinuous Lagrange finite element of degree 2 on a interval") == 0)
    return new ffc_08_dof_map_0();
  if (strcmp(signature, "FFC dof map for Discontinuous Lagrange finite element of degree 0 on a triangle") == 0)
    return new ffc_09_dof_map_0();
  if (strcmp(signature, "FFC dof map for Discontinuous Lagrange finite element of degree 1 on a triangle") == 0)
    return new ffc_10_dof_map_0();
  if (strcmp(signature, "FFC dof map for Discontinuous Lagrange finite element of degree 2 on a triangle") == 0)
    return new ffc_11_dof_map_0();
  if (strcmp(signature, "FFC dof map for Discontinuous Lagrange finite element of degree 0 on a tetrahedron") == 0)
    return new ffc_12_dof_map_0();
  if (strcmp(signature, "FFC dof map for Discontinuous Lagrange finite element of degree 1 on a tetrahedron") == 0)
    return new ffc_13_dof_map_0();
  if (strcmp(signature, "FFC dof map for Discontinuous Lagrange finite element of degree 2 on a tetrahedron") == 0)
    return new ffc_14_dof_map_0();
  if (strcmp(signature, "FFC dof map for Mixed finite element: [Lagrange finite element of degree 1 on a triangle, Lagrange finite element of degree 1 on a triangle]") == 0)
    return new ffc_15_dof_map_0();
  if (strcmp(signature, "FFC dof map for Mixed finite element: [Lagrange finite element of degree 2 on a triangle, Lagrange finite element of degree 2 on a triangle]") == 0)
    return new ffc_16_dof_map_0();
  if (strcmp(signature, "FFC dof map for Mixed finite element: [Lagrange finite element of degree 1 on a tetrahedron, Lagrange finite element of degree 1 on a tetrahedron, Lagrange finite element of degree 1 on a tetrahedron]") == 0)
    return new ffc_17_dof_map_0();
  if (strcmp(signature, "FFC dof map for Mixed finite element: [Lagrange finite element of degree 2 on a tetrahedron, Lagrange finite element of degree 2 on a tetrahedron, Lagrange finite element of degree 2 on a tetrahedron]") == 0)
    return new ffc_18_dof_map_0();
  if (strcmp(signature, "FFC dof map for Mixed finite element: [Discontinuous Lagrange finite element of degree 0 on a triangle, Discontinuous Lagrange finite element of degree 0 on a triangle]") == 0)
    return new ffc_19_dof_map_0();
  if (strcmp(signature, "FFC dof map for Mixed finite element: [Discontinuous Lagrange finite element of degree 1 on a triangle, Discontinuous Lagrange finite element of degree 1 on a triangle]") == 0)
    return new ffc_20_dof_map_0();
  if (strcmp(signature, "FFC dof map for Mixed finite element: [Discontinuous Lagrange finite element of degree 0 on a tetrahedron, Discontinuous Lagrange finite element of degree 0 on a tetrahedron, Discontinuous Lagrange finite element of degree 0 on a tetrahedron]") == 0)
    return new ffc_21_dof_map_0();
  if (strcmp(signature, "FFC dof map for Mixed finite element: [Discontinuous Lagrange finite element of degree 1 on a tetrahedron, Discontinuous Lagrange finite element of degree 1 on a tetrahedron, Discontinuous Lagrange finite element of degree 1 on a tetrahedron]") == 0)
    return new ffc_22_dof_map_0();
  if (strcmp(signature, "FFC dof map for Mixed finite element: [Discontinuous Lagrange finite element of degree 2 on a tetrahedron, Discontinuous Lagrange finite element of degree 2 on a tetrahedron, Discontinuous Lagrange finite element of degree 2 on a tetrahedron]") == 0)
    return new ffc_23_dof_map_0();
  if (strcmp(signature, "FFC dof map for Brezzi-Douglas-Marini finite element of degree 1 on a triangle") == 0)
    return new ffc_24_dof_map_0();
  return 0;
}

#endif
