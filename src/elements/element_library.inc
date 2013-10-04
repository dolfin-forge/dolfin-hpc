// Automatically generated code mapping element and dof map signatures
// to the corresponding ufc::finite_element and ufc::dof_map classes

#include <dolfin/config/dolfin_config.h>
#include <cstring>

#include <dolfin/fem/DofMap.h>

#if ENABLE_UFL
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

#endif

