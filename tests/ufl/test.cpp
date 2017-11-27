#include "../tests.h"

#ifdef HAVE_CHECK

// Include all test cases
#include "Algebra/UFL_Algebra.h"
#include "Argument/UFL_Argument.h"
#include "array/UFL_array.h"
#include "CellSurfaceArea/UFL_CellSurfaceArea.h"
#include "Cell/UFL_Cell.h"
#include "CellVolume/UFL_CellVolume.h"
#include "Circumradius/UFL_Circumradius.h"
#include "Class/UFL_Class.h"
#include "Coefficient/UFL_Coefficient.h"
#include "Conditional/UFL_Conditional.h"
#include "Differentiation/UFL_Differentitation.h"
#include "Domain/UFL_Domain.h"
#include "ElementList/UFL_ElementList.h"
#include "EnrichedElement/UFL_EnrichedElement.h"
#include "Equation/UFL_Equation.h"
#include "Expression/UFL_Expression.h"
#include "FacetArea/UFL_FacetArea.h"
#include "FacetNormal/UFL_FacetNormal.h"
#include "Family/UFL_Family.h"
#include "FiniteElementSpace/UFL_FiniteElementSpace.h"
#include "FiniteElement/UFL_FiniteElement.h"
#include "Form/UFL_Form.h"
#include "GeometricQuantity/UFL_GeometricQuantitiy.h"
#include "Indexed/UFL_Indexed.h"
#include "IndexSum/UFL_IndexSum.h"
#include "Index/UFL_Index.h"
#include "Integral/UFL_Integral.h"
#include "List/UFL_List.h"
#include "MixedElement/UFL_MixedElement.h"
#include "Object/UFL_Object.h"
#include "QuadratureScheme/UFL_QuadratureScheme.h"
#include "repr/UFL_repr.h"
#include "RestrictedElement/UFL_RestrictedElement.h"
#include "Space/UFL_Space.h"
#include "SpatialCoordinate/UFL_SpatialCoordinate.h"
#include "TensorElement/UFL_TensorElement.h"
#include "Tensors/UFL_Tensors.h"
#include "tuple/UFL_tuple.h"
#include "Tuple/UFL_Tuple.h"
#include "type/UFL_type.h"
#include "Variable/UFL_Variable.h"
#include "VectorElement/UFL_VectorElement.h"

DOLFIN_SUITE_BEGIN(suite, "ufl")
{
  DOLFIN_TCASE_CREATE("UFL_Algebra");
  DOLFIN_TCASE_ADD(test_UFL_Algebra);

  DOLFIN_TCASE_CREATE("UFL_Argument");
  DOLFIN_TCASE_ADD(test_UFL_Argument);
 
  DOLFIN_TCASE_CREATE("UFL_array");
  DOLFIN_TCASE_ADD(test_UFL_array);

  DOLFIN_TCASE_CREATE("UFL_Cell");
  DOLFIN_TCASE_ADD(test_UFL_Cell);

  DOLFIN_TCASE_CREATE("UFL_Cell Surface area");
  DOLFIN_TCASE_ADD(test_UFL_CellSurfaceArea);

  DOLFIN_TCASE_CREATE("UFL_Cell Volume");
  DOLFIN_TCASE_ADD(test_UFL_CellVolume);

  DOLFIN_TCASE_CREATE("UFL_Circumradius");
  DOLFIN_TCASE_ADD(test_UFL_Circumradius);

  DOLFIN_TCASE_CREATE("UFL_Class");
  DOLFIN_TCASE_ADD(test_UFL_Class);

  DOLFIN_TCASE_CREATE("UFL_Coefficient");
  DOLFIN_TCASE_ADD(test_UFL_Coefficient);
  
  DOLFIN_TCASE_CREATE("UFL_Conditional");
  DOLFIN_TCASE_ADD(test_UFL_Conditional);

  DOLFIN_TCASE_CREATE("UFL_Differentiation");
  DOLFIN_TCASE_ADD(test_UFL_Differentiation);

  DOLFIN_TCASE_CREATE("UFL_Domain");
  DOLFIN_TCASE_ADD(test_UFL_Domain);

  DOLFIN_TCASE_CREATE("UFL_ElementList");
  DOLFIN_TCASE_ADD(test_UFL_ElementList);

  DOLFIN_TCASE_CREATE("UFL_EnrichedElement");
  DOLFIN_TCASE_ADD(test_UFL_EnrichedElement);

  DOLFIN_TCASE_CREATE("UFL_Equation");
  DOLFIN_TCASE_ADD(test_UFL_Equation);

  DOLFIN_TCASE_CREATE("UFL_Expression");
  DOLFIN_TCASE_ADD(test_UFL_Expression);

  DOLFIN_TCASE_CREATE("UFL_FacetArea");
  DOLFIN_TCASE_ADD(test_UFL_FacetArea);

  DOLFIN_TCASE_CREATE("UFL_FacetNormal");
  DOLFIN_TCASE_ADD(test_UFL_FacetNormal);

  DOLFIN_TCASE_CREATE("UFL_Family");
  DOLFIN_TCASE_ADD(test_UFL_Family);

  DOLFIN_TCASE_CREATE("UFL_FiniteElement");
  DOLFIN_TCASE_ADD(test_UFL_FiniteElement);

  DOLFIN_TCASE_CREATE("UFL_FiniteElementSpace");
  DOLFIN_TCASE_ADD(test_UFL_FiniteElementSpace);

  DOLFIN_TCASE_CREATE("UFL_Form");
  DOLFIN_TCASE_ADD(test_UFL_Form);

  DOLFIN_TCASE_CREATE("UFL_GeometricQuantity");
  DOLFIN_TCASE_ADD(test_UFL_GeometricQuantity);

  DOLFIN_TCASE_CREATE("UFL_Index");
  DOLFIN_TCASE_ADD(test_UFL_Index);

  DOLFIN_TCASE_CREATE("UFL_Indexed");
  DOLFIN_TCASE_ADD(test_UFL_Indexed);

  DOLFIN_TCASE_CREATE("UFL_IndexSum");
  DOLFIN_TCASE_ADD(test_UFL_IndexSum);

  DOLFIN_TCASE_CREATE("UFL_Integral");
  DOLFIN_TCASE_ADD(test_UFL_Integral);
  
  DOLFIN_TCASE_CREATE("UFL_List");
  DOLFIN_TCASE_ADD(test_UFL_List);

  DOLFIN_TCASE_CREATE("UFL_MixedElement");
  DOLFIN_TCASE_ADD(test_UFL_MixedElement);

  DOLFIN_TCASE_CREATE("UFL_Object");
  DOLFIN_TCASE_ADD(test_UFL_Object);
  
  DOLFIN_TCASE_CREATE("UFL_QuadratureScheme");
  DOLFIN_TCASE_ADD(test_UFL_QuadratureScheme);

  DOLFIN_TCASE_CREATE("UFL_repr");
  DOLFIN_TCASE_ADD(test_UFL_repr);

  DOLFIN_TCASE_CREATE("UFL_RestrictedElement");
  DOLFIN_TCASE_ADD(test_UFL_RestrictedElement);

  DOLFIN_TCASE_CREATE("UFL_Space");
  DOLFIN_TCASE_ADD(test_UFL_Space);

  DOLFIN_TCASE_CREATE("UFL_SpatialCoordinate");
  DOLFIN_TCASE_ADD(test_UFL_SpatialCoordinate);

  DOLFIN_TCASE_CREATE("UFL_TensorElement");
  DOLFIN_TCASE_ADD(test_UFL_TensorElement);

  DOLFIN_TCASE_CREATE("UFL_Tensors");
  DOLFIN_TCASE_ADD(test_UFL_Tensors);

  DOLFIN_TCASE_CREATE("UFL_tuple");
  DOLFIN_TCASE_ADD(test_UFL_tuple);

  DOLFIN_TCASE_CREATE("UFL_Tuple");
  DOLFIN_TCASE_ADD(test_UFL_Tuple);

  DOLFIN_TCASE_CREATE("UFL_type");
  DOLFIN_TCASE_ADD(test_UFL_type_int);
  DOLFIN_TCASE_ADD(test_UFL_type_real);
  DOLFIN_TCASE_ADD(test_UFL_type_string);
  
  DOLFIN_TCASE_CREATE("UFL_Variable");
  DOLFIN_TCASE_ADD(test_UFL_Variable);
  
  DOLFIN_TCASE_CREATE("UFL_VectorElement");
  DOLFIN_TCASE_ADD(test_UFL_VectorElement);
}
DOLFIN_SUITE_END
//-----------------------------------------------------------------------------
DOLFIN_CHECK_SUITE("dolfin/ufl", suite)
//-----------------------------------------------------------------------------

#endif
