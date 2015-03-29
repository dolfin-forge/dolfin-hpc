#include <dolfin/config/dolfin_config.h>
#include <dolfin/main/init.h>

#include <iostream>
#include <iomanip>

#ifdef HAVE_CHECK

#include <check.h>

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
#include "FiniteElementBase/UFL_FiniteElementBase.h"
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

int argc;
char **argv;

void setup()
{
  //  dolfin_init(argc, argv);
}

void teardown()
{
  //  dolfin_finalize();
}

Suite *ufl_suite()
{
  TCase *tc;
  Suite *s;

  s = suite_create("UFL");

  tcase_add_test(tc, test_UFL_Algebra); 
  tc = tcase_create("UFL Algebra"); 
  suite_add_tcase (s, tc);
  tcase_add_checked_fixture (tc, setup, teardown);

  tc = tcase_create("UFL Argument"); 
  tcase_add_test(tc, test_UFL_Argument);
  suite_add_tcase (s, tc);
  tcase_add_checked_fixture (tc, setup, teardown);
 
  tc = tcase_create("UFL array"); 
  tcase_add_test(tc, test_UFL_array);
  suite_add_tcase (s, tc);
  tcase_add_checked_fixture (tc, setup, teardown);

  tc = tcase_create("UFL Cell"); 
  tcase_add_test(tc, test_UFL_Cell);
  suite_add_tcase (s, tc);
  tcase_add_checked_fixture (tc, setup, teardown);

  tc = tcase_create("UFL Cell Surface area"); 
  tcase_add_test(tc, test_UFL_CellSurfaceArea);
  suite_add_tcase (s, tc);
  tcase_add_checked_fixture (tc, setup, teardown);

  tc = tcase_create("UFL Cell Volume"); 
  tcase_add_test(tc, test_UFL_CellVolume);
  suite_add_tcase (s, tc);
  tcase_add_checked_fixture (tc, setup, teardown);;

  tc = tcase_create("UFL Circumradius"); 
  tcase_add_test(tc, test_UFL_Circumradius);
  suite_add_tcase (s, tc);
  tcase_add_checked_fixture (tc, setup, teardown);;

  tc = tcase_create("UFL Class");
  tcase_add_test(tc, test_UFL_Class);
  suite_add_tcase (s, tc);
  tcase_add_checked_fixture (tc, setup, teardown);;

  tc = tcase_create("UFL Coefficient");
  tcase_add_test(tc, test_UFL_Coefficient);
  suite_add_tcase (s, tc);
  tcase_add_checked_fixture (tc, setup, teardown);;
  
  tc = tcase_create("UFL Conditional");
  tcase_add_test(tc, test_UFL_Conditional);
  suite_add_tcase (s, tc);
  tcase_add_checked_fixture (tc, setup, teardown);;

  tc = tcase_create("UFL Differentiation");
  tcase_add_test(tc, test_UFL_Differentiation);
  suite_add_tcase (s, tc);
  tcase_add_checked_fixture (tc, setup, teardown);;

  tc = tcase_create("UFL Domain");
  tcase_add_test(tc, test_UFL_Domain);
  suite_add_tcase (s, tc);
  tcase_add_checked_fixture (tc, setup, teardown);;

  tc = tcase_create("UFL ElementList");
  tcase_add_test(tc, test_UFL_ElementList);
  suite_add_tcase (s, tc);
  tcase_add_checked_fixture (tc, setup, teardown);;

  tc = tcase_create("UFL EnrichedElement");
  tcase_add_test(tc, test_UFL_EnrichedElement); 
  suite_add_tcase (s, tc);
  tcase_add_checked_fixture (tc, setup, teardown);;

  tc = tcase_create("UFL Equation");
  tcase_add_test(tc, test_UFL_Equation); 
  suite_add_tcase (s, tc);
  tcase_add_checked_fixture (tc, setup, teardown);;

  tc = tcase_create("UFL Expression");
  tcase_add_test(tc, test_UFL_Expression); 
  suite_add_tcase (s, tc);
  tcase_add_checked_fixture (tc, setup, teardown);;

  tc = tcase_create("UFL FacetArea");
  tcase_add_test(tc, test_UFL_FacetArea); 
  suite_add_tcase (s, tc);
  tcase_add_checked_fixture (tc, setup, teardown);;

  tc = tcase_create("UFL FacetNormal");
  tcase_add_test(tc, test_UFL_FacetNormal); 
  suite_add_tcase (s, tc);
  tcase_add_checked_fixture (tc, setup, teardown);;

  tc = tcase_create("UFL Family");
  tcase_add_test(tc, test_UFL_Family); 
  suite_add_tcase (s, tc);
  tcase_add_checked_fixture (tc, setup, teardown);;

  tc = tcase_create("UFL FiniteElement");
  tcase_add_test(tc, test_UFL_FiniteElement);
  suite_add_tcase (s, tc);
  tcase_add_checked_fixture (tc, setup, teardown);;

  tc = tcase_create("UFL FiniteElementBase");
  tcase_add_test(tc, test_UFL_FiniteElementBase);
  suite_add_tcase (s, tc);
  tcase_add_checked_fixture (tc, setup, teardown);;

  tc = tcase_create("UFL Form");
  tcase_add_test(tc, test_UFL_Form);
  suite_add_tcase (s, tc);
  tcase_add_checked_fixture (tc, setup, teardown);;

  tc = tcase_create("UFL GeometricQuantity");
  tcase_add_test(tc, test_UFL_GeometricQuantity);
  suite_add_tcase (s, tc);
  tcase_add_checked_fixture (tc, setup, teardown);;

  tc = tcase_create("UFL Index");
  tcase_add_test(tc, test_UFL_Index);
  suite_add_tcase (s, tc);
  tcase_add_checked_fixture (tc, setup, teardown);;

  tc = tcase_create("UFL Indexed");
  tcase_add_test(tc, test_UFL_Indexed);
  suite_add_tcase (s, tc);
  tcase_add_checked_fixture (tc, setup, teardown);;

  tc = tcase_create("UFL IndexSum");
  tcase_add_test(tc, test_UFL_IndexSum);
  suite_add_tcase (s, tc);
  tcase_add_checked_fixture (tc, setup, teardown);;

  tc = tcase_create("UFL Integral");
  tcase_add_test(tc, test_UFL_Integral);
  suite_add_tcase (s, tc);
  tcase_add_checked_fixture (tc, setup, teardown);;
  
  tc = tcase_create("UFL List");
  tcase_add_test(tc, test_UFL_List);
  suite_add_tcase (s, tc);
  tcase_add_checked_fixture (tc, setup, teardown);;

  tc =  tcase_create("UFL MixedElement");
  tcase_add_test(tc, test_UFL_MixedElement);
  suite_add_tcase (s, tc);
  tcase_add_checked_fixture (tc, setup, teardown);;

  tc = tcase_create("UFL Object");
  tcase_add_test(tc, test_UFL_Object);
  suite_add_tcase (s, tc);
  tcase_add_checked_fixture (tc, setup, teardown);;
  
  tc = tcase_create("UFL QuadratureScheme");
  tcase_add_test(tc, test_UFL_QuadratureScheme);
  suite_add_tcase (s, tc);
  tcase_add_checked_fixture (tc, setup, teardown);;

  tc = tcase_create("UFL repr");
  tcase_add_test(tc, test_UFL_repr);
  suite_add_tcase (s, tc);
  tcase_add_checked_fixture (tc, setup, teardown);;

  tc = tcase_create("UFL RestrictedElement");
  tcase_add_test(tc, test_UFL_RestrictedElement);
  suite_add_tcase (s, tc);
  tcase_add_checked_fixture (tc, setup, teardown);;

  tc = tcase_create("UFL Space");
  tcase_add_test(tc, test_UFL_Space);
  suite_add_tcase (s, tc);
  tcase_add_checked_fixture (tc, setup, teardown);;

  tc = tcase_create("UFL SpatialCoordinate");
  tcase_add_test(tc, test_UFL_SpatialCoordinate);
  suite_add_tcase (s, tc);
  tcase_add_checked_fixture (tc, setup, teardown);;

  tc = tcase_create("UFL TensorElement");
  tcase_add_test(tc, test_UFL_TensorElement);
  suite_add_tcase (s, tc);
  tcase_add_checked_fixture (tc, setup, teardown);;

  tc = tcase_create("UFL Tensors");
  tcase_add_test(tc, test_UFL_Tensors);
  suite_add_tcase (s, tc);
  tcase_add_checked_fixture (tc, setup, teardown);;

  tc = tcase_create("UFL tuple");
  tcase_add_test(tc, test_UFL_tuple);
  suite_add_tcase (s, tc);
  tcase_add_checked_fixture (tc, setup, teardown);;

  tc = tcase_create("UFL Tuple");
  tcase_add_test(tc, test_UFL_Tuple);
  suite_add_tcase (s, tc);
  tcase_add_checked_fixture (tc, setup, teardown);;

  tc = tcase_create("UFL type");
  tcase_add_test(tc, test_UFL_type_int);
  tcase_add_test(tc, test_UFL_type_real);
  tcase_add_test(tc, test_UFL_type_string);
  suite_add_tcase (s, tc);
  tcase_add_checked_fixture (tc, setup, teardown);;
  
  tc = tcase_create("UFL Variable");
  tcase_add_test(tc, test_UFL_Variable);
  suite_add_tcase (s, tc);
  tcase_add_checked_fixture (tc, setup, teardown);;
  
  tc = tcase_create("UFL VectorElement");
  tcase_add_test(tc, test_UFL_VectorElement);
  suite_add_tcase (s, tc);
  tcase_add_checked_fixture (tc, setup, teardown);;

  return s;
}

int main(void)
{
  int number_failed;
  Suite* s = ufl_suite();
  SRunner* sr = srunner_create(s);

  srunner_run_all(sr, CK_NORMAL);
  number_failed = srunner_ntests_failed(sr);
  srunner_free(sr);

  return (number_failed == 0) ? 0 : 1;
}

#else

int main(void)
{
  fprintf(stderr, "*** Check is required for dolfin/ufl tests ***\n");
  return 0;
}

#endif
