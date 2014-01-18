#include <dolfin/elements/FE.h>

namespace FE
{

//-----------------------------------------------------------------------------
Array<std::string> const init_FunctionSpaceNames()
{
  Array<std::string> ret;
  ret.push_back(FE::FINITE_ELEMENT);
  ret.push_back(FE::VECTOR_ELEMENT);
  ret.push_back(FE::MIXED_ELEMENT);
  return ret;
}

Array<std::string> const FE::FunctionSpace::Names = init_FunctionSpaceNames();

//-----------------------------------------------------------------------------
std::map<FE::FunctionSpace::Type, std::string> const init_TypesOfFunctionSpaces()
{
  std::map<FE::FunctionSpace::Type, std::string> m;
  m.insert(std::pair<FE::FunctionSpace::Type, std::string>(FE::FunctionSpace::FiniteElement, FE::FINITE_ELEMENT));
  m.insert(std::pair<FE::FunctionSpace::Type, std::string>(FE::FunctionSpace::VectorElement, FE::VECTOR_ELEMENT));
  m.insert(std::pair<FE::FunctionSpace::Type, std::string>(FE::FunctionSpace::MixedElement, FE::MIXED_ELEMENT));
  return m;
}

std::map<FE::FunctionSpace::Type, std::string> const FE::FunctionSpace::TypesOfFunctionSpaces =  init_TypesOfFunctionSpaces();

//-----------------------------------------------------------------------------
std::map<std::string, FE::FunctionSpace::Type> const init_NamesOfFunctionSpaces()
{
  std::map<std::string, FE::FunctionSpace::Type> m;
  m.insert(std::pair<std::string, FE::FunctionSpace::Type>(FE::FINITE_ELEMENT, FE::FunctionSpace::FiniteElement));
  m.insert(std::pair<std::string, FE::FunctionSpace::Type>(FE::VECTOR_ELEMENT, FE::FunctionSpace::VectorElement));
  m.insert(std::pair<std::string, FE::FunctionSpace::Type>(FE::MIXED_ELEMENT, FE::FunctionSpace::MixedElement));
  return m;
}

std::map<std::string, FE::FunctionSpace::Type> const FE::FunctionSpace::NamesOfFunctionSpaces =  init_NamesOfFunctionSpaces();

//-----------------------------------------------------------------------------
Array<std::string> const init_FamilyNames()
{
  Array<std::string> ret;
  ret.push_back(FE::LAGRANGE);
  ret.push_back(FE::DG);
  ret.push_back(FE::BDM);
  return ret;
}

Array<std::string> const FE::Family::Names = init_FamilyNames();

//-----------------------------------------------------------------------------
std::map<FE::Family::Type, std::string> const init_TypesOfFamilies()
{
  std::map<FE::Family::Type, std::string> m;
  m.insert(std::pair<FE::Family::Type, std::string>(FE::Family::Lagrange, FE::LAGRANGE));
  m.insert(std::pair<FE::Family::Type, std::string>(FE::Family::DG, FE::DG));
  m.insert(std::pair<FE::Family::Type, std::string>(FE::Family::BDM, FE::BDM));
  return m;
}

std::map<FE::Family::Type, std::string> const FE::Family::TypesOfFamilies =  init_TypesOfFamilies();

//-----------------------------------------------------------------------------
std::map<std::string, FE::Family::Type> const init_NamesOfFamilies()
{
  std::map<std::string, FE::Family::Type> m;
  m.insert(std::pair<std::string, FE::Family::Type>(FE::LAGRANGE, FE::Family::Lagrange));
  m.insert(std::pair<std::string, FE::Family::Type>(FE::DG, FE::Family::DG));
  m.insert(std::pair<std::string, FE::Family::Type>(FE::BDM, FE::Family::BDM));
  return m;
}

std::map<std::string, FE::Family::Type> const FE::Family::NamesOfFamilies =  init_NamesOfFamilies();

//-----------------------------------------------------------------------------
Array<std::string> const init_CellNames()
{
  Array<std::string> ret;
  ret.push_back(FE::INTERVAL);
  ret.push_back(FE::TRIANGLE);
  ret.push_back(FE::QUADRILATERAL);
  ret.push_back(FE::TETRAHEDRON);
  ret.push_back(FE::HEXAHEDRON);
  return ret;
}

Array<std::string> const FE::Cell::Names = init_CellNames();

//-----------------------------------------------------------------------------
std::map<FE::Cell::Type, std::string> const init_TypesOfCells()
{
  std::map<FE::Cell::Type, std::string> m;
  m.insert(std::pair<FE::Cell::Type, std::string>(FE::Cell::interval, FE::INTERVAL));
  m.insert(std::pair<FE::Cell::Type, std::string>(FE::Cell::triangle, FE::TRIANGLE));
  m.insert(std::pair<FE::Cell::Type, std::string>(FE::Cell::quadrilateral, FE::QUADRILATERAL));
  m.insert(std::pair<FE::Cell::Type, std::string>(FE::Cell::tetrahedron, FE::TETRAHEDRON));
  m.insert(std::pair<FE::Cell::Type, std::string>(FE::Cell::hexahedron, FE::HEXAHEDRON));
  return m;
}

std::map<FE::Cell::Type, std::string> const FE::Cell::TypesOfCells =  init_TypesOfCells();

//-----------------------------------------------------------------------------
std::map<std::string, FE::Cell::Type> const init_NamesOfCells()
{
  std::map<std::string, FE::Cell::Type> m;
  m.insert(std::pair<std::string, FE::Cell::Type>(FE::INTERVAL, FE::Cell::interval));
  m.insert(std::pair<std::string, FE::Cell::Type>(FE::TRIANGLE, FE::Cell::triangle));
  m.insert(std::pair<std::string, FE::Cell::Type>(FE::QUADRILATERAL, FE::Cell::quadrilateral));
  m.insert(std::pair<std::string, FE::Cell::Type>(FE::TETRAHEDRON, FE::Cell::tetrahedron));
  m.insert(std::pair<std::string, FE::Cell::Type>(FE::HEXAHEDRON, FE::Cell::hexahedron));
  return m;
}

std::map<std::string, FE::Cell::Type> const FE::Cell::NamesOfCells =  init_NamesOfCells();

//-----------------------------------------------------------------------------
std::string const get_signature(FE::FunctionSpace::Type const space,
                            FE::Family::Type const family,
                            FE::Cell::Type const cell,
                            uint const space_dim,
                            uint const degree,
                            uint const value_dim)
{
  std::stringstream ss;

  // FiniteElement('Lagrange', Cell('interval', Space(1)), 1, None)

  ss << FE::FunctionSpace::type2string(space) << "(";
  ss << "'" << FE::Family::type2string(family) << "', ";
  ss << "Cell('" << FE::Cell::type2string(cell) << "', ";
  ss << "Space(" << space_dim << ")), ";
  ss << degree << ", ";
  if (space == FE::FunctionSpace::VectorElement)
  {
    ss << value_dim << ", ";
  }
  ss << ", None)";

  return ss.str();
}

//-----------------------------------------------------------------------------
attributes const get_attributes(std::string const signature)
{
  std::string type;
  std::string family;
  std::string shape;
  uint space = 1;
  uint degree = -1;
  uint value = 1;

  //
  std::string s(signature);
  size_t pos = s.find("(", 0);
  type = s.substr(0, pos);

  //FIXME: Not mixed element aware
  if (type != FE::MIXED_ELEMENT)
  {
    s.erase(0, pos + 1);
    size_t element = 0;
    size_t t0 = 0;
    size_t t1 = 0;
    std::string tok;
    while ((t1 = s.find(",", t0)) != std::string::npos)
    {
      tok = s.substr(t0, t1 - t0);
      if (element == 0)
      {
        family = s.substr(t0 + 1, t1 - t0 - 2);
      }
      if (element == 1)
      {
        shape = s.substr(t0 + 7, t1 - t0 - 8);
      }
      if (element == 2)
      {
        std::stringstream ss;
        ss << s.substr(t0 + 7, t1 - t0 - 9);
        ss >> space;
      }
      if (element == 3)
      {
        std::stringstream ss;
        ss << s.substr(t0 + 1, t1 - t0 - 1);
        ss >> degree;
      }
      if (element == 4)
      {
        std::stringstream ss;
        ss << s.substr(t0 + 1, t1 - t0 - 1);
        ss >> value;
      }
      ++element;
      t0 = t1 + 1;
    }
  }
  else
  {
    family = FE::MIXED_ELEMENT;
  }
  return FE::attributes(type, family, shape, space, degree,
                        value);
}

}
