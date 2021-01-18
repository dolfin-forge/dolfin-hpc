// Copyright (C) 2017 Aurelien Larcher.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_ELEMENTS_H
#define __DOLFIN_ELEMENTS_H

#include <string>

namespace dolfin
{

// we only really use CG and DG at the moment, but here is a complete list of
// element families ufc/ffc can generate
enum ElementFamilyType
{
  ARG,
  AW,
  BDFM,
  BDM,
  CR,
  DG,
  HER,
  CG,
  MTW,
  MOR,
  N1curl,
  N2curl,
  RT,
  BQ,
  B,
  Q,
  R,
  U,
  Mixed,
  Vector,
  Tensor,
  Enriched,
  Restricted
};

struct Elements
{
  struct cg
  {
    static ElementFamilyType const type = ElementFamilyType::CG;
  };

  struct dg
  {
    static ElementFamilyType const type = ElementFamilyType::DG;
  };
};

inline bool operator==( std::string const &       family_str,
                        ElementFamilyType const & family_type )
{
  bool ret = false;
  switch ( family_type )
  {
    case CG:
      ret = family_str == "Lagrange";
      break;
    case DG:
      ret = family_str == "Discontinuous Lagrange";
      break;
    case ARG:
      ret = family_str == "Argyris";
      break;
    case AW:
      ret = family_str == "Arnold-Winther";
      break;
    case BDFM:
      ret = family_str == "Brezzi-Douglas-Fortin-Marini";
      break;
    case BDM:
      ret = family_str == "Brezzi-Douglas-Marini";
      break;
    case CR:
      ret = family_str == "Crouzeix-Raviart";
      break;
    case HER:
      ret = family_str == "Hermite";
      break;
    case MTW:
      ret = family_str == "Mardal-Tai-Winther";
      break;
    case MOR:
      ret = family_str == "Morley";
      break;
    case N1curl:
      ret = family_str == "Nedelec 1st kind H(curl)";
      break;
    case N2curl:
      ret = family_str == "Nedelec 2nd kind H(curl)";
      break;
    case RT:
      ret = family_str == "Raviart-Thomas";
      break;
    case BQ:
      ret = family_str == "Boundary Quadrature";
      break;
    case B:
      ret = family_str == "Bubble";
      break;
    case R:
      ret = family_str == "Real";
      break;
    case U:
      ret = family_str == "Undefined";
    // case Q:
    //   ret = family_str == ;
    //   break;
    // case Mixed:
    //   ret = family_str == ;
    //   break;
    // case Vector:
    //   ret = family_str == ;
    //   break;
    // case Tensor:
    //   ret = family_str == ;
    //   break;
    // case Enriched:
    //   ret = family_str == ;
    //   break;
    // case Restricted:
    //   ret = family_str == ;
    //   break;
    default:
      ret = false;
      break;
  }

  return ret;
}

inline bool operator==( ElementFamilyType const & family_type,
                        std::string const &       family_str )
{
  return family_str == family_type;
}

inline bool operator==( char const *              family_str,
                        ElementFamilyType const & family_type )
{
  return std::string( family_str ) == family_type;
}

inline bool operator==( ElementFamilyType const & family_type,
                        char const *              family_str )
{
  return std::string( family_str ) == family_type;
}

} /* namespace dolfin */

#endif /* __DOLFIN_ELEMENTS_H */
