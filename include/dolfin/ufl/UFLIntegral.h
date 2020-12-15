// Copyright (C) 2014 Bärbel Janssen.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_UFL_INTEGRAL_H
#define __DOLFIN_UFL_INTEGRAL_H

#include <dolfin/ufl/UFLClass.h>
#include <dolfin/ufl/UFLExpression.h>

namespace ufl
{

/**
 *  DOCUMENTATION:
 *
 *  @class  MeasureData
 *
 *  @brief  A domain description for different kinds of data that could be
 *  passed to a Measure.
 */

class MeasureData : public type<std::string>
{

public:

  enum Type
  {
    None
  };

  ///
  MeasureData(Type const& t);

  ///
  MeasureData(repr_t const& repr);

  ///
  ~MeasureData() override;

  ///
  Type type() const;

  ///
  void display() const override;

private:

  MeasureData::Type const type_;

  //--- STATIC ----------------------------------------------------------------

  ///
  static std::string type_repr(MeasureData::Type const& t);

  ///
  static MeasureData::Type repr_type(repr_t const& repr);
  typedef dolfin::_ordered_map<repr_t, MeasureData::Type> MappingReprToType;
  typedef dolfin::_ordered_map<MeasureData::Type, repr_t> MappingTypeToRepr;
  typedef std::pair<repr_t, MeasureData::Type> MappingReprToTypeItem;
  typedef std::pair<MeasureData::Type, repr_t> MappingTypeToReprItem;

  static MappingReprToType const& MappingToType()
  {
    static MappingReprToType const MappingToType =
        __init_mapping_repr_to_type();
    return MappingToType;
  }

  static MappingTypeToRepr const& MappingToRepr()
  {
    static MappingTypeToRepr const MappingToRepr =
        __init_mapping_type_to_repr();
    return MappingToRepr;
  }

  static MappingReprToType __init_mapping_repr_to_type();
  static MappingTypeToRepr __init_mapping_type_to_repr();
};

/**
 *  DOCUMENTATION:
 *
 *  @class  MeasureDomain
 *
 *  @brief  A domain description for different kinds of measures.
 */

class MeasureDomain : public type<std::string>
{

public:

  enum Type
  {
    None, cell, exterior_facet, interior_facet, macro_cell, surface
  };

  ///
  MeasureDomain(Type const& t);

  ///
  MeasureDomain(repr_t const& repr);

  ///
  ~MeasureDomain() override;

  ///
  Type type() const;

  ///
  void display() const override;

private:

  MeasureDomain::Type const type_;

  //--- STATIC ----------------------------------------------------------------

public:
  ///
  static std::string type_repr(MeasureDomain::Type const& t);

private:
  ///
  static MeasureDomain::Type repr_type(repr_t const& repr);
  typedef dolfin::_ordered_map<repr_t, MeasureDomain::Type> MappingReprToType;
  typedef dolfin::_ordered_map<MeasureDomain::Type, repr_t> MappingTypeToRepr;
  typedef std::pair<repr_t, MeasureDomain::Type> MappingReprToTypeItem;
  typedef std::pair<MeasureDomain::Type, repr_t> MappingTypeToReprItem;

  static MappingReprToType MappingToType()
  {
    static MappingReprToType const MappingToType =
        __init_mapping_repr_to_type();
    return MappingToType;
  }

  static MappingTypeToRepr MappingToRepr()
  {
    static MappingTypeToRepr const MappingToRepr =
        __init_mapping_type_to_repr();
    return MappingToRepr;
  }

  static MappingReprToType __init_mapping_repr_to_type();
  static MappingTypeToRepr __init_mapping_type_to_repr();
};

/**
 *  DOCUMENTATION:
 *
 *  @class  Measure
 *
 *  @brief  Provides an interface complying with UFL Measure.
 */

class Expression;

class Measure : public Class
{
public:

  ///
  Measure(MeasureDomain::Type const& measure_type,
          MeasureData const& measure_data, dolfin::uint& measure_id);

  ///
  Measure(repr_t const& repr);

//      Measure(MeasureDomain::Type const& measure_type,
//          MeasureData const& meta_data,
//          MeasureData const& measure_data,
//          dolfin::uint measure_id);

///
  ~Measure() override;

  ///
  std::vector<Class const *> const operands(std::string const& name) const;

  ///
  std::vector<std::vector<Class const *> > const level_operands(
      std::vector<std::vector<Class const *> > const& operands) const;

  ///Evaluate the expression tree at the given quadrature_points
  virtual std::vector<std::vector<std::vector<dolfin::real> > > const evaluate(
      dolfin::uint n,
      std::vector<std::vector<std::vector<dolfin::real> > > const& tensor,
      dolfin::UFCReferenceCell const& ref_cell,
      std::vector<dolfin::real*> const& q_points,
      const double * const * coordinates) const;

  //--- INTERFACE -------------------------------------------------------------

  static Measure const * create(Object::repr_t const& repr);

  MeasureDomain::Type measure_type() const;

//      MeasureData const& measure_data() const;

  MeasureData const& meta_data() const;

  dolfin::uint measure_domain_id() const;

  //--- INTERFACE inherited from UFLClass -------------------------------------

  /// __repr__
  repr_t const& repr() const override;

  /// __str__
  std::string const& str() const override;

  ///
  void display() const override;

private:
  MeasureDomain const measure_domain_;
  type<dolfin::uint> const measure_id_;
  MeasureData const meta_data_;

  repr_t const repr_;
  std::string const str_;
};

class Integral : public Class
{
public:

  ///
  Integral(Expression const& integrand, const Measure& measure);

  ///
  Integral(repr_t const & repr);

  ///
  ~Integral() override;

  ///
  std::vector<Class const *> const operands(std::string const& name) const;

  ///
  std::vector<std::vector<Class const *> > const level_operands(
      std::vector<std::vector<Class const *> > const& operands) const;

  ///Evaluate the expression tree at the given quadrature_points
  virtual std::vector<std::vector<std::vector<dolfin::real> > > const evaluate(
      dolfin::uint n,
      std::vector<std::vector<std::vector<dolfin::real> > > const& tensor,
      dolfin::UFCReferenceCell const& ref_cell,
      std::vector<dolfin::real*> const& q_points,
      const double * const * coordinates) const;

  //--- INTERFACE -------------------------------------------------------------

  static Integral const * create(Object::repr_t const& repr);

  std::vector<Expression const *> const& integrand() const;

  Measure const& measure() const;

  //--- INTERFACE inherited from UFLClass -------------------------------------

  /// __repr__
  repr_t const& repr() const override;

  /// __str__
  std::string const& str() const override;

  ///
  void display() const override;

private:

  std::vector<Expression const *> fill_expressions(
      std::vector<repr_t> const& reprs);
  std::vector<Expression const *> fill_expressions(Expression const& e);

  std::vector<Expression const *> const expressions_;
  Measure const measure_;

  repr_t const repr_;
  std::string const str_;
};
} /* namespace ufl */
#endif /* __DOLFIN_UFL_INTEGRAL_H */
