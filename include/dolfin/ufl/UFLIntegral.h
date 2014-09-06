// Copyright (C) 2014 Bärbel Janssen.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:
// Last changed:

#ifndef __UFL_INTEGRAL_H_
#define __UFL_INTEGRAL_H_

#include <dolfin/ufl/UFLClass.h>
#include <dolfin/ufl/UFLExpression.h>

namespace ufl
{

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
        None,
        cell,
        exterior_facet,
        interior_facet,
        macro_cell,
        surface
      };

      ///
      MeasureDomain(Type const& t);
      
      ///
      MeasureDomain(repr_t const& repr);

      ///
      ~MeasureDomain();

      ///
      Type const type() const;

      ///
      void display() const;

    private:

      MeasureDomain::Type const type_;

      //--- STATIC ----------------------------------------------------------------

      ///
      static std::string const type_repr(MeasureDomain::Type const& t);

      ///
      static MeasureDomain::Type const repr_type(repr_t const& repr);
      typedef std::map<repr_t const, MeasureDomain::Type> MappingReprToType;
      typedef std::map<MeasureDomain::Type, repr_t const> MappingTypeToRepr;
      typedef std::pair<repr_t const, MeasureDomain::Type> MappingReprToTypeItem;
      typedef std::pair<MeasureDomain::Type, repr_t const> MappingTypeToReprItem;

      static MappingReprToType const MappingToType()
      {
        static MappingReprToType const MappingToType = __init_mapping_repr_to_type();
        return MappingToType;
      }

      static MappingTypeToRepr const MappingToRepr()
      {
        static MappingTypeToRepr const MappingToRepr = __init_mapping_type_to_repr();
        return MappingToRepr;
      }

      static MappingReprToType const __init_mapping_repr_to_type();
      static MappingTypeToRepr const __init_mapping_type_to_repr();
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
//          MeasureData const& measure_data, 
          dolfin::uint measure_id);

      ///
      Measure(repr_t const & repr);
                                                                                                    
//      Measure(MeasureDomain::Type const& measure_type, 
//          MeasureData const& meta_data, 
//          MeasureData const& measure_data, 
//          dolfin::uint measure_id);

      ///
      ~Measure();

      //--- INTERFACE -------------------------------------------------------------

//      Class const& reconstruct(const Class& measure);

      MeasureDomain::Type const& measure_type() const;

//      MeasureData const& measure_data() const;

//      MeasureData const& meta_data() const;

      dolfin::uint measure_id() const;

      //--- INTERFACE inherited from UFLClass -------------------------------------

      /// __repr__
      repr_t const repr() const;

      /// __str__
      std::string const str() const;

      ///
      void display() const;

    private:
      MeasureDomain::Type const measure_;
      type<dolfin::uint> const measure_id_;

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
      ~Integral();

      //--- INTERFACE -------------------------------------------------------------

      Expression const& integrand() const;

      Measure const& measure() const;

//      Class const& reconstruct(const Class& integral);

      //--- INTERFACE inherited from UFLClass -------------------------------------

      /// __repr__
      repr_t const repr() const;

      /// __str__
      std::string const str() const;

      ///
      void display() const;

    private:

      Expression const& integrand_;
      Measure const measure_;

      repr_t const repr_;
      std::string const str_;
  };
} /* namespace ufl */
#endif /* __UFL_INTEGRAL_H_ */
