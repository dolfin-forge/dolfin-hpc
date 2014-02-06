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
 *  @class  Integral
 *
 *  @brief  Provides an interface complying with UFL Integral.
 */

  class Expression;

  class Measure : public Class
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
      Measure(Type const& measure_type, 
//          MeasureData const& measure_data, 
          dolfin::uint measure_id);
                                                                                                    
//      Measure(Type const& measure_type, 
//          MeasureData const& meta_data, 
//          MeasureData const& measure_data, 
//          dolfin::uint measure_id);

      ///
      ~Measure();

      //--- INTERFACE -------------------------------------------------------------

//      Class const& reconstruct(const Class& measure);

      Type const& measure_type() const;

//      MeasureData const& measure_data() const;

//      MeasureData const& meta_data() const;

      dolfin::uint const& measure_id() const;

      //--- INTERFACE inherited from UFLClass -------------------------------------

      /// __repr__
      repr_t const repr() const;

      /// __str__
      std::string const str() const;

      ///
      void display() const;

      ///
      Measure const* create(repr_t const & repr) const;

    private:
      Measure::Type const measure_;
      dolfin::uint const measure_id_;
      mutable repr_t repr_;
      mutable std::string str_;
  };

  class Integral : public Class
  {
    public:

      ///
      Integral(Expression const& integrand, const Measure& measure);

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

      ///
      Integral const* create(repr_t const & repr) const;

    private:

      Expression const& integrand_;
      Measure const measure_;

      mutable repr_t repr_;
      mutable std::string str_;
  };
} /* namespace ufl */
#endif /* __UFL_INTEGRAL_H_ */
