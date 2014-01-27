// Copyright (C) 2014 Bärbel Janssen.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:
// Last changed:

#ifndef __UFL_INTEGRAL_H_
#define __UFL_INTEGRAL_H_


namespace ufl
{

/**
 *  DOCUMENTATION:
 *
 *  @class  Integral
 *
 *  @brief  Provides an interface complying with UFL Integral.
 */

  class Integral : public Expression
  {
    public:

      ///
      Integral(const Expression& integrand, const Measure& measure);

      ///
      ~Integral();

      //--- INTERFACE -------------------------------------------------------------

      const Expression& integrand();

      const Measure& measure();

      const Expression& reconstruct(const Expression& integral)
      //--- INTERFACE inherited from UFLClass -------------------------------------

      /// __repr__
      std::string const repr() const;

      /// __str__
      std::string const str() const;

    private:
      const Expression integrand;
      const Measure measure;

      std::string const repr_;
      std::string const str_;
  };

  class Measure : public Expression
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
      Measure(const Type& measure_type, 
          const MeasureData& measure_data, 
          unsigned int measure_id);
                                                                                                    
      Measure(const Type& measure_type, 
          const MeasureData& meta_data, 
          const MeasureData& measure_data, 
          unsigned int measure_id);

      ///
      ~Measure();

      //--- INTERFACE -------------------------------------------------------------

      const Expression& reconstruct(const Expression& measure)

      const Type& measure_type();

      const MeasureData& measure_data();

      const MeasureData& meta_data();

      const unsigned int& measure_id();

      //--- INTERFACE inherited from UFLClass -------------------------------------

      /// __repr__
      std::string const repr() const;

      /// __str__
      std::string const str() const;

    private:
      std::string const repr_;
      std::string const str_;
  };
} /* namespace ufl */
#endif /* __UFL_INTEGRAL_H_ */
