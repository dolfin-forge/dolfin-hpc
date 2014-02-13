// Copyright (C) 2014 Bärbel Janssen.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:
// Last changed:

#ifndef __UFL_FORM_H_
#define __UFL_FORM_H_

#include <dolfin/ufl/UFLList.h>

namespace ufl
{

/**
 *  DOCUMENTATION:
 *
 *  @class  Form
 *
 *  @brief  Provides an interface complying with UFL Form.
 */

  class Form : public Class
  {
    public:

      ///
      Form(List const& list);

      ///
      Form(repr_t const & repr);

      ///
      ~Form();

      //--- INTERFACE -------------------------------------------------------------

      Cell const& cell() const;

      std::vector<Integral> const& integrals(Measure::Type const& measure_type=Measure::None) const;

      std::vector<Measure> const measures(Measure::Type const& measure_type=Measure::None) const;

      std::vector<Measure::Type> const domains(Measure::Type const& measure_type=Measure::None) const;

      std::vector<Integral> const& cell_integrals() const;

      std::vector<Integral> const& exterior_facet_integrals() const;
      
      std::vector<Integral> const& interior_facet_integrals() const;
      
      std::vector<Integral> const& macro_cell_integrals() const;
      
      std::vector<Integral> const& surface_integrals() const;
      
//      FormData const& form_data() const;
      
//      FormData const& compute_form_data() const;

      bool const is_preprocessed() const;
      
      //--- INTERFACE inherited from UFLClass -------------------------------------

      repr_t const repr() const;

      /// __str__
      std::string const str() const;
      
      ///
      void display() const;

    private:
      List const list_;

//      const FormData form_data;

      repr_t const repr_;
      std::string const str_;

      bool const is_preprocessed_;
  };
  
  

} /* namespace ufl */
#endif /* __UFL_FORM_H_ */
