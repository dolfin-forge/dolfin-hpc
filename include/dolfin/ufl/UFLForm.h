// Copyright (C) 2014 Bärbel Janssen.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:
// Last changed:

#ifndef __UFL_FORM_H_
#define __UFL_FORM_H_


namespace ufl
{

/**
 *  DOCUMENTATION:
 *
 *  @class  Form
 *
 *  @brief  Provides an interface complying with UFL Form.
 */

  class Form : public Expression
  {
    public:

      ///
      Form(const std::vector<Integral>& integrals);

      ///
      Form(const Integral& integral);

      ///
      ~Form();

      //--- INTERFACE -------------------------------------------------------------

      const Cell& cell() const;

      const std::vector<Integral>& integrals() const;

      const std::vector<Measure>& measures() const;

      const std::vector<Domain::Type>& domains() const;

      const std::vector<Integral>& cell_integrals() const;

      const std::vector<Integral>& exterior_facet_integrals() const;
      
      const std::vector<Integral>& interior_facet_integrals() const;
      
      const std::vector<Integral>& macro_cell_integrals() const;
      
      const std::vector<Integral>& surface_integrals() const;
      
      const FormData& form_data() const;
      
      const FormData& compute_form_data() const;

      bool is_preprocessed() const;
      
      //--- INTERFACE inherited from UFLClass -------------------------------------

      /// __repr__
      std::string const repr() const;

      /// __str__
      std::string const str() const;

    private:
      const std::vector<Integral> integrals;

      const FormData form_data;

      std::string const repr_;
      std::string const str_;

      bool is_preprocessed;
  };
} /* namespace ufl */
#endif /* __UFL_FORM_H_ */
