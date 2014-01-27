// Copyright (C) 2014 Bärbel Janssen.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:
// Last changed:

#ifndef __UFL_EXPRESSION_H_
#define __UFL_EXPRESSION_H_


namespace ufl
{

/**
 *  DOCUMENTATION:
 *
 *  @class  Expression
 *
 *  @brief  Provides an interface complying with UFL Expression.
 */

  class Expression : public Class
  {
    public:

      ///
      Expression();
    
      ///
      ~Expression();

      //--- INTERFACE -------------------------------------------------------------

      Expression& reconstruct(const Expression &expression);

      operands();
      
      shape();
      
      rank();
      
      cell();
      
      geometric_dimension();
      
      is_cellwise_constant();
      
      evaluate(const Coordinate& x, const Mapping& mapping, const IndexValues& index_values);
      
      free_indices();
      
      index_dimensions();
      
      hash();


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

  class Operator : public Expression
  {
    public:

      Operator (const Expression& expression);
      ~Operator ();

      const Expression& reconstruct(const Expression& expression);
      bool is_cellwise_constant() const;
  
  };
} /* namespace ufl */
#endif /* __UFL_EXPRESSION_H_ */
