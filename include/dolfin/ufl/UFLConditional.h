// Copyright (C) 2014 Bärbel Janssen.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_UFL_CONDITIONAL_H
#define __DOLFIN_UFL_CONDITIONAL_H

#include <dolfin/ufl/UFLClass.h>
#include <dolfin/ufl/UFLExpression.h>

namespace ufl
{

/**
 *  DOCUMENTATION:
 *
 *  @class  Condition
 *
 *  @brief  Provides an interface complying with UFL Condition.
 */

class Condition : public Class
{

public:

  ///
  ~Condition() override;

  ///
  static Condition const * create(Object::repr_t const& repr);

  //--- INTERFACE -------------------------------------------------------------

//      ///
//      free_indices() const;

//      ///
//      index_dimensions() const;

//      ///
//      shape() const;

//      ///
//      evaluate(self, x, mapping, component, index_values):

//--- INTERFACE inherited from UFLClass -------------------------------------

//      /// __repr__
//      repr_t const& repr() const;
//
//      /// __str__
//      std::string const& str() const;
//
//      ///
//      void display() const;

protected:

  ///
  Condition(std::string const& name);

  ///
  Condition(std::string const& name, repr_t const& repr);

//      repr_t const repr_;
//      std::string const str_;
};

/**
 *  DOCUMENTATION:
 *
 *  @class  BinaryCondition
 *
 *  @brief  Provides an interface complying with UFL BinaryCondition.
 */

class BinaryCondition : public Condition
{

public:

  ///
  BinaryCondition(std::string const& name, Expression const& left_expression,
                  Expression const& right_expression);

  ///
  BinaryCondition(std::string const& name, repr_t const& repr);

  ///
  ~BinaryCondition() override;

  //--- INTERFACE inherited from UFLClass -------------------------------------

  /// __repr__
  repr_t const& repr() const override;

  /// __str__
  std::string const& str() const override;

  ///
  void display() const override;

private:

  std::vector<Expression const *> const fill_expressions(
      std::vector<repr_t> const& reprs);
  std::vector<Expression const *> const fill_expressions(Expression const& e1,
                                                         Expression const& e2);
  std::vector<Expression const *> const expressions_;

  repr_t const repr_;
  std::string const str_;
};

/**
 *  DOCUMENTATION:
 *
 *  @class  EQ
 *
 *  @brief  Provides an interface complying with UFL EQ.
 */

class EQ : public BinaryCondition
{

public:

  ///
  EQ(Expression const& left_expression, Expression const& right_expression);

  ///
  EQ(repr_t const& repr);

  ///
  ~EQ() override;

//      ///
//      evaluate(self, x, mapping, component, index_values):

};

/**
 *  DOCUMENTATION:
 *
 *  @class  NE
 *
 *  @brief  Provides an interface complying with UFL NE.
 */

class NE : public BinaryCondition
{

public:

  ///
  NE(Expression const& left_expression, Expression const& right_expression);

  ///
  NE(repr_t const& repr);

  ///
  ~NE() override;

//      ///
//      evaluate(self, x, mapping, component, index_values):

};

/**
 *  DOCUMENTATION:
 *
 *  @class  LE
 *
 *  @brief  Provides an interface complying with UFL LE.
 */

class LE : public BinaryCondition
{

public:

  ///
  LE(Expression const& left_expression, Expression const& right_expression);

  ///
  LE(repr_t const& repr);

  ///
  ~LE() override;

//      ///
//      evaluate(self, x, mapping, component, index_values):

};

/**
 *  DOCUMENTATION:
 *
 *  @class  GE
 *
 *  @brief  Provides an interface complying with UFL GE.
 */

class GE : public BinaryCondition
{

public:

  ///
  GE(Expression const& left_expression, Expression const& right_expression);

  ///
  GE(repr_t const& repr);

  ///
  ~GE() override;

//      ///
//      evaluate(self, x, mapping, component, index_values):

};

/**
 *  DOCUMENTATION:
 *
 *  @class  LT
 *
 *  @brief  Provides an interface complying with UFL LT.
 */

class LT : public BinaryCondition
{

public:

  ///
  LT(Expression const& left_expression, Expression const& right_expression);

  ///
  LT(repr_t const& repr);

  ///
  ~LT() override;

//      ///
//      evaluate(self, x, mapping, component, index_values):

};

/**
 *  DOCUMENTATION:
 *
 *  @class  GT
 *
 *  @brief  Provides an interface complying with UFL GT.
 */

class GT : public BinaryCondition
{

public:

  ///
  GT(Expression const& left_expression, Expression const& right_expression);

  ///
  GT(repr_t const& repr);

  ///
  ~GT() override;

//      ///
//      evaluate(self, x, mapping, component, index_values):

};

/**
 *  DOCUMENTATION:
 *
 *  @class  AndCondition
 *
 *  @brief  Provides an interface complying with UFL AndCondition.
 */

class AndCondition : public BinaryCondition
{

public:

  ///
  AndCondition(Expression const& left_expression,
               Expression const& right_expression);

  ///
  AndCondition(repr_t const& repr);

  ///
  ~AndCondition() override;

//      ///
//      evaluate(self, x, mapping, component, index_values):

};

/**
 *  DOCUMENTATION:
 *
 *  @class  OrCondition
 *
 *  @brief  Provides an interface complying with UFL OrCondition.
 */

class OrCondition : public BinaryCondition
{

public:

  ///
  OrCondition(Expression const& left_expression,
              Expression const& right_expression);

  ///
  OrCondition(repr_t const& repr);

  ///
  ~OrCondition() override;

//      ///
//      evaluate(self, x, mapping, component, index_values):

};

/**
 *  DOCUMENTATION:
 *
 *  @class  NotCondition
 *
 *  @brief  Provides an interface complying with UFL NotCondition.
 */

class NotCondition : public Condition
{

public:

  ///
  NotCondition(Expression const& e);

  ///
  NotCondition(repr_t const& repr);

  ///
  ~NotCondition() override;

//      ///
//      evaluate(self, x, mapping, component, index_values):

//--- INTERFACE inherited from UFLClass -------------------------------------

  /// __repr__
  repr_t const& repr() const override;

  /// __str__
  std::string const& str() const override;

  ///
  void display() const override;

private:

//      Expression const expression_;
  std::vector<Expression const *> const fill_expressions(
      std::vector<repr_t> const& reprs);
  std::vector<Expression const *> const fill_expressions(Expression const& e);
  std::vector<Expression const *> const expressions_;

  repr_t const repr_;
  std::string const str_;
};

/**
 *  DOCUMENTATION:
 *
 *  @class  Conditional
 *
 *  @brief  Provides an interface complying with UFL Conditional.
 */

class Conditional : public Class
{

public:

  ///
  Conditional(Condition const& c, Expression const& s1, Expression const& s2);

  ///
  Conditional(repr_t const& repr);

  ///
  ~Conditional() override;

  //--- INTERFACE -------------------------------------------------------------

//      ///
  std::vector<Expression const *> const& operands() const;

//      ///
//      free_indices() const;

//      ///
//      index_dimensions() const;

//      ///
//      shape() const;

//      ///
//      evaluate(self, x, mapping, component, index_values):

//--- INTERFACE inherited from UFLClass -------------------------------------

  /// __repr__
  repr_t const& repr() const override;

  /// __str__
  std::string const& str() const override;

  ///
  void display() const override;

private:

  Condition const * c_;
  std::vector<Expression const *> const fill_expressions(
      std::vector<repr_t> const& reprs);
  std::vector<Expression const *> const fill_expressions(Expression const& e1,
                                                         Expression const& e2);
  std::vector<Expression const *> const expressions_;
//      Expression const e1_;
//      Expression const e2_;

  repr_t const repr_;
  std::string const str_;
};

} /* namespace ufl */
#endif /* __DOLFIN_UFL_CONDITIONAL_H */
