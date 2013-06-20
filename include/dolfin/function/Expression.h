#ifndef __EXPRESSION_H
#define __EXPRESSION_H

#include <dolfin/common/types.h>

namespace dolfin
{

/// Derived classes should implement the copy constructor if the allocate
/// new objects.

class Expression
{
public:

	/// Create user-defined function
	Expression()
	{
	}

	/// Return the rank of the value space
	virtual uint rank() const = 0;

	/// Return the dimension of the value space for axis i
	virtual uint dim(uint i) const = 0;

	/// Evaluate function at given point
	virtual void eval(real* values, const real* x) const = 0;

protected:

	/// Destructor
	virtual ~Expression()
	{
	}

private:

};

class RealReference: public Expression
{
public:

	/// Create user-defined function
	RealReference(real const& r) :
			Expression(),
			r_(r)
	{
	}

	/// Destructor
	~RealReference()
	{
	}

	/// Return the rank of the value space
	virtual uint rank() const
	{
		return 0;
	}

	/// Return the dimension of the value space for axis i
	virtual uint dim(uint i) const
	{
		return 1;
	}

	/// Evaluate function at given point
	virtual void eval(real* values, const real* x) const
	{
		values[0] = r_;
	}

private:
	real const& r_;

};

}

#endif
