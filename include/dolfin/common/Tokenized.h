// Copyright (C) 2015 Aurelien Larcher.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_TOKENIZED_H
#define __DOLFIN_TOKENIZED_H

#include <string>

namespace dolfin
{

/**
 *  @class    Tokenized
 *
 *  @brief    Implements an interface for tracking the internal state of an
 *            instance with a token.
 */

class Tokenized
{

public:
  ///
  virtual int token() const = 0;

protected:
  ///
  Tokenized()
  {
  }

  ///
  virtual ~Tokenized()
  {
  }

  ///
  virtual void update_token() = 0;
};

} /* namespace dolfin */

#endif /* __DOLFIN_TOKENIZED_H */
