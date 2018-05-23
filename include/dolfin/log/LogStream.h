// Copyright (C) 2018 Aurelien Larcher
// Licensed under the GNU LGPL Version 2.1.
//
#ifndef __DOLFIN_LOG_STREAM_H
#define __DOLFIN_LOG_STREAM_H

#include <dolfin/common/types.h>
#include <dolfin/common/assert.h>

#include <cstdio>
#include <ostream>
#include <iomanip>
#include <string>
#include <cstring>
#include <stdarg.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
struct __log_pipe
{
  __log_pipe(std::ostream & src) :
    src_(src),
    buf_(src.rdbuf())
  {
  }

  void operator()(std::ostream & dst)
  {
    buf_ = src_.rdbuf(dst.rdbuf());
  }

  ~__log_pipe()
  {
    src_.rdbuf(buf_);
  }

private:

    std::ostream &   src_;
    std::streambuf * buf_;

};
//-----------------------------------------------------------------------------
template<class T>
struct __sink
{
  typedef T stream;

  static inline void init(stream& ss)
  {
    if (!ss.good()) ss.clear();
    if (!ss.good()) std::perror("Invalid stream");
  }
};
//-----------------------------------------------------------------------------
template<class Sink>
struct __logstream : protected std::streambuf, public std::ostream
{
  typedef typename Sink::stream stream;

  __logstream(stream * ss) :
      std::streambuf(),
      std::ostream(this), //XXX
      ss_(ss),
      sb_(ss->rdbuf()),
      W_c_('\n'),
      W_i_(0)
  {
    init(ss);
  }

  ~__logstream()
  {
  }

  ///
  inline void init(stream * ss)
  {
    if (ss) Sink::init(*ss);
    std::fill_n(W_t_, sizeof(W_t_), ' ');
  }

  /// Pass through streambuf to honour indent, just declare supported types
  inline __logstream& operator<<(short          x) { std::ostream::operator<<(x); return *this; }
  inline __logstream& operator<<(int            x) { std::ostream::operator<<(x); return *this; }
  inline __logstream& operator<<(long           x) { std::ostream::operator<<(x); return *this; }
  inline __logstream& operator<<(unsigned short x) { std::ostream::operator<<(x); return *this; }
  inline __logstream& operator<<(unsigned int   x) { std::ostream::operator<<(x); return *this; }
  inline __logstream& operator<<(unsigned long  x) { std::ostream::operator<<(x); return *this; }
  inline __logstream& operator<<(float          x) { std::ostream::operator<<(x); return *this; }
  inline __logstream& operator<<(double         x) { std::ostream::operator<<(x); return *this; }
  inline __logstream& operator<<(long double    x) { std::ostream::operator<<(x); return *this; }

  /// Declare char inserter explicitly
  inline __logstream& operator<<(char  x) { std::operator<<(*this, x); return *this;   }

  /// Pair output
  template<typename T, typename V>
  inline __logstream& operator<<(std::pair<T, V> x)
  {
    return (*this) << "( " << x.first << ", " << x.second << " )";
  }

  /// Output n char
  __logstream& nputc(uint n, char c)
  {
    for (uint i = 0; i < n; ++i) put(c);
    return *this;
  }

  /// Formatted output
  __logstream& format(char const * msg, va_list aptr,
                      char const * pre = NULL, char const * suf = NULL,
                      size_t * n = NULL)
  {
    if (n)
    {
      // tellp does not work with any stream, switch to string buffer
      sb_ = os_.rdbuf(); os_.clear(); os_.seekp(0);
    }
    if (pre) (*this) << pre;
    if (msg)
    {
      char const *c0, *c1 = msg;
      char fc = std::ostream::fill();
      int  pn = std::ostream::precision();
      int  wn = 0;
      for (;;)
      {
        c0 = c1;
        while (*(c1 = msg++) != '%') { if (*c1 == '\0') goto ret; put(*c1); }
        bool sh = false; char ln = '\0';
        // Backup ioflags
        std::ios_base::fmtflags ff(flags());
fmt:
        /**
          * Reduced C89 format: %[flags][width][.precision][length]specifier
          *
          * flags:
          *   - left justify
          *   + show positive sign
          *   # format modifier: showbase (integer), showpoint (floating-point)
          *   0 zero padding
          *
          * width:
          *   [0-9] numeric value
          *   *     read from va_arg
          *
          * precision:
          *   [0-9] numeric value
          *   *     read from va_arg
          *
          * length:
          *   h short integer
          *   l long integer
          *   L long double
          *
          * specifiers:
          *   d Signed decimal integer
          *   i Signed decimal integer
          *   u Unsigned decimal integer
          *   o Unsigned octal
          *   x Unsigned hexadecimal integer
          *   X Unsigned hexadecimal integer (uppercase)
          *   f Decimal floating point, lowercase
          *   F Decimal floating point, uppercase
          *   e Scientific notation (mantissa/exponent), lowercase
          *   E Scientific notation (mantissa/exponent), uppercase
          *   g Use %f
          *   G Use %F
          *   a Hexadecimal floating point, lowercase
          *   A Hexadecimal floating point, uppercase
          *   c Character a
          *   s String of characters
          *   p Pointer address
          *   % Percent character
          */
        switch(*(c1 = msg++))
        {
          // flags
          case '-': // left-justify
            std::ostream::setf(std::ios_base::left, std::ios::adjustfield);
            goto fmt;
          case '+': // force sign
            std::ostream::setf(std::ios_base::showpos);
            goto fmt;
          case ' ': // space in place of sign
            //XXX: not supported
            goto fmt;
          case '#': // format modifier
            sh = true;
            goto fmt;
          case '0': // left padding
            fc = std::ostream::fill('0');
            goto fmt;
          // width
          case '1':
          case '2':
          case '3':
          case '4':
          case '5':
          case '6':
          case '7':
          case '8':
          case '9':
            for (wn = 0;(*c1 >= '0') && (*c1 <= '9'); ++c1)
            {
              wn = wn * 10 + (*c1 - '0');
            }
            (*this) << std::setw(wn);
            msg = c1;
            goto fmt;
          case '*':
            (*this) << std::setw(va_arg(aptr, int));
            goto fmt;
          // precision
          case '.':
            if (*(++c1) == '*')
            {
              wn = va_arg(aptr, int); ++c1;
            }
            else
            {
              for (wn = 0; (*c1 >= '0') && (*c1 <= '9'); ++c1)
              {
                wn = wn * 10 + (*c1 - '0');
              }
            }
            pn = std::ostream::precision(wn);
            msg = c1;
            goto fmt;
          // length
          case 'h':
          case 'l':
          case 'L':
            ln = *c1;
            goto fmt;
          // specifier
          case 'd':
          case 'i':
            std::ostream::setf(std::ios_base::dec, std::ios_base::basefield);
            if (*c1 == 'l')
              std::ostream::operator<<(va_arg(aptr, long int));
            else
              std::ostream::operator<<(va_arg(aptr, int));
            break;
          case 'u':
            std::ostream::setf(std::ios_base::dec, std::ios_base::basefield);
fmtu:
            if (sh) { std::ostream::setf(std::ios_base::showbase); }
            if (*c1 == 'l')
              std::ostream::operator<<(va_arg(aptr, unsigned long int));
            else
              std::ostream::operator<<(va_arg(aptr, unsigned int));
            break;
          case 'o':
            std::ostream::setf(std::ios_base::oct, std::ios_base::basefield);
            goto fmtu;
          case 'X':
            std::ostream::setf(std::ios_base::uppercase);
            /* FALLTHROUGH */
          case 'x':
            std::ostream::setf(std::ios_base::hex, std::ios_base::basefield);
            goto fmtu;
          case 'F':
          case 'G':
            std::ostream::setf(std::ios_base::uppercase);
            /* FALLTHROUGH */
          case 'f':
          case 'g':
            std::ostream::setf(std::ios_base::fixed, std::ios::floatfield);
fmtf:
            if (sh) std::ostream::setf(std::ios_base::showpoint);
            if (*c1 == 'L')
              std::ostream::operator<<(va_arg(aptr, long double));
            else
              std::ostream::operator<<(va_arg(aptr, double));
            break;
          case 'E':
            std::ostream::setf(std::ios_base::uppercase);
            /* FALLTHROUGH */
          case 'e':
            setf(std::ios_base::scientific, std::ios::floatfield);
            goto fmtf;
          case 'A':
            std::ostream::setf(std::ios_base::uppercase);
            /* FALLTHROUGH */
          case 'a':
            std::ostream::setf(std::ios_base::scientific, std::ios::floatfield);
            std::ostream::setf(std::ios_base::hex, std::ios_base::basefield);
            goto fmtf;
          case 'c':
            // Use non-member char inserter
            (*this) << static_cast<char>(va_arg(aptr, int));
            break;
          case 's':
            // Use non-member string inserter
            (*this) << static_cast<char *>(va_arg(aptr, char*));
            break;
          case 'p':
            std::ostream::operator<<(va_arg(aptr, void*));
            break;
          case '%':
            std::ostream::operator<<('%');
            break;
          // default
          default:
            break;
        }

        // Reset
        fc = std::ostream::fill(fc);
        pn = std::ostream::precision(pn);
        flags(ff);
      }
    }
ret:
    if (suf) (*this) << suf;
    if (n)
    {
      sb_ = ss_->rdbuf();
      *n = sb_->sputn(os_.str().c_str(), os_.str().size());
    }
    return *this;
  }

  /// Default implemented and filtered using the streambuf
  //  __logstream& operator<<(char const* x);
  //  __logstream& operator<<(std::string x);

  /// Manipulators
  __logstream& operator<<(__logstream& (*fp)(__logstream&)) { return fp(*this); }

  /// Manipulator adapters for stream implementation
  __logstream& operator<<(stream& (*fp)(stream&)) { fp(*this); return *this; }

  /// Manipulator adapters for ios_base
  __logstream& operator<<(std::ios_base& (*fp)(std::ios_base&)) { fp(*this); return *this; }

  __logstream& operator++() { if (W_i_ < INDENTMAX) ++W_i_; return *this; }
  __logstream& operator--() { if (W_i_ > INDENTMIN) --W_i_; return *this; }

  int verbose() { return W_v_; }

  int verbose(int n)
  {
    std::swap(W_v_, n);
    if (W_v_ < 0)
      clear(std::iostream::badbit);
    else
      clear(std::iostream::goodbit);
    return n;
  }

  int silence()
  {
    return verbose(-1);
  }

private:

  ///
  stream * const         ss_;
  std::streambuf *       sb_;
  std::ostringstream     os_;

  static int const LINEWIDTH = 256;
  static int const INDENTTAB = 2;
  static int const INDENTMIN = 0;
  static int const INDENTMAX = 128;

  std::streambuf::int_type overflow(std::streambuf::int_type c)
  {
    dolfin_assert(sb_);
    if (std::streambuf::traits_type::eq_int_type(std::streambuf::traits_type::eof(), c))
    {
      return std::streambuf::traits_type::not_eof(c);
    }
    if (W_i_)
    {
      switch (W_c_)
      {
        case '\n':
        case '\r':
          sb_->sputn(W_t_, INDENTTAB * W_i_);
          break;
        default:
          break;
      }
    }
    sb_->sputc(c);
    W_c_ = c;
    return c;
  }

  char W_c_;
  int  W_i_;
  int  W_v_;
  char W_t_[LINEWIDTH];

};
//-----------------------------------------------------------------------------
typedef __logstream<__sink<std::ostream> > LogStream;
//-----------------------------------------------------------------------------
/// Pair output
template<typename T, typename V>
inline std::ostream& operator<<(std::ostream& ss, std::pair<T, V> x)
{
  return ss << "( " << x.first << ", " << x.second << " )";
}
//-----------------------------------------------------------------------------

} /* namespace dolfin */

#endif /* __DOLFIN_LOG_STREAM_H */
