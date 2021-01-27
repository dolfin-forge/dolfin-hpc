// Copyright (C) 2018 Aurelien Larcher
// Licensed under the GNU LGPL Version 2.1.
#ifndef __DOLFIN_LOG_STREAM_H
#define __DOLFIN_LOG_STREAM_H

#include <dolfin/common/types.h>
#include <dolfin/common/assert.h>
#include <dolfin/common/maybe_unused.h>

#include <cstdio>
#include <sstream>
#include <ostream>
#include <iomanip>
#include <string>
#include <cstdarg>

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
  using stream = T;

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
  using stream = typename Sink::stream;

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

  ~__logstream() override = default;

  ///
  inline void init(stream * ss)
  {
    if (ss) Sink::init(*ss);
    std::fill_n(W_t_, sizeof(W_t_), ' ');
  }

  /// Pass through streambuf to honour indent, just declare supported types
  inline auto operator<<(short          x) -> __logstream& { std::ostream::operator<<(x); return *this; }
  inline auto operator<<(int            x) -> __logstream& { std::ostream::operator<<(x); return *this; }
  inline auto operator<<(long           x) -> __logstream& { std::ostream::operator<<(x); return *this; }
  inline auto operator<<(unsigned short x) -> __logstream& { std::ostream::operator<<(x); return *this; }
  inline auto operator<<(unsigned int   x) -> __logstream& { std::ostream::operator<<(x); return *this; }
  inline auto operator<<(unsigned long  x) -> __logstream& { std::ostream::operator<<(x); return *this; }
  inline auto operator<<(unsigned long long x) -> __logstream& { std::ostream::operator<<(x); return *this; }
  inline auto operator<<(float          x) -> __logstream& { std::ostream::operator<<(x); return *this; }
  inline auto operator<<(double         x) -> __logstream& { std::ostream::operator<<(x); return *this; }
  inline auto operator<<(long double    x) -> __logstream& { std::ostream::operator<<(x); return *this; }

  /// Declare char inserter explicitly
  inline auto operator<<(char  x) -> __logstream& { std::operator<<(*this, x); return *this;   }

  /// Pair output
  template<typename T, typename V>
  inline auto operator<<(std::pair<T, V> x) -> __logstream&
  {
    return (*this) << "( " << x.first << ", " << x.second << " )";
  }

  /// Output n char
  auto nputc(size_t n, char c) -> __logstream&
  {
    for (size_t i = 0; i < n; ++i) put(c);
    return *this;
  }

  /// Formatted output
  auto format(char const * msg, va_list aptr,
                      char const * pre = nullptr, char const * suf = nullptr,
                      size_t * n = nullptr) -> __logstream&
  {
    if (n)
    {
      // tellp does not work with any stream, switch to string buffer
      os_.clear(); os_.seekp(0);
      //sb_->rdbuf(os_.rdbuf());
    }
    if (pre)
      (*this) << pre;
    if (msg)
    {
      char const *c1 = msg;
      char fc = std::ostream::fill();
      int  pn = std::ostream::precision();
      int  wn = 0;
      for (;;)
      {
        while (*(c1 = msg++) != '%') {
          if (*c1 == '\0')
            goto ret;
          put(*c1);
        }
        bool sh = false;
        char ln = '\0';
        MAYBE_UNUSED(ln);
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
    if (suf)
      (*this) << suf;
    ss_->flush();
    if (n)
    {
      //sb_->rdbuf(ss_->rdbuf());
      *n = sb_->sputn(os_.str().c_str(), os_.str().size());
    }
    return *this;
  }

  /// Default implemented and filtered using the streambuf
  //  __logstream& operator<<(char const* x);
  //  __logstream& operator<<(std::string x);

  /// Manipulators
  auto operator<<(__logstream& (*fp)(__logstream&)) -> __logstream& { return fp(*this); }

  /// Manipulator adapters for stream implementation
  auto operator<<(stream& (*fp)(stream&)) -> __logstream& { fp(*this); return *this; }

  /// Manipulator adapters for ios_base
  auto operator<<(std::ios_base& (*fp)(std::ios_base&)) -> __logstream& { fp(*this); return *this; }

  auto operator++() -> __logstream& { if (W_i_ < INDENTMAX) ++W_i_; return *this; }
  auto operator--() -> __logstream& { if (W_i_ > INDENTMIN) --W_i_; return *this; }

  auto verbose() -> int { return W_v_; }

  auto verbose(int n) -> int
  {
    std::swap(W_v_, n);
    if (W_v_ < 0)
      clear(std::iostream::badbit);
    else
      clear(std::iostream::goodbit);
    return n;
  }

  auto silence() -> int
  {
    return verbose(-1);
  }

  auto indent() -> int { return W_i_; }

  auto indentwidth() -> int { return W_i_ * INDENTTAB; }

  auto indenttab() -> int { return INDENTTAB; }

private:

  ///
  stream * const         ss_;
  std::streambuf *       sb_;
  std::ostringstream     os_;

  static constexpr size_t LINEWIDTH = 256;
  static constexpr size_t INDENTTAB = 2;
  static constexpr size_t INDENTMIN = 0;
  static constexpr size_t INDENTMAX = 128;

  auto overflow(std::streambuf::int_type c) -> std::streambuf::int_type override
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
  size_t  W_i_;
  int  W_v_;
  char W_t_[LINEWIDTH];

};
//-----------------------------------------------------------------------------
using LogStream = __logstream<__sink<std::ostream> >;
//-----------------------------------------------------------------------------
/// Pair output
template<typename T, typename V>
inline auto operator<<(std::ostream& ss, std::pair<T, V> x) -> std::ostream&
{
  return ss << "( " << x.first << ", " << x.second << " )";
}
//-----------------------------------------------------------------------------

} /* namespace dolfin */

#endif /* __DOLFIN_LOG_STREAM_H */
