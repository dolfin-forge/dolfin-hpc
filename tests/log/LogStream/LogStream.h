#include <dolfin_tests.h>

#ifdef HAVE_CHECK

#include <dolfin/log/log.h>
#include <dolfin/log/LogStream.h>

using namespace dolfin;

//-----------------------------------------------------------------------------
template<class T>
void check_log_stream(std::string, T& sb)
{
  __logstream<__sink<T> >  os(&sb);
  {
    os << '1';
    os << '2';
    os << '3';
    os << std::endl;
    os << '1' << '2' << '3';
  }
  os << std::endl;
  {
    short a = 1; os << a;
    short b = 2; os << b;
    short c = 3; os << c;
    os << std::endl;
    os << a << b << c;
  }
  os << std::endl;
  {
    int a = 1; os << a;
    int b = 2; os << b;
    int c = 3; os << c;
    os << std::endl;
    os << a << b << c;
  }
  os << std::endl;
  {
    long a = 1; os << a;
    long b = 2; os << b;
    long c = 3; os << c;
    os << std::endl;
    os << a << b << c;
  }
  os << std::endl;
  {
    uint a = 1; os << a;
    uint b = 2; os << b;
    uint c = 3; os << c;
    os << std::endl;
    os << a << b << c;
  }
  os << std::endl;
  {
    uidx a = 1; os << a;
    uidx b = 2; os << b;
    uidx c = 3; os << c;
    os << std::endl;
    os << a << b << c;
  }
  os << std::endl;
  {
    os << std::fixed;
    float a = 1.; os << a;
    float b = 2.; os << b;
    float c = 3.; os << c;
    os << std::endl;
    os << std::fixed;
    os << a << b << c;
    os << std::endl;
    os << std::scientific;
    os << a << b << c;
  }
  os << std::endl;
  {
    os << std::fixed;
    real a = 1.; os << a;
    real b = 2.; os << b;
    real c = 3.; os << c;
    os << std::endl;
    os << std::fixed;
    os << a << b << c;
    os << std::endl;
    os << std::scientific;
    os << a << b << c;
  }
  os << std::endl;
  {
    os << "charstar0";
    os << "charstar1";
    os << "charstar2";
  }
  os << std::endl;
}
//-----------------------------------------------------------------------------
DOLFIN_START_TEST( test_LogStream )
  {
    //---
    {
      check_log_stream<std::ostream>("cout", std::cout);
    }
    //---
    {
        dolfin::cout << "<html>\n";
      ++dolfin::cout << "<body>\n";
      ++dolfin::cout << "<content>\n";
      --dolfin::cout << "</body>\n";
      --dolfin::cout << "</html>\n";
    }
    //---
    {
      ++dolfin::cout;
      message("Characters: %c %c", 'a', 65);
      message("Decimals: %d %ld", 1977, 650000L);
      message("Preceding with blanks: %10d", 1977);
      message("Preceding with zeros: %010d", 1977);
      message("Some different radices: %d %x %o %#x %#o", 100, 100, 100, 100, 100);
      message("Floats: %4.2f %+.0e %E", 3.1416, 3.1416, 3.1416);
      message("Width trick: %*d", 5, 10);
      message("%s", "A string");
      char c = 0;
      message("Pointer value: %p", &c);
      --dolfin::cout;
    }
    //---
    {
      message(std::string("message"));
      message("message");
      message("message: %d %u", -1, 2);
      message("message: %f %e", -1.f, 2.f);
    }
    //---
    {
      warning(std::string("warning"));
      warning("warning");
    }
    //---
    {
      try {error(std::string("error"));} catch (...) {}
      try {error("error");} catch (...) {}
    }
    //---
    {
      begin("block0");
      message("a0");
      begin("block1");
      message("a1");
      begin("block3");
      message("a3\nb3\nc3");
      end();
      begin("block5");
      message("a5\nb5\nc5");
      end();
      message("b1");
      end();
      begin("block2");
      message("a2");
      begin("block4");
      message("a4\nb4\nc4");
      end();
      begin("block6");
      message("a6\nb6\nc6");
      end();
      message("b2");
      end();
      message("b0");
      end();
      //
      begin(std::string("blockS"));
      end();
    }
    //---
    {
      section("section0");
      cout << "This is the beginning of the hierarchy\n";
      section("section1");
      section("section1.1");
      end();
      end();
      section("section2");
      end();
      end();
      section(std::string("sectionS"));
      end();
    }
    //---
    {
      header("header0");
      header(std::string("headerS"));
    }
    //---
    {
      mark("");
      mark("time0");
    }
  }
DOLFIN_END_TEST
//-----------------------------------------------------------------------------

#endif
