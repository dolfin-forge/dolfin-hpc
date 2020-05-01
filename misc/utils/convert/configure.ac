AC_PREREQ([2.61])
AC_INIT([DOLFINconvert], [1.0.0])
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE([foreign])
AM_MAINTAINER_MODE
AC_LANG(C++)
AC_ARG_ENABLE(debug, AC_HELP_STRING([--enable-debug],
             [Compile with debug symbols]),
             [enable_debug=${enableval}],[enable_debug=no])

AM_PATH_XML2(,,AC_MSG_ERROR([libxml2 not found]))
CPPFLAGS="$CPPFLAGS $XML_CPPFLAGS"
LIBS="$XML_LIBS $LIBS"

AC_CONFIG_FILES([Makefile\
		 src/Makefile])

AC_OUTPUT
