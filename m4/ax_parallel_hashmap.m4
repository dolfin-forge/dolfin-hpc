#
# ----------------------------------------------------------------------------
# "THE BEER-WARE LICENSE" (Revision 42):
# <julian@ingridcloud.com> wrote this file. As long as you retain this notice
# you can do whatever you want with this stuff. If we meet some day, and you
# think this stuff is worth it, you can buy me a beer in return Julian Hornich
# ----------------------------------------------------------------------------
#

AC_DEFUN([AX_PHMAP],[
	AC_ARG_WITH([phmap],
	AS_HELP_STRING([--with-phmap=DIR],
	[Directory for parallel hashmap: parallel_hashmap/phmap.h]),
	[
	if test -d "$withval"; then
		ac_phmap_path="$withval/parallel_hashmap";
		PHMAP_CPPFLAGS="-I$ac_phmap_path"
	fi
	],)

	if test -d "$ac_phmap_path"; then
		CPPFLAGS_SAVED="$CPPFLAGS"
		CPPFLAGS="$PHMAP_CPPFLAGS $CPPFLAGS"
		export CPPFLAGS
	fi

	have_phmap_h=no
	if test "x${with_phmap}" != xno; then
		AC_CHECK_HEADER([parallel_hashmap/phmap.h],[have_phmap_h=yes],[have_phmap_h=no])
		if test x"${have_phmap_h}" = xyes; then
			AC_DEFINE(HAVE_PARALLEL_HASH_MAP,[1],[Define if parallel hash map is present])
		else
			if test -d "$ac_phmap_path"; then
				CPPFLAGS="$CPPFLAGS_SAVED"
			fi
		fi
	fi
])
