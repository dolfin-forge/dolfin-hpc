AC_DEFUN([AX_CXX_HEADER_NONTR1_UNORDERED_SET], [
  AC_CACHE_CHECK(for unordered_set,
  ax_cv_cxx_unordered_set,
  [AC_LANG_SAVE
  AC_LANG_CPLUSPLUS
  AC_TRY_COMPILE([#include <unordered_set>], [using std::unordered_set;],
  ax_cv_cxx_unordered_set=yes, ax_cv_cxx_unordered_set=no)
  AC_LANG_RESTORE
  ])
  if test "$ax_cv_cxx_unordered_set" = yes; then
    AC_DEFINE(HAVE_UNORDERED_SET,[1],[Define if unordered_set is present. ])
  fi
])
