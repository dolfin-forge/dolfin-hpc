AC_DEFUN([AX_CXX_HEADER_NONTR1_UNORDERED_MAP], [
  AC_CACHE_CHECK(for unordered_map,
  ax_cv_cxx_unordered_map,
  [AC_LANG_SAVE
  AC_LANG_CPLUSPLUS
  AC_TRY_COMPILE([#include <unordered_map>], [using std::unordered_map;],
  ax_cv_cxx_unordered_map=yes, ax_cv_cxx_unordered_map=no)
  AC_LANG_RESTORE
  ])
  if test "$ax_cv_cxx_unordered_map" = yes; then
    AC_DEFINE(HAVE_UNORDERED_MAP,[1],[Define if unordered_map is present. ])
  fi
])
