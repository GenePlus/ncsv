# Process this file with autoconf to produce a configure script.
AC_PREREQ([2.69])
AC_INIT([sv_helper], [1.0], [liuchang@geneplus.org.cn])
AM_INIT_AUTOMAKE([foreign])
AC_CONFIG_HEADER([config.h])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_RANLIB

# Check for libraries
AC_SEARCH_LIBS([gzopen],[z],,[AC_MSG_ERROR([libz not found, please install zlib (http://www.zlib.net/)])])
AC_SEARCH_LIBS([clock_gettime], [rt], [AC_DEFINE([HAVE_CLOCK_GETTIME], [1], [clock_getttime found])], )


# Checks for header files.
AC_LANG([C++])
AC_CHECK_HEADERS([ \
		   unordered_map tr1/unordered_map ext/hash_map \
		   unordered_set tr1/unordered_set ext/hash_set \
])

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_INT32_T
AC_TYPE_INT64_T
AC_TYPE_SIZE_T
AC_TYPE_UINT16_T
AC_TYPE_UINT32_T
AC_TYPE_UINT64_T
AC_TYPE_UINT8_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_CHECK_FUNCS([clock_gettime floor gettimeofday memset pow])

AC_SUBST(CXXFLAGS, "$CXXFLAGS -fpic")
AC_SUBST(CFLAGS, "$CFLAGS -fpic")

AC_CONFIG_FILES([Makefile
                 SGA/Algorithm/Makefile
                 SGA/Bigraph/Makefile
                 SGA/SGA/Makefile
                 SGA/SQG/Makefile
                 SGA/StringGraph/Makefile
                 SGA/SuffixTools/Makefile
                 SGA/Util/Makefile])
AC_CONFIG_SUBDIRS([SeqLib])

AC_OUTPUT
