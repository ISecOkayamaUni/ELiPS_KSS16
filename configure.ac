#Autoconf pre-req version is 2.69
AC_PREREQ(2.69)
#AC_INIT(FULL-PACKAGE-NAME, VERSION, BUG-REPORT-ADDRESS)
AC_INIT([elipsbnbls12], [1.0.1])
#Automake init
AM_INIT_AUTOMAKE(elipsbnbls12, 1.0.1, no-define)
AC_CONFIG_SRCDIR([src/bn_inits.c])
AC_CONFIG_MACRO_DIRS([m4])
#Initialize libtool with min ver 2.4
LT_INIT
LT_PREREQ([2.4])

#Avoid warning of README
AM_INIT_AUTOMAKE([foreign 1.13 subdir-objects])



# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_CPP
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_MAKE_SET
AC_PROG_LIBTOOL


# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([stdlib.h string.h gmp.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_TYPE_SIZE_T
AC_CHECK_SIZEOF(char)
AC_CHECK_SIZEOF(short int)
AC_CHECK_SIZEOF(int)
AC_CHECK_SIZEOF(long)
AC_CHECK_SIZEOF(long long)


# Checks for libraries.
# Checks for library functions.
AC_FUNC_MALLOC
AC_CHECK_LIB(gmp, __gmpz_add)
AC_CONFIG_HEADERS([config.h])

#make files
AC_CONFIG_FILES([Makefile
                 src/Makefile])
AC_OUTPUT
