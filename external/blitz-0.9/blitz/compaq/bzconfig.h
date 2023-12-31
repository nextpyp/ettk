#ifndef _BLITZ_COMPAQ_BZCONFIG_H
#define _BLITZ_COMPAQ_BZCONFIG_H 1
 
/* blitz/compaq/bzconfig.h. Generated automatically at end of configure. */
/* blitz/config.h.  Generated by configure.  */
/* blitz/config.h.in.  Generated from configure.ac by autoheader.  */


/******************************************************************************
 * config.h           Compiler language support flags
 *
 * This file was generated automatically when running the configure script.
 * You should rerun configure each time you switch compilers, install new
 * standard libraries, or change compiler versions.
 *
 */



/* define if bool is a built-in type */
/* #undef BZ_HAVE_BOOL */

/* define if the compiler has <climits> header */
/* #undef BZ_HAVE_CLIMITS */

/* define if the compiler has complex<T> */
/* #undef BZ_HAVE_COMPLEX */

/* define if the compiler has standard complex<T> functions */
/* #undef BZ_HAVE_COMPLEX_FCNS */

/* define if the compiler has complex math functions */
/* #undef BZ_HAVE_COMPLEX_MATH1 */

/* define if the compiler has more complex math functions */
/* #undef BZ_HAVE_COMPLEX_MATH2 */

/* define if complex math functions are in namespace std */
/* #undef BZ_HAVE_COMPLEX_MATH_IN_NAMESPACE_STD */

/* define if the compiler supports const_cast<> */
/* #undef BZ_HAVE_CONST_CAST */

/* Define to 1 if you have the <cstring> header file. */
#ifndef BZ_HAVE_CSTRING 
#define BZ_HAVE_CSTRING  1 
#endif

/* define if the compiler supports default template parameters */
/* #undef BZ_HAVE_DEFAULT_TEMPLATE_PARAMETERS */

/* Define to 1 if you have the <dlfcn.h> header file. */
#ifndef BZ_HAVE_DLFCN_H 
#define BZ_HAVE_DLFCN_H  1 
#endif

/* define if the compiler supports dynamic_cast<> */
/* #undef BZ_HAVE_DYNAMIC_CAST */

/* define if the compiler handle computations inside an enum */
/* #undef BZ_HAVE_ENUM_COMPUTATIONS */

/* define if the compiler handles (int) casts in enum computations */
/* #undef BZ_HAVE_ENUM_COMPUTATIONS_WITH_CAST */

/* define if the compiler supports exceptions */
/* #undef BZ_HAVE_EXCEPTIONS */

/* define if the compiler supports the explicit keyword */
/* #undef BZ_HAVE_EXPLICIT */

/* define if the compiler supports explicit template function qualification */
/* #undef BZ_HAVE_EXPLICIT_TEMPLATE_FUNCTION_QUALIFICATION */

/* define if the compiler recognizes the full specialization syntax */
/* #undef BZ_HAVE_FULL_SPECIALIZATION_SYNTAX */

/* define if the compiler supports function templates with non-type parameters
   */
/* #undef BZ_HAVE_FUNCTION_NONTYPE_PARAMETERS */

/* define if the compiler supports IEEE math library */
/* #undef BZ_HAVE_IEEE_MATH */

/* Define to 1 if you have the <inttypes.h> header file. */
#ifndef BZ_HAVE_INTTYPES_H 
#define BZ_HAVE_INTTYPES_H  1 
#endif

/* Define to 1 if you have the `m' library (-lm). */
/* #undef BZ_HAVE_LIBM */

/* define if the compiler supports member constants */
/* #undef BZ_HAVE_MEMBER_CONSTANTS */

/* define if the compiler supports member templates */
/* #undef BZ_HAVE_MEMBER_TEMPLATES */

/* define if the compiler supports member templates outside the class
   declaration */
/* #undef BZ_HAVE_MEMBER_TEMPLATES_OUTSIDE_CLASS */

/* Define to 1 if you have the <memory.h> header file. */
#ifndef BZ_HAVE_MEMORY_H 
#define BZ_HAVE_MEMORY_H  1 
#endif

/* define if the compiler supports the mutable keyword */
/* #undef BZ_HAVE_MUTABLE */

/* define if the compiler implements namespaces */
/* #undef BZ_HAVE_NAMESPACES */

/* define if the compiler supports the Numerical C Extensions Group restrict
   keyword */
/* #undef BZ_HAVE_NCEG_RESTRICT */

/* define if the compiler supports the __restrict__ keyword */
/* #undef BZ_HAVE_NCEG_RESTRICT_EGCS */

/* define if the compiler has numeric_limits<T> */
/* #undef BZ_HAVE_NUMERIC_LIMITS */

/* define if the compiler accepts the old for scoping rules */
/* #undef BZ_HAVE_OLD_FOR_SCOPING */

/* define if the compiler supports partial ordering */
/* #undef BZ_HAVE_PARTIAL_ORDERING */

/* define if the compiler supports partial specialization */
/* #undef BZ_HAVE_PARTIAL_SPECIALIZATION */

/* define if the compiler supports reinterpret_cast<> */
/* #undef BZ_HAVE_REINTERPRET_CAST */

/* define if the compiler supports Run-Time Type Identification */
/* #undef BZ_HAVE_RTTI */

/* define if the compiler has getrusage() function */
/* #undef BZ_HAVE_RUSAGE */

/* define if the compiler supports static_cast<> */
/* #undef BZ_HAVE_STATIC_CAST */

/* define if the compiler supports ISO C++ standard library */
/* #undef BZ_HAVE_STD */

/* Define to 1 if you have the <stdint.h> header file. */
#ifndef BZ_HAVE_STDINT_H 
#define BZ_HAVE_STDINT_H  1 
#endif

/* Define to 1 if you have the <stdlib.h> header file. */
#ifndef BZ_HAVE_STDLIB_H 
#define BZ_HAVE_STDLIB_H  1 
#endif

/* define if the compiler supports Standard Template Library */
/* #undef BZ_HAVE_STL */

/* Define to 1 if you have the <strings.h> header file. */
#ifndef BZ_HAVE_STRINGS_H 
#define BZ_HAVE_STRINGS_H  1 
#endif

/* Define to 1 if you have the <string.h> header file. */
#ifndef BZ_HAVE_STRING_H 
#define BZ_HAVE_STRING_H  1 
#endif

/* define if the compiler supports System V math library */
/* #undef BZ_HAVE_SYSTEM_V_MATH */

/* Define to 1 if you have the <sys/stat.h> header file. */
#ifndef BZ_HAVE_SYS_STAT_H 
#define BZ_HAVE_SYS_STAT_H  1 
#endif

/* Define to 1 if you have the <sys/types.h> header file. */
#ifndef BZ_HAVE_SYS_TYPES_H 
#define BZ_HAVE_SYS_TYPES_H  1 
#endif

/* define if the compiler supports basic templates */
/* #undef BZ_HAVE_TEMPLATES */

/* define if the compiler supports templates as template arguments */
/* #undef BZ_HAVE_TEMPLATES_AS_TEMPLATE_ARGUMENTS */

/* define if the compiler supports use of the template keyword as a qualifier
   */
/* #undef BZ_HAVE_TEMPLATE_KEYWORD_QUALIFIER */

/* define if the compiler supports template-qualified base class specifiers */
/* #undef BZ_HAVE_TEMPLATE_QUALIFIED_BASE_CLASS */

/* define if the compiler supports template-qualified return types */
/* #undef BZ_HAVE_TEMPLATE_QUALIFIED_RETURN_TYPE */

/* define if the compiler supports function matching with argument types which
   are template scope-qualified */
/* #undef BZ_HAVE_TEMPLATE_SCOPED_ARGUMENT_MATCHING */

/* define if the compiler recognizes typename */
/* #undef BZ_HAVE_TYPENAME */

/* define if the compiler supports the vector type promotion mechanism */
/* #undef BZ_HAVE_TYPE_PROMOTION */

/* Define to 1 if you have the <unistd.h> header file. */
#ifndef BZ_HAVE_UNISTD_H 
#define BZ_HAVE_UNISTD_H  1 
#endif

/* define if the compiler supports numeric traits promotions */
/* #undef BZ_HAVE_USE_NUMTRAIT */

/* define if the compiler has valarray<T> */
/* #undef BZ_HAVE_VALARRAY */

/* define if the compiler has isnan function in namespace std */
/* #undef BZ_ISNAN_IN_NAMESPACE_STD */

/* define if the compiler has C math abs(integer types) in namespace std */
/* #undef BZ_MATH_ABSINT_IN_NAMESPACE_STD */

/* define if the compiler has C math functions in namespace std */
/* #undef BZ_MATH_FN_IN_NAMESPACE_STD */

/* Name of package */
#ifndef BZ_PACKAGE 
#define BZ_PACKAGE  "blitz" 
#endif

/* Define to the address where bug reports for this package should be sent. */
#ifndef BZ_PACKAGE_BUGREPORT 
#define BZ_PACKAGE_BUGREPORT  "blitz-bugs@oonumerics.org" 
#endif

/* Define to the full name of this package. */
#ifndef BZ_PACKAGE_NAME 
#define BZ_PACKAGE_NAME  "blitz" 
#endif

/* Define to the full name and version of this package. */
#ifndef BZ_PACKAGE_STRING 
#define BZ_PACKAGE_STRING  "blitz 0.9" 
#endif

/* Define to the one symbol short name of this package. */
#ifndef BZ_PACKAGE_TARNAME 
#define BZ_PACKAGE_TARNAME  "blitz" 
#endif

/* Define to the version of this package. */
#ifndef BZ_PACKAGE_VERSION 
#define BZ_PACKAGE_VERSION  "0.9" 
#endif

/* Define to 1 if you have the ANSI C header files. */
#ifndef BZ_STDC_HEADERS 
#define BZ_STDC_HEADERS  1 
#endif

/* Enable Blitz thread-safety features */
/* #undef BZ_THREADSAFE */

/* Version number of package */
#ifndef BZ_VERSION 
#define BZ_VERSION  "0.9" 
#endif

/* CXX */
#ifndef BZ__compiler_name 
#define BZ__compiler_name  "mpicxx" 
#endif

/* CXXFLAGS */
#ifndef BZ__compiler_options 
#define BZ__compiler_options  "-std ansi -D__USE_STD_IOSTREAM -DBZ_ENABLE_XOPEN_SOURCE -D_OSF_SOURCE -ieee -model ansi -accept restrict_keyword -nousing_std" 
#endif

/* date */
#ifndef BZ__config_date 
#define BZ__config_date  "Tue Sep 22 20:56:51 EDT 2015" 
#endif

/* uname -a */
#ifndef BZ__os_name 
#define BZ__os_name  "Linux biowulf2.nih.gov 2.6.32-573.3.1.el6.x86_64 #1 SMP Mon Aug 10 09:44:54 EDT 2015 x86_64 x86_64 x86_64 GNU/Linux" 
#endif

/* target */
#ifndef BZ__platform 
#define BZ__platform  "x86_64-unknown-linux-gnu" 
#endif
 
/* once: _BLITZ_COMPAQ_BZCONFIG_H */
#endif
