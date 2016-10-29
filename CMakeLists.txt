CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
PROJECT(POSEST)

# levmar library
SET(LEVMAR_INCLUDE_DIR "/usr/local/levmar-2.6" CACHE PATH "Path to LEVMAR library header")
#ADD_SUBDIRECTORY(${LEVMAR_INCLUDE_DIR})
INCLUDE_DIRECTORIES(${LEVMAR_INCLUDE_DIR})
MESSAGE(STATUS "LEVMAR_INCLUDE_DIR = ${LEVMAR_INCLUDE_DIR}")

IF(MSVC)
  ADD_DEFINITIONS(/arch:SSE2)
  # get rid of CRT warnings
  ADD_DEFINITIONS(-D_CRT_SECURE_NO_WARNINGS)
ENDIF(MSVC)

INCLUDE_DIRECTORIES(mlsl)

ADD_LIBRARY(posest buckets.c lqs.c ransac.c rngs.c prosac.c deal.c p3p.c p4pf.c planep4p.c svd3.c polysolve.c
  poseproj.c posest.c align.c lhm.c sam.c compiler.h lqs.h ransac.h prosac.h p3p.h p4pf.h planep4p.h svd3.h
  polysolve.h poseproj.h posest.h util.h rngs.h sam.h
  mlsl/mlsl.c mlsl/mt19937ar.c mlsl/redblack.c mlsl/sobolseq.c
)

OPTION(BUILD_DEMOS "Build demo programs?" TRUE)

# demo program
IF(BUILD_DEMOS)
  # paths to levmar & lapack/blas libraries
  SET(LEVMAR_BINARY_DIR "${LEVMAR_INCLUDE_DIR}/w32" CACHE PATH "Path to levmar library")
  MESSAGE(STATUS "LEVMAR_BINARY_DIR = ${LEVMAR_BINARY_DIR}")

  SET(LAPACKBLAS_DIR "/usr/lib" CACHE PATH "Path to lapack/blas libraries")
  MESSAGE(STATUS "LAPACKBLAS_DIR = ${LAPACKBLAS_DIR}")


  # actual names for the lapack/blas/f2c libraries
  SET(LAPACKBLAS_LIB_NAMES "lapack;blas" CACHE STRING "The name of the lapack & blas libraries")
  SET(F2C_LIB_NAME f2c CACHE STRING "The name of the f2c or F77/I77 library")
  # f2c is sometimes equivalent to libF77 & libI77
  #SET(F2C_LIB_NAME "libF77;libI77" CACHE STRING "The name of the f2c or F77/I77 library")

  SET(LIBS posest levmar)
  SET(LIBS ${LIBS} ${LAPACKBLAS_LIB_NAMES} ${F2C_LIB_NAME})

  LINK_DIRECTORIES(${CMAKE_BINARY_DIR}) # location of the posest library
  LINK_DIRECTORIES(${LEVMAR_BINARY_DIR}) # location of the levmar library
  LINK_DIRECTORIES(${LAPACKBLAS_DIR})

  ADD_EXECUTABLE(posest_demo posest_demo.c)
  TARGET_LINK_LIBRARIES(posest_demo ${LIBS})

  ADD_EXECUTABLE(binocposest_demo binocposest_demo.c)
  TARGET_LINK_LIBRARIES(binocposest_demo ${LIBS})

  MESSAGE(STATUS "demos will be linked against ${LIBS}")

# make sure that the libraries are built before the demos
  ADD_DEPENDENCIES(posest_demo posest levmar)
  ADD_DEPENDENCIES(binocposest_demo posest levmar)
ENDIF(BUILD_DEMOS)
