cmake_minimum_required(VERSION 3.13)

include(CheckIncludeFile)

# set the project name
project(CAFE VERSION 1.1)

set (CMAKE_PREFIX_PATH      "${PROJECT_SOURCE_DIR}/cmake")

set(OPTIMIZER_HIGH_PRECISION 1e-6)
set(OPTIMIZER_LOW_PRECISION 1e-3)
set(PHASED_OPTIMIZER_PHASE1_ATTEMPTS 4)
set(NUM_OPTIMIZER_INITIALIZATION_ATTEMPTS 1)
set(LAMBDA_PERTURBATION_STEP_SIZE 1)
set(LOG_OFFSET 1.0)

if (CMAKE_HOST_WIN32)
	set(NOMINMAX 1)
endif()

find_package(MKL QUIET)
find_package(GSL QUIET)
find_package(OpenBLAS QUIET)

CHECK_INCLUDE_FILE(getopt.h HAVE_GETOPT_H)

configure_file(config.h.in config.h)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(OPTIMIZER_STRATEGY "NelderMead" CACHE STRING "Optimizer strategy (LBFGS, Nelder-Mead, some others)")
set(DISCRETIZATION_RANGE "200" CACHE STRING "How finely to divide the continuous range")
set(MAX_STACK_FAMILY_SIZE "1000" CACHE STRING "Max family size that can be handled efficiently")

option(MODEL_GENE_EXPRESSION_LOGS "Perform calculations using log values rather than linear values" ON)
option(USE_MAX_PROBABILITY "Evaluate probability as maximum value rather than cumulative" OFF)
option(INCLUDE_GIT_IN_VERSION_INFO "Write information about the last GIT commit into version info" OFF)
option(MAX_STACK_FAMILY_SIZE "Max family size that can be handled efficiently" 1000)
 

if (MODEL_GENE_EXPRESSION_LOGS)
  add_definitions(-DMODEL_GENE_EXPRESSION_LOGS)
endif()

if (USE_MAX_PROBABILITY)
  add_definitions(-DUSE_MAX_PROBABILITY)
endif()

add_definitions(-DOPTIMIZER_STRATEGY=${OPTIMIZER_STRATEGY})
add_definitions(-DDISCRETIZATION_RANGE=${DISCRETIZATION_RANGE})
add_definitions(-DMAX_STACK_FAMILY_SIZE=${MAX_STACK_FAMILY_SIZE})


if (INCLUDE_GIT_IN_VERSION_INFO)
  # Define the two required variables before including
  # the source code for watching a git repository.
  set(PRE_CONFIGURE_FILE "git_version.cpp.in")
  set(POST_CONFIGURE_FILE "${CMAKE_CURRENT_BINARY_DIR}/git_version.cpp")
  include(git_watcher.cmake)

  # Create a library out of the compiled post-configure file.
  add_definitions(-DINCLUDE_GIT_IN_VERSION_INFO)
  add_library(git STATIC ${POST_CONFIGURE_FILE})
  target_include_directories(git PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
  add_dependencies(git check_git)
endif()

# add the executable
add_executable(cafe5 main.cpp)
add_executable(tests EXCLUDE_FROM_ALL test.cpp)

target_compile_definitions(cafe5 PRIVATE DOCTEST_CONFIG_DISABLE)
target_compile_definitions(tests PRIVATE SILENT ELPP_NO_CHECK_MACROS)
target_compile_options(tests PRIVATE -Wno-range-loop-construct)
add_subdirectory(src)

if (INCLUDE_GIT_IN_VERSION_INFO)
  target_link_libraries(cafe5 PRIVATE git)
  target_link_libraries(tests PRIVATE git)
endif()

find_package(ZLIB)
target_link_libraries(cafe5 PRIVATE ZLIB::ZLIB)
target_link_libraries(tests PRIVATE ZLIB::ZLIB)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(cafe5 PRIVATE OpenMP::OpenMP_CXX)
    target_link_libraries(tests PRIVATE OpenMP::OpenMP_CXX)
endif()

if(Boost_FOUND)
    target_link_libraries(cafe5 PRIVATE Boost::program_options)
    target_link_libraries(tests PRIVATE Boost::program_options)
endif()

if(MKL_FOUND)
    target_include_directories(cafe5 PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_INCLUDE_DIRECTORIES>)
    target_include_directories(tests PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_INCLUDE_DIRECTORIES>)

    target_compile_options(cafe5 PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_COMPILE_OPTIONS>)
    target_compile_options(tests PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_COMPILE_OPTIONS>)

    target_link_libraries(cafe5 PUBLIC $<LINK_ONLY:MKL::MKL>)
    target_link_libraries(tests PUBLIC $<LINK_ONLY:MKL::MKL>)
    add_compile_definitions(BLAS_AVAILABLE)
elseif(OpenBLAS_FOUND)
    target_include_directories(cafe5 PUBLIC ${OpenBLAS_INCLUDE_DIRS})
    target_include_directories(tests PUBLIC ${OpenBLAS_INCLUDE_DIRS})

    target_link_libraries(cafe5 PUBLIC ${OpenBLAS_LIBRARIES})
    target_link_libraries(tests PUBLIC ${OpenBLAS_LIBRARIES})
    add_compile_definitions(BLAS_AVAILABLE)
elseif (GSL_FOUND)
    target_include_directories(cafe5 PUBLIC ${GSL_INCLUDE_DIRS})
    target_include_directories(tests PUBLIC ${GSL_INCLUDE_DIRS})

    target_link_libraries(cafe5 PUBLIC ${GSL_LIBRARIES})
    target_link_libraries(tests PUBLIC ${GSL_LIBRARIES})
    add_compile_definitions(BLAS_AVAILABLE)
endif()

# GCC or Clang
add_definitions(-include config.h)
target_compile_options(tests PRIVATE -g -rdynamic -Wall -Werror)


install(TARGETS cafe5 DESTINATION bin)
install(DIRECTORY examples DESTINATION .)
install(DIRECTORY docs/manual DESTINATION .)
install(FILES README.md LICENSE DESTINATION .)

include(CPack)

