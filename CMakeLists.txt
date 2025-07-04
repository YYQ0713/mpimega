cmake_minimum_required(VERSION 3.18)
project(mpimeta)
set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_CXX_STANDARD 17)

add_subdirectory(thirdparty/ga-5.8.2)
include_directories(thirdparty/ga-5.8.2/global/src)
include_directories(thirdparty/ga-5.8.2/ma)
include_directories(${CMAKE_BINARY_DIR}/thirdparty/ga-5.8.2/gaf2c)

add_subdirectory(thirdparty/fmt-11.2.0)

set(TBB_TEST OFF CACHE BOOL "Disable TBB tests")
set(TBB_EXAMPLES OFF CACHE BOOL "Disable TBB examples")
add_subdirectory(thirdparty/oneTBB-2022.2.0)

#if (CMAKE_VERSION VERSION_LESS "3.1")
#    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
#else ()
#endif ()

option(ZLIB_ENABLE_TESTS "" OFF)
option(ZLIB_COMPAT "" ON)

add_subdirectory(thirdparty/zlib-ng-2.2.4)

set(CMAKE_C_COMPILER mpicc)
set(CMAKE_CXX_COMPILER mpicxx)

option(COVERAGE "Generate coverage report" OFF)
option(STATIC_BUILD "Build static executation" OFF)
option(SANITIZER "Enable sanitizers" OFF)
option(TSAN "Enable thread sanitizers" OFF)

include_directories(src)
FILE(GLOB_RECURSE ASMBL_SOURCE "src/assembly/*.cpp")
FILE(GLOB_RECURSE SDBG_SOURCE "src/sdbg/*.cpp")
FILE(GLOB_RECURSE CX1_SOURCE "src/sorting/*.cpp")
FILE(GLOB_RECURSE SEQ_SOURCE "src/sequence/*.cpp")
FILE(GLOB_RECURSE MPI_SOURCE "src/mpienv/*.cpp")
FILE(GLOB_RECURSE LCASM_SOURCE "src/localasm/*.cpp")
FILE(GLOB_RECURSE IDBA_SOURCE "src/idba/*.cpp")

LIST(APPEND OTHER_SOURCE
		src/main_assemble.cpp
        src/main.cpp
        src/main_buildlib.cpp
        src/main_sdbg_build.cpp
		src/main_local_assemble.cpp
        src/utils/options_description.cpp
        )

if (STATIC_BUILD)
    set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
endif (STATIC_BUILD)

#find_package(ZLIB REQUIRED)
find_package(OpenMP REQUIRED)


set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DXXH_INLINE_ALL -ftemplate-depth=3000")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-unused-function")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprefetch-loop-arrays -funroll-loops")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D__XROOT__='\"${CMAKE_SOURCE_DIR}/src\"'")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D__XFILE__='\"$(subst ${CMAKE_SOURCE_DIR}/,,$(abspath $<))\"'")

set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${ZLIB_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")

set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
set(CMAKE_CXX_FLAGS_DEBUG "-g -ggdb -O1 -D_LIBCPP_DEBUG -D_GLIBCXX_DEBUG")


if (COVERAGE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O0 --coverage")
    set(COV_PY "coverage run")
endif (COVERAGE)

if (SANITIZER)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address -fsanitize=leak -fsanitize=undefined")
endif (SANITIZER)

if (TSAN)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=thread")
endif (TSAN)

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}: ${CMAKE_CXX_FLAGS}")

add_executable(mpimeta_core ${OTHER_SOURCE} ${SDBG_SOURCE} ${SEQ_SOURCE} ${CX1_SOURCE} ${MPI_SOURCE} ${ASMBL_SOURCE} ${LCASM_SOURCE} ${IDBA_SOURCE})
add_executable(mpimeta_core_popcnt ${OTHER_SOURCE} ${SDBG_SOURCE} ${SEQ_SOURCE} ${CX1_SOURCE} ${MPI_SOURCE} ${ASMBL_SOURCE} ${LCASM_SOURCE} ${IDBA_SOURCE})
add_executable(mpimeta_core_no_hw_accel ${OTHER_SOURCE} ${SDBG_SOURCE} ${SEQ_SOURCE} ${CX1_SOURCE} ${MPI_SOURCE} ${ASMBL_SOURCE} ${LCASM_SOURCE} ${IDBA_SOURCE})

set_target_properties(mpimeta_core PROPERTIES COMPILE_FLAGS "-mbmi2 -DUSE_BMI2 -mpopcnt")
set_target_properties(mpimeta_core_popcnt PROPERTIES COMPILE_FLAGS "-mpopcnt")

#target_link_libraries(mpimeta_core ${ZLIB_LIBRARIES})
#target_link_libraries(mpimeta_core_popcnt ${ZLIB_LIBRARIES})
#target_link_libraries(mpimeta_core_no_hw_accel ${ZLIB_LIBRARIES})

target_link_libraries(mpimeta_core 
		ga
		ga++
		m
		armci
		comex
		fmt::fmt
		zlib
		TBB::tbb)
target_link_libraries(mpimeta_core_popcnt 
		ga
		ga++
		m
		armci
		comex
		fmt::fmt
		zlib
		TBB::tbb)
target_link_libraries(mpimeta_core_no_hw_accel 
		ga
		ga++
		armci
		m
		comex
		fmt::fmt
		zlib
		TBB::tbb)

target_include_directories(mpimeta_core PRIVATE thirdparty/zlib-ng-2.2.4)
target_include_directories(mpimeta_core_popcnt PRIVATE thirdparty/zlib-ng-2.2.4)
target_include_directories(mpimeta_core_no_hw_accel PRIVATE thirdparty/zlib-ng-2.2.4)
target_include_directories(mpimeta_core PRIVATE thirdparty/oneTBB-2022.2.0/include)
target_include_directories(mpimeta_core_popcnt PRIVATE thirdparty/oneTBB-2022.2.0/include)
target_include_directories(mpimeta_core_no_hw_accel PRIVATE thirdparty/oneTBB-2022.2.0/include)
add_custom_target(
        mpimeta ALL
        COMMAND cp ${CMAKE_SOURCE_DIR}/src/mpimeta .
)

add_dependencies(mpimeta mpimeta_core mpimeta_core_popcnt mpimeta_core_no_hw_accel)