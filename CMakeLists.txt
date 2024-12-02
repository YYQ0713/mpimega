cmake_minimum_required(VERSION 2.8)
project(mpimeta)
set(CMAKE_VERBOSE_MAKEFILE ON)


set(CMAKE_C_COMPILER mpicc)
set(CMAKE_CXX_COMPILER mpicxx)


if (CMAKE_VERSION VERSION_LESS "3.1")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
else ()
    set(CMAKE_CXX_STANDARD 11)
endif ()

option(COVERAGE "Generate coverage report" OFF)
option(STATIC_BUILD "Build static executation" OFF)
option(SANITIZER "Enable sanitizers" OFF)
option(TSAN "Enable thread sanitizers" OFF)

include_directories(src)
FILE(GLOB_RECURSE SDBG_SOURCE "src/sdbg/*.cpp")
FILE(GLOB_RECURSE CX1_SOURCE "src/sorting/*.cpp")
FILE(GLOB_RECURSE SEQ_SOURCE "src/sequence/*.cpp")
FILE(GLOB_RECURSE MPI_SOURCE "src/mpienv/*.cpp")

LIST(APPEND OTHER_SOURCE
        src/main.cpp
        src/main_buildlib.cpp
        src/main_sdbg_build.cpp
        src/utils/options_description.cpp
        )


if (STATIC_BUILD)
    set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
endif (STATIC_BUILD)

find_package(ZLIB REQUIRED)
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

add_executable(mpimeta_core ${OTHER_SOURCE} ${SDBG_SOURCE} ${SEQ_SOURCE} ${CX1_SOURCE} ${MPI_SOURCE})
add_executable(mpimeta_core_popcnt ${OTHER_SOURCE} ${SDBG_SOURCE} ${SEQ_SOURCE} ${CX1_SOURCE} ${MPI_SOURCE})
add_executable(mpimeta_core_no_hw_accel ${OTHER_SOURCE} ${SDBG_SOURCE} ${SEQ_SOURCE} ${CX1_SOURCE} ${MPI_SOURCE})

set_target_properties(mpimeta_core PROPERTIES COMPILE_FLAGS "-mbmi2 -DUSE_BMI2 -mpopcnt")
set_target_properties(mpimeta_core_popcnt PROPERTIES COMPILE_FLAGS "-mpopcnt")


target_link_libraries(mpimeta_core ${ZLIB_LIBRARIES})
target_link_libraries(mpimeta_core_popcnt ${ZLIB_LIBRARIES})
target_link_libraries(mpimeta_core_no_hw_accel ${ZLIB_LIBRARIES})


add_custom_target(
        mpimeta ALL
        COMMAND cp ${CMAKE_SOURCE_DIR}/src/mpimeta .
)

add_dependencies(mpimeta mpimeta_core mpimeta_core_popcnt mpimeta_core_no_hw_accel)