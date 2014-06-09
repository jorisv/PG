cmake_minimum_required(VERSION 2.8)

INCLUDE(cmake/base.cmake)
INCLUDE(cmake/boost.cmake)
INCLUDE(cmake/cpack.cmake)
INCLUDE(cmake/eigen.cmake)

SET(PROJECT_NAME PG)
SET(PROJECT_DESCRIPTION "...")
SET(PROJECT_URL "https://github.com/jorisv/PG")

# Disable -Werror on Unix for now.
SET(CXX_DISABLE_WERROR True)

SETUP_PROJECT()

option(PYTHON_BINDING "Generate python binding." ON)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -std=c++0x -pedantic")

#########################
# External dependencies
SEARCH_FOR_EIGEN()
SET(BOOST_REQUIRED 1.48)
SEARCH_FOR_BOOST()

ADD_REQUIRED_DEPENDENCY("sch-core")
ADD_REQUIRED_DEPENDENCY("SpaceVecAlg")
ADD_REQUIRED_DEPENDENCY("RBDyn")
ADD_REQUIRED_DEPENDENCY("roboptim-core > 2.0")

#########################

add_definitions(-DROBOPTIM_DO_NOT_CHECK_ALLOCATION)


add_subdirectory(src)
add_subdirectory(tests)

if(${PYTHON_BINDING})
  add_subdirectory(binding/python)
endif()

# Add dependency towards the library in the pkg-config file.
PKG_CONFIG_APPEND_LIBS(PG)

SETUP_PROJECT_FINALIZE()
SETUP_PROJECT_CPACK()

