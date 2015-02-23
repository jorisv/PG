include(../../cmake/python.cmake)

# Look for Python 2.7
SET(Python_ADDITIONAL_VERSIONS 2.7)
FINDPYTHON()

# PG

# create the package in build dir for testing purpose
set(LIBRARY_OUTPUT_PATH ${CMAKE_CURRENT_BINARY_DIR}/pg)
configure_file(__init__.py ${CMAKE_CURRENT_BINARY_DIR}/pg/__init__.py COPYONLY)

set(OUTPUT_BINDING ${CMAKE_CURRENT_BINARY_DIR}/pg.cpp)

# generate python binding code
add_custom_command (
  OUTPUT ${OUTPUT_BINDING}
  COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/generate.py ${OUTPUT_BINDING}
  DEPENDS generate.py
)

# build the library
set(SOURCES ${OUTPUT_BINDING})
include_directories(.)
include_directories(../../src)
include_directories(SYSTEM ${PYTHON_INCLUDE_DIRS})

add_library(_pg SHARED ${SOURCES})
target_link_libraries(_pg PG)
PKG_CONFIG_USE_DEPENDENCY(_pg sch-core)
PKG_CONFIG_USE_DEPENDENCY(_pg SpaceVecAlg)
PKG_CONFIG_USE_DEPENDENCY(_pg RBDyn)
PKG_CONFIG_USE_DEPENDENCY(_pg roboptim-core)

set_target_properties(_pg PROPERTIES PREFIX "")
set_target_properties(_pg PROPERTIES INSTALL_RPATH
    "${CMAKE_INSTALL_RPATH};${ROBOPTIM_CORE_PLUGIN_IPOPT_LIBRARY_DIRS}")

# install rules
install(TARGETS _pg DESTINATION "${PYTHON_SITELIB}/pg")
PYTHON_INSTALL_BUILD(pg __init__.py "${PYTHON_SITELIB}")
