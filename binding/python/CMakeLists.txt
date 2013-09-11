find_package(PythonInterp REQUIRED)
find_package(PythonLibs REQUIRED)

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
include_directories(${PYTHON_INCLUDE_DIRS})

add_library(_pg SHARED ${SOURCES})
PKG_CONFIG_USE_DEPENDENCY(_pg scd)
PKG_CONFIG_USE_DEPENDENCY(_pg SpaceVecAlg)
PKG_CONFIG_USE_DEPENDENCY(_pg RBDyn)
PKG_CONFIG_USE_DEPENDENCY(_pg roboptim-core)
PKG_CONFIG_USE_DEPENDENCY(_pg roboptim-core-plugin-ipopt)

set_target_properties(_pg PROPERTIES PREFIX "")
set_target_properties(_pg PROPERTIES INSTALL_RPATH
    "${CMAKE_INSTALL_RPATH};${ROBOPTIM_CORE_PLUGIN_IPOPT_LIBRARY_DIRS}")


# install rules
# find python version (not implemented in cmake 2.8.5)
execute_process(COMMAND "${PYTHON_EXECUTABLE}" -V ERROR_VARIABLE _VERSION
                OUTPUT_QUIET ERROR_STRIP_TRAILING_WHITESPACE)
string(REPLACE "Python " "" PYTHON_VERSION_STRING "${_VERSION}")
string(REGEX REPLACE "^([0-9]+)\\.[0-9]+\\.[0-9]+.*" "\\1"
       PYTHON_VERSION_MAJOR "${PYTHON_VERSION_STRING}")
string(REGEX REPLACE "^[0-9]+\\.([0-9])+\\.[0-9]+.*" "\\1"
       PYTHON_VERSION_MINOR "${PYTHON_VERSION_STRING}")
set(INSTALL_PATH
  "lib/python${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}/dist-packages/pg/")
install(TARGETS _pg DESTINATION ${INSTALL_PATH})
install(FILES __init__.py DESTINATION ${INSTALL_PATH})
