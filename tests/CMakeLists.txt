INCLUDE_DIRECTORIES(BEFORE ${Boost_INCLUDE_DIR})

ENABLE_TESTING()

set(Boost_USE_STATIC_LIBS ON)
set(BOOST_COMPONENTS unit_test_framework)
SEARCH_FOR_BOOST()

include_directories("${PROJECT_SOURCE_DIR}/src")
include_directories(${Boost_INCLUDE_DIRS})

set(HEADERS Z12Arm.h XYZ12Arm.h)

macro(addUnitTest name)
  add_executable(${name} ${name}.cpp ${HEADERS})
  target_link_libraries(${name} ${Boost_LIBRARIES} PG)
  PKG_CONFIG_USE_DEPENDENCY(${name} sch-core)
  PKG_CONFIG_USE_DEPENDENCY(${name} SpaceVecAlg)
  PKG_CONFIG_USE_DEPENDENCY(${name} RBDyn)
  PKG_CONFIG_USE_DEPENDENCY(${name} roboptim-core)
  add_test(${name}Unit ${name})
endmacro(addUnitTest)

addUnitTest("PGTest")
addUnitTest("DiffTest")
