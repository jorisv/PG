INCLUDE_DIRECTORIES(BEFORE ${Boost_INCLUDE_DIR})

ENABLE_TESTING()

set(Boost_USE_STATIC_LIBS ON)
set(BOOST_COMPONENTS unit_test_framework)
SEARCH_FOR_BOOST()

include_directories("${PROJECT_SOURCE_DIR}/src")
include_directories(${Boost_INCLUDE_DIRS})

add_executable(PGTest PGTest.cpp)
target_link_libraries(PGTest ${Boost_LIBRARIES})
PKG_CONFIG_USE_DEPENDENCY(PGTest scd)
PKG_CONFIG_USE_DEPENDENCY(PGTest SpaceVecAlg)
PKG_CONFIG_USE_DEPENDENCY(PGTest RBDyn)
PKG_CONFIG_USE_DEPENDENCY(PGTest roboptim-core)
PKG_CONFIG_USE_DEPENDENCY(PGTest roboptim-core-plugin-ipopt)
add_test(PGTestUnit PGTest)

