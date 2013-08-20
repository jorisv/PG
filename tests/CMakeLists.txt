ENABLE_TESTING()

set(Boost_USE_STATIC_LIBS ON)
find_package(Boost COMPONENTS unit_test_framework REQUIRED)

include_directories("${PROJECT_SOURCE_DIR}/src")
include_directories(${Boost_INCLUDE_DIRS})

add_executable(AlgoTest AlgoTest.cpp)
target_link_libraries(AlgoTest ${Boost_LIBRARIES} PG)
PKG_CONFIG_USE_DEPENDENCY(AlgoTest SpaceVecAlg)
PKG_CONFIG_USE_DEPENDENCY(AlgoTest RBDyn)
add_test(AlgosTestUnit AlgoTest)

