ENABLE_TESTING()

set(Boost_USE_STATIC_LIBS ON)
find_package(Boost COMPONENTS unit_test_framework REQUIRED)

include_directories("${PROJECT_SOURCE_DIR}/src")
include_directories(${Boost_INCLUDE_DIRS})

#add_executable(QPSolverTest QPSolverTest.cpp)
#target_link_libraries(QPSolverTest ${Boost_LIBRARIES} PG)
#PKG_CONFIG_USE_DEPENDENCY(QPSolverTest SpaceVecAlg)
#PKG_CONFIG_USE_DEPENDENCY(QPSolverTest RBDyn)
#add_test(QPSolverTestUnit QPSolverTest)

