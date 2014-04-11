INCLUDE_DIRECTORIES(BEFORE ${Boost_INCLUDE_DIR})

set(SOURCES PGData.cpp
            StdCostFunc.cpp
            FixedContactConstr.cpp StaticStabilityConstr.cpp
            PositiveForceConstr.cpp FrictionConeConstr.cpp
            PlanarSurfaceConstr.cpp CollisionConstr.cpp
            PostureGenerator.cpp)
set(HEADERS PGData.h
            ConfigStruct.h
            StdCostFunc.h
            FixedContactConstr.h StaticStabilityConstr.h
            PositiveForceConstr.h FrictionConeConstr.h
            PlanarSurfaceConstr.h TorqueConstr.h
            CollisionConstr.h EllipseContactConstr.h
            IterationCallback.h
            PostureGenerator.h)

add_library(PG SHARED ${SOURCES} ${HEADERS})
PKG_CONFIG_USE_DEPENDENCY(PG scd)
PKG_CONFIG_USE_DEPENDENCY(PG SpaceVecAlg)
PKG_CONFIG_USE_DEPENDENCY(PG RBDyn)
PKG_CONFIG_USE_DEPENDENCY(PG roboptim-core)

set(INSTALL_PATH include/PG)
install(FILES ${HEADERS} DESTINATION ${INSTALL_PATH})
install(TARGETS PG DESTINATION "lib")

