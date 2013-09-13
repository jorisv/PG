set(SOURCES lib.cpp)
set(HEADERS FK.h ID.h
            PGData.h EigenAutoDiffScalar.h AutoDiffFunction.h
            ConfigStruct.h
            StdCostFunc.h
            FixedContactConstr.h StaticStabilityConstr.h
            PositiveForceConstr.h FrictionConeConstr.h
            PlanarSurfaceConstr.h TorqueConstr.h
            CollisionConstr.h
            IpoptIntermediateCallback.h
            PostureGenerator.h)

# !!! PG is an empty library !!!
# !!! Its aim is to test header file and have an IDE integration !!!
# !!! Don't link against it !!!
add_library(PG SHARED ${SOURCES} ${HEADERS})


PKG_CONFIG_USE_DEPENDENCY(PG scd)
PKG_CONFIG_USE_DEPENDENCY(PG SpaceVecAlg)
PKG_CONFIG_USE_DEPENDENCY(PG RBDyn)
PKG_CONFIG_USE_DEPENDENCY(PG roboptim-core)
PKG_CONFIG_USE_DEPENDENCY(PG roboptim-core-plugin-ipopt)


set(INSTALL_PATH include/PG)
install(FILES ${HEADERS} DESTINATION ${INSTALL_PATH})

