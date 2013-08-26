set(SOURCES lib.cpp)
set(HEADERS FK.h ID.h
            PGData.h EigenAutoDiffScalar.h AutoDiffFunction.h
            StdCostFunc.h
            FixedContactConstr.h StaticStabilityConstr.h
            PositiveForceConstr.h FrictionConeConstr.h
            PostureGenerator.h)

add_library(PG SHARED ${SOURCES} ${HEADERS})


PKG_CONFIG_USE_DEPENDENCY(PG SpaceVecAlg)
PKG_CONFIG_USE_DEPENDENCY(PG RBDyn)
PKG_CONFIG_USE_DEPENDENCY(PG roboptim-core)
PKG_CONFIG_USE_DEPENDENCY(PG roboptim-core-plugin-ipopt)


set(INSTALL_PATH include/PG)
install(FILES ${HEADERS} DESTINATION ${INSTALL_PATH})
install(TARGETS PG DESTINATION "lib")

