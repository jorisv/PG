set(SOURCES lib.cpp)
set(HEADERS FK.h PGData.h EigenAutoDiffScalar.h
            StdCostFunc.h
            FixedContactConstr.h
            PostureGenerator.h)

add_library(PG SHARED ${SOURCES} ${HEADERS})


PKG_CONFIG_USE_DEPENDENCY(PG SpaceVecAlg)
PKG_CONFIG_USE_DEPENDENCY(PG RBDyn)
PKG_CONFIG_USE_DEPENDENCY(PG roboptim-core)
PKG_CONFIG_USE_DEPENDENCY(PG roboptim-core-plugin-ipopt)


set(INSTALL_PATH include/PG)
install(FILES ${HEADERS} DESTINATION ${INSTALL_PATH})
install(TARGETS PG DESTINATION "lib")

