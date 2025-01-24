# DO NOT EDIT THIS FILE! THIS FILE IS AUTOMATICALLY GENERATED BY xml2cmake
CMAKE_MINIMUM_REQUIRED(VERSION 3.10 FATAL_ERROR)
SET(VISIT_PLUGIN_DIR "/Users/benwibking/.visit/3.4.2/darwin-arm64/plugins")
INCLUDE("/Applications/VisIt.app/Contents/Resources/3.4.2/darwin-arm64/include/PluginVsInstall.cmake")
INCLUDE("VisItLibraryDependencies.cmake")

PROJECT(openpmd_database)

SET(COMMON_SOURCES
openpmdPluginInfo.C
openpmdCommonPluginInfo.C
)

SET(LIBI_SOURCES
openpmdPluginInfo.C
)

SET(LIBM_SOURCES
openpmdMDServerPluginInfo.C
${COMMON_SOURCES}
avtopenpmdFileFormat.C
dataLayoutTransform.C
)

SET(LIBE_SOURCES
openpmdEnginePluginInfo.C
${COMMON_SOURCES}
avtopenpmdFileFormat.C
dataLayoutTransform.C
)

INCLUDE_DIRECTORIES(
${CMAKE_CURRENT_SOURCE_DIR}
${VISIT_DATABASE_INCLUDES}
)

LINK_DIRECTORIES(${VISIT_LIBRARY_DIR})

# OpenPMD-api superbuild
set(openPMD_USE_PYTHON OFF CACHE BOOL "")
set(openPMD_BUILD_TESTING OFF CACHE BOOL "")
set(openPMD_BUILD_EXAMPLES OFF CACHE BOOL "")
set(openPMD_BUILD_CLI_TOOLS OFF CACHE BOOL "")
set(openPMD_INSTALL OFF CACHE BOOL "")
add_subdirectory(${openpmd_database_SOURCE_DIR}/extern/openPMD-api ${openpmd_database_BINARY_DIR}/openPMD-api)
include_directories(${OpenPMD_INCLUDE_DIRS_RET})
link_libraries(openPMD::openPMD)

# mdspan
add_subdirectory(${openpmd_database_SOURCE_DIR}/extern/mdspan ${openpmd_database_BINARY_DIR}/mdspan)
include_directories(${mdspan_INCLUDE_DIRS_RET})
link_libraries(std::mdspan)

ADD_LIBRARY(IopenpmdDatabase ${LIBI_SOURCES})
TARGET_LINK_LIBRARIES(IopenpmdDatabase visitcommon)
SET(INSTALLTARGETS IopenpmdDatabase)

IF(NOT VISIT_ENGINE_ONLY AND NOT VISIT_DBIO_ONLY)
    ADD_LIBRARY(MopenpmdDatabase ${LIBM_SOURCES})
    TARGET_LINK_LIBRARIES(MopenpmdDatabase visitcommon avtdbatts avtdatabase_ser)
    ADD_TARGET_DEFINITIONS(MopenpmdDatabase MDSERVER)
    SET(INSTALLTARGETS ${INSTALLTARGETS} MopenpmdDatabase)
ENDIF(NOT VISIT_ENGINE_ONLY AND NOT VISIT_DBIO_ONLY)

ADD_LIBRARY(EopenpmdDatabase_ser ${LIBE_SOURCES})
TARGET_LINK_LIBRARIES(EopenpmdDatabase_ser visitcommon avtpipeline_ser avtdatabase_ser)
ADD_TARGET_DEFINITIONS(EopenpmdDatabase_ser ENGINE)
SET(INSTALLTARGETS ${INSTALLTARGETS} EopenpmdDatabase_ser)

IF(VISIT_PARALLEL)
    ADD_PARALLEL_LIBRARY(EopenpmdDatabase_par ${LIBE_SOURCES})
    TARGET_LINK_LIBRARIES(EopenpmdDatabase_par visitcommon avtpipeline_par avtdatabase_par)
    ADD_TARGET_DEFINITIONS(EopenpmdDatabase_par ENGINE)
    SET(INSTALLTARGETS ${INSTALLTARGETS} EopenpmdDatabase_par)
ENDIF(VISIT_PARALLEL)

VISIT_INSTALL_DATABASE_PLUGINS(${INSTALLTARGETS})

