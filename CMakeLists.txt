###############    CMakeLists for DM_CALC build    #####################
#  Written by Quentin Riffard (qriffard@lbl.gov)                       #
#  To create an Xcode project:                                         #
#         cmake -G Xcode .                                             #
########################################################################

###############             cmake version          #####################
cmake_minimum_required(VERSION 3.2.3)

###############    project name and version        #####################
project(DmCalc VERSION 1.0.0 LANGUAGES CXX)

###############            Build type              #####################
set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_BUILD_TYPE RelWithDebInfo)


# If the install directory not specified, use the build directory
if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
   MESSAGE(STATUS "Using default installation directory in build directory")
 set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}")
endif ()

################################################################################
################################################################################


################################################################################
# Find all the packages
################################################################################
################################################################################
# ROOT
find_package(ROOT COMPONENTS MathMore Matrix Tree)
IF(NOT ${ROOT_FOUND})
	message("Root not found... Let's try an other way ;)")
	include($ENV{ROOTSYS}/etc/root/cmake/FindROOT.cmake)
	include_directories(${ROOT_INCLUDE_DIRS})
	link_directories( ${ROOT_LIBRARY_DIR})
ELSE(NOT ${ROOT_FOUND})
	 add_definitions(${ROOT_CXX_FLAGS})
	 include(${ROOT_USE_FILE})
ENDIF()


###############             Project files           #####################

file(GLOB sources_Common ${PROJECT_SOURCE_DIR}/Common/*.cc)
file(GLOB headers_Common ${PROJECT_SOURCE_DIR}/Common/*.hh)

file(GLOB sources_models ${PROJECT_SOURCE_DIR}/Models/*/*.cc)
file(GLOB headers_models ${PROJECT_SOURCE_DIR}/Models/*/*.hh)

file(GLOB sources_det ${PROJECT_SOURCE_DIR}/Detectors/*.cc)
file(GLOB headers_det ${PROJECT_SOURCE_DIR}/Detectors/*.hh)



###############      build the library     #####################
add_library(DmCalc SHARED ${sources_Common} ${sources_models} ${sources_bkg} ${sources_det} ${headers_det})

target_link_libraries(DmCalc ${ROOT_LIBRARIES})

target_include_directories( DmCalc
	          PUBLIC
                        $<INSTALL_INTERFACE:include>
                  PRIVATE
		        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Detectors>
		        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Models/Backgrounds>
		        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Models/Signals>
		        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Common>)

##############################################
# Installation instructions

# Destinamtion dir
set(INSTALL_BIN_DIR "bin")
set(INSTALL_INCL_DIR "include")
set(INSTALL_LIB_DIR "lib")
set(INSTALL_CMAKE_DIR "cmake")


# Target properties
set_target_properties(DmCalc PROPERTIES PUBLIC_HEADER "${headers_Common};${headers_models};${headers_bkg};${headers_det}")

set_target_properties(DmCalc PROPERTIES
			     ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
			     LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")

install(TARGETS DmCalc
	EXPORT DmCalc-targets
	LIBRARY DESTINATION "${INSTALL_LIB_DIR}"
	PUBLIC_HEADER DESTINATION "${INSTALL_INCL_DIR}")


#Export the targets to a script
install(EXPORT DmCalc-targets
 	FILE DmCalcTargets.cmake
	DESTINATION ${INSTALL_CMAKE_DIR}
)


include(CMakePackageConfigHelpers)

configure_package_config_file(${CMAKE_CURRENT_SOURCE_DIR}/cmake/DmCalcConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/DmCalcConfig.cmake
    INSTALL_DESTINATION ${INSTALL_CMAKE_DIR}
)


install(FILES
    ${CMAKE_CURRENT_LIST_DIR}/cmake/FindDmCalc.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/DmCalcConfig.cmake
    DESTINATION ${INSTALL_CMAKE_DIR}
)

install(DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/DataBase
	DESTINATION .)


export(PACKAGE DmCalc)

############### Build the example ###################
file( GLOB  APP_SOURCES RELATIVE "${CMAKE_BINARY_DIR}" "Examples/*.cc" )


foreach( testsourcefile ${APP_SOURCES} )
    get_filename_component(testsourcefilename ${testsourcefile} NAME)
    string( REPLACE ".cc" ".bin" testname ${testsourcefilename} )

    add_executable(${testname} ${testsourcefile})
    target_include_directories(${testname}
                  PUBLIC
                        $<INSTALL_INTERFACE:include>
                  PRIVATE
		        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Models/Backgrounds>
		        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Models/Signals>
		        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Common>
                )

    target_link_libraries(${testname}
                PUBLIC
                        ${ROOT_LIBRARIES}
                PRIVATE
			DmCalc )
     set_target_properties(${testname} PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

     install(TARGETS ${testname} DESTINATION ${INSTALL_BIN_DIR})


endforeach( testsourcefile ${APP_SOURCES} )

