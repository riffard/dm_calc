###############  CMakeLists for MimacTools build   #####################
#  Written by Quentin Riffard (riffard@lpsc.in2p3.fr)                  #

#  To create an Xcode project:                                         #
#         cmake -G Xcode .                                             #

###############             cmake version          #####################
cmake_minimum_required(VERSION 2.6)

###############              Project name          #####################
project(dm_calc)

###############        Root include/lib dir        #####################
find_package(ROOT COMPONENTS MathCore)
IF(NOT ${ROOT_FOUND})
	message("Root not found... Let's try an other way ;)")
	include($ENV{ROOTSYS}/etc/root/cmake/FindROOT.cmake)
ENDIF()

include_directories(${ROOT_INCLUDE_DIRS})

link_directories( ${ROOT_LIBRARY_DIR})

###############        Root include/lib dir        #####################

message("$ENV{Boost_INCLUDE_DIR}")

find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})
link_directories( ${Boost_LIBRARY_DIR})


###############                #####################
if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
   MESSAGE(STATUS "Using default installation directory in build directory")
 set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}")
endif ()

set(INSTALL_LIB_DIR "lib" CACHE PATH "Installation directory for libraries")
set(INSTALL_INCLUDE_DIR "include" CACHE PATH
  "Installation directory for header files")
set(INSTALL_CMAKE_DIR "${DEF_INSTALL_CMAKE_DIR}" CACHE PATH
  "Installation directory for CMake files")

foreach(p LIB INCLUDE CMAKE)
  set(var INSTALL_${p}_DIR)
  if(NOT IS_ABSOLUTE "${${var}}")
    set(${var} "${CMAKE_INSTALL_PREFIX}/${${var}}")
  endif()
endforeach()




###############            Minimal flags           #####################
IF(${UNIX})
  set (CMAKE_CXX_FLAGS "-std=c++11 -g3 -Wall -O0")
ENDIF(${UNIX})

IF(${APPLE})
  set (CMAKE_CXX_FLAGS "-std=c++0x -g3 -Wall -fPIC -fno-common -O0")
ENDIF(${APPLE})

if(NOT DEFINED CMAKE_MACOSX_RPATH)
  set(CMAKE_MACOSX_RPATH 0)
endif()

#######   Check the compiler and set the compile and link flags  #######
set(CMAKE_BUILD_TYPE Debug)




###############             Root flags            #####################
if (CMAKE_SYSTEM_NAME MATCHES Darwin)
   EXEC_PROGRAM("sw_vers -productVersion | cut -d . -f 1-2" OUTPUT_VARIABLE MAC_OS_VERSION)
   MESSAGE("--- Found a Mac OS X System ${MAC_OS_VERSION}")
   if (CMAKE_COMPILER_IS_GNUCXX)
      MESSAGE("--- Found GNU compiler collection")

      STRING(COMPARE EQUAL "10.5" "${MAC_OS_VERSION}" MAC_OS_10_5)
      IF(MAC_OS_10_5)
        SET(CMAKE_CXX_FLAGS "-m64")
        SET(CMAKE_Fortran_FLAGS "-m64")
      ENDIF(MAC_OS_10_5)

      SET(CMAKE_SHARED_LIBRARY_CREATE_C_FLAGS "${CMAKE_SHARED_LIBRARY_CREATE_C_FLAGS} -flat_namespace -single_module -undefined dynamic_lookup")
      SET(CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS "${CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS} -flat_namespace -single_module -undefined dynamic_lookup")

      # Select flags.
      set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g")
      set(CMAKE_CXX_FLAGS_RELEASE        "-O2")
      set(CMAKE_CXX_FLAGS_DEBUG          "-g -O2 -fno-reorder-blocks -fno-schedule-insns -fno-inline")
      set(CMAKE_CXX_FLAGS_DEBUGFULL      "-g3 -fno-inline -Wnon-virtual-dtor -Wno-long-long -ansi -Wundef -Wcast-align -Wchar-subscripts -Wall -W -Wpointer-arith -Wformat-security -fno-exceptions -fno-check-new -fno-common")
      set(CMAKE_CXX_FLAGS_PROFILE        "-g3 -fno-inline -ftest-coverage -fprofile-arcs")
      set(CMAKE_C_FLAGS_RELWITHDEBINFO   "-O2 -g")
      set(CMAKE_C_FLAGS_RELEASE          "-O2")
      set(CMAKE_C_FLAGS_DEBUG            "-g -O2 -fno-reorder-blocks -fno-schedule-insns -fno-inline")
      set(CMAKE_C_FLAGS_DEBUGFULL        "-g3 -fno-inline -Wno-long-long -std=iso9899:1990 -Wundef -Wcast-align -Werror-implicit-function-declaration -Wchar-subscripts -Wall -W -Wpointer-arith -Wwrite-strings -Wformat-security -Wmissing-format-attribute -fno-common")
      set(CMAKE_C_FLAGS_PROFILE          "-g3 -fno-inline -ftest-coverage -fprofile-arcs")

      endif (CMAKE_COMPILER_IS_GNUCXX)

endif (CMAKE_SYSTEM_NAME MATCHES Darwin)


###############             Project files           #####################

file(GLOB sources_common ${PROJECT_SOURCE_DIR}/common/src/*.cc)
file(GLOB headers_common ${PROJECT_SOURCE_DIR}/common/include/*.hh)

file(GLOB sources_models ${PROJECT_SOURCE_DIR}/Models/*/*.cc)
file(GLOB headers_models ${PROJECT_SOURCE_DIR}/Models/*/*.hh)

file(GLOB sources_det ${PROJECT_SOURCE_DIR}/Detectors/src/*.cc)
file(GLOB headers_det ${PROJECT_SOURCE_DIR}/Detectors/include/*.hh)




set(all_header "${headers_common};${headers_models};${headers_bkg};${headers_det}")

foreach (header ${all_header})
	GET_FILENAME_COMPONENT(dir_path ${header} PATH)
	include_directories(${dir_path})
endforeach(header ${all_header})


###############      library extension  #####################
IF(${UNIX})
  SET(LIBRARY_PROPERTIES ${LIBRARY_PROPERTIES} SUFFIX ".so" )
ENDIF(${UNIX})

IF(${APPLE})
  SET(LIBRARY_PROPERTIES ${LIBRARY_PROPERTIES} VERSION "${LIB_VERSION}" SOVERSION "${MAJOR_VERSION}" SUFFIX ".dylib"  )
ENDIF(${APPLE})



###############          Output directory          #####################
#set(LIBRARY_OUTPUT_PATH "${CMAKE_BINARY_DIR}/lib")
#set(EXECUTABLE_OUTPUT_PATH "${CMAKE_BINARY_DIR}/bin")



###############      build the library     #####################
add_library(dm_calc SHARED ${sources_common} ${sources_models} ${sources_bkg} ${sources_det} ${headers_det})

target_link_libraries(dm_calc ${ROOT_LIBRARIES} ${BOOST_LIBRARIES} MathMore)

LINK_DIRECTORIES(${CMAKE_BINARY_DIR}/lib)

set_target_properties(dm_calc PROPERTIES PUBLIC_HEADER "${headers_common};${headers_models};${headers_bkg};${headers_det}")


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
		        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/common>
                )



    target_link_libraries(${testname}
                PUBLIC
                        ${ROOT_LIBRARIES}
                PRIVATE
		dm_calc )

endforeach( testsourcefile ${APP_SOURCES} )


############### install the library ###################
#install(TARGETS dm_calc DESTINATION ${PROJECT_SOURCE_DIR}/lib)
install(TARGETS dm_calc
  # IMPORTANT: Add the foo library to the "export-set"
  EXPORT dm_calcTargets
  LIBRARY DESTINATION "${INSTALL_LIB_DIR}"
  PUBLIC_HEADER DESTINATION "${INSTALL_INCLUDE_DIR}")

export(PACKAGE dm_calc)
export(TARGETS dm_calc
  FILE "${PROJECT_BINARY_DIR}/dm_calcTargets.cmake")



file(RELATIVE_PATH REL_INCLUDE_DIR "${INSTALL_CMAKE_DIR}"
   "${INSTALL_INCLUDE_DIR}")
# ... for the build tree
set(CONF_INCLUDE_DIRS "${PROJECT_SOURCE_DIR}" "${PROJECT_BINARY_DIR}")
configure_file(dm_calcConfig.cmake.in
  "${PROJECT_BINARY_DIR}/dm_calcConfig.cmake" @ONLY)
# ... for the install tree
set(CONF_INCLUDE_DIRS "\${dm_calc_CMAKE_DIR}/${REL_INCLUDE_DIR}")
configure_file(dm_calcConfig.cmake.in
  "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/dm_calcConfig.cmake" @ONLY)
# ... for both
configure_file(dm_calcConfigVersion.cmake.in
  "${PROJECT_BINARY_DIR}/dm_calcConfigVersion.cmake" @ONLY)
 
# Install the dm_calcConfig.cmake and dm_calcConfigVersion.cmake
install(FILES
  "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/dm_calcConfig.cmake"
  "${PROJECT_BINARY_DIR}/dm_calcConfigVersion.cmake"
  DESTINATION "${INSTALL_CMAKE_DIR}" COMPONENT dev)
 
# Install the export set for use with the install-tree
install(EXPORT dm_calcTargets DESTINATION
  "${INSTALL_CMAKE_DIR}" COMPONENT dev)

