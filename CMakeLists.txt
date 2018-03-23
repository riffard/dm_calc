###############  CMakeLists for MimacTools build   #####################
#  Written by Quentin Riffard (riffard@lpsc.in2p3.fr)                  #

#  To create an Xcode project:                                         #
#         cmake -G Xcode .                                             #

###############             cmake version          #####################
cmake_minimum_required(VERSION 2.6)

###############              Project name          #####################
project(DM_CALC)

###############        Root include/lib dir        #####################
find_package(ROOT REQUIRED COMPONENTS RooStats MathCore)

include_directories(${ROOT_INCLUDE_DIRS})

link_directories( ${ROOT_LIBRARY_DIR})

###############        Root include/lib dir        #####################
find_package(Boost)

include_directories(${BOOST_INCLUDE_DIRS})

link_directories( ${BOOST_LIBRARY_DIR})







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

file(GLOB sources_model_wimp ${PROJECT_SOURCE_DIR}/model_wimp/src/*.cc)
file(GLOB headers_model_wimp ${PROJECT_SOURCE_DIR}/model_wimp/include/*.hh)

file(GLOB sources_bkg ${PROJECT_SOURCE_DIR}/background/src/*.cc)
file(GLOB headers_bkg ${PROJECT_SOURCE_DIR}/background/include/*.hh)

file(GLOB sources_det ${PROJECT_SOURCE_DIR}/Detectors/src/*.cc)
file(GLOB headers_det ${PROJECT_SOURCE_DIR}/Detectors/include/*.hh)


set(all_header "${headers_common};${headers_model_wimp};${headers_bkg};${headers_det}")

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
set(LIBRARY_OUTPUT_PATH "${CMAKE_BINARY_DIR}/lib")
set(EXECUTABLE_OUTPUT_PATH "${CMAKE_BINARY_DIR}/bin")



###############      build the library     #####################
add_library(dm_calc SHARED ${sources_common} ${sources_model_wimp} ${sources_bkg} ${sources_det} ${headers_det})

target_link_libraries(dm_calc ${ROOT_LIBRARIES} ${BOOST_LIBRARIES} MathMore)

LINK_DIRECTORIES(${CMAKE_BINARY_DIR}/lib)

############### install the library ###################
install(TARGETS dm_calc DESTINATION ${PROJECT_SOURCE_DIR}/lib)