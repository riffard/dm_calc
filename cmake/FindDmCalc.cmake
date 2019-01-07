#FindDmCalc.cmake
# 
# Finds the rapidjson library
#
# This will define the following variables
#
#    DmCalc_FOUND
#    DmCalc_INCLUDE_DIRS

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(DmCalc
    REQUIRED_VARS DmCalc_INCLUDE_DIR)


message("Cacauette !!!!!!!!!!")

if(DmCalc_FOUND)
    get_filename_component(DmCalc_INCLUDE_DIRS ${DmCalc_INCLUDE_DIR} DIRECTORY)
endif()

