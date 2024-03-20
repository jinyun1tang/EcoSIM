add_subdirectory(utils)
add_subdirectory(minimath)
add_subdirectory(modelconfig)
add_subdirectory(mesh)
add_subdirectory(modelpars)
add_subdirectory(ecosim_datatype)
add_subdirectory(APIData)
add_subdirectory(disturbances)
add_subdirectory(microbial_bgc)
add_subdirectory(plant_bgc)
add_subdirectory(geochem)
add_subdirectory(HydroTherm)
add_subdirectory(APIs)
add_subdirectory(ioutils)
add_subdirectory(ecosim_mods)
add_subdirectory(balances)
add_subdirectory(ATSUtils)
add_subdirectory(modelforc)
add_subdirectory(transport)
add_subdirectory(main)

#message("ECOSIM_LIBRARIES ${ECOSIM_LIBRARIES}")
# Collect all the intermediate libraries
#set(LIBRARY_LIST subdir1_lib subdir2_lib subdir3_lib)

set(ECOSIM_LIBRARIES ${ECOSIM_LIBRARIES} PARENT_SCOPE)

# Combine the intermediate libraries into a single dynamic library
#add_library(ecosimff SHARED)
