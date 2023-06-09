add_subdirectory(PhysData)
add_subdirectory(SnowPhys)
<<<<<<< HEAD
add_subdirectory(SurfPhys)
add_subdirectory(SoilPhys)
set(ECOSIM_LIBRARIES ${ECOSIM_LIBRARIES} PARENT_SCOPE)
add_subdirectory(SurfPhys)
add_subdirectory(SoilPhys)
set(ECOSIM_LIBRARIES PhysData SnowPhys;SurfPhys;SoilPhys;${ECOSIM_LIBRARIES} PARENT_SCOPE)
=======
add_subdirectory(SoilPhys)
set(ECOSIM_LIBRARIES PhysData SnowPhys;SoilPhys;${ECOSIM_LIBRARIES} PARENT_SCOPE)
>>>>>>> 0ab9258 (separated the snow hydrothermal code from soil)
