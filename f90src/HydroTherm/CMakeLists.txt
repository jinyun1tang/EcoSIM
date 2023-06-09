add_subdirectory(PhysData)
add_subdirectory(SnowPhys)
add_subdirectory(SoilPhys)
set(ECOSIM_LIBRARIES PhysData SnowPhys;SoilPhys;${ECOSIM_LIBRARIES} PARENT_SCOPE)
