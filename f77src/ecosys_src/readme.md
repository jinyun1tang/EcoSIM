# ecosys_src

|File                 | Description                                                                         |
|---------------------|-------------------------------------------------------------------------------------|
|day.f                |Reinitialize daily variables                                                         |
|erosion.f            |Soil erosion model                                                                   |
|exec.f               |Balance check of C, N, P, water and heat                                             |
|extract.f            |Aggregate all exchange fluxes and send them to redist                                |
|grosub.f             |Calculate all plant transformations                                                  |
|hfunc.f              |Plant phenology calculation                                                          |
|hour1.f              |Reinitialize hourly variables                                                        |
|nitro.f              |Microbial C, N, and P transformations                                                |
|redist.f             |Update soil C, N, P, water, heat and solute fluxes                                   |
|solute.f             |Solute model                                                                         |
|starte.f             |Initialize soil chemistry variables                                                  |
|startq.f             |Initialize plant variables                                                           |
|starts.f             |Initialize soil variables                                                            |
|stomate.f            |Stomatal conductance at maximum turgor                                               |
|trnsfr.f             |3D fluxes of non-salt solutes and gases                                              |
|trnsfrs.f            |3D fluxes of salt solutes                                                            |
|uptake.f             |C, N, P, water and heat exchange between root and soil, and atmosphere               |
|visual.f             |Output for model check                                                               |
|watsub.f             |Soil hydrothermal processes                                                          |
|wthr.f               |Process weather variables                                                            |
|solutepar.h          |Parameters for the solute model                                                      |
