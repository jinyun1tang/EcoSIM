module TranspNoSaltDataMod
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use GridConsts
  use TracerIDMod
  use EcoSIMConfig, only : jcplx => jcplxc
implicit none
  CHARACTER(LEN=*), Private,PARAMETER :: MOD_FILENAME=&
  __FILE__
  real(r8), parameter :: tiny_p=1.e-3_r8
  real(r8), parameter :: XFRS=0.05_r8                             !maximum fraction of total soil volume to act as macropore for exchange with micropore
  real(r8), Parameter :: VFLWX=0.5_r8                             !default outflow in water for micropore-macropore exchange

  real(r8), allocatable :: trcg_AquaADV_Snow2Soil_flxM(:,:,:)     !gas advection flux from snow to top soil at iteration M [g d-2]
  real(r8), allocatable :: trcn_AquaADV_Snow2Band_flxM(:,:,:)     !nutrient advection flux from snow to top band soil at iteration M [g d-2]
  real(r8), allocatable :: trcn_AquaADV_Snow2Soil_flxM(:,:,:)     !nutrient advection flux from snow to top soil at iteration M [g d-2]
  real(r8), allocatable :: Gas_AdvDif_FlxMM_2DH(:,:,:)            !2D horizontal gas fluxes
  real(r8), allocatable :: DOM_transpFlxM_2DH(:,:,:,:)            !2D horizontal DOM fluxes [g d-2]
  real(r8), allocatable :: trcs_transpFlxM_2DH(:,:,:)             !2D horizontal solute fluxes [g d-2]
  real(r8), allocatable :: DOM_MicP2_vr(:,:,:,:,:)                !copy of DOM in micropore for transport iteration [g d-2]
  real(r8), allocatable :: DOM_MicpTranspFlxM_3D(:,:,:,:,:,:)     !3D DOM transport through micropores at iteration M [g d-2]
  real(r8), allocatable :: DOM_MacpTranspFlxM_3D(:,:,:,:,:,:)     !3D DOM transport through macropores at iteration M [g d-2]
  real(r8), allocatable :: DOM_Transp2Micp_flxM_vr(:,:,:,:,:)     !DOM transport to micropore at iteration M [g d-2]
  real(r8), allocatable :: DOM_Transp2Macp_flxM_vr(:,:,:,:,:)     !DOM transport to macropore at iteration M [g d-2]  
  real(r8), allocatable :: RDOM_CumEcoProd_vr(:,:,:,:,:)
  real(r8), allocatable :: WaterFlow2SoilMM_3D(:,:,:,:)           !3D water flow used for fast iteration 
  real(r8), allocatable :: GasDiff2Litr_flx_col(:,:,:)            !diffusive gas flux to litter [g d-2 h-1]
  real(r8), allocatable :: RGasAtmDisol2SoilM_col(:,:,:)          !Gas dissolution from atmosphere into soil at iteration M
  real(r8), allocatable :: RGasAtmDisol2LitrM_col(:,:,:)          !Gas dissolution from atmosphere into litter at iteration M


  real(r8), allocatable :: DOM_FloXSurRunoff_PotFlxM(:,:,:,:)      !Potential DOM flux through surface runoff at iteration M [g d-2]
  real(r8), allocatable :: RBGCSinkGasMM_vr(:,:,:,:)               !BGC sink for gaseous tracers at gas iteration MM
  real(r8), allocatable :: RBGCSinkSoluteM_vr(:,:,:,:)             !BGC sink for solute tracers at heat iteration M
  real(r8), allocatable :: trcg_AquaAdv_flxM_snvr(:,:,:,:)         !gas advection by water advection in snow at iteration M [g d-2]                   
  real(r8), allocatable :: trcn_AquaAdv_flxM_snvr(:,:,:,:)         !solute advection by water advection in snow at iteration M [g d-2]
  real(r8), allocatable :: DOM_Flo2LitrM(:,:,:,:)                  !DOM flow into litter layer at iteration M [gC d-2]
  real(r8), allocatable :: DOM_Flo2TopSoilM(:,:,:,:)               !DOM flow into topsoil layer at iteration M [gC d-2]
  real(r8), allocatable :: trcg_Precip2LitrM_col(:,:,:)                !wet precipiation gas flux to litter at iteration M [g d-2]
  real(r8), allocatable :: trcn_Precip2LitrM_col(:,:,:)                !wet precipiation nutrient flux to litter at iteration M [g d-2]
  real(r8), allocatable :: trcs_Precip2MicpM_col(:,:,:)                !wet precipiation nutrient flux to soil micropore at iteration M [g d-2]
  real(r8), allocatable :: DOM_Mac2MicPore_flxM_vr(:,:,:,:,:)      !DOM exchange from macropore and micropore at iteration M [g d-2]

  real(r8), allocatable ::  DOM_SurfRunoff_flxM(:,:,:,:)           !DOM incoming flux to grid from surface runoff at iteration M [g d-2]

  real(r8), allocatable ::  trcg_SnowDrift_flxM(:,:,:)             !gas incoming flux to grid from snow drift at iteration M [g d-2]
  real(r8), allocatable ::  trcn_SnowDrift_flxM(:,:,:)             !nutrient incoming flux to grid from snow drift at iteration M [g d-2]
  real(r8), allocatable ::  trcg_SurfRunoff_flxM(:,:,:)            !Gas incoming flux to grid from surface runoff at iteration M [g d-2]
  real(r8), allocatable ::  trcn_SurfRunoff_flxM(:,:,:)            !Nutrient incoming flux to grid from surface runoff at iteration M [g d-2]
  real(r8), allocatable ::  trcs_Transp2Micp_flxM_vr(:,:,:,:)      !total 3D micropore flux at iteration M [g d-2]
  real(r8), allocatable ::  trcs_Transp2Macp_flxM_vr(:,:,:,:)      !total 3D macropore flux at iteration M [g d-2]
  real(r8), allocatable ::  trcsol_Irrig_flxM_vr(:,:,:,:)          !solute flux due to irrigation
  real(r8), allocatable ::  GasDifctScaledMM_vr(:,:,:,:)           !time scaled gas diffusivity
  real(r8), allocatable ::  SoluteDifusivitytscaledM_vr(:,:,:,:)   !time scaled solute diffusivity
  real(r8), allocatable ::  trcg_VLWatMicP_vr(:,:,:,:)             !solubility-scaled effective volume for dissolved gas species [m3 d-2]

  real(r8), allocatable ::  tScalReductVLsoiAirPMM_vr(:,:,:)          !time scaled air-filled reduction volume [m3 h /d2 ]
  real(r8), allocatable ::  VLsoiAirPMA_vr(:,:,:)                      !Effective air-filled soil pore volume occupied by NH4 [m3]
  real(r8), allocatable ::  VLsoiAirPMB_vr(:,:,:)                      !Effective air-filled soil pore volume occupied by band-NH4 [m3]
  real(r8), allocatable ::  VLWatMicPMA_vr(:,:,:)                      !Effective water-filled soil micropore volume occupied by NH4 [m3]
  real(r8), allocatable ::  VLWatMicPMB_vr(:,:,:)                      !Effective air-filled soil micropore volume occupied by band-NH4 [m3]
  real(r8), allocatable ::  VLWatMicPXA_vr(:,:,:)                      !maximum non-band NH4 mass in micropore [gN /d2]
  real(r8), allocatable ::  VLWatMicPXB_vr(:,:,:)                      !maximum band NH4 mass in micropore [gN /d2]
  real(r8), allocatable ::  PARGas_CefMM(:,:,:)                        !surface conductance for gas diffusion []
  real(r8), allocatable ::  Gas_WetDepo2Litr_col(:,:,:)                !wet deposition of gas to litter [g d-2 h-1]
  real(r8), allocatable ::  trcg_Ebu_flxM_vr(:,:,:,:)                !
  real(r8), allocatable ::  trcg_FloXSurRunoff_PotFlxM(:,:,:)        !Potential gas flux through surface runoff at iteration M [g d-2]
  real(r8), allocatable ::  trcn_FloXSurRunoff_PotFlxM(:,:,:)           !Potential nutrient flux through surface runoff at iteration M [g d-2]
  real(r8), allocatable ::  DOM_MacP2_vr(:,:,:,:,:)                     !copy of DOM in macropore for transport iteration [g d-2]
  real(r8), allocatable ::  DOM_diffusivitytscaledM_vr(:,:,:,:)      !time scaled DOM diffusivity at iteration M 

  real(r8), allocatable :: trcg_solsml2_snvr(:,:,:,:)                   !copy of dissolved gas tracer in snow for transport [g d-2]
  real(r8), allocatable :: trcn_solsml2_snvr(:,:,:,:)                   !copy of solute tracer in snow for transport [g d-2]
  real(r8), allocatable ::  trcg_Aqua_flxM_snvr(:,:,:,:)                !Aquatic gas flow into snow layer at iteration M [g d-2]
  real(r8), allocatable ::  trcn_Aqua_flxM_snvr(:,:,:,:)                !Aquatic nutrient flow into snow layer at iteration M [g d-2]
  real(r8), allocatable ::  DOM_FloXSurRof_flxM_2DH(:,:,:,:,:,:)        !2D flow of DOM through surface runoff at iteration M [g d-2]
  real(r8), allocatable ::  trcg_SurRof_flxM_2DH(:,:,:,:,:)             !2D flow of gas through surface runoff at iteration M [g d-2]
  real(r8), allocatable ::  trcn_SurRof_flxM_2DH(:,:,:,:,:)             !2D flow of nutrient through surface runoff at iteration M [g d-2]
  real(r8), allocatable :: trcg_gasml2_vr(:,:,:,:)                       !copy of gas tracer for transport [g d-2]
  real(r8), allocatable :: trcs_solml2_vr(:,:,:,:)                      !copy of solute tracer in micropore for transport [g d-2]
  real(r8), allocatable :: trcs_soHml2_vr(:,:,:,:)                      !copy of solute tracer in macropore for transport [g d-2]
  real(r8), allocatable ::  trcn_SnowDrift_flxM_2DH(:,:,:,:,:)           !nutrient flux through snow drift at iteration M [g d-2]
  real(r8), allocatable ::  trcg_SnowDrift_flxM_2DH(:,:,:,:,:)           !gas flux through snow drift at iteration M  [g d-2]                  
  real(r8), allocatable ::  RGas_Disol_FlxMM_vr(:,:,:,:)                !gas dissolution (>0 into aqueous phase) at iteration MM [g d-2]
  real(r8), allocatable ::  Gas_AdvDif_FlxMM_3D(:,:,:,:,:)              !3D gas flux advection + diffusion at iteration MM [g d-2]
  real(r8), allocatable ::  Gas_AdvDif_FlxMM_vr(:,:,:,:)                !total 3D gas flux advection + diffusion at iteration MM [g d-2]
  real(r8), allocatable :: trcs_MacpTranspFlxM_3D(:,:,:,:,:)            !solute 3D macropore flux at iteration M [g d-2]
  real(r8), allocatable :: trcs_MicpTranspFlxM_3D(:,:,:,:,:)            !solute 3D micropore flux at iteration M [g d-2]
  real(r8), allocatable :: trcs_Mac2MicPore_flxM_vr(:,:,:,:)            !Mac to micropore nutrient flux at iteration M [g d-2]
  real(r8), allocatable :: trcg_AquaADV_Snow2Litr_flxM(:,:,:)           !aqueous gas flow from snow to litr [g d-2]
  real(r8), allocatable :: trcn_AquaADV_Snow2Litr_flxM(:,:,:)           !aqueous nutrient flow from snow to litr [g d-2]

  real(r8), allocatable :: Gas_litr2Soil_flx_col(:,:,:)                 !gas flux from litter to soil [g d-2 h-1]
  real(r8), allocatable :: Gas_litr2Soil_flxM_col(:,:,:)                !gas flux from litter to soil at iteration M [g d-2]
  real(r8), allocatable :: GasHydroLoss_litr_flx_col(:,:,:)             !gas flux loss through surface runoff [g d-2 h-1]
  real(r8), allocatable :: GasDiff2Soil_flx_col(:,:,:)                  !gas diffusion flux into soil [g d-2 h-1]
  real(r8), allocatable :: Gas_WetDepo2Soil_col(:,:,:)                  !gas wet deposition flux into soil [g d-2 h-1]
  real(r8), allocatable :: RGasNetProdSoil_col(:,:,:)                   !gas production in soil [g d-2 h-1]
  real(r8), allocatable :: trc_topsoil_flx_col(:,:,:)                   !net tracer flux at the top of soil [g d-2 -h1]
  real(r8), allocatable :: TranspNetSoil_flx_col(:,:,:)                 !total tracer flux due to transport in soil [g d-2 h-1]
  real(r8), allocatable :: TranspNetSoil_flx2_col(:,:,:)                !total tracer flux due to transport in soil at iteration M [g d-2 h-1]
  real(r8), allocatable :: TranspNetSoil_fast_flxM_col(:,:,:)           !total tracer flux due to fast transport in soil at iteration M [g d-2 h-1]
  real(r8), allocatable :: TranspNetSoil_slow_flxM_col(:,:,:)           !total tracer flux due to slow transport in soil at iteration M [g d-2 h-1]
  real(r8), allocatable :: transp_diff_slow_vr(:,:,:,:)                 !error in computing the transport profile [g d-2 h-1]
  real(r8), allocatable :: Gas_WetDepo2Snow_col(:,:,:)                  !wet deposition to snow [g d-2 h-1] 
  real(r8), allocatable :: Gas_Snowloss_flx_col(:,:,:)                      !tracer loss from snow [g d-2 h-1]
  real(r8), allocatable :: trcg_mass_begf(:,:,:)                        !initial tracer mass check in fast transport
  real(r8), allocatable :: trcg_mass_begs(:,:,:)                        !initial tracer mass check in slow transport
  real(r8), allocatable :: trcg_mass_snow_begs(:,:,:)                   !tracer in snow during slow iteration [g d-2]
  real(r8), allocatable :: trcg_mass_litr_begs(:,:,:)                   !tracer in litter during slow iteration [g d-2]
  real(r8), allocatable :: trcg_mass_soil_begs(:,:,:)                   !tracer in soil during slow iteration [g d-2]
  real(r8), allocatable :: trcg_netflx2soil_col(:,:,:)                  !net tracer flux into soil during gas transport
  real(r8), allocatable :: TranspNetSoil_fast_flx_col(:,:,:)            !
  real(r8), allocatable :: RGas_Disol_FlxM_vr(:,:,:,:)
  real(r8), allocatable :: RGasSinkScalar_vr(:,:,:,:)                   !BGC sink scalar for numerical stability, [none]
  real(r8), allocatable :: DOM_mass_begs(:,:,:,:)                       !DOM  mass check in slow transport
  real(r8), allocatable :: trcs_mass_beg(:,:,:)                         !initial tracer mass before doing transport
  real(r8), allocatable :: trcs_mass_litr(:,:,:)
  real(r8), allocatable :: trcs_mass_snow(:,:,:)
  real(r8), allocatable :: trcs_mass_soil(:,:,:)  
  real(r8), allocatable :: DOM_Hydroloss_col(:,:,:,:)                  !total dom loss through hydrological pathways
  real(r8), allocatable :: TranspNetDOM_flx_col(:,:,:,:)               !net dom flux during slow transport
  real(r8), allocatable :: TranspNetDOM_flx2_col(:,:,:,:)              !net dom flux during slow transport  
  real(r8), allocatable :: TranspNetDOM_flxM_col(:,:,:,:)              !net dom flux during iteration of slow transport
  real(r8), allocatable :: DOM_mass3_col(:,:,:,:)                      !total dom mass for conservation check
  real(r8), allocatable :: Gas_WetDeposit_slow_flx_col(:,:,:)          !wet deposition in one iteration of slow transport [g d-2]
  real(r8), allocatable :: GasDiff2Surf_slow_flx_col(:,:,:)            !surface diffusion in one iteration of slow transport [g d-2]
  real(r8), allocatable :: trcs_NetProd_slow_flx_col(:,:,:)            !net tracer production flux during slow transport [g d-2]
  real(r8), allocatable :: trcs_hydrloss_slow_flx_col(:,:,:)           !tracer loss through snow drift, surf/subsurf runoff, and drainage in slow transport [g d-2]
  real(r8), allocatable :: Gas_WetDepo2Snow_slow_flx_col(:,:,:)
  real(r8), allocatable :: Gas_Snowloss_slow_flx_col(:,:,:)
  real(r8), allocatable :: trcs_NetFlow2Litr_slow_flx_col(:,:,:)       !total flux into litter layer during one iteration of slow transport [g d-2]
  real(r8), allocatable :: trcs_netflow2soil_slow_flx_col(:,:,:)       !net tracer flow to soil in one iteration of slow transport [g d-2]
!----------------------------------------------------------------------

contains
  subroutine InitTransfrData
  use EcoSiMParDataMod, only : micpar
  implicit none
  integer :: NumOfLitrCmplxs

  NumOfLitrCmplxs=micpar%NumOfLitrCmplxs
  allocate(trcs_netflow2soil_slow_flx_col(ids_beg:ids_end,JY,JX)); trcs_netflow2soil_slow_flx_col=0._r8
  allocate(trcs_NetFlow2Litr_slow_flx_col(idg_beg:idg_NH3,JY,JX)); trcs_NetFlow2Litr_slow_flx_col=0._r8
  allocate(trcg_mass_snow_begs(idg_beg:idg_NH3,JY,JX)); trcg_mass_snow_begs=0._r8
  allocate(trcg_mass_litr_begs(idg_beg:idg_NH3,JY,JX)); trcg_mass_litr_begs=0._r8
  allocate(trcg_mass_soil_begs(idg_beg:idg_NH3,JY,JX)); trcg_mass_soil_begs=0._r8
  allocate(trcs_hydrloss_slow_flx_col(ids_beg:ids_end,JY,JX)); trcs_hydrloss_slow_flx_col=0._r8
  allocate(trcs_NetProd_slow_flx_col(ids_beg:ids_end,JY,JX)); trcs_NetProd_slow_flx_col=0._r8
  allocate(GasDiff2Surf_slow_flx_col(idg_beg:idg_NH3,JY,JX)); GasDiff2Surf_slow_flx_col=0._r8
  allocate(Gas_WetDeposit_slow_flx_col(idg_beg:idg_NH3,JY,JX)); Gas_WetDeposit_slow_flx_col=0._r8
  allocate(DOM_mass3_col(idom_beg:idom_end,jcplx,JY,JX)); DOM_mass3_col=0._r8
  allocate(TranspNetDOM_flx_col(idom_beg:idom_end,jcplx,JY,JX));TranspNetDOM_flx_col=0._r8
  allocate(TranspNetDOM_flx2_col(idom_beg:idom_end,jcplx,JY,JX));TranspNetDOM_flx2_col=0._r8
  allocate(Gas_WetDepo2Snow_slow_flx_col(idg_beg:idg_NH3,JY,JX));Gas_WetDepo2Snow_slow_flx_col=0._r8
  allocate(Gas_Snowloss_slow_flx_col(idg_beg:idg_NH3,JY,JX));Gas_Snowloss_slow_flx_col=0._r8

  allocate(TranspNetDOM_flxM_col(idom_beg:idom_end,jcplx,JY,JX));TranspNetDOM_flxM_col=0._r8
  allocate(DOM_Hydroloss_col(idom_beg:idom_end,jcplx,JY,JX)); DOM_Hydroloss_col=0._r8
  allocate(DOM_mass_begs(idom_beg:idom_end,jcplx,JY,JX));DOM_mass_begs=0._r8
  allocate(trcg_mass_begs(idg_beg:idg_NH3,JY,JX)); trcg_mass_begs=0._r8
  allocate(RGasSinkScalar_vr(idg_beg:idg_NH3,0:JZ,JY,JX));RGasSinkScalar_vr=1._r8
  allocate(RGas_Disol_FlxM_vr(idg_beg:idg_NH3,JZ,JY,JX));RGas_Disol_FlxM_vr=0._r8
  allocate(TranspNetSoil_flx2_col(idg_beg:idg_NH3,JY,JX)); TranspNetSoil_flx2_col=0._r8
  allocate(trcg_mass_begf(idg_beg:idg_NH3,JY,JX)); trcg_mass_begf=0._r8
  allocate(Gas_WetDepo2Snow_col(idg_beg:idg_NH3,JY,JX)); Gas_WetDepo2Snow_col=0._r8
  allocate(Gas_Snowloss_flx_col(idg_beg:idg_NH3,JY,JX)); Gas_Snowloss_flx_col=0._r8
  allocate(transp_diff_slow_vr(idg_beg:idg_NH3,JZ,JY,JX)); transp_diff_slow_vr=0._r8
  allocate(TranspNetSoil_fast_flx_col(idg_beg:idg_NH3,JY,JX)); TranspNetSoil_fast_flx_col=0._r8
  allocate(trc_topsoil_flx_col(idg_beg:idg_NH3,JY,JX)); trc_topsoil_flx_col=0._r8
  allocate(GasDiff2Soil_flx_col(idg_beg:idg_NH3,JY,JX));GasDiff2Soil_flx_col=0._r8
  allocate(Gas_WetDepo2Soil_col(idg_beg:idg_NH3,JY,JX));Gas_WetDepo2Soil_col=0._r8
  allocate(RGasNetProdSoil_col(idg_beg:idg_NH3,JY,JX)); RGasNetProdSoil_col=0._r8
  allocate(GasHydroLoss_litr_flx_col(idg_beg:idg_NH3,JY,JX)); GasHydroLoss_litr_flx_col=0._r8
  allocate(Gas_litr2Soil_flxM_col(idg_beg:idg_end,JY,JX)); Gas_litr2Soil_flxM_col=0._r8
  allocate(GasDiff2Litr_flx_col(idg_beg:idg_NH3,JY,JX)); GasDiff2Litr_flx_col=0._r8
  allocate(trcs_mass_litr(idg_beg:idg_NH3,JY,JX));trcs_mass_litr=0._r8
  allocate(trcs_mass_snow(idg_beg:idg_NH3,JY,JX));trcs_mass_snow=0._r8
  allocate(trcs_mass_soil(idg_beg:idg_NH3,JY,JX));trcs_mass_soil=0._r8
  allocate(Gas_WetDepo2Litr_col(idg_beg:idg_NH3,JY,JX)); Gas_WetDepo2Litr_col=0._r8
  allocate(Gas_litr2Soil_flx_col(idg_beg:idg_NH3,JY,JX)); Gas_litr2Soil_flx_col=0._r8
  allocate(trcs_mass_beg(ids_beg:ids_end,JY,JX)); trcs_mass_beg=0._r8
  allocate(trcsol_Irrig_flxM_vr(ids_beg:ids_end,JZ,JY,JX)); trcsol_Irrig_flxM_vr=0._r8  !solute flux due to irrigation
  allocate(RGasAtmDisol2SoilM_col(idg_beg:idg_end,JY,JX));RGasAtmDisol2SoilM_col=0._r8
  allocate(RGasAtmDisol2LitrM_col(idg_beg:idg_NH3,JY,JX));RGasAtmDisol2LitrM_col=0._r8

  allocate(WaterFlow2SoilMM_3D(3,JD,JV,JH)); WaterFlow2SoilMM_3D=0._r8
  allocate(trcg_AquaADV_Snow2Litr_flxM(idg_beg:idg_NH3,JY,JX));   trcg_AquaADV_Snow2Litr_flxM=0._r8
  allocate(trcn_AquaADV_Snow2Litr_flxM(ids_nut_beg:ids_nuts_end,JY,JX)); trcn_AquaADV_Snow2Litr_flxM=0._r8
  allocate(DOM_FloXSurRunoff_PotFlxM(idom_beg:idom_end,1:jcplx,JV,JH));    DOM_FloXSurRunoff_PotFlxM=0._r8
  allocate(DOM_MicP2_vr(idom_beg:idom_end,1:jcplx,0:JZ,JY,JX));    DOM_MicP2_vr=0._r8
  allocate(RDOM_CumEcoProd_vr(idom_beg:idom_end,1:jcplx,0:JZ,JY,JX));  RDOM_CumEcoProd_vr=0._r8
  allocate(DOM_MicpTranspFlxM_3D(idom_beg:idom_end,1:jcplx,3,0:JD,JV,JH));DOM_MicpTranspFlxM_3D=0._r8
  allocate(DOM_MacpTranspFlxM_3D(idom_beg:idom_end,1:jcplx,3,JD,JV,JH));  DOM_MacpTranspFlxM_3D=0._r8
  allocate(DOM_Transp2Micp_flxM_vr(idom_beg:idom_end,1:jcplx,JZ,JY,JX));    DOM_Transp2Micp_flxM_vr=0._r8

  allocate(trcg_AquaAdv_flxM_snvr(idg_beg:idg_NH3,JS,JY,JX));   trcg_AquaAdv_flxM_snvr=0._r8
  allocate(trcn_AquaAdv_flxM_snvr(ids_nut_beg:ids_nuts_end,JS,JY,JX));   trcn_AquaAdv_flxM_snvr=0._r8
  allocate(DOM_Flo2LitrM(idom_beg:idom_end,1:NumOfLitrCmplxs,JY,JX));  DOM_Flo2LitrM=0._r8
  allocate(DOM_Flo2TopSoilM(idom_beg:idom_end,1:NumOfLitrCmplxs,JY,JX));  DOM_Flo2TopSoilM=0._r8
  allocate(trcg_Precip2LitrM_col(ids_beg:ids_end,JY,JX));      trcg_Precip2LitrM_col=0._r8
  allocate(trcn_Precip2LitrM_col(ids_nut_beg:ids_nuts_end,JY,JX));      trcn_Precip2LitrM_col=0._r8
  allocate(trcs_Precip2MicpM_col(ids_beg:ids_end,JY,JX));      trcs_Precip2MicpM_col=0._r8
  allocate(DOM_Mac2MicPore_flxM_vr(idom_beg:idom_end,1:jcplx,JZ,JY,JX));DOM_Mac2MicPore_flxM_vr=0._r8
  allocate(DOM_SurfRunoff_flxM(idom_beg:idom_end,1:jcplx,JY,JX)); DOM_SurfRunoff_flxM=0._r8

  allocate(trcg_SnowDrift_flxM(idg_beg:idg_NH3,JY,JX));      trcg_SnowDrift_flxM=0._r8
  allocate(trcn_SnowDrift_flxM(ids_nut_beg:ids_nuts_end,JY,JX)); trcn_SnowDrift_flxM=0._r8
  allocate(trcg_SurfRunoff_flxM(idg_beg:idg_NH3,JY,JX)); trcg_SurfRunoff_flxM=0._r8
  allocate(trcn_SurfRunoff_flxM(ids_nut_beg:ids_nuts_end,JY,JX));trcn_SurfRunoff_flxM=0._r8

  allocate(GasDifctScaledMM_vr(idg_beg:idg_end,JZ,JY,JX)); GasDifctScaledMM_vr=0._r8
  allocate(SoluteDifusivitytscaledM_vr(ids_beg:ids_end,0:JZ,JY,JX));SoluteDifusivitytscaledM_vr=0._r8
  allocate(trcg_VLWatMicP_vr(idg_beg:idg_end,0:JZ,JY,JX)); trcg_VLWatMicP_vr=0._r8

  allocate(trcs_transpFlxM_2DH(ids_beg:ids_end,JY,JX)); trcs_transpFlxM_2DH=0._r8
  allocate(DOM_transpFlxM_2DH(idom_beg:idom_end,jcplx,JY,JX));DOM_transpFlxM_2DH=0._r8
  allocate(Gas_AdvDif_FlxMM_2DH(idg_beg:idg_NH3,JY,JX)); Gas_AdvDif_FlxMM_2DH=0._r8
  allocate(tScalReductVLsoiAirPMM_vr(JZ,JY,JX));     tScalReductVLsoiAirPMM_vr=0._r8
  allocate(VLsoiAirPMA_vr(JZ,JY,JX));   VLsoiAirPMA_vr=0._r8
  allocate(VLsoiAirPMB_vr(JZ,JY,JX));   VLsoiAirPMB_vr=0._r8
  allocate(VLWatMicPMA_vr(JZ,JY,JX));   VLWatMicPMA_vr=0._r8
  allocate(VLWatMicPMB_vr(JZ,JY,JX));   VLWatMicPMB_vr=0._r8
  allocate(VLWatMicPXA_vr(0:JZ,JY,JX)); VLWatMicPXA_vr=0._r8
  allocate(VLWatMicPXB_vr(JZ,JY,JX));   VLWatMicPXB_vr=0._r8
  allocate(PARGas_CefMM(idg_beg:idg_NH3,JY,JX));      PARGas_CefMM=0._r8

  allocate(RBGCSinkGasMM_vr(idg_beg:idg_NH3,0:JZ,JY,JX));RBGCSinkGasMM_vr=0._r8
  allocate(RBGCSinkSoluteM_vr(ids_nuts_beg:ids_nuts_end,0:JZ,JY,JX));RBGCSinkSoluteM_vr=0._r8

  allocate(trcg_Ebu_flxM_vr(idg_beg:idg_end,JZ,JY,JX));   trcg_Ebu_flxM_vr                        = 0._r8
  allocate(trcg_FloXSurRunoff_PotFlxM(idg_beg:idg_end,JV,JH));     trcg_FloXSurRunoff_PotFlxM           = 0._r8
  allocate(trcn_FloXSurRunoff_PotFlxM(ids_nut_beg:ids_nuts_end,JV,JH));     trcn_FloXSurRunoff_PotFlxM  = 0._r8
  allocate(DOM_MacP2_vr(idom_beg:idom_end,1:jcplx,JZ,JY,JX));DOM_MacP2_vr                     = 0._r8
  allocate(DOM_Transp2Macp_flxM_vr(idom_beg:idom_end,1:jcplx,JZ,JY,JX));DOM_Transp2Macp_flxM_vr = 0._r8

  allocate(DOM_diffusivitytscaledM_vr(idom_beg:idom_end,0:JZ,JY,JX)); DOM_diffusivitytscaledM_vr = 0._r8

  allocate(trcg_solsml2_snvr(idg_beg:idg_NH3,JS,JY,JX));trcg_solsml2_snvr = 0._r8
  allocate(trcn_solsml2_snvr(ids_nut_beg:ids_nuts_end,JS,JY,JX));trcn_solsml2_snvr  = 0._r8

  allocate(trcg_Aqua_flxM_snvr(idg_beg:idg_NH3,JS,JY,JX));   trcg_Aqua_flxM_snvr = 0._r8
  allocate(trcn_Aqua_flxM_snvr(ids_nut_beg:ids_nuts_end,JS,JY,JX));trcn_Aqua_flxM_snvr     = 0._r8

  allocate(DOM_FloXSurRof_flxM_2DH(idom_beg:idom_end,1:jcplx,2,2,JV,JH));DOM_FloXSurRof_flxM_2DH = 0._r8
  allocate(trcg_SurRof_flxM_2DH(idg_beg:idg_NH3,2,2,JV,JH));  trcg_SurRof_flxM_2DH       = 0._r8

  allocate(trcg_gasml2_vr(idg_beg:idg_NH3,0:JZ,JY,JX)); trcg_gasml2_vr(:,:,:,:) = 0._r8
  allocate(trcs_solml2_vr(ids_beg:ids_end,0:JZ,JY,JX)); trcs_solml2_vr(:,:,:,:) = 0._r8
  allocate(trcs_soHml2_vr(ids_beg:ids_end,0:JZ,JY,JX)); trcs_soHml2_vr(:,:,:,:) = 0._r8

  allocate(trcg_SnowDrift_flxM_2DH(idg_beg:idg_NH3,2,2,JV,JH));    trcg_SnowDrift_flxM_2DH           = 0._r8
  allocate(trcn_SurRof_flxM_2DH(ids_nut_beg:ids_nuts_end,2,2,JV,JH));  trcn_SurRof_flxM_2DH = 0._r8
  allocate(Gas_AdvDif_FlxMM_3D(idg_beg:idg_NH3,3,JD,JV,JH));Gas_AdvDif_FlxMM_3D                             = 0._r8
  allocate(Gas_AdvDif_FlxMM_vr(idg_beg:idg_NH3,JZ,JY,JH));Gas_AdvDif_FlxMM_vr                     = 0._r8

  allocate(RGas_Disol_FlxMM_vr(idg_beg:idg_end,0:JZ,JY,JX)); RGas_Disol_FlxMM_vr = 0._r8

  allocate(trcn_SnowDrift_flxM_2DH(ids_nut_beg:ids_nuts_end,2,2,JV,JH)); trcn_SnowDrift_flxM_2DH = 0._r8
  allocate(trcs_MacpTranspFlxM_3D(ids_beg:ids_end,3,JD,JV,JH));trcs_MacpTranspFlxM_3D       = 0._r8
  allocate(trcs_MicpTranspFlxM_3D(ids_beg:ids_end,3,0:JD,JV,JH));trcs_MicpTranspFlxM_3D     = 0._r8
  allocate(trcs_Mac2MicPore_flxM_vr(ids_beg:ids_end,JZ,JY,JX));trcs_Mac2MicPore_flxM_vr     = 0._r8
  allocate(trcs_Transp2Micp_flxM_vr(ids_beg:ids_end,JZ,JY,JX));trcs_Transp2Micp_flxM_vr = 0._r8
  allocate(trcs_Transp2Macp_flxM_vr(ids_beg:ids_end,JZ,JY,JX));trcs_Transp2Macp_flxM_vr = 0._r8
  allocate(trcg_AquaADV_Snow2Soil_flxM(idg_beg:idg_NH3,JY,JX));trcg_AquaADV_Snow2Soil_flxM   = 0._r8
  allocate(trcn_AquaADV_Snow2Band_flxM(ids_nutb_beg:ids_nutb_end,JY,JX));trcn_AquaADV_Snow2Band_flxM  = 0._r8
  allocate(trcn_AquaADV_Snow2Soil_flxM(ids_nut_beg:ids_nuts_end,JY,JX));trcn_AquaADV_Snow2Soil_flxM   = 0._r8
  allocate(TranspNetSoil_flx_col(idg_beg:idg_NH3,JY,JX)); TranspNetSoil_flx_col=0._r8
  allocate(TranspNetSoil_fast_flxM_col(idg_beg:idg_NH3,JY,JX)); TranspNetSoil_fast_flxM_col=0._r8
  allocate(TranspNetSoil_slow_flxM_col(idg_beg:idg_end,JY,JX)); TranspNetSoil_slow_flxM_col=0._r8
  allocate(trcg_netflx2soil_col(idg_beg:idg_NH3,JY,JX)); trcg_netflx2soil_col=0._r8

  end subroutine InitTransfrData

!----------------------------------------------------------------------
  subroutine DestructTransfrData
  use abortutils, only : destroy
  implicit none

  call destroy(trcs_netflow2soil_slow_flx_col)
  call destroy(Gas_WetDepo2Snow_slow_flx_col)
  call destroy(Gas_Snowloss_slow_flx_col)
  call destroy(trcs_NetFlow2Litr_slow_flx_col)
  call destroy(trcg_mass_snow_begs)
  call destroy(trcg_mass_litr_begs)
  call destroy(trcg_mass_soil_begs)
  call destroy(trcs_hydrloss_slow_flx_col)
  call destroy(trcs_NetProd_slow_flx_col)
  call destroy(GasDiff2Surf_slow_flx_col)
  call destroy(Gas_WetDeposit_slow_flx_col)
  call destroy(DOM_mass3_col)
  call destroy(TranspNetDOM_flx2_col)
  call destroy(TranspNetDOM_flx_col)
  call destroy(TranspNetDOM_flxM_col)
  call destroy(DOM_Hydroloss_col)
  call destroy(DOM_mass_begs)
  call destroy(trcg_mass_begs)
  call destroy(RGasSinkScalar_vr)
  call destroy(RGas_Disol_FlxM_vr)
  call destroy(TranspNetSoil_flx2_col)
  call destroy(TranspNetSoil_fast_flx_col)  
  call destroy(trcg_netflx2soil_col)
  call destroy(trcg_mass_begf)
  call destroy(transp_diff_slow_vr)
  call destroy(TranspNetSoil_flx_col)
  call destroy(trc_topsoil_flx_col)
  call destroy(GasDiff2Soil_flx_col)
  call destroy(Gas_WetDepo2Soil_col)
  call destroy(RGasNetProdSoil_col)
  call destroy(GasHydroLoss_litr_flx_col)
  call destroy(Gas_litr2Soil_flxM_col)
  call destroy(Gas_litr2Soil_flx_col)
  call destroy(trcs_mass_litr)
  call destroy(trcs_mass_snow)
  call destroy(trcs_mass_soil)
  call destroy(GasDiff2Litr_flx_col)
  call destroy(Gas_WetDepo2Litr_col)
  call destroy(trcs_mass_beg)
  call destroy(RGasAtmDisol2SoilM_col)
  call destroy(RGasAtmDisol2LitrM_col)
  call destroy(trcs_transpFlxM_2DH)
  call destroy(DOM_transpFlxM_2DH)
  call destroy(Gas_AdvDif_FlxMM_2DH)
  call destroy(WaterFlow2SoilMM_3D)
  call destroy(trcg_AquaADV_Snow2Litr_flxM)
  call destroy(trcn_AquaADV_Snow2Litr_flxM)
  call destroy(GasDifctScaledMM_vr)
  call destroy(SoluteDifusivitytscaledM_vr)
  call destroy(DOM_MicP2_vr)
  call destroy(trcsol_Irrig_flxM_vr)
  call destroy(RDOM_CumEcoProd_vr)
  call destroy(DOM_MicpTranspFlxM_3D)
  call destroy(DOM_MacpTranspFlxM_3D)
  call destroy(RBGCSinkSoluteM_vr)
  call destroy(DOM_FloXSurRunoff_PotFlxM)
  call destroy(DOM_FloXSurRof_flxM_2DH)
  call destroy(DOM_Flo2LitrM)
  call destroy(DOM_Flo2TopSoilM)
  call destroy(trcg_AquaADV_Snow2Soil_flxM)
  call destroy(trcs_Precip2MicpM_col)
  call destroy(trcn_Precip2LitrM_col)
  call destroy(DOM_Mac2MicPore_flxM_vr)
  call destroy(trcg_AquaAdv_flxM_snvr)
  call destroy(trcn_AquaAdv_flxM_snvr)
  call destroy(trcg_Ebu_flxM_vr)
  call destroy(RBGCSinkGasMM_vr)
  call destroy(trcg_FloXSurRunoff_PotFlxM)
  call destroy(DOM_SurfRunoff_flxM)
  call destroy(trcg_SurfRunoff_flxM)
  call destroy(trcn_SurfRunoff_flxM)
  call destroy(trcg_SnowDrift_flxM)
  call destroy(trcn_SnowDrift_flxM)
  call destroy(trcs_Transp2Micp_flxM_vr)
  call destroy(trcg_VLWatMicP_vr)
  call destroy(trcn_FloXSurRunoff_PotFlxM)
  call destroy(tScalReductVLsoiAirPMM_vr)
  call destroy(VLsoiAirPMA_vr)
  call destroy(VLsoiAirPMB_vr)
  call destroy(VLWatMicPMA_vr)
  call destroy(VLWatMicPMB_vr)
  call destroy(VLWatMicPXA_vr)
  call destroy(VLWatMicPXB_vr)
  call destroy(PARGas_CefMM)
  call destroy(DOM_diffusivitytscaledM_vr)
  call destroy(Gas_AdvDif_FlxMM_3D)
  call destroy(trcg_solsml2_snvr)
  call destroy(trcn_solsml2_snvr)
  call destroy(trcn_SurRof_flxM_2DH)
  call destroy(trcg_Aqua_flxM_snvr)
  call destroy(trcn_Aqua_flxM_snvr)
  call destroy(trcg_Precip2LitrM_col)
  call destroy(trcg_SnowDrift_flxM_2DH)
  call destroy(trcs_Transp2Macp_flxM_vr)
  call destroy(Gas_WetDepo2Snow_col)
  call destroy(Gas_Snowloss_flx_col)
  call destroy(DOM_MacP2_vr)
  call destroy(Gas_AdvDif_FlxMM_vr)
  call destroy(trcg_gasml2_vr)
  call destroy(trcs_solml2_vr)
  call destroy(trcs_soHml2_vr)
  call destroy(trcn_AquaADV_Snow2Band_flxM)
  call destroy(trcn_AquaADV_Snow2Soil_flxM)
  call destroy(trcg_SurRof_flxM_2DH)
  call destroy(RGas_Disol_FlxMM_vr)
  call destroy(trcn_SnowDrift_flxM_2DH)
  call destroy(trcs_MacpTranspFlxM_3D)
  call destroy(trcs_MicpTranspFlxM_3D)
  call destroy(trcs_Mac2MicPore_flxM_vr)
  call destroy(TranspNetSoil_flx2_col)
  call destroy(DOM_Transp2Macp_flxM_vr)
  call destroy(DOM_Transp2Micp_flxM_vr)
  call destroy(TranspNetSoil_slow_flxM_col)
  call destroy(TranspNetSoil_fast_flxM_col)
  end subroutine DestructTransfrData

end module TranspNoSaltDataMod
