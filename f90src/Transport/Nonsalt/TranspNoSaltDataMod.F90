module TranspNoSaltDataMod
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use GridConsts
  use TracerIDMod
  use EcoSIMConfig, only : jcplx => jcplxc
implicit none
  CHARACTER(LEN=*), Private,PARAMETER :: MOD_FILENAME=&
  __FILE__
  real(r8), parameter :: XFRS=0.05_r8
  real(r8), Parameter :: VFLWX=0.5_r8           !default outflow in water for micropore-macropore exchange

  real(r8), allocatable :: trcg_AquaADV_Snow2Soil_flxM(:)     !gas advection flux from snow to top soil at iteration M [g d-2]
  real(r8), allocatable :: trcn_AquaADV_Snow2Band_flxM(:)     !nutrient advection flux from snow to top band soil at iteration M [g d-2]
  real(r8), allocatable :: trcn_AquaADV_Snow2Soil_flxM(:)     !nutrient advection flux from snow to top soil at iteration M [g d-2]
  
  real(r8), allocatable :: DOM_MicP2(:,:,:,:,:)               !copy of DOM in micropore for transport iteration [g d-2]
  real(r8), allocatable :: DOM_MicpTranspFlxM_3D(:,:,:,:,:,:) !3D DOM transport through micropores at iteration M [g d-2]
  real(r8), allocatable :: DOM_MacpTranspFlxM_3D(:,:,:,:,:,:) !3D DOM transport through macropores at iteration M [g d-2]
  real(r8), allocatable :: DOM_Transp2Micp_flxM_vr(:,:,:,:,:) !DOM transport to micropore at iteration M [g d-2]
  real(r8), allocatable :: DOM_Transp2Macp_flxM_vr(:,:,:,:,:) !DOM transport to macropore at iteration M [g d-2]  
  real(r8), allocatable :: RDOM_CumEcoProd_vr(:,:,:,:,:)

  real(r8), allocatable :: DOM_FloXSurRunoff_flxM(:,:,:,:)    !DOM flux through surface runoff at iteration M [g d-2]
  real(r8), allocatable :: RBGCSinkGasMM_vr(:,:,:,:)               !BGC sink for gaseous tracers at gas iteration MM
  real(r8), allocatable :: RBGCSinkSoluteM_vr(:,:,:,:)             !BGC sink for solute tracers at heat iteration M
  real(r8), allocatable :: trcg_AquaAdv_flxM_snvr(:,:,:,:)         !gas advection by water advection in snow at iteration M [g d-2]                   
  real(r8), allocatable :: trcn_AquaAdv_flxM_snvr(:,:,:,:)         !solute advection by water advection in snow at iteration M [g d-2]
  real(r8), allocatable :: DOM_Flo2LitrM(:,:,:,:)                  !DOM flow into litter layer at iteration M [gC d-2]
  real(r8), allocatable :: DOM_Flo2TopSoilM(:,:,:,:)               !DOM flow into topsoil layer at iteration M [gC d-2]
  real(r8), allocatable :: trcg_RFL0(:,:,:)                        !
  real(r8), allocatable :: trcn_RFL0(:,:,:)                        !
  real(r8), allocatable :: trcs_RFL1(:,:,:)                        !
  real(r8), allocatable :: DOM_Mac2MicPore_flxM_vr(:,:,:,:,:)      !DOM exchange from macropore and micropore

  real(r8), allocatable ::  POSGL2(:,:,:)                          !
  real(r8), allocatable ::  O2AquaDiffusvity2(:,:,:)               !
  real(r8), allocatable ::  DOM_SurfRunoff_flxM(:,:,:,:)           !DOM incoming flux to grid from surface runoff at iteration M [g d-2]

  real(r8), allocatable ::  trcg_SnowDrift_flxM(:,:,:)             !gas incoming flux to grid from snow drift at iteration M [g d-2]
  real(r8), allocatable ::  trcn_SnowDrift_flxM(:,:,:)             !nutrient incoming flux to grid from snow drift at iteration M [g d-2]
  real(r8), allocatable ::  trcg_SurfRunoff_flxM(:,:,:)            !Gas incoming flux to grid from surface runoff at iteration M [g d-2]
  real(r8), allocatable ::  trcn_SurfRunoff_flxM(:,:,:)            !Nutrient incoming flux to grid from surface runoff at iteration M [g d-2]
  real(r8), allocatable ::  trcs_Transp2Micp_flxM_vr(:,:,:,:)      !total 3D micropore flux at iteration M [g d-2]
  real(r8), allocatable ::  trcs_Transp2Macp_flxM_vr(:,:,:,:)      !total 3D macropore flux at iteration M [g d-2]

  real(r8), allocatable ::  GasDifctScaled_vr(:,:,:,:)            !time scaled gas diffusivity
  real(r8), allocatable ::  SoluteDifusivitytscaled_vr(:,:,:,:)
  real(r8), allocatable ::  DifuscG_vr(:,:,:,:,:)
  real(r8), allocatable ::  trcg_VLWatMicP_vr(:,:,:,:)                      !

  real(r8), allocatable ::  DH2GG(:,:,:,:)                     !
  real(r8), allocatable ::  RHGFXS(:,:,:)                      !
  real(r8), allocatable ::  THGFLG(:,:,:)                      !
  real(r8), allocatable ::  CumReductVLsoiAirPM(:,:,:)                        !
  real(r8), allocatable ::  THETHL(:,:,:)                      !
  real(r8), allocatable ::  VLsoiAirPMA(:,:,:)                      !
  real(r8), allocatable ::  VLsoiAirPMB(:,:,:)                      !
  real(r8), allocatable ::  VLWatMicPMA_vr(:,:,:)                      !
  real(r8), allocatable ::  VLWatMicPMB_vr(:,:,:)                      !
  real(r8), allocatable ::  VLWatMicPXA_vr(:,:,:)                      !maximum non-band NH4 mass in micropore [gN /d2]
  real(r8), allocatable ::  VLWatMicPXB_vr(:,:,:)                      !maximum band NH4 mass in micropore [gN /d2]
  real(r8), allocatable ::  PARG_cef(:,:,:)                        !

  real(r8), allocatable ::  trcg_Ebu_flxM_vr(:,:,:,:)                      !
  real(r8), allocatable ::  trcg_FloXSurRunoff_flxM(:,:,:)           !gas flux through runoff at iteration M [g d-2]
  real(r8), allocatable ::  trcn_FloXSurRunoff_flxM(:,:,:)           !nutrient flux through runoff at iteration M [g d-2]
  real(r8), allocatable ::  DOM_MacP2(:,:,:,:,:)                     !copy of DOM in macropore for transport iteration [g d-2]
  real(r8), allocatable ::  DOMdiffusivity2_vr(:,:,:,:)              !

  real(r8), allocatable :: trcg_solsml2_snvr(:,:,:,:)
  real(r8), allocatable :: trcn_solsml2_snvr(:,:,:,:)

  real(r8), allocatable ::  trcg_TBLS_snvr(:,:,:,:)                      !
  real(r8), allocatable ::  trcn_TBLS_snvr(:,:,:,:)
  real(r8), allocatable ::  DOM_FloXSurRof_flxM_2DH(:,:,:,:,:,:)        !2D flow of DOM through surface runoff at iteration M [g d-2]
  real(r8), allocatable ::  trcg_FloXSurRof_flxM_2DH(:,:,:,:,:)         !2D flow of gas through surface runoff at iteration M [g d-2]
  real(r8), allocatable ::  trcn_FloXSurRof_flxM_2DH(:,:,:,:,:)         !2D flow of nutrient through surface runoff at iteration M [g d-2]

  real(r8), allocatable :: trc_gasml2_vr(:,:,:,:)         !copy of gas tracer for transport [g d-2]
  real(r8), allocatable :: trc_solml2_vr(:,:,:,:)         !copy of solute tracer in micropore for transport [g d-2]
  real(r8), allocatable :: trc_soHml2_vr(:,:,:,:)         !copy of solute tracer in macropore for transport [g d-2]

  real(r8), allocatable ::  trcn_SnowDrift_flxM_2D(:,:,:,:)                      !nutrient flux through snow drift at iteration M [g d-2]
  real(r8), allocatable ::  trcg_SnowDrift_flxM_2D(:,:,:,:)                      !gas flux through snow drift at iteration M  [g d-2]
                  !
  real(r8), allocatable ::  RGas_Disol_FlxMM_vr(:,:,:,:)         !gas dissolution (>0 into aqueous phase)
  real(r8), allocatable ::  RGasADFlxMM_3D(:,:,:,:,:)            !3D gas flux advection + diffusion
  real(r8), allocatable ::  Gas_AdvDif_FlxMM_vr(:,:,:,:)         !total 3D gas flux advection + diffusion at iteration MM [g d-2]

  real(r8), allocatable :: trcs_MacpTranspFlxM_3D(:,:,:,:,:)     !solute 3D macropore flux at iteration M [g d-2]
  real(r8), allocatable :: trcs_MicpTranspFlxM_3D(:,:,:,:,:)     !solute 3D micropore flux at iteration M [g d-2]
  real(r8), allocatable :: trcs_Mac2MicPore_flxM_vr(:,:,:,:)     !Mac to micropore nutrient flux at iteration M [g d-2]
!----------------------------------------------------------------------

contains
  subroutine InitTransfrData
  use EcoSiMParDataMod, only : micpar
  implicit none
  integer :: NumOfLitrCmplxs
  NumOfLitrCmplxs=micpar%NumOfLitrCmplxs

  allocate(DOM_FloXSurRunoff_flxM(idom_beg:idom_end,1:jcplx,JV,JH));    DOM_FloXSurRunoff_flxM=0._r8
  allocate(DOM_MicP2(idom_beg:idom_end,1:jcplx,0:JZ,JY,JX));    DOM_MicP2=0._r8
  allocate(RDOM_CumEcoProd_vr(idom_beg:idom_end,1:jcplx,0:JZ,JY,JX));  RDOM_CumEcoProd_vr=0._r8
  allocate(DOM_MicpTranspFlxM_3D(idom_beg:idom_end,1:jcplx,3,0:JD,JV,JH));DOM_MicpTranspFlxM_3D=0._r8
  allocate(DOM_MacpTranspFlxM_3D(idom_beg:idom_end,1:jcplx,3,JD,JV,JH));  DOM_MacpTranspFlxM_3D=0._r8
  allocate(DOM_Transp2Micp_flxM_vr(idom_beg:idom_end,1:jcplx,JZ,JY,JX));    DOM_Transp2Micp_flxM_vr=0._r8

  allocate(trcg_AquaAdv_flxM_snvr(idg_beg:idg_NH3,JS,JY,JX));   trcg_AquaAdv_flxM_snvr=0._r8
  allocate(trcn_AquaAdv_flxM_snvr(ids_nut_beg:ids_nuts_end,JS,JY,JX));   trcn_AquaAdv_flxM_snvr=0._r8
  allocate(DOM_Flo2LitrM(idom_beg:idom_end,1:NumOfLitrCmplxs,JY,JX));  DOM_Flo2LitrM=0._r8
  allocate(DOM_Flo2TopSoilM(idom_beg:idom_end,1:NumOfLitrCmplxs,JY,JX));  DOM_Flo2TopSoilM=0._r8
  allocate(trcg_RFL0(idg_beg:idg_NH3,JY,JX));      trcg_RFL0=0._r8
  allocate(trcn_RFL0(ids_nut_beg:ids_nuts_end,JY,JX));      trcn_RFL0=0._r8
  allocate(trcs_RFL1(ids_beg:ids_end,JY,JX));      trcs_RFL1=0._r8
  allocate(DOM_Mac2MicPore_flxM_vr(idom_beg:idom_end,1:jcplx,JZ,JY,JX));DOM_Mac2MicPore_flxM_vr=0._r8
  allocate(POSGL2(0:JZ,JY,JX)); POSGL2=0._r8
  allocate(O2AquaDiffusvity2(0:JZ,JY,JX)); O2AquaDiffusvity2=0._r8

  allocate(DOM_SurfRunoff_flxM(idom_beg:idom_end,1:jcplx,JY,JX)); DOM_SurfRunoff_flxM=0._r8

  allocate(trcg_SnowDrift_flxM(idg_beg:idg_NH3,JY,JX));      trcg_SnowDrift_flxM=0._r8
  allocate(trcn_SnowDrift_flxM(ids_nut_beg:ids_nuts_end,JY,JX)); trcn_SnowDrift_flxM=0._r8
  allocate(trcg_SurfRunoff_flxM(idg_beg:idg_NH3,JY,JX)); trcg_SurfRunoff_flxM=0._r8
  allocate(trcn_SurfRunoff_flxM(ids_nut_beg:ids_nuts_end,JY,JX));trcn_SurfRunoff_flxM=0._r8


  allocate(GasDifctScaled_vr(idg_beg:idg_end,JZ,JY,JX)); GasDifctScaled_vr=0._r8
  allocate(SoluteDifusivitytscaled_vr(ids_beg:ids_end,0:JZ,JY,JX));SoluteDifusivitytscaled_vr=0._r8
  allocate(DifuscG_vr(idg_beg:idg_end,3,JZ,JY,JX)); DifuscG_vr=0._r8

  allocate(trcg_VLWatMicP_vr(idg_beg:idg_end,0:JZ,JY,JX)); trcg_VLWatMicP_vr=0._r8

  allocate(DH2GG(3,JZ,JY,JX));  DH2GG=0._r8
  allocate(RHGFXS(JZ,JY,JX));   RHGFXS=0._r8
  allocate(THGFLG(JZ,JY,JX));   THGFLG=0._r8
  allocate(CumReductVLsoiAirPM(JZ,JY,JX));     CumReductVLsoiAirPM=0._r8
  allocate(THETHL(JZ,JY,JX));   THETHL=0._r8
  allocate(VLsoiAirPMA(JZ,JY,JX));   VLsoiAirPMA=0._r8
  allocate(VLsoiAirPMB(JZ,JY,JX));   VLsoiAirPMB=0._r8
  allocate(VLWatMicPMA_vr(JZ,JY,JX));   VLWatMicPMA_vr=0._r8
  allocate(VLWatMicPMB_vr(JZ,JY,JX));   VLWatMicPMB_vr=0._r8
  allocate(VLWatMicPXA_vr(0:JZ,JY,JX)); VLWatMicPXA_vr=0._r8
  allocate(VLWatMicPXB_vr(JZ,JY,JX));   VLWatMicPXB_vr=0._r8
  allocate(PARG_cef(idg_beg:idg_NH3,JY,JX));      PARG_cef=0._r8

  allocate(RBGCSinkGasMM_vr(idg_beg:idg_end,0:JZ,JY,JX));RBGCSinkGasMM_vr=0._r8
  allocate(RBGCSinkSoluteM_vr(ids_nuts_beg:ids_nuts_end,0:JZ,JY,JX));RBGCSinkSoluteM_vr=0._r8

  allocate(trcg_Ebu_flxM_vr(idg_beg:idg_end,JZ,JY,JX));   trcg_Ebu_flxM_vr                        = 0._r8
  allocate(trcg_FloXSurRunoff_flxM(idg_beg:idg_end,JV,JH));     trcg_FloXSurRunoff_flxM           = 0._r8
  allocate(trcn_FloXSurRunoff_flxM(ids_nut_beg:ids_nuts_end,JV,JH));     trcn_FloXSurRunoff_flxM  = 0._r8
  allocate(DOM_MacP2(idom_beg:idom_end,1:jcplx,JZ,JY,JX));DOM_MacP2                     = 0._r8
  allocate(DOM_Transp2Macp_flxM_vr(idom_beg:idom_end,1:jcplx,JZ,JY,JX));DOM_Transp2Macp_flxM_vr = 0._r8

  allocate(DOMdiffusivity2_vr(idom_beg:idom_end,0:JZ,JY,JX)); DOMdiffusivity2_vr = 0._r8

  allocate(trcg_solsml2_snvr(idg_beg:idg_NH3,JS,JY,JX));trcg_solsml2_snvr = 0._r8
  allocate(trcn_solsml2_snvr(ids_nut_beg:ids_nuts_end,JS,JY,JX));trcn_solsml2_snvr  = 0._r8

  allocate(trcg_TBLS_snvr(idg_beg:idg_NH3,JS,JY,JX));   trcg_TBLS_snvr = 0._r8
  allocate(trcn_TBLS_snvr(ids_nut_beg:ids_nuts_end,JS,JY,JX));trcn_TBLS_snvr     = 0._r8

  allocate(DOM_FloXSurRof_flxM_2DH(idom_beg:idom_end,1:jcplx,2,2,JV,JH));DOM_FloXSurRof_flxM_2DH = 0._r8
  allocate(trcg_FloXSurRof_flxM_2DH(idg_beg:idg_NH3,2,2,JV,JH));  trcg_FloXSurRof_flxM_2DH       = 0._r8

  allocate(trc_gasml2_vr(idg_beg:idg_end,0:JZ,JY,JX)); trc_gasml2_vr(:,:,:,:) = 0._r8
  allocate(trc_solml2_vr(ids_beg:ids_end,0:JZ,JY,JX)); trc_solml2_vr(:,:,:,:) = 0._r8
  allocate(trc_soHml2_vr(ids_beg:ids_end,0:JZ,JY,JX)); trc_soHml2_vr(:,:,:,:) = 0._r8

  allocate(trcg_SnowDrift_flxM_2D(idg_beg:idg_NH3,2,JV,JH));    trcg_SnowDrift_flxM_2D                    = 0._r8
  allocate(trcn_FloXSurRof_flxM_2DH(ids_nut_beg:ids_nuts_end,2,2,JV,JH));  trcn_FloXSurRof_flxM_2DH = 0._r8
  allocate(RGasADFlxMM_3D(idg_beg:idg_NH3,3,JD,JV,JH));RGasADFlxMM_3D                             = 0._r8
  allocate(Gas_AdvDif_FlxMM_vr(idg_beg:idg_NH3,JZ,JY,JH));Gas_AdvDif_FlxMM_vr                     = 0._r8

  allocate(RGas_Disol_FlxMM_vr(idg_beg:idg_end,0:JZ,JY,JX)); RGas_Disol_FlxMM_vr = 0._r8

  allocate(trcn_SnowDrift_flxM_2D(ids_nut_beg:ids_nuts_end,2,JV,JH)); trcn_SnowDrift_flxM_2D = 0._r8


  allocate(trcs_MacpTranspFlxM_3D(ids_beg:ids_end,3,JD,JV,JH));trcs_MacpTranspFlxM_3D       = 0._r8
  allocate(trcs_MicpTranspFlxM_3D(ids_beg:ids_end,3,0:JD,JV,JH));trcs_MicpTranspFlxM_3D     = 0._r8
  allocate(trcs_Mac2MicPore_flxM_vr(ids_beg:ids_end,JZ,JY,JX));trcs_Mac2MicPore_flxM_vr     = 0._r8
  allocate(trcs_Transp2Micp_flxM_vr(ids_beg:ids_end,JZ,JY,JX));trcs_Transp2Micp_flxM_vr = 0._r8
  allocate(trcs_Transp2Macp_flxM_vr(ids_beg:ids_end,JZ,JY,JX));trcs_Transp2Macp_flxM_vr = 0._r8
  allocate(trcg_AquaADV_Snow2Soil_flxM(idg_beg:idg_end));trcg_AquaADV_Snow2Soil_flxM                      = 0._r8
  allocate(trcn_AquaADV_Snow2Band_flxM(ids_nutb_beg:ids_nutb_end));trcn_AquaADV_Snow2Band_flxM  = 0._r8
  allocate(trcn_AquaADV_Snow2Soil_flxM(ids_nut_beg:ids_nuts_end));trcn_AquaADV_Snow2Soil_flxM   = 0._r8

  end subroutine InitTransfrData

!----------------------------------------------------------------------
  subroutine DestructTransfrData
  use abortutils, only : destroy
  implicit none

  call destroy(GasDifctScaled_vr)
  call destroy(DifuscG_vr)
  call destroy(SoluteDifusivitytscaled_vr)
  call destroy(DOM_MicP2)

  call destroy(RDOM_CumEcoProd_vr)
  call destroy(DOM_MicpTranspFlxM_3D)
  call destroy(DOM_MacpTranspFlxM_3D)
  call destroy(RBGCSinkSoluteM_vr)
  call destroy(DOM_FloXSurRunoff_flxM)
  call destroy(DOM_FloXSurRof_flxM_2DH)
  call destroy(DOM_Flo2LitrM)
  call destroy(DOM_Flo2TopSoilM)
  call destroy(trcg_AquaADV_Snow2Soil_flxM)
  call destroy(trcs_RFL1)
  call destroy(trcn_RFL0)
  call destroy(DOM_Mac2MicPore_flxM_vr)
  call destroy(trcg_AquaAdv_flxM_snvr)
  call destroy(trcn_AquaAdv_flxM_snvr)
  call destroy(trcg_Ebu_flxM_vr)
  call destroy(RBGCSinkGasMM_vr)
  call destroy(POSGL2)
  call destroy(O2AquaDiffusvity2)
  call destroy(trcg_FloXSurRunoff_flxM)
  call destroy(DOM_SurfRunoff_flxM)
  call destroy(trcg_SurfRunoff_flxM)
  call destroy(trcn_SurfRunoff_flxM)
  call destroy(trcg_SnowDrift_flxM)
  call destroy(trcn_SnowDrift_flxM)
  call destroy(trcs_Transp2Micp_flxM_vr)
  call destroy(trcg_VLWatMicP_vr)
  call destroy(trcn_FloXSurRunoff_flxM)
  call destroy(DH2GG)
  call destroy(RHGFXS)
  call destroy(THGFLG)
  call destroy(CumReductVLsoiAirPM)
  call destroy(THETHL)
  call destroy(VLsoiAirPMA)
  call destroy(VLsoiAirPMB)
  call destroy(VLWatMicPMA_vr)
  call destroy(VLWatMicPMB_vr)
  call destroy(VLWatMicPXA_vr)
  call destroy(VLWatMicPXB_vr)
  call destroy(PARG_cef)
  call destroy(DOMdiffusivity2_vr)
  call destroy(RGasADFlxMM_3D)
  call destroy(trcg_solsml2_snvr)
  call destroy(trcn_solsml2_snvr)
  call destroy(trcn_FloXSurRof_flxM_2DH)
  call destroy(trcg_TBLS_snvr)
  call destroy(trcn_TBLS_snvr)
  call destroy(trcg_RFL0)
  call destroy(trcg_SnowDrift_flxM_2D)
  call destroy(trcs_Transp2Macp_flxM_vr)

  call destroy(DOM_MacP2)
  call destroy(Gas_AdvDif_FlxMM_vr)
  call destroy(trc_gasml2_vr)
  call destroy(trc_solml2_vr)
  call destroy(trc_soHml2_vr)
  call destroy(trcn_AquaADV_Snow2Band_flxM)
  call destroy(trcn_AquaADV_Snow2Soil_flxM)
  call destroy(trcg_FloXSurRof_flxM_2DH)
  call destroy(RGas_Disol_FlxMM_vr)
  call destroy(trcn_SnowDrift_flxM_2D)
  call destroy(trcs_MacpTranspFlxM_3D)
  call destroy(trcs_MicpTranspFlxM_3D)
  call destroy(trcs_Mac2MicPore_flxM_vr)

  call destroy(DOM_Transp2Macp_flxM_vr)
  call destroy(DOM_Transp2Micp_flxM_vr)
  end subroutine DestructTransfrData

end module TranspNoSaltDataMod
