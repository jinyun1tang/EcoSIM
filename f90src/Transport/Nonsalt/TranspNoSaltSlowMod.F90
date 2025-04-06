module TranspNoSaltSlowMod
  use data_kind_mod,     only: r8 => DAT_KIND_R8
  use abortutils,        only: destroy,  endrun
  use TracerPropMod,     only: MolecularWeight
  use PlantDataRateType, only: trcs_deadroot2soil_vr
  use minimathmod    
  use IrrigationDataType
  use DebugToolMod
  use GridConsts
  use EcoSIMSolverPar
  use TranspNoSaltDataMod
  use AqueChemDatatype
  use GridDataType
  use SOMDataType
  use SoilBGCDataType
  use SurfSoilDataType
  use SurfLitterDataType
  use SnowDataType
  use EcoSIMCtrlDataType
  use SoilWaterDataType
  use LandSurfDataType
  use SoilPropertyDataType
  use ChemTranspDataType
  use ClimForcDataType
  USE EcoSimConst
  USE SoilHeatDataType
  use TracerIDMod
implicit none

  private
  CHARACTER(LEN=*), PARAMETER :: MOD_FILENAME=&
  __FILE__

  public :: TransptSlowNoSaltM
  public :: BubbleEffluxM

  contains
!------------------------------------------------------------------------------------------
  subroutine TransptSlowNoSaltM(I,J,M,NHE,NHW,NVS,NVN)
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: M
  integer, intent(in) :: NHE,NHW,NVS,NVN
  character(len=*), parameter :: subname='TransptSlowNoSaltM'

  real(r8) :: pscal(ids_beg:ids_end), dpscal(ids_beg:ids_end)
  real(r8) :: pscal1_dom(idom_beg:idom_end),dpscal_dom(idom_beg:idom_end)
  real(r8) :: dpscal_max
  integer :: ids,idom

  call PrintInfo('beg '//subname)
  dpscal=1._r8;dpscal_dom=1._r8

  call StageSlowTranspIterationM(I,J,M,NHE,NHW,NVS,NVN)

  call SnowSoluteVertDrainM(I,J,M,NHE,NHW,NVS,NVN)
  dpscal_max=1._r8
  DO while(dpscal_max>1.e-2_r8)
    pscal=1.001_r8
    call ZeroFluxesM(I,J,M,NHE,NHW,NVS,NVN)  

    call SoluteFluxSurfaceM(I,J,M,NHE,NHW,NVS,NVN)

    call TracerXBoundariesM(I,J,M,NHW,NHE,NVN,NVS)

    call TracerXGridsM(I,J,M,NHE,NHW,NVS,NVN)

    call SlowUpdateStateVars(I,J,M,NHW,NHE,NVN,NVS,dpscal,dpscal_dom,pscal,pscal1_dom)

    call AccumSlowFluxesM(I,J,M,NHE,NHW,NVS,NVN,dpscal,dpscal_dom,pscal,pscal1_dom)

    call copyStateVars(NHE,NHW,NVS,NVN)

    dpscal_max=0._r8
    DO ids=ids_beg,ids_end
      dpscal(ids)=dpscal(ids)*(1._r8-pscal(ids))
      dpscal_max=AMAX1(dpscal_max,dpscal(ids))
    enddo

    DO idom=idom_beg,idom_end
      dpscal_dom(idom)=dpscal_dom(idom)*(1._r8-pscal1_dom(idom))
      dpscal_max=AMAX1(dpscal_max,dpscal_dom(idom))
    ENDDO
  ENDDO

  call PrintInfo('end '//subname)
  end subroutine TransptSlowNoSaltM
!------------------------------------------------------------------------------------------
  subroutine SlowUpdateStateVars(I,J,M,NHW,NHE,NVN,NVS,dpscal,dpscal_dom,pscal1,pscal1_dom)
  implicit none
  integer, intent(in) :: I,J,M,NHW,NHE,NVN,NVS
  real(r8),intent(in) :: dpscal(ids_beg:ids_end)
  real(r8),intent(in) :: dpscal_dom(idom_beg:idom_end)
  real(r8),optional, intent(inout) :: pscal1(ids_beg:ids_end)
  real(r8),optional, intent(inout) :: pscal1_dom(idom_beg:idom_end)

  if(present(pscal1) .and. present(pscal1_dom))then
    pscal1     = 1.001_r8
    pscal1_dom = 1.001_r8
    call UpdateSnowTracersM(I,J,M,NHW,NHE,NVN,NVS,dpscal, pscal1)

    call UpdateSurfTracerM(I,J,M,NHW,NHE,NVN,NVS,dpscal,dpscal_dom,pscal1,pscal1_dom)

    call UpdateSubSurfTracerM(I,J,M,NHW,NHE,NVN,NVS,dpscal,dpscal_dom,pscal1,pscal1_dom)
    pscal1     = pscal1*0.9999_r8
    pscal1_dom = pscal1_dom*0.9999_r8
  else
    call UpdateSnowTracersM(I,J,M,NHW,NHE,NVN,NVS,dpscal)

    call UpdateSurfTracerM(I,J,M,NHW,NHE,NVN,NVS,dpscal,dpscal_dom)

    call UpdateSubSurfTracerM(I,J,M,NHW,NHE,NVN,NVS,dpscal,dpscal_dom)
  endif
  end subroutine SlowUpdateStateVars
!------------------------------------------------------------------------------------------

  subroutine copyStateVars(NHE,NHW,NVS,NVN)
  implicit none
  integer, intent(in) :: NHE,NHW,NVS,NVN
  character(len=*), parameter :: subname='copyStateVars'
  integer :: NY,NX,idg,ids,L,K,idom
  
  call PrintInfo('beg '//subname)
  DO NY=NHW,NHE
    DO  NX=NVN,NVS
      !to litter
      DO  K=1,jcplx
        DO idom=idom_beg,idom_end
          DOM_MicP2_vr(idom,K,0,NY,NX)=DOM_MicP_vr(idom,K,0,NY,NX)
        ENDDO
      ENDDO

      !exclude NH3B
      DO idg=idg_beg,idg_NH3
        trcs_solml2_vr(idg,0,NY,NX)=trcs_solml_vr(idg,0,NY,NX)
      ENDDO
      
      DO ids=ids_nut_beg,ids_nuts_end
        trcs_solml2_vr(ids,0,NY,NX)=trcs_solml_vr(ids,0,NY,NX)
      ENDDO

      DO L=NU(NY,NX),NL(NY,NX)    
        DO ids=ids_beg,ids_end
          trcs_solml2_vr(ids,L,NY,NX) = trcs_solml_vr(ids,L,NY,NX)
          trcs_soHml2_vr(ids,L,NY,NX) = trcs_soHml_vr(ids,L,NY,NX)
        ENDDO

        DO  K=1,jcplx
          DO idom=idom_beg,idom_end
            DOM_MicP2_vr(idom,K,L,NY,NX)=DOM_MicP_vr(idom,K,L,NY,NX)
            DOM_MacP2_vr(idom,K,L,NY,NX)=DOM_MacP_vr(idom,K,L,NY,NX)
          ENDDO
        ENDDO
      ENDDO

      DO L=1,JS
        trcg_solsml2_snvr(idg_beg:idg_NH3,L,NY,NX)          = trcg_solsml_snvr(idg_beg:idg_NH3,L,NY,NX)
        trcn_solsml2_snvr(ids_nut_beg:ids_nuts_end,L,NY,NX) = trcn_solsml_snvr(ids_nut_beg:ids_nuts_end,L,NY,NX)
      enddo
    ENDDO
  ENDDO    
  call PrintInfo('end '//subname)
  end subroutine copyStateVars
!------------------------------------------------------------------------------------------

  subroutine UpdateSubSurfTracerM(I,J,M,NHW,NHE,NVN,NVS,dpscal,dpscal_dom,pscal1,pscal1_dom)
  implicit none
  integer, intent(in) :: I,J,M,NHW,NHE,NVN,NVS
  real(r8), intent(in):: dpscal(ids_beg:ids_end)
  real(r8), intent(in):: dpscal_dom(idom_beg:idom_end)
  real(r8),optional,intent(inout) :: pscal1(ids_beg:ids_end)
  real(r8),optional,intent(inout) :: pscal1_dom(idom_beg:idom_end)

  character(len=*), parameter :: subname = 'UpdateSubSurfTracerM'
  integer :: NY,NX,L,K,idg,ids,idom
  real(r8) :: pscal(ids_beg:ids_end)
  real(r8) :: pscal_dom(idom_beg:idom_end),flux

  call PrintInfo('beg '//subname)
  if(present(pscal1))then
    pscal=pscal1
  else
    pscal=1.001_r8
  endif
  if(present(pscal1_dom))then
    pscal_dom=pscal1_dom
  else
    pscal_dom=1.001_r8
  endif

  DO NY=NHW,NHE
    DO  NX=NVN,NVS
      DO L=NU(NY,NX)+1,NL(NY,NX)    
        IF(VLSoilPoreMicP_vr(L,NY,NX).GT.ZEROS2(NY,NX))THEN        

          DO idg=idg_beg,idg_end     
            if(dpscal(idg)>tiny_p)then
              flux=trcsol_Irrig_flxM_vr(ids,L,NY,NX)+trcs_Transp2Micp_flxM_vr(ids,L,NY,NX) + &
                trcs_Mac2MicPore_flxM_vr(ids,L,NY,NX)      
              flux=flux*dpscal(idg)        
              call get_flux_scalar(trcs_solml2_vr(idg,L,NY,NX),flux,trcs_solml_vr(idg,L,NY,NX),pscal(idg))
            endif
          ENDDO

          DO ids=idg_NH3B+1,ids_end
            if(dpscal(ids)>tiny_p)then
              flux=trcsol_Irrig_flxM_vr(ids,L,NY,NX)+trcs_Transp2Micp_flxM_vr(ids,L,NY,NX) + &
                trcs_Mac2MicPore_flxM_vr(ids,L,NY,NX)
              flux=flux*dpscal(ids)  
              call get_flux_scalar(trcs_solml2_vr(ids,L,NY,NX),flux,trcs_solml_vr(ids,L,NY,NX),pscal(ids))            
            endif
          ENDDO

          DO  K=1,jcplx
            DO idom=idom_beg,idom_end
              if(dpscal_dom(idom)>tiny_p)then
                flux=DOM_Transp2Micp_flxM_vr(idom,K,L,NY,NX) + DOM_Mac2MicPore_flxM_vr(idom,K,L,NY,NX)
                flux=flux*dpscal_dom(idom)
                call get_flux_scalar(DOM_MicP2_vr(idom,K,L,NY,NX), flux, DOM_MicP_vr(idom,K,L,NY,NX),pscal_dom(idom))

                flux=DOM_Transp2Macp_flxM_vr(idom,K,L,NY,NX) - DOM_Mac2MicPore_flxM_vr(idom,K,L,NY,NX)
                flux=flux*dpscal_dom(idom)
                call get_flux_scalar(DOM_MacP2_vr(idom,K,L,NY,NX),flux, DOM_MacP_vr(idom,K,L,NY,NX),pscal_dom(idom))
              endif
            ENDDO
          ENDDO

          DO ids=ids_beg,ids_end
            if(dpscal(ids)>tiny_p)then
              flux=trcs_Transp2Macp_flxM_vr(ids,L,NY,NX) - trcs_Mac2MicPore_flxM_vr(ids,L,NY,NX)
              flux=flux*dpscal(ids)
              call get_flux_scalar(trcs_soHml2_vr(ids,L,NY,NX),flux, trcs_soHml_vr(ids,L,NY,NX),pscal(ids))            
            endif
          enddo
        ENDIF
      ENDDO  
    ENDDO
  ENDDO  
  if(present(pscal1))pscal1=pscal
  if(present(pscal1_dom))pscal1_dom=pscal_dom
  call PrintInfo('end '//subname)
  end subroutine UpdateSubSurfTracerM
!------------------------------------------------------------------------------------------
  subroutine AccumSlowFluxesM(I,J,M,NHE,NHW,NVS,NVN,dpscal,dpcsal_dom,pscal1,pscal1_dom)
  !
  !DESCRIPTION
  !Accumulate fluxes for mass conservation check
  !Lateral runoff, dissolution (aka dry gas deposition), 
  !wet deposition through rain and irrigation, 
  !snow fluxes to litter and soil
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: M
  integer, intent(in) :: NHE,NHW,NVS,NVN
  real(r8),intent(in) :: dpscal(ids_beg:ids_end)
  real(r8), intent(in):: pscal1(ids_beg:ids_end)
  real(r8), intent(in):: dpcsal_dom(idom_beg:idom_end)
  real(r8), intent(in):: pscal1_dom(idom_beg:idom_end)

  character(len=*), parameter :: subname='AccumSlowFluxesM'
  integer :: idom,K,ids,NY,NX,idg,idn,L
  real(r8):: ppscal(ids_beg:ids_end),pscal_max
  real(r8):: ppscal_dom(idom_beg:idom_end)

  call PrintInfo('beg '//subname)
  pscal_max=0._r8
  DO ids=ids_beg,ids_end
    ppscal(ids) = dpscal(ids)*pscal1(ids)
    pscal_max=AMAX1(pscal1(ids),pscal_max)
  ENDDO

  do idom=idom_beg,idom_end    
    ppscal_dom(idom)=dpcsal_dom(idom)*pscal1_dom(idom)
    pscal_max=AMAX1(pscal1_dom(idom),pscal_max)
  enddo

  if(pscal_max<0.999_r8)then
    call SlowUpdateStateVars(I,J,M,NHW,NHE,NVN,NVS,ppscal,ppscal_dom)  
  endif

  DO NX=NHW,NHE
    DO  NY=NVN,NVS
      DO idg=idg_beg,idg_NH3      
        GasDiff2Surf_flx_col(idg,NY,NX) = GasDiff2Surf_flx_col(idg,NY,NX)+ppscal(idg) &
          *(RGasAtmDisol2LitrM_col(idg,NY,NX)+RGasAtmDisol2SoilM_col(idg,NY,NX))
      ENDDO

      RCH4PhysexchPrev_vr(0,NY,NX)= RCH4PhysexchPrev_vr(0,NY,NX) + ppscal(idg_CH4)*(RGasAtmDisol2LitrM_col(idg_CH4,NY,NX) &
        -trcg_Precip2LitrM(idg_CH4,NY,NX))

      RO2AquaSourcePrev_vr(0,NY,NX)   =RO2AquaSourcePrev_vr(0,NY,NX)+ppscal(idg_O2)*(RGasAtmDisol2LitrM_col(idg_O2,NY,NX) &
        -trcg_Precip2LitrM(idg_O2,NY,NX)) 

      RCH4PhysexchPrev_vr(NU(NY,NX),NY,NX)= RCH4PhysexchPrev_vr(NU(NY,NX),NY,NX)+ppscal(idg_CH4)*(RGasAtmDisol2SoilM_col(idg_CH4,NY,NX) &
        - trcs_Precip2MicpM(idg_CH4,NY,NX))

      RO2AquaSourcePrev_vr(NU(NY,NX),NY,NX)   =RO2AquaSourcePrev_vr(NU(NY,NX),NY,NX)+ppscal(idg_O2)*(RGasAtmDisol2SoilM_col(idg_O2,NY,NX) &
        - trcs_Precip2MicpM(idg_O2,NY,NX))

      idg=idg_NH3B
      GasDiff2Surf_flx_col(idg,NY,NX) = GasDiff2Surf_flx_col(idg,NY,NX)+RGasAtmDisol2SoilM_col(idg,NY,NX)*ppscal(idg)              
      
      DO idg=idg_beg,idg_NH3        
        GasHydroLossFlx_col(idg,NY,NX)=GasHydroLossFlx_col(idg,NY,NX)+trcg_SurfRunoff_flxM(idg,NY,NX)*ppscal(idg)
      ENDDO
      
      do idg=idg_beg,idg_NH3
        trcg_AquaADV_Snow2Litr_flx(idg,NY,NX)=trcg_AquaADV_Snow2Litr_flx(idg,NY,NX)+trcg_AquaADV_Snow2Litr_flxM(idg,NY,NX)*ppscal(idg)           
      ENDDO

      do ids=ids_nut_beg,ids_nuts_end
        trcn_AquaADV_Snow2Litr_flx(ids,NY,NX)=trcn_AquaADV_Snow2Litr_flx(ids,NY,NX)+trcn_AquaADV_Snow2Litr_flxM(ids,NY,NX)*ppscal(ids)
      enddo

      DO  K=1,jcplx
        do idom=idom_beg,idom_end
          DOM_transpFlx_2DH(idom,K,NY,NX) = DOM_transpFlx_2DH(idom,K,NY,NX)+DOM_transpFlxM_2DH(idom,K,NY,NX)*ppscal_dom(idom)

          DOM_draing_col(idom,NY,NX)    = DOM_draing_col(idom,NY,NX)+ ppscal_dom(idom) &
            * (DOM_MicpTranspFlxM_3D(idom,K,3,NL(NY,NX),NY,NX)+DOM_MacpTranspFlxM_3D(idom,K,3,NL(NY,NX),NY,NX))
        ENDDO
      ENDDO

      DO ids=ids_beg,ids_end
        trcs_transpFlx_2DH(ids,NY,NX) = trcs_transpFlx_2DH(ids,NY,NX)+trcs_transpFlxM_2DH(ids,NY,NX)*ppscal(ids)
        trcs_draing_col(ids,NY,NX)    = trcs_draing_col(ids,NY,NX) + ppscal(ids) &
          *(trcs_MicpTranspFlxM_3D(ids,3,NL(NY,NX),NY,NX)+trcs_MacpTranspFlxM_3D(ids,3,NL(NY,NX),NY,NX))
      ENDDO

      DO idg=idg_beg,idg_NH3
        trcg_SnowDrift_flx_col(idg,NY,NX)=trcg_SnowDrift_flx_col(idg,NY,NX)+trcg_SnowDrift_flxM(idg,NY,NX)*ppscal(idg)
        GasHydroLossFlx_col(idg,NY,NX)=GasHydroLossFlx_col(idg,NY,NX)+trcg_SnowDrift_flxM(idg,NY,NX)*ppscal(idg)
      ENDDO  

      DO idn=ids_nut_beg,ids_nuts_end
        trcn_SnowDrift_flx_col(idn,NY,NX) =trcn_SnowDrift_flx_col(idn,NY,NX) + trcn_SnowDrift_flxM(idn,NY,NX)*ppscal(idn)
      ENDDO  

      DO L=NU(NY,NX),NL(NY,NX)
        RO2AquaSourcePrev_vr(L,NY,NX)   =RO2AquaSourcePrev_vr(L,NY,NX) +(trcs_Transp2Micp_flxM_vr(idg_O2,L,NY,NX) + &
          trcs_Mac2MicPore_flxM_vr(idg_O2,L,NY,NX)+trcsol_Irrig_flxM_vr(idg_O2,L,NY,NX))*ppscal(idg_O2)
        
        RCH4PhysexchPrev_vr(L,NY,NX) = RCH4PhysexchPrev_vr(L,NY,NX)+(trcs_Transp2Micp_flxM_vr(idg_CH4,L,NY,NX) + &
          trcs_Mac2MicPore_flxM_vr(idg_CH4,L,NY,NX)+trcsol_Irrig_flxM_vr(idg_CH4,L,NY,NX))*ppscal(idg_CH4)
      ENDDO
    ENDDO
  ENDDO
  call PrintInfo('end '//subname)
  end subroutine AccumSlowFluxesM
!------------------------------------------------------------------------------------------
  subroutine StageSlowTranspIterationM(I,J,M,NHE,NHW,NVS,NVN)
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: M
  integer, intent(in) :: NHE,NHW,NVS,NVN
  character(len=*), parameter :: subname='StageSlowTranspIterationM'
  integer :: NY,NX,L, K

  call PrintInfo('beg '//subname)

  DO NX=NHW,NHE
    DO  NY=NVN,NVS
      DO L=NU(NY,NX),NL(NY,NX)
        VLWatMicPMA_vr(L,NY,NX) = VLWatMicPM_vr(M,L,NY,NX)*trcs_VLN_vr(ids_NH4,L,NY,NX)
        VLWatMicPMB_vr(L,NY,NX) = VLWatMicPM_vr(M,L,NY,NX)*trcs_VLN_vr(ids_NH4B,L,NY,NX)
        VLWatMicPXA_vr(L,NY,NX) = natomw*VLWatMicPMA_vr(L,NY,NX)
        VLWatMicPXB_vr(L,NY,NX) = natomw*VLWatMicPMB_vr(L,NY,NX)
        VLsoiAirPMA_vr(L,NY,NX) = VLsoiAirPM_vr(M,L,NY,NX)*trcs_VLN_vr(ids_NH4,L,NY,NX)
        VLsoiAirPMB_vr(L,NY,NX) = VLsoiAirPM_vr(M,L,NY,NX)*trcs_VLN_vr(ids_NH4B,L,NY,NX)
      ENDDO
    ENDDO  
  ENDDO
  call PrintInfo('end '//subname)
  end subroutine StageSlowTranspIterationM
!------------------------------------------------------------------------------------------
  subroutine ZeroFluxesM(I,J,M,NHE,NHW,NVS,NVN)  

  implicit none
  integer, intent(in) :: I,J,M,NHE,NHW,NVS,NVN

  integer :: L,NY,NX,K

  DOM_MacpTranspFlxM_3D   = 0._r8
  DOM_MicpTranspFlxM_3D   = 0._r8
  DOM_Mac2MicPore_flxM_vr = 0._r8

  trcs_MacpTranspFlxM_3D   = 0._r8
  trcs_MicpTranspFlxM_3D   = 0._r8
  trcs_Mac2MicPore_flxM_vr = 0._r8

  DOM_transpFlxM_2DH  = 0._r8
  trcs_transpFlxM_2DH = 0._r8

  DO NX=NHW,NHE
    DO  NY=NVN,NVS
      DO  K=1,jcplx
        DOM_SurfRunoff_flxM(idom_beg:idom_end,K,NY,NX)=0.0_r8
      ENDDO

      trcg_SurfRunoff_flxM(idg_beg:idg_NH3,NY,NX)          = 0.0_r8
      trcn_SurfRunoff_flxM(ids_nut_beg:ids_nuts_end,NY,NX) = 0.0_r8
      trcg_SnowDrift_flxM(idg_beg:idg_NH3,NY,NX)          = 0.0_r8
      trcn_SnowDrift_flxM(ids_nut_beg:ids_nuts_end,NY,NX) = 0.0_r8

      !   INITIALIZE SNOWPACK NET FLUX ACCUMULATORS
      DO  L=1,JS
        trcg_Aqua_flxM_snvr(idg_beg:idg_NH3,L,NY,NX)          = 0._r8
        trcn_Aqua_flxM_snvr(ids_nut_beg:ids_nuts_end,L,NY,NX) = 0._r8
      ENDDO

      DO L=NU(NY,NX),NL(NY,NX)    
        DO  K=1,jcplx
          DOM_Transp2Micp_flxM_vr(idom_beg:idom_end,K,L,NY,NX) = 0.0_r8
          DOM_Transp2Macp_flxM_vr(idom_beg:idom_end,K,L,NY,NX)= 0.0_r8
        ENDDO
        trcs_Transp2Micp_flxM_vr(ids_beg:ids_end,L,NY,NX) = 0.0_r8
        trcs_Transp2Macp_flxM_vr(ids_beg:ids_end,L,NY,NX) = 0._r8
      ENDDO  
    ENDDO
  ENDDO
  end subroutine ZeroFluxesM
!------------------------------------------------------------------------------------------
  subroutine TracerXGridsM(I,J,M,NHE,NHW,NVS,NVN)
  !
  !Description
  !Tracer exchange across the internal grids
  implicit none
  integer, intent(in) :: I,J,M
  integer, intent(in) :: NHE,NHW,NVS,NVN  
  character(len=*), parameter :: subname='TracerXGridsM'  
  integer :: NY,NX,L,N
  integer :: N3,N2,N1,LL
  integer :: N6,N5,N4
  LOGICAL :: LBOUND

  call PrintInfo('beg '//subname)
  DO NX=NHW,NHE
    DO  NY=NVN,NVS
      D125: DO L=1,NL(NY,NX)       
        N1=NX;N2=NY;N3=L    !source grid, avoid edges

        IF(N3.GE.NUM(N2,N1))call MicMacPoresSoluteExchM(M,N1,N2,N3)

        !identify the dest grid, which can be at the edges.                      
        D120: DO N=FlowDirIndicator_col(N2,N1),3
          IF(N.EQ.iWestEastDirection)THEN            
            IF(NX.EQ.NHE)CYCLE  !skip eastern boundary
            N4 = NX+1;N5 = NY;N6 = L
          ELSEIF(N.EQ.iNorthSouthDirection)THEN                        
            IF(NY.EQ.NVS)CYCLE  !skip southern boundary
            N4 = NX;N5 = NY+1;N6 = L
          ELSEIF(N.EQ.iVerticalDirection)THEN                        
            IF(L.EQ.NL(NY,NX))CYCLE  !skip bottom boundary
            N4 = NX;N5 = NY;N6 = L+1
          ENDIF
          
          DO LL=N6,NL(NY,NX)
            IF(VLSoilPoreMicP_vr(LL,N5,N4).GT.ZEROS2(N5,N4))THEN
              N6=LL
              exit
            ENDIF
          ENDDO
          
          IF(VLSoilPoreMicP_vr(N3,N2,N1).GT.ZEROS2(N2,N1))THEN
            IF(N3.GE.NUM(N2,N1) .AND. N6.GE.NUM(N5,N4) .AND. N3.LE.NL(N2,N1) .AND. N6.LE.NL(N5,N4))THEN
              !
              call SoluteAdvDifusTranspM(I,J,M,N,N1,N2,N3,N4,N5,N6)
              !
              IF(FlowDirIndicator_col(N2,N1).NE.3 .OR. N.EQ.iVerticalDirection)THEN
                call NetTracerFlowXSoilPoresM(NY,NX,N,M,N1,N2,N3,N4,N5,N6)
              ENDIF
            ENDIF                
          ENDIF
        ENDDO D120  
      ENDDO D125
    ENDDO
  ENDDO  
  call PrintInfo('end '//subname)    
  end subroutine TracerXGridsM

!------------------------------------------------------------------------------------------

  subroutine SoluteFluxSurfaceM(I,J,M,NHE,NHW,NVS,NVN)
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: M
  integer, intent(in) :: NHE,NHW,NVS,NVN
  character(len=*), parameter :: subname='SoluteFluxSurfaceM'

  integer  :: NY,NX,K
  real(r8) :: FLWRM1
  real(r8) :: trcs_cl_litr(ids_beg:ids_end)
  real(r8) :: trcs_cl_soil(ids_beg:ids_end)
  real(r8) :: solute_adv_Lit2Soil_flxM(ids_beg:ids_end)
  real(r8) :: trcs_Dif_Litr2Soil_flxM(ids_beg:ids_end)
  real(r8) :: CDOM_MicP_litr(idom_beg:idom_end,1:jcplx)
  real(r8) :: CDOM_MicP_soil(idom_beg:idom_end,1:jcplx)
  real(r8) :: DOM_Adv_Litr2Soil_flxM(idom_beg:idom_end,1:jcplx)
  real(r8) :: DOM_Difus_Litr2Soil_flxM(idom_beg:idom_end,1:jcplx)


  call PrintInfo('beg '//subname)

  DO NX=NHW,NHE
    DO  NY=NVN,NVS
      !-----------------------------------------------------------------
      !Vertical 1D exchange
      call LitterAtmosExchangeM(I,J,M,NY,NX,CDOM_MicP_litr,trcs_cl_litr)
      !
      call SoilAtmosExchangeM(I,J,M,NY,NX,CDOM_MicP_soil,trcs_cl_soil)
      !
      call Lit2SoilTracerAdvectM(I,J,M,NY,NX,FLWRM1,solute_adv_Lit2Soil_flxM,DOM_Adv_Litr2Soil_flxM)
      !
      call LitterSoilTracerDiffusionM(I,J,M,NY,NX,CDOM_MicP_litr,CDOM_MicP_soil,FLWRM1,trcs_cl_litr,trcs_cl_soil,&
        trcs_Dif_Litr2Soil_flxM,DOM_Difus_Litr2Soil_flxM)

      !-----------------------------------------------------------------        
      !
      call SurfLayerNetTracerFluxM(I,J,M,NY,NX,trcs_Dif_Litr2Soil_flxM,solute_adv_Lit2Soil_flxM,&
        DOM_Adv_Litr2Soil_flxM,DOM_Difus_Litr2Soil_flxM)

      !
      !     SOLUTE TRANSPORT FROM WATER OVERLAND FLOW
      !     IN 'WATSUB' AND FROM SOLUTE CONCENTRATIONS
      !     IN SOIL SURFACE LAYER
      !
      call SurfXGridTracerTransptM(M,NY,NX,NHE,NHW,NVS,NVN)
    ENDDO
  ENDDO
  call PrintInfo('end '//subname)
  end subroutine SoluteFluxSurfaceM

!------------------------------------------------------------------------------------------

  subroutine LitterAtmosExchangeM(I,J,M,NY,NX,CDOM_MicP_litr,trcs_cl_litr)
  implicit none
  integer, intent(in) :: I,J,NY,NX,M
  real(r8),intent(out) :: trcs_cl_litr(ids_beg:ids_end)
  real(r8),intent(out) :: CDOM_MicP_litr(idom_beg:idom_end,1:jcplx)

  character(len=*), parameter :: subname = 'LitterAtmosExchangeM'
  integer :: K,ids,idg,idom
  real(r8) :: trc_gasq(idg_beg:idg_NH3)  !equilibrium tracer concentation
  real(r8) :: DFGcc(idg_beg:idg_NH3)
  real(r8) :: DLYR0,TORT0

  call PrintInfo('beg '//subname)

  IF(VGeomLayer_vr(0,NY,NX).GT.ZEROS2(NY,NX) .AND. VLWatMicPM_vr(M,0,NY,NX).GT.ZEROS2(NY,NX))THEN

    DLYR0 = AMAX1(ZERO2,DLYR_3D(3,0,NY,NX)) !vertical layer thickness
    TORT0 = TortMicPM_vr(M,0,NY,NX)*AREA(3,NU(NY,NX),NY,NX)/(0.5_r8*DLYR0)*FracSurfByLitR_col(NY,NX)

    DO idg=idg_beg,idg_NH3
      DFGcc(idg)=SoluteDifusivitytscaledM_vr(idg,0,NY,NX)*TORT0
    ENDDO

    DO  K=1,jcplx
      do idom=idom_beg,idom_end
        CDOM_MicP_litr(idom,K)=AZMAX1(DOM_MicP2_vr(idom,K,0,NY,NX)/VLWatMicPM_vr(M,0,NY,NX))
      enddo
    ENDDO

    !temporary solute concentration
    DO ids=ids_beg,ids_end
      trcs_cl_litr(ids)=AZMAX1(trcs_solml2_vr(ids,0,NY,NX)/VLWatMicPM_vr(M,0,NY,NX))
    ENDDO
!
!     EQUILIBRIUM CONCENTRATIONS AT RESIDUE SURFACE AT WHICH
!     AQUEOUS DIFFUSION THROUGH RESIDUE SURFACE LAYER = GASEOUS
!     DIFFUSION THROUGH ATMOSPHERE BOUNDARY LAYER CALCULATED
!     FROM AQUEOUS DIFFUSIVITY AND BOUNDARY LAYER CONDUCTANCE
!
!     PARR=boundary layer conductance above litter surface from watsub.f
!
    DO idg=idg_beg,idg_NH3
!
!     RGasAtmDisol2LitrM_col: dissolution into litter water, based on two-point gradient flow
!
      trc_gasq(idg)=(PARR_col(NY,NX)*GasSolbility_vr(idg,0,NY,NX)*AtmGasCgperm3(idg,NY,NX) &
        +DFGcc(idg)*trcs_cl_litr(idg))/(DFGcc(idg)+PARR_col(NY,NX))

      RGasAtmDisol2LitrM_col(idg,NY,NX)=(trc_gasq(idg)-trcs_cl_litr(idg))*AMIN1(VLWatMicPM_vr(M,0,NY,NX),DFGcc(idg))

    ENDDO

  ELSE
    RGasAtmDisol2LitrM_col(idg_beg:idg_NH3,NY,NX)=0.0_r8
  ENDIF
  call PrintInfo('end '//subname)
  end subroutine LitterAtmosExchangeM

!------------------------------------------------------------------------------------------
  subroutine SoilAtmosExchangeM(I,J,M,NY,NX,CDOM_MicP_soil,trcs_cl_soil)
  !
  !Description
  !do soil-atmosphere tracer exchange 
  implicit none
  integer,  intent(in) :: I,J  
  integer,  intent(in) :: M,NY,NX
  real(r8), intent(out) :: trcs_cl_soil(ids_beg:ids_end)
  real(r8), intent(out) :: CDOM_MicP_soil(idom_beg:idom_end,1:jcplx)

  character(len=*), parameter :: subname='SoilAtmosExchangeM'
  real(r8) :: DLYR1,TORT1,VLWatMicPOA,VLWatMicPOB,VLWatMicPPA,VLWatMicPPB
  real(r8) :: DiffusivitySolutEff
  real(r8) :: trc_gsolc,trc_gsolc2
  integer  :: K,idg,idom
!
!     SURFACE EXCHANGE OF AQUEOUS CO2, CH4, O2, N2, NH3
!     THROUGH VOLATILIZATION-DISSOLUTION FROM AQUEOUS
!     DIFFUSIVITIES IN SURFACE SOIL LAYER
!
  call PrintInfo('beg '//subname)

  IF(VLWatMicPM_vr(M,NU(NY,NX),NY,NX).GT.ZEROS2(NY,NX))THEN
    VLWatMicPOA = VLWatMicPM_vr(M,NU(NY,NX),NY,NX)*trcs_VLN_vr(ids_NO3,NU(NY,NX),NY,NX)
    VLWatMicPOB = VLWatMicPM_vr(M,NU(NY,NX),NY,NX)*trcs_VLN_vr(ids_NO3B,NU(NY,NX),NY,NX)
    VLWatMicPPA = VLWatMicPM_vr(M,NU(NY,NX),NY,NX)*trcs_VLN_vr(ids_H1PO4,NU(NY,NX),NY,NX)
    VLWatMicPPB = VLWatMicPM_vr(M,NU(NY,NX),NY,NX)*trcs_VLN_vr(ids_H1PO4B,NU(NY,NX),NY,NX)
    DLYR1       = AMAX1(ZERO2,DLYR_3D(3,NU(NY,NX),NY,NX))
    TORT1       = TortMicPM_vr(M,NU(NY,NX),NY,NX)*AREA(3,NU(NY,NX),NY,NX)/(0.5_r8*DLYR1)

    D8910: DO  K=1,jcplx
      do idom=idom_beg,idom_end
        CDOM_MicP_soil(idom,K)=AZMAX1(DOM_MicP2_vr(idom,K,NU(NY,NX),NY,NX)/VLWatMicPM_vr(M,NU(NY,NX),NY,NX))
      enddo
    ENDDO D8910

   !excldue  NH3 and NH3B
    DO idg=idg_beg,idg_NH3-1
      trcs_cl_soil(idg)=trcs_solml2_vr(idg,NU(NY,NX),NY,NX)/VLWatMicPM_vr(M,NU(NY,NX),NY,NX)
    ENDDO

    !deal with inorganic nutrients, including NH3 and NH3B
    IF(VLWatMicPMA_vr(NU(NY,NX),NY,NX).GT.ZEROS2(NY,NX))THEN
      trcs_cl_soil(idg_NH3) = AZMAX1(trcs_solml2_vr(idg_NH3,NU(NY,NX),NY,NX)/VLWatMicPMA_vr(NU(NY,NX),NY,NX))
      trcs_cl_soil(ids_NH4) = AZMAX1(trcs_solml2_vr(ids_NH4,NU(NY,NX),NY,NX)/VLWatMicPMA_vr(NU(NY,NX),NY,NX))
    ELSE
      trcs_cl_soil(idg_NH3) = 0.0_r8
      trcs_cl_soil(ids_NH4) = 0.0_r8
    ENDIF

    IF(VLWatMicPOA.GT.ZEROS2(NY,NX))THEN
      trcs_cl_soil(ids_NO3) = AZMAX1(trcs_solml2_vr(ids_NO3,NU(NY,NX),NY,NX)/VLWatMicPOA)
      trcs_cl_soil(ids_NO2) = AZMAX1(trcs_solml2_vr(ids_NO2,NU(NY,NX),NY,NX)/VLWatMicPOA)
    ELSE
      trcs_cl_soil(ids_NO3) = 0.0_r8
      trcs_cl_soil(ids_NO2) = 0.0_r8
    ENDIF

    IF(VLWatMicPPA.GT.ZEROS2(NY,NX))THEN
      trcs_cl_soil(ids_H1PO4) = AZMAX1(trcs_solml2_vr(ids_H1PO4,NU(NY,NX),NY,NX)/VLWatMicPPA)
      trcs_cl_soil(ids_H2PO4) = AZMAX1(trcs_solml2_vr(ids_H2PO4,NU(NY,NX),NY,NX)/VLWatMicPPA)
    ELSE
      trcs_cl_soil(ids_H1PO4) = 0.0_r8
      trcs_cl_soil(ids_H2PO4) = 0.0_r8
    ENDIF
    IF(VLWatMicPMB_vr(NU(NY,NX),NY,NX).GT.ZEROS2(NY,NX))THEN
      trcs_cl_soil(idg_NH3B) = AZMAX1(trcs_solml2_vr(idg_NH3B,NU(NY,NX),NY,NX)/VLWatMicPMB_vr(NU(NY,NX),NY,NX))
      trcs_cl_soil(ids_NH4B) = AZMAX1(trcs_solml2_vr(ids_NH4B,NU(NY,NX),NY,NX)/VLWatMicPMB_vr(NU(NY,NX),NY,NX))
    ELSE
      trcs_cl_soil(idg_NH3B) = trcs_cl_soil(idg_NH3)
      trcs_cl_soil(ids_NH4B) = trcs_cl_soil(ids_NH4)
    ENDIF

    IF(VLWatMicPOB.GT.ZEROS2(NY,NX))THEN
      trcs_cl_soil(ids_NO3B) = AZMAX1(trcs_solml2_vr(ids_NO3B,NU(NY,NX),NY,NX)/VLWatMicPOB)
      trcs_cl_soil(ids_NO2B)  = AZMAX1(trcs_solml2_vr(ids_NO2B,NU(NY,NX),NY,NX)/VLWatMicPOB)
    ELSE
      trcs_cl_soil(ids_NO3B) = trcs_cl_soil(ids_NO3)
      trcs_cl_soil(ids_NO2B) = trcs_cl_soil(ids_NO2)
    ENDIF

    IF(VLWatMicPPB.GT.ZEROS2(NY,NX))THEN
      trcs_cl_soil(ids_H1PO4B) = AZMAX1(trcs_solml2_vr(ids_H1PO4B,NU(NY,NX),NY,NX)/VLWatMicPPB)
      trcs_cl_soil(ids_H2PO4B) = AZMAX1(trcs_solml2_vr(ids_H2PO4B,NU(NY,NX),NY,NX)/VLWatMicPPB)
    ELSE
      trcs_cl_soil(ids_H1PO4B) = trcs_cl_soil(ids_H1PO4)
      trcs_cl_soil(ids_H2PO4B) = trcs_cl_soil(ids_H2PO4)
    ENDIF
    !
    !     EQUILIBRIUM CONCENTRATIONS AT SOIL SURFACE AT WHICH
    !     AQUEOUS DIFFUSION THROUGH SOIL SURFACE LAYER = GASEOUS
    !     DIFFUSION THROUGH ATMOSPHERE BOUNDARY LAYER CALCULATED
    !     FROM AQUEOUS DIFFUSIVITY AND BOUNDARY LAYER CONDUCTANCE

    !     SURFACE VOLATILIZATION-DISSOLUTION FROM DIFFERENCES
    !     BETWEEN ATMOSPHERIC AND SOIL SURFACE EQUILIBRIUM
    !     CONCENTRATIONS
    ! include NH3B

    DO idg=idg_beg,idg_end
      DiffusivitySolutEff = SoluteDifusivitytscaledM_vr(idg,NU(NY,NX),NY,NX)*TORT1
      trc_gsolc           = (CondGasXSnowM_col(M,NY,NX)*AtmGasCgperm3(idg,NY,NX)*GasSolbility_vr(idg,NU(NY,NX),NY,NX) &
        +DiffusivitySolutEff*trcs_cl_soil(idg))/(DiffusivitySolutEff+CondGasXSnowM_col(M,NY,NX))

      RGasAtmDisol2SoilM_col(idg,NY,NX)=(trc_gsolc-trcs_cl_soil(idg)) &
        *AMIN1(VLWatMicPM_vr(M,NU(NY,NX),NY,NX)*trcs_VLN_vr(idg,NU(NY,NX),NY,NX),DiffusivitySolutEff)

    ENDDO

  ELSE
    RGasAtmDisol2SoilM_col(idg_beg:idg_end,NY,NX)      = 0.0_r8
  ENDIF

  call PrintInfo('end '//subname)
  end subroutine SoilAtmosExchangeM

!------------------------------------------------------------------------------------------

  subroutine Lit2SoilTracerAdvectM(I,J,M,NY,NX,FLWRM1,solute_adv_Lit2Soil_flxM,DOM_Adv_Litr2Soil_flxM)
  !
  ! Description:
  ! Litter and top soil tracer exchange by (upstream) advection
  !
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: M, NY, NX
  real(r8),intent(out) :: solute_adv_Lit2Soil_flxM(ids_beg:ids_end)
  real(r8),intent(out) :: FLWRM1    !water flow from litter to soil at iteration M
  real(r8),intent(out) :: DOM_Adv_Litr2Soil_flxM(idom_beg:idom_end,1:jcplx)

  character(len=*), parameter :: subname='Lit2SoilTracerAdvectM'
  REAL(R8) :: VFLW
  integer :: K,idn,idom,idg,idn1

  call PrintInfo('beg '//subname)

  FLWRM1=WatFLoLitr2SoilM_col(M,NY,NX)
!
!
  !water flow into litter
  IF(FLWRM1.GT.0.0_r8)THEN
    IF(VLWatMicPM_vr(M,0,NY,NX).GT.ZEROS2(NY,NX))THEN
      VFLW=AZMAX1(AMIN1(VFLWX,FLWRM1/VLWatMicPM_vr(M,0,NY,NX)))
    ELSE
      VFLW=VFLWX
    ENDIF

    DO  K=1,jcplx
      DO idom=idom_beg,idom_end
        DOM_Adv_Litr2Soil_flxM(idom,K)=VFLW*AZMAX1(DOM_MicP2_vr(idom,K,0,NY,NX))
      ENDDO
    ENDDO

    !all volatiles
    DO idg=idg_beg,idg_end
      solute_adv_Lit2Soil_flxM(idg)=VFLW*AZMAX1(trcs_solml2_vr(idg,0,NY,NX))
    enddo

    !from NH4 to H2PO4
    DO idn=ids_nut_beg,ids_nuts_end
      solute_adv_Lit2Soil_flxM(idn)=VFLW*AZMAX1(trcs_solml2_vr(idn,0,NY,NX))*trcs_VLN_vr(idn,NU(NY,NX),NY,NX)
    ENDDO

    !from NH4B to ids_H2PO4B
    DO idn=0,ids_nuts     
      idn1=ids_NH4B+idn
      solute_adv_Lit2Soil_flxM(idn1)=VFLW*AZMAX1(trcs_solml2_vr(idn+ids_NH4,0,NY,NX))*trcs_VLN_vr(idn1,NU(NY,NX),NY,NX)
    ENDDO

    solute_adv_Lit2Soil_flxM(idg_NH3B)=VFLW*AZMAX1(trcs_solml2_vr(idg_NH3,0,NY,NX))*trcs_VLN_vr(idg_NH3B,NU(NY,NX),NY,NX)

    ! flow from soil to litter
  ELSE
    IF(VLWatMicPM_vr(M,NU(NY,NX),NY,NX).GT.ZEROS2(NY,NX))THEN
      VFLW=AZMIN1(AMAX1(-VFLWX,FLWRM1/VLWatMicPM_vr(M,NU(NY,NX),NY,NX)))
    ELSE
      VFLW=-VFLWX
    ENDIF
    DO K=1,jcplx
      DO idom=idom_beg,idom_end
        DOM_Adv_Litr2Soil_flxM(idom,K)=VFLW*AZMAX1(DOM_MicP2_vr(idom,K,NU(NY,NX),NY,NX))
      ENDDO
    ENDDO

    DO idn=ids_beg,ids_end
      solute_adv_Lit2Soil_flxM(idn)=VFLW*AZMAX1(trcs_solml2_vr(idn,NU(NY,NX),NY,NX))
    ENDDO
  ENDIF
  call PrintInfo('end '//subname)
  end subroutine Lit2SoilTracerAdvectM
!------------------------------------------------------------------------------------------

  subroutine LitterSoilTracerDiffusionM(I,J,M,NY,NX,CDOM_MicP_litr,CDOM_MicP_soil,FLWRM1,trcs_cl_litr,trcs_cl_soil,&
    trcs_Dif_Litr2Soil_flxM,DOM_Difus_Litr2Soil_flxM)
  !
  !Description
  ! Diffusion flux between litter and soil  
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: M, NY, NX
  real(r8), intent(in) :: FLWRM1                         !water flow from litter to topsoil
  real(r8), intent(in) :: trcs_cl_litr(ids_beg:ids_end)  !tracer concentation in litter
  real(r8), intent(in) :: trcs_cl_soil(ids_beg:ids_end)  !tracer concentation in topsoil
  real(r8), intent(in) :: CDOM_MicP_litr(idom_beg:idom_end,1:jcplx)  !DOM concentation in litter
  real(r8), intent(in) :: CDOM_MicP_soil(idom_beg:idom_end,1:jcplx)  !DOM concentation in topsoil
  real(r8), intent(out):: trcs_Dif_Litr2Soil_flxM(ids_beg:ids_end)   !tracer diffusion from litter to soil
  real(r8), intent(out):: DOM_Difus_Litr2Soil_flxM(idom_beg:idom_end,1:jcplx)  !tracer diffusion from litter to soil

  character(len=*), parameter :: subname='LitterSoilTracerDiffusionM'
  real(r8) :: TORT0,TORT1
  real(r8) :: DLYR0,DLYR1
  real(r8) :: DIFDOM(idom_beg:idom_end)

  real(r8) :: DIFDOM0
  real(r8) :: DIFDOM1
  real(r8) :: DISPN,DIF0,DIF1
  real(r8) :: SDifc(ids_beg:ids_end)
  integer  :: K,idn,idg,idom
!
!     VOLT,DLYR,AREA=soil surface volume, thickness, area
!     VLWatMicPM=micropore water-filled porosity from watsub.f
!
  call PrintInfo('beg '//subname)

  IF((VGeomLayer_vr(0,NY,NX).GT.ZEROS2(NY,NX) .AND. VLWatMicPM_vr(M,0,NY,NX).GT.ZEROS2(NY,NX)) &
    .AND. (VLWatMicPM_vr(M,NU(NY,NX),NY,NX).GT.ZEROS2(NY,NX)))THEN
!
!     DIFFUSIVITIES IN RESIDUE AND SOIL SURFACE

    DLYR0 = AMAX1(ZERO2,DLYR_3D(3,0,NY,NX))
    TORT0 = TortMicPM_vr(M,0,NY,NX)/DLYR0*FracSurfByLitR_col(NY,NX)
    DLYR1 = AMAX1(ZERO2,DLYR_3D(3,NU(NY,NX),NY,NX))
    TORT1 = TortMicPM_vr(M,NU(NY,NX),NY,NX)/DLYR1
    DISPN = DISP_3D(3,NU(NY,NX),NY,NX)*AMIN1(VFLWX,ABS(FLWRM1/AREA(3,NU(NY,NX),NY,NX)))

    DO idom=idom_beg,idom_end
      DIFDOM0      = (DOM_diffusivitytscaledM_vr(idom,0,NY,NX)*TORT0+DISPN)
      DIFDOM1      = (DOM_diffusivitytscaledM_vr(idom,NU(NY,NX),NY,NX)*TORT1+DISPN)
      DIFDOM(idom) = DIFDOM0*DIFDOM1/(DIFDOM0+DIFDOM1)*AREA(3,NU(NY,NX),NY,NX)
    ENDDO

    DO idn=ids_beg,ids_end
      DIF0        = (SoluteDifusivitytscaledM_vr(idn,0,NY,NX)*TORT0+DISPN)
      DIF1        = (SoluteDifusivitytscaledM_vr(idn,NU(NY,NX),NY,NX)*TORT1+DISPN)
      SDifc(idn) = DIF0*DIF1/(DIF0+DIF1)*AREA(3,NU(NY,NX),NY,NX)
    ENDDO
!
    DO  K=1,jcplx
      DO idom=idom_beg,idom_end
        DOM_Difus_Litr2Soil_flxM(idom,K)=DIFDOM(idom)*(CDOM_MicP_litr(idom,K)-CDOM_MicP_soil(idom,K))
      ENDDO
    ENDDO

    !exclude NH3B and NH3
    DO idg=idg_beg,idg_NH3-1
        trcs_Dif_Litr2Soil_flxM(idg)=SDifc(idg)*(trcs_cl_litr(idg)-trcs_cl_soil(idg))
    ENDDO

    !include NH3B and NH3
    DO idn=ids_nuts_beg,ids_nuts_end
      trcs_Dif_Litr2Soil_flxM(idn)=SDifc(idn)*(trcs_cl_litr(idn)-trcs_cl_soil(idn))*trcs_VLN_vr(idn,NU(NY,NX),NY,NX)
    ENDDO
  ELSE
    DO  K=1,jcplx
      DOM_Difus_Litr2Soil_flxM(idom_beg:idom_end,K)=0.0_r8
    ENDDO
    trcs_Dif_Litr2Soil_flxM(ids_beg:ids_end)=0._r8
  ENDIF
  call PrintInfo('end '//subname)
  end subroutine LitterSoilTracerDiffusionM
!------------------------------------------------------------------------------------------

  subroutine SurfLayerNetTracerFluxM(I,J,M,NY,NX,trcs_Dif_Litr2Soil_flxM,solute_adv_Lit2Soil_flxM,&
    DOM_Adv_Litr2Soil_flxM,DOM_Difus_Litr2Soil_flxM)
  use EcoSiMParDataMod, only : micpar
  implicit none
  integer, intent(in) :: I,J,M
  integer, intent(in) :: NY, NX
  real(r8), intent(in) :: trcs_Dif_Litr2Soil_flxM(ids_beg:ids_end)             !diffusion flux
  real(r8), intent(in) :: solute_adv_Lit2Soil_flxM(ids_beg:ids_end)      !advection flux
  real(r8), intent(in) :: DOM_Adv_Litr2Soil_flxM(idom_beg:idom_end,1:jcplx)     !diffusion advection flux
  real(r8), intent(in) :: DOM_Difus_Litr2Soil_flxM(idom_beg:idom_end,1:jcplx)   !DOM diffusion flux
  character(len=*), parameter :: subname ='SurfLayerNetTracerFluxM'   
  integer :: K,idg,idom,ids
  real(r8):: dom_flx

  call PrintInfo('beg '//subname)

  !litter layer
  DO K=1,micpar%NumOfLitrCmplxs
    DO idom=idom_beg,idom_end
      dom_flx=DOM_Adv_Litr2Soil_flxM(idom,K)+DOM_Difus_Litr2Soil_flxM(idom,K)
      DOM_MicpTranspFlxM_3D(idom,K,3,0,NY,NX)         = DOM_Flo2LitrM(idom,K,NY,NX)-dom_flx
      DOM_MicpTranspFlxM_3D(idom,K,3,NU(NY,NX),NY,NX) = DOM_Flo2TopSoilM(idom,K,NY,NX)+dom_flx
    ENDDO
  ENDDO

  DO idg=idg_beg,idg_NH3
    trcs_MicpTranspFlxM_3D(idg,3,0,NY,NX)=trcg_Precip2LitrM(idg,NY,NX)+trcg_AquaADV_Snow2Litr_flxM(idg,NY,NX) &
      -solute_adv_Lit2Soil_flxM(idg)-trcs_Dif_Litr2Soil_flxM(idg)
  ENDDO

  do ids=ids_nut_beg,ids_nuts_end
    trcs_MicpTranspFlxM_3D(ids,3,0,NY,NX)=trcn_Precip2LitrM(ids,NY,NX)+trcn_AquaADV_Snow2Litr_flxM(ids,NY,NX) &
      -solute_adv_Lit2Soil_flxM(ids)-trcs_Dif_Litr2Soil_flxM(ids)
  enddo

  trcs_MicpTranspFlxM_3D(idg_NH3,3,0,NY,NX)   = trcs_MicpTranspFlxM_3D(idg_NH3,3,0,NY,NX)-solute_adv_Lit2Soil_flxM(idg_NH3B)-trcs_Dif_Litr2Soil_flxM(idg_NH3B)
  trcs_MicpTranspFlxM_3D(ids_NH4,3,0,NY,NX)   = trcs_MicpTranspFlxM_3D(ids_NH4,3,0,NY,NX)-solute_adv_Lit2Soil_flxM(ids_NH4B)-trcs_Dif_Litr2Soil_flxM(ids_NH4B)
  trcs_MicpTranspFlxM_3D(ids_NO3,3,0,NY,NX)   = trcs_MicpTranspFlxM_3D(ids_NO3,3,0,NY,NX)-solute_adv_Lit2Soil_flxM(ids_NO3B)-trcs_Dif_Litr2Soil_flxM(ids_NO3B)
  trcs_MicpTranspFlxM_3D(ids_NO2,3,0,NY,NX)   = trcs_MicpTranspFlxM_3D(ids_NO2,3,0,NY,NX)-solute_adv_Lit2Soil_flxM(ids_NO2B)-trcs_Dif_Litr2Soil_flxM(ids_NO2B)
  trcs_MicpTranspFlxM_3D(ids_H1PO4,3,0,NY,NX) = trcs_MicpTranspFlxM_3D(ids_H1PO4,3,0,NY,NX)-solute_adv_Lit2Soil_flxM(ids_H1PO4B)-trcs_Dif_Litr2Soil_flxM(ids_H1PO4B)
  trcs_MicpTranspFlxM_3D(ids_H2PO4,3,0,NY,NX) = trcs_MicpTranspFlxM_3D(ids_H2PO4,3,0,NY,NX)-solute_adv_Lit2Soil_flxM(ids_H2PO4B)-trcs_Dif_Litr2Soil_flxM(ids_H2PO4B)

  !top soil
  DO idg=idg_beg,idg_NH3
    trcs_MicpTranspFlxM_3D(idg,3,NU(NY,NX),NY,NX)=trcs_Precip2MicpM(idg,NY,NX)+trcg_AquaADV_Snow2Soil_flxM(idg,NY,NX) &
      +solute_adv_Lit2Soil_flxM(idg)+trcs_Dif_Litr2Soil_flxM(idg)
  ENDDO

  do ids=ids_nut_beg,ids_nuts_end
    trcs_MicpTranspFlxM_3D(ids,3,NU(NY,NX),NY,NX)=trcs_Precip2MicpM(ids,NY,NX)+trcn_AquaADV_Snow2Soil_flxM(ids,NY,NX) &
      +solute_adv_Lit2Soil_flxM(ids)+trcs_Dif_Litr2Soil_flxM(ids)
  enddo

  do ids=ids_nutb_beg,ids_nutb_end
    trcs_MicpTranspFlxM_3D(ids,3,NU(NY,NX),NY,NX)=trcs_Precip2MicpM(ids,NY,NX)+trcn_AquaADV_Snow2Band_flxM(ids,NY,NX) &
      +solute_adv_Lit2Soil_flxM(ids)+trcs_Dif_Litr2Soil_flxM(ids)
  enddo

  call PrintInfo('end '//subname)
  end subroutine SurfLayerNetTracerFluxM

!------------------------------------------------------------------------------------------

  subroutine SurfXGridTracerTransptM(M,NY,NX,NHE,NHW,NVS,NVN)
  !
  !Description:
  !Do inter-grid lateral transport at the surface
  !
  implicit none
  integer, intent(in) :: M, NY, NX, NHE, NHW, NVS, NVN
  character(len=*), parameter :: subname='SurfXGridTracerTransptM'

  real(r8) :: FQRM,VFLW
  integer :: K,N,NN,N1,N2,N4,N5,N4B,N5B,idg,idn
  integer :: idom
!
  call PrintInfo('beg '//subname)

  N1=NX;N2=NY

  !DO SURFACE RUNOFF  
  IF(SurfRunoffPotentM_col(M,N2,N1).GT.ZEROS(N2,N1))THEN
    !Obtain potential tracer flux through surface runoff, include both inner grids and boundaries  
    IF(VLWatMicPM_vr(M,0,N2,N1).GT.ZEROS2(N2,N1))THEN
      VFLW=AMIN1(VFLWX,SurfRunoffPotentM_col(M,N2,N1)/VLWatMicPM_vr(M,0,N2,N1))
    ELSE
      VFLW=VFLWX
    ENDIF

    DO  K=1,jcplx
      DO idom=idom_beg,idom_end
        DOM_FloXSurRunoff_PotFlxM(idom,K,N2,N1)=VFLW*AZMAX1(DOM_MicP2_vr(idom,K,0,N2,N1))
      enddo
    ENDDO

    DO idg=idg_beg,idg_NH3
      trcg_FloXSurRunoff_PotFlxM(idg,N2,N1)=VFLW*AZMAX1(trcs_solml2_vr(idg,0,N2,N1))
    ENDDO

    DO idn=ids_nut_beg,ids_nuts_end
      trcn_FloXSurRunoff_PotFlxM(idn,N2,N1)=VFLW*AZMAX1(trcs_solml2_vr(idn,0,N2,N1))
    ENDDO
    !-----------------------------------------------------------------
    !     LOCATE INTERNAL BOUNDARIES BETWEEN ADJACENT GRID CELLS
    !
    D4311: DO N=1,2
      D4306: DO NN=1,2

        IF(N.EQ.iWestEastDirection)THEN  !skip boundaries          
          IF((NX.EQ.NHE .AND. NN.EQ.iFront) .OR. (NX.EQ.NHW .AND. NN.EQ.iBehind))THEN
            cycle
          ELSE
            N4  = NX+1;N5  = NY
            N4B = NX-1;N5B = NY
          ENDIF
        ELSEIF(N.EQ.iNorthSouthDirection)THEN  !skip boundaries          
          IF((NY.EQ.NVS .AND. NN.EQ.iFront) .OR. (NY.EQ.NVN .AND. NN.EQ.iBehind))THEN
            cycle
          ELSE
            N4  = NX;N5  = NY+1
            N4B = NX;N5B = NY-1
          ENDIF
        ENDIF      
        !
        ! IF OVERLAND FLOW IS FROM CURRENT TO ADJACENT GRID CELL
        !
        IF(NN.EQ.iFront)THEN
          FQRM=QflxSurfRunoffM_2DH(M,N,iBehind,N5,N4)/SurfRunoffPotentM_col(M,N2,N1)

          IF(FQRM.GT.ZEROS(N2,N1))then
            DO  K=1,jcplx
              do idom=idom_beg,idom_end
                DOM_FloXSurRof_flxM_2DH(idom,K,N,iBehind,N5,N4)=DOM_FloXSurRunoff_PotFlxM(idom,K,N2,N1)*FQRM
              enddo
            ENDDO

            DO idg=idg_beg,idg_NH3
              trcg_SurRof_flxM_2DH(idg,N,iBehind,N5,N4)=trcg_FloXSurRunoff_PotFlxM(idg,N2,N1)*FQRM
            ENDDO

            DO idn=ids_nut_beg,ids_nuts_end
              trcn_SurRof_flxM_2DH(idn,N,iBehind,N5,N4)=trcn_FloXSurRunoff_PotFlxM(idn,N2,N1)*FQRM
            ENDDO
            !------------------------------------------------
            ! ACCUMULATE HOURLY FLUXES FOR USE IN REDIST.F
            DO  K=1,jcplx
              do idom=idom_beg,idom_end
                DOM_FloXSurRunoff_2DH(idom,K,N,iBehind,N5,N4)=DOM_FloXSurRunoff_2DH(idom,K,N,iBehind,N5,N4) &
                  +DOM_FloXSurRof_flxM_2DH(idom,K,N,iBehind,N5,N4)
              enddo
            ENDDO

            DO idg=idg_beg,idg_NH3
              trcg_FloXSurRunoff_2D(idg,N,iBehind,N5,N4)=trcg_FloXSurRunoff_2D(idg,N,iBehind,N5,N4) &
                +trcg_SurRof_flxM_2DH(idg,N,iBehind,N5,N4)
            ENDDO

            DO idn=ids_nut_beg,ids_nuts_end
              trcn_FloXSurRunoff_2D(idn,N,iBehind,N5,N4)=trcn_FloXSurRunoff_2D(idn,N,iBehind,N5,N4) &
                +trcn_SurRof_flxM_2DH(idn,N,iBehind,N5,N4)
            ENDDO
          ELSE
            DO K=1,jcplx
              DOM_FloXSurRof_flxM_2DH(idom_beg:idom_end,K,N,iBehind,N5,N4)=0.0_r8
            ENDDO

            trcg_SurRof_flxM_2DH(idg_beg:idg_NH3,N,iBehind,N5,N4)          = 0.0_r8
            trcn_SurRof_flxM_2DH(ids_nut_beg:ids_nuts_end,N,iBehind,N5,N4) = 0.0_r8
          ENDIF

          !  IF OVERLAND FLOW IS FROM CURRENT TO ADJACENT GRID CELL
        ELSEIF(NN.EQ.iBehind)THEN        
          !legitimate grid
          IF(N4B.GT.0 .AND. N5B.GT.0)THEN

            FQRM=QflxSurfRunoffM_2DH(M,N,iFront,N5B,N4B)/SurfRunoffPotentM_col(M,N2,N1)
            IF(FQRM.GT.ZEROS(N2,N1))then
              DO  K=1,jcplx
                do idom=idom_beg,idom_end
                  DOM_FloXSurRof_flxM_2DH(idom,K,N,iFront,N5B,N4B)=DOM_FloXSurRunoff_PotFlxM(idom,K,N2,N1)*FQRM
                enddo
              ENDDO

              DO idg=idg_beg,idg_NH3
                trcg_SurRof_flxM_2DH(idg,N,iFront,N5B,N4B)=trcg_FloXSurRunoff_PotFlxM(idg,N2,N1)*FQRM
              ENDDO

              DO idn=ids_nut_beg,ids_nuts_end
                trcn_SurRof_flxM_2DH(idn,N,iFront,N5B,N4B)=trcn_FloXSurRunoff_PotFlxM(idn,N2,N1)*FQRM
              ENDDO
              !-----------------
              !accumulate fluxes
              DO K=1,jcplx
                do idom=idom_beg,idom_end
                  DOM_FloXSurRunoff_2DH(idom,K,N,iFront,N5B,N4B)=DOM_FloXSurRunoff_2DH(idom,K,N,iFront,N5B,N4B) &
                    +DOM_FloXSurRof_flxM_2DH(idom,K,N,iFront,N5B,N4B)
                enddo
              ENDDO

              DO idg=idg_beg,idg_NH3
                trcg_FloXSurRunoff_2D(idg,N,iFront,N5B,N4B)=trcg_FloXSurRunoff_2D(idg,N,iFront,N5B,N4B) &
                  +trcg_SurRof_flxM_2DH(idg,N,iFront,N5B,N4B)
              ENDDO

              DO idn=ids_nut_beg,ids_nuts_end
                trcn_FloXSurRunoff_2D(idn,N,iFront,N5B,N4B)=trcn_FloXSurRunoff_2D(idn,N,iFront,N5B,N4B) &
                  +trcn_SurRof_flxM_2DH(idn,N,iFront,N5B,N4B)
              ENDDO
            ELSE
              DO  K=1,jcplx
                DOM_FloXSurRof_flxM_2DH(idom_beg:idom_end,K,N,iFront,N5B,N4B)=0.0_r8
              ENDDO
              trcg_SurRof_flxM_2DH(idg_beg:idg_end,N,iFront,N5B,N4B)          = 0.0_r8
              trcn_SurRof_flxM_2DH(ids_nut_beg:ids_nuts_end,N,iFront,N5B,N4B) = 0.0_r8
            ENDIF
          ENDIF  
        ENDIF        
      ENDDO  D4306
    ENDDO D4311
  ELSE
    DO K=1,jcplx
      DOM_FloXSurRunoff_PotFlxM(idom_beg:idom_end,K,N2,N1)=0.0_r8
    ENDDO
    trcg_FloXSurRunoff_PotFlxM(idg_beg:idg_NH3,N2,N1)          = 0.0_r8
    trcn_FloXSurRunoff_PotFlxM(ids_nut_beg:ids_nuts_end,N2,N1) = 0.0_r8
  ENDIF
  !----------------------------------
  !DO SNOW REDISTRIBUTION
  !
  D4310: DO N=1,2
    D4305: DO NN=1,2
      !
      IF(N.EQ.iWestEastDirection)THEN
        !skip boundaries
        IF((NX.EQ.NHE .AND. NN.EQ.iFront) .OR. (NX.EQ.NHW .AND. NN.EQ.iBehind))THEN
          cycle
        ELSE
          N4  = NX+1;N5  = NY
          N4B = NX-1;N5B = NY
        ENDIF
      ELSEIF(N.EQ.iNorthSouthDirection)THEN
        !skip boundaries
        IF((NY.EQ.NVS .AND. NN.EQ.iFront) .OR. (NY.EQ.NVN .AND. NN.EQ.iBehind))THEN
          cycle
        ELSE
          N4  = NX;N5  = NY+1
          N4B = NX;N5B = NY-1
        ENDIF
      ENDIF

      IF(NN.EQ.iFront)THEN
        !
        !   IF NO SNOW DRIFT THEN NO TRANSPORT
        !
        IF(ABS(DrySnoFlxByRedistM_2DH(M,N,iBehind,N5,N4)).LE.ZEROS2(N2,N1))THEN
          trcg_SnowDrift_flxM_2D(idg_beg:idg_NH3,N,iBehind,N5,N4)          = 0.0_r8
          trcn_SnowDrift_flxM_2D(ids_nut_beg:ids_nuts_end,N,iBehind,N5,N4) = 0.0_r8
          VFLW                                                             = 0._r8
          !
          !  IF DRIFT IS FROM CURRENT TO ADJACENT GRID CELL
          !
        ELSEIF(DrySnoFlxByRedistM_2DH(M,N,iBehind,N5,N4).GT.ZEROS2(N2,N1))THEN
          IF(VcumSnoDWI_col(N2,N1).GT.ZEROS2(N2,N1))THEN
            VFLW=AZMAX1(AMIN1(VFLWX,DrySnoFlxByRedistM_2DH(M,N,iBehind,N5,N4)/VcumSnoDWI_col(N2,N1)))
          ELSE
            VFLW=VFLWX
          ENDIF

          DO idg=idg_beg,idg_NH3
            trcg_SnowDrift_flxM_2D(idg,N,iBehind,N5,N4) = VFLW*AZMAX1(trcg_solsml2_snvr(idg,1,N2,N1))
          ENDDO

          DO idn=ids_nut_beg,ids_nuts_end
            trcn_SnowDrift_flxM_2D(idn,N,iBehind,N5,N4) = VFLW*AZMAX1(trcn_solsml2_snvr(idn,1,N5,N4))
          ENDDO
        ENDIF  
        !check legitimate behind grid  
      ELSEIF(NN.EQ.iBehind .AND. N5B.GT.0 .AND. N4B.GT.0)THEN   

        IF(ABS(DrySnoFlxByRedistM_2DH(M,N,iFront,N5B,N4B)).LE.ZEROS2(N2,N1))THEN
          trcg_SnowDrift_flxM_2D(idg_beg:idg_NH3,N,iFront,N5B,N4B)          = 0.0_r8
          trcn_SnowDrift_flxM_2D(ids_nut_beg:ids_nuts_end,N,iFront,N5B,N4B) = 0.0_r8
          VFLW                                                              = 0._r8
          !
          !     IF DRIFT IS FROM CURRENT TO ADJACENT GRID CELL
          !
        ELSEIF(DrySnoFlxByRedistM_2DH(M,N,iFront,N5B,N4B).GT.ZEROS2(N2,N1))THEN
          IF(VcumSnoDWI_col(N2,N1).GT.ZEROS2(N2,N1))THEN
            VFLW=AZMAX1(AMIN1(VFLWX,DrySnoFlxByRedistM_2DH(M,N,iFront,N5B,N4B)/VcumSnoDWI_col(N2,N1)))
          ELSE
            VFLW=VFLWX
          ENDIF
          DO idg=idg_beg,idg_NH3
            trcg_SnowDrift_flxM_2D(idg,N,iFront,N5B,N4B) = VFLW*AZMAX1(trcg_solsml2_snvr(idg,1,N2,N1))
          ENDDO

          DO idn=ids_nut_beg,ids_nuts_end
            trcn_SnowDrift_flxM_2D(idn,N,iFront,N5B,N4B) = VFLW*AZMAX1(trcn_solsml2_snvr(idn,1,N2,N1))
          ENDDO          
        ENDIF
      ENDIF  
      
    ENDDO D4305
  ENDDO D4310
  call PrintInfo('end '//subname)
  end subroutine SurfXGridTracerTransptM
!------------------------------------------------------------------------------------------

  subroutine SnowSoluteVertDrainM(I,J,M,NHE,NHW,NVS,NVN)
  !
  !Description
  !Vertical drainage of solute in and out of snow.
  !out of snow fluxes are added to litter, soil and nutrient band, accordingly. 
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: M,NHE,NHW,NVS,NVN
  integer :: NY, NX

  character(len=*), parameter :: subname='SnowSoluteVertDrainM'
  logical :: LCHKBOTML  !flag for dealing with the bottom snow layer
  integer :: L,L2
  real(r8) :: VFLWS,VFLWW,VFLWR
  real(r8) :: VFLWNOB,VFLWNO3,VFLWNHB
  real(r8) :: VFLWNH4,VFLWPO4,VFLWPOB
  integer  :: idg,idn
!
  call PrintInfo('beg '//subname)

  DO NX=NHW,NHE
    DO  NY=NVN,NVS

      !initialize with zeros
      trcg_AquaADV_Snow2Litr_flxM(idg_beg:idg_NH3,NY,NX)           = 0.0_r8
      trcn_AquaADV_Snow2Litr_flxM(ids_nut_beg:ids_nuts_end,NY,NX)  = 0.0_r8
      trcg_AquaADV_Snow2Soil_flxM(idg_beg:idg_end,NY,NX)           = 0.0_r8
      trcn_AquaADV_Snow2Soil_flxM(ids_nut_beg:ids_nuts_end,NY,NX)  = 0.0_r8
      trcn_AquaADV_Snow2Band_flxM(ids_nutb_beg:ids_nutb_end,NY,NX) = 0.0_r8

      LCHKBOTML=.false.
      DO  L=1,JS

        !check having meaningful snow mass
        IF(VLSnowHeatCapM_snvr(M,L,NY,NX).GT.VLHeatCapSnowMin_col(NY,NX)*1.e-4_r8)THEN
          !next layer index
          L2=MIN(JS,L+1)
          !L is inner layers
          IF(L.LT.JS .AND. VLSnowHeatCapM_snvr(M,L2,NY,NX).GT.VLHeatCapSnowMin_col(NY,NX))THEN
            !some flowout
            IF(VLWatSnow_snvr(L,NY,NX).GT.ZEROS2(NY,NX))THEN
              VFLWW=AZMAX1(AMIN1(1.0_r8,WatFlowInSnowM_snvr(M,L2,NY,NX)/VLWatSnow_snvr(L,NY,NX)))
              !no water is layer L, all water is gone to L2
            ELSE          
              VFLWW=1.0_r8
            ENDIF

            DO idg=idg_beg,idg_NH3
              trcg_AquaAdv_flxM_snvr(idg,L2,NY,NX) = trcg_solsml2_snvr(idg,L,NY,NX)*VFLWW
              trcg_AquaAdv_flx_snvr(idg,L2,NY,NX)  = trcg_AquaAdv_flx_snvr(idg,L2,NY,NX)+trcg_AquaAdv_flxM_snvr(idg,L2,NY,NX)
            ENDDO

            DO idn=ids_nut_beg,ids_nuts_end
              trcn_AquaAdv_flxM_snvr(idn,L2,NY,NX) = trcn_solsml2_snvr(idn,L,NY,NX)*VFLWW
              !accumulate flux
              trcn_AquaAdv_flx_snvr(idn,L2,NY,NX)  = trcn_AquaAdv_flx_snvr(idn,L2,NY,NX)+trcn_AquaAdv_flxM_snvr(idn,L2,NY,NX)
            ENDDO
            !L == JS: bottom layer or L<LS and next layer has no significant snow
          ELSE
            !nothing to add to L2
            IF(L.LT.JS)THEN
              trcg_AquaAdv_flxM_snvr(idg_beg:idg_NH3,L2,NY,NX)          = 0.0_r8
              trcn_AquaAdv_flxM_snvr(ids_nut_beg:ids_nuts_end,L2,NY,NX) = 0.0_r8
            ENDIF
    !
    !     SNOWPACK SOLUTE DISCHARGE TO SURFACE LITTER, SOIL SURFACE
    !
            IF(.not.LCHKBOTML)THEN
              IF(VLWatSnow_snvr(L,NY,NX).GT.ZEROS2(NY,NX))THEN
                VFLWR=AZMAX1(AMIN1(1.0_r8,WatFlowSno2LitRM_col(M,NY,NX)/VLWatSnow_snvr(L,NY,NX)))
                VFLWS=AZMAX1(AMIN1(1.0_r8,(WatFlowSno2MicPM_col(M,NY,NX)+WatFlowSno2MacPM_col(M,NY,NX))/VLWatSnow_snvr(L,NY,NX)))
              ELSE
                VFLWR=FracSurfByLitR_col(NY,NX)
                VFLWS=FracSurfBareSoil_col(NY,NX)
              ENDIF
              VFLWNH4 = VFLWS*trcs_VLN_vr(ids_NH4,NU(NY,NX),NY,NX)
              VFLWNHB = VFLWS*trcs_VLN_vr(ids_NH4B,NU(NY,NX),NY,NX)
              VFLWNO3 = VFLWS*trcs_VLN_vr(ids_NO3,NU(NY,NX),NY,NX)
              VFLWNOB = VFLWS*trcs_VLN_vr(ids_NO3B,NU(NY,NX),NY,NX)
              VFLWPO4 = VFLWS*trcs_VLN_vr(ids_H1PO4,NU(NY,NX),NY,NX)
              VFLWPOB = VFLWS*trcs_VLN_vr(ids_H1PO4B,NU(NY,NX),NY,NX)
              !snow flow to litter layer
              DO idg=idg_beg,idg_NH3
                trcg_AquaADV_Snow2Litr_flxM(idg,NY,NX)=trcg_solsml2_snvr(idg,L,NY,NX)*VFLWR
                trcg_AquaADV_Snow2Litr_flx(idg,NY,NX)=trcg_AquaADV_Snow2Litr_flx(idg,NY,NX)+trcg_AquaADV_Snow2Litr_flxM(idg,NY,NX)
              ENDDO

              DO idn=ids_nut_beg,ids_nuts_end
                trcn_AquaADV_Snow2Litr_flxM(idn,NY,NX)=trcn_solsml2_snvr(idn,L,NY,NX)*VFLWR
                !accumulate flux
                trcn_AquaADV_Snow2Litr_flx(idn,NY,NX)=trcn_AquaADV_Snow2Litr_flx(idn,NY,NX)+trcn_AquaADV_Snow2Litr_flxM(idn,NY,NX)
              ENDDO
              !-----------------------
              !snow flow to topsoil 
              DO idg=idg_beg,idg_NH3-1
                trcg_AquaADV_Snow2Soil_flxM(idg,NY,NX)  = trcg_solsml2_snvr(idg,L,NY,NX)*VFLWS
              ENDDO

              trcg_AquaADV_Snow2Soil_flxM(idg_NH3,NY,NX)  = trcg_solsml2_snvr(idg_NH3,L,NY,NX)*VFLWNH4
              trcg_AquaADV_Snow2Soil_flxM(idg_NH3B,NY,NX) = trcg_solsml2_snvr(idg_NH3,L,NY,NX)*VFLWNHB

              trcn_AquaADV_Snow2Soil_flxM(ids_NO2,NY,NX)    = 0._r8
              trcn_AquaADV_Snow2Soil_flxM(ids_NH4,NY,NX)    = trcn_solsml2_snvr(ids_NH4,L,NY,NX)*VFLWNH4
              trcn_AquaADV_Snow2Soil_flxM(ids_NO3,NY,NX)    = trcn_solsml2_snvr(ids_NO3,L,NY,NX)*VFLWNO3
              trcn_AquaADV_Snow2Soil_flxM(ids_H1PO4,NY,NX)  = trcn_solsml2_snvr(ids_H1PO4,L,NY,NX)*VFLWPO4
              trcn_AquaADV_Snow2Soil_flxM(ids_H2PO4,NY,NX)  = trcn_solsml2_snvr(ids_H2PO4,L,NY,NX)*VFLWPO4

              trcn_AquaADV_Snow2Band_flxM(ids_NO2B,NY,NX)   = 0._r8
              trcn_AquaADV_Snow2Band_flxM(ids_NH4B,NY,NX)   = trcn_solsml2_snvr(ids_NH4,L,NY,NX)*VFLWNHB
              trcn_AquaADV_Snow2Band_flxM(ids_NO3B,NY,NX)   = trcn_solsml2_snvr(ids_NO3,L,NY,NX)*VFLWNOB
              trcn_AquaADV_Snow2Band_flxM(ids_H1PO4B,NY,NX) = trcn_solsml2_snvr(ids_H1PO4,L,NY,NX)*VFLWPOB
              trcn_AquaADV_Snow2Band_flxM(ids_H2PO4B,NY,NX) = trcn_solsml2_snvr(ids_H2PO4,L,NY,NX)*VFLWPOB

              !accumulate flux
              DO idg=idg_beg,idg_end
                trcg_AquaADV_Snow2Soil_flx(idg,NY,NX) = trcg_AquaADV_Snow2Soil_flx(idg,NY,NX)+trcg_AquaADV_Snow2Soil_flxM(idg,NY,NX)
              enddo
              
              DO idn=ids_nut_beg,ids_nuts_end
                trcn_AquaADV_Snow2Soil_flx(idn,NY,NX) = trcn_AquaADV_Snow2Soil_flx(idn,NY,NX)+trcn_AquaADV_Snow2Soil_flxM(idn,NY,NX)
              ENDDO

              DO idn=ids_nutb_beg,ids_nutb_end
                trcn_AquaADV_Snow2Band_flx(idn,NY,NX) = trcn_AquaADV_Snow2Band_flx(idn,NY,NX)+trcn_AquaADV_Snow2Band_flxM(idn,NY,NX)
              ENDDO
              LCHKBOTML = .true.
            ENDIF
          ENDIF
        ENDIF
      ENDDO
    ENDDO
  ENDDO
  call PrintInfo('end '//subname)
  end subroutine SnowSoluteVertDrainM

! ----------------------------------------------------------------------
  subroutine SoluteAdvDifusTranspM(I,J,M,N,N1,N2,N3,N4,N5,N6)
  implicit none
  integer, intent(in) :: I,J,M,N
  integer, intent(in) :: N1,N2,N3  !source grid
  integer, intent(in) :: N4,N5,N6  !dest grid
  character(len=*), parameter :: subname='SoluteAdvDifusTranspM'
  real(r8)  :: DOM_Difus_Mac2Micp_flxM(idom_beg:idom_end,1:jcplx)  

  integer  :: K,ids,idom

  call PrintInfo('beg '//subname)

  !------------------------------------------------
  !     SOLUTE TRANSPORT IN MICROPORES
  !
  call MicroporeSoluteAdvectionM(M,N,N1,N2,N3,N4,N5,N6)
  !
  call MicroporeSoluteDiffusionM(M,N,N1,N2,N3,N4,N5,N6)
  !------------------------------------------------
  !     SOLUTE TRANSPORT IN MACROPORES
  !
  call MacroporeSoluteAdvectionM(M,N,N1,N2,N3,N4,N5,N6)
  !
  call MacroporeSoluteDispersionM(M,N,N1,N2,N3,N4,N5,N6)

  call PrintInfo('end '//subname)
  end subroutine SoluteAdvDifusTranspM

! ----------------------------------------------------------------------
  subroutine MicroporeSoluteDiffusionM(M,N,N1,N2,N3,N4,N5,N6)
  implicit none
  integer , intent(in) :: M,N,N1,N2,N3,N4,N5,N6
  real(r8) :: DOM_Difus_Mac2Micp_flxM(idom_beg:idom_end,1:jcplx)

  character(len=*), parameter :: subname='MicroporeSoluteDiffusionM'
  real(r8) :: THETW1_vr(JZ,JY,JX)  !soil saturation   
  real(r8) :: VLWatMicPOA,VLWatMicPOB,VLWatMicPPA,VLWatMicPPB
  real(r8) :: VLWatMicP2A,VLWatMicP2B,VLWatMicP3A,VLWatMicP3B,VLWatMicP4A,VLWatMicP4B
  real(r8) :: trcsolc1(ids_beg:ids_end)
  real(r8) :: trcsolc2(ids_beg:ids_end)
  real(r8) :: CDOM_MicP1(idom_beg:idom_end,1:jcplx)
  real(r8) :: CDOM_MicP2_vr(idom_beg:idom_end,1:jcplx)
  real(r8) :: SDifc(ids_beg:ids_end),SDifFlx(ids_beg:ids_end)
  real(r8) :: DISPN,DOMDifc(idom_beg:idom_end)
  real(r8) :: DLYR1,DLYR2,TORTL
  integer  :: K,ids,idg,idom

  call PrintInfo('beg '//subname)
  THETW1_vr(N3,N2,N1)=AZMAX1(safe_adb(VLWatMicPM_vr(M,N3,N2,N1),VLSoilMicP_vr(N3,N2,N1)))
  THETW1_vr(N6,N5,N4)=AZMAX1(safe_adb(VLWatMicPM_vr(M,N6,N5,N4),VLSoilMicP_vr(N6,N5,N4)))

  IF(THETW1_vr(N3,N2,N1).GT.SoilWatAirDry_vr(N3,N2,N1) .AND. THETW1_vr(N6,N5,N4).GT.SoilWatAirDry_vr(N6,N5,N4) &
    .AND. VLWatMicPM_vr(M,N3,N2,N1).GT.ZEROS2(N2,N1) .AND. VLWatMicPM_vr(M,N6,N5,N4).GT.ZEROS2(N5,N4))THEN

    !------------------------------------------------
    !     MICROPORE CONCENTRATIONS FROM WATER-FILLED POROSITY
    !     IN CURRENT AND ADJACENT GRID CELLS
    !
    D9810: DO K=1,jcplx
      do idom=idom_beg,idom_end
        !source
        CDOM_MicP1(idom,K)=AZMAX1(DOM_MicP2_vr(idom,K,N3,N2,N1)/VLWatMicPM_vr(M,N3,N2,N1))
        !dest
        CDOM_MicP2_vr(idom,K)=AZMAX1(DOM_MicP2_vr(idom,K,N6,N5,N4)/VLWatMicPM_vr(M,N6,N5,N4))
      enddo
    ENDDO D9810

    DO idg=idg_beg,idg_NH3-1
      trcsolc1(idg)=AZMAX1(trcs_solml2_vr(idg,N3,N2,N1)/VLWatMicPM_vr(M,N3,N2,N1))
    ENDDO

    VLWatMicP4A = VLWatMicPM_vr(M,N3,N2,N1)*trcs_VLN_vr(ids_NH4,N3,N2,N1)
    IF(VLWatMicP4A.GT.ZEROS2(N2,N1))THEN
      trcsolc1(ids_NH4) = AMIN1(AZMAX1(trcs_solml2_vr(ids_NH4,N3,N2,N1)/VLWatMicP4A),tracerSolc_max(ids_NH4))
      trcsolc1(idg_NH3) = AMIN1(AZMAX1(trcs_solml2_vr(idg_NH3,N3,N2,N1)/VLWatMicP4A),tracerSolc_max(idg_NH3))
    ELSE
      trcsolc1(ids_NH4)=0.0_r8
      trcsolc1(idg_NH3)=0.0_r8
    ENDIF

    VLWatMicP3A = VLWatMicPM_vr(M,N3,N2,N1)*trcs_VLN_vr(ids_NO3,N3,N2,N1)
    IF(VLWatMicP3A.GT.ZEROS2(N2,N1))THEN
      trcsolc1(ids_NO3) = AMIN1(AZMAX1(trcs_solml2_vr(ids_NO3,N3,N2,N1)/VLWatMicP3A),tracerSolc_max(ids_NO3))
      trcsolc1(ids_NO2) = AZMAX1(trcs_solml2_vr(ids_NO2,N3,N2,N1)/VLWatMicP3A)      
    ELSE
      trcsolc1(ids_NO3)=0.0_r8
      trcsolc1(ids_NO2)=0.0_r8
    ENDIF

    VLWatMicP2A = VLWatMicPM_vr(M,N3,N2,N1)*trcs_VLN_vr(ids_H1PO4,N3,N2,N1)
    IF(VLWatMicP2A.GT.ZEROS2(N2,N1))THEN
      trcsolc1(ids_H1PO4) = AZMAX1(trcs_solml2_vr(ids_H1PO4,N3,N2,N1)/VLWatMicP2A)
      trcsolc1(ids_H2PO4) = AZMAX1(trcs_solml2_vr(ids_H2PO4,N3,N2,N1)/VLWatMicP2A)
    ELSE
      trcsolc1(ids_H1PO4)=0.0_r8
      trcsolc1(ids_H2PO4)=0.0_r8
    ENDIF

    VLWatMicP4B = VLWatMicPM_vr(M,N3,N2,N1)*trcs_VLN_vr(ids_NH4B,N3,N2,N1)
    IF(VLWatMicP4B.GT.ZEROS2(N2,N1))THEN
      trcsolc1(ids_NH4B) = AMIN1(AZMAX1(trcs_solml2_vr(ids_NH4B,N3,N2,N1)/VLWatMicP4B),tracerSolc_max(ids_NH4B))
      trcsolc1(idg_NH3B) = AMIN1(AZMAX1(trcs_solml2_vr(idg_NH3B,N3,N2,N1)/VLWatMicP4B),tracerSolc_max(idg_NH3B))
    ELSE
      trcsolc1(ids_NH4B)=0.0_r8
      trcsolc1(idg_NH3B)=0.0_r8
    ENDIF

    VLWatMicP3B = VLWatMicPM_vr(M,N3,N2,N1)*trcs_VLN_vr(ids_NO3B,N3,N2,N1)
    IF(VLWatMicP3B.GT.ZEROS2(N2,N1))THEN
      trcsolc1(ids_NO3B) = AMIN1(AZMAX1(trcs_solml2_vr(ids_NO3B,N3,N2,N1)/VLWatMicP3B),tracerSolc_max(ids_NO3B))
      trcsolc1(ids_NO2B) = AZMAX1(trcs_solml2_vr(ids_NO2B,N3,N2,N1)/VLWatMicP3B)
    ELSE
      trcsolc1(ids_NO3B) = trcsolc1(ids_NO3)
      trcsolc1(ids_NO2B) = trcsolc1(ids_NO2)
    ENDIF

    VLWatMicP2B = VLWatMicPM_vr(M,N3,N2,N1)*trcs_VLN_vr(ids_H1PO4B,N3,N2,N1)
    IF(VLWatMicP2B.GT.ZEROS2(N2,N1))THEN
      trcsolc1(ids_H1PO4B) = AZMAX1(trcs_solml2_vr(ids_H1PO4B,N3,N2,N1)/VLWatMicP2B)
      trcsolc1(ids_H2PO4B) = AZMAX1(trcs_solml2_vr(ids_H2PO4B,N3,N2,N1)/VLWatMicP2B)
    ELSE
      trcsolc1(ids_H1PO4B)=trcsolc1(ids_H1PO4)
      trcsolc1(ids_H2PO4B)=trcsolc1(ids_H2PO4)
    ENDIF

    DO idg=idg_beg,idg_NH3-1
      trcsolc2(idg)=AZMAX1(trcs_solml2_vr(idg,N6,N5,N4)/VLWatMicPM_vr(M,N6,N5,N4))
    ENDDO

    IF(VLWatMicPMA_vr(N6,N5,N4).GT.ZEROS2(N5,N4))THEN
      !gN/m3 H2O, the maximum solubility is about 9.7gN/100g water at 25oC    
      trcsolc2(idg_NH3) = AMIN1(AZMAX1(trcs_solml2_vr(idg_NH3,N6,N5,N4)/VLWatMicPMA_vr(N6,N5,N4)),tracerSolc_max(idg_NH3))
      trcsolc2(ids_NH4) = AMIN1(AZMAX1(trcs_solml2_vr(ids_NH4,N6,N5,N4)/VLWatMicPMA_vr(N6,N5,N4)),tracerSolc_max(ids_NH4))
    ELSE
      trcsolc2(idg_NH3) = 0.0_r8
      trcsolc2(ids_NH4) = 0.0_r8
    ENDIF

    VLWatMicPOA = VLWatMicPM_vr(M,N6,N5,N4)*trcs_VLN_vr(ids_NO3,N6,N5,N4)
    IF(VLWatMicPOA.GT.ZEROS2(N5,N4))THEN
      trcsolc2(ids_NO3) = AZMAX1(trcs_solml2_vr(ids_NO3,N6,N5,N4)/VLWatMicPOA)
      trcsolc2(ids_NO2) = AZMAX1(trcs_solml2_vr(ids_NO2,N6,N5,N4)/VLWatMicPOA)
    ELSE
      trcsolc2(ids_NO3) = 0.0_r8
      trcsolc2(ids_NO2) = 0.0_r8
    ENDIF

    VLWatMicPPA = VLWatMicPM_vr(M,N6,N5,N4)*trcs_VLN_vr(ids_H1PO4,N6,N5,N4)
    IF(VLWatMicPPA.GT.ZEROS2(N5,N4))THEN
      trcsolc2(ids_H1PO4) = AZMAX1(trcs_solml2_vr(ids_H1PO4,N6,N5,N4)/VLWatMicPPA)
      trcsolc2(ids_H2PO4) = AZMAX1(trcs_solml2_vr(ids_H2PO4,N6,N5,N4)/VLWatMicPPA)
    ELSE
      trcsolc2(ids_H1PO4) = 0.0_r8
      trcsolc2(ids_H2PO4) = 0.0_r8
    ENDIF

    IF(VLWatMicPMB_vr(N6,N5,N4).GT.ZEROS2(N5,N4))THEN
      trcsolc2(idg_NH3B) = AMIN1(AZMAX1(trcs_solml2_vr(idg_NH3B,N6,N5,N4)/VLWatMicPMB_vr(N6,N5,N4)),tracerSolc_max(idg_NH3B))
      trcsolc2(ids_NH4B) = AMIN1(AZMAX1(trcs_solml2_vr(ids_NH4B,N6,N5,N4)/VLWatMicPMB_vr(N6,N5,N4)),tracerSolc_max(ids_NH4B))
    ELSE
      trcsolc2(idg_NH3B) = trcsolc2(idg_NH3)
      trcsolc2(ids_NH4B) = trcsolc2(ids_NH4)
    ENDIF

    VLWatMicPOB = VLWatMicPM_vr(M,N6,N5,N4)*trcs_VLN_vr(ids_NO3B,N6,N5,N4)
    IF(VLWatMicPOB.GT.ZEROS2(N5,N4))THEN
      trcsolc2(ids_NO3B) = AMIN1(AZMAX1(trcs_solml2_vr(ids_NO3B,N6,N5,N4)/VLWatMicPOB),tracerSolc_max(ids_NO3B))
      trcsolc2(ids_NO2B) = AZMAX1(trcs_solml2_vr(ids_NO2B,N6,N5,N4)/VLWatMicPOB)
    ELSE
      trcsolc2(ids_NO3B) = trcsolc2(ids_NO3)
      trcsolc2(ids_NO2B) = trcsolc2(ids_NO2)
    ENDIF

    VLWatMicPPB = VLWatMicPM_vr(M,N6,N5,N4)*trcs_VLN_vr(ids_H1PO4B,N6,N5,N4)
    IF(VLWatMicPPB.GT.ZEROS2(N5,N4))THEN
      trcsolc2(ids_H1PO4B) = AZMAX1(trcs_solml2_vr(ids_H1PO4B,N6,N5,N4)/VLWatMicPPB)
      trcsolc2(ids_H2PO4B) = AZMAX1(trcs_solml2_vr(ids_H2PO4B,N6,N5,N4)/VLWatMicPPB)
    ELSE
      trcsolc2(ids_H1PO4B) = trcsolc2(ids_H1PO4)
      trcsolc2(ids_H2PO4B) = trcsolc2(ids_H2PO4)
    ENDIF
    !------------------------------------------------
    !     DIFFUSIVITIES IN CURRENT AND ADJACENT GRID CELL MICROPORES
    !
    DLYR1 = AMAX1(ZERO2,DLYR_3D(N,N3,N2,N1))
    DLYR2 = AMAX1(ZERO2,DLYR_3D(N,N6,N5,N4))
    TORTL = (TortMicPM_vr(M,N3,N2,N1)*DLYR1+TortMicPM_vr(M,N6,N5,N4)*DLYR2)/(DLYR1+DLYR2)
    DISPN = DISP_3D(N,N6,N5,N4)*AMIN1(VFLWX,ABS(WaterFlow2MicPM_3D(M,N,N6,N5,N4)/AREA(N,N6,N5,N4)))

    DO idom=idom_beg,idom_end
      DOMDifc(idom)=(DOM_diffusivitytscaledM_vr(idom,N6,N5,N4)*TORTL+DISPN)*XDPTH_3D(N,N6,N5,N4)
    ENDDO

    DO ids=ids_beg,ids_end
      SDifc(ids)=(SoluteDifusivitytscaledM_vr(ids,N6,N5,N4)*TORTL+DISPN)*XDPTH_3D(N,N6,N5,N4)
    ENDDO
    !
    D9805: DO K=1,jcplx
      do idom=idom_beg,idom_end
        DOM_Difus_Mac2Micp_flxM(idom,K)=DOMDifc(idom)*(CDOM_MicP1(idom,K)-CDOM_MicP2_vr(idom,K))        
      enddo
    ENDDO D9805

    DO idg=idg_beg,idg_NH3-1
      SDifFlx(idg)=SDifc(idg)*(trcsolc1(idg)-trcsolc2(idg))
    ENDDO

    DO ids=ids_nuts_beg,ids_nuts_end
      SDifFlx(ids)=SDifc(ids)*(trcsolc1(ids)-trcsolc2(ids))*AMIN1(trcs_VLN_vr(ids,N3,N2,N1),trcs_VLN_vr(ids,N6,N5,N4))      
    ENDDO
  ELSE
    D9905: DO K=1,jcplx
      DOM_Difus_Mac2Micp_flxM(idom_beg:idom_end,K)=0.0_r8
    ENDDO D9905
    SDifFlx(ids_beg:ids_end)=0._r8
  ENDIF

  DO ids=ids_beg,ids_end
    trcs_MicpTranspFlxM_3D(ids,N,N6,N5,N4)=trcs_MicpTranspFlxM_3D(ids,N,N6,N5,N4)+SDifFlx(ids)

    if(abs(trcs_MicpTranspFlxM_3D(ids,N,N6,N5,N4))>1.e10)then
      write(*,*)ids,N,N6,N5,N4,trcs_MicpTranspFlxM_3D(ids,N,N6,N5,N4),SDifFlx(ids),trcs_names(ids)   
      write(*,*)SDifc(ids),trcsolc1(ids),trcsolc2(ids),AMIN1(trcs_VLN_vr(ids,N3,N2,N1),trcs_VLN_vr(ids,N6,N5,N4))   
      call endrun(trim(mod_filename)//' at line',__LINE__)
    endif
  ENDDO

  DO K=1,jcplx
    do idom=idom_beg,idom_end
      DOM_MicpTranspFlxM_3D(idom,K,N,N6,N5,N4) = DOM_MicpTranspFlxM_3D(idom,K,N,N6,N5,N4)+DOM_Difus_Mac2Micp_flxM(idom,K)
    ENDDO  
  ENDDO
  call PrintInfo('end '//subname)
  end subroutine MicroporeSoluteDiffusionM

! ----------------------------------------------------------------------
  subroutine MacroporeSoluteAdvectionM(M,N,N1,N2,N3,N4,N5,N6)
  implicit none
  integer , intent(in) :: M,N,N1,N2,N3,N4,N5,N6
  real(r8) :: DOM_Adv2MacP_flxM(idom_beg:idom_end,1:jcplx)
  character(len=*), parameter :: subname='MacroporeSoluteAdvectionM'
  real(r8) ::  trcs_Adv2MacP_flxM(ids_beg:ids_end)
  integer  :: K,idg,ids,idom
  real(r8) :: VFLW

  call PrintInfo('beg '//subname)
  !     WaterFlow2MacPM_3D=water flux through soil macropore from watsub.f
  !
  IF(WaterFlow2MacPM_3D(M,N,N6,N5,N4).GT.0.0_r8)THEN
    !
    !     IF MACROPORE WATER FLUX FROM 'WATSUB' IS FROM CURRENT TO
    !     ADJACENT GRID CELL THEN CONVECTIVE TRANSPORT IS THE PRODUCT
    !     OF WATER FLUX AND MACROPORE SOLUTE CONCENTRATIONS IN CURRENT
    !     GRID CELL
    !
    IF(VLWatMacPM_vr(M,N3,N2,N1).GT.ZEROS2(N2,N1))THEN
      VFLW=AZMAX1(AMIN1(VFLWX,WaterFlow2MacPM_3D(M,N,N6,N5,N4)/VLWatMacPM_vr(M,N3,N2,N1)))
    ELSE
      VFLW=VFLWX
    ENDIF
    !
    !     ACCOUNT FOR MACROPORE-MICROPORE EXCHANGE
    !
    IF(N.EQ.iVerticalDirection .AND. VLMacP_vr(N6,N5,N4).GT.VLWatMacPM_vr(M,N6,N5,N4))THEN
      D9800: DO K=1,jcplx
        do idom=idom_beg,idom_end
          DOM_Adv2MacP_flxM(idom,K)=VFLW*AZMAX1((DOM_MacP2_vr(idom,K,N3,N2,N1) &
            -AZMIN1(DOM_Mac2MicPore_flxM_vr(idom,K,NU(N2,N1),N2,N1))))
        enddo
      ENDDO D9800

      DO idg=idg_beg,idg_NH3-1
         trcs_Adv2MacP_flxM(idg)=VFLW*AZMAX1((trcs_soHml2_vr(idg,N3,N2,N1)-AZMIN1(trcs_Mac2MicPore_flxM_vr(idg,NU(N2,N1),N2,N1))))
      ENDDO

      DO ids=ids_nuts_beg,ids_nuts_end
         trcs_Adv2MacP_flxM(ids)=VFLW*AZMAX1((trcs_soHml2_vr(ids,N3,N2,N1) &
          -AZMIN1(trcs_Mac2MicPore_flxM_vr(ids,NU(N2,N1),N2,N1)*trcs_VLN_vr(ids,N3,N2,N1)))) &
          *trcs_VLN_vr(ids,N6,N5,N4)
      ENDDO
!
!     OTHERWISE
!
    ELSE
      D9850: DO K=1,jcplx
        do idom=idom_beg,idom_end
          DOM_Adv2MacP_flxM(idom,K)=VFLW*AZMAX1(DOM_MacP2_vr(idom,K,N3,N2,N1))
        enddo
      ENDDO D9850
!exclude NH3 and NH3B
      DO idg=idg_beg,idg_end-2
         trcs_Adv2MacP_flxM(idg)=VFLW*AZMAX1(trcs_soHml2_vr(idg,N3,N2,N1))
      ENDDO

      DO ids=ids_nuts_beg,ids_nuts_end
         trcs_Adv2MacP_flxM(ids)=VFLW*AZMAX1(trcs_soHml2_vr(ids,N3,N2,N1))*trcs_VLN_vr(ids,N6,N5,N4)
      ENDDO
    ENDIF
!
!     IF MACROPORE WATER FLUX FROM 'WATSUB' IS FROM ADJACENT TO
!     CURRENT GRID CELL THEN CONVECTIVE TRANSPORT IS THE PRODUCT
!     OF WATER FLUX AND MACROPORE SOLUTE CONCENTRATIONS IN ADJACENT
!     GRID CELL
!
  ELSEIF(WaterFlow2MacPM_3D(M,N,N6,N5,N4).LT.0.0_r8)THEN
    IF(VLWatMacPM_vr(M,N6,N5,N4).GT.ZEROS2(N5,N4))THEN
      VFLW=AZMIN1(AMAX1(-VFLWX,WaterFlow2MacPM_3D(M,N,N6,N5,N4)/VLWatMacPM_vr(M,N6,N5,N4)))
    ELSE
      VFLW=-VFLWX
    ENDIF

    D9665: DO K=1,jcplx
      do idom=idom_beg,idom_end
        DOM_Adv2MacP_flxM(idom,K)=VFLW*AZMAX1(DOM_MacP2_vr(idom,K,N6,N5,N4))
      enddo
    ENDDO D9665

    DO idg=idg_beg,idg_end-2
       trcs_Adv2MacP_flxM(idg)=VFLW*AZMAX1(trcs_soHml2_vr(idg,N6,N5,N4))
    ENDDO

    DO ids=ids_nuts_beg,ids_nuts_end
       trcs_Adv2MacP_flxM(ids)=VFLW*AZMAX1(trcs_soHml2_vr(ids,N6,N5,N4))*trcs_VLN_vr(ids,N6,N5,N4)
    ENDDO
  ELSE
!
!     NO MACROPORE FLUX
!
    D9795: DO K=1,jcplx
      DOM_Adv2MacP_flxM(idom_beg:idom_end,K)=0.0_r8
    ENDDO D9795
     trcs_Adv2MacP_flxM(ids_beg:ids_end)=0.0_r8
  ENDIF

  DO ids=ids_beg,ids_end
    trcs_MacpTranspFlxM_3D(ids,N,N6,N5,N4)= trcs_MacpTranspFlxM_3D(ids,N,N6,N5,N4)+trcs_Adv2MacP_flxM(ids)
  ENDDO

  DO K=1,jcplx
    do idom=idom_beg,idom_end
      DOM_MacpTranspFlxM_3D(idom,K,N,N6,N5,N4) =DOM_MacpTranspFlxM_3D(idom,K,N,N6,N5,N4)+ DOM_Adv2MacP_flxM(idom,K)
    ENDDO
  ENDDO      
  call PrintInfo('end '//subname)
  end subroutine MacroporeSoluteAdvectionM

! ----------------------------------------------------------------------
  subroutine MacroporeSoluteDispersionM(M,N,N1,N2,N3,N4,N5,N6)
  implicit none
  integer , intent(in) :: M,N,N1,N2,N3,N4,N5,N6

  character(len=*), parameter :: subname='MacroporeSoluteDispersionM'
  real(r8) :: DOM_Difus_Macp_flxM(idom_beg:idom_end,1:jcplx)
  real(r8) :: VOLH2A,VOLH2B,VOLH3A,VOLH3B,VOLH4A,VOLH4B
  real(r8) :: VOLHMA,VOLHMB,VOLHOA,VOLHOB,VOLHPA,VOLHPB
  real(r8) :: trcs_coH1(ids_beg:ids_end)
  real(r8) :: trcs_coH2(ids_beg:ids_end)
  real(r8) :: SDifc(ids_beg:ids_end),TORTL
  real(r8) :: DOMDifc(idom_beg:idom_end)
  real(r8) :: DISPN,DLYR1,DLYR2
  real(r8) :: SDifHFlx(ids_beg:ids_end)
  integer  :: K,ids,idg,idom
  real(r8) :: CDOM_MacP1(idom_beg:idom_end,1:jcplx)
  real(r8) :: CDOM_MacP2_vr(idom_beg:idom_end,1:jcplx)

  call PrintInfo('beg '//subname)

  IF(VLWatMacPM_vr(M,N3,N2,N1).GT.SoilWatAirDry_vr(N3,N2,N1)*VLMacP_vr(N3,N2,N1) &
    .AND. VLWatMacPM_vr(M,N6,N5,N4).GT.SoilWatAirDry_vr(N6,N5,N4)*VLMacP_vr(N6,N5,N4))THEN
    !------------------------------------------------------------------------------------
    !     MACROPORE CONCENTRATIONS IN CURRENT AND ADJACENT GRID CELLS
    !
    D9790: DO K=1,jcplx
      do idom=idom_beg,idom_end
        CDOM_MacP1(idom,K)=AZMAX1(DOM_MacP2_vr(idom,K,N3,N2,N1)/VLWatMacPM_vr(M,N3,N2,N1))
        CDOM_MacP2_vr(idom,K)=AZMAX1(DOM_MacP2_vr(idom,K,N6,N5,N4)/VLWatMacPM_vr(M,N6,N5,N4))
      enddo
    ENDDO D9790

    !exclude NH3 and NH3B
    DO idg=idg_beg,idg_NH3-1
      trcs_coH1(idg)=AZMAX1(trcs_soHml2_vr(idg,N3,N2,N1)/VLWatMacPM_vr(M,N3,N2,N1))
    ENDDO

    VOLH4A = VLWatMacPM_vr(M,N3,N2,N1)*trcs_VLN_vr(ids_NH4,N3,N2,N1)
    IF(VOLH4A.GT.ZEROS2(N2,N1))THEN
      trcs_coH1(ids_NH4)=AZMAX1(trcs_soHml2_vr(ids_NH4,N3,N2,N1)/VOLH4A)
      trcs_coH1(idg_NH3)=AZMAX1(trcs_soHml2_vr(idg_NH3,N3,N2,N1)/VOLH4A)
    ELSE
      trcs_coH1(ids_NH4)=0.0_r8
      trcs_coH1(idg_NH3)=0.0_r8
    ENDIF

    VOLH3A = VLWatMacPM_vr(M,N3,N2,N1)*trcs_VLN_vr(ids_NO3,N3,N2,N1)
    IF(VOLH3A.GT.ZEROS2(N2,N1))THEN
      trcs_coH1(ids_NO3)=AZMAX1(trcs_soHml2_vr(ids_NO3,N3,N2,N1)/VOLH3A)
      trcs_coH1(ids_NO2)=AZMAX1(trcs_soHml2_vr(ids_NO2,N3,N2,N1)/VOLH3A)
    ELSE
      trcs_coH1(ids_NO3)=0.0_r8
      trcs_coH1(ids_NO2)=0.0_r8
    ENDIF

    VOLH2A = VLWatMacPM_vr(M,N3,N2,N1)*trcs_VLN_vr(ids_H1PO4,N3,N2,N1)
    IF(VOLH2A.GT.ZEROS2(N2,N1))THEN
      trcs_coH1(ids_H1PO4)=AZMAX1(trcs_soHml2_vr(ids_H1PO4,N3,N2,N1)/VOLH2A)
      trcs_coH1(ids_H2PO4)=AZMAX1(trcs_soHml2_vr(ids_H2PO4,N3,N2,N1)/VOLH2A)
    ELSE
      trcs_coH1(ids_H1PO4)=0.0_r8
      trcs_coH1(ids_H2PO4)=0.0_r8
    ENDIF

    VOLH4B = VLWatMacPM_vr(M,N3,N2,N1)*trcs_VLN_vr(ids_NH4B,N3,N2,N1)
    IF(VOLH4B.GT.ZEROS2(N2,N1))THEN
      trcs_coH1(ids_NH4B)=AZMAX1(trcs_soHml2_vr(ids_NH4B,N3,N2,N1)/VOLH4B)
      trcs_coH1(idg_NH3B)=AZMAX1(trcs_soHml2_vr(idg_NH3B,N3,N2,N1)/VOLH4B)
    ELSE
      trcs_coH1(ids_NH4B)=trcs_coH1(ids_NH4)
      trcs_coH1(idg_NH3B)=trcs_coH1(idg_NH3)
    ENDIF

    VOLH3B = VLWatMacPM_vr(M,N3,N2,N1)*trcs_VLN_vr(ids_NO3B,N3,N2,N1)
    IF(VOLH3B.GT.ZEROS2(N2,N1))THEN
      trcs_coH1(ids_NO3B)=AZMAX1(trcs_soHml2_vr(ids_NO3B,N3,N2,N1)/VOLH3B)
      trcs_coH1(ids_NO2B)=AZMAX1(trcs_soHml2_vr(ids_NO2B,N3,N2,N1)/VOLH3B)
    ELSE
      trcs_coH1(ids_NO3B)=trcs_coH1(ids_NO3)
      trcs_coH1(ids_NO2B)=trcs_coH1(ids_NO2)
    ENDIF

    VOLH2B = VLWatMacPM_vr(M,N3,N2,N1)*trcs_VLN_vr(ids_H1PO4B,N3,N2,N1)
    IF(VOLH2B.GT.ZEROS2(N2,N1))THEN
      trcs_coH1(ids_H1PO4B)=AZMAX1(trcs_soHml2_vr(ids_H1PO4B,N3,N2,N1)/VOLH2B)
      trcs_coH1(ids_H2PO4B)=AZMAX1(trcs_soHml2_vr(ids_H2PO4B,N3,N2,N1)/VOLH2B)
    ELSE
      trcs_coH1(ids_H1PO4B)=trcs_coH1(ids_H1PO4)
      trcs_coH1(ids_H2PO4B)=trcs_coH1(ids_H2PO4)
    ENDIF

   !excldue NH3 and NH3B
    DO idg=idg_beg,idg_NH3-1
      trcs_coH2(idg)=AZMAX1(trcs_soHml2_vr(idg,N6,N5,N4)/VLWatMacPM_vr(M,N6,N5,N4))
    ENDDO

    VOLHMA=VLWatMacPM_vr(M,N6,N5,N4)*trcs_VLN_vr(ids_NH4,N6,N5,N4)
    IF(VOLHMA.GT.ZEROS2(N5,N4))THEN
      trcs_coH2(ids_NH4)=AZMAX1(trcs_soHml2_vr(ids_NH4,N6,N5,N4)/VOLHMA)
      trcs_coH2(idg_NH3)=AZMAX1(trcs_soHml2_vr(idg_NH3,N6,N5,N4)/VOLHMA)
    ELSE
      trcs_coH2(ids_NH4)=0.0_r8
      trcs_coH2(idg_NH3)=0.0_r8
    ENDIF

    VOLHOA = VLWatMacPM_vr(M,N6,N5,N4)*trcs_VLN_vr(ids_NO3,N6,N5,N4)
    IF(VOLHOA.GT.ZEROS2(N5,N4))THEN
      trcs_coH2(ids_NO3)=AZMAX1(trcs_soHml2_vr(ids_NO3,N6,N5,N4)/VOLHOA)
      trcs_coH2(ids_NO2)=AZMAX1(trcs_soHml2_vr(ids_NO2,N6,N5,N4)/VOLHOA)
    ELSE
      trcs_coH2(ids_NO3)=0.0_r8
      trcs_coH2(ids_NO2)=0.0_r8
    ENDIF

    VOLHPA = VLWatMacPM_vr(M,N6,N5,N4)*trcs_VLN_vr(ids_H1PO4,N6,N5,N4)
    IF(VOLHPA.GT.ZEROS2(N5,N4))THEN
      trcs_coH2(ids_H1PO4)=AZMAX1(trcs_soHml2_vr(ids_H1PO4,N6,N5,N4)/VOLHPA)
      trcs_coH2(ids_H2PO4)=AZMAX1(trcs_soHml2_vr(ids_H2PO4,N6,N5,N4)/VOLHPA)
    ELSE
      trcs_coH2(ids_H1PO4)=0.0_r8
      trcs_coH2(ids_H2PO4)=0.0_r8
    ENDIF

    VOLHMB=VLWatMacPM_vr(M,N6,N5,N4)*trcs_VLN_vr(ids_NH4B,N6,N5,N4)
    IF(VOLHMB.GT.ZEROS2(N5,N4))THEN
      trcs_coH2(ids_NH4B)=AZMAX1(trcs_soHml2_vr(ids_NH4B,N6,N5,N4)/VOLHMB)
      trcs_coH2(idg_NH3B)=AZMAX1(trcs_soHml2_vr(idg_NH3B,N6,N5,N4)/VOLHMB)
    ELSE
      trcs_coH2(ids_NH4B)=trcs_coH2(ids_NH4)
      trcs_coH2(idg_NH3B)=trcs_coH2(idg_NH3)
    ENDIF

    VOLHOB = VLWatMacPM_vr(M,N6,N5,N4)*trcs_VLN_vr(ids_NO3B,N6,N5,N4)
    IF(VOLHOB.GT.ZEROS2(N5,N4))THEN
      trcs_coH2(ids_NO3B)=AZMAX1(trcs_soHml2_vr(ids_NO3B,N6,N5,N4)/VOLHOB)
      trcs_coH2(ids_NO2B)=AZMAX1(trcs_soHml2_vr(ids_NO2B,N6,N5,N4)/VOLHOB)
    ELSE
      trcs_coH2(ids_NO3B)=trcs_coH2(ids_NO3)
      trcs_coH2(ids_NO2B)=trcs_coH2(ids_NO2)
    ENDIF

    VOLHPB = VLWatMacPM_vr(M,N6,N5,N4)*trcs_VLN_vr(ids_H1PO4B,N6,N5,N4)
    IF(VOLHPB.GT.ZEROS2(N5,N4))THEN
      trcs_coH2(ids_H1PO4B)=AZMAX1(trcs_soHml2_vr(ids_H1PO4B,N6,N5,N4)/VOLHPB)
      trcs_coH2(ids_H2PO4B)=AZMAX1(trcs_soHml2_vr(ids_H2PO4B,N6,N5,N4)/VOLHPB)
    ELSE
      trcs_coH2(ids_H1PO4B)=trcs_coH2(ids_H1PO4)
      trcs_coH2(ids_H2PO4B)=trcs_coH2(ids_H2PO4)
    ENDIF

    !------------------------------------------------------------------------------------
    !
    !     DIFFUSIVITIES IN CURRENT AND ADJACENT GRID CELL MACROPORES
    !
    DLYR1 = AMAX1(ZERO2,DLYR_3D(N,N3,N2,N1))
    DLYR2 = AMAX1(ZERO2,DLYR_3D(N,N6,N5,N4))
    TORTL = (TortMacPM_vr(M,N3,N2,N1)*DLYR1+TortMacPM_vr(M,N6,N5,N4)*DLYR2)/(DLYR1+DLYR2)
    DISPN = DISP_3D(N,N6,N5,N4)*AMIN1(VFLWX,ABS(WaterFlow2MacPM_3D(M,N,N6,N5,N4)/AREA(N,N6,N5,N4)))

    DO idom=idom_beg,idom_end
      DOMDifc(idom)=(DOM_diffusivitytscaledM_vr(idom,N6,N5,N4)*TORTL+DISPN)*XDPTH_3D(N,N6,N5,N4)
    ENDDO

    DO ids=ids_beg,ids_end
      SDifc(ids)=(SoluteDifusivitytscaledM_vr(ids,N6,N5,N4)*TORTL+DISPN)*XDPTH_3D(N,N6,N5,N4)
    ENDDO
    !
    !     DIFFUSIVE FLUXES BETWEEN CURRENT AND ADJACENT GRID CELL
    !     MACROPORES
    !
    D9785: DO K=1,jcplx
      do idom=idom_beg,idom_end
        DOM_Difus_Macp_flxM(idom,K)=DOMDifc(idom)*(CDOM_MacP1(idom,K)-CDOM_MacP2_vr(idom,K))
      enddo
    ENDDO D9785

    ! exclude NH3 and NH3B
    DO idg=idg_beg,idg_NH3-1
      SDifHFlx(idg)=SDifc(idg)*(trcs_coH1(idg)-trcs_coH2(idg))
    ENDDO

    DO ids=ids_nuts_beg,ids_end
      SDifHFlx(ids)=SDifc(ids)*(trcs_coH1(ids)-trcs_coH2(ids)) &
        *AMIN1(trcs_VLN_vr(ids,N3,N2,N1),trcs_VLN_vr(ids,N6,N5,N4))
    ENDDO
  ELSE
    D9780: DO K=1,jcplx
      DOM_Difus_Macp_flxM(idom_beg:idom_end,K)=0.0_r8
    ENDDO D9780
    SDifHFlx(ids_beg:ids_end)=0.0_r8
  ENDIF

  DO ids = ids_beg,ids_end
    trcs_MacpTranspFlxM_3D(ids,N,N6,N5,N4)=trcs_MacpTranspFlxM_3D(ids,N,N6,N5,N4)+SDifHFlx(ids)
  ENDDO

  DO K=1,jcplx
    do idom=idom_beg,idom_end
      DOM_MacpTranspFlxM_3D(idom,K,N,N6,N5,N4) = DOM_MacpTranspFlxM_3D(idom,K,N,N6,N5,N4)+DOM_Difus_Macp_flxM(idom,K)
    enddo
  ENDDO

  call PrintInfo('end '//subname)
  end subroutine MacroporeSoluteDispersionM

!----------------------------------------------------------------------

  subroutine MicroporeSoluteAdvectionM(M,N,N1,N2,N3,N4,N5,N6)
  implicit none
  integer, intent(in) :: M,N,N1,N2,N3,N4,N5,N6

  character(len=*), parameter :: subname='MicroporeSoluteAdvectionM'
  real(r8) :: DOM_Adv2MicP_flx(idom_beg:idom_end,1:jcplx)
  real(r8) :: trcs_Adv2MicP_flx(ids_beg:ids_end)
  integer :: K,ids,idom
  real(r8) :: VFLW
   
  call PrintInfo('beg '//subname)
  IF(WaterFlow2MicPM_3D(M,N,N6,N5,N4).GT.0.0_r8)THEN
    IF(VLWatMicPM_vr(M,N3,N2,N1).GT.ZEROS2(N2,N1))THEN
      VFLW=AZMAX1(AMIN1(VFLWX,WaterFlow2MicPM_3D(M,N,N6,N5,N4)/VLWatMicPM_vr(M,N3,N2,N1)))
    ELSE
      VFLW=VFLWX
    ENDIF
    DO  K=1,jcplx
      do idom=idom_beg,idom_end
        DOM_Adv2MicP_flx(idom,K)=VFLW*AZMAX1(DOM_MicP2_vr(idom,K,N3,N2,N1))
      enddo
    ENDDO

    DO ids=ids_beg,ids_end
      trcs_Adv2MicP_flx(ids)=VFLW*AZMAX1(trcs_solml2_vr(ids,N3,N2,N1))
    ENDDO
    !
  ELSE
    IF(VLWatMicPM_vr(M,N6,N5,N4).GT.ZEROS2(N5,N4))THEN
      VFLW=AZMIN1(AMAX1(-VFLWX,WaterFlow2MicPM_3D(M,N,N6,N5,N4)/VLWatMicPM_vr(M,N6,N5,N4)))
    ELSE
      VFLW=-VFLWX
    ENDIF
    D9815: DO K=1,jcplx
      do idom=idom_beg,idom_end
        DOM_Adv2MicP_flx(idom,K)=VFLW*AZMAX1(DOM_MicP2_vr(idom,K,N6,N5,N4))
      enddo
    ENDDO D9815
    DO ids=ids_beg,ids_end
      trcs_Adv2MicP_flx(ids)=VFLW*AZMAX1(trcs_solml2_vr(ids,N6,N5,N4))
    ENDDO
  ENDIF

  DO ids=ids_beg,ids_end
    trcs_MicpTranspFlxM_3D(ids,N,N6,N5,N4)=trcs_MicpTranspFlxM_3D(ids,N,N6,N5,N4)+trcs_Adv2MicP_flx(ids)

    if(abs(trcs_MicpTranspFlxM_3D(ids,N,N6,N5,N4))>1.e10)then
      write(*,*)ids,N,N6,N5,N4,trcs_Adv2MicP_flx(ids)
      write(*,*)VFLW,trcs_solml2_vr(ids,N6,N5,N4),trcs_solml2_vr(ids,N3,N2,N1),trcs_names(ids)
      write(*,*)trcs_solml_vr(ids,N6,N5,N4),trcs_solml_vr(ids,N3,N2,N1)
      call endrun(trim(mod_filename)//' at line',__LINE__)
    endif
  ENDDO

  DO K=1,jcplx
    do idom=idom_beg,idom_end
      DOM_MicpTranspFlxM_3D(idom,K,N,N6,N5,N4) = DOM_MicpTranspFlxM_3D(idom,K,N,N6,N5,N4)+DOM_Adv2MicP_flx(idom,K)
    enddo
  ENDDO
  call PrintInfo('end '//subname)
  end subroutine MicroporeSoluteAdvectionM

! ----------------------------------------------------------------------
  subroutine MicMacPoresSoluteExchM(M,N1,N2,N3)
  !
  !Description:
  !     MACROPORE-MICROPORE SOLUTE EXCHANGE WITHIN SOIL
  !     LAYER FROM WATER EXCHANGE IN 'WATSUB' AND
  !     FROM MACROPORE OR MICROPORE SOLUTE CONCENTRATIONS
  !
  implicit none
  integer, intent(in) :: M,N1,N2,N3
  character(len=*), parameter :: subname='MicMacPoresSoluteExchM'
  integer :: K,ids,idg,idom

  call PrintInfo('beg '//subname)

  do idg=idg_beg,idg_NH3-1
    trcg_VLWatMicP_vr(idg,N3,N2,N1)=VLWatMicPM_vr(M,N3,N2,N1)*GasSolbility_vr(idg,N3,N2,N1)
  enddo

  trcg_VLWatMicP_vr(idg_NH3,N3,N2,N1)=VLWatMicPMA_vr(N3,N2,N1)*GasSolbility_vr(idg_NH3,N3,N2,N1)
  trcg_VLWatMicP_vr(idg_NH3B,N3,N2,N1)=VLWatMicPMB_vr(N3,N2,N1)*GasSolbility_vr(idg_NH3,N3,N2,N1)

  call MicMacPoresSoluteAdvExchM(M,N1,N2,N3)
 
  call MicMacPoresSoluteDifExchM(M,N1,N2,N3)
  
  call PrintInfo('end '//subname)

  end subroutine MicMacPoresSoluteExchM

! ----------------------------------------------------------------------
  subroutine MicMacPoresSoluteAdvExchM(M,N4,N5,N6)
  !
  !Description
  !Do micropore-macropore tracer exchange by upstream advection
  implicit none
  integer, intent(in) :: M,N4,N5,N6

  character(len=*), parameter :: subname='MicMacPoresSoluteAdvExchM'
  real(r8) :: DOM_Adv2MicP_flx(idom_beg:idom_end,1:jcplx)
  real(r8) :: VFLW    !fraction in flow
  real(r8) :: trcs_Adv2MicP_flx(ids_beg:ids_end)
  integer :: K,idg,ids,idom

  call PrintInfo('beg '//subname)
!
!     MACROPORE TO MICROPORE TRANSFER
!
  IF(FWatExMacP2MicPM_vr(M,N6,N5,N4).GT.0.0_r8)THEN
    IF(VLWatMacPM_vr(M,N6,N5,N4).GT.ZEROS2(N5,N4))THEN
      VFLW=AZMAX1(AMIN1(VFLWX,FWatExMacP2MicPM_vr(M,N6,N5,N4)/VLWatMacPM_vr(M,N6,N5,N4)))
    ELSE
      VFLW=VFLWX
    ENDIF

    D9970: DO K=1,jcplx
      do idom=idom_beg,idom_end
        DOM_Adv2MicP_flx(idom,K)=VFLW*AZMAX1(DOM_MacP2_vr(idom,K,N6,N5,N4))
      enddo
    ENDDO D9970

    DO idg=idg_beg,idg_NH3-1
      trcs_Adv2MicP_flx(idg)=VFLW*AZMAX1(trcs_soHml2_vr(idg,N6,N5,N4))
    ENDDO

    DO ids=ids_nuts_beg,ids_nuts_end
      trcs_Adv2MicP_flx(ids)=VFLW*AZMAX1(trcs_soHml2_vr(ids,N6,N5,N4))*trcs_VLN_vr(ids,N6,N5,N4)
    ENDDO
!
!     MICROPORE TO MACROPORE TRANSFER
!
  ELSEIF(FWatExMacP2MicPM_vr(M,N6,N5,N4).LT.0.0_r8)THEN
    IF(VLWatMicPM_vr(M,N6,N5,N4).GT.ZEROS2(N5,N4))THEN
      VFLW=AZMIN1(AMAX1(-VFLWX,FWatExMacP2MicPM_vr(M,N6,N5,N4)/VLWatMicPM_vr(M,N6,N5,N4)))
    ELSE
      VFLW=-VFLWX
    ENDIF

    D9965: DO K=1,jcplx
      do idom=idom_beg,idom_end
        DOM_Adv2MicP_flx(idom,K)=VFLW*AZMAX1(DOM_MicP2_vr(idom,K,N6,N5,N4))
      enddo
    ENDDO D9965
    !
    DO idg=idg_beg,idg_NH3-1
      trcs_Adv2MicP_flx(idg)=VFLW*AZMAX1(trcs_solml2_vr(idg,N6,N5,N4))
    ENDDO

    !include NH3 and NH3B below
    DO ids=ids_nuts_beg,ids_nuts_end
      trcs_Adv2MicP_flx(ids)=VFLW*AZMAX1(trcs_solml2_vr(ids,N6,N5,N4))*trcs_VLN_vr(ids,N6,N5,N4)
    ENDDO
    !
    !     NO MACROPORE TO MICROPORE TRANSFER
    !
  ELSE
    D9960: DO K=1,jcplx
      DOM_Adv2MicP_flx(idom_beg:idom_end,K)=0.0_r8
    ENDDO D9960
    trcs_Adv2MicP_flx(ids_beg:ids_end)=0.0_r8
  ENDIF

  DO ids=ids_beg,ids_end
    trcs_Mac2MicPore_flxM_vr(ids,N6,N5,N4)=trcs_Mac2MicPore_flxM_vr(ids,N6,N5,N4)+trcs_Adv2MicP_flx(ids)
  ENDDO
  !
  !     TOTAL CONVECTIVE TRANSFER BETWEEN MACROPOES AND MICROPORES
  !
  DO  K=1,jcplx
    do idom=idom_beg,idom_end
      DOM_Mac2MicPore_flxM_vr(idom,K,N6,N5,N4)=DOM_Mac2MicPore_flxM_vr(idom,K,N6,N5,N4)+DOM_Adv2MicP_flx(idom,K)
    enddo
  enddo
  call PrintInfo('end '//subname)
  end subroutine MicMacPoresSoluteAdvExchM

! ----------------------------------------------------------------------

  subroutine MicMacPoresSoluteDifExchM(M,N4,N5,N6)
  !
  !Description
  !Do macropore-micropore tracer exchange by diffusion
  !
  implicit none
  integer, intent(in) :: M,N4,N5,N6

  character(len=*), parameter :: subname='MicMacPoresSoluteDifExchM'
  real(r8) :: trcs_Difus_Mac2Micp_flxM(ids_beg:ids_end)
  real(r8) :: VLWatMacPS,VOLWT
  integer  :: K,ids,idg,idom
  real(r8) :: DOM_Difus_Mac2Micp_flxM(idom_beg:idom_end,1:jcplx)

  call PrintInfo('beg '//subname)
  IF(VLWatMacPM_vr(M,N6,N5,N4).GT.ZEROS2(N5,N4))THEN    
    VLWatMacPS = AMIN1(XFRS*VGeomLayer_vr(N6,N5,N4),VLWatMacPM_vr(M,N6,N5,N4))
    VOLWT      = VLWatMicPM_vr(M,N6,N5,N4)+VLWatMacPS

    D9955: DO K=1,jcplx
      do idom=idom_beg,idom_end
        DOM_Difus_Mac2Micp_flxM(idom,K)=dts_HeatWatTP*(AZMAX1(DOM_MacP2_vr(idom,K,N6,N5,N4)) &
          *VLWatMicPM_vr(M,N6,N5,N4)-AZMAX1(DOM_MicP2_vr(idom,K,N6,N5,N4))*VLWatMacPS)/VOLWT
      enddo
    ENDDO D9955

    DO idg=idg_beg,idg_NH3-1
      trcs_Difus_Mac2Micp_flxM(idg)=dts_HeatWatTP*(AZMAX1(trcs_soHml2_vr(idg,N6,N5,N4))*VLWatMicPM_vr(M,N6,N5,N4) &
        -AZMAX1(trcs_solml2_vr(idg,N6,N5,N4))*VLWatMacPS)/VOLWT
    ENDDO

    DO ids=ids_nuts_beg,ids_nuts_end
      trcs_Difus_Mac2Micp_flxM(ids)=dts_HeatWatTP*(AZMAX1(trcs_soHml2_vr(ids,N6,N5,N4))*VLWatMicPM_vr(M,N6,N5,N4) &
        -AZMAX1(trcs_solml2_vr(ids,N6,N5,N4))*VLWatMacPS)/VOLWT &
        *trcs_VLN_vr(ids,N6,N5,N4)
    ENDDO
  ELSE
    D9975: DO K=1,jcplx
      DOM_Difus_Mac2Micp_flxM(idom_beg:idom_end,K)=0.0_r8
    ENDDO D9975
    trcs_Difus_Mac2Micp_flxM(ids_beg:ids_end)=0.0_r8
  ENDIF

  DO ids=ids_beg,ids_end
    trcs_Mac2MicPore_flxM_vr(ids,N6,N5,N4)=trcs_Mac2MicPore_flxM_vr(ids,N6,N5,N4)+trcs_Difus_Mac2Micp_flxM(ids)
  ENDDO
  !
  ! TOTAL CONVECTIVE +DIFFUSIVE TRANSFER BETWEEN MACROPOES AND MICROPORES
  !
  DO  K=1,jcplx
    do idom=idom_beg,idom_end
      DOM_Mac2MicPore_flxM_vr(idom,K,N6,N5,N4)=DOM_Mac2MicPore_flxM_vr(idom,K,N6,N5,N4)+DOM_Difus_Mac2Micp_flxM(idom,K)
    enddo
  enddo
  call PrintInfo('end '//subname)
  end subroutine MicMacPoresSoluteDifExchM

! ----------------------------------------------------------------------
  subroutine BubbleEffluxM(I,J,M,NHE,NHW,NVS,NVN)
  implicit none
  integer, intent(in) :: I,J,M,NHE,NHW,NVS,NVN

  integer  :: N1,N2,N3
  integer  :: iDisableEbu                        !bubbling flag

  character(len=*), parameter :: subname='BubbleEffluxM'  
  real(r8) :: THETW1
  real(r8) :: GasMassSolubility(idg_beg:idg_end)
  real(r8) :: trcg_VOLG(idg_beg:idg_end)   !gas concentation corresponding to each volatile [mol d-2]
  integer  :: idg
  real(r8) :: VTATM,VTGAS,DVTGAS
  real(r8) :: FracEbu                      !fraction removed through ebullition
  real(r8) :: dPond
  integer  :: LG

  !
  call PrintInfo('beg '//subname)
  LG=0
  DO N1=NHW,NHE
    DO N2=NVN,NVS
      iDisableEbu=ifalse
      dPond=0._r8
      DO N3=NUM(N2,N1),NL(N2,N1)  !sweep downward untill reaching layer where ebullition is disabled
        !
        THETW1=AZMAX1(safe_adb(VLWatMicPM_vr(M,N3,N2,N1),VLSoilMicP_vr(N3,N2,N1)))
        !
        IF(THETW1.GT.SoilWatAirDry_vr(N3,N2,N1) & ! has significant water 
          .AND. iDisableEbu.EQ.ifalse)THEN        ! ebullition is allowed

          do idg=idg_beg,idg_NH3
            GasMassSolubility(idg) =MolecularWeight(idg)*GasSolbility_vr(idg,N3,N2,N1)  !conver into carbon g /mol
          enddo        
          GasMassSolubility(idg_NH3B)=GasMassSolubility(idg_NH3)
          !
          !     GASEOUS EQUIVALENT PARTIAL CONCENTRATIONS
          !
          DO idg=idg_beg,idg_end
            trcg_VOLG(idg)=AZMAX1(trcs_solml2_vr(idg,N3,N2,N1))/GasMassSolubility(idg)  !mol/d2 
          ENDDO 
          !
          !     GASEOUS EQUIVALENT ATMOSPHERIC CONCENTRATION
          !
          !     VTATM=molar gas concentration at atmospheric pressure
          !     VTGAS=total molar gas concentration, should the gaseous phase be included as well?     
          !    
          if(iPondFlag_col(N2,N1))then 
            if(SoilBulkDensity_vr(N3,N2,N1).LE.ZERO)then
              !still in the ponding water, \rho*g*h,  1.e3*10
              dPond=(CumDepz2LayBottom_vr(N3,N2,N1)-CumDepz2LayBottom_vr(NUM(N2,N1)-1,N2,N1))*10._r8
            else
              !in the soil
              dPond=(CumDepz2LayBottom_vr(iPondBotLev_col(N2,N1),N2,N1)-CumDepz2LayBottom_vr(NUM(N2,N1)-1,N2,N1))*10._r8
            endif  
          endif  
          VTATM = AZMAX1((PBOT_col(N2,N1)+dPond)*VLWatMicPM_vr(M,N3,N2,N1)/(RGasC*TKS_vr(N3,N2,N1)))*1.E3_r8  !mol gas/ d2
          VTGAS = sum(trcg_VOLG(idg_beg:idg_end))
          !
          !     PROPORTIONAL REMOVAL OF EXCESS AQUEOUS GASES
          !
          IF(VTGAS.GT.VTATM)THEN
            DVTGAS  = 0.5_r8*(VTATM-VTGAS)         !<0, bubbling occurs
            FracEbu = AZMIN1(DVTGAS/VTGAS)
            DO idg=idg_beg,idg_end
              trcg_Ebu_flxM_vr(idg,N3,N2,N1)=FracEbu*trcg_VOLG(idg)*GasMassSolubility(idg)  

              trcs_solml2_vr(idg,N3,N2,N1)=trcs_solml2_vr(idg,N3,N2,N1)+trcg_Ebu_flxM_vr(idg,N3,N2,N1)
              if(LG==0)then
                !accumulate bubbling flux
                if(idg==idg_NH3B)then
                  trcg_ebu_flx_col(idg_NH3,N2,N1) = trcg_ebu_flx_col(idg_NH3,N2,N1)+trcg_Ebu_flxM_vr(idg,N3,N2,N1)
                else
                  trcg_ebu_flx_col(idg,N2,N1) = trcg_ebu_flx_col(idg,N2,N1)+trcg_Ebu_flxM_vr(idg,N3,N2,N1)                
                endif
              else
                 trcg_gasml_vr(idg,N3,N2,N1)=trcg_gasml_vr(idg,N3,N2,N1)-trcg_Ebu_flxM_vr(idg,N3,N2,N1)
              endif              
            ENDDO
            RO2AquaSourcePrev_vr(N3,N2,N1) = RO2AquaSourcePrev_vr(N3,N2,N1)+trcg_Ebu_flxM_vr(idg_O2,N3,N2,N1)
            RCH4PhysexchPrev_vr(N3,N2,N1)  = RCH4PhysexchPrev_vr(N3,N2,N1)+trcg_Ebu_flxM_vr(idg_CH4,N3,N2,N1)
          ELSE
            LG=N3
            trcg_Ebu_flxM_vr(idg_beg:idg_end,N3,N2,N1)=0.0_r8
          ENDIF
        ELSE
          iDisableEbu=itrue
          trcg_Ebu_flxM_vr(idg_beg:idg_end,N3,N2,N1)=0.0_r8
        ENDIF
      ENDDO  
    ENDDO
  ENDDO  


  call PrintInfo('end '//subname)
  end subroutine BubbleEffluxM
! ----------------------------------------------------------------------
  subroutine TracerThruRofSnowXYM(M,N,N1,N2,N4,N5,N4B,N5B)
  !
  !Tracer flux through surface runoff and 
  !snow drift
  implicit none
  integer, intent(in) :: M,N
  integer, intent(in) :: N1,N2,N4,N5,N4B,N5B

  character(len=*), parameter :: subname='TracerThruRofSnowXYM'
  integer :: NN,K,idg,ids,idom

  call PrintInfo('beg '//subname)
  DO NN=1,2
    DO  K=1,jcplx
      do idom=idom_beg,idom_end
        DOM_SurfRunoff_flxM(idom,K,N2,N1)=DOM_SurfRunoff_flxM(idom,K,N2,N1)+DOM_FloXSurRof_flxM_2DH(idom,K,N,NN,N2,N1)
      enddo
    enddo

    DO idg=idg_beg,idg_NH3
      trcg_SurfRunoff_flxM(idg,N2,N1)=trcg_SurfRunoff_flxM(idg,N2,N1)+trcg_SurRof_flxM_2DH(idg,N,NN,N2,N1)
    ENDDO
    
    DO ids=ids_nut_beg,ids_nuts_end
      trcn_SurfRunoff_flxM(ids,N2,N1)=trcn_SurfRunoff_flxM(ids,N2,N1)+trcn_SurRof_flxM_2DH(ids,N,NN,N2,N1)
    ENDDO

    IF(IFLBM_2DH(M,N,NN,N5,N4).EQ.0)THEN
      D9551: DO  K=1,jcplx
        do idom=idom_beg,idom_end
          DOM_SurfRunoff_flxM(idom,K,N2,N1)=DOM_SurfRunoff_flxM(idom,K,N2,N1)-DOM_FloXSurRof_flxM_2DH(idom,K,N,NN,N5,N4)
        enddo
      enddo D9551

      DO idg=idg_beg,idg_NH3
        trcg_SurfRunoff_flxM(idg,N2,N1)=trcg_SurfRunoff_flxM(idg,N2,N1)-trcg_SurRof_flxM_2DH(idg,N,NN,N5,N4)
      ENDDO

      DO ids=ids_nut_beg,ids_nuts_end
        trcn_SurfRunoff_flxM(ids,N2,N1)=trcn_SurfRunoff_flxM(ids,N2,N1)-trcn_SurRof_flxM_2DH(ids,N,NN,N5,N4)
      ENDDO
    ENDIF

    IF(N4B.GT.0 .AND. N5B.GT.0 .AND. NN.EQ.iFront)THEN
      DO  K=1,jcplx
        do idom=idom_beg,idom_end
          DOM_SurfRunoff_flxM(idom,K,N2,N1)=DOM_SurfRunoff_flxM(idom,K,N2,N1)-DOM_FloXSurRof_flxM_2DH(idom,K,N,NN,N5B,N4B)
        enddo
      enddo
      DO idg=idg_beg,idg_NH3
        trcg_SurfRunoff_flxM(idg,N2,N1)=trcg_SurfRunoff_flxM(idg,N2,N1)-trcg_SurRof_flxM_2DH(idg,N,NN,N5B,N4B)
      ENDDO

      DO ids=ids_nut_beg,ids_nuts_end
        trcn_SurfRunoff_flxM(ids,N2,N1)=trcn_SurfRunoff_flxM(ids,N2,N1)-trcn_SurRof_flxM_2DH(ids,N,NN,N5B,N4B)
      ENDDO
    ENDIF
  enddo
  !---------------
  !     NET OVERLAND SOLUTE FLUX IN SNOW
  !
  DO idg=idg_beg,idg_NH3
    trcg_SnowDrift_flxM(idg,N2,N1)=trcg_SnowDrift_flxM(idg,N2,N1)+trcg_SnowDrift_flxM_2D(idg,N,NN,N2,N1)
  ENDDO

  do ids=ids_nut_beg,ids_nuts_end
    trcn_SnowDrift_flxM(ids,N2,N1)=trcn_SnowDrift_flxM(ids,N2,N1)+trcn_SnowDrift_flxM_2D(ids,N,NN,N2,N1)
  enddo

  if(IFLBSM_2DH(M,N,NN,N5,N4).EQ.0)then
    DO idg=idg_beg,idg_NH3
      trcg_SnowDrift_flxM(idg,N2,N1)=trcg_SnowDrift_flxM(idg,N2,N1)-trcg_SnowDrift_flxM_2D(idg,N,NN,N5,N4)
    ENDDO

    do ids=ids_nut_beg,ids_nuts_end
      trcn_SnowDrift_flxM(ids,N2,N1)=trcn_SnowDrift_flxM(ids,N2,N1)-trcn_SnowDrift_flxM_2D(ids,N,NN,N5,N4)
    enddo
  endif

  IF(N4B.GT.0.AND.N5B.GT.0.AND.NN.EQ.iFront)THEN
    DO idg=idg_beg,idg_NH3
      trcg_SnowDrift_flxM(idg,N2,N1)=trcg_SnowDrift_flxM(idg,N2,N1)-trcg_SnowDrift_flxM_2D(idg,N,NN,N5B,N4B)
    ENDDO

    do ids=ids_nut_beg,ids_nuts_end
      trcn_SnowDrift_flxM(ids,N2,N1)=trcn_SnowDrift_flxM(ids,N2,N1)-trcn_SnowDrift_flxM_2D(ids,N,NN,N5B,N4B)
    enddo
  ENDIF
  call PrintInfo('end '//subname)
  end subroutine TracerThruRofSnowXYM

!------------------------------------------------------------------------------------------

  subroutine TracerFluxThruSnowZM(I,J,M,N1,N2)
  !
  !tracer flux from snow to litter and soil
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: M,N1,N2

  character(len=*), parameter :: subname='TracerFluxThruSnowZM'
  integer :: LS,LS2
  integer :: idg,idn

  call PrintInfo('beg '//subname)
  DO  LS=1,JS

    IF(VLSnowHeatCapM_snvr(M,LS,N2,N1).GT.VLHeatCapSnowMin_col(N2,N1)*1.e-4_r8)THEN
      LS2=MIN(JS,LS+1)

      ! NEXT LAYER IS IN THE SNOWPACK   
      IF(LS.LT.JS .AND. VLSnowHeatCapM_snvr(M,LS2,N2,N1).GT.VLHeatCapSnowMin_col(N2,N1))THEN
        DO idg=idg_beg,idg_NH3
          trcg_Aqua_flxM_snvr(idg,LS,N2,N1)=trcg_Aqua_flxM_snvr(idg,LS,N2,N1) &
            +trcg_AquaAdv_flxM_snvr(idg,LS,N2,N1)-trcg_AquaAdv_flxM_snvr(idg,LS2,N2,N1)
        ENDDO

        DO idn=ids_nut_beg,ids_nuts_end
          trcn_Aqua_flxM_snvr(idn,LS,N2,N1)=trcn_Aqua_flxM_snvr(idn,LS,N2,N1) &
            +trcn_AquaAdv_flxM_snvr(idn,LS,N2,N1)-trcn_AquaAdv_flxM_snvr(idn,LS2,N2,N1)
        ENDDO

        ! NEXT LAYER IS THE LITTER AND SOIL SURFACE
      ELSE

        ! exclude NH3 and NH3B
        DO idg=idg_beg,idg_NH3
          trcg_Aqua_flxM_snvr(idg,LS,N2,N1)=trcg_Aqua_flxM_snvr(idg,LS,N2,N1) &
            +trcg_AquaAdv_flxM_snvr(idg,LS,N2,N1)   &
            -trcg_AquaADV_Snow2Litr_flxM(idg,N2,N1) &
            -trcg_AquaADV_Snow2Soil_flxM(idg,N2,N1)
        ENDDO

        trcg_Aqua_flxM_snvr(idg_NH3,LS,N2,N1)=trcg_Aqua_flxM_snvr(idg_NH3,LS,N2,N1) &
          -trcg_AquaADV_Snow2Soil_flxM(idg_NH3B,N2,N1)

        !from NH4 to H2PO4
        DO idn=0,ids_nuts
          trcn_Aqua_flxM_snvr(ids_NH4+idn,LS,N2,N1)=trcn_Aqua_flxM_snvr(ids_NH4+idn,LS,N2,N1) &
            +trcn_AquaAdv_flxM_snvr(ids_NH4+idn,LS,N2,N1)   &
            -trcn_AquaADV_Snow2Litr_flxM(ids_NH4+idn,N2,N1) &
            -trcn_AquaADV_Snow2Soil_flxM(ids_NH4+idn,N2,N1) &
            -trcn_AquaADV_Snow2Band_flxM(ids_NH4B+idn,N2,N1)
        ENDDO
      ENDIF
    ENDIF
  enddo
  call PrintInfo('end '//subname)
  end subroutine TracerFluxThruSnowZM
! ----------------------------------------------------------------------

  subroutine XBoundaryTracerTranspXYZM(M,N,NN,M1,M2,M3,M4,M5,M6)
  !
  !Description:
  !tracer advection across the boundaries in 3 different directions
  !
  implicit none
  integer, intent(in) :: M,N,NN
  integer, intent(in) :: M1,M2,M3
  integer, intent(in) :: M4,M5,M6    !note when NN==iBehind, (M3,M2,M1)==(M6,M5,M4)

  character(len=*), parameter :: subname='XBoundaryTracerTranspXYZM'
  integer :: K,idg,idn,ids,idom
  real(r8) :: VFLW

  call PrintInfo('beg '//subname)
  !------------------
  !Micropore
  !out of cell  (M3,M2,M1), horizontal flow + bottom, exclude top
  !
  IF(NN.EQ.iFront .AND. WaterFlow2MicPM_3D(M,N,M6,M5,M4).GT.0.0_r8 &  
    .OR. (NN.EQ.iBehind .AND. WaterFlow2MicPM_3D(M,N,M6,M5,M4).LT.0.0_r8))THEN
    IF(VLWatMicPM_vr(M,M3,M2,M1).GT.ZEROS2(M2,M1))THEN
      VFLW=AMAX1(-VFLWX,AMIN1(VFLWX,WaterFlow2MicPM_3D(M,N,M6,M5,M4)/VLWatMicPM_vr(M,M3,M2,M1)))
    ELSE
      VFLW=0.0_r8
    ENDIF

    DO  K=1,jcplx
      do idom=idom_beg,idom_end
        DOM_MicpTranspFlxM_3D(idom,K,N,M6,M5,M4)=VFLW*AZMAX1(DOM_MicP2_vr(idom,K,M3,M2,M1))
      enddo
    enddo

    DO idn=ids_beg,ids_end
      trcs_MicpTranspFlxM_3D(idn,N,M6,M5,M4)=VFLW*AZMAX1(trcs_solml2_vr(idn,M3,M2,M1))*trcs_VLN_vr(idn,M3,M2,M1)

      if(abs(trcs_MicpTranspFlxM_3D(idn,N,M6,M5,M4))>1.e3)then
        write(*,*)VFLW,trcs_solml2_vr(idn,M3,M2,M1),trcs_VLN_vr(idn,M3,M2,M1)
        call endrun(trim(mod_filename)//' at line',__LINE__)                
      endif
    ENDDO
    !
    !  SOLUTE GAIN WITH SUBSURFACE MICROPORE WATER GAIN
    !coming into cell (N3,N2,N1) from the boundaries
  ELSE
    DO  K=1,jcplx
      DOM_MicpTranspFlxM_3D(idom_beg:idom_end,K,N,M6,M5,M4)=0.0_r8
    enddo
    trcs_MicpTranspFlxM_3D(idg_beg:idg_NH3-1,N,M6,M5,M4)=0.0_r8

    !Assuming the nutrient tracer concentations follow that from irrigation
    !This seems not that reasonable when irrigation is not included as part of WaterFlow2MicPM_3D
    !double check  
    DO ids=ids_nuts_beg,ids_nuts_end
      trcs_MicpTranspFlxM_3D(ids,N,M6,M5,M4)=WaterFlow2MicPM_3D(M,N,M6,M5,M4)*trcn_irrig_vr(ids,M3,M2,M1)*trcs_VLN_vr(ids,M3,M2,M1)
      if(abs(trcs_MicpTranspFlxM_3D(ids,N,M6,M5,M4))>1.e3)then
        write(*,*)WaterFlow2MicPM_3D(M,N,M6,M5,M4), &
          trcn_irrig_vr(ids,M3,M2,M1),trcs_VLN_vr(ids,M3,M2,M1)
        call endrun(trim(mod_filename)//' at line',__LINE__)                
      endif        
    ENDDO
  ENDIF

  !---------------------------
  !     SOLUTE LOSS WITH SUBSURFACE MACROPORE WATER LOSS
  !
  IF(NN.EQ.iFront .AND. WaterFlow2MacPM_3D(M,N,M6,M5,M4).GT.0.0_r8 &
    .OR. (NN.EQ.iBehind .AND. WaterFlow2MacPM_3D(M,N,M6,M5,M4).LT.0.0_r8))THEN
    IF(VLWatMacPM_vr(M,M3,M2,M1).GT.ZEROS2(M2,M1))THEN
      VFLW=AMAX1(-VFLWX,AMIN1(VFLWX,WaterFlow2MacPM_3D(M,N,M6,M5,M4)/VLWatMacPM_vr(M,M3,M2,M1)))
    ELSE
      VFLW=0.0_r8
    ENDIF

    DO  K=1,jcplx
      do idom=idom_beg,idom_end
        DOM_MacpTranspFlxM_3D(idom,K,N,M6,M5,M4)=VFLW*AZMAX1(DOM_MacP2_vr(idom,K,M3,M2,M1))
      enddo
    enddo

    DO idn=ids_beg,ids_end
      trcs_MacpTranspFlxM_3D(idn,N,M6,M5,M4)=VFLW*AZMAX1(trcs_soHml2_vr(idn,M3,M2,M1))*trcs_VLN_vr(idn,M3,M2,M1)
    ENDDO
    !
    !     NO SOLUTE GAIN IN SUBSURFACE MACROPORES
    !
  ELSE
    DO  K=1,jcplx
      DOM_MacpTranspFlxM_3D(idom_beg:idom_end,K,N,M6,M5,M4)=0.0_r8
    enddo
    trcs_MacpTranspFlxM_3D(ids_beg:ids_end,N,M6,M5,M4)=0.0_r8
  ENDIF
  call PrintInfo('end '//subname)
  end subroutine XBoundaryTracerTranspXYZM
!------------------------------------------------------------------------------------------

  subroutine XSurfBoundaryTracerRunoffXYM(M,N,NN,N1,N2,M4,M5,RCHQF)

  !Description
  !
  implicit none
  integer, intent(in) :: M,N,NN,N1,N2,M4,M5
  real(r8), intent(in) :: RCHQF
  character(len=*), parameter :: subname='XSurfBoundaryTracerRunoffXYM'
  real(r8) :: FQRM
  integer :: K,idg,idn,ids,idom
  logical :: LNOFlow

  call PrintInfo('beg '//subname)
  LNOFlow=.not.XGridRunoffFlag_2DH(NN,N,N2,N1) .OR. isclose(RCHQF,0.0_r8) .OR. SurfRunoffPotentM_col(M,N2,N1).LE.ZEROS(N2,N1)

  IF(LNOFlow)THEN
    DO  K=1,jcplx
      DOM_FloXSurRof_flxM_2DH(idom_beg:idom_end,K,N,NN,M5,M4)=0.0_r8
    ENDDO
    trcg_SurRof_flxM_2DH(idg_beg:idg_NH3,N,NN,M5,M4)          = 0.0_r8
    trcn_SurRof_flxM_2DH(ids_nut_beg:ids_nuts_end,N,NN,M5,M4) = 0.0_r8
  ELSE
    !
    ! SOLUTE LOSS FROM RUNOFF DEPENDING ON ASPECT
    ! AND BOUNDARY CONDITIONS SET IN SITE FILE
    ! Lose tracer to external environment
    IF((NN.EQ.iFront .AND. QflxSurfRunoffM_2DH(M,N,NN,M5,M4).GT.ZEROS(N2,N1)) &          
      .OR. (NN.EQ.iBehind .AND. QflxSurfRunoffM_2DH(M,N,NN,M5,M4).LT.ZEROS(N2,N1)))THEN  

      FQRM=QflxSurfRunoffM_2DH(M,N,NN,M5,M4)/SurfRunoffPotentM_col(M,N2,N1)
      DO  K=1,jcplx
        DO idom=idom_beg,idom_end
          DOM_FloXSurRof_flxM_2DH(idom,K,N,NN,M5,M4)=DOM_FloXSurRunoff_PotFlxM(idom,K,N2,N1)*FQRM
        ENDDO
      enddo

      DO idg=idg_beg,idg_NH3
        trcg_SurRof_flxM_2DH(idg,N,NN,M5,M4)=trcg_FloXSurRunoff_PotFlxM(idg,N2,N1)*FQRM
      ENDDO

      DO ids=ids_nut_beg,ids_nuts_end
        trcn_SurRof_flxM_2DH(ids,N,NN,M5,M4)=trcn_FloXSurRunoff_PotFlxM(ids,N2,N1)*FQRM
      ENDDO
      !
      !     SOLUTE GAIN FROM RUNON DEPENDING ON ASPECT
      !     AND BOUNDARY CONDITIONS SET IN SITE FILE
      !
    ELSEIF((NN.EQ.iBehind .AND. QflxSurfRunoffM_2DH(M,N,NN,M5,M4).GT.ZEROS(N2,N1)) &
      .OR.(NN.EQ.iFront .AND. QflxSurfRunoffM_2DH(M,N,NN,M5,M4).LT.ZEROS(N2,N1)))THEN

      !For incoming flux, assume the DOM concentation is zero 
      DO  K=1,jcplx
        DOM_FloXSurRof_flxM_2DH(idom_beg:idom_end,K,N,NN,M5,M4)=0.0_r8
      enddo

      DO idg=idg_beg,idg_NH3
        trcg_SurRof_flxM_2DH(idg,N,NN,M5,M4) = QflxSurfRunoffM_2DH(M,N,NN,M5,M4)*trcg_rain_mole_conc_col(idg,M5,M4)*MolecularWeight(idg)
      ENDDO

      trcn_SurRof_flxM_2DH(ids_nut_beg:ids_nuts_end,N,NN,M5,M4)=0.0_r8

    ELSE
      DO  K=1,jcplx
        DOM_FloXSurRof_flxM_2DH(idom_beg:idom_end,K,N,NN,M5,M4)=0.0_r8
      enddo
      trcg_SurRof_flxM_2DH(idg_beg:idg_NH3,N,NN,M5,M4)          = 0.0_r8
      trcn_SurRof_flxM_2DH(ids_nut_beg:ids_nuts_end,N,NN,M5,M4) = 0.0_r8
    ENDIF

  ENDIF
  !
  !     BOUNDARY SNOW FLUX
  ! No snow-drift across the boundaries
  IF(NN.EQ.iFront)THEN
    trcg_SnowDrift_flxM_2D(idg_beg:idg_NH3,N,iBehind,M5,M4)          = 0.0_r8
    trcn_SnowDrift_flxM_2D(ids_nut_beg:ids_nuts_end,N,iBehind,M5,M4) = 0.0_r8
  ENDIF
  call PrintInfo('end '//subname)
  end subroutine XSurfBoundaryTracerRunoffXYM
!------------------------------------------------------------------------------------------

  subroutine TracerXBoundariesM(I,J,M,NHW,NHE,NVN,NVS)
  !
  !Description:
  !Do transport across the boundaries, upper surface, subsurface lateral,
  ! and bottom
  !
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: M,NHW, NHE, NVN, NVS

  character(len=*), parameter :: subname='TracerXBoundariesM'
  integer :: NY,NX,L
  integer :: NN,N
  integer :: N1,N2,N3,N4,N5,N6,N4B,N5B  !inner grids
  integer :: M1,M2,M3,M4,M5,M6          !boundary exchange
  real(r8) :: RCHQF

  call PrintInfo('beg '//subname)
!
  DO  NX=NHW,NHE
    DO  NY=NVN,NVS
      D9585: DO L=NU(NY,NX),NL(NY,NX)
        N1=NX;N2=NY;N3=L
        M1 = NX;M2 = NY;M3 = L             !source grid        
        D9580: DO  N=FlowDirIndicator_col(NY,NX),3
          D9575: DO  NN=1,2
            IF(N.EQ.iWestEastDirection)THEN
              N4  = NX+1; N5  = NY          
              N4B = NX-1;N5B = NY;N6  = L
              IF(NN.EQ.iFront)THEN                   !eastward                
                IF(NX.EQ.NHE)THEN                    !eastern boundary
                  M4 = NX+1;M5 = NY;M6 = L           !dest
                  RCHQF  = RechargEastSurf(M2,M1)
                ELSE
                  cycle
                ENDIF
              ELSEIF(NN.EQ.iBehind)THEN             !west                
                IF(NX.EQ.NHW)THEN                   !western boundary
                  M4 = NX;M5 = NY;M6 = L            !dest
                  RCHQF  = RechargWestSurf(M5,M4)
                ELSE
                  cycle
                ENDIF
              ENDIF
            ELSEIF(N.EQ.iNorthSouthDirection)THEN              
              N4  = NX;N5  = NY+1  
              N4B = NX;N5B = NY-1; N6  = L
              IF(NN.EQ.iFront)THEN  
                IF(NY.EQ.NVS)THEN                     ! southern boundary
                  M4 = NX;M5 = NY+1;M6 = L            !target grid
                  RCHQF  = RechargSouthSurf(M2,M1)
                ELSE
                  cycle
                ENDIF
              ELSEIF(NN.EQ.iBehind)THEN              !north                
                IF(NY.EQ.NVN)THEN                    !northern boundary
                  M4 = NX;M5 = NY;M6 = L             !target
                  RCHQF  = RechargNorthSurf(M5,M4)
                ELSE
                  cycle
                ENDIF
              ENDIF
            ELSEIF(N.EQ.iVerticalDirection)THEN      !vertical
              N4 = NX;N5 = NY;N6 = L+1               !target
              IF(NN.EQ.iFront)THEN                
                IF(L.EQ.NL(NY,NX))THEN               !lower boundary
                  M4 = NX;M5 = NY;M6 = L+1           !target grid
                ELSE
                  cycle
                ENDIF
              ELSEIF(NN.EQ.iBehind)THEN              
                cycle
              ENDIF
            ENDIF            
            !
            call XBoundaryTracerFlowM(L,N,NN,M,N1,N2,N3,M1,M2,M3,M4,M5,M6,RCHQF)
            !
          ENDDO D9575
          !
          IF(L.EQ.NUM(N2,N1))THEN
            IF(N.NE.iVerticalDirection)THEN
              !horizontal flow
              call TracerThruRofSnowXYM(M,N,N1,N2,N4,N5,N4B,N5B)
              !
            ELSEIF(N.EQ.iVerticalDirection)THEN
              ! vertical direction
              call TracerFluxThruSnowZM(I,J,M,N1,N2)
            ENDIF
          ENDIF

        ENDDO D9580
      ENDDO D9585

    ENDDO
  ENDDO
  call PrintInfo('end '//subname)
  end subroutine TracerXBoundariesM

! ----------------------------------------------------------------------
  subroutine XBoundaryTracerFlowM(L,N,NN,M,N1,N2,N3,M1,M2,M3,M4,M5,M6,RCHQF)
  implicit none
  !
  !Do cross boundary tracer transport, lateral+bottom, exclude top 
  !
  integer, intent(in) :: L,N, NN, M 
  integer, intent(in) :: N1, N2,N3
  integer, intent(in) :: M1, M2,M3,M4,M5,M6
  real(r8), intent(in):: RCHQF
  character(len=*), parameter :: subname='XBoundaryTracerFlowM'
  real(r8) :: FLGM,FQRM,VFLW
  integer :: K,idg

! begin_execution
  call PrintInfo('beg '//subname)

  !Surface/litter horizontal tracer flow
  IF(L.EQ.NUM(M2,M1) .AND. N.NE.iVerticalDirection)THEN        
    call XSurfBoundaryTracerRunoffXYM(M,N,NN,N1,N2,M4,M5,RCHQF)
  ENDIF
!
  IF(VLSoilPoreMicP_vr(N3,N2,N1).GT.ZEROS2(N2,N1))THEN
    IF(FlowDirIndicator_col(M2,M1).NE.3 .OR. N.EQ.iVerticalDirection)THEN        
      call XBoundaryTracerTranspXYZM(M,N,NN,M1,M2,M3,M4,M5,M6)
    ENDIF
  ENDIF

  call PrintInfo('end '//subname)
  end subroutine XBoundaryTracerFlowM

!------------------------------------------------------------------------------------------

  subroutine NetTracerFlowXSoilPoresM(NY,NX,N,M,N1,N2,N3,N4,N5,N6)
  !
  !Description
  !Vertical or lateral flux between grid grids
  !does not include litter layer
  implicit none

  integer, intent(in) :: NY,NX,N,M
  integer, intent(in) :: N1,N2,N3
  integer, intent(in) :: N4,N5,N6
  character(len=*), parameter :: subname='NetTracerFlowXSoilPoresM'
  integer :: K,LL,ids,idg,idom

  call PrintInfo('beg '//subname)
  
  DO  K=1,jcplx
    do idom=idom_beg,idom_end
      DOM_Transp2Micp_flxM_vr(idom,K,N3,N2,N1)=DOM_Transp2Micp_flxM_vr(idom,K,N3,N2,N1) &
        +DOM_MicpTranspFlxM_3D(idom,K,N,N3,N2,N1)-DOM_MicpTranspFlxM_3D(idom,K,N,N6,N5,N4)
      DOM_Transp2Macp_flxM_vr(idom,K,N3,N2,N1)=DOM_Transp2Macp_flxM_vr(idom,K,N3,N2,N1) &
        +DOM_MacpTranspFlxM_3D(idom,K,N,N3,N2,N1)-DOM_MacpTranspFlxM_3D(idom,K,N,N6,N5,N4)

      if(N.NE.iVerticalDirection)then    
        DOM_transpFlxM_2DH(idom,K,N2,N1)=DOM_transpFlxM_2DH(idom,K,N2,N1) &
          +DOM_MicpTranspFlxM_3D(idom,K,N,N3,N2,N1)-DOM_MicpTranspFlxM_3D(idom,K,N,N6,N5,N4) &
          +DOM_MacpTranspFlxM_3D(idom,K,N,N3,N2,N1)-DOM_MacpTranspFlxM_3D(idom,K,N,N6,N5,N4)
      endif    
    enddo
  enddo
  !
  !incoming flux (N3,N2,N1), outgoing flux (N6,N5,N4), N is direction
  !
  DO ids=ids_beg,ids_end
    trcs_Transp2Micp_flxM_vr(ids,N3,N2,N1)=trcs_Transp2Micp_flxM_vr(ids,N3,N2,N1) &
      +trcs_MicpTranspFlxM_3D(ids,N,N3,N2,N1)-trcs_MicpTranspFlxM_3D(ids,N,N6,N5,N4)

    trcs_Transp2Macp_flxM_vr(ids,N3,N2,N1)=trcs_Transp2Macp_flxM_vr(ids,N3,N2,N1) &
      +trcs_MacpTranspFlxM_3D(ids,N,N3,N2,N1)-trcs_MacpTranspFlxM_3D(ids,N,N6,N5,N4)

    if(N.NE.iVerticalDirection)then
       trcs_transpFlxM_2DH(ids,N2,N1)=trcs_transpFlxM_2DH(ids,N2,N1) &
         +trcs_MicpTranspFlxM_3D(ids,N,N3,N2,N1)-trcs_MicpTranspFlxM_3D(ids,N,N6,N5,N4) &
         +trcs_MacpTranspFlxM_3D(ids,N,N3,N2,N1)-trcs_MacpTranspFlxM_3D(ids,N,N6,N5,N4) 
    endif

    if(abs(trcs_Transp2Micp_flxM_vr(ids,N3,N2,N1))>1.e10)then
      write(*,*)N,N3,N2,N1,N6,N5,N4,trcs_names(ids)
      write(*,*)trcs_MicpTranspFlxM_3D(ids,N,N3,N2,N1),-trcs_MicpTranspFlxM_3D(ids,N,N6,N5,N4)
      call endrun(trim(mod_filename)//' at line',__LINE__)                
    endif  
  ENDDO

  call PrintInfo('end '//subname)

  end subroutine NetTracerFlowXSoilPoresM

!------------------------------------------------------------------------------------------

  subroutine UpdateSurfTracerM(I,J,M,NHW,NHE,NVN,NVS,dpscal,dpscal_dom,pscal1,pscal1_dom)
  !
  !Description
  !
  implicit none
  integer, intent(in) :: I,J,M,NHW,NHE,NVN,NVS
  real(r8),intent(in) :: dpscal(ids_beg:ids_end)
  real(r8),intent(in) :: dpscal_dom(idom_beg:idom_end)
  real(r8), optional, intent(inout) :: pscal1(ids_beg:ids_end)
  real(r8), optional, intent(inout) :: pscal1_dom(idom_beg:idom_end)

  character(len=*), parameter :: subname='UpdateSurfTracerM'

  real(r8) :: pscal(ids_beg:ids_end)
  real(r8) :: pscal_dom(idom_beg:idom_end)
  integer :: K,idg,ids,idom,L,NY,NX
  real(r8) :: flux

  call PrintInfo('beg '//subname)
  if(present(pscal1))then
    pscal=pscal1
  else
    pscal=1.001_r8
  endif  

  if(present(pscal1_dom))then
    pscal_dom=pscal1_dom
  else
    pscal_dom=pscal1_dom
  endif  

  DO NY=NHW,NHE
    DO  NX=NVN,NVS  
      !--------------------
      !to litter
      DO  K=1,jcplx
        DO idom=idom_beg,idom_end
          if(dpscal_dom(idom)>tiny_p)then
            flux=DOM_SurfRunoff_flxM(idom,K,NY,NX)+DOM_MicpTranspFlxM_3D(idom,K,3,0,NY,NX)
            flux=flux*dpscal_dom(idom)
            call get_flux_scalar(DOM_MicP2_vr(idom,K,0,NY,NX),flux,DOM_MicP_vr(idom,K,0,NY,NX),pscal_dom(idom))
          endif
        ENDDO
      ENDDO

     !exclude NH3B
      DO idg=idg_beg,idg_NH3-1
        if(dpscal(idg)>tiny_p)then
          flux = RGasAtmDisol2LitrM_col(idg,NY,NX)+trcg_SurfRunoff_flxM(idg,NY,NX)+trcs_MicpTranspFlxM_3D(ids,3,0,NY,NX)
          flux=flux*dpscal(idg)
          call get_flux_scalar(trcs_solml2_vr(idg,0,NY,NX),flux,trcs_solml_vr(idg,0,NY,NX),pscal(idg))
        endif
      ENDDO

      idg=idg_NH3
      if(dpscal(idg)>tiny_p)then
        flux = RGasAtmDisol2LitrM_col(idg,NY,NX)+trcg_SurfRunoff_flxM(idg,NY,NX)+trcs_MicpTranspFlxM_3D(idg,3,0,NY,NX) &
          -RBGCSinkSoluteM_vr(idg_NH3,0,NY,NX)
        flux=flux*dpscal(idg)
        call get_flux_scalar(trcs_solml2_vr(idg,0,NY,NX),flux,trcs_solml_vr(idg,0,NY,NX),pscal(idg))
      endif

      DO ids=ids_nut_beg,ids_nuts_end
        if(dpscal(ids)>tiny_p)then
          flux=trcs_MicpTranspFlxM_3D(ids,3,0,NY,NX)+trcn_SurfRunoff_flxM(ids,NY,NX)-RBGCSinkSoluteM_vr(ids,0,NY,NX)
          flux=flux*dpscal(ids)
          call get_flux_scalar(trcs_solml2_vr(ids,0,NY,NX),flux,trcs_solml_vr(ids,0,NY,NX),pscal(ids))
        endif
      ENDDO

      !-------------------------
      !topsoil
      L=NU(NY,NX)
      IF(VLSoilPoreMicP_vr(L,NY,NX).GT.ZEROS2(NY,NX))THEN            
        DO idg=idg_beg,idg_NH3-1
          if(dpscal(idg)>tiny_p)then
            flux=RGasAtmDisol2SoilM_col(idg,NY,NX)+trcsol_Irrig_flxM_vr(ids,L,NY,NX) & 
              +trcs_Transp2Micp_flxM_vr(ids,L,NY,NX) + trcs_Mac2MicPore_flxM_vr(ids,L,NY,NX)       
            flux=flux*dpscal(idg)    
            call get_flux_scalar(trcs_solml2_vr(idg,L,NY,NX),flux,trcs_solml_vr(idg,L,NY,NX),pscal(idg))
          endif
        ENDDO

        DO idg=idg_NH3,idg_NH3B
          if(dpscal(idg)>tiny_p)then
            flux=RGasAtmDisol2SoilM_col(idg,NY,NX)+trcsol_Irrig_flxM_vr(ids,L,NY,NX) & 
              +trcs_Transp2Micp_flxM_vr(ids,L,NY,NX) + trcs_Mac2MicPore_flxM_vr(ids,L,NY,NX)-RBGCSinkSoluteM_vr(ids,L,NY,NX)       
            flux=flux*dpscal(idg)    
            call get_flux_scalar(trcs_solml2_vr(idg,L,NY,NX),flux,trcs_solml_vr(idg,L,NY,NX),pscal(idg))
          endif
        ENDDO

        DO ids=ids_NH4B,ids_nuts_end
          if(dpscal(ids)>tiny_p)then
            flux = trcsol_Irrig_flxM_vr(ids,L,NY,NX)+trcs_Transp2Micp_flxM_vr(ids,L,NY,NX) + trcs_Mac2MicPore_flxM_vr(ids,L,NY,NX) &
              -RBGCSinkSoluteM_vr(ids,L,NY,NX)
            flux = flux*dpscal(ids)
            call get_flux_scalar(trcs_solml2_vr(ids,L,NY,NX),flux,trcs_solml_vr(ids,L,NY,NX),pscal(ids))
          endif
        ENDDO

        DO ids=ids_beg,ids_end
          if(dpscal(ids)>tiny_p)then
            flux = trcs_Transp2Macp_flxM_vr(ids,L,NY,NX) - trcs_Mac2MicPore_flxM_vr(ids,L,NY,NX)
            flux = flux*dpscal(ids)
            call get_flux_scalar(trcs_soHml2_vr(ids,L,NY,NX), flux, trcs_soHml_vr(ids,L,NY,NX), pscal(ids))
          endif
        ENDDO

        DO  K=1,jcplx
          DO idom=idom_beg,idom_end
            if(dpscal_dom(idom)>tiny_p)then
              flux=DOM_Transp2Micp_flxM_vr(idom,K,L,NY,NX) + DOM_Mac2MicPore_flxM_vr(idom,K,L,NY,NX)
              flux=flux*dpscal_dom(idom)
              call get_flux_scalar(DOM_MicP2_vr(idom,K,L,NY,NX),flux, DOM_MicP_vr(idom,K,L,NY,NX),pscal_dom(idom))

              flux=DOM_Transp2Macp_flxM_vr(idom,K,L,NY,NX) - DOM_Mac2MicPore_flxM_vr(idom,K,L,NY,NX)
              flux=flux*dpscal_dom(idom)
              call get_flux_scalar(DOM_MacP2_vr(idom,K,L,NY,NX), flux, DOM_MacP_vr(idom,K,L,NY,NX),pscal_dom(idom))
            endif
          ENDDO
        ENDDO
      ENDIF  
    ENDDO
  ENDDO
  if(present(pscal1))pscal1=pscal
  call PrintInfo('end '//subname)
  end subroutine UpdateSurfTracerM

!------------------------------------------------------------------------------------------

  subroutine UpdateSnowTracersM(I,J,M,NHW,NHE,NVN,NVS,dpscal,pscal1)
  implicit none
  integer, intent(in) :: I,J,M
  integer, intent(in) :: NHW,NHE,NVN,NVS
  real(r8),intent(in) :: dpscal(ids_beg:ids_end)
  real(r8),optional,intent(inout):: pscal1(ids_beg:ids_end)
  integer :: L,idg,idn,NY,NX
  real(r8) :: b,flux
  real(r8) :: pscal(ids_beg:ids_end)
!
!
  if(present(pscal1))then
    pscal=pscal1
  else
    pscal=1.001_r8
  endif  

  DO NY=NHW,NHE
    DO  NX=NVN,NVS  

      DO idg=idg_beg,idg_NH3
        if(dpscal(idg)>tiny_p)then
          flux=trcg_SnowDrift_flxM(idg,NY,NX)+trcg_Aqua_flxM_snvr(idg,1,NY,NX)
          flux=flux*dpscal(idg)
          call get_flux_scalar(trcg_solsml2_snvr(idg,1,NY,NX),flux,trcg_solsml_snvr(idg,1,NY,NX),pscal(idg))
        endif
      ENDDO

      DO idn=ids_nut_beg,ids_nuts_end
        if(dpscal(idn)>tiny_p)then
          flux=trcn_SnowDrift_flxM(idn,NY,NX)+trcn_Aqua_flxM_snvr(idn,1,NY,NX)
          flux=flux*dpscal(idn)
          call get_flux_scalar(trcn_solsml2_snvr(idn,1,NY,NX),flux,trcn_solsml_snvr(idn,1,NY,NX),pscal(idn))
        endif
      ENDDO

      DO  L=2,JS
        DO idg=idg_beg,idg_NH3
          if(trcg_solsml2_snvr(idg,L,NY,NX)>0._r8 .and. dpscal(idg)>tiny_p) then
            flux=trcg_Aqua_flxM_snvr(idg,L,NY,NX)*dpscal(idg)
            call get_flux_scalar(trcg_solsml2_snvr(idg,L,NY,NX),flux, trcg_solsml_snvr(idg,L,NY,NX),pscal(idg))
          endif  
        ENDDO

        DO idn=ids_nut_beg,ids_nuts_end
          if(trcn_solsml2_snvr(idn,L,NY,NX)>0._r8 .and. dpscal(idn)>tiny_p) then
            flux=trcn_Aqua_flxM_snvr(idn,L,NY,NX)*dpscal(idn)
            call get_flux_scalar(trcn_solsml2_snvr(idn,L,NY,NX),flux,trcn_solsml2_snvr(idn,L,NY,NX),pscal(idn))
          endif  
        ENDDO
      ENDDO
    ENDDO
  ENDDO
  if(present(pscal1))pscal1=pscal

  end subroutine UpdateSnowTracersM

end module TranspNoSaltSlowMod
