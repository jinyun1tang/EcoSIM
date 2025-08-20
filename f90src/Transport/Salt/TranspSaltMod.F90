module TranspSaltMod
!!
! Description:
!
  use data_kind_mod,    only: r8 => DAT_KIND_R8
  use minimathmod  
  use DebugToolMod
  use SOMDataType
  use GridConsts
  use FlagDataType
  use SoilPhysDataType
  use EcoSIMSolverPar
  use EcoSIMCtrlDataType
  use ClimForcDataType
  use FertilizerDataType
  use SoilWaterDataType
  use SurfLitterDataType
  use SnowDataType
  use SurfSoilDataType
  use ChemTranspDataType
  use LandSurfDataType
  use SoilBGCDataType
  use RootDataType
  use AqueChemDatatype
  use SoilPropertyDataType
  use IrrigationDataType
  use PlantDataRateType
  use GridDataType
  use IngridTranspMod
  use TranspSaltDataMod
  implicit none

  private
  CHARACTER(LEN=*), PARAMETER :: MOD_FILENAME=&
  __FILE__

  public :: TranspSalt
  contains

!----------------------------------------------------------------------
  SUBROUTINE TranspSalt(I,J,NHW,NHE,NVN,NVS)
!
!     Description:
!
!     THIS SUBROUTINE CALCULATES 3-DIMENSIONAL FLUXES OF ALL SOIL
!     SALT SOLUTES
!
  implicit none
  integer, intent(in) :: I, J
  integer, intent(in) :: NHW,NHE,NVN,NVS

  integer :: NX,NY,L,M,ids
  real(r8)  :: dpscal(idsalt_beg:idsaltb_end)
  real(r8)  :: pscal(idsalt_beg:idsaltb_end)
  real(r8) :: dpscal_max

!     execution begins here
!
!     TIME STEPS FOR SOLUTE FLUX CALCULATIONS
!
  call InitSaltTransportModel(I,NHW,NHE,NVN,NVS)
!
!     TIME STEP USED IN GAS AND SOLUTE FLUX CALCULATIONS
!
  D30: DO M=1,NPH
   !
    pscal=1._r8;dpscal=1._r8
    dpscal_max=1._r8

    do while(dpscal_max>1.e-2_r8)
      call ZeroFluxArraysM(NHW,NHE,NVN,NVS)

      call GetSaltTranspFlxM(I,M,NHW,NHE,NVN,NVS)
      !
      call UpdateSaltTranspM(M,NHW,NHE,NVN,NVS,dpscal,pscal)

      call AccumSaltFluxesM(M,NHW,NHE,NVN,NVS,dpscal,pscal)

      call copySaltStates(NHW,NHE,NVN,NVS)

      dpscal_max=0._r8
      DO ids=idsalt_beg,idsaltb_end
        dpscal(ids)=dpscal(ids)*(1._r8-pscal(ids))
        dpscal_max=AMAX1(dpscal_max,dpscal(ids))
      enddo

    enddo
  ENDDO D30
  RETURN

  END subroutine TranspSalt
!------------------------------------------------------------------------------------------

  subroutine AccumSaltFluxesM(M,NHW,NHE,NVN,NVS,dpscal,pscal)
  implicit none
  integer, intent(in) :: M,NHW,NHE,NVN,NVS
  real(r8), intent(in) :: dpscal(idsalt_beg:idsaltb_end)
  real(r8), intent(in):: pscal(idsalt_beg:idsaltb_end)

  real(r8) :: pscal_max
  real(r8) :: ppscal(idsalt_beg:idsaltb_end)
  integer  :: idsalt,NY,NX,ids

  pscal_max=0._r8
  DO ids=idsalt_beg,idsaltb_end
    ppscal(ids) = dpscal(ids)*pscal(ids)
    pscal_max=AMAX1(pscal(ids),pscal_max)
  ENDDO

  if(pscal_max<0.999_r8)then
    call UpdateSaltTranspM(M,NHW,NHE,NVN,NVS,ppscal)
  endif

  DO NX=NHW,NHE
    DO NY=NVN,NVS     
      DO idsalt=idsalt_beg,idsalt_end
        trcSalt_AquaADV_Snow2Litr_flx(idsalt,NY,NX)  = trcSalt_AquaADV_Snow2Litr_flx(idsalt,NY,NX)+ppscal(idsalt) &
          * trcSalt_AquaADV_Snow2Litr_flxM(idsalt,NY,NX)

        trcSalt_SnowDrift_flx_col(idsalt,NY,NX) = trcSalt_SnowDrift_flx_col(idsalt,NY,NX) + ppscal(idsalt) &
          * trcSalt_SnowDrift_flxM_col(idsalt,NY,NX)  
      ENDDO

      DO idsalt=idsalt_beg,idsaltb_end
        trcSalt_AquaADV_Snow2Soil_flx(idsalt,NY,NX)  = trcSalt_AquaADV_Snow2Soil_flx(idsalt,NY,NX)+ppscal(idsalt) &
          *trcSalt_AquaADV_Snow2Soil_flxM(idsalt,NY,NX)
      ENDDO

    ENDDO
  ENDDO  

  
  end subroutine AccumSaltFluxesM
!------------------------------------------------------------------------------------------

  subroutine ZeroAtmosSoluteFlux(NY,NX)
!
!     Description:
!
  implicit none
  integer, intent(in) :: NY,NX
  integer :: idsalt

  trcSalt_Precip2LitrM(:,NY,NX) = 0._r8

  trcSalt_Precip2MicpM(:,NY,NX) = 0._r8

  end subroutine ZeroAtmosSoluteFlux

!------------------------------------------------------------------------------------------

  subroutine AtmosSoluteFluxToTopsoil(I,NY,NX)
!
!     Description:
!
  implicit none
  integer, intent(in) :: I,NY,NX
  integer :: idsalt,ids
  character(len=*), parameter :: subname='AtmosSoluteFluxToTopsoil'


!     begin_execution
  call PrintInfo('beg '//subname)
  !
  DO idsalt=idsalt_beg,idsalt_end    
    trcSalt_Precip2LitrM(idsalt,NY,NX)   =(Rain2LitRSurf_col(NY,NX)*trcsalt_rain_mole_conc_col(idsalt,NY,NX) &
      +Irrig2LitRSurf_col(NY,NX)*trcsalt_irrig_mole_conc_col(idsalt,I,NY,NX))*dts_HeatWatTP  
  ENDDO


  DO idsalt=idsalt_beg,idsalt_KSO4
    trcSalt_Precip2MicpM(idsalt,NY,NX)=(Rain2SoilSurf_col(NY,NX)*trcsalt_rain_mole_conc_col(idsalt,NY,NX) &
      +Irrig2SoilSurf_col(NY,NX)*trcsalt_irrig_mole_conc_col(idsalt,I,NY,NX))*dts_HeatWatTP
  ENDDO
  
  DO idsalt=idsalt_H0PO4,idsalt_MgHPO4
    ids=idsalt-idsalt_H0PO4+idsalt_H0PO4B  
    !soil 
    trcSalt_Precip2MicpM(idsalt,NY,NX)=(Rain2SoilSurf_col(NY,NX)*trcsalt_rain_mole_conc_col(idsalt,NY,NX) &
      +Irrig2SoilSurf_col(NY,NX)*trcsalt_irrig_mole_conc_col(idsalt,I,NY,NX))*trcs_VLN_vr(ids_H1PO4,NU_col(NY,NX),NY,NX)*dts_HeatWatTP
    !band   
    trcSalt_Precip2MicpM(ids,NY,NX)=(Rain2SoilSurf_col(NY,NX)*trcsalt_rain_mole_conc_col(idsalt,NY,NX) &
      +Irrig2SoilSurf_col(NY,NX)*trcsalt_irrig_mole_conc_col(idsalt,I,NY,NX))*trcs_VLN_vr(ids_H1PO4B,NU_col(NY,NX),NY,NX)*dts_HeatWatTP    
  ENDDO
  call PrintInfo('end '//subname)
  end subroutine AtmosSoluteFluxToTopsoil

!------------------------------------------------------------------------------------------

  subroutine GetSubHourFlux(I,NY,NX)
!
!     Description:
!
  implicit none
  integer, intent(in) :: I,NY,NX
  integer :: L,idsalt,ids
  real(r8) :: FLWU


!     begin_execution
!
!     dts_HeatWatTP=1/no. of cycles h-1 for water, heat and solute flux calculations

  D10: DO L=NU_col(NY,NX),NL_col(NY,NX)
    DO idsalt=idsalt_beg,idsaltb_end
      trcSalt_RGeoChem_flxM_vr(idsalt,L,NY,NX)=-trcSalt_RGeoChem_flx_vr(idsalt,L,NY,NX)*dts_HeatWatTP
    ENDDO

    trcSalt_RGeoChem_flxM_vr(idsalt_Hp,L,NY,NX)=trcSalt_RGeoChem_flxM_vr(idsalt_Hp,L,NY,NX)-(RProd_Hp_vr(L,NY,NX))*dts_HeatWatTP
    !
    ! SOLUTE FLUXES FROM SUBSURFACE IRRIGATION
    ! placeholder: salt uptake by plant roots can be applied here. 
    !
    FLWU=TWaterPlantRoot2Soil_vr(L,NY,NX)*dts_HeatWatTP

    DO idsalt=idsalt_beg,idsalt_KSO4
      trcSalt_Irrig_vr(idsalt,L,NY,NX)=FWatIrrigate2MicP_vr(L,NY,NX)*trcsalt_irrig_mole_conc_col(idsalt,I,NY,NX)
    ENDDO

    DO idsalt=idsalt_H0PO4,idsalt_MgHPO4
      ids=idsalt-idsalt_H0PO4+idsalt_H0PO4B
      trcSalt_Irrig_vr(idsalt,L,NY,NX) = FWatIrrigate2MicP_vr(L,NY,NX)*trcsalt_irrig_mole_conc_col(idsalt,I,NY,NX)*trcs_VLN_vr(ids_H1PO4,L,NY,NX)
      trcSalt_Irrig_vr(ids,L,NY,NX)    = FWatIrrigate2MicP_vr(L,NY,NX)*trcsalt_irrig_mole_conc_col(idsalt,I,NY,NX)*trcs_VLN_vr(ids_H1PO4B,L,NY,NX)
    ENDDO
    !
    !     SUB-HOURLY SOLUTE FLUXES FROM SUBSURFACE IRRIGATION
    !
    DO idsalt=idsalt_beg,idsaltb_end
      trcSalt_Irrig_flxM_vr(idsalt,L,NY,NX)=trcSalt_Irrig_vr(idsalt,L,NY,NX)*dts_HeatWatTP
    ENDDO
    !
    !     SOLUTE DIFFUSIVITIES AT SUB-HOURLY TIME STEP
    !
    SoluteDifusivitytscaledM_vr(L,NY,NX)                 = SoluteDifusvty_vr(ids_H1PO4,L,NY,NX)*dts_HeatWatTP
    AquaIonDifusivty2_vr(idsalt_beg:idsalt_mend,L,NY,NX) = AquaIonDifusivty_vr(idsalt_beg:idsalt_mend,L,NY,NX)*dts_HeatWatTP

  ENDDO D10
  end subroutine GetSubHourFlux

!------------------------------------------------------------------------------------------
  subroutine copySaltStates(NHW,NHE,NVN,NVS)

  implicit none
  integer, intent(in) :: NHW,NHE,NVN,NVS
  integer :: NY,NX
  integer :: L,idsalt

  DO NX=NHW,NHE
    DO NY=NVN,NVS
      D20: DO L=1,JS
        DO idsalt=idsalt_beg,idsalt_end
          trcSalt_ml2_snvr(idsalt,L,NY,NX)=trcSalt_ml_snvr(idsalt,L,NY,NX)
        ENDDO
      ENDDO D20

      L=0
      DO idsalt=idsalt_beg,idsaltb_end
        trcSalt_solml2_vr(idsalt,L,NY,NX)=trcSalt_solml_vr(idsalt,L,NY,NX)
      ENDDO
      DO L=NU_col(NY,NX),NL_col(NY,NX)
        DO idsalt=idsalt_beg,idsaltb_end
          trcSalt_solml2_vr(idsalt,L,NY,NX)=trcSalt_solml_vr(idsalt,L,NY,NX)
          trcSalt_soHml2_vr(idsalt,L,NY,NX)=trcSalt_soHml_vr(idsalt,L,NY,NX)
        ENDDO
      ENDDO
    ENDDO
  ENDDO
  end subroutine copySaltStates
!------------------------------------------------------------------------------------------

  subroutine InitSaltTransportModel(I,NHW,NHE,NVN,NVS)
!
!     Description:
!
  implicit none
  integer, intent(in) :: I,NHW,NHE,NVN,NVS
  character(len=*), parameter :: subname='InitSaltTransportModel'
  integer :: NY,NX

  !     begin_execution
  call PrintInfo('beg '//subname)

  D9995: DO NX=NHW,NHE
    D9990: DO NY=NVN,NVS
      !
      call ZeroAtmosSoluteFlux(NY,NX)      

      IF(SnoFalPrec_col(NY,NX).GT.0.0.OR.(RainFalPrec_col(NY,NX).GT.0.0 &
        .AND.VLSnowHeatCapM_snvr(1,1,NY,NX).GT.VLHeatCapSnowMin_col(NY,NX)))THEN

        ! there is snowpack
        call SaltFall2Snowpack(I,NY,NX)
        !
        !     HOURLY SOLUTE FLUXES FROM ATMOSPHERE TO SOIL SURFACE
        !     IN RAINFALL AND IRRIGATION ACCORDING TO CONCENTRATIONS
        !     ENTERED IN WEATHER AND IRRIGATION FILES
        !
      ELSEIF((PrecAtm_col(NY,NX).GT.0.0_r8 .OR. IrrigSurface_col(NY,NX).GT.0.0_r8) & !there is significant precipiation
        .AND. VLSnowHeatCapM_snvr(1,1,NY,NX).LE.VLHeatCapSnowMin_col(NY,NX))THEN     !there is no significant snowpack
        
        call AtmosSoluteFluxToTopsoil(I,NY,NX)
        !
        !     NO SOLUTE FLUXES FROM ATMOSPHERE
        !
      ENDIF
      !
      !     GAS AND SOLUTE FLUXES AT SUB-HOURLY FLUX TIME STEP
      !     ENTERED IN SITE FILE
      !
      call GetSubHourFlux(I,NY,NX)

    ENDDO D9990
  ENDDO D9995

  call copySaltStates(NHW,NHE,NVN,NVS)  
  call PrintInfo('end '//subname)
  end subroutine InitSaltTransportModel


!------------------------------------------------------------------------------------------

  subroutine UpdateSaltTranspM(M,NHW,NHE,NVN,NVS,dpscal,pscal)
!
!     Description:
!
  implicit none
  integer , intent(in) :: M,NHW,NHE,NVN,NVS
  real(r8), intent(in) :: dpscal(idsalt_beg:idsaltb_end)
  real(r8), optional, intent(out):: pscal(idsalt_beg:idsaltb_end)

  character(len=*), parameter :: subname='UpdateSaltTranspM'
  integer :: L,idsalt,NY,NX
  real(r8) :: pscal1(idsalt_beg:idsaltb_end)
  real(r8) :: flux
  !
  !     begin_execution
  !
  call PrintInfo('beg '//subname)

  pscal1=1.001_r8
  D9695: DO NX=NHW,NHE
    D9690: DO NY=NVN,NVS
      !snow
      DO idsalt=idsalt_beg,idsalt_end
        flux=trcSalt_SnowDrift_flxM_col(idsalt,NY,NX)+trcSalt_Aqua_flxM_snvr(idsalt,1,NY,NX)
        flux=flux*dpscal(idsalt)
        call get_flux_scalar(trcSalt_ml2_snvr(idsalt,1,NY,NX),flux,trcSalt_ml_snvr(idsalt,1,NY,NX),pscal(idsalt))
      ENDDO

      D9670: DO L=2,JS
        DO idsalt=idsalt_beg,idsalt_end
          flux=trcSalt_Aqua_flxM_snvr(idsalt,L,NY,NX)
          flux=flux*dpscal(idsalt)
          call get_flux_scalar(trcSalt_ml2_snvr(idsalt,L,NY,NX),flux,trcSalt_ml_snvr(idsalt,L,NY,NX),pscal(idsalt))
        ENDDO
      ENDDO D9670
      !litter
      DO idsalt=idsalt_beg,idsalt_end        
        flux=trcSalt_FloXSurRof_flxM_col(idsalt,NY,NX)+trcSalt_MicpTranspFlxM_3D(idsalt,3,0,NY,NX)
        flux=flux*dpscal(idsalt)
        call get_flux_scalar(trcSalt_solml2_vr(idsalt,0,NY,NX),flux,trcSalt_solml_vr(idsalt,0,NY,NX),pscal(idsalt))
      ENDDO

      D9685: DO L=NU_col(NY,NX),NL_col(NY,NX)
        IF(VLSoilPoreMicP_vr(L,NY,NX).GT.ZEROS2(NY,NX))THEN
          DO idsalt=idsalt_beg,idsaltb_end            
            flux=trcSalt_Transp2Micp_flxM_vr(idsalt,L,NY,NX)+trcSalt_Mac2MicPore_flxM_vr(idsalt,L,NY,NX) &
              +trcSalt_Irrig_flxM_vr(idsalt,L,NY,NX)-trcSalt_RGeoChem_flxM_vr(idsalt,L,NY,NX)
            flux=flux*dpscal(idsalt)
            call get_flux_scalar(trcSalt_solml2_vr(idsalt,L,NY,NX),flux,trcSalt_solml2_vr(idsalt,L,NY,NX),pscal(idsalt))

            flux=trcSalt_Transp2Macp_flxM_vr(idsalt,L,NY,NX)-trcSalt_Mac2MicPore_flxM_vr(idsalt,L,NY,NX)
            flux=flux*dpscal(idsalt)
            call get_flux_scalar(trcSalt_soHml2_vr(idsalt,L,NY,NX),flux,trcSalt_soHml_vr(idsalt,L,NY,NX),pscal(idsalt))
            
          ENDDO
        ENDIF
      ENDDO D9685
    ENDDO D9690
  ENDDO D9695  
  if(present(pscal))pscal=pscal1
  call PrintInfo('end '//subname)
  end subroutine UpdateSaltTranspM

!------------------------------------------------------------------------------------------

  subroutine ZeroFluxArraysM(NHW,NHE,NVN,NVS)

  implicit none
  integer , intent(in) :: NHW,NHE,NVN,NVS

  integer :: NY,NX
  integer :: L, idsalt

  DO NX=NHW,NHE
    DO  NY=NVN,NVS
!
!     INITIALIZE SNOWPACK NET FLUX ACCUMULATORS
!
      D9855: DO L=1,JS
        DO idsalt=idsalt_beg,idsalt_end
          trcSalt_Aqua_flxM_snvr(idsalt,L,NY,NX)=0.0
        ENDDO
      ENDDO D9855
!
!     INITIALIZE SOIL SOLUTE NET FLUX ACCUMULATORS
!
      D9885: DO L=NU_col(NY,NX),NL_col(NY,NX)
        DO idsalt=idsalt_beg,idsalt_end      
          trcSalt_Transp2Micp_flxM_vr(idsalt,L,NY,NX)=0.0_r8
          trcSalt_Transp2Macp_flxM_vr(idsalt,L,NY,NX)=0.0_r8
        enddo
      ENDDO D9885      
    ENDDO
  ENDDO
  end subroutine ZeroFluxArraysM


end module TranspSaltMod
