module TranspSaltMod
!!
! Description:
!
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use minimathmod, only : isclose,AZMAX1
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

  real(r8) :: RCHQF,RCHGFU,RCHGFT
  real(r8) :: XN
  real(r8),allocatable ::  trcSalt_TQR(:,:,:)                         !

  real(r8),allocatable ::  trcSalt_TQ(:,:,:)                         !
  real(r8),allocatable ::  trcSalt_Flo2MicP_vr(:,:,:,:)                      !
  real(r8),allocatable ::  trcSalt_Flo2MacP_vr(:,:,:,:)                      !
  real(r8),allocatable ::  trcSalt_RFLZ(:,:,:,:)                      !
!----------------------------------------------------------------------



  public :: TranspSalt
  public :: initTranspSalt
  public :: destructTranspSalt
  contains
!----------------------------------------------------------------------
  subroutine DestructTranspSalt
  use abortutils, only : destroy
  implicit none

  call destroy(trcSalt_TQR)
  call destroy(trcSalt_Flo2MicP_vr)


  end subroutine DestructTranspSalt

!----------------------------------------------------------------------

  subroutine initTranspSalt()
  implicit none

  allocate(trcSalt_TQR(idsalt_beg:idsalt_end,JY,JX));       trcSalt_TQR=0._r8
  allocate(trcSalt_Flo2MicP_vr(idsalt_beg:idsaltb_end,JZ,JY,JX));   trcSalt_Flo2MicP_vr=0._r8
  allocate(trcSalt_Flo2MacP_vr(idsalt_beg:idsaltb_end,JZ,JY,JX));   trcSalt_Flo2MacP_vr=0._r8
  allocate(trcSalt_RFLZ(idsalt_beg:idsaltb_end,JZ,JY,JX));   trcSalt_RFLZ=0._r8
  end subroutine InitTranspSalt
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

  integer :: NX,NY,L,M
!     execution begins here

!
!     TIME STEPS FOR SOLUTE FLUX CALCULATIONS
!
  call SaltModelSoluteFlux(I,NHW,NHE,NVN,NVS)
!
!     TIME STEP USED IN GAS AND SOLUTE FLUX CALCULATIONS
!
  D30: DO M=1,NPH
!
    call InitFluxArrays(NHW,NHE,NVN,NVS)

    call SaltModelSoluteHydroFlux(I,M,NHW,NHE,NVN,NVS)
!
!     BOUNDARY SOLUTE AND GAS FLUXES
!
    call SaltModelInternalFlux(M,NHW,NHE,NVN,NVS)
!
!     UPDATE STATE VARIABLES FROM TOTAL FLUXES CALCULATED ABOVE
!
    D9695: DO NX=NHW,NHE
      D9690: DO NY=NVN,NVS
!
!     STATE VARIABLES FOR SOLUTES IN MICROPORES AND MACROPORES IN
!     SOIL SURFACE LAYER FROM OVERLAND FLOW
!
        call UpdateSoluteInSnow(NY,NX)
!
!     STATE VARIABLES FOR SOLUTES IN SURFACE RESIDUE FROM OVERLAND
!     FLOW AND SURFACE FLUX
        call UpdateSoluteInResidue(NY,NX)
!
!     STATE VARIABLES FOR GASES AND FOR SOLUTES IN MICROPORES AND
!     MACROPORES IN SOIL LAYERS FROM SUBSURFACE FLOW, EQUILIBRIUM
!     REACTIONS IN SOLUTE
        call UpdateSoluteInMicMacpores(NY,NX)
      ENDDO D9690
    ENDDO D9695
  ENDDO D30
  RETURN

  END subroutine TranspSalt
!------------------------------------------------------------------------------------------

  subroutine ZeroAtmosSoluteFlux(NY,NX)
!
!     Description:
!
  implicit none
  integer, intent(in) :: NY,NX
  integer :: idsalt
!     begin_execution
  DO idsalt=idsalt_beg,idsalt_end
    trcSaltFlo2SnowLay(idsalt,1,NY,NX)=0.0
    trcSalt3DFlo2Cell_3D(idsalt,3,0,NY,NX)=0.0
  ENDDO

  DO idsalt=idsalt_beg,idsaltb_end
    trcSalt3DFlo2Cell_3D(idsalt,3,NU(NY,NX),NY,NX)=0.0
  ENDDO
  end subroutine ZeroAtmosSoluteFlux
!------------------------------------------------------------------------------------------

  subroutine AtmosSoluteFluxToSnowpack(I,NY,NX)
!
!     Description:
!
  implicit none
  integer, intent(in) :: I,NY,NX
  integer :: idsalt

!     begin_execution

!
!     HOURLY SOLUTE FLUXES FROM ATMOSPHERE TO SOIL SURFACE
!     IN RAINFALL AND IRRIGATION ARE ZERO IF SNOWPACK IS PRESENT
!
!     X*FLS,X*FLB=hourly solute flux to micropores in non-band,band
!
  DO idsalt=idsalt_beg,idsalt_end
    trcSalt3DFlo2Cell_3D(idsalt,3,0,NY,NX)=0.0
    trcSaltFlo2SnowLay(idsalt,1,NY,NX)=Rain2SoilSurf_col(NY,NX)*trcsalt_rain_conc(idsalt,NY,NX) &
      +Irrig2SoilSurf_col(NY,NX)*trcsalt_irrig_conc(idsalt,I,NY,NX)    
  ENDDO

  DO idsalt=idsalt_beg,idsaltb_end
    trcSalt3DFlo2Cell_3D(idsalt,3,NU(NY,NX),NY,NX)=0.0
  ENDDO
  end subroutine AtmosSoluteFluxToSnowpack
!------------------------------------------------------------------------------------------

  subroutine AtmosSoluteFluxToTopsoil(I,NY,NX)
!
!     Description:
!
  implicit none
  integer, intent(in) :: I,NY,NX
  integer :: idsalt,ids
!     begin_execution
!
!     HOURLY SOLUTE FLUXES FROM ATMOSPHERE TO SNOWPACK
!     IN SNOWFALL AND IRRIGATION IS ZERO IF SNOWPACK IS ABSENT
!
!     PrecAtm_col,PRECI=snow+rain,irrigation
!     X*BLS=hourly solute flux to snowpack
!     X*FLS,X*FLB=hourly solute flux to surface litter,soil surface micropore non-band,band
!     Rain2LitRSurf_col,Irrig2LitRSurf=water flux to surface litter from rain,irrigation
!     FLQGQ,FLQGI=water flux to soil surface from rain,irrigation
!     C*R,C*Q=precipitation,irrigation solute concentrations
!     salt code: *HY*=H+,*OH*=OH-,*AL*=Al3+,*FE*=Fe3+,*CA*=Ca2+,*MG*=Mg2+
!          :*NA*=Na+,*KA*=K+,*SO4*=SO42-,*CL*=Cl-,*CO3*=CO32-,*HCO3*=HCO3-
!          :*CO2*=CO2,*ALO1*=AlOH2-,*ALOH2=AlOH2-,*ALOH3*=AlOH3
!          :*ALOH4*=AlOH4+,*ALS*=AlSO4+,*FEO1*=FeOH2-,*FEOH2=F3OH2-
!          :*FEOH3*=FeOH3,*FEOH4*=FeOH4+,*FES*=FeSO4+,*CAO*=CaOH
!          :*CAC*=CaCO3,*CAH*=CaHCO3-,*CAS*=CaSO4,*MGO*=MgOH,*MGC*=MgCO3
!          :*MHG*=MgHCO3-,*MGS*=MgSO4,*NAC*=NaCO3-,*NAS*=NaSO4-,*KAS*=KSO4-
!     phosphorus code: *H0P*=PO43-,*H3P*=H3PO4,*F1P*=FeHPO42-,*F2P*=F1H2PO4-
!          :*C0P*=CaPO4-,*C1P*=CaHPO4,*C2P*=CaH4P2O8+,*M1P*=MgHPO4,*COO*=COOH-
!          :*1=non-band,*B=band
!
  DO idsalt=idsalt_beg,idsalt_end
    trcSaltFlo2SnowLay(idsalt,1,NY,NX)=0.0_r8
    trcSalt3DFlo2Cell_3D(idsalt,3,0,NY,NX)=Rain2LitRSurf_col(NY,NX)*trcsalt_rain_conc(idsalt,NY,NX) &
      +Irrig2LitRSurf(NY,NX)*trcsalt_irrig_conc(idsalt,I,NY,NX)
  ENDDO

  DO idsalt=idsalt_beg,idsalt_KSO4
    trcSalt3DFlo2Cell_3D(idsalt,3,NU(NY,NX),NY,NX)=Rain2SoilSurf_col(NY,NX)*trcsalt_rain_conc(idsalt,NY,NX) &
      +Irrig2SoilSurf_col(NY,NX)*trcsalt_irrig_conc(idsalt,I,NY,NX)
  ENDDO
  
  DO idsalt=idsalt_H0PO4,idsalt_MgHPO4
    ids=idsalt-idsalt_H0PO4+idsalt_H0PO4B  
    trcSalt3DFlo2Cell_3D(idsalt,3,NU(NY,NX),NY,NX)=(Rain2SoilSurf_col(NY,NX)*trcsalt_rain_conc(idsalt,NY,NX) &
      +Irrig2SoilSurf_col(NY,NX)*trcsalt_irrig_conc(idsalt,I,NY,NX))*trcs_VLN_vr(ids_H1PO4,NU(NY,NX),NY,NX)
    trcSalt3DFlo2Cell_3D(ids,3,NU(NY,NX),NY,NX)=(Rain2SoilSurf_col(NY,NX)*trcsalt_rain_conc(idsalt,NY,NX) &
      +Irrig2SoilSurf_col(NY,NX)*trcsalt_irrig_conc(idsalt,I,NY,NX))*trcs_VLN_vr(ids_H1PO4B,NU(NY,NX),NY,NX)    
  ENDDO

  end subroutine AtmosSoluteFluxToTopsoil
!------------------------------------------------------------------------------------------

  subroutine InitSolutesInSnowpack(NY,NX)
!
!     Description:
!
  implicit none
  integer, intent(in) :: NY,NX

  integer :: L,NTA
!     begin_execution
!
!     Z*W=solute content in snowpacl
!     salt code: *HY*=H+,*OH*=OH-,*AL*=Al3+,*FE*=Fe3+,*CA*=Ca2+,*MG*=Mg2+
!          :*NA*=Na+,*KA*=K+,*SO4*=SO42-,*CL*=Cl-,*CO3*=CO32-,*HCO3*=HCO3-
!          :*CO2*=CO2,*ALO1*=AlOH2-,*ALOH2=AlOH2-,*ALOH3*=AlOH3
!          :*ALOH4*=AlOH4+,*ALS*=AlSO4+,*FEO1*=FeOH2-,*FEOH2=F3OH2-
!          :*FEOH3*=FeOH3,*FEOH4*=FeOH4+,*FES*=FeSO4+,*CAO*=CaOH
!          :*CAC*=CaCO3,*CAH*=CaHCO3-,*CAS*=CaSO4,*MGO*=MgOH,*MGC*=MgCO3
!          :*MHG*=MgHCO3-,*MGS*=MgSO4,*NAC*=NaCO3-,*NAS*=NaSO4-,*KAS*=KSO4-
!     phosphorus code: *H0P*=PO43-,*H3P*=H3PO4,*F1P*=FeHPO42-,*F2P*=F1H2PO4-
!          :*C0P*=CaPO4-,*C1P*=CaHPO4,*C2P*=CaH4P2O8+,*M1P*=MgHPO4,*COO*=COOH-C
  D20: DO L=1,JS
    DO NTA=idsalt_beg,idsalt_end
      trc_Saltml2_snvr(NTA,L,NY,NX)=trc_Saltml_snvr(NTA,L,NY,NX)
    ENDDO
  ENDDO D20
  end subroutine InitSolutesInSnowpack
!------------------------------------------------------------------------------------------

  subroutine GetSubHourFlux(NY,NX)
!
!     Description:
!
  implicit none
  integer, intent(in) :: NY,NX

  integer :: idsalt
!     begin_execution
!
!
!     dts_HeatWatTP=1/no. of cycles h-1 for water, heat and solute flux calculations
!     R*BLS,R*FL0,R*FL1,R*FL2=solute flux to snowpack,surface litter,soil surface non-band,band

!
  DO idsalt=idsalt_beg,idsalt_end
    trcSaltAdv2SowLay(idsalt,1,NY,NX) = trcSaltFlo2SnowLay(idsalt,1,NY,NX)*dts_HeatWatTP
    trcSalt_RFL0(idsalt,NY,NX)        = trcSalt3DFlo2Cell_3D(idsalt,3,0,NY,NX)*dts_HeatWatTP
  ENDDO

  DO idsalt=idsalt_beg,idsaltb_end
    trcSalt_RFL1(idsalt,NY,NX)=trcSalt3DFlo2Cell_3D(idsalt,3,NU(NY,NX),NY,NX)*dts_HeatWatTP
  ENDDO

  end subroutine GetSubHourFlux
!------------------------------------------------------------------------------------------

  subroutine GetSubHourlyFluxByLayer(I,NY,NX)
!
!     Description:
!
  implicit none
  integer, intent(in) :: I,NY,NX

  integer :: L,idsalt,ids
  real(r8) :: FLWU(JZ,JY,JX)


!     begin_execution
!
!     dts_HeatWatTP=1/no. of cycles h-1 for water, heat and solute flux calculations
!     RZ*2=solute flux at time step for flux calculations
!     TR*=solute transformations from solute.f
!
  D10: DO L=NU(NY,NX),NL(NY,NX)
    DO idsalt=idsalt_beg,idsaltb_end
      trcSalt_solml2R(idsalt,L,NY,NX)=-trcSalt_TR_vr(idsalt,L,NY,NX)*dts_HeatWatTP
    ENDDO

    trcSalt_solml2R(idsalt_Hp,L,NY,NX)=trcSalt_solml2R(idsalt_Hp,L,NY,NX)-(RProd_Hp_vr(L,NY,NX))*dts_HeatWatTP
!
!     SOLUTE FLUXES FROM SUBSURFACE IRRIGATION
!
!     FLU=subsurface water flux from watsub.f
!     R*FLU,R*FBU=subsurface solute flux in non-band,band
!     C*Q=irrigation solute concentrations
!     VLNH4,VLNO3,VLPO4=non-band NH4,NO3,PO4 volume fraction
!     VLNHB,VLNOB,VLPOB=band NH4,NO3,PO4 volume fraction
!
    FLWU(L,NY,NX)=TPlantRootH2OLoss_vr(L,NY,NX)*dts_HeatWatTP

    DO idsalt=idsalt_beg,idsalt_KSO4
      trcSalt_Irrig_vr(idsalt,L,NY,NX)=FWatIrrigate2MicP_vr(L,NY,NX)*trcsalt_irrig_conc(idsalt,I,NY,NX)
    ENDDO

    DO idsalt=idsalt_H0PO4,idsalt_MgHPO4
      ids=idsalt-idsalt_H0PO4+idsalt_H0PO4B
      trcSalt_Irrig_vr(idsalt,L,NY,NX)=FWatIrrigate2MicP_vr(L,NY,NX)*trcsalt_irrig_conc(idsalt,I,NY,NX)*trcs_VLN_vr(ids_H1PO4,L,NY,NX)
      trcSalt_Irrig_vr(ids,L,NY,NX)=FWatIrrigate2MicP_vr(L,NY,NX)*trcsalt_irrig_conc(idsalt,I,NY,NX)*trcs_VLN_vr(ids_H1PO4B,L,NY,NX)
    ENDDO

!
!     SUB-HOURLY SOLUTE FLUXES FROM SUBSURFACE IRRIGATION
!
!     R*FLZ,R*FBZ=subsurface solute flux in non-band,band
!     dts_HeatWatTP=1/no. of cycles h-1 for water, heat and solute flux calculations
!
    DO idsalt=idsalt_beg,idsaltb_end
      trcSalt_RFLZ(idsalt,L,NY,NX)=trcSalt_Irrig_vr(idsalt,L,NY,NX)*dts_HeatWatTP
    ENDDO
!
!     SOLUTE DIFFUSIVITIES AT SUB-HOURLY TIME STEP
!
!     dts_HeatWatTP=1/no. of cycles h-1 for water, heat and solute flux calculations
!     *SGL*=solute diffusivity from hour1.f
!     solute code:PO=PO4,AL=Al,FE=Fe,HY=H,CA=Ca,GM=Mg,AN=Na,AK=KOH=OH
!                :SO=SO4,CL=Cl,C3=CO3,HC=HCO3
    POSGL2(L,NY,NX)=SoluteDifusvty_vr(ids_H1PO4,L,NY,NX)*dts_HeatWatTP

    AquaIonDifusivty2_vr(idsalt_beg:idsalt_mend,L,NY,NX)=AquaIonDifusivty_vr(idsalt_beg:idsalt_mend,L,NY,NX)*dts_HeatWatTP
!
!     STATE VARIABLES FOR SOLUTES USED IN 'TranspSalt'
!     TO STORE SUB-HOURLY CHANGES DURING FLUX CALCULATIONS
!     INCLUDING TRANSFORMATIONS FROM REACTIONS IN SOLUTE.F
!
!
    DO idsalt=idsalt_beg,idsaltb_end
      trcSalt_solml2(idsalt,L,NY,NX)=trcSalt_solml_vr(idsalt,L,NY,NX)
      trcSalt_soHml2(idsalt,L,NY,NX)=trcSalt_soHml_vr(idsalt,L,NY,NX)
    ENDDO

  ENDDO D10
  end subroutine GetSubHourlyFluxByLayer
!------------------------------------------------------------------------------------------

  subroutine SaltModelSoluteFlux(I,NHW,NHE,NVN,NVS)
!
!     Description:
!
  implicit none
  integer, intent(in) :: I,NHW,NHE,NVN,NVS

  integer :: NY,NX
!     begin_execution

  D9995: DO NX=NHW,NHE
    D9990: DO NY=NVN,NVS
!
!     HOURLY SOLUTE FLUXES FROM ATMOSPHERE TO SNOWPACK
!     IN SNOWFALL AND IRRIGATION ACCORDING TO CONCENTRATIONS
!     ENTERED IN WEATHER AND IRRIGATION FILES
!
!
      IF(SnoFalPrec_col(NY,NX).GT.0.0.OR.(RainFalPrec_col(NY,NX).GT.0.0 &
        .AND.VLSnowHeatCapM_snvr(1,1,NY,NX).GT.VLHeatCapSnowMin_col(NY,NX)))THEN
  !     there is snowpack

        call AtmosSoluteFluxToSnowpack(I,NY,NX)
  !
  !     HOURLY SOLUTE FLUXES FROM ATMOSPHERE TO SOIL SURFACE
  !     IN RAINFALL AND IRRIGATION ACCORDING TO CONCENTRATIONS
  !     ENTERED IN WEATHER AND IRRIGATION FILES
  !
      ELSEIF((PrecAtm_col(NY,NX).GT.0.0_r8.OR.IrrigSurface_col(NY,NX).GT.0.0_r8) &
        .AND.VLSnowHeatCapM_snvr(1,1,NY,NX).LE.VLHeatCapSnowMin_col(NY,NX))THEN
  !     there is no significant snowpack, but precipitation
        call AtmosSoluteFluxToTopsoil(I,NY,NX)
  !
  !     NO SOLUTE FLUXES FROM ATMOSPHERE
  !
      ELSE
        !no precipitation
        call ZeroAtmosSoluteFlux(NY,NX)
      ENDIF
!
!     GAS AND SOLUTE FLUXES AT SUB-HOURLY FLUX TIME STEP
!     ENTERED IN SITE FILE
      call GetSubHourFlux(NY,NX)
!
!     INITIAL SOLUTES IN SNOWPACK
      call InitSolutesInSnowpack(NY,NX)
!
!     SOLUTE FLUXES FROM SOLUTE.F
      call GetSubHourlyFluxByLayer(I,NY,NX)

    ENDDO D9990
  ENDDO D9995
  end subroutine SaltModelSoluteFlux
!------------------------------------------------------------------------------------------

  subroutine InitFluxAccumulatorsInSnowpack(NY,NX)
!
!     Description:
!
  implicit none
  integer, intent(in) :: NY,NX
  integer :: L,idsalt

  D9855: DO L=1,JS
    DO idsalt=idsalt_beg,idsalt_end
      trcSalt_TBLS(idsalt,L,NY,NX)=0.0
    ENDDO
  ENDDO D9855
  end subroutine InitFluxAccumulatorsInSnowpack

!------------------------------------------------------------------------------------------

  subroutine ZeroBoundarySnowFlux(N,M5,M4)
!
!     Description:
!
  implicit none
  integer, intent(in) :: N,M5,M4
!     begin_execution
  trcSalt_RQ(idsalt_beg:idsalt_end,N,M5,M4)=0.0

  end subroutine ZeroBoundarySnowFlux
!------------------------------------------------------------------------------------------

  subroutine ZeroRechargeSoluteFlux(N,NN,M5,M4)
!
!     Description:
!
  implicit none
  integer, intent(in) :: N,NN,M5,M4
!     begin_execution

  trcSalt_FloXSurRof_flxM_2DH(idsalt_beg:idsalt_end,N,NN,M5,M4)=0.0

  end subroutine ZeroRechargeSoluteFlux
!------------------------------------------------------------------------------------------

  subroutine SoluteExportThruBoundary(N1,N2,M,N,NN,M5,M4)
!
!     Description:
!
  implicit none
  integer, intent(in) :: N1,N2,M,N,NN,M5,M4
  real(r8) :: FQRM
  integer :: idsalt
!     begin_execution

  FQRM=QflxSurfRunoffM_2DH(M,N,NN,M5,M4)/SurfRunoffWatFluxM_2DH(M,N2,N1)

  DO idsalt=idsalt_beg,idsalt_end
    trcSalt_FloXSurRof_flxM_2DH(idsalt,N,NN,M5,M4)=trcSalt_FloXSurRof_flxM(idsalt,N2,N1)*FQRM
  ENDDO
!
!     ACCUMULATE HOURLY FLUXES FOR USE IN REDIST.F
!
!     X*QRS=hourly solute in runoff
!     RQR*=solute in runoff
!
  DO idsalt=idsalt_beg,idsalt_end
    trc_salt_rof_bounds(idsalt,N,NN,M5,M4)=trc_salt_rof_bounds(idsalt,N,NN,M5,M4)+trcSalt_FloXSurRof_flxM_2DH(idsalt,N,NN,M5,M4)
  ENDDO
  end subroutine SoluteExportThruBoundary
!------------------------------------------------------------------------------------------

  subroutine ZeroSoluteInfluxThruBoundary(N,NN,M5,M4)
!
!     Description:
!
  implicit none
  integer, intent(in) :: N,NN,M5,M4
!     begin_execution

  trcSalt_FloXSurRof_flxM_2DH(idsalt_beg:idsalt_end,N,NN,M5,M4)=0.0

  end subroutine ZeroSoluteInfluxThruBoundary
!------------------------------------------------------------------------------------------

  subroutine SoluteGainSubsurfMicropore(M3,M2,M1,M,N,M6,M5,M4)
!
!     Description:
!
  implicit none
  integer, intent(in) :: M,N
  integer, intent(in) :: M3,M2,M1
  integer, intent(in) :: M6,M5,M4
  integer :: idsalt, ids
  real(r8) :: VFLOW
!     begin_execution

  VFLOW=WaterFlow2MicPM_3D(M,N,M6,M5,M4)
  DO idsalt=idsalt_beg,idsalt_KSO4
    trcSalt3DFlo2CellM(idsalt,N,M6,M5,M4)=VFLOW*trcsalt_subirrig_conc(idsalt,M3,M2,M1)
  ENDDO
  DO idsalt=idsalt_H0PO4,idsalt_MgHPO4
    ids=idsalt-idsalt_H0PO4+idsalt_H0PO4B
    trcSalt3DFlo2CellM(idsalt,N,M6,M5,M4)=VFLOW*trcsalt_subirrig_conc(idsalt,M3,M2,M1)*trcs_VLN_vr(ids_H1PO4,M3,M2,M1)
    trcSalt3DFlo2CellM(ids,N,M6,M5,M4)=VFLOW*trcsalt_subirrig_conc(idsalt,M3,M2,M1)*trcs_VLN_vr(ids_H1PO4B,M3,M2,M1)
  ENDDO
  end subroutine SoluteGainSubsurfMicropore
!------------------------------------------------------------------------------------------

  subroutine SoluteLossSubsurfMicropore(M,N,M1,M2,M3,M4,M5,M6)
!
!     Description:
!
  implicit none
  integer, intent(in) :: M,N,M1,M2,M3,M4,M5,M6
  real(r8) :: VFLW
  integer :: idsalt,ids
!     begin_execution

  IF(VLWatMicPM_vr(M,M3,M2,M1).GT.ZEROS2(M2,M1))THEN
    VFLW=AMAX1(-VFLWX,AMIN1(VFLWX,WaterFlow2MicPM_3D(M,N,M6,M5,M4) &
      /VLWatMicPM_vr(M,M3,M2,M1)))
  ELSE
    VFLW=0.0_r8
  ENDIF
  DO idsalt=idsalt_beg,idsalt_KSO4
    trcSalt3DFlo2CellM(idsalt,N,M6,M5,M4)=VFLW*AZMAX1(trcSalt_solml2(idsalt,M3,M2,M1))
  ENDDO

  DO idsalt=idsalt_H0PO4,idsalt_MgHPO4
    ids=idsalt-idsalt_H0PO4+idsalt_H0PO4B
    trcSalt3DFlo2CellM(idsalt,N,M6,M5,M4)=VFLW*AZMAX1(trcSalt_solml2(idsalt,M3,M2,M1))*trcs_VLN_vr(ids_H1PO4,M3,M2,M1)
    trcSalt3DFlo2CellM(ids,N,M6,M5,M4)=VFLW*AZMAX1(trcSalt_solml2(ids,M3,M2,M1))*trcs_VLN_vr(ids_H1PO4B,M3,M2,M1)
  ENDDO
  end subroutine SoluteLossSubsurfMicropore
!------------------------------------------------------------------------------------------

  subroutine SoluteLossSubsurfMacropore(M,N,M1,M2,M3,M4,M5,M6)
!
!     Description:
!
  implicit none
  integer, intent(in) :: M,N,M1,M2,M3,M4,M5,M6
  real(r8) :: VFLW
  integer :: idsalt,ids
!     begin_execution

  IF(VLWatMacPM_vr(M,M3,M2,M1).GT.ZEROS2(M2,M1))THEN
    VFLW=AMAX1(-VFLWX,AMIN1(VFLWX,WaterFlow2MacPM_3D(M,N,M6,M5,M4)/VLWatMacPM_vr(M,M3,M2,M1)))

    DO idsalt=idsalt_beg,idsalt_KSO4
      trcSalt_RFHS(idsalt,N,M6,M5,M4)=VFLW*AZMAX1(trcSalt_soHml2(idsalt,M3,M2,M1))
    ENDDO

    DO idsalt=idsalt_H0PO4,idsalt_MgHPO4
      ids=idsalt-idsalt_H0PO4+idsalt_H0PO4B
      trcSalt_RFHS(idsalt,N,M6,M5,M4)=VFLW*AZMAX1(trcSalt_soHml2(idsalt,M3,M2,M1)) &
        *trcs_VLN_vr(ids_H1PO4,M3,M2,M1)
      trcSalt_RFHS(ids,N,M6,M5,M4)=VFLW*AZMAX1(trcSalt_soHml2(ids,M3,M2,M1)) &
        *trcs_VLN_vr(ids_H1PO4B,M3,M2,M1)
    ENDDO

  ELSE
    VFLW=0.0_r8
    trcSalt_RFHS(idsalt_beg:idsaltb_end,N,M6,M5,M4)=0._r8
  ENDIF

  end subroutine SoluteLossSubsurfMacropore

!------------------------------------------------------------------------------------------

  subroutine AccumFluxMacMicPores(N,M6,M5,M4)
!
!     Description:
!
  implicit none
  integer, intent(in) :: N,M6,M5,M4
  integer :: idsalt
!     begin_execution
!
!     X*FLS,X*FLW,X*FLB=hourly solute flux in non-band,band micropores
!     X*FHS,X*FHW,X*FHB=hourly solute flux in non-band,band macropores
!     R*FLS,R*FLW,R*FLB=solute flux in non-band,band micropores
!     R*FHS,R*FHW,R*FHB=solute flux in non-band,band macropores
!
  DO idsalt=idsalt_beg,idsaltb_end
    trcSalt3DFlo2Cell_3D(idsalt,N,M6,M5,M4)=trcSalt3DFlo2Cell_3D(idsalt,N,M6,M5,M4)+trcSalt3DFlo2CellM(idsalt,N,M6,M5,M4)
    trcSalt_XFHS_3D(idsalt,N,M6,M5,M4)=trcSalt_XFHS_3D(idsalt,N,M6,M5,M4)+trcSalt_RFHS(idsalt,N,M6,M5,M4)
  ENDDO
  end subroutine AccumFluxMacMicPores
!------------------------------------------------------------------------------------------

  subroutine NetOverloadFluxInWater(M,N,N1,N2,N4,N5,N4B,N5B)
!
!     Description:
!
  integer, intent(in) :: M,N,N1,N2,N4,N5,N4B,N5B
  integer :: NN,idsalt
!     begin_execution
!
!     TQR*=net overland solute flux
!     RQR*=overland solute flux
!     salt code: *HY*=H+,*OH*=OH-,*AL*=Al3+,*FE*=Fe3+,*CA*=Ca2+,*MG*=Mg2+
!          :*NA*=Na+,*KA*=K+,*SO4*=SO42-,*CL*=Cl-,*CO3*=CO32-,*HCO3*=HCO3-
!          :*CO2*=CO2,*ALO1*=AlOH2-,*ALOH2=AlOH2-,*ALOH3*=AlOH3
!          :*ALOH4*=AlOH4+,*ALS*=AlSO4+,*FEO1*=FeOH2-,*FEOH2=F3OH2-
!          :*FEOH3*=FeOH3,*FEOH4*=FeOH4+,*FES*=FeSO4+,*CAO*=CaOH
!          :*CAC*=CaCO3,*CAH*=CaHCO3-,*CAS*=CaSO4,*MGO*=MgOH,*MGC*=MgCO3
!          :*MHG*=MgHCO3-,*MGS*=MgSO4,*NAC*=NaCO3-,*NAS*=NaSO4-,*KAS*=KSO4-
!     phosphorus code: *H0P*=PO43-,*H3P*=H3PO4,*F1P*=FeHPO42-,*F2P*=F1H2PO4-
!          :*C0P*=CaPO4-,*C1P*=CaHPO4,*C2P*=CaH4P2O8+,*M1P*=MgHPO4,*COO*=COOH-
!          :*1=non-band,*B=band
  D1202: DO NN=1,2
    DO idsalt=idsalt_beg,idsalt_end
      trcSalt_TQR(idsalt,N2,N1)=trcSalt_TQR(idsalt,N2,N1)+trcSalt_FloXSurRof_flxM_2DH(idsalt,N,NN,N2,N1)
    ENDDO

    IF(IFLBM(M,N,NN,N5,N4).EQ.0)THEN
      DO idsalt=idsalt_beg,idsalt_end
        trcSalt_TQR(idsalt,N2,N1)=trcSalt_TQR(idsalt,N2,N1)-trcSalt_FloXSurRof_flxM_2DH(idsalt,N,NN,N5,N4)
      ENDDO

    ENDIF
    IF(N4B.GT.0.AND.N5B.GT.0.AND.NN.EQ.iOutflow)THEN
      DO idsalt=idsalt_beg,idsalt_end
        trcSalt_TQR(idsalt,N2,N1)=trcSalt_TQR(idsalt,N2,N1)-trcSalt_FloXSurRof_flxM_2DH(idsalt,N,NN,N5B,N4B)
      ENDDO
    ENDIF
  ENDDO D1202
  end subroutine NetOverloadFluxInWater
!------------------------------------------------------------------------------------------

  subroutine NetOverloadFLuxInSnow(N,N1,N2,N4,N5)
!
!     Description:
!
  implicit none
  integer, intent(in) :: N,N1,N2,N4,N5

  integer :: idsalt
!     begin_execution
!     TQS*=net solute flux in snow transfer
!     RQS*=solute flux in snow transfer
!
  DO idsalt=idsalt_beg,idsalt_end
    trcSalt_TQ(idsalt,N2,N1)=trcSalt_TQ(idsalt,N2,N1)+trcSalt_RQ(idsalt,N,N2,N1)-trcSalt_RQ(idsalt,N,N5,N4)
  ENDDO
  end subroutine NetOverloadFLuxInSnow
!------------------------------------------------------------------------------------------

  subroutine NetFluxInSnowpack(M,NY,NX,N1,N2)
!
!     Description:
!
  integer, intent(in) :: M,NY,NX,N1,N2
  integer :: LS,LS2,idsalt
!     begin_execution

  D1205: DO LS=1,JS
    IF(VLSnowHeatCapM_snvr(M,LS,NY,NX).GT.VLHeatCapSnowMin_col(NY,NX))THEN
      LS2=MIN(JS,LS+1)
!
!     IF LOWER LAYER IS IN THE SNOWPACK
  !
      IF(LS.LT.JS.AND.VLSnowHeatCapM_snvr(M,LS2,N2,N1).GT.VLHeatCapSnowMin_col(N2,N1))THEN
        DO idsalt=idsalt_beg,idsalt_end
          trcSalt_TBLS(idsalt,LS,NY,NX)=trcSalt_TBLS(idsalt,LS,NY,NX) &
            +trcSaltAdv2SowLay(idsalt,LS,NY,NX)-trcSaltAdv2SowLay(idsalt,LS2,NY,NX)
        ENDDO
      ELSE
  !
  !     IF LOWER LAYER IS THE LITTER AND SOIL SURFACE
  !
        DO idsalt=idsalt_beg,idsalt_end
          trcSalt_TBLS(idsalt,LS,NY,NX)=trcSalt_TBLS(idsalt,LS,NY,NX)+trcSaltAdv2SowLay(idsalt,LS,NY,NX) &
            -trcSalt3DFlo2CellM(idsalt,3,0,N2,N1)-trcSalt3DFlo2CellM(idsalt,3,NUM(N2,N1),N2,N1)
        ENDDO

        DO idsalt=0,idsalt_nuts
          trcSalt_TBLS(idsalt_H0PO4+idsalt,LS,NY,NX)=trcSalt_TBLS(idsalt_H0PO4+idsalt,LS,NY,NX) &
            -trcSalt3DFlo2CellM(idsalt_H0PO4B+idsalt,3,NUM(N2,N1),N2,N1)
        ENDDO
      ENDIF
    ENDIF
  ENDDO D1205
  end subroutine NetFluxInSnowpack
!------------------------------------------------------------------------------------------

  subroutine TotFluxInMacMicPores(N,N1,N2,N3,N4,N5,NY,NX,N6)
!
!     Description:
!
  implicit none
  integer, intent(in) :: N,N1,N2,N3,N4,N5,NY,NX
  integer, intent(inout) :: N6
  integer :: LL,idsalt
!     begin_execution
!
!     T*FLS=net convective + diffusive solute flux through micropores
!     R*FLS=convective + diffusive solute flux through micropores
!     R*FLW,R*FLB=convective + diffusive solute flux through micropores in non-band,band
!     T*FHS=net convective + diffusive solute flux through macropores
!     R*FHS=convective + diffusive solute flux through macropores
!     R*FHW,R*FHB=convective + diffusive solute flux through macropores in non-band,band


  IF(FlowDirIndicator(N2,N1).NE.3.OR.N.EQ.3)THEN
    D1200: DO LL=N6,NL(NY,NX)
      IF(VLSoilPoreMicP_vr(LL,N2,N1).GT.ZEROS2(N2,N1))THEN
        N6=LL
        exit
      ENDIF
    ENDDO D1200

    IF(VLSoilPoreMicP_vr(N3,N2,N1).GT.ZEROS2(N2,N1))THEN
      DO idsalt=idsalt_beg,idsaltb_end
        trcSalt_Flo2MicP_vr(idsalt,N3,N2,N1)=trcSalt_Flo2MicP_vr(idsalt,N3,N2,N1) &
          +trcSalt3DFlo2CellM(idsalt,N,N3,N2,N1)-trcSalt3DFlo2CellM(idsalt,N,N6,N5,N4)
        trcSalt_Flo2MacP_vr(idsalt,N3,N2,N1)=trcSalt_Flo2MacP_vr(idsalt,N3,N2,N1) &
          +trcSalt_RFHS(idsalt,N,N3,N2,N1)-trcSalt_RFHS(idsalt,N,N6,N5,N4)
      ENDDO

    ELSE
      trcSalt_Flo2MicP_vr(idsalt_beg:idsaltb_end,N3,N2,N1)=0.0
      trcSalt_Flo2MacP_vr(idsalt_beg:idsaltb_end,N3,N2,N1)=0.0
    ENDIF
  ENDIF
  end subroutine TotFluxInMacMicPores
!------------------------------------------------------------------------------------------

  subroutine SaltModelInternalFlux(M,NHW,NHE,NVN,NVS)
!
!     Description:
!
  implicit none
  integer, intent(in) :: M,NHW,NHE,NVN,NVS
  integer :: NY,NX,L,N1,N2,N3,N4,N5,N6,N4B,N5B,N,NN
  integer :: M1,M2,M3,M4,M5,M6
!     begin_execution
!     N3,N2,N1=L,NY,NX of source grid cell
!     M6,M5,M4=L,NY,NX of destination grid cell
!

  D9595: DO NX=NHW,NHE
    D9590: DO NY=NVN,NVS
      D9585: DO L=NU(NY,NX),NL(NY,NX)
        !source grid
        N1=NX;N2=NY;N3=L
!
!     LOCATE ALL EXTERNAL BOUNDARIES AND SET BOUNDARY CONDITIONS
!     ENTERED IN 'READS'
!
        !locate dest grid
        D9580: DO N=FlowDirIndicator(NY,NX),3
          D9575: DO NN=1,2
            IF(N.EQ.iEastWestDirection)THEN
              N4=NX+1
              N5=NY
              N4B=NX-1
              N5B=NY
              N6=L
              IF(NN.EQ.iOutflow)THEN
                IF(NX.EQ.NHE)THEN
                  M1 = NX
                  M2 = NY
                  M3 = L
                  M4 = NX+1
                  M5 = NY
                  M6 = L
                  XN=-1.0_r8  !out of the boundary
                  RCHQF=RechargEastSurf(M2,M1)
                  RCHGFU=RechargEastSubSurf(M2,M1)
                  RCHGFT=RechargRateEastWTBL(M2,M1)
                ELSE
                  cycle
                ENDIF
              ELSEIF(NN.EQ.iInflow)THEN
                IF(NX.EQ.NHW)THEN
                  M1 = NX
                  M2 = NY
                  M3 = L
                  M4 = NX
                  M5 = NY
                  M6 = L
                  XN = 1.0  !in from the boundary
                  RCHQF=RechargWestSurf(M5,M4)
                  RCHGFU=RechargWestSubSurf(M5,M4)
                  RCHGFT=RechargRateWestWTBL(M5,M4)
                ELSE
                  CYCLE
                ENDIF
              ENDIF
            ELSEIF(N.EQ.iNorthSouthDirection)THEN
              N4=NX
              N5=NY+1
              N4B=NX
              N5B=NY-1
              N6=L
              IF(NN.EQ.iOutflow)THEN
                IF(NY.EQ.NVS)THEN
                  M1=NX
                  M2=NY
                  M3=L
                  M4=NX
                  M5=NY+1
                  M6=L
                  XN=-1.0_r8  !out of the boundary
                  RCHQF=RechargSouthSurf(M2,M1)
                  RCHGFU=RechargSouthSubSurf(M2,M1)
                  RCHGFT=RechargRateSouthWTBL(M2,M1)
                ELSE
                  CYCLE
                ENDIF
              ELSEIF(NN.EQ.iInflow)THEN
                IF(NY.EQ.NVN)THEN
                  M1=NX
                  M2=NY
                  M3=L
                  M4=NX
                  M5=NY
                  M6=L
                  XN=1.0_r8  !in from the boundary
                  RCHQF=RechargNorthSurf(M5,M4)
                  RCHGFU=RechargNorthSubSurf(M5,M4)
                  RCHGFT=RechargRateNorthWTBL(M5,M4)
                ELSE
                  CYCLE
                ENDIF
              ENDIF
            ELSEIF(N.EQ.iVerticalDirection)THEN
              N1 = NX
              N2 = NY
              N3 = L
              N4 = NX
              N5 = NY
              N6 = L+1
              IF(NN.EQ.iOutflow)THEN
                IF(L.EQ.NL(NY,NX))THEN
                  M1 = NX
                  M2 = NY
                  M3 = L
                  M4 = NX
                  M5 = NY
                  M6 = L+1
                  XN=-1.0_r8  !out of the boundary
                ELSE
                  CYCLE
                ENDIF
              ELSEIF(NN.EQ.iInflow)THEN
                CYCLE
              ENDIF
            ENDIF
!
!     SURFACE SOLUTE TRANSPORT FROM BOUNDARY SURFACE
!     RUNOFF IN WATSUB AND CONCENTRATIONS IN THE SURFACE SOIL LAYER
!
!
!     SURFACE SOLUTE TRANSPORT FROM BOUNDARY SURFACE
!     RUNOFF IN 'WATSUB' AND CONCENTRATIONS IN THE SURFACE SOIL LAYER
!
!     QRM =runoff from watsub.f
!     RQR*=solute in runoff
!
            IF(L.EQ.NUM(M2,M1).AND.N.NE.3)THEN
              IF(.not.XGridRunoffFlag(NN,N,N2,N1).OR.isclose(RCHQF,0.0_r8) &
                .OR.SurfRunoffWatFluxM_2DH(M,N2,N1).LE.ZEROS(N2,N1))THEN
                call ZeroRechargeSoluteFlux(N,NN,M5,M4)
              ELSE
!
!     SOLUTE LOSS FROM RUNOFF DEPENDING ON ASPECT
!     AND BOUNDARY CONDITIONS SET IN SITE FILE
!
                IF((NN.EQ.iOutflow.AND.QflxSurfRunoffM_2DH(M,N,NN,M5,M4).GT.ZEROS(N2,N1)) &
                  .OR.(NN.EQ.iInflow.AND.QflxSurfRunoffM_2DH(M,N,NN,M5,M4).LT.ZEROS(N2,N1)))THEN

                  call SoluteExportThruBoundary(N1,N2,M,N,NN,M5,M4)
!
!     SOLUTE GAIN FROM RUNON DEPENDING ON ASPECT
!     AND BOUNDARY CONDITIONS SET IN SITE FILE
!
                ELSEIF((NN.EQ.iInflow.AND.QflxSurfRunoffM_2DH(M,N,NN,M5,M4).GT.ZEROS(N2,N1)) &
                  .OR.(NN.EQ.iOutflow.AND.QflxSurfRunoffM_2DH(M,N,NN,M5,M4).LT.ZEROS(N2,N1)))THEN
                  call ZeroSoluteInfluxThruBoundary(N,NN,M5,M4)
                ELSE
                  call ZeroSoluteInfluxThruBoundary(N,NN,M5,M4)
                ENDIF
              ENDIF
!
!     BOUNDARY SNOW FLUX
!
              IF(NN.EQ.iOutflow)THEN
                call ZeroBoundarySnowFlux(N,M5,M4)
              ENDIF
            ENDIF
!
!     SOLUTE LOSS WITH SUBSURFACE MICROPORE WATER LOSS
!
!     WaterFlow2MicPM_3D=water flux through soil micropore from watsub.f
!     VLWatMicPM=micropore water-filled porosity from watsub.f
!     R*FLS=convective solute flux through micropores
!     R*FLW,R*FLB=convective solute flux through micropores in non-band,band
!
            IF(VLSoilPoreMicP_vr(N3,N2,N1).GT.ZEROS(NY,NX))THEN

              IF(FlowDirIndicator(M2,M1).NE.3.OR.N.EQ.3)THEN
                IF(NN.EQ.iOutflow.AND.WaterFlow2MicPM_3D(M,N,M6,M5,M4).GT.0.0 &
                  .OR.NN.EQ.iInflow.AND.WaterFlow2MicPM_3D(M,N,M6,M5,M4).LT.0.0)THEN

                  call SoluteLossSubsurfMicropore(M,N,M1,M2,M3,M4,M5,M6)

                ELSE
!     SOLUTE GAIN WITH SUBSURFACE MICROPORE WATER GAIN
!
                  call SoluteGainSubsurfMicropore(M3,M2,M1,M,N,M6,M5,M4)
                ENDIF
!
!     SOLUTE LOSS WITH SUBSURFACE MACROPORE WATER LOSS
!
!     WaterFlow2MacPM_3D=water flux through soil macropore from watsub.f
!     VLWatMacPM=macropore water-filled porosity from watsub.f
!     RFH*S=solute diffusive flux through macropore

                IF(NN.EQ.iOutflow.AND.WaterFlow2MacPM_3D(M,N,M6,M5,M4).GT.0.0 &
                  .OR.NN.EQ.iInflow.AND.WaterFlow2MacPM_3D(M,N,M6,M5,M4).LT.0.0)THEN

                  call SoluteLossSubsurfMacropore(M,N,M1,M2,M3,M4,M5,M6)

                ELSE
!
!     NO SOLUTE GAIN IN SUBSURFACE MACROPORES
!
                  trcSalt_RFHS(idsalt_beg:idsaltb_end,N,M6,M5,M4)=0.0_r8
                ENDIF
!
!     ACCUMULATE HOURLY FLUXES FOR USE IN REDIST.F
                call AccumFluxMacMicPores(N,M6,M5,M4)
              ENDIF
            ENDIF
          ENDDO D9575
!
!     TOTAL SOLUTE FLUXES IN EACH GRID CELL
!
          IF(L.EQ.NUM(N2,N1))THEN
            IF(N.NE.3)THEN
!
!     NET OVERLAND SOLUTE FLUX IN WATER
!
              call NetOverloadFluxInWater(M,N,N1,N2,N4,N5,N4B,N5B)
!
!     NET OVERLAND SOLUTE FLUX IN SNOW
              call NetOverloadFLuxInSnow(N,N1,N2,N4,N5)
!
            ELSEIF(N.EQ.3)THEN
!     NET SOLUTE FLUX IN SNOWPACK
!
!     VLSnowHeatCapM,VLHeatCapSnowMin_col=current,minimum volumetric heat capacity of snowpack
!     T*BLS=net solute flux in snowpack
!     R*BLS=solute flux in snowpack
              call NetFluxInSnowpack(M,NY,NX,N1,N2)
            ENDIF
          ENDIF
!
!     TOTAL SOLUTE FLUX IN MICROPORES AND MACROPORES
!
          call TotFluxInMacMicPores(N,N1,N2,N3,N4,N5,NY,NX,N6)
        ENDDO D9580
      ENDDO D9585
    ENDDO D9590
  ENDDO D9595
  end subroutine SaltModelInternalFlux
!------------------------------------------------------------------------------------------

  subroutine UpdateSoluteInSnow(NY,NX)
!
!     Description:
!
  integer, intent(in) :: NY,NX
  integer :: L,NTA
!     begin_execution
!     *W2=solute content of snowpack
!     TQS*=net overland solute flux in snow
!     T*BLS=net solute flux in snowpack
!     salt code: *HY*=H+,*OH*=OH-,*AL*=Al3+,*FE*=Fe3+,*CA*=Ca2+,*MG*=Mg2+
!          :*NA*=Na+,*KA*=K+,*SO4*=SO42-,*CL*=Cl-,*CO3*=CO32-,*HCO3*=HCO3-
!          :*CO2*=CO2,*ALO1*=AlOH2-,*ALOH2=AlOH2-,*ALOH3*=AlOH3
!          :*ALOH4*=AlOH4+,*ALS*=AlSO4+,*FEO1*=FeOH2-,*FEOH2=F3OH2-
!          :*FEOH3*=FeOH3,*FEOH4*=FeOH4+,*FES*=FeSO4+,*CAO*=CaOH
!          :*CAC*=CaCO3,*CAH*=CaHCO3-,*CAS*=CaSO4,*MGO*=MgOH,*MGC*=MgCO3
!          :*MHG*=MgHCO3-,*MGS*=MgSO4,*NAC*=NaCO3-,*NAS*=NaSO4-,*KAS*=KSO4-
!     phosphorus code: *H0P*=PO43-,*H3P*=H3PO4,*F1P*=FeHPO42-,*F2P*=F1H2PO4-
!          :*C0P*=CaPO4-,*C1P*=CaHPO4,*C2P*=CaH4P2O8+,*M1P*=MgHPO4,*COO*=COOH-
!          :*1=non-band,*B=band
!
  DO NTA=idsalt_beg,idsalt_end
    trc_Saltml2_snvr(NTA,1,NY,NX)=trc_Saltml2_snvr(NTA,1,NY,NX)+trcSalt_TQ(NTA,NY,NX)
  ENDDO
  D9670: DO L=1,JS
    DO NTA=idsalt_beg,idsalt_end
      trc_Saltml2_snvr(NTA,L,NY,NX)=trc_Saltml2_snvr(NTA,L,NY,NX)+trcSalt_TBLS(NTA,L,NY,NX)
    ENDDO
  ENDDO D9670
  end subroutine UpdateSoluteInSnow
!------------------------------------------------------------------------------------------

  subroutine UpdateSoluteInResidue(NY,NX)
!
!     Description:
!
  integer, intent(in) :: NY,NX

  integer :: idsalt
!     begin_execution
!
!     *S2=litter solute content
!     R*DFR=gas exchange between atmosphere and surface litter water
!     R*DFS=gas exchange between atmosphere and soil surface water
!     R*FLS=convective + diffusive solute flux into litter,soil surface
!     R*FLW,R*FLB=convective + diffusive solute flux into litter from non-band,band
!     TQR*=net overland solute flux

!
  DO idsalt=idsalt_beg,idsalt_end
    trcSalt_solml2(idsalt,0,NY,NX)=trcSalt_solml2(idsalt,0,NY,NX) &
      +trcSalt_TQR(idsalt,NY,NX)+trcSalt3DFlo2CellM(idsalt,3,0,NY,NX)
  ENDDO
  end subroutine UpdateSoluteInResidue
!------------------------------------------------------------------------------------------

  subroutine UpdateSoluteInMicMacpores(NY,NX)
!
!     Description:
!
  implicit none
  integer, intent(in) :: NY,NX
  integer :: L,idsalt
!     begin_execution
!
  D9685: DO L=NU(NY,NX),NL(NY,NX)
    IF(VLSoilPoreMicP_vr(L,NY,NX).GT.ZEROS2(NY,NX))THEN
      DO idsalt=idsalt_beg,idsaltb_end
        trcSalt_solml2(idsalt,L,NY,NX)=trcSalt_solml2(idsalt,L,NY,NX) &
          +trcSalt_Flo2MicP_vr(idsalt,L,NY,NX)+trcSalt_RFXS(idsalt,L,NY,NX) &
          +trcSalt_RFLZ(idsalt,L,NY,NX)
        trcSalt_soHml2(idsalt,L,NY,NX)=trcSalt_soHml2(idsalt,L,NY,NX) &
          +trcSalt_Flo2MacP_vr(idsalt,L,NY,NX)-trcSalt_RFXS(idsalt,L,NY,NX)
      ENDDO

    ENDIF
  ENDDO D9685
  end subroutine UpdateSoluteInMicMacpores

!------------------------------------------------------------------------------------------

  subroutine InitFluxArrays(NHW,NHE,NVN,NVS)

  implicit none
  integer, intent(in) :: NHW,NHE,NVN,NVS

  integer :: NY,NX
  DO NX=NHW,NHE
    DO  NY=NVN,NVS
!
!     INITIALIZE SOLUTE RUNOFF NET FLUX ACCUMULATORS
      call InitFluxAccumlatorsByRunoff(NY,NX)
!
!     INITIALIZE SNOWPACK NET FLUX ACCUMULATORS
      call InitFluxAccumulatorsInSnowpack(NY,NX)
!
!     INITIALIZE SOIL SOLUTE NET FLUX ACCUMULATORS
      call InitFluxAccumulatorsInSoil(NY,NX)
    ENDDO
  ENDDO
  end subroutine InitFluxArrays


!------------------------------------------------------------------------------------------

  subroutine InitFluxAccumlatorsByRunoff(NY,NX)
!
!     Description:
!
  implicit none
  integer, intent(in) :: NY,NX
!     begin_execution
!
  trcSalt_TQR(idsalt_beg:idsalt_end,NY,NX)=0.0

  trcSalt_TQ(idsalt_beg:idsalt_end,NY,NX)=0.0
  end subroutine InitFluxAccumlatorsByRunoff
!------------------------------------------------------------------------------------------

  subroutine InitFluxAccumulatorsInSoil(NY,NX)
!
!     Description:
!
  implicit none
  integer, intent(in) :: NY,NX
  integer :: L

  D9885: DO L=NU(NY,NX),NL(NY,NX)
    trcSalt_Flo2MicP_vr(idsalt_beg:idsaltb_end,L,NY,NX)=0.0_r8
    trcSalt_Flo2MacP_vr(idsalt_beg:idsaltb_end,L,NY,NX)=0.0_r8
  ENDDO D9885
  end subroutine InitFluxAccumulatorsInSoil

end module TranspSaltMod
