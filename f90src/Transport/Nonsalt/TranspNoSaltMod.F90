module TranspNoSaltMod
!!
! Description:
!
  use data_kind_mod,    only: r8 => DAT_KIND_R8
  use abortutils,       only: destroy,endrun
  USE MiniMathMod,      ONLY: AZMAX1, fixnegmass, flux_mass_limiter,AZERO
  use TracerPropMod,    only: MolecularWeight
  use EcoSiMParDataMod, only: micpar
  use SnowTransportMod, only: TracerFall2Snowpack
  use TranspNoSaltFastMod
  use InitNoSaltTransportMod
  use TranspNoSaltSlowMod
  use DebugToolMod
  use SOMDataType
  use ChemTranspDataType
  use GridConsts
  use EcoSIMSolverPar
  use SoilPhysDataType
  use SoilHeatDatatype
  use SoilWaterDataType
  use EcoSIMCtrlDataType
  use SoilBGCDataType
  use ClimForcDataType
  use FertilizerDataType
  use EcosimConst
  use SurfLitterDataType
  use SnowDataType
  use SurfSoilDataType
  use LandSurfDataType
  use RootDataType
  use AqueChemDatatype
  use SoilPropertyDataType
  use PlantDataRateType
  use GridDataType
  use TranspNoSaltDataMod
  use IrrigationDataType
  implicit none

  private
  CHARACTER(LEN=*), PARAMETER :: MOD_FILENAME=&
  __FILE__

  real(r8), parameter :: tinyval=1.e-13_r8
  public :: TranspNoSalt
  public :: InitTranspNoSalt,DestructTranspNoSalt
  contains

  subroutine InitTranspNoSalt()
  implicit none

  call InitTransfrData

  end subroutine InitTranspNoSalt
!------------------------------------------------------------------------------------------
  subroutine DestructTranspNoSalt
  implicit none

  call DestructTransfrData

  end subroutine DestructTranspNoSalt

!------------------------------------------------------------------------------------------

  SUBROUTINE TranspNoSalt(I,J,NHW,NHE,NVN,NVS)
!
!     THIS SUBROUTINE CALCULATES 3-DIMENSIONAL FLUXES OF ALL SOIL
!     NON-SALT SOLUTES AND GASES
!
  implicit none
  integer, intent(in) :: I, J
  integer, intent(in) :: NHW,NHE,NVN,NVS

  character(len=*), parameter :: subname='TranspNoSalt'
  integer :: MX,MM
  integer :: M

!     execution begins here
  call PrintInfo('beg '//subname)
!
!
  call InitTranspNoSaltModel(I,J,NHW,NHE,NVN,NVS)

!
! TIME STEP USED IN GAS AND SOLUTE FLUX CALCULATIONS
! NPH=NPX,NPT=NPY
! NPH=no. of cycles per hour for water, heat and solute flux calculations
! NPG=NPG=NPH*NPT,number of cycles per hour for gas flux calculations, 
! dt_GasCyc=1.0_r8/NPT, number of cycles for gas flux calculations for 1 iteration in NPH,

  DO M=1,NPH

    !transport without reaction rates, update trcs_solml2_vr
    call TransptSlowNoSaltM(I,J,M,NHE,NHW,NVS,NVN)

    !DO gas dissolution-volatization, and solute uptake in each iteration MM        
    call PassSlow2FastIterationM(I,J,M,NHE,NHW,NVS,NVN)

    DO MM=1,NPT      
      call TransptFastNoSaltMM(I,J,M,NHE,NHW,NVS,NVN)
    ENDDO
    !Do bubbling
    call BubbleEffluxM(I,J,M,NHE,NHW,NVS,NVN)
  ENDDO

  call PrintInfo('end '//subname)
  
  END subroutine TranspNoSalt

!------------------------------------------------------------------------------------------
  subroutine PassSlow2FastIterationM(I,J,M,NHE,NHW,NVS,NVN)
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: M
  integer, intent(in) :: NHE,NHW,NVS,NVN
  character(len=*), parameter :: subname='PassSlow2FastIterationM'
  integer :: NY,NX,L,N
  real(r8) :: PARGM

  call PrintInfo('beg '//subname)

  DO NX=NHW,NHE
    DO  NY=NVN,NVS
      !
      ! GASEOUS BOUNDARY LAYER CONDUCTANCES
      !
      ! PARG=boundary layer conductance above soil surface from watsub.f
      !
      PARGM                       = CondGasXSnowM_col(M,NY,NX)*dt_GasCyc
      PARGas_CefMM(idg_CO2,NY,NX) = PARGM*0.74_r8
      PARGas_CefMM(idg_CH4,NY,NX) = PARGM*1.04_r8
      PARGas_CefMM(idg_O2,NY,NX)  = PARGM*0.83_r8
      PARGas_CefMM(idg_N2,NY,NX)  = PARGM*0.86_r8
      PARGas_CefMM(idg_N2O,NY,NX) = PARGM*0.74_r8
      PARGas_CefMM(idg_NH3,NY,NX) = PARGM*1.02_r8
      PARGas_CefMM(idg_H2,NY,NX)  = PARGM*2.08_r8
      PARGas_CefMM(idg_AR,NY,NX)  = PARGM*0.72_r8


      RBGCSinkGasMM_vr(idg_O2,0,NY,NX)=RO2UptkSoilM_vr(M,0,NY,NX)*dt_GasCyc
      L=NU(NY,NX)
      WaterFlow2SoilMM_3D(3,L,NY,NX)     = (WaterFlow2MicPM_3D(M,3,L,NY,NX)+WaterFlow2MacPM_3D(M,3,L,NY,NX))*dt_GasCyc
      tScalReductVLsoiAirPMM_vr(L,NY,NX) = ReductVLsoiAirPM_vr(M,L,NY,NX)*dt_GasCyc

      VLWatMicPMA_vr(L,NY,NX)          = VLWatMicPM_vr(M,L,NY,NX)*trcs_VLN_vr(ids_NH4,L,NY,NX)
      VLWatMicPMB_vr(L,NY,NX)          = VLWatMicPM_vr(M,L,NY,NX)*trcs_VLN_vr(ids_NH4B,L,NY,NX)
      VLWatMicPXA_vr(L,NY,NX)          = natomw*VLWatMicPMA_vr(L,NY,NX)
      VLWatMicPXB_vr(L,NY,NX)          = natomw*VLWatMicPMB_vr(L,NY,NX)
      VLsoiAirPMA_vr(L,NY,NX)          = VLsoiAirPM_vr(M,L,NY,NX)*trcs_VLN_vr(ids_NH4,L,NY,NX)
      VLsoiAirPMB_vr(L,NY,NX)          = VLsoiAirPM_vr(M,L,NY,NX)*trcs_VLN_vr(ids_NH4B,L,NY,NX)
      RBGCSinkGasMM_vr(idg_O2,L,NY,NX) = RO2UptkSoilM_vr(M,L,NY,NX)*dt_GasCyc-trcs_deadroot2soil_vr(idg_O2,L,NY,NX)*dts_gas

      DO L=NU(NY,NX)+1,NL(NY,NX)
        tScalReductVLsoiAirPMM_vr(L,NY,NX) = ReductVLsoiAirPM_vr(M,L,NY,NX)*dt_GasCyc
        DO N=FlowDirIndicator_col(NY,NX),3
          WaterFlow2SoilMM_3D(N,L,NY,NX)=(WaterFlow2MicPM_3D(M,N,L,NY,NX)+WaterFlow2MacPM_3D(M,N,L,NY,NX))*dt_GasCyc
        ENDDO

        VLWatMicPMA_vr(L,NY,NX) = VLWatMicPM_vr(M,L,NY,NX)*trcs_VLN_vr(ids_NH4,L,NY,NX)
        VLWatMicPMB_vr(L,NY,NX) = VLWatMicPM_vr(M,L,NY,NX)*trcs_VLN_vr(ids_NH4B,L,NY,NX)
        VLWatMicPXA_vr(L,NY,NX) = natomw*VLWatMicPMA_vr(L,NY,NX)
        VLWatMicPXB_vr(L,NY,NX) = natomw*VLWatMicPMB_vr(L,NY,NX)

        VLsoiAirPMA_vr(L,NY,NX)          = VLsoiAirPM_vr(M,L,NY,NX)*trcs_VLN_vr(ids_NH4,L,NY,NX)
        VLsoiAirPMB_vr(L,NY,NX)          = VLsoiAirPM_vr(M,L,NY,NX)*trcs_VLN_vr(ids_NH4B,L,NY,NX)
        RBGCSinkGasMM_vr(idg_O2,L,NY,NX) = RO2UptkSoilM_vr(M,L,NY,NX)*dt_GasCyc-trcs_deadroot2soil_vr(idg_O2,L,NY,NX)*dts_gas
      ENDDO
    ENDDO
  ENDDO

  call PrintInfo('end '//subname)
  end subroutine PassSlow2FastIterationM

end module TranspNoSaltMod
