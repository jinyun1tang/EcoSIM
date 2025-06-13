module TranspNoSaltMod
!!
! Description:
!
  use data_kind_mod,    only: r8 => DAT_KIND_R8
  use abortutils,       only: destroy, endrun
  USE MiniMathMod,      ONLY: AZMAX1,  fixnegmass, flux_mass_limiter, AZERO
  use TracerPropMod,    only: MolecularWeight
  use EcoSiMParDataMod, only: micpar
  use EcoSIMCtrlMod,    only: iVerbLevel
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

  subroutine EnterMassCheck(NHW,NHE,NVN,NVS)
  implicit none
  integer, intent(in) :: NHW,NHE,NVN,NVS
  integer :: NY,NX,idg,L

  DO NX=NHW,NHE
    DO  NY=NVN,NVS

      trcs_mass_beg(:,NY,NX)=0._r8
      DO idg=idg_beg,idg_NH3

        trcs_mass_litr(idg,NY,NX)=trcs_solml_vr(idg,0,NY,NX)
        trcs_mass_snow(idg,NY,NX)=0._r8
        trcs_mass_soil(idg,NY,NX)=0._r8
        
        DO L=1,nsnol_col(NY,NX)
          trcs_mass_snow(idg,NY,NX)=trcs_mass_snow(idg,NY,NX)+trcg_solsml_snvr(idg,L,NY,NX)
        ENDDO

        trcs_mass_beg(idg,NY,NX)=trcs_mass_litr(idg,NY,NX)+trcs_mass_snow(idg,NY,NX)
        
        DO L=NU_col(NY,NX),NL_col(NY,NX)
          trcs_mass_beg(idg,NY,NX)=trcs_mass_beg(idg,NY,NX)+trcs_solml_vr(idg,L,NY,NX)+trcg_gasml_vr(idg,L,NY,NX) &
            + trcs_soHml_vr(idg,L,NY,NX)
          trcs_mass_soil(idg,NY,NX)=trcs_mass_soil(idg,NY,NX)+trcs_solml_vr(idg,L,NY,NX)+trcg_gasml_vr(idg,L,NY,NX) &
            + trcs_soHml_vr(idg,L,NY,NX)  
        ENDDO        
      ENDDO

      idg=idg_NH3
      trcs_mass_litr(idg,NY,NX)=trcs_mass_litr(idg,NY,NX)+trcs_solml_vr(idg_NH3B,0,NY,NX)
      trcs_mass_beg(idg,NY,NX)=trcs_mass_beg(idg,NY,NX)+trcs_solml_vr(idg_NH3B,0,NY,NX)
      DO L=NU_col(NY,NX),NL_col(NY,NX)
        trcs_mass_beg(idg,NY,NX)=trcs_mass_beg(idg,NY,NX)+trcs_solml_vr(idg_NH3B,L,NY,NX) &
          + trcs_soHml_vr(idg_NH3B,L,NY,NX)
        trcs_mass_soil(idg,NY,NX)=trcs_mass_soil(idg,NY,NX)+trcs_solml_vr(idg_NH3B,L,NY,NX) &
          + trcs_soHml_vr(idg_NH3B,L,NY,NX)  
      ENDDO

    ENDDO
  ENDDO  
  end subroutine EnterMassCheck  
!------------------------------------------------------------------------------------------

  subroutine ExitMassCheck(I,J,NHW,NHE,NVN,NVS,M)
  implicit none
  integer, intent(in) :: I,J,NHW,NHE,NVN,NVS
  integer, optional, intent(in) :: M
  integer :: NY,NX,idg,L  
  real(r8) :: errmass,mass_litr,mass_snow,mass_soil
  real(r8) :: netflx2soil
  real(r8) :: trcs_mass_now(idg_beg:idg_NH3)

  DO NX=NHW,NHE
    DO  NY=NVN,NVS

      trcs_mass_now(:)=0._r8
      DO idg=idg_beg,idg_NH3

        SurfGasEmiss_flx_col(idg,NY,NX) =  trcg_ebu_flx_col(idg,NY,NX) &
          + GasDiff2Surf_flx_col(idg,NY,NX)+Gas_WetDeposit_flx_col(idg,NY,NX) 

        mass_litr=trcs_solml_vr(idg,0,NY,NX)
        mass_snow=0._r8
        DO L=1,nsnol_col(NY,NX)
          mass_snow=mass_snow+trcg_solsml_snvr(idg,L,NY,NX)
        ENDDO        

        mass_soil=0._r8
        DO L=NU_col(NY,NX),NL_col(NY,NX)
          mass_soil=mass_soil+trcg_gasml_vr(idg,L,NY,NX)+trcs_solml_vr(idg,L,NY,NX)+trcs_soHml_vr(idg,L,NY,NX)  
        ENDDO

        if(idg==idg_NH3)then
          !also include NH3B
          mass_litr=mass_litr+trcs_solml_vr(idg_NH3B,0,NY,NX)
          DO L=NU_col(NY,NX),NL_col(NY,NX)
            mass_soil=mass_soil+trcs_solml_vr(idg_NH3B,L,NY,NX)+trcs_soHml_vr(idg_NH3B,L,NY,NX)
          ENDDO
        endif

        trcs_mass_now(idg)=mass_snow+mass_litr+mass_soil
        errmass=trcs_mass_beg(idg,NY,NX)-trcs_mass_now(idg)  &
          +SurfGasEmiss_flx_col(idg,NY,NX)+GasHydroLoss_flx_col(idg,NY,NX) &
          +RGasNetProd_col(idg,NY,NX)-trcs_Soil2plant_uptake_col(idg,NY,NX)

        if(idg==idg_NH3)then
          errmass=errmass-trcs_Soil2plant_uptake_col(idg_NH3B,NY,NX)
        endif
        if(abs(errmass)>1.e-5)then
          if(iVerbLevel==1 .or. abs(errmass)>1.e-4)then
            write(121,*)('-',L=1,50)
            if(present(M))then
              write(121,*)(I*1000+J)*10+M,trcs_names(idg),'total'
            else
              write(121,*)(I*1000+J)*10,trcs_names(idg),'total'
            endif
            write(121,*)'errmass',errmass, 'NY NX=',NY,NX
            write(121,*)'mass beg, end, delta',trcs_mass_beg(idg,NY,NX),trcs_mass_now(idg),trcs_mass_now(idg)-trcs_mass_beg(idg,NY,NX)
            write(121,*)'beg sno,litr,soil',trcs_mass_snow(idg,NY,NX),trcs_mass_litr(idg,NY,NX),trcs_mass_soil(idg,NY,NX)
            write(121,*)'mass sno,litr,soil',mass_snow,mass_litr,mass_soil
            write(121,*)'Delta mass sno, litr,soil',mass_snow-trcs_mass_snow(idg,NY,NX),mass_litr-trcs_mass_litr(idg,NY,NX), &
              mass_soil-trcs_mass_soil(idg,NY,NX)
            write(121,*)'mass snow',trcg_solsml_snvr(idg,1:nsnol_col(NY,NX),NY,NX)  
            write(121,*)'dif/wetdepo to litr =',GasDiff2Litr_flx_col(idg,NY,NX),Gas_WetDepo2Litr_col(idg,NY,NX)
            write(121,*)'gas litr2soil       =',Gas_litr2Soil_flx_col(idg,NY,NX)  
            write(121,*)'snow2litr           =',trcg_AquaADV_Snow2Litr_flx(idg,NY,NX)
            write(121,*)'netflx2litr         =',trcg_AquaADV_Snow2Litr_flx(idg,NY,NX)+GasDiff2Litr_flx_col(idg,NY,NX) &
              +Gas_WetDepo2Litr_col(idg,NY,NX)-Gas_litr2Soil_flx_col(idg,NY,NX)+GasHydroLoss_litr_flx_col(idg,NY,NX) 

            netflx2soil=GasDiff2Soil_flx_col(idg,NY,NX)+Gas_WetDepo2Soil_col(idg,NY,NX) &
              +Gas_litr2Soil_flx_col(idg,NY,NX)+trcg_ebu_flx_col(idg,NY,NX) &
              +trcg_AquaADV_Snow2Soil_flx(idg,NY,NX)+trcs_SubsurTransp_flx_2DH(idg,NY,NX)  &
              +trcs_irrig_flx_col(idg,NY,NX)-trcs_drainage_flx_col(idg,NY,NX)

            if(idg==idg_NH3)then
              netflx2soil=netflx2soil+trcn_AquaADV_Snow2Band_flx(idg_NH3B,NY,NX)+trcs_SubsurTransp_flx_2DH(idg_NH3B,NY,NX) &
                +trcs_irrig_flx_col(idg_NH3B,NY,NX)-trcs_drainage_flx_col(idg_NH3B,NY,NX)
              write(121,*)'topsoil, irrig flx = ',trc_topsoil_flx_col(idg,NY,NX),trcs_irrig_flx_col(idg,NY,NX) &
                +trcs_irrig_flx_col(idg_NH3B,NY,NX) 
            else 
              write(121,*)'topsoil, irrig flx = ',trc_topsoil_flx_col(idg,NY,NX),trcs_irrig_flx_col(idg,NY,NX)               
            endif    
            write(121,*)'netflx2soil        =',netflx2soil, TranspNetSoil_flx_col(idg,NY,NX)+trcg_ebu_flx_col(idg,NY,NX)
              
            write(121,*)'err litr           =', trcg_AquaADV_Snow2Litr_flx(idg,NY,NX)+GasDiff2Litr_flx_col(idg,NY,NX) &
              +Gas_WetDepo2Litr_col(idg,NY,NX)-Gas_litr2Soil_flx_col(idg,NY,NX)+GasHydroLoss_litr_flx_col(idg,NY,NX) &
              -mass_litr+trcs_mass_litr(idg,NY,NX) 
            write(121,*)'err soil           =',GasDiff2Soil_flx_col(idg,NY,NX) &
              +Gas_WetDepo2Soil_col(idg,NY,nX)+Gas_litr2Soil_flx_col(idg,NY,NX)+trcg_AquaADV_Snow2Soil_flx(idg,NY,NX)&
              -trcs_drainage_flx_col(idg,NY,NX)+trcs_SubsurTransp_flx_2DH(idg,NY,NX)+trcg_ebu_flx_col(idg,NY,NX) &
              - mass_soil+trcs_mass_soil(idg,NY,NX)&
              + trcs_irrig_flx_col(idg,NY,NX)+RGasNetProdSoil_col(idg,NY,NX)

            write(121,*)'nsnol_col           =',nsnol_col(NY,NX)
            write(121,*)'err snow            =',Gas_WetDepo2Snow_col(idg,NY,NX)-Gas_Snowloss_flx_col(idg,NY,NX)-mass_snow+trcs_mass_snow(idg,NY,NX)
            write(121,*)'snow depo loss drift=',Gas_WetDepo2Snow_col(idg,NY,NX),Gas_Snowloss_flx_col(idg,NY,NX),trcg_SnowDrift_flx_col(idg,NY,NX)
            write(121,*)'ntransp2soil        =',TranspNetSoil_flx_col(idg,NY,NX),TranspNetSoil_flx2_col(idg,NY,NX),TranspNetSoil_flx_col(idg,NY,NX)-TranspNetSoil_flx2_col(idg,NY,NX)              
            write(121,*)'drain               =',trcs_drainage_flx_col(idg,NY,NX)  
            write(121,*)'ebu                 =',trcg_ebu_flx_col(idg,NY,NX)
            write(121,*)'dif                 =',GasDiff2Surf_flx_col(idg,NY,NX)
            write(121,*)'soil2roots          =',trcs_Soil2plant_uptake_col(idg,NY,NX)
            write(121,*)'--------------------------'
            write(121,*)'surfemis            =',SurfGasEmiss_flx_col(idg,NY,NX)          
            write(121,*)'hydloss             =',GasHydroLoss_flx_col(idg,NY,NX)
            write(121,*)'netpro              =',RGasNetProd_col(idg,NY,NX),RGasNetProdSoil_col(idg,NY,NX)          
            write(121,*)'--------------------------'          
            write(121,*)'latloss             =',trcs_SubsurTransp_flx_2DH(idg,NY,NX)
            write(121,*)'wetdep              =',Gas_WetDeposit_flx_col(idg,NY,NX)          
            write(121,*)'sol                 =',trcs_solml_vr(idg,NU_col(NY,NX):NL_col(NY,NX),NY,NX)
            write(121,*)'gas                 =',trcg_gasml_vr(idg,NU_col(NY,NX):NL_col(NY,NX),NY,NX)          
            write(121,*)'Hml                 =',trcs_soHml_vr(idg,NU_col(NY,NX):NL_col(NY,NX),NY,NX)  
            if(idg==idg_NH3)then
              write(121,*)'sol_NH3B            =',trcs_solml_vr(idg_NH3B,NU_col(NY,NX):NL_col(NY,NX),NY,NX)
              write(121,*)'Hml_NH3B            =',trcs_soHml_vr(idg_NH3B,NU_col(NY,NX):NL_col(NY,NX),NY,NX)  
            endif
            write(121,*)'transp',NU_col(NY,NX),NL_col(NY,NX)
            write(121,*)transp_diff_slow_vr(idg,NU_col(NY,NX):NL_col(NY,NX),NY,NX)
          endif
          if(abs(errmass)>1.e-4) &
          call endrun(trim(mod_filename)//' at line',__LINE__)
        endif
      ENDDO
    ENDDO
  ENDDO  
  end subroutine ExitMassCheck

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

  call EnterMassCheck(NHW,NHE,NVN,NVS)
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
      call TransptFastNoSaltMM(I,J,M,MM,NHE,NHW,NVS,NVN)
    ENDDO
    
    !Do bubbling
    call BubbleEffluxM(I,J,M,NHE,NHW,NVS,NVN)

  ENDDO

  call BackCopyStateVars(I,J,NHE,NHW,NVS,NVN)
      
  call ExitMassCheck(I,J,NHW,NHE,NVN,NVS)

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

      RBGCSinkGasMM_vr(idg_O2,0,NY,NX) = RO2UptkSoilM_vr(M,0,NY,NX)*dt_GasCyc
      VLWatMicPXA_vr(0,NY,NX)          = natomw*VLWatMicPM_vr(M,0,NY,NX)

      L=NU_col(NY,NX)
      WaterFlow2SoilMM_3D(3,L,NY,NX)     = (WaterFlow2MicPM_3D(M,3,L,NY,NX)+WaterFlow2MacPM_3D(M,3,L,NY,NX))*dt_GasCyc
      tScalReductVLsoiAirPMM_vr(L,NY,NX) = ReductVLsoiAirPM_vr(M,L,NY,NX)*dt_GasCyc

      VLWatMicPMA_vr(L,NY,NX)          = VLWatMicPM_vr(M,L,NY,NX)*trcs_VLN_vr(ids_NH4,L,NY,NX)
      VLWatMicPMB_vr(L,NY,NX)          = VLWatMicPM_vr(M,L,NY,NX)*trcs_VLN_vr(ids_NH4B,L,NY,NX)
      VLWatMicPXA_vr(L,NY,NX)          = natomw*VLWatMicPMA_vr(L,NY,NX)
      VLWatMicPXB_vr(L,NY,NX)          = natomw*VLWatMicPMB_vr(L,NY,NX)
      VLsoiAirPMA_vr(L,NY,NX)          = VLsoiAirPM_vr(M,L,NY,NX)*trcs_VLN_vr(ids_NH4,L,NY,NX)
      VLsoiAirPMB_vr(L,NY,NX)          = VLsoiAirPM_vr(M,L,NY,NX)*trcs_VLN_vr(ids_NH4B,L,NY,NX)
      RBGCSinkGasMM_vr(idg_O2,L,NY,NX) = RO2UptkSoilM_vr(M,L,NY,NX)*dt_GasCyc-trcs_deadroot2soil_vr(idg_O2,L,NY,NX)*dts_gas

      DO L=NU_col(NY,NX)+1,NL_col(NY,NX)
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
