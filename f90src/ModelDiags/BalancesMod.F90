module BalancesMod
  use data_kind_mod,     only: r8 => DAT_KIND_R8
  use CanopyDataType,    only: QVegET_col
  use GridDataType,      only: NU, NL
  use EcoSimConst,       only: DENSICE
  use abortutils,        only: endrun
  use EcoSIMCtrlMod,     only: fixWaterLevel
  use AqueChemDatatype,  only: trcg_mass_cumerr_col
  use PlantMgmtDataType, only: iDayPlanting_pft, iDayPlantHarvest_pft
  use RootDataType
  use SoilBGCDataType
  use GridDataType
  use SurfLitterDataType
  use CanopyDataType
  use SnowDataType
  use BalanceCheckDataType
  use SoilWaterDataType
  use ClimForcDataType
  use SurfSoilDataType
  use EcosimBGCFluxType
  use SoilHeatDataType
  use PlantDataRateType
  use EcoSimSumDataType
implicit none
  private
  character(len=*), parameter :: mod_filename = &
  __FILE__
  public :: BegCheckBalances
  public :: EndCheckBalances
  public :: SummarizeTracerMass
  public :: SummarizeSnowMass
  public :: SummarizeTracers
contains

  subroutine BegCheckBalances(I,J,NHW,NHE,NVN,NVS)  
  !
  !Description
  !Prepare for next mass balance check
  implicit none
  integer, intent(in) :: I,J,NHW,NHE,NVN,NVS
  integer :: NY,NX,idg
  
  DO  NX=NHW,NHE
    DO  NY=NVN,NVS
      WaterErr_col(NY,NX)           = WatMass_col(NY,NX)
      HeatErr_col(NY,NX)            = HeatStore_col(NY,NX)
      SnowEngyBeg_col(NY,NX)        = SnowEngyEnd_col(NY,NX)
      CanopyWaterMassBeg_col(NY,NX) = CanopyWaterMassEnd_col(NY,NX)
      SnowMassBeg_col(NY,NX)        = SnowMassEnd_col(NY,NX)
      LitWatMassBeg_col(NY,NX)      = LitWatMassEnd_col(NY,NX)
      SoilWatMassBeg_col(NY,NX)     = SoilWatMassEnd_col(NY,NX)

      DO idg=idg_beg,idg_end
        trcg_TotalMass_beg_col(idg,NY,NX) = trcg_TotalMass_col(idg,NY,NX)
        trcg_soilMass_beg_col(idg,NY,NX) = trcg_soilMass_col(idg,NY,NX)
      enddo
      DO idg=idg_beg,idg_NH3      
        trcg_rootMass_beg_col(idg,NY,NX)  = trcg_rootMass_col(idg,NY,NX)  
        trcg_snowMass_beg_col(idg,NY,NX)  = trcg_snowMass_col(idg,NY,NX)
      enddo  
!      if(I==140 .and. J>=20)then
!        write(116,*)'==============================================='
!        write(116,*)trcg_rootMass_beg_col(idg_N2,NY,NX),I*1000+J,'beg'
!      endif
    ENDDO
  ENDDO  
  end subroutine BegCheckBalances
!------------------------------------------------------------------------------------------
  subroutine SummarizeStorage(I,J,NHW,NHE,NVN,NVS)  

  implicit none
  integer, intent(in) :: I,J,NHW,NHE,NVN,NVS
  integer :: NY,NX

  WatMassStore_lnd = 0._r8
  HeatStore_lnd    = 0._r8
  DO NX = NHW, NHE
    DO  NY=NVN,NVS

      call SumUpWaterStorage(NY,NX)

      call SumUpHeatStorage(NY,NX)

    ENDDO
  ENDDO
  end subroutine SummarizeStorage
!------------------------------------------------------------------------------------------

  subroutine SummarizeSnowMass(NY,NX)
  implicit none
  integer, intent(in) :: NY,NX
  integer :: L
  SnowMassEnd_col(NY,NX)=0._r8
  DO  L=1,JS
    if(AMAX1(VLDrySnoWE_snvr(L,NY,NX),VLWatSnow_snvr(L,NY,NX),VLIceSnow_snvr(L,NY,NX))<=0._r8)cycle
    SnowMassEnd_col(NY,NX) = SnowMassEnd_col(NY,NX)+VLDrySnoWE_snvr(L,NY,NX)+VLWatSnow_snvr(L,NY,NX)+VLIceSnow_snvr(L,NY,NX)*DENSICE
  ENDDO

  end subroutine SummarizeSnowMass
!------------------------------------------------------------------------------------------
  subroutine SumUpWaterStorage(NY,NX)

  implicit none
  integer, intent(in) :: NY,NX
  integer :: L

  !canopy
  CanopyWaterMassEnd_col(NY,NX) = WatHeldOnCanopy_col(NY,NX)+CanopyWat_col(NY,NX)
  !snow  
  SnowMassEnd_col(NY,NX)=0._r8
  DO  L=1,JS
    if(AMAX1(VLDrySnoWE_snvr(L,NY,NX),VLWatSnow_snvr(L,NY,NX),VLIceSnow_snvr(L,NY,NX))<=0._r8)cycle
    SnowMassEnd_col(NY,NX) = SnowMassEnd_col(NY,NX)+VLDrySnoWE_snvr(L,NY,NX)+VLWatSnow_snvr(L,NY,NX)+VLIceSnow_snvr(L,NY,NX)*DENSICE
  ENDDO
  !litter
  LitWatMassEnd_col(NY,NX) = VLWatMicP_vr(0,NY,NX)+VLiceMicP_vr(0,NY,NX)*DENSICE

  WatMass_col(NY,NX)       = CanopyWaterMassEnd_col(NY,NX)+SnowMassEnd_col(NY,NX)+LitWatMassEnd_col(NY,NX)
  SoilWatMassEnd_col(NY,NX)= 0._r8
  !add water in soil 
  DO L=NU(NY,NX),NL(NY,NX)
    SoilWatMassEnd_col(NY,NX) = SoilWatMassEnd_col(NY,NX)+ VLWatMicP_vr(L,NY,NX)+VLWatMacP_vr(L,NY,NX)+(VLiceMicP_vr(L,NY,NX)+VLiceMacP_vr(L,NY,NX))*DENSICE
  ENDDO

  WatMass_col(NY,NX)=WatMass_col(NY,NX)+SoilWatMassEnd_col(NY,NX)

  WatMassStore_lnd = WatMassStore_lnd+WatMass_col(NY,NX)

  end subroutine SumUpWaterStorage

!------------------------------------------------------------------------------------------
  subroutine SumUpHeatStorage(NY,NX)

  implicit none
  integer, intent(in) :: NY,NX
  integer :: L

  HeatStore_col(NY,NX) = 0._r8
  !add snow
  DO  L=1,JS
    if(VLHeatCapSnow_snvr(L,NY,NX)<=0._r8)cycle
    HeatStore_col(NY,NX) = HeatStore_col(NY,NX)+ VLHeatCapSnow_snvr(L,NY,NX)*TKSnow_snvr(L,NY,NX)    
  ENDDO
  SnowEngyEnd_col(NY,NX) = HeatStore_col(NY,NX)
  !add canopy + litter
  HeatStore_col(NY,NX)   = HeatStore_col(NY,NX)+CanopyHeatStor_col(NY,NX) + TKS_vr(0,NY,NX)*VHeatCapacity_vr(0,NY,NX)

  !add soil
  DO L=NU(NY,NX),NL(NY,NX)
    HeatStore_col(NY,NX) = HeatStore_col(NY,NX)+ VHeatCapacity_vr(L,NY,NX)*TKS_vr(L,NY,NX)
  ENDDO

  HeatStore_lnd = HeatStore_lnd+HeatStore_col(NY,NX)

  end subroutine SumUpHeatStorage
!------------------------------------------------------------------------------------------

  subroutine EndCheckBalances(I,J,NHW,NHE,NVN,NVS)  
  !
  !Description
  !Do mass balance check
  implicit none
  integer, intent(in) :: I,J,NHW,NHE,NVN,NVS
  integer :: NY,NX
  real(r8), parameter :: err_h2o=1.e-5_r8
  real(r8), parameter :: err_engy=1.e-6_r8
  real(r8) :: WaterErr_test
  real(r8) :: HeatErr_test
  real(r8) :: literH2Oerr_test
  real(r8) :: canopyH2Oerr_test
  real(r8) :: SnowMassErr_test
  real(r8) :: SoilWatErr_test
  real(r8) :: precipErr_test
  real(r8) :: prec2expSErr_test
  real(r8) :: prec2SnoErr_test
  real(r8) :: tracer_mass_err
  real(r8) :: tracer_rootmass_err
  real(r8) :: tracer_snowmass_err
  integer :: ii,idg

  call SummarizeStorage(I,J,NHW,NHE,NVN,NVS)

  call SummarizeTracers(I,J,NHW,NHE,NVN,NVS)

  DO  NX=NHW,NHE
    DO  NY=NVN,NVS

      if(fixWaterLevel)then
        SoilWatErr_test=SoilWatMassBeg_col(NY,NX)-SoilWatMassEnd_col(NY,NX)      
      else
        SoilWatErr_test=SoilWatMassBeg_col(NY,NX)-SoilWatMassEnd_col(NY,NX)+Qinflx2Soil_col(NY,NX) &
          -QDrain_col(NY,NX)-QDischar_col(NY,NX)+TPlantRootH2OUptake_col(NY,NX)+QWatIntLaterFlow_col(NY,NX)
      endif
      precipErr_test = RainPrecThrufall_col(NY,NX)-Rain2LitR_col(NY,NX)-Rain2Soil_col(NY,NX)-RainPrec2Sno_col(NY,NX)
      prec2expSErr_test=Rain2ExposedSurf_col(NY,NX)-Rain2LitR_col(NY,NX)-Rain2Soil_col(NY,NX)
      prec2SnoErr_test= Prec2Snow_col(NY,NX)-RainPrec2Sno_col(NY,NX)-SnoFalPrec_col(NY,NX)  

      literH2Oerr_test = LitWatMassBeg_col(NY,NX)-LitWatMassEnd_col(NY,NX)+QRunSurf_col(NY,NX)+WatFLo2LitR_col(NY,NX)

      SnowMassErr_test = SnowMassBeg_col(NY,NX)-SnowMassEnd_col(NY,NX)+Prec2Snow_col(NY,NX)-QSnowH2Oloss_col(NY,NX)

      canopyH2Oerr_test=CanopyWaterMassBeg_col(NY,NX)-CanopyWaterMassEnd_col(NY,NX)+PrecIntceptByCanopy_col(NY,NX) &
        +QVegET_col(NY,NX)-TPlantRootH2OUptake_col(NY,NX)-QCanopyWat2Dist_col(NY,NX)

      if(fixWaterLevel)then
        WaterErr_test = WaterErr_col(NY,NX)-WatMass_col(NY,NX)
      else      
        WaterErr_test = WaterErr_col(NY,NX)-WatMass_col(NY,NX)+PrecAtm_col(NY,NX)+Irrigation_col(NY,NX)+QWatIntLaterFlow_col(NY,NX) &
          +RainLitr_col(NY,NX)+VapXAir2GSurf_col(NY,NX)+QVegET_col(NY,NX)+QRunSurf_col(NY,NX) &
          -QDrain_col(NY,NX)-QDischar_col(NY,NX)+TPlantRootH2OUptake_col(NY,NX)-QCanopyWat2Dist_col(NY,NX)
      endif
      HeatErr_test = HeatErr_col(NY,NX)-HeatStore_col(NY,NX)+THeatRootRelease_col(NY,NX) &
        +HeatSource_col(NY,NX)+Eco_NetRad_col(NY,NX)+Eco_Heat_Latent_col(NY,NX)+Eco_Heat_Sens_col(NY,NX)&
        +PrecHeat_col(NY,NX)+THeatSoiThaw_col(NY,NX)+THeatSnowThaw_col(NY,NX)+HeatRunSurf_col(NY,NX) &
        -HeatDrain_col(NY,NX)-HeatDischar_col(NY,NX)-HeatCanopy2Dist_col(NY,NX)

      if(abs(WaterErr_test)>err_h2o)then
        write(110,*)('=',ii=1,50)
        write(110,*)I*1000+J,'NY,NX ',NY,NX
        write(110,*)'init H2O         =',WaterErr_col(NY,NX)
        write(110,*)'final H2O        =',WatMass_col(NY,NX)
        write(110,*)'delta H2O        =',WatMass_col(NY,NX)-WaterErr_col(NY,NX)
        write(110,*)'prec, irrig      =',PrecAtm_col(NY,NX),Irrigation_col(NY,NX) 
        write(110,*)'litfall H2O      =',RainLitr_col(NY,NX)
        write(110,*)'surf evap        =',VapXAir2GSurf_col(NY,NX)    
        write(110,*)'PrecIntceptByCan =',PrecIntceptByCanopy_col(NY,NX)   
        write(110,*)'plant evap       =',QVegET_col(NY,NX)
        write(110,*)'run on           =',QRunSurf_col(NY,NX)
        write(110,*)'Qinflx2Soil_col  =',Qinflx2Soil_col(NY,NX)
        write(110,*)'QWatIntLaterflow =',QWatIntLaterFlow_col(NY,NX)
        write(110,*)'drain            =',QDrain_col(NY,NX)
        write(110,*)'discharge        =',QDischar_col(NY,NX)
        write(110,*)'root uptake      =',TPlantRootH2OUptake_col(NY,NX)        
        write(110,*)'SnowMassBeg,end  =',SnowMassBeg_col(NY,NX),SnowMassEnd_col(NY,NX)
        write(110,*)'CanopyWatbeg,end =',CanopyWaterMassBeg_col(NY,NX),CanopyWaterMassEnd_col(NY,NX)
        write(110,*)'WatHeldOnCanopy  =',WatHeldOnCanopy_col(NY,NX)
        write(110,*)'CanopyWat_col    =',CanopyWat_col(NY,NX)
        write(110,*)'SoilWatMassBegEnd=',SoilWatMassBeg_col(NY,NX),SoilWatMassEnd_col(NY,NX),NU(NY,NX)
        write(110,*)'snofall          =',Prec2Snow_col(NY,NX),SnoFalPrec_col(NY,NX)
        write(110,*)'nsnol_col        =',nsnol_col(NY,NX)
        write(110,*)'snowloss         =',QSnowH2Oloss_col(NY,NX)
        write(110,*)'Snowxfer         =',QSnoWatXfer2Soil_col(NY,NX)+QSnoIceXfer2Soil_col(NY,NX)*DENSICE
        write(110,*)('-',ii=1,50)
        write(110,*)'precipErr_test   =',precipErr_test,RainPrecThrufall_col(NY,NX),RainPrec2Sno_col(NY,NX)
        write(110,*)'prec2expSErr_test=',prec2expSErr_test
        write(110,*)'prec2SnoErr_test =',prec2SnoErr_test        
        write(110,*)'canopyH2Oerr_test=',canopyH2Oerr_test
        write(110,*)'SnowMassErr_test =',SnowMassErr_test        
        write(110,*)'literH2Oerr_test =',literH2Oerr_test,LitWatMassBeg_col(NY,NX),LitWatMassEnd_col(NY,NX)
        write(110,*)'SoilWatErr_test  =',SoilWatErr_test
        if(abs(SoilWatErr_test)>1.e-4_r8) &
        call endrun('H2O error test failure in '//trim(mod_filename)//' at line',__LINE__)
      endif

      !Turn off the tracer mass conservation check momentarily, due to complication of 
      !grid change. It will be turned on when a safe strategy will be figured out later.
      
      DO idg=idg_beg,idg_NH3-1        
        tracer_mass_err = trcg_TotalMass_beg_col(idg,NY,NX)+SurfGasEmisFlx_col(idg,NY,NX)+GasHydroLossFlx_col(idg,NY,NX) &
          +RGasNetProd_col(idg,NY,NX)-trcg_TotalMass_col(idg,NY,NX)

        trcg_mass_cumerr_col(idg,NY,NX)=trcg_mass_cumerr_col(idg,NY,NX)+ tracer_mass_err 

        tracer_rootmass_err = trcg_rootMass_beg_col(idg,NY,NX)-trcg_rootMass_col(idg,NY,NX) &
          +trcg_air2root_flx_col(idg,NY,NX)-trcs_deadroot2soil_col(idg,NY,NX)+TRootGasLossDisturb_col(idg,NY,NX)

        if(idg==idg_O2)then
          tracer_rootmass_err = tracer_rootmass_err-RUptkRootO2_col(NY,NX)
        elseif(idg==idg_CO2)then
          tracer_rootmass_err = tracer_rootmass_err+RootCO2Emis2Root_col(NY,NX)
        else 
          tracer_rootmass_err = tracer_rootmass_err+trcs_plant_uptake_col(idg,NY,NX)
        endif
        
        tracer_snowmass_err=trcg_snowMass_beg_col(idg,NY,NX)-trcg_snowMass_col(idg,NY,NX) + &
          trcg_AquaAdv_flx_snvr(idg,1,NY,NX)-trcg_snowMassloss_col(idg,NY,NX)
        if(AMAX1(abs(tracer_mass_err),abs(tracer_rootmass_err))>1.e-5_r8)then
          write(111,*)('-',ii=1,50)
          write(111,*)I*1000+J,NY,NX,trcs_names(idg),iDayPlantHarvest_pft(1,NY,NX),iDayPlanting_pft(1,NY,NX)
          write(111,*)'beg end trc mass=',trcg_TotalMass_beg_col(idg,NY,NX),trcg_TotalMass_col(idg,NY,NX),&
            trcg_TotalMass_beg_col(idg,NY,NX)-trcg_TotalMass_col(idg,NY,NX)
          write(111,*)'mass_err         =',tracer_mass_err
          write(111,*)'surf emis        =',SurfGasEmisFlx_col(idg,NY,NX)
          write(111,*)'gasdif,ebu       =',GasDiff2Surf_flx_col(idg,NY,NX),trcg_ebu_flx_col(idg,NY,NX)
          write(111,*)'phenoflx         =',TRootGasLossDisturb_col(idg,NY,NX)
          write(111,*)'wedepo           =',Gas_WetDeposition_col(idg,NY,NX)          
          write(111,*)'---'
          write(111,*)'GasHydroloss     =',GasHydroLossFlx_col(idg,NY,NX)
          if(idg==idg_N2)then
            write(111,*)'RGasNetProd,rNfix,mup=',RGasNetProd_col(idg,NY,NX),RootN2Fix_col(NY,NX),trcs_RMicbUptake_col(idg,NY,NX)
          elseif(idg==idg_CO2)then
            write(111,*)'RGasNetProd     =',RGasNetProd_col(idg,NY,NX),-RootCO2AutorPrev_col(NY,NX)  
          else
            write(111,*)'RGasNetProd,micup =',RGasNetProd_col(idg,NY,NX),trcs_RMicbUptake_col(idg,NY,NX)
          endif
          
          write(111,*)'------------------'
          write(111,*)'snow beg_end mass=',trcg_snowMass_beg_col(idg,NY,NX),trcg_snowMass_col(idg,NY,NX)
          write(111,*)'snowloss2SoilLitr=',trcg_snowMassloss_col(idg,NY,NX)
          write(111,*)'tracersnowfall   =',trcg_AquaAdv_flx_snvr(idg,1,NY,NX)
          write(111,*)'snowmass err     =',tracer_snowmass_err
          write(111,*)'------------------'
          write(111,*)'soil beg_end mass=',trcg_soilMass_beg_col(idg,NY,NX), trcg_soilMass_col(idg,NY,NX),&
            trcg_soilMass_beg_col(idg,NY,NX)-trcg_soilMass_col(idg,NY,NX)
          write(111,*)'soil2atm         =',trcg_ebu_flx_col(idg,NY,NX)+GasDiff2Surf_flx_col(idg,NY,NX)+Gas_WetDeposition_col(idg,NY,NX)  
          
          if(idg==idg_CO2)then
            write(111,*)'ar2soil          =',RootCO2Ar2Soil_col(NY,NX)
          endif    
          write(111,*)'soil loss        =',trcg_ebu_flx_col(idg,NY,NX)+GasDiff2Surf_flx_col(idg,NY,NX)+Gas_WetDeposition_col(idg,NY,NX) &
              -trcs_plant_uptake_col(idg,NY,NX) 
          
          write(111,*)'=================='
          write(111,*)'root beg_end mass=',trcg_rootMass_beg_col(idg,NY,NX),trcg_rootMass_col(idg,NY,NX),&
            trcg_rootMass_beg_col(idg,NY,NX)-trcg_rootMass_col(idg,NY,NX)
          write(111,*)'root trcmass_err =',tracer_rootmass_err
          write(111,*)'pltair2root      =',trcg_air2root_flx_col(idg,NY,NX)     

          if(idg==idg_O2)then
            write(111,*)'plt_uptake       =',-RUptkRootO2_col(NY,NX)
          elseif(idg==idg_CO2)then
            write(111,*)'tplt2root,soi2root =',RootCO2Emis2Root_col(NY,NX),trcs_plant_uptake_col(idg_CO2,NY,NX)
            write(111,*)'ar2root            =',RootCO2Ar2Root_col(NY,NX)
          else
            write(111,*)'plt_soiluptake     =',trcs_plant_uptake_col(idg,NY,NX)
          endif
          write(111,*)'deadroot2soil    =',trcs_deadroot2soil_col(idg,NY,NX)

!          if(abs(tracer_mass_err)>1.e-1_r8) &
!            call endrun('tracer'//trcs_names(idg)//' error test failure in '//trim(mod_filename)//' at line',__LINE__)
        endif
      enddo      

      cycle      
      if(abs(HeatErr_test)>err_engy)then
        write(110,*)('-',ii=1,50)
        write(110,*)I+J/24.,'NY,NX,Engy conservation error',NY,NX,HeatErr_test
        write(110,*)'Init Engy=',HeatErr_col(NY,NX)
        write(110,*)'Final Engy=',HeatStore_col(NY,NX)
        write(110,*)'Delta Engy=',HeatStore_col(NY,NX)-HeatErr_col(NY,NX)
        write(110,*)'Delta snow engy=',SnowEngyEnd_col(NY,NX)-SnowEngyBeg_col(NY,NX)
        write(110,*)'Heat source=',HeatSource_col(NY,NX)
        write(110,*)'Net Rad=',Eco_NetRad_col(NY,NX)
        write(110,*)'Latent Heat=',Eco_Heat_Latent_col(NY,NX)
        write(110,*)'Sensible Heat=',Eco_Heat_Sens_col(NY,NX)
        write(110,*)'Heat prec=',PrecHeat_col(NY,NX)
        write(110,*)'Heat phase change, soil, snow=',THeatSoiThaw_col(NY,NX),THeatSnowThaw_col(NY,NX)
        write(110,*)'Heat runoff=',HeatRunSurf_col(NY,NX)
        write(110,*)'Heat drain=',HeatDrain_col(NY,NX)
        write(110,*)'Heat disch=',HeatDischar_col(NY,NX)
        call endrun('Engy error test failure in '//trim(mod_filename)//' at line',__LINE__)
      endif

    ENDDO
  ENDDO
  end subroutine EndCheckBalances 


!------------------------------------------------------------------------------------------
  subroutine SummarizeTracers(I,J,NHW,NHE,NVN,NVS)
  !
  !Description
  !sum up mass of tracers
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: NHW,NHE,NVN,NVS  
  integer :: NY,NX
  integer :: idg,L
  real(r8) :: trcg_snow(idg_beg:idg_NH3)
  real(r8) :: trcg_root(idg_beg:idg_NH3)
  real(r8) :: trcg_litr(idg_beg:idg_NH3)
  real(r8) :: trcg_soil(idg_beg:idg_end)

  DO  NX=NHW,NHE
    DO  NY=NVN,NVS
      trcg_TotalMass_col(idg_beg:idg_end,NY,NX)=0._r8
      trcg_snow=0._r8
      trcg_root=0._r8
      trcg_soil=0._r8
      trcg_litr=0._r8
      DO idg=idg_beg,idg_NH3
        trcg_litr(idg)=trcs_solml_vr(idg,0,NY,NX)
      enddo
      !note NL may change due to soil relayering, which can be tested using
      !by comparing the values of NL and NLI      
      DO L=NUI(NY,NX),NL(NY,NX)

        DO idg=idg_beg,idg_NH3
          trcg_soil(idg)=trcg_soil(idg) + trcg_gasml_vr(idg,L,NY,NX)
          trcg_root(idg)=trcg_root(idg) + trcg_root_vr(idg,L,NY,NX)                  
        ENDDO

        DO idg=idg_beg,idg_end  
          trcg_soil(idg)=trcg_soil(idg) + trcs_solml_vr(idg,L,NY,NX) + trcs_soHml_vr(idg,L,NY,NX)
        ENDDO
      ENDDO

      DO idg=idg_beg,idg_NH3  
        trcg_soilMass_col(idg,NY,NX)=trcg_soil(idg)+trcg_litr(idg)
      ENDDO    
      trcg_soilMass_col(idg_NH3B,NY,NX)=trcg_soil(idg_NH3B)

      !Because idg_NH3B does not exist in snow
      !add now      
      DO L=1,nsnol_col(NY,NX)          
        DO idg=idg_beg,idg_NH3
          trcg_snow(idg)=trcg_snow(idg)+trcg_solsml_snvr(idg,L,NY,NX)
        ENDDO
      ENDDO    
      DO idg=idg_beg,idg_NH3  
        trcg_snowMass_col(idg,NY,NX)=trcg_snow(idg)
      ENDDO    
      DO idg=idg_beg,idg_NH3
        trcg_rootMass_col(idg,NY,NX) = trcg_root(idg)
      ENDDO
      
      DO idg=idg_beg,idg_NH3
        trcg_TotalMass_col(idg,NY,NX)=trcg_snowMass_col(idg,NY,NX)+trcg_rootMass_col(idg,NY,NX)+trcg_soilMass_col(idg,NY,NX)
      ENDDO
      trcg_TotalMass_col(idg_NH3B,NY,NX)=trcg_soilMass_col(idg_NH3B,NY,NX)
      cycle
      if(I==2 .and. J>=5 .and. NY==1 .and. NX==5)then 
        idg=idg_CO2
        write(115,*)I*1000+J,'nsnol_col     trcname     snow   litr  root  soil'
        write(115,*)nsnol_col(NY,NX),trcs_names(idg),trcg_snow(idg),trcg_litr(idg),trcg_root(idg),trcg_soilMass_col(idg,NY,NX),NY,NX
        DO L=NUI(NY,NX),NL(NY,NX)
          write(115,*)L,trcg_gasml_vr(idg,L,NY,NX),trcs_solml_vr(idg,L,NY,NX), trcs_soHml_vr(idg,L,NY,NX)
        ENDDO  
      endif

    ENDDO
  ENDDO
  end subroutine SummarizeTracers

!------------------------------------------------------------------------------------------

  subroutine SummarizeTracerMass(I,J,NHW,NHE,NVN,NVS)

  implicit none 
  integer, intent(in) :: I,J,NHW,NHE,NVN,NVS
  integer :: NY,NX

  call SummarizeStorage(I,J,NHW,NHE,NVN,NVS)

  call SummarizeTracers(I,J,NHW,NHE,NVN,NVS)

  end subroutine SummarizeTracerMass

end module BalancesMod
