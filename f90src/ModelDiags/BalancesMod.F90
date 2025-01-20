module BalancesMod
  use data_kind_mod,  only: r8 => DAT_KIND_R8
  use CanopyDataType, only: QVegET_col
  use GridDataType, only : NU,NL
  use EcoSimConst, only : DENSICE
  use abortutils , only : endrun
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
  public :: SumUpStorage
contains

  subroutine BegCheckBalances(I,J,NHW,NHE,NVN,NVS)  
  !
  !Description
  !Prepare for next mass balance check
  implicit none
  integer, intent(in) :: I,J,NHW,NHE,NVN,NVS
  integer :: NY,NX
  
  DO  NX=NHW,NHE
    DO  NY=NVN,NVS
      WaterErr_col(NY,NX)           = WatMass_col(NY,NX)
      HeatErr_col(NY,NX)            = HeatStore_col(NY,NX)
      SnowEngyBeg_col(NY,NX)        = SnowEngyEnd_col(NY,NX)
      CanopyWaterMassBeg_col(NY,NX) = CanopyWaterMassEnd_col(NY,NX)
      SnowMassBeg_col(NY,NX)        = SnowMassEnd_col(NY,NX)
      LitWatMassBeg_col(NY,NX)      = LitWatMassEnd_col(NY,NX)
      SoilWatMassBeg_col(NY,NX)     = SoilWatMassEnd_col(NY,NX)
    ENDDO
  ENDDO  
  end subroutine BegCheckBalances 
!------------------------------------------------------------------------------------------
  subroutine SumUpStorage(I,J,NHW,NHE,NVN,NVS)  

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
  end subroutine SumUpStorage
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
  real(r8), parameter :: err_h2o=1.e-6_r8
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
  integer :: ii

  call SumUpStorage(I,J,NHW,NHE,NVN,NVS)

  DO  NX=NHW,NHE
    DO  NY=NVN,NVS

      SoilWatErr_test=SoilWatMassBeg_col(NY,NX)-SoilWatMassEnd_col(NY,NX)+Qinflx2Soil_col(NY,NX) &
        -QDrain_col(NY,NX)-QDischar_col(NY,NX)+TPlantRootH2OUptake_col(NY,NX)+QWatIntLaterFlow_col(NY,NX)

      precipErr_test = RainPrecThrufall_col(NY,NX)-Rain2LitR_col(NY,NX)-Rain2Soil_col(NY,NX)-RainPrec2Sno_col(NY,NX)
      prec2expSErr_test=Rain2ExposedSurf_col(NY,NX)-Rain2LitR_col(NY,NX)-Rain2Soil_col(NY,NX)
      prec2SnoErr_test= Prec2Snow_col(NY,NX)-RainPrec2Sno_col(NY,NX)-SnoFalPrec_col(NY,NX)  

      literH2Oerr_test = LitWatMassBeg_col(NY,NX)-LitWatMassEnd_col(NY,NX)+QRunSurf_col(NY,NX)+WatFLo2LitR_col(NY,NX)

      SnowMassErr_test = SnowMassBeg_col(NY,NX)-SnowMassEnd_col(NY,NX)+Prec2Snow_col(NY,NX)-QSnowH2Oloss_col(NY,NX)

      canopyH2Oerr_test=CanopyWaterMassBeg_col(NY,NX)-CanopyWaterMassEnd_col(NY,NX)+PrecIntceptByCanopy_col(NY,NX) &
        +QVegET_col(NY,NX)-TPlantRootH2OUptake_col(NY,NX)-QCanopyWat2Dist_col(NY,NX)

      WaterErr_test = WaterErr_col(NY,NX)-WatMass_col(NY,NX)+PrecAtm_col(NY,NX)+Irrigation_col(NY,NX)+QWatIntLaterFlow_col(NY,NX) &
        +RainLitr_col(NY,NX)+VapXAir2GSurf_col(NY,NX)+QVegET_col(NY,NX)+QRunSurf_col(NY,NX) &
        -QDrain_col(NY,NX)-QDischar_col(NY,NX)+TPlantRootH2OUptake_col(NY,NX)-QCanopyWat2Dist_col(NY,NX)

      HeatErr_test = HeatErr_col(NY,NX)-HeatStore_col(NY,NX)+THeatRootRelease_col(NY,NX) &
        +HeatSource_col(NY,NX)+Eco_NetRad_col(NY,NX)+Eco_Heat_Latent_col(NY,NX)+Eco_Heat_Sens_col(NY,NX)&
        +PrecHeat_col(NY,NX)+THeatSoiThaw_col(NY,NX)+THeatSnowThaw_col(NY,NX)+HeatRunSurf_col(NY,NX) &
        -HeatDrain_col(NY,NX)-HeatDischar_col(NY,NX)-HeatCanopy2Dist_col(NY,NX)

      if(abs(WaterErr_test)>err_h2o)then
        write(110,*)I+J/24.,'NY,NX, H2O conservation error',NY,NX,WaterErr_test,NU(NY,NX)
        write(110,*)'init H2O         =',WaterErr_col(NY,NX)
        write(110,*)'final H2O        =',WatMass_col(NY,NX)
        write(110,*)'delta H2O        =',WatMass_col(NY,NX)-WaterErr_col(NY,NX)
        write(110,*)'precipErr_test   =',precipErr_test,RainPrecThrufall_col(NY,NX),RainPrec2Sno_col(NY,NX)
        write(110,*)'prec2expSErr_test=',prec2expSErr_test
        write(110,*)'prec2SnoErr_test =',prec2SnoErr_test
        write(110,*)'prec, irrig      =',PrecAtm_col(NY,NX),Irrigation_col(NY,NX) 
        write(110,*)'litfall H2O      =',RainLitr_col(NY,NX)
        write(110,*)'surf evap        =',VapXAir2GSurf_col(NY,NX)    
        write(110,*)'PrecIntceptByCan =',PrecIntceptByCanopy_col(NY,NX)   
        write(110,*)'plant evap       =',QVegET_col(NY,NX)
        write(110,*)'root uptake      =',TPlantRootH2OUptake_col(NY,NX)
        write(110,*)'run on           =',QRunSurf_col(NY,NX)
        write(110,*)'drain            =',QDrain_col(NY,NX)
        write(110,*)'discharge        =',QDischar_col(NY,NX)
        write(110,*)'SnowMassBeg,end  =',SnowMassBeg_col(NY,NX),SnowMassEnd_col(NY,NX)
        write(110,*)'CanopyWatbeg,end =',CanopyWaterMassBeg_col(NY,NX),CanopyWaterMassEnd_col(NY,NX)
        write(110,*)'WatHeldOnCanopy  =',WatHeldOnCanopy_col(NY,NX)
        write(110,*)'CanopyWat_col    =',CanopyWat_col(NY,NX)
        write(110,*)'canopyH2Oerr_test=',canopyH2Oerr_test
        write(110,*)'SnowMassErr_test =',SnowMassErr_test        
        write(110,*)'literH2Oerr_test =',literH2Oerr_test,LitWatMassBeg_col(NY,NX),LitWatMassEnd_col(NY,NX)
        write(110,*)'SoilWatErr_test  =',SoilWatErr_test
        write(110,*)'SoilWatMassBegEnd=',SoilWatMassBeg_col(NY,NX),SoilWatMassEnd_col(NY,NX)
        write(110,*)'snofall          =',Prec2Snow_col(NY,NX),SnoFalPrec_col(NY,NX)
        write(110,*)'snowloss         =',QSnowH2Oloss_col(NY,NX)
        write(110,*)'Snowxfer         =',QSnoWatXfer2Soil_col(NY,NX)+QSnoIceXfer2Soil_col(NY,NX)*DENSICE
        write(110,*)('-',ii=1,50)
        if(abs(SoilWatErr_test)>1.e-4_r8) &
        call endrun('H2O error test failure in '//trim(mod_filename)//' at line',__LINE__)
      endif

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

end module BalancesMod
