module MethanotrophMod
  use data_kind_mod,        only: r8 => DAT_KIND_R8
  use MicFLuxTypeMod,       only: micfluxtype
  use MicStateTraitTypeMod, only: micsttype
  use MicForcTypeMod,       only: micforctype
  use EcoSiMParDataMod,     only: micpar
  use DebugToolMod,         only: PrintInfo
  use minimathmod,          only: AZMAX1, safe_adb
  use EcoSIMSolverPar,      only: dts_gas, NPH, NPT
  use EcoSimConst
  use NitroPars
  use MicrobeDiagTypes
  use MicrobMathFuncMod,    only: CalcRespMaint, StageAutotroph

  implicit none

  private
  public :: AeroMethanotrophCatabolism, AMONC10Catabolism, AMOANME2dCatabolism

  save
  character(len=*), parameter :: mod_filename = &
  __FILE__

  contains

!------------------------------------------------------------------------------------------

  subroutine AeroMethanotrophCatabolism(I,J,N,RMOMK,TOMEAutoKC,micfor,micstt,naqfdiag,nmicf,nmics,micflx,nmicdiag)
  !
  !CH4 uptake is paritioned into catabolic reaction RVOXP and 1st step anabolic reaction, where
  !the latter is partitioned into gross respiration (growth + maintenance respiration ) and biomass synthesis
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: N
  real(r8), intent(in) :: RMOMK(2)
  real(r8),intent(in) :: TOMEAutoKC
  type(micforctype), intent(in) :: micfor
  type(micsttype), intent(in) :: micstt
  type(Cumlate_Flux_Diag_type), intent(in) :: naqfdiag
  type(Microbe_Flux_type), intent(inout) :: nmicf
  type(Microbe_State_type), intent(inout):: nmics
  type(micfluxtype), intent(inout) :: micflx
  type(Microbe_Diag_type), intent(inout) :: nmicdiag
  character(len=*), parameter :: subname='AeroMethanotrophCatabolism'
  integer  :: M,MM,NGL
  real(r8)  :: RGOMP     !methane oxidized into CH2O
  real(r8)  :: RVOXP
  real(r8)  :: RVOXPA    !methane oxidation
  real(r8)  :: RVOXPB
  real(r8) :: CH4G1,CH4S1,CCH4S1
  real(r8) :: RCH4L1,RCH4F1,RCH4S1
  real(r8) :: RGOMP1,RCHDF
  real(r8) :: VMAX1,VOLWCH
  real(r8) :: RVOXP1
  real(r8) :: VMAX,FCH4X
  real(r8) :: pscal
  REAL(R8) :: VOLWPM

  associate(                                                &
    GrowthEnvScalAutor     => nmics%GrowthEnvScalAutor,     &
    FBiomNutStoiScalAutor  => nmics%FBiomNutStoiScalAutor,  &
    FSBSTAutor             => nmicdiag%FSBSTAutor,          &
    OMActAutor             => nmics%OMActAutor,             &
    FracOMActAutor         => nmics%FracOMActAutor,         &
    RO2Dmnd4GrossRespAutor => nmicf%RO2Dmnd4GrossRespAutor, &
    RespGrossAutor         => nmicf%RespGrossAutor,         &
    ECHZAutor              => nmicf%ECHZAutor,              &
    RSMetaOxidSoilAutor    => nmicf%RSMetaOxidSoilAutor,    &
    RSMetaOxidBandAutor    => nmicf%RSMetaOxidBandAutor,    &
    CCH4E                  => micfor%CCH4E,                 &
    VLsoiAirPM             => micfor%VLsoiAirPM,            &
    VLWatMicPM             => micfor%VLWatMicPM,            &
    ZEROS2                 => micfor%ZEROS2,                &
    ZEROS                  => micfor%ZEROS,                 &
    THETPM                 => micfor%THETPM,                &
    DiffusivitySolutEff    => micfor%DiffusivitySolutEff,   &
    litrm                  => micfor%litrm,                 &
    JGniA                  => micpar%JGniA,                 &
    JGnfA                  => micpar%JGnfA,                 &
    CCH4G                  => micstt%CCH4G,                 &
    CH4S                   => micstt%CH4S,                  &
    RCH4MetaDmndAutor      => micflx%RCH4MetaDmndAutor,     &
    RCH4MetaDmndAutorPrev  => micflx%RCH4MetaDmndAutorPrev, &
    CH4AquaSolubility      => micstt%CH4AquaSolubility,     &
    RCH4EcoDmndPrev        => micfor%RCH4EcoDmndPrev,       &
    RCH4PhysexchPrev       => micfor%RCH4PhysexchPrev,      &
    RCH4GasXchangePrev     => micfor%RCH4GasXchangePrev,    &
    RO2MetaDmndAutor       => micflx%RO2MetaDmndAutor       &
  )
!     begin_execution
!
!     CH4 OXIDATION FROM SPECIFIC OXIDATION RATE, ENERGY YIELD,
!     ACTIVE OXIDIZER BIOMASS, TEMPERATURE, AQUEOUS CO2 AND
!     CH4 CONCENTRATIONS IN BAND AND NON-BAND SOIL ZONES
!
!     ECHZ=growth respiration efficiency
!     VMAX=potential oxidation
!     RCH4PhysexchPrev=total aqueous CH4 exchange from previous hour
!     RCH4GasXchangePrev=total gaseous CH4 exchange from previous hour
!     tCH4ProdAceto+tCH4ProdH2=total CH4 generated from methanogenesis
!     dts_gas=1.0/(NPH*NPT)
!     CH4G1,CH4S1=CH4 gaseous, aqueous amounts
!     CCH4E,CCH4G=CH4 gas concentration in atmosphere, soil
!     VLsoiAirPM,VLWatMicPM=air,water-filled porosity
!     CH4AquaSolubility=CH4 aqueous solubility
!     CCK4=Km for CH4 uptake
!     ECHO=efficiency CO2 conversion to biomass
!     RGOMP1=substrate-limited CH4 oxidation
!     RCHDF=gaseous-aqueous CH4 exchange
!     DiffusivitySolutEff=rate constant for gaseous-aqueous exchange
!
  call PrintInfo('beg '//subname)
  DO NGL=micpar%JGniA(N),micpar%JGnfA(N)
    IF(OMActAutor(NGL).LE.0.0_r8)cycle
    call StageAutotroph(NGL,N,TOMEAutoKC,micfor,nmics,nmicdiag)

    call CalcRespMaint(I,J,NGL,RMOMK,micfor,micstt,micflx,nmicf,nmics)

    IF(RCH4EcoDmndPrev.GT.ZEROS)THEN
      FCH4X=AMAX1(FMN,RCH4MetaDmndAutorPrev(NGL)/RCH4EcoDmndPrev)
    ELSE
      FCH4X=AMAX1(FMN,FracOMActAutor(NGL))
    ENDIF

    ECHZAutor(NGL) = EH4X
    VMAX           = GrowthEnvScalAutor(NGL)*FBiomNutStoiScalAutor(NGL)*OMActAutor(NGL)*VMXCH4OxiAero
    RCH4L1         = RCH4PhysexchPrev*dts_gas*FCH4X
    RCH4F1         = RCH4GasXchangePrev*dts_gas*FCH4X
    RCH4S1         = (naqfdiag%tCH4ProdAceto+naqfdiag%tCH4ProdH2)*dts_gas*FCH4X

    IF(litrm)THEN
      !surface residue layer
      CH4G1 = CCH4E*VLsoiAirPM(1)*FCH4X
      VMAX1 = AZMAX1(AMIN1(VMAX*dts_gas,CH4G1))  !apparent vmax for uptake
    ELSE
      CH4G1 = CCH4G*VLsoiAirPM(1)*FCH4X
      VMAX1 = VMAX*dts_gas  !apparent vmax for uptake
    ENDIF
    CH4S1 = CH4S*FCH4X
    RVOXP = 0.0_r8
    RGOMP = 0.0_r8
    !
    !     CH4 DISSOLUTION FROM GASEOUS PHASE SOLVED IN SHORTER TIME STEP
    !     TO MAINTAIN AQUEOUS CH4 CONCENTRATION DURING OXIDATION
    ! for aerobic methanotrophs, CH4 is oxiized to CO2 for energy, and also
    ! to intracellular C for respiration. (Grant 1999), the C yield is approximated
    ! as the energy required for turning CH4 into organic C, and the energy released
    ! from turnining CH4 into CO2.
    ! the catabolic reaction is (molar-basis)
    !  CH4 + 2O2 -> CO2 + 2H2O
    !  mass basis becomes (2*32/12=5.33)
    !  CH4 + 5.33 O2 -> CO2 + 2H2O, RVOXP
    !
    !  CH2O+ O2 -> CO2,  32/12=2.667

    D320: DO M=1,NPH
      IF(VLWatMicPM(M).GT.ZEROS2)THEN
        VOLWCH = VLWatMicPM(M)*CH4AquaSolubility
        VOLWPM = VOLWCH+VLsoiAirPM(M)

        !CH4 uptake by aerobic oxidation
        !RCH4F1: net gaseous CH4 flux from transport into the layer
        !RCH4L1: net aqueous CH4 flux from transport into the layer
        !RCH4S1: aquoues CH4 production from methanogenesis in the layer
        D325: DO MM=1,NPT
          CH4G1           = AZMAX1(CH4G1+RCH4F1)
          CH4S1           = AZMAX1(CH4S1+RCH4L1+RCH4S1)
          CCH4S1          = safe_adb(CH4S1,VLWatMicPM(M))
          FSBSTAutor(NGL) = CCH4S1/(CCH4S1+CCK4)

          !CH4->CH2O->CO2
          !RVOXP1 CH4 used to produce energy
          !ECHO: carbon yield,
          !ECHZAutor(NGL): respiration ratio
          !CH4 is oxidized to generate energy RVOXP1, which supports the production of RGOMP1
          RVOXP1=AMIN1(CH4S1*0.9999_r8/(1.0_r8+ECHO*ECHZAutor(NGL)),VMAX1*FSBSTAutor(NGL))

          !the respiration yield of CH2O
          RGOMP1 = RVOXP1*ECHO*ECHZAutor(NGL)
          CH4S1  = CH4S1-RVOXP1-RGOMP1
          !dissolution-vaporization
          IF(THETPM(M).GT.AirFillPore_Min)THEN
            RCHDF=DiffusivitySolutEff(M)*(AMAX1(ZEROS,CH4G1)*VOLWCH-CH4S1*VLsoiAirPM(M))/VOLWPM
            RCHDF=AMAX1(AMIN1(CH4G1,RCHDF),-CH4S1)
          ELSE
            RCHDF=0.0_r8
          ENDIF
          CH4G1 = CH4G1-RCHDF
          CH4S1 = CH4S1+RCHDF
          RVOXP = RVOXP+RVOXP1
          RGOMP = RGOMP+RGOMP1
        ENDDO D325
      ENDIF
    ENDDO D320

    RVOXPA = AZMAX1(RVOXP) !CH4 oxidized to support production for RGOMP, which is used to compute growth + (growth/maint resp)
    RVOXPB = 0.0_r8
    !
    !     O2 DEMAND FROM CH4 OXIDATION
    ! CH4 taken up is partitioned into RVOXP and RGOMP/ECHZAutor(NGL)
    ! RGOMP/ECHZAutor(NGL) is partitioned into growth + maintenance respiraiton + growth
    ! with growth = growth respiration*(1/ECHZAutor(NGL)-1).
    ! thus the uptake for growth is growth respiration/ECHZAutor(NGL),
    ! so the total upatke is RVOXP + maintenance + uptake for growth
    ! note maintenance + growth_resp/ECHZAutor(NGL) = gross_resp/ECHZAutor(NGL)+(1-/ECHZAutor(NGL))*maintenance
    ! < gross_resp/ECHZAutor(NGL), meaning the model does not have exact stoichiometry balance.
    !     RO2Dmnd4RespHeter=O2 demand from respiration
    !     ROXYP=O2 demand from respiration + CH4 oxidation
    !CH4+O2  -> CH2O + H2O, RGOXP*ECHO
    !CH2O+O2 -> CO2 + 2H2O, RGOMP as CO2
    RO2Dmnd4GrossRespAutor(NGL) = 2.667_r8*RGOMP
    RO2MetaDmndAutor(NGL)       = RO2Dmnd4GrossRespAutor(NGL)+5.333_r8*RVOXP
    RCH4MetaDmndAutor(NGL)      = RVOXPA
    RespGrossAutor(NGL)         = RGOMP
    RSMetaOxidSoilAutor(NGL)    = RVOXPA
    RSMetaOxidBandAutor(NGL)    = 0._r8
  ENDDO
  call PrintInfo('end '//subname)
  end associate
  end subroutine AeroMethanotrophCatabolism

!------------------------------------------------------------------------------------------
  subroutine AMONC10Catabolism(I,J,N,RMOMK,TOMEAutoKC,micfor,micstt,naqfdiag,nmicf,nmics,micflx,nmicdiag)
  !
  !Description:
  !  Catabolic reaction:
  !8/3NO2(-)+8/3H(+)+CH4(aq) -> CO2(aq)+4/3N2(aq)+10/3H2O, N/C=14*8/(3*12)=3.111
  !Gibbs free energy: -1050 kJ (molCH4)-1, -87.5 kJ (molC)-1

  !Wei et al. (2022) The denitrifying anaerobic methane oxidation process and microorganisms in the environments: A review
  !when fixing carbon for biomass, it takes CO2 from the environment

  implicit none
  integer,  intent(in)  :: I,J,N
  real(r8), intent(in) :: RMOMK(2)
  real(r8), intent(in)  :: TOMEAutoKC
  type(micforctype), intent(in) :: micfor
  type(micsttype), intent(inout) :: micstt
  type(micfluxtype), intent(inout) :: micflx
  type(Cumlate_Flux_Diag_type),INTENT(INOUT) :: naqfdiag
  type(Microbe_State_type), intent(inout) :: nmics
  type(Microbe_Flux_type), intent(inout) :: nmicf
  type(Microbe_Diag_type), intent(inout) :: nmicdiag
  character(len=*), parameter :: subname='AMONC10Catabolism'
  real(r8) :: FNO2S,FNO2B,VMX2S,VMX2B,XCH4,FNO2,FNB2,XCO2,FCH4X
  real(r8) :: RNNO2,RNNOB,VMAX,RGOMP,GCH4X,GCH4O,RVOXP,SCAL,OXYI
  real(r8) :: GCNCX= 87.5_r8
  real(r8) :: GANCX= 13.17_r8
  real(r8) :: GRNCX  = 50.33_r8  !free energy yields of redox reaction for respiration, [kJ gC-1]
  real(r8) :: ENC10,ECH0,RNDTP              !growth respiraiton efficiency for AMO NC10, [-]
  integer :: NGL

  associate(                                              &
   GrowthEnvScalAutor    => nmics%GrowthEnvScalAutor,     &
   FBiomNutStoiScalAutor => nmics%FBiomNutStoiScalAutor,  &
   FracOMActAutor        => nmics%FracOMActAutor,         &
   RNO2EcoUptkSoilPrev   => micfor%RNO2EcoUptkSoilPrev,   &
   RNO2XupAutorPrev      => micflx%RNO2XupAutorPrev,      &
   RNO2XupAutorBandPrev  => micflx%RNO2XupAutorBandPrev,  &
   RNO2EcoUptkBandPrev   => micfor%RNO2EcoUptkBandPrev,   &
   RNO2XupAutor          => micflx%RNO2XupAutor,          &
   RNO2XupAutorBand      => micflx%RNO2XupAutorBand,      &
   FracNO2XupAutor       => nmics%FracNO2XupAutor,        &
   FSBSTAutor            => nmicdiag%FSBSTAutor,          &
   OMActAutor            => nmics%OMActAutor,             &
   RCH4MetaDmndAutor     => micflx%RCH4MetaDmndAutor,     &
   RCH4MetaDmndAutorPrev => micflx%RCH4MetaDmndAutorPrev, &
   RCH4EcoDmndPrev       => micfor%RCH4EcoDmndPrev,       &
   RNOxReduxAutorSoil    => nmicf%RNOxReduxAutorSoil,     & !NO2(-) reduction by N2
   RNOxReduxAutorBand    => nmicf%RNOxReduxAutorBand,     & !NO2(-) reduction by N2
   RSMetaOxidSoilAutor   => nmicf%RSMetaOxidSoilAutor,    &
   COXYS                 => micstt%COXYS,                 &
   ZEROS                 => micfor%ZEROS,                 &
   CCO2S                 => micstt%CCO2S,                 &
   RespGrossAutor        => nmicf%RespGrossAutor,         &
   RCO2ProdAutor         => nmicf%RCO2ProdAutor,          &
   ECHZAutor             => nmicf%ECHZAutor,              &
   VLNO3                 => micfor%VLNO3,                 &
   VLNOB                 => micfor%VLNOB,                 &
   ZERO                  => micfor%ZERO,                  &
   JGniA                 => micpar%JGniA,                 &
   JGnfA                 => micpar%JGnfA,                 &
   TKS                   => micfor%TKS,                   &
   CCH4S                 => micstt%CCH4S,                 &
   ZNO2B                 => micstt%ZNO2B,                 &
   ZNO2S                 => micstt%ZNO2S,                 &
   CNO2B                 => micstt%CNO2B,                 &
   CNO2S                 => micstt%CNO2S,                 &
   CH4S                  => micstt%CH4S                   &
  )
  call PrintInfo('beg '//subname)
  FNO2S = VLNO3
  FNO2B = VLNOB
  ENC10 = GRNCX/EOMH
  OXYI  = 1.0_r8-1.0_r8/(1.0_r8+EXP(1.0_r8*AMAX1(-COXYS+2.5_r8,-50._r8)))
  XCO2  = CCO2S/(CCO2S+CCKM)*OXYI

  DO NGL  = micpar%JGniA(N), micpar%JGnfA(N)
    IF(OMActAutor(NGL).LE.0.0_r8)cycle
    call StageAutotroph(NGL,N,TOMEAutoKC,micfor,nmics,nmicdiag)

    call CalcRespMaint(I,J,NGL,RMOMK,micfor,micstt,micflx,nmicf,nmics)

    IF(RCH4EcoDmndPrev.GT.ZEROS)THEN
      FCH4X=AMAX1(FMN,RCH4MetaDmndAutorPrev(NGL)/RCH4EcoDmndPrev)
    ELSE
      FCH4X=AMAX1(FMN,FracOMActAutor(NGL))
    ENDIF

    XCH4=CCH4S/(CCH4S+CCK4)
    IF(RNO2EcoUptkSoilPrev.GT.ZEROS)THEN
      FNO2=AMAX1(FMN,RNO2XupAutorPrev(NGL)/RNO2EcoUptkSoilPrev)
    ELSE
      FNO2=AMAX1(FMN,FracNO2XupAutor(NGL)*VLNO3)
    ENDIF

    IF(RNO2EcoUptkBandPrev.GT.ZEROS)THEN
      FNB2=AMAX1(FMN,RNO2XupAutorBandPrev(NGL)/RNO2EcoUptkBandPrev)
    ELSE
      FNB2=AMAX1(FMN,FracNO2XupAutor(NGL)*VLNOB)
    ENDIF
    naqfdiag%TFNO2X=naqfdiag%TFNO2X+FNO2
    naqfdiag%TFNO2B=naqfdiag%TFNO2B+FNB2

    VMAX=GrowthEnvScalAutor(NGL)*FBiomNutStoiScalAutor(NGL)*XCH4*OMActAutor(NGL)*VMX2AMONC10*XCO2
    IF(CNO2S.GT.ZERO)THEN
      VMX2S=VMAX*FNO2S*CNO2S/(CNO2S+Z2KM)
    ELSE
      VMX2S=0.0_r8
    ENDIF

    !band-soil
    IF(CNO2B.GT.ZERO)THEN
      VMX2B=VMAX*FNO2B*CNO2B/(CNO2B+Z3KM)
    ELSE
      VMX2B=0.0_r8
    ENDIF

    RNNO2          = AZMAX1(AMIN1(VMX2S,FNO2*ZNO2S))
    RNNOB          = AZMAX1(AMIN1(VMX2B,FNB2*ZNO2B))
    GCH4X          = RGASC*1.E-3_r8*TKS*LOG((AMAX1(1.0E-08_r8,CCH4S)/12._r8))/12._r8
    ECHZAutor(NGL) = AMAX1(EO2X,AMIN1(1.0_r8,1.0_r8/(1.0_r8+2._r8*AZMAX1((GCNCX+GCH4X))/EOMH)))
    RNDTP= RNNOB+RNNO2

    IF(RNDTP>3.111_r8*CH4S*FCH4X)THEN
      scal = 3.111_r8*CH4S*FCH4X/RNDTP
      RNNOB=RNNOB*SCAL
      RNNO2=RNNO2*SCAL
      RNDTP= RNNOB+RNNO2
    ENDIF
    RGOMP=RNDTP/3.111_r8
    RespGrossAutor(NGL)      = RGOMP
    RSMetaOxidSoilAutor(NGL) = RGOMP !CH4 oxidized
    RCH4MetaDmndAutor(NGL)   = RGOMP !CH4 demanded

    RNO2XupAutor(NGL)        = RNNO2
    RNO2XupAutorBand(NGL)    = RNNOB
    RNOxReduxAutorSoil(NGL)  = RNNO2
    RNOxReduxAutorBand(NGL)  = RNNOB
    RCO2ProdAutor(NGL)       = RGOMP  !BIOMASS C will be CO2 from the environment
  ENDDO
  call PrintInfo('end '//subname)
  end associate

  end subroutine AMONC10Catabolism
!------------------------------------------------------------------------------------------
  subroutine AMOANME2dCatabolism(I,J,N,RMOMK,TOMEAutoKC,VOLWZ,micfor,micstt,naqfdiag,nmicf,nmics,micflx,nmicdiag)
  !
  !Description:
  !ANME-2d uses CH4 for both energy and carbon biomass, specifically, CH4 is first oxidized to CO2 to produce energy,
  !then some CO2 is re-assimilated for biomass.
  !
  !  Catabolic reaction to produce respiration:
  !4NO3(-)+CH4(aq)->4NO2(-)+HCO3(-)+H2O+H(+), 4*14/12=4.667
  !Gibbs free energy: -510.0kJ/molC = 41.7 kJ/gC

  !Ref: Bhattarai et al. (2019), Physiology and Distribution of Archaeal Methanotrophs That Couple Anaerobic Oxidation of Methane with Sulfate Reduction

  implicit none
  integer,  intent(in)  :: I,J,N
  real(r8), intent(in) :: RMOMK(2)
  REAL(R8), INTENT(IN)  :: VOLWZ
  real(r8), intent(in)  :: TOMEAutoKC
  type(micforctype), intent(in) :: micfor
  type(micsttype), intent(inout) :: micstt
  type(micfluxtype), intent(inout) :: micflx
  type(Cumlate_Flux_Diag_type),INTENT(INOUT) :: naqfdiag
  type(Microbe_State_type), intent(inout) :: nmics
  type(Microbe_Flux_type), intent(inout) :: nmicf
  type(Microbe_Diag_type), intent(inout) :: nmicdiag
  character(len=*), parameter :: subname='AMOANME2dCatabolism'
  real(r8) :: FNO3S,FNO3B,VMXDXS,VMXDXB,VMXDXT,FVMXDX,XCH4,scal
  real(r8) :: XCO2,FCH4X,OXYI
  real(r8) :: RGOMP,FNB3X,FNO3X,VMAX,RNNO3,RNN3B,RVOXP,GCH4X,GCH4O
  integer :: NGL

  associate(                                              &
   GrowthEnvScalAutor    => nmics%GrowthEnvScalAutor,     &
   FBiomNutStoiScalAutor => nmics%FBiomNutStoiScalAutor,  &
   FSBSTAutor            => nmicdiag%FSBSTAutor,          &
   FracOMActAutor        => nmics%FracOMActAutor,         &
   OMActAutor            => nmics%OMActAutor,             &
   RespGrossAutor        => nmicf%RespGrossAutor,         &
   RCO2ProdAutor         => nmicf%RCO2ProdAutor,          &
   RNO3XupAutor          => micflx%RNO3XupAutor    ,      &
   RNO3XupAutorBand      => micflx%RNO3XupAutorBand,      &
   RNO3XupAutorPrev      => micflx%RNO3XupAutorPrev,      &
   RNO3XupAutorBandPrev  => micflx%RNO3XupAutorBandPrev,  &
   RNO3EcoDmndSoilPrev   => micfor%RNO3EcoDmndSoilPrev,   &
   RNO3EcoDmndBandPrev   => micfor%RNO3EcoDmndBandPrev,   &
   RSMetaOxidSoilAutor   => nmicf%RSMetaOxidSoilAutor,    &
   RCH4MetaDmndAutor     => micflx%RCH4MetaDmndAutor,     &
   RCH4MetaDmndAutorPrev => micflx%RCH4MetaDmndAutorPrev, &
   RCH4EcoDmndPrev       => micfor%RCH4EcoDmndPrev,       &
   RNOxReduxAutorSoil    => nmicf%RNOxReduxAutorSoil,     & !NO3(-) reduction by NO2(-)
   RNOxReduxAutorBand    => nmicf%RNOxReduxAutorBand,     & !NO3(-) reduction by NO2(-)
   ECHZAutor             => nmicf%ECHZAutor,              & !respiraiton efficiency
   COXYS                 => micstt%COXYS,                 &
   CCO2S                 => micstt%CCO2S,                 &
   JGniA                 => micpar%JGniA,                 &
   JGnfA                 => micpar%JGnfA,                 &
   TKS                   => micfor%TKS,                   &
   ZEROS                 => micfor%ZEROS,                 &
   ZEROS2                => micfor%ZEROS2,                &
   VLNO3                 => micfor%VLNO3,                 &
   VLNOB                 => micfor%VLNOB,                 &
   ZERO                  => micfor%ZERO,                  &
   CCH4S                 => micstt%CCH4S,                 &
   CH4S                  => micstt%CH4S,                  &
   CNO3B                 => micstt%CNO3B,                 &
   CNO3S                 => micstt%CNO3S,                 &
   ZNO3B                 => micstt%ZNO3B,                 &
   ZNO3S                 => micstt%ZNO3S                  &
  )
  call PrintInfo('beg '//subname)
  FNO3S = VLNO3
  FNO3B = VLNOB
  OXYI  = 1.0_r8-1.0_r8/(1.0_r8+EXP(1.0_r8*AMAX1(-COXYS+2.5_r8,-50._r8)))
  DO NGL  = micpar%JGniA(N), micpar%JGnfA(N)
    IF(OMActAutor(NGL).LE.0.0_r8)cycle
    call StageAutotroph(NGL,N,TOMEAutoKC,micfor,nmics,nmicdiag)

    call CalcRespMaint(I,J,NGL,RMOMK,micfor,micstt,micflx,nmicf,nmics)

    IF(RCH4EcoDmndPrev.GT.ZEROS)THEN
      FCH4X=AMAX1(FMN,RCH4MetaDmndAutorPrev(NGL)/RCH4EcoDmndPrev)
    ELSE
      FCH4X=AMAX1(FMN,FracOMActAutor(NGL))
    ENDIF

    IF(RNO3EcoDmndSoilPrev.GT.ZEROS)THEN
      FNO3X=AMAX1(FMN,RNO3XupAutorPrev(NGL)/RNO3EcoDmndSoilPrev)
    ELSE
      FNO3X=AMAX1(FMN,FracOMActAutor(NGL)*VLNO3)
    ENDIF

    IF(RNO3EcoDmndBandPrev.GT.ZEROS)THEN
      FNB3X=AMAX1(FMN,RNO3XupAutorBandPrev(NGL)/RNO3EcoDmndBandPrev)
    ELSE
      FNB3X=AMAX1(FMN,FracOMActAutor(NGL)*VLNOB)
    ENDIF

    naqfdiag%TFNO3X = naqfdiag%TFNO3X+FNO3X
    naqfdiag%TFNO3B = naqfdiag%TFNO3B+FNB3X

    XCH4=CCH4S/(CCH4S+CCK4)
    VMAX=GrowthEnvScalAutor(NGL)*FBiomNutStoiScalAutor(NGL)*XCH4*OMActAutor(NGL)*VMX3AMO2D*OXYI
    IF(CNO3S.GT.ZERO)THEN
      VMXDXS=VMAX*FNO3S*CNO3S/(CNO3S+Z3KM)
    ELSE
      VMXDXS=0.0_r8
    ENDIF
    !band-soil
    IF(CNO3B.GT.ZERO)THEN
      VMXDXB=VMAX*FNO3B*CNO3B/(CNO3B+Z3KM)
    ELSE
      VMXDXB=0.0_r8
    ENDIF
    VMXDXT=VMXDXS+VMXDXB

    !product inhibition
    IF(VOLWZ.GT.ZEROS2)THEN
      FVMXDX=1.0_r8/(1.0_r8+VMXDXT/(VMKI*VOLWZ))
    ELSE
      FVMXDX=0.0_r8
    ENDIF
    VMXDXS = VMXDXS*FVMXDX
    VMXDXB = VMXDXB*FVMXDX

    RNNO3               = AZMAX1(AMIN1(VMXDXS,FNO3X*ZNO3S))
    RNN3B               = AZMAX1(AMIN1(VMXDXB,FNB3X*ZNO3B))
    RVOXP               = RNNO3+RNN3B !partitioned between catabolic reaction and anabolic reaction
    GCH4X               = RGASC*1.E-3_r8*TKS*LOG((AMAX1(1.0E-08_r8,CCH4S)/12._r8))
    GCH4O               = GCH4X/12._r8
    ECHZAutor(NGL)      = AMAX1(EO2X,AMIN1(1.0_r8,1.0_r8/(1.0_r8+2._r8*AZMAX1((GN3CX+GCH4O))/EOMH)))

    !follow the approach by Hydrogenotrophic methanogen
    !catabolic energy measured by C, 4NO3(-)+CH4->4NO2(-)+HCO3(-)+H2O+H(+), 4*14/12=4.667
    !RGOMP is the C-eqv gross respiraiton, or CH4 oxidized to derive energy for anabolic reaction
    !assuming electrons/reducing power is from the catabolic reaction, so no more NO3(-) is used in converting CH4 to C biomass.
    !
    IF(RVOXP/4.667_r8>CH4S*FCH4X)THEN
      SCAL  = CH4S*FCH4X*4.667_r8/RVOXP
      RGOMP = CH4S*FCH4X
      RNNO3 = RNNO3*scal
      RNN3B = RNN3B*scal
    ELSE
      RGOMP                  = RVOXP/4.667_r8    !total energy generated from NO3(-) reduction,
    ENDIF
    RespGrossAutor(NGL)      = RGOMP
    RSMetaOxidSoilAutor(NGL) = RGOMP !total CH4 oxidized
    RCH4MetaDmndAutor(NGL)   = RGOMP

    RNO3XupAutor(NGL)        = RNNO3
    RNO3XupAutorBand(NGL)    = RNN3B
    RNOxReduxAutorSoil(NGL)  = RNNO3
    RNOxReduxAutorBand(NGL)  = RNN3B
    RCO2ProdAutor(NGL)       = RGOMP    !CO2 produced from catabolic reaction, some of it will be reassimilated for biomass
  ENDDO
  call PrintInfo('end '//subname)
  end associate
  end subroutine AMOANME2dCatabolism
end module MethanotrophMod
