module PlantNonstElmDynMod
  use minimathmod  , only : safe_adb,AZMAX1,AZMIN1
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use ElmIDMod
  use PlantBGCPars
  use PlantAPIData  
  implicit none

  private

! end_include_section

  character(len=*), private, parameter :: mod_filename = &
  __FILE__

  public :: PlantNonstElmTransfer
  public :: SeasonStoreShootTransfer
  public :: StalkRsrvShootNonstTransfer
  public :: StalkRsrvRootNonstTransfer
  contains
  ![header]
!----------------------------------------------------------------------------------------------------
  subroutine WithinBranchElmTransfer(I,J,NZ)

  implicit none
  integer, intent(in) :: I,J,NZ
  real(r8) :: TotNonstElm_loc(NumPlantChemElms)
  real(r8) :: NonstElm_loc(NumPlantChemElms,NumCanopyLayers1)  
  real(r8) :: TwoCompMassC,NonstElmGradt
  real(r8) :: StalkRsrvGradt  
  real(r8) :: TotStalkRsrv_loc(NumPlantChemElms)  
  real(r8) :: TotStalkMassC  
  real(r8) :: sumchk1,sumchk2  
  integer  :: NB,NE
  real(r8) :: XFRE(NumPlantChemElms)
  real(r8) :: mass_inital(NumPlantChemElms)
  real(r8) :: mass_finale(NumPlantChemElms)
  real(r8) :: LeafPetoMassC_brch(NumCanopyLayers1)
  associate(                                                       &
    iPlantPhenolPattern_pft => plt_pheno%iPlantPhenolPattern_pft  ,& !input  :plant growth habit: annual or perennial,[-]
    iPlantCalendar_brch     => plt_pheno%iPlantCalendar_brch      ,& !input  :plant growth stage, [-]
    Hours2LeafOut_brch      => plt_pheno%Hours2LeafOut_brch       ,& !input  :counter for mobilizing nonstructural C during spring leafout/dehardening, [h]
    SapwoodBiomassC_brch  => plt_biom%SapwoodBiomassC_brch    ,& !input  :branch live stalk C, [gC d-2]
    LeafPetolBiomassC_brch  => plt_biom%LeafPetolBiomassC_brch    ,& !input  :plant branch leaf + sheath C, [g d-2]
    ZERO4Groth_pft          => plt_biom%ZERO4Groth_pft            ,& !input  :threshold zero for plang growth calculation, [-]
    iPlantBranchState_brch  => plt_pheno%iPlantBranchState_brch   ,& !input  :flag to detect branch death, [-]
    NumOfBranches_pft       => plt_morph%NumOfBranches_pft        ,& !input  :number of branches,[-]
    StalkRsrvElms_brch      => plt_biom%StalkRsrvElms_brch        ,& !inoput :branch reserve element mass, [g d-2]
    CanopyNonstElms_brch    => plt_biom%CanopyNonstElms_brch       & !inoput :branch nonstructural element, [g d-2]
  )
    !     TRANSFER NON-STRUCTURAL C,N,P AMONG BRANCH LEAVES
    !     FROM NON-STRUCTURAL C,N,P CONCENTRATION DIFFERENCES
    !     WHEN SEASONAL STORAGE C IS NOT BEING MOBILIZED
    !
    !     iPlantBranchState_brch=branch living flag: 0=alive,1=dead
    !     Hours2LeafOut_brch=hourly leafout counter
    !     HourReq2InitSStor4LeafOut=number of hours required to initiate remobilization of storage C for leafout
    !     LeafPetolBiomassC_brch=leaf+petiole mass
    !     CPOOL,ZPOOL,PPOOL=non-structural C,N,P mass in branch
    !
  TwoCompMassC=0._r8
  TotNonstElm_loc(1:NumPlantChemElms)=0._r8
  D300: DO NB=1,NumOfBranches_pft(NZ)
    IF(iPlantBranchState_brch(NB,NZ).EQ.iLive)THEN
      IF(Hours2LeafOut_brch(NB,NZ).GT.HourReq2InitSStor4LeafOut(iPlantPhenolPattern_pft(NZ)))THEN
        LeafPetoMassC_brch(NB)=AZMAX1(LeafPetolBiomassC_brch(NB,NZ))
        DO NE=1,NumPlantChemElms
          NonstElm_loc(NE,NB)=AZMAX1(CanopyNonstElms_brch(NE,NB,NZ))
          TotNonstElm_loc(NE)=TotNonstElm_loc(NE)+NonstElm_loc(NE,NB)
        ENDDO
        TwoCompMassC=TwoCompMassC+LeafPetoMassC_brch(NB)
      ENDIF
    ENDIF
  ENDDO D300

  !!Nonst check
  DO NE=1,NumPlantChemElms
    mass_inital(NE)=SUM(CanopyNonstElms_brch(NE,1:NumOfBranches_pft(NZ),NZ))
  enddo
    !
  D305: DO NB=1,NumOfBranches_pft(NZ)
    IF(iPlantBranchState_brch(NB,NZ).EQ.iLive)THEN
        !leaf out criterion met
      IF(Hours2LeafOut_brch(NB,NZ).GT.HourReq2InitSStor4LeafOut(iPlantPhenolPattern_pft(NZ)))THEN
        IF(TwoCompMassC.GT.ZERO4Groth_pft(NZ) .AND. TotNonstElm_loc(ielmc).GT.ZERO4Groth_pft(NZ))THEN
          NonstElmGradt=TotNonstElm_loc(ielmc)*LeafPetoMassC_brch(NB)-NonstElm_loc(ielmc,NB)*TwoCompMassC
          XFRE(ielmc)=0.01_r8*NonstElmGradt/TwoCompMassC
          DO NE=2,NumPlantChemElms
            NonstElmGradt=TotNonstElm_loc(NE)*NonstElm_loc(ielmc,NB)-NonstElm_loc(NE,NB)*TotNonstElm_loc(ielmc)
            XFRE(NE)=0.01_r8*NonstElmGradt/TotNonstElm_loc(ielmc)
          ENDDO
          DO NE=1,NumPlantChemElms
            CanopyNonstElms_brch(NE,NB,NZ)=CanopyNonstElms_brch(NE,NB,NZ)+XFRE(NE)
          ENDDO
        ENDIF
      ENDIF
    ENDIF
  ENDDO D305

  DO NE=1,NumPlantChemElms
    mass_finale(NE)=SUM(CanopyNonstElms_brch(NE,1:NumOfBranches_pft(NZ),NZ))
  enddo
  
!=============================================================================
!     TRANSFER NON-STRUCTURAL C,N,P AMONG BRANCH STALK RESERVES
!     FROM NON-STRUCTURAL C,N,P CONCENTRATION DIFFERENCES
!
!     iPlantBranchState_brch=branch living flag: 0=alive,1=dead
!     SapwoodBiomassC_brch=stalk sapwood mass
!     WTRSVB,WTRSBN,WTRSBP=stalk reserve C,N,P mass
!     iPlantCalendar_brch(ipltcal_BeginSeedFill,=start of grain filling and setting max seed size
! the algorithm below is different from that for within canopy transfer between different branches
! why?
  DO NE=1,NumPlantChemElms
    mass_inital(NE)=sum(StalkRsrvElms_brch(NE,1:NumOfBranches_pft(NZ),NZ))
  ENDDO

  TotStalkMassC=0._r8
  TotStalkRsrv_loc(1:NumPlantChemElms)=0._r8
  D330: DO NB=1,NumOfBranches_pft(NZ)
    IF(iPlantBranchState_brch(NB,NZ).EQ.iLive)THEN
      IF(iPlantCalendar_brch(ipltcal_BeginSeedFill,NB,NZ).NE.0)THEN
        TotStalkMassC=TotStalkMassC+SapwoodBiomassC_brch(NB,NZ)
        DO NE=1,NumPlantChemElms
          TotStalkRsrv_loc(NE)=TotStalkRsrv_loc(NE)+StalkRsrvElms_brch(NE,NB,NZ)
        ENDDO
      ENDIF
    ENDIF
  ENDDO D330

  sumchk1=TotStalkRsrv_loc(ielmc)
  sumchk2=0._r8
  IF(TotStalkMassC.GT.ZERO4Groth_pft(NZ).AND.TotStalkRsrv_loc(ielmc).GT.ZERO4Groth_pft(NZ))THEN
    D335: DO NB=1,NumOfBranches_pft(NZ)
      IF(iPlantBranchState_brch(NB,NZ).EQ.iLive)THEN
        IF(iPlantCalendar_brch(ipltcal_BeginSeedFill,NB,NZ).NE.0)THEN
          StalkRsrvGradt=TotStalkRsrv_loc(ielmc)*SapwoodBiomassC_brch(NB,NZ)-StalkRsrvElms_brch(ielmc,NB,NZ)*TotStalkMassC
          XFRE(ielmc)=0.1_r8*StalkRsrvGradt/TotStalkMassC            
          StalkRsrvElms_brch(ielmc,NB,NZ)=StalkRsrvElms_brch(ielmc,NB,NZ)+XFRE(ielmc)
          sumchk2=sumchk2+StalkRsrvElms_brch(ielmc,NB,NZ)
          !based on stoichiometry gradient
          DO NE=2,NumPlantChemElms
            StalkRsrvGradt=TotStalkRsrv_loc(NE)*StalkRsrvElms_brch(ielmc,NB,NZ) &
              -StalkRsrvElms_brch(NE,NB,NZ)*TotStalkRsrv_loc(ielmc)
            XFRE(NE)=0.1_r8*StalkRsrvGradt/TotStalkRsrv_loc(ielmc)
            StalkRsrvElms_brch(NE,NB,NZ)=StalkRsrvElms_brch(NE,NB,NZ)+XFRE(NE)
          ENDDO
        ENDIF
      ENDIF
    ENDDO D335
  ENDIF
  DO NE=1,NumPlantChemElms
    mass_finale(NE)=sum(StalkRsrvElms_brch(NE,1:NumOfBranches_pft(NZ),NZ))
  ENDDO
  end associate   
  end subroutine WithinBranchElmTransfer    

!----------------------------------------------------------------------------------------------------
  subroutine RootMycoNonstTransfer(I,J,NZ)
  implicit none
  integer, intent(in) :: i,J,nz

  integer   :: NE,N,L
  real(r8) :: mass_inital(NumPlantChemElms)
  real(r8) :: mass_finale(NumPlantChemElms)
  real(r8) :: XFRE(NumPlantChemElms)  
  real(r8) :: WTRTD1,WTRTD2,CPOOLT
  real(r8) :: TwoCompMassC,NonstElmGradt
  associate(                                                    &
    NU                     => plt_site%NU                      ,& !input  :current soil surface layer number, [-]
    ZERO4Groth_pft         => plt_biom%ZERO4Groth_pft          ,& !input  :threshold zero for plang growth calculation, [-]
    NIXBotRootLayer_pft    => plt_morph%NIXBotRootLayer_pft    ,& !input  :maximum soil layer number for all root axes, [-]
    PopuRootMycoC_pvr      => plt_biom% PopuRootMycoC_pvr      ,& !input  :root layer C, [gC d-2]
    ZERO4LeafVar_pft       => plt_biom%ZERO4LeafVar_pft        ,& !input  :threshold zero for leaf calculation, [-]
    Myco_pft               => plt_morph%Myco_pft               ,& !input  :mycorrhizal type (no or yes),[-]
    RootMycoNonstElms_rpvr => plt_biom%RootMycoNonstElms_rpvr   & !inoput :root layer nonstructural element, [g d-2]
  )

!     TRANSFER NON-STRUCTURAL C,N,P BWTWEEN ROOT AND MYCORRHIZAE
!     IN EACH ROOTED SOIL LAYER FROM NON-STRUCTURAL C,N,P
!     CONCENTRATION DIFFERENCES
!
!     MY=mycorrhizal:1=no,2=yes
!     CPOOLR,ZPOOLR,PPOOLR=non-structural C,N,P mass in 1:root,2:mycorrhizae
!     WTRTD=1:root,2:mycorrhizal C mass
!     FSNK=min ratio of branch or mycorrhizae to root for calculating C transfer
!     FMYC=rate constant for root-mycorrhizal C,N,P exchange (h-1)
!  
  DO NE=1,NumPlantChemElms
    mass_inital(NE)=sum(RootMycoNonstElms_rpvr(NE,1:Myco_pft(NZ),NU:NIXBotRootLayer_pft(NZ),NZ))
  ENDDO

  !this enables extension to other mycorrhizae
  DO N=2,Myco_pft(NZ)
    D425: DO L=NU,NIXBotRootLayer_pft(NZ)
      IF(RootMycoNonstElms_rpvr(ielmc,ipltroot,L,NZ).GT.ZERO4Groth_pft(NZ) &
        .AND. PopuRootMycoC_pvr(ipltroot,L,NZ).GT.ZERO4LeafVar_pft(NZ))THEN
        !root
        WTRTD1=PopuRootMycoC_pvr(ipltroot,L,NZ)
        WTRTD2=AMIN1(PopuRootMycoC_pvr(ipltroot,L,NZ),AMAX1(FSNK &
        *PopuRootMycoC_pvr(ipltroot,L,NZ), PopuRootMycoC_pvr(N,L,NZ)))
        TwoCompMassC=WTRTD1+WTRTD2
        IF(TwoCompMassC.GT.ZERO4Groth_pft(NZ))THEN
          NonstElmGradt=(RootMycoNonstElms_rpvr(ielmc,ipltroot,L,NZ)*WTRTD2 &
            - RootMycoNonstElms_rpvr(ielmc,N,L,NZ)*WTRTD1)/TwoCompMassC
          XFRE(ielmc)=FMYC*NonstElmGradt
          RootMycoNonstElms_rpvr(ielmc,ipltroot,L,NZ)=RootMycoNonstElms_rpvr(ielmc,ipltroot,L,NZ)-XFRE(ielmc)
          RootMycoNonstElms_rpvr(ielmc,N,L,NZ)=RootMycoNonstElms_rpvr(ielmc,N,L,NZ)+XFRE(ielmc)

          CPOOLT=RootMycoNonstElms_rpvr(ielmc,ipltroot,L,NZ)+RootMycoNonstElms_rpvr(ielmc,N,L,NZ)
          IF(CPOOLT.GT.ZERO4Groth_pft(NZ))THEN
            !exchange based on stoichiometry gradient
            DO NE=2,NumPlantChemElms
              NonstElmGradt=(RootMycoNonstElms_rpvr(NE,ipltroot,L,NZ)*RootMycoNonstElms_rpvr(ielmc,N,L,NZ) &
                - RootMycoNonstElms_rpvr(NE,N,L,NZ)*RootMycoNonstElms_rpvr(ielmc,ipltroot,L,NZ))/CPOOLT
              XFRE(NE)=FMYC*NonstElmGradt
              RootMycoNonstElms_rpvr(NE,ipltroot,L,NZ)=RootMycoNonstElms_rpvr(NE,ipltroot,L,NZ)-XFRE(NE)
              RootMycoNonstElms_rpvr(NE,N,L,NZ)=RootMycoNonstElms_rpvr(NE,N,L,NZ)+XFRE(NE)
            ENDDO
          ENDIF
        ENDIF
      ENDIF
    ENDDO D425
  ENDDO

  DO NE=1,NumPlantChemElms
    mass_finale(NE)=sum(RootMycoNonstElms_rpvr(NE,1:Myco_pft(NZ),NU:NIXBotRootLayer_pft(NZ),NZ))
  ENDDO
  end associate
  end subroutine RootMycoNonstTransfer

!----------------------------------------------------------------------------------------------------
  subroutine SeasonStoreRootNonstTransfer(I,J,NZ)

  !for perennial plants only
  implicit none
  integer, intent(in) :: I,J,NZ
  integer  :: N,L,NE
  real(r8) :: CNL,CPL  
  real(r8) :: XFREX(NumPlantChemElms)
  real(r8) :: XFRE(NumPlantChemElms)  
  real(r8) :: mass_inital(NumPlantChemElms)
  real(r8) :: mass_finale(NumPlantChemElms)

  associate(                                                           &
    MaxSoiL4Root_pft          => plt_morph%MaxSoiL4Root_pft           ,& !input  :maximum soil layer number for all root axes,[-]
    NU                        => plt_site%NU                          ,& !input  :current soil surface layer number, [-]
    Myco_pft                  => plt_morph%Myco_pft                   ,& !input  :mycorrhizal type (no or yes),[-]
    ZERO                      => plt_site%ZERO                        ,& !input  :threshold zero for numerical stability, [-]
    iPlantTurnoverPattern_pft => plt_pheno%iPlantTurnoverPattern_pft  ,& !input  :phenologically-driven above-ground turnover: all, foliar only, none,[-]
    RootNonstructElmConc_rpvr => plt_biom%RootNonstructElmConc_rpvr   ,& !input  :root layer nonstructural C concentration, [g g-1]
    SeasonalNonstElms_pft     => plt_biom%SeasonalNonstElms_pft       ,& !inoput :plant stored nonstructural element at current step, [g d-2]
    RootMycoNonstElms_rpvr    => plt_biom%RootMycoNonstElms_rpvr       & !inoput :root layer nonstructural element, [g d-2]
  )
    DO NE=1,NumPlantChemElms
      mass_inital(NE)=sum(RootMycoNonstElms_rpvr(NE,1:Myco_pft(NZ),NU:MaxSoiL4Root_pft(NZ),NZ))+SeasonalNonstElms_pft(NE,NZ)
    ENDDO

    D5545: DO N=1,Myco_pft(NZ)
      D5550: DO L=NU,MaxSoiL4Root_pft(NZ)
        IF(RootNonstructElmConc_rpvr(ielmc,N,L,NZ).GT.ZERO)THEN
          CNL=RootNonstructElmConc_rpvr(ielmc,N,L,NZ)/(RootNonstructElmConc_rpvr(ielmc,N,L,NZ) &
            +RootNonstructElmConc_rpvr(ielmn,N,L,NZ)/CNKI)
          CPL=RootNonstructElmConc_rpvr(ielmc,N,L,NZ)/(RootNonstructElmConc_rpvr(ielmc,N,L,NZ) &
            +RootNonstructElmConc_rpvr(ielmp,N,L,NZ)/CPKI)
        ELSE
          CNL=0._r8
          CPL=0._r8
        ENDIF
        XFREX(ielmc) = RateK4RootSeaStorNonstEXfer(iPlantTurnoverPattern_pft(NZ))*AZMAX1(RootMycoNonstElms_rpvr(ielmc,N,L,NZ))
        XFREX(ielmn) = RateK4RootSeaStorNonstEXfer(iPlantTurnoverPattern_pft(NZ))*AZMAX1(RootMycoNonstElms_rpvr(ielmn,N,L,NZ))*(1.0_r8+CNL)
        XFREX(ielmp) = RateK4RootSeaStorNonstEXfer(iPlantTurnoverPattern_pft(NZ))*AZMAX1(RootMycoNonstElms_rpvr(ielmp,N,L,NZ))*(1.0_r8+CPL)
        XFRE(ielmc)  = AMIN1(XFREX(ielmc),XFREX(ielmn)/CNMN,XFREX(ielmp)/CPMN)
        XFRE(ielmn)  = AMIN1(XFREX(ielmn),XFRE(ielmc)*CNMX,XFREX(ielmp)*CNMX/CPMN*0.5_r8)
        XFRE(ielmp)  = AMIN1(XFREX(ielmp),XFRE(ielmc)*CPMX,XFREX(ielmn)*CPMX/CNMN*0.5_r8)

        DO NE=1,NumPlantChemElms
          RootMycoNonstElms_rpvr(NE,N,L,NZ)=RootMycoNonstElms_rpvr(NE,N,L,NZ)-XFRE(NE)
          SeasonalNonstElms_pft(NE,NZ)=SeasonalNonstElms_pft(NE,NZ)+XFRE(NE)
        ENDDO
      ENDDO D5550      
    ENDDO D5545
    DO NE=1,NumPlantChemElms
      mass_finale(NE)=sum(RootMycoNonstElms_rpvr(NE,1:Myco_pft(NZ),NU:MaxSoiL4Root_pft(NZ),NZ))+SeasonalNonstElms_pft(NE,NZ)
    ENDDO

  end associate
  end subroutine SeasonStoreRootNonstTransfer  

!----------------------------------------------------------------------------------------------------
  subroutine ShootRootElmTransfer(I,J,NZ,PTRT,RootSinkC_vr,Root1stSink_pvr,&
    Root2ndSink_pvr,RootSinkC)
  implicit none
  integer, intent(in) :: I,J,NZ
  real(r8), intent(in):: PTRT
  real(r8), INTENT(IN) :: RootSinkC_vr(pltpar%jroots,JZ1)
  real(r8), intent(in) :: Root1stSink_pvr(pltpar%jroots,JZ1,pltpar%MaxNumRootAxes)
  real(r8), intent(in) :: Root2ndSink_pvr(pltpar%jroots,JZ1,pltpar%MaxNumRootAxes)
  real(r8), intent(in) :: RootSinkC(pltpar%jroots)

  integer :: L,NB,N,NR,NE
  real(r8) :: ZPOOLS,ZPOOLT
  real(r8) :: ZPOOLB
  real(r8) :: ZPOOLD,EPOOLD
  REAL(R8) :: RootSinkWeight_pft(JZ1),BranchSinkWeight_pft(JP1)
  real(r8) :: CPOOLT
  real(r8) :: NonstElmRootE,NonstElmBrchE
  real(r8) :: NonstElmGradt
  real(r8) :: CPOOLB,CPOOLS
  real(r8) :: FWTC
  real(r8) :: FWTS
  real(r8) :: PPOOLB
  real(r8) :: PPOOLT
  real(r8) :: PTSHTR
  real(r8) :: PPOOLS
  real(r8) :: TwoCompMassC
  real(r8) :: WTRTLX
  real(r8) :: WTLSBX,WTLSBB
  real(r8) :: WTRTLR
  real(r8) :: XFRE(NumPlantChemElms)
  real(r8) :: mass_inital(NumPlantChemElms)
  real(r8) :: mass_finale(NumPlantChemElms)

  associate(                                                                   &
    NU                            => plt_site%NU                              ,& !input  :current soil surface layer number, [-]
    Myco_pft                      => plt_morph%Myco_pft                       ,& !input  :mycorrhizal type (no or yes),[-]
    RootMyco2ndStrutElms_rpvr     => plt_biom%RootMyco2ndStrutElms_rpvr       ,& !input  :root layer element secondary axes, [g d-2]
    FracShootLeafElmAlloc2Litr    => plt_allom%FracShootLeafElmAlloc2Litr     ,& !input  :woody element allocation, [-]
    FracRootElmAlloc2Litr         => plt_allom%FracRootElmAlloc2Litr          ,& !input  :C woody fraction in root,[-]
    RootMyco1stElm_raxs           => plt_biom%RootMyco1stElm_raxs             ,& !input  :root C primary axes, [g d-2]
    RootMyco1stStrutElms_rpvr     => plt_biom%RootMyco1stStrutElms_rpvr       ,& !input  :root layer element primary axes, [g d-2]
    LeafPetolBiomassC_brch        => plt_biom%LeafPetolBiomassC_brch          ,& !input  :plant branch leaf + sheath C, [g d-2]
    RootElms_pft                  => plt_biom%RootElms_pft                    ,& !input  :plant root element mass, [g d-2]
    ZERO4Groth_pft                => plt_biom%ZERO4Groth_pft                  ,& !input  :threshold zero for plang growth calculation, [-]
    iPlantBranchState_brch        => plt_pheno%iPlantBranchState_brch         ,& !input  :flag to detect branch death, [-]
    iPlantPhenolPattern_pft       => plt_pheno%iPlantPhenolPattern_pft        ,& !input  :plant growth habit: annual or perennial,[-]
    ShutRutNonstElmntConducts_pft => plt_pheno%ShutRutNonstElmntConducts_pft  ,& !input  :shoot-root rate constant for nonstructural C exchange, [h-1]
    RootCO2Autor_pvr              => plt_rbgc%RootCO2Autor_pvr                ,& !input  :root respiration constrained by O2, [g d-2 h-1]
    k_fine_litr                   => pltpar%k_fine_litr                       ,& !input  :fine litter complex id
    NIXBotRootLayer_rpft          => plt_morph%NIXBotRootLayer_rpft           ,& !input  :maximum soil layer number for root axes, [-]
    MaxSoiL4Root_pft              => plt_morph%MaxSoiL4Root_pft               ,& !input  :maximum soil layer number for all root axes,[-]
    NumRootAxes_pft               => plt_morph%NumRootAxes_pft                ,& !input  :root primary axis number,[-]
    NumOfBranches_pft             => plt_morph%NumOfBranches_pft              ,& !input  :number of branches,[-]
    RootMycoActiveBiomC_pvr       => plt_biom%RootMycoActiveBiomC_pvr         ,& !inoput :root layer structural C, [gC d-2]
    PopuRootMycoC_pvr             => plt_biom% PopuRootMycoC_pvr              ,& !inoput :root layer C, [gC d-2]
    RootMycoNonstElms_rpvr        => plt_biom%RootMycoNonstElms_rpvr          ,& !inoput :root layer nonstructural element, [g d-2]
    CanopyNonstElms_brch          => plt_biom%CanopyNonstElms_brch            ,& !inoput :branch nonstructural element, [g d-2]
    CanopyLeafShethC_pft          => plt_biom%CanopyLeafShethC_pft            ,& !inoput :canopy leaf + sheath C, [g d-2]
    ECO_ER_col                    => plt_bgcr%ECO_ER_col                      ,& !inoput :ecosystem respiration, [g d-2 h-1]
    Eco_AutoR_CumYr_col           => plt_bgcr%Eco_AutoR_CumYr_col              & !inoput :ecosystem autotrophic respiration, [g d-2 h-1]
  )
!     ROOT AND NODULE TOTALS
!
!     RootMycoActiveBiomC_pvr,WTRTD=active,actual root C mass
!     WTRT1,WTRT2=primary,secondary root C mass in soil layer
!     GrossResp_pft=total PFT respiration
!     RootCO2Autor_pvr=total root respiration
!     ECO_ER_col=ecosystem respiration
!     Eco_AutoR_CumYr_col=total autotrophic respiration
!
  D5445: DO N=1,Myco_pft(NZ)
    D5450: DO L=NU,MaxSoiL4Root_pft(NZ)
      RootMycoActiveBiomC_pvr(N,L,NZ) = 0._r8
      PopuRootMycoC_pvr(N,L,NZ)       = 0._r8
      D5460: DO NR=1,NumRootAxes_pft(NZ)
        RootMycoActiveBiomC_pvr(N,L,NZ) = RootMycoActiveBiomC_pvr(N,L,NZ)+RootMyco2ndStrutElms_rpvr(ielmc,N,L,NR,NZ)
        PopuRootMycoC_pvr(N,L,NZ)       = PopuRootMycoC_pvr(N,L,NZ)+RootMyco2ndStrutElms_rpvr(ielmc,N,L,NR,NZ) &
          +RootMyco1stStrutElms_rpvr(ielmc,N,L,NR,NZ)
      ENDDO D5460
      ECO_ER_col          = ECO_ER_col+RootCO2Autor_pvr(N,L,NZ)
      Eco_AutoR_CumYr_col = Eco_AutoR_CumYr_col+RootCO2Autor_pvr(N,L,NZ)
    ENDDO D5450

    DO  NR=1,NumRootAxes_pft(NZ)
      RootMycoActiveBiomC_pvr(N,NIXBotRootLayer_rpft(NR,NZ),NZ)=RootMycoActiveBiomC_pvr(N,NIXBotRootLayer_rpft(NR,NZ),NZ)&
        +RootMyco1stElm_raxs(ielmc,N,NR,NZ)
    ENDDO
  ENDDO D5445

!     TRANSFER NON-STRUCTURAL C,N,P BETWEEN ROOT AND SHOOT
!
!     SINK STRENGTH OF ROOTS IN EACH SOIL LAYER AS A FRACTION
!     OF TOTAL SINK STRENGTH OF ROOTS IN ALL SOIL LAYERS
!
!     iPlantPhenolPattern_pft=growth habit:0=annual,1=perennial from PFT file
!     WTLS,WTRT=total PFT leaf+petiole,root C mass
!     FWTC,FWTS,FWTR=canopy,root system,root layer sink weighting factor
!     RootSinkC_vr,RootSinkC=root layer,root system sink strength
!
!     IF(iPlantPhenolPattern_pft(NZ).EQ.iplt_perennial)THEN
!  why 2/3 here?
  IF(CanopyLeafShethC_pft(NZ).GT.ZERO4Groth_pft(NZ))THEN
    FWTC=AMIN1(1.0_r8,0.667_r8*RootElms_pft(ielmc,NZ)/CanopyLeafShethC_pft(NZ))
  ELSE
    FWTC=1.0_r8
  ENDIF
  IF(RootElms_pft(ielmc,NZ).GT.ZERO4Groth_pft(NZ))THEN
    FWTS=AMIN1(1.0_r8,CanopyLeafShethC_pft(NZ)/(0.667_r8*RootElms_pft(ielmc,NZ)))
  ELSE
    FWTS=1.0_r8
  ENDIF
!     ELSE
!     FWTC=1.0_r8
!     FWTS=1.0_r8
!     ENDIF
  D290: DO L=NU,MaxSoiL4Root_pft(NZ)
    IF(RootSinkC(ipltroot).GT.ZERO4Groth_pft(NZ))THEN
      RootSinkWeight_pft(L)=AZMAX1(RootSinkC_vr(ipltroot,L)/RootSinkC(ipltroot))
    ELSE
      RootSinkWeight_pft(L)=1.0_r8
    ENDIF
  ENDDO D290
!     RATE CONSTANT FOR TRANSFER IS SET FROM INPUT IN 'READQ'
!     BUT IS NOT USED FOR ANNUALS DURING GRAIN FILL
!
!     WTLS,WTLSB=total,branch PFT leaf+petiole C mass
!
  CanopyLeafShethC_pft(NZ)=0._r8
  D309: DO NB=1,NumOfBranches_pft(NZ)
    CanopyLeafShethC_pft(NZ)=CanopyLeafShethC_pft(NZ)+LeafPetolBiomassC_brch(NB,NZ)
  ENDDO D309
    
!
!     SINK STRENGTH OF BRANCHES IN EACH CANOPY AS A FRACTION
!     OF TOTAL SINK STRENGTH OF THE CANOPY
!
!     iPlantBranchState_brch=branch living flag: 0=alive,1=dead
!     iPlantPhenolPattern_pft=growth habit:0=annual,1=perennial from PFT file
!     iPlantCalendar_brch(ipltcal_SetSeedNumber,=end date for setting final seed number
!     BranchSinkWeight_pft=branch sink weighting factor
!     ShutRutNonstElmntConducts_pft=rate constant for equilibrating shoot-root nonstructural C concn from PFT file
!     PTRT=allocation to leaf+petiole used to modify ShutRutNonstElmntConducts_pft,in annuals
!     FWTC,FWTS,FWTR=canopy,root system,root layer sink weighting factor
!     FWOOD,FWOODN,FWOODP=C,N,P woody fraction in root:0=woody,1=non-woody
!     FWODB=C woody fraction in branch:0=woody,1=non-woody
!     FSNK=min ratio of branch or mycorrhizae to root for calculating C transfer
!     CPOOL,ZPOOL,PPOOL=non-structural C,N,P mass in branch
!     CPOOLR,ZPOOLR,PPOOLR=non-structural C,N,P mass in root
!
    DO NE=1,NumPlantChemElms
      mass_inital(NE)=sum(RootMycoNonstElms_rpvr(NE,1:Myco_pft(NZ),NU:MaxSoiL4Root_pft(NZ),NZ)) &
        +SUM(CanopyNonstElms_brch(NE,1:NumOfBranches_pft(NZ),NZ))
    ENDDO

  D310: DO NB=1,NumOfBranches_pft(NZ)
    IF(iPlantBranchState_brch(NB,NZ).EQ.iLive)THEN
      IF(CanopyLeafShethC_pft(NZ).GT.ZERO4Groth_pft(NZ))THEN
        BranchSinkWeight_pft(NB)=AZMAX1(LeafPetolBiomassC_brch(NB,NZ)/CanopyLeafShethC_pft(NZ))
      ELSE
        BranchSinkWeight_pft(NB)=1.0_r8
      ENDIF
      IF(iPlantPhenolPattern_pft(NZ).EQ.iplt_annual)THEN
        PTSHTR=ShutRutNonstElmntConducts_pft(NZ)*PTRT**0.167_r8
      ELSE
        PTSHTR=ShutRutNonstElmntConducts_pft(NZ)
      ENDIF

      D415: DO L=NU,MaxSoiL4Root_pft(NZ)
        WTLSBX       = LeafPetolBiomassC_brch(NB,NZ)*FracShootLeafElmAlloc2Litr(ielmc,k_fine_litr)*RootSinkWeight_pft(L)*FWTC
        WTRTLX       = RootMycoActiveBiomC_pvr(ipltroot,L,NZ)*FracRootElmAlloc2Litr(ielmc,k_fine_litr)*BranchSinkWeight_pft(NB)*FWTS
        WTLSBB       = AZMAX1(WTLSBX,FSNK*WTRTLX)
        WTRTLR       = AZMAX1(WTRTLX,FSNK*WTLSBX)
        TwoCompMassC = WTLSBB+WTRTLR
        IF(TwoCompMassC.GT.ZERO4Groth_pft(NZ))THEN
          CPOOLB                                      = AZMAX1(CanopyNonstElms_brch(ielmc,NB,NZ)*RootSinkWeight_pft(L))
          CPOOLS                                      = AZMAX1(RootMycoNonstElms_rpvr(ielmc,ipltroot,L,NZ)*BranchSinkWeight_pft(NB))
          NonstElmGradt                               = (CPOOLB*WTRTLR-CPOOLS*WTLSBB)/TwoCompMassC
          XFRE(ielmc)                                 = PTSHTR*NonstElmGradt
          CPOOLT                                      = CPOOLS+CPOOLB
          CanopyNonstElms_brch(ielmc,NB,NZ)           = CanopyNonstElms_brch(ielmc,NB,NZ)-XFRE(ielmc)
          RootMycoNonstElms_rpvr(ielmc,ipltroot,L,NZ) = RootMycoNonstElms_rpvr(ielmc,ipltroot,L,NZ)+XFRE(ielmc)

          !N & P tranfer based on stoichiometry ratio
          IF(CPOOLT.GT.ZERO4Groth_pft(NZ))THEN
            DO NE=2,NumPlantChemElms
              NonstElmBrchE = AZMAX1(CanopyNonstElms_brch(NE,NB,NZ)*RootSinkWeight_pft(L))
              NonstElmRootE = AZMAX1(RootMycoNonstElms_rpvr(NE,ipltroot,L,NZ)*BranchSinkWeight_pft(NB))
              NonstElmGradt = (NonstElmBrchE*CPOOLS-NonstElmRootE*CPOOLB)/CPOOLT
              XFRE(NE)      = PTSHTR*NonstElmGradt
              IF(XFRE(NE)>0._r8)then
                XFRE(NE)=AMIN1(CanopyNonstElms_brch(NE,NB,NZ),XFRE(NE))
              ELSE  
                XFRE(NE)=AMAX1(-RootMycoNonstElms_rpvr(NE,ipltroot,L,NZ),XFRE(NE))
              ENDIF
              CanopyNonstElms_brch(NE,NB,NZ)           = CanopyNonstElms_brch(NE,NB,NZ)-XFRE(NE)
              RootMycoNonstElms_rpvr(NE,ipltroot,L,NZ) = RootMycoNonstElms_rpvr(NE,ipltroot,L,NZ)+XFRE(NE)
            ENDDO
          ENDIF  
        ENDIF
      ENDDO D415
    ENDIF
  ENDDO D310
  DO NE=1,NumPlantChemElms
    mass_finale(NE)=sum(RootMycoNonstElms_rpvr(NE,1:Myco_pft(NZ),NU:MaxSoiL4Root_pft(NZ),NZ)) &
      +SUM(CanopyNonstElms_brch(NE,1:NumOfBranches_pft(NZ),NZ))
  ENDDO
  end associate
  end subroutine ShootRootElmTransfer  

!----------------------------------------------------------------------------------------------------
  subroutine PlantNonstElmTransfer(I,J,NZ,PTRT,RootSinkC_vr,Root1stSink_pvr,&
    Root2ndSink_pvr,RootSinkC,BegRemoblize)
  !
  !DESCRIPTION
  !transfer of nonstructural C/N/P 
  !
  implicit none
  integer,  intent(in) :: I,J,NZ
  integer,  intent(in) :: BegRemoblize
  real(r8), intent(in):: PTRT
  real(r8), INTENT(IN) :: RootSinkC_vr(pltpar%jroots,JZ1)
  real(r8), intent(in) :: Root1stSink_pvr(pltpar%jroots,JZ1,pltpar%MaxNumRootAxes)
  real(r8), intent(in) :: Root2ndSink_pvr(pltpar%jroots,JZ1,pltpar%MaxNumRootAxes)
  real(r8), intent(in) :: RootSinkC(pltpar%jroots)

!     begin_execution
  associate(                                                       &
    iPlantPhenolPattern_pft => plt_pheno%iPlantPhenolPattern_pft  ,& !input  :plant growth habit: annual or perennial,[-]
    Myco_pft                => plt_morph%Myco_pft                 ,& !input  :mycorrhizal type (no or yes),[-]
    NumOfBranches_pft       => plt_morph%NumOfBranches_pft         & !input  :number of branches,[-]
  )


  IF(NumOfBranches_pft(NZ).GT.1)THEN
    call WithinBranchElmTransfer(I,J,NZ)
  ENDIF

  IF(Myco_pft(NZ).GT.ipltroot)then
    call RootMycoNonstTransfer(I,J,NZ)
  endif    

!     TRANSFER ROOT NON-STRUCTURAL C,N,P TO SEASONAL STORAGE
!     IN PERENNIALS
!
  IF(BegRemoblize.EQ.itrue .AND. iPlantPhenolPattern_pft(NZ).NE.iplt_annual)THEN
    call SeasonStoreRootNonstTransfer(I,J,NZ)
  ENDIF


  call ShootRootElmTransfer(I,J,NZ,PTRT,RootSinkC_vr,Root1stSink_pvr,Root2ndSink_pvr,RootSinkC)

  end associate
  end subroutine PlantNonstElmTransfer

!----------------------------------------------------------------------------------------------------
  subroutine SeasonStoreShootTransfer(I,J,NB,NZ)

  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: NB,NZ
  integer :: ne
  real(r8) :: CWTRSV
  REAL(R8) :: CWTRSN,CWTRSP
  real(r8) :: CNR,CPR
  real(r8) :: CNL,CPL
  real(r8) :: XFREX(1:NumPlantChemElms)
  real(r8) :: XFRE(1:NumPlantChemElms)

  associate(                                                           &
    LeafPetoNonstElmConc_brch => plt_biom%LeafPetoNonstElmConc_brch   ,& !input  :branch nonstructural C concentration, [g d-2]
    ZERO4Groth_pft            => plt_biom%ZERO4Groth_pft              ,& !input  :threshold zero for plang growth calculation, [-]
    iPlantTurnoverPattern_pft => plt_pheno%iPlantTurnoverPattern_pft  ,& !input  :phenologically-driven above-ground turnover: all, foliar only, none,[-]
    SapwoodBiomassC_brch    => plt_biom%SapwoodBiomassC_brch      ,& !input  :branch live stalk C, [gC d-2]
    SeasonalNonstElms_pft     => plt_biom%SeasonalNonstElms_pft       ,& !inoput :plant stored nonstructural element at current step, [g d-2]
    CanopyNonstElms_brch      => plt_biom%CanopyNonstElms_brch        ,& !inoput :branch nonstructural element, [g d-2]
    StalkRsrvElms_brch        => plt_biom%StalkRsrvElms_brch           & !inoput :branch reserve element mass, [g d-2]
  )

  IF(SapwoodBiomassC_brch(NB,NZ).GT.ZERO4Groth_pft(NZ).AND. StalkRsrvElms_brch(ielmc,NB,NZ).GT.ZERO4Groth_pft(NZ))THEN
    CWTRSV = AZMAX1(StalkRsrvElms_brch(ielmc,NB,NZ)/SapwoodBiomassC_brch(NB,NZ))
    CWTRSN = AZMAX1(StalkRsrvElms_brch(ielmn,NB,NZ)/SapwoodBiomassC_brch(NB,NZ))
    CWTRSP = AZMAX1(StalkRsrvElms_brch(ielmp,NB,NZ)/SapwoodBiomassC_brch(NB,NZ))
    CNR    = CWTRSV/(CWTRSV+CWTRSN/CNKI)
    CPR    = CWTRSV/(CWTRSV+CWTRSP/CPKI)
  ELSE
    CNR=0._r8
    CPR=0._r8
  ENDIF
  XFREX(ielmc) = RateK4ShootSeaStoreNonstEXfer(iPlantTurnoverPattern_pft(NZ))*AZMAX1(StalkRsrvElms_brch(ielmc,NB,NZ))
  XFREX(ielmn) = RateK4ShootSeaStoreNonstEXfer(iPlantTurnoverPattern_pft(NZ))*AZMAX1(StalkRsrvElms_brch(ielmn,NB,NZ))*(1.0_r8+CNR)
  XFREX(ielmp) = RateK4ShootSeaStoreNonstEXfer(iPlantTurnoverPattern_pft(NZ))*AZMAX1(StalkRsrvElms_brch(ielmp,NB,NZ))*(1.0_r8+CPR)
  XFRE(ielmc)  = AMIN1(XFREX(ielmc),XFREX(ielmn)/CNMN,XFREX(ielmp)/CPMN)
  XFRE(ielmn)  = AMIN1(XFREX(ielmn),XFRE(ielmc)*CNMX,XFREX(ielmp)*CNMX/CPMN*0.5_r8)
  XFRE(ielmp)  = AMIN1(XFREX(ielmp),XFRE(ielmc)*CPMX,XFREX(ielmn)*CPMX/CNMN*0.5_r8)

  DO NE=1,NumPlantChemElms
    StalkRsrvElms_brch(NE,NB,NZ) = StalkRsrvElms_brch(NE,NB,NZ)-XFRE(NE)
    SeasonalNonstElms_pft(NE,NZ) = SeasonalNonstElms_pft(NE,NZ)+XFRE(NE)
  ENDDO
  IF(LeafPetoNonstElmConc_brch(ielmc,NB,NZ).GT.ZERO4Groth_pft(NZ))THEN
    CNL=LeafPetoNonstElmConc_brch(ielmc,NB,NZ)/(LeafPetoNonstElmConc_brch(ielmc,NB,NZ) &
      +LeafPetoNonstElmConc_brch(ielmn,NB,NZ)/CNKI)
    CPL=LeafPetoNonstElmConc_brch(ielmc,NB,NZ)/(LeafPetoNonstElmConc_brch(ielmc,NB,NZ) &
      +LeafPetoNonstElmConc_brch(ielmp,NB,NZ)/CPKI)
  ELSE
    CNL = 0._r8
    CPL = 0._r8
  ENDIF
  XFREX(ielmc)=RateK4ShootSeaStoreNonstEXfer(iPlantTurnoverPattern_pft(NZ))*AZMAX1(CanopyNonstElms_brch(ielmc,NB,NZ))
  XFREX(ielmn)=RateK4ShootSeaStoreNonstEXfer(iPlantTurnoverPattern_pft(NZ))*AZMAX1(CanopyNonstElms_brch(ielmn,NB,NZ))*(1.0_r8+CNL)
  XFREX(ielmp)=RateK4ShootSeaStoreNonstEXfer(iPlantTurnoverPattern_pft(NZ))*AZMAX1(CanopyNonstElms_brch(ielmp,NB,NZ))*(1.0_r8+CPL)
  XFRE(ielmc)=AMIN1(XFREX(ielmc),XFREX(ielmn)/CNMN,XFREX(ielmp)/CPMN)
  XFRE(ielmn)=AMIN1(XFREX(ielmn),XFRE(ielmc)*CNMX,XFREX(ielmp)*CNMX/CPMN*0.5_r8)
  XFRE(ielmp)=AMIN1(XFREX(ielmp),XFRE(ielmc)*CPMX,XFREX(ielmn)*CPMX/CNMN*0.5_r8)
  DO NE=1,NumPlantChemElms
    CanopyNonstElms_brch(NE,NB,NZ) = CanopyNonstElms_brch(NE,NB,NZ)-XFRE(NE)
    SeasonalNonstElms_pft(NE,NZ)   = SeasonalNonstElms_pft(NE,NZ)+XFRE(NE)
  ENDDO
  end associate
  end subroutine SeasonStoreShootTransfer

!----------------------------------------------------------------------------------------------------
  subroutine StalkRsrvShootNonstTransfer(I,J,NB,NZ)
  implicit none
  integer, intent(in) :: I,J,NB,NZ
  REAL(R8) :: ShootBiomC_brch
  real(r8) :: CPOOLT,CPOOLD
  real(r8) :: NonstGradt
  real(r8) :: XFRE(1:NumPlantChemElms)
  integer :: NE

  associate(                                                       &
    ZERO4Groth_pft          => plt_biom%ZERO4Groth_pft            ,& !input  :threshold zero for plang growth calculation, [-]
    LeafPetolBiomassC_brch  => plt_biom%LeafPetolBiomassC_brch    ,& !input  :plant branch leaf + sheath C, [g d-2]
    iPlantPhenolPattern_pft => plt_pheno%iPlantPhenolPattern_pft  ,& !input  :plant growth habit: annual or perennial,[-]
    SapwoodBiomassC_brch  => plt_biom%SapwoodBiomassC_brch    ,& !input  :branch live stalk C, [gC d-2]
    StalkRsrvElms_brch      => plt_biom%StalkRsrvElms_brch        ,& !inoput :branch reserve element mass, [g d-2]
    CanopyNonstElms_brch    => plt_biom%CanopyNonstElms_brch       & !inoput :branch nonstructural element, [g d-2]
  )

  ShootBiomC_brch=LeafPetolBiomassC_brch(NB,NZ)+SapwoodBiomassC_brch(NB,NZ)
  CPOOLT=CanopyNonstElms_brch(ielmc,NB,NZ)+StalkRsrvElms_brch(ielmc,NB,NZ)
  IF(ShootBiomC_brch.GT.ZERO4Groth_pft(NZ))THEN
    CPOOLD=(CanopyNonstElms_brch(ielmc,NB,NZ)*SapwoodBiomassC_brch(NB,NZ) &
      -StalkRsrvElms_brch(ielmc,NB,NZ)*LeafPetolBiomassC_brch(NB,NZ))/ShootBiomC_brch
    XFRE(ielmc)=FXFY(iPlantPhenolPattern_pft(NZ))*CPOOLD
    CanopyNonstElms_brch(ielmc,NB,NZ)=CanopyNonstElms_brch(ielmc,NB,NZ)-XFRE(ielmc)
    StalkRsrvElms_brch(ielmc,NB,NZ)=StalkRsrvElms_brch(ielmc,NB,NZ)+XFRE(ielmc)
  ENDIF

  IF(CPOOLT.GT.ZERO4Groth_pft(NZ))THEN
    DO NE=2,NumPlantChemElms
      NonstGradt=(CanopyNonstElms_brch(NE,NB,NZ)*StalkRsrvElms_brch(ielmc,NB,NZ) &
        -StalkRsrvElms_brch(NE,NB,NZ)*CanopyNonstElms_brch(ielmc,NB,NZ))/CPOOLT
      XFRE(NE)=FXFZ(iPlantPhenolPattern_pft(NZ))*NonstGradt
      CanopyNonstElms_brch(NE,NB,NZ)=CanopyNonstElms_brch(NE,NB,NZ)-XFRE(NE)
      StalkRsrvElms_brch(NE,NB,NZ)=StalkRsrvElms_brch(NE,NB,NZ)+XFRE(NE)
    ENDDO
  ENDIF
  end associate
  end subroutine StalkRsrvShootNonstTransfer

!----------------------------------------------------------------------------------------------------
  subroutine StalkRsrvRootNonstTransfer(I,J,NB,NZ)
  implicit none
  integer, intent(in) :: I,J,NB,NZ

  integer :: L,NE
  real(r8) :: WTPLTX,WTRTRX,CPOOLD,CPOOLT
  real(r8) :: XFRE(1:NumPlantChemElms)
  real(r8) :: NonstGradt

  associate(                                                       &
    NU                      => plt_site%NU                        ,& !input  :current soil surface layer number, [-]
    ZEROS2                  => plt_site%ZEROS2                    ,& !input  :threshold zero for numerical stability,[-]
    k_woody_litr            => pltpar%k_woody_litr                ,& !input  :woody litter complex id
    SapwoodBiomassC_brch  => plt_biom%SapwoodBiomassC_brch    ,& !input  :branch live stalk C, [gC d-2]
    VLSoilPoreMicP_vr       => plt_soilchem%VLSoilPoreMicP_vr     ,& !input  :volume of soil layer, [m3 d-2]
    RootMycoActiveBiomC_pvr => plt_biom%RootMycoActiveBiomC_pvr   ,& !input  :root layer structural C, [gC d-2]
    FracRootElmAlloc2Litr   => plt_allom%FracRootElmAlloc2Litr    ,& !input  :C woody fraction in root,[-]
    ZERO4Groth_pft          => plt_biom%ZERO4Groth_pft            ,& !input  :threshold zero for plang growth calculation, [-]
    iPlantPhenolPattern_pft => plt_pheno%iPlantPhenolPattern_pft  ,& !input  :plant growth habit: annual or perennial,[-]
    MaxSoiL4Root_pft        => plt_morph%MaxSoiL4Root_pft         ,& !input  :maximum soil layer number for all root axes,[-]
    StalkRsrvElms_brch      => plt_biom%StalkRsrvElms_brch        ,& !inoput :branch reserve element mass, [g d-2]
    RootMycoNonstElms_rpvr  => plt_biom%RootMycoNonstElms_rpvr     & !inoput :root layer nonstructural element, [g d-2]
  )
  D2050: DO L=NU,MaxSoiL4Root_pft(NZ)
    IF(VLSoilPoreMicP_vr(L).GT.ZEROS2)THEN
      WTRTRX = AMAX1(ZERO4Groth_pft(NZ),RootMycoActiveBiomC_pvr(ipltroot,L,NZ)*FracRootElmAlloc2Litr(ielmc,k_woody_litr))
      WTPLTX = WTRTRX+SapwoodBiomassC_brch(NB,NZ)
      IF(WTPLTX.GT.ZERO4Groth_pft(NZ))THEN
        CPOOLD=(RootMycoNonstElms_rpvr(ielmc,ipltroot,L,NZ)*SapwoodBiomassC_brch(NB,NZ) &
          -StalkRsrvElms_brch(ielmc,NB,NZ)*WTRTRX)/WTPLTX
        XFRE(ielmc)                                 = AZMAX1(FXFY(iPlantPhenolPattern_pft(NZ))*CPOOLD)
        RootMycoNonstElms_rpvr(ielmc,ipltroot,L,NZ) = RootMycoNonstElms_rpvr(ielmc,ipltroot,L,NZ)-XFRE(ielmc)
        StalkRsrvElms_brch(ielmc,NB,NZ)             = StalkRsrvElms_brch(ielmc,NB,NZ)+XFRE(ielmc)
        CPOOLT                                      = RootMycoNonstElms_rpvr(ielmc,ipltroot,L,NZ)+StalkRsrvElms_brch(ielmc,NB,NZ)
        IF(CPOOLT.GT.ZERO4Groth_pft(NZ))THEN
          DO NE=2,NumPlantChemElms
            NonstGradt=(RootMycoNonstElms_rpvr(NE,ipltroot,L,NZ)*StalkRsrvElms_brch(ielmc,NB,NZ)-&
              StalkRsrvElms_brch(NE,NB,NZ)*RootMycoNonstElms_rpvr(ielmc,ipltroot,L,NZ))/CPOOLT
            XFRE(NE)                                 = AZMAX1(FXFZ(iPlantPhenolPattern_pft(NZ))*NonstGradt)
            RootMycoNonstElms_rpvr(NE,ipltroot,L,NZ) = RootMycoNonstElms_rpvr(NE,ipltroot,L,NZ)-XFRE(NE)
            StalkRsrvElms_brch(NE,NB,NZ)             = StalkRsrvElms_brch(NE,NB,NZ)+XFRE(NE)
          ENDDO
        ENDIF
      ENDIF
    ENDIF
  ENDDO D2050
  end associate
  end subroutine StalkRsrvRootNonstTransfer
  ![tail]
end module PlantNonstElmDynMod