module PlantNonstElmDynMod
  use minimathmod,   only: safe_adb, AZMAX1, AZMIN1
  use data_kind_mod, only: r8 => DAT_KIND_R8
  use DebugToolMod,  only: PrintInfo
  use PlantMathFuncMod, only : ExchFluxLimiter
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
  public :: RepleteLowSeaStorByRoot
  public :: RepleteSeaStoreByStalk
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
    SapwoodBiomassC_brch    => plt_biom%SapwoodBiomassC_brch      ,& !input  :branch live stalk C, [gC d-2]
    CanopyLeafSheathC_brch  => plt_biom%CanopyLeafSheathC_brch    ,& !input  :plant branch leaf + sheath C, [g d-2]
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
    !     CanopyLeafSheathC_brch=leaf+petiole mass
    !     CPOOL,ZPOOL,PPOOL=non-structural C,N,P mass in branch
    !
  TwoCompMassC=0._r8
  TotNonstElm_loc(1:NumPlantChemElms)=0._r8
  D300: DO NB=1,NumOfBranches_pft(NZ)
    IF(iPlantBranchState_brch(NB,NZ).EQ.iLive)THEN
      IF(Hours2LeafOut_brch(NB,NZ).GT.HourReq2InitSStor4LeafOut(iPlantPhenolPattern_pft(NZ)))THEN
        LeafPetoMassC_brch(NB)=AZMAX1(CanopyLeafSheathC_brch(NB,NZ))
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
            NonstElmGradt = TotNonstElm_loc(NE)*NonstElm_loc(ielmc,NB)-NonstElm_loc(NE,NB)*TotNonstElm_loc(ielmc)
            XFRE(NE)      = 0.01_r8*NonstElmGradt/TotNonstElm_loc(ielmc)
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

  TotStalkMassC                        = 0._r8
  TotStalkRsrv_loc(1:NumPlantChemElms) = 0._r8
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
          StalkRsrvGradt                  = TotStalkRsrv_loc(ielmc)*SapwoodBiomassC_brch(NB,NZ)-StalkRsrvElms_brch(ielmc,NB,NZ)*TotStalkMassC
          XFRE(ielmc)                     = 0.1_r8*StalkRsrvGradt/TotStalkMassC
          StalkRsrvElms_brch(ielmc,NB,NZ) = StalkRsrvElms_brch(ielmc,NB,NZ)+XFRE(ielmc)
          sumchk2                         = sumchk2+StalkRsrvElms_brch(ielmc,NB,NZ)
          !based on stoichiometry gradient
          DO NE=2,NumPlantChemElms
            StalkRsrvGradt               = TotStalkRsrv_loc(NE)*StalkRsrvElms_brch(ielmc,NB,NZ)-StalkRsrvElms_brch(NE,NB,NZ)*TotStalkRsrv_loc(ielmc)
            XFRE(NE)                     = 0.1_r8*StalkRsrvGradt/TotStalkRsrv_loc(ielmc)
            StalkRsrvElms_brch(NE,NB,NZ) = StalkRsrvElms_brch(NE,NB,NZ)+XFRE(NE)
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
  real(r8) :: XFRE
  real(r8) :: WTRTD1,WTRTD2,CPOOLT
  real(r8) :: TwoCompMassC,NonstElmGradt
  character(len=*), parameter :: subname='RootMycoNonstTransfer'

  real(r8), parameter :: scalp=0.9999_r8
  associate(                                                    &
    NU                     => plt_site%NU                      ,& !input  :current soil surface layer number, [-]
    ZERO4Groth_pft         => plt_biom%ZERO4Groth_pft          ,& !input  :threshold zero for plang growth calculation, [-]
    NMaxRootBotLayer_pft   => plt_morph%NMaxRootBotLayer_pft   ,& !input  :maximum soil layer number for all root axes, [-]
    PopuRootMycoC_pvr      => plt_biom% PopuRootMycoC_pvr      ,& !input  :root layer C, [gC d-2]
    ZERO4LeafVar_pft       => plt_biom%ZERO4LeafVar_pft        ,& !input  :threshold zero for leaf calculation, [-]
    Myco_pft               => plt_morph%Myco_pft               ,& !input  :mycorrhizal type (no or yes),[-]
    RootMycoNonstElms_rpvr => plt_biom%RootMycoNonstElms_rpvr   & !inoput :root layer nonstructural element, [g d-2]
  )
  call PrintInfo('beg '//subname)
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
    mass_inital(NE)=sum(RootMycoNonstElms_rpvr(NE,1:Myco_pft(NZ),NU:NMaxRootBotLayer_pft(NZ),NZ))
  ENDDO

  !this enables extension to other mycorrhizae
  DO N=2,Myco_pft(NZ)
    D425: DO L=NU,NMaxRootBotLayer_pft(NZ)
      IF(RootMycoNonstElms_rpvr(ielmc,ipltroot,L,NZ).GT.ZERO4Groth_pft(NZ) &
        .AND. PopuRootMycoC_pvr(ipltroot,L,NZ).GT.ZERO4LeafVar_pft(NZ))THEN
        !root
        WTRTD1=PopuRootMycoC_pvr(ipltroot,L,NZ)
        WTRTD2=AMIN1(PopuRootMycoC_pvr(ipltroot,L,NZ),AMAX1(FSNK*PopuRootMycoC_pvr(ipltroot,L,NZ), PopuRootMycoC_pvr(N,L,NZ)))
        TwoCompMassC=WTRTD1+WTRTD2
        
        IF(TwoCompMassC.GT.ZERO4Groth_pft(NZ))THEN
          NonstElmGradt=(RootMycoNonstElms_rpvr(ielmc,ipltroot,L,NZ)*WTRTD2 &
            - RootMycoNonstElms_rpvr(ielmc,N,L,NZ)*WTRTD1)/TwoCompMassC
          XFRE = FMYC*NonstElmGradt

          call ExchFluxLimiter(RootMycoNonstElms_rpvr(ielmc,ipltroot,L,NZ),RootMycoNonstElms_rpvr(ielmc,N,L,NZ),XFRE)
          RootMycoNonstElms_rpvr(ielmc,ipltroot,L,NZ) = RootMycoNonstElms_rpvr(ielmc,ipltroot,L,NZ)-XFRE
          RootMycoNonstElms_rpvr(ielmc,N,L,NZ)        = RootMycoNonstElms_rpvr(ielmc,N,L,NZ)+XFRE

          CPOOLT=RootMycoNonstElms_rpvr(ielmc,ipltroot,L,NZ)+RootMycoNonstElms_rpvr(ielmc,N,L,NZ)
          IF(CPOOLT.GT.ZERO4Groth_pft(NZ))THEN
            !exchange based on stoichiometry gradient
            DO NE=2,NumPlantChemElms
              NonstElmGradt=(RootMycoNonstElms_rpvr(NE,ipltroot,L,NZ)*RootMycoNonstElms_rpvr(ielmc,N,L,NZ) &
                - RootMycoNonstElms_rpvr(NE,N,L,NZ)*RootMycoNonstElms_rpvr(ielmc,ipltroot,L,NZ))/CPOOLT
              XFRE = FMYC*NonstElmGradt
              
              call ExchFluxLimiter(RootMycoNonstElms_rpvr(NE,ipltroot,L,NZ),RootMycoNonstElms_rpvr(NE,N,L,NZ),XFRE)
              RootMycoNonstElms_rpvr(NE,ipltroot,L,NZ) = RootMycoNonstElms_rpvr(NE,ipltroot,L,NZ)-XFRE
              RootMycoNonstElms_rpvr(NE,N,L,NZ)        = RootMycoNonstElms_rpvr(NE,N,L,NZ)+XFRE
            ENDDO
          ENDIF
        ENDIF
      ENDIF
    ENDDO D425
  ENDDO

  DO NE=1,NumPlantChemElms
    mass_finale(NE)=sum(RootMycoNonstElms_rpvr(NE,1:Myco_pft(NZ),NU:NMaxRootBotLayer_pft(NZ),NZ))
  ENDDO
  call PrintInfo('end '//subname)
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
          CNL=RootNonstructElmConc_rpvr(ielmc,N,L,NZ)/(RootNonstructElmConc_rpvr(ielmc,N,L,NZ)+RootNonstructElmConc_rpvr(ielmn,N,L,NZ)/CNKI)
          CPL=RootNonstructElmConc_rpvr(ielmc,N,L,NZ)/(RootNonstructElmConc_rpvr(ielmc,N,L,NZ)+RootNonstructElmConc_rpvr(ielmp,N,L,NZ)/CPKI)
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
          call ExchFluxLimiter(RootMycoNonstElms_rpvr(NE,N,L,NZ),SeasonalNonstElms_pft(NE,NZ),XFRE(NE))
          RootMycoNonstElms_rpvr(NE,N,L,NZ) = RootMycoNonstElms_rpvr(NE,N,L,NZ)-XFRE(NE)
          SeasonalNonstElms_pft(NE,NZ)      = SeasonalNonstElms_pft(NE,NZ)+XFRE(NE)
        ENDDO
      ENDDO D5550      
    ENDDO D5545
    DO NE=1,NumPlantChemElms
      mass_finale(NE)=sum(RootMycoNonstElms_rpvr(NE,1:Myco_pft(NZ),NU:MaxSoiL4Root_pft(NZ),NZ))+SeasonalNonstElms_pft(NE,NZ)
    ENDDO

  end associate
  end subroutine SeasonStoreRootNonstTransfer  

!----------------------------------------------------------------------------------------------------
  subroutine ShootRootElmTransfer(I,J,NZ,GrothPART2LeafPetole,RootSinkC_vr,RootSinkC)
  !
  !The shoot-root exchange of nonstrucal element is done as follows:
  !1)sum up the strengths over all branches, and all roots LAYERS
  !2)compute weighting factor for each branch and each root layer
  !3)do exchange between one branch and all different root layers  
  use PlantMathFuncMod, only : get_zero_turg_ccpolt,update_osmo_turg_pressure
  implicit none
  integer, intent(in) :: I,J,NZ
  real(r8), intent(in):: GrothPART2LeafPetole
  real(r8), INTENT(IN) :: RootSinkC_vr(pltpar%jroots,JZ1)
  real(r8), intent(in) :: RootSinkC(pltpar%jroots)
  character(len=*), parameter :: subname='ShootRootElmTransfer'
  integer :: L,NB,N,NR,NE,L1
  real(r8) :: ZPOOLS,ZPOOLT
  real(r8) :: ZPOOLB
  real(r8) :: ZPOOLD,EPOOLD
  REAL(R8) :: RootSinkWeight_vr(JZ1),BranchSinkWeight_pft(JP1)
  real(r8) :: CPOOLT
  real(r8) :: NonstElmRootE,NonstElmBrchE
  real(r8) :: NonstElmGradt
  real(r8) :: CPOOLB,CPOOLS
  real(r8) :: FWTC, FWTS !canopy, root sink weighting factor
  real(r8) :: PPOOLB
  real(r8) :: PPOOLT
  real(r8) :: PTSHTR
  real(r8) :: PPOOLS
  real(r8) :: TwoCompMassC
  real(r8) :: WTRTLX
  real(r8) :: WTLSBX,WTLSBB
  real(r8) :: WTRTLR
  real(r8) :: XFRE
  real(r8) :: ZTOL
  real(r8) :: mass_inital(NumPlantChemElms)
  real(r8) :: mass_finale(NumPlantChemElms)
  real(r8) :: RootDepzMean,PSIOsmo,PSITurg
  real(r8) :: CPOLT(JZ1),CCPOLT
  associate(                                                                   &
    NU                            => plt_site%NU                              ,& !input  :current soil surface layer number, [-]
    Myco_pft                      => plt_morph%Myco_pft                       ,& !input  :mycorrhizal type (no or yes),[-]
    RootMyco2ndStrutElms_rpvr     => plt_biom%RootMyco2ndStrutElms_rpvr       ,& !input  :root layer element secondary axes, [g d-2]
    FracShootElmAllocm            => plt_allom%FracShootElmAllocm             ,& !input  :woody element allocation, [-]
    FracRootElmAllocm             => plt_allom%FracRootElmAllocm              ,& !input  :C woody fraction in root,[-]
    RootMyco1stElm_raxs           => plt_biom%RootMyco1stElm_raxs             ,& !input  :root C primary axes, [g d-2]
    Root1stDepz_raxes             => plt_morph%Root1stDepz_raxes              ,& !input  :root layer depth, [m]    
    RootMyco1stStrutElms_rpvr     => plt_biom%RootMyco1stStrutElms_rpvr       ,& !input  :root layer element primary axes, [g d-2]
    CanopyLeafSheathC_brch        => plt_biom%CanopyLeafSheathC_brch          ,& !input  :plant branch leaf + sheath C, [g d-2]
    CumSoilThickness_vr           => plt_site%CumSoilThickness_vr             ,& !input  :depth to bottom of soil layer from surface of grid cell, [m]    
    RootElms_pft                  => plt_biom%RootElms_pft                    ,& !input  :plant root element mass, [g d-2]
    PlantPopulation_pft           => plt_site%PlantPopulation_pft             ,& !input  :plant population, [d-2]    
    ZERO4Groth_pft                => plt_biom%ZERO4Groth_pft                  ,& !input  :threshold zero for plang growth calculation, [-]
    iPlantBranchState_brch        => plt_pheno%iPlantBranchState_brch         ,& !input  :flag to detect branch death, [-]
    iPlantPhenolPattern_pft       => plt_pheno%iPlantPhenolPattern_pft        ,& !input  :plant growth habit: annual or perennial,[-]
    ShootRootNonstElmConduts_pft  => plt_pheno%ShootRootNonstElmConduts_pft   ,& !input  :shoot-root rate constant for nonstructural C exchange, [h-1]
    RootCO2Autor_pvr              => plt_rbgc%RootCO2Autor_pvr                ,& !input  :root respiration constrained by O2, [g d-2 h-1]
    PSIRoot_pvr                   => plt_ew%PSIRoot_pvr                       ,& !input  :root total water potential, [Mpa]    
    PSIRootTurg_vr                => plt_ew%PSIRootTurg_vr                    ,& !input  :root turgor water potential, [Mpa]    
    PSIRootOSMO_vr                => plt_ew%PSIRootOSMO_vr                    ,& !input  :root osmotic water potential, [Mpa]    
    TKS_vr                        => plt_ew%TKS_vr                            ,& !input  :mean annual soil temperature, [K]    
    DLYR3                         => plt_site%DLYR3                           ,& !input  :vertical thickness of soil layer, [m]    
    OrganOsmoPsi0pt_pft           => plt_ew%OrganOsmoPsi0pt_pft               ,& !input  :Organ osmotic potential when canopy water potential = 0 MPa, [MPa]    
    k_fine_comp                   => pltpar%k_fine_comp                       ,& !input  :fine litter complex id
    flag2ndGrowth_pvr             => plt_morph%flag2ndGrowth_pvr              ,& !input  :flag for secondary growth of primary roots, [-]        
    NRoot1stTipLay_raxes          => plt_morph%NRoot1stTipLay_raxes           ,& !input  :maximum soil layer number for root axes, [-]
    MaxSoiL4Root_pft              => plt_morph%MaxSoiL4Root_pft               ,& !input  :maximum soil layer number for all root axes,[-]
    NumPrimeRootAxes_pft          => plt_morph%NumPrimeRootAxes_pft           ,& !input  :root primary axis number,[-]
    NMaxRootBotLayer_pft          => plt_morph%NMaxRootBotLayer_pft           ,& !input  :maximum soil layer number for all root axes, [-]    
    NumOfBranches_pft             => plt_morph%NumOfBranches_pft              ,& !input  :number of branches,[-]
    fRootTube_rpvr                 => plt_morph%fRootTube_rpvr                  ,& !input  :fraction of root for transport,[-]
    NGTopRootLayer_pft            => plt_morph%NGTopRootLayer_pft             ,& !input  :soil layer at planting depth, [-]    
    RootMycoActiveBiomC_pvr       => plt_biom%RootMycoActiveBiomC_pvr         ,& !inoput :root layer structural C, [gC d-2]
    PopuRootMycoC_pvr             => plt_biom%PopuRootMycoC_pvr               ,& !inoput :root layer C, [gC d-2]
    RootMycoNonstElms_rpvr        => plt_biom%RootMycoNonstElms_rpvr          ,& !inoput :root layer nonstructural element, [g d-2]
    CanopyNonstElms_brch          => plt_biom%CanopyNonstElms_brch            ,& !inoput :branch nonstructural element, [g d-2]
    CanopyLeafSheathC_pft         => plt_biom%CanopyLeafSheathC_pft           ,& !inoput :canopy leaf + sheath C, [g d-2]
    RootSinkWeight_pvr            => plt_morph%RootSinkWeight_pvr             ,& !output :Root nonst element sink profile, [d-2]
    ECO_ER_col                    => plt_bgcr%ECO_ER_col                      ,& !inoput :ecosystem respiration, [g d-2 h-1]
    ShootRootXferElm_pft          => plt_bgcr%ShootRootXferElm_pft            ,& !inoput :shoot-root nonstructural element transfer, [ g d-2 h-1]
    Eco_AutoR_CumYr_col           => plt_bgcr%Eco_AutoR_CumYr_col              & !inoput :ecosystem autotrophic respiration, [g d-2 h-1]
  )

  call PrintInfo('beg '//subname)
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
      D5460: DO NR=1,NumPrimeRootAxes_pft(NZ)
        RootMycoActiveBiomC_pvr(N,L,NZ) = RootMycoActiveBiomC_pvr(N,L,NZ)+RootMyco2ndStrutElms_rpvr(ielmc,N,L,NR,NZ)
        PopuRootMycoC_pvr(N,L,NZ)       = PopuRootMycoC_pvr(N,L,NZ)+RootMyco2ndStrutElms_rpvr(ielmc,N,L,NR,NZ)           
      ENDDO D5460
      ECO_ER_col          = ECO_ER_col+RootCO2Autor_pvr(N,L,NZ)
      Eco_AutoR_CumYr_col = Eco_AutoR_CumYr_col+RootCO2Autor_pvr(N,L,NZ)
    ENDDO D5450    
  ENDDO D5445
  
  !recognizing that non-elongation zone mostly serves as strorage and highway for resource transport to the tip.
  DO  NR=1,NumPrimeRootAxes_pft(NZ)
    L1=NRoot1stTipLay_raxes(NR,NZ)    
    DO L=NU,MaxSoiL4Root_pft(NZ)    
      if(flag2ndGrowth_pvr(L,NR,NZ))then
        RootMycoActiveBiomC_pvr(ipltroot,L1,NZ)=RootMycoActiveBiomC_pvr(ipltroot,L1,NZ)+RootMyco1stStrutElms_rpvr(ielmc,L,NR,NZ) &
          *fRootTube_rpvr(L,NR,NZ)
      ELSE
        RootMycoActiveBiomC_pvr(ipltroot,L1,NZ)=RootMycoActiveBiomC_pvr(ipltroot,L1,NZ)+RootMyco1stStrutElms_rpvr(ielmc,L,NR,NZ)  
      ENDIF
    ENDDO
  ENDDO

  DO L=NU,MaxSoiL4Root_pft(NZ)  
    DO NR=1,NumPrimeRootAxes_pft(NZ)
      PopuRootMycoC_pvr(ipltroot,L,NZ)  = PopuRootMycoC_pvr(ipltroot,L,NZ)+RootMyco1stStrutElms_rpvr(ielmc,L,NR,NZ)
    ENDDO
  ENDDO

  !     TRANSFER NON-STRUCTURAL C,N,P BETWEEN ROOT AND SHOOT
  !
  !     SINK STRENGTH OF ROOTS IN EACH SOIL LAYER AS A FRACTION
  !     OF TOTAL SINK STRENGTH OF ROOTS IN ALL SOIL LAYERS
  !
  !     iPlantPhenolPattern_pft=growth habit:0=annual,1=perennial from PFT file
  !     WTLS,WTRT=total PFT leaf+petiole,root C mass
  !     FWTC,FWTS,FWTR=canopy,root system,root layer 
  !     RootSinkC_vr,RootSinkC=root layer,root system sink strength
  !
  !     IF(iPlantPhenolPattern_pft(NZ).EQ.iplt_perennial)THEN
  !canopy sink weighting factor
  !  why 2/3 here?
  IF(CanopyLeafSheathC_pft(NZ).GT.ZERO4Groth_pft(NZ))THEN
    FWTC=AMIN1(1.0_r8,0.667_r8*RootElms_pft(ielmc,NZ)/CanopyLeafSheathC_pft(NZ))
  ELSE
    FWTC=1.0_r8
  ENDIF

  !root sink weighting factor
  IF(RootElms_pft(ielmc,NZ).GT.ZERO4Groth_pft(NZ))THEN
    FWTS=AMIN1(1.0_r8,CanopyLeafSheathC_pft(NZ)/(0.667_r8*RootElms_pft(ielmc,NZ)))
  ELSE
    FWTS=1.0_r8
  ENDIF

  IF(RootSinkC(ipltroot).GT.ZERO4Groth_pft(NZ))THEN  
    DO L=NU,MaxSoiL4Root_pft(NZ)  
      RootSinkWeight_pvr(L,NZ)=AZMAX1(RootSinkC_vr(ipltroot,L)/RootSinkC(ipltroot))
    ENDDO    
  ELSE
    ZTOL=CumSoilThickness_vr(MaxSoiL4Root_pft(NZ))-CumSoilThickness_vr(NU-1)  
    DO L=NU,MaxSoiL4Root_pft(NZ)  
      RootSinkWeight_pvr(L,NZ)=DLYR3(L)/ZTOL
    ENDDO
  ENDIF

  !     RATE CONSTANT FOR TRANSFER IS SET FROM INPUT IN 'READQ'
  !     BUT IS NOT USED FOR ANNUALS DURING GRAIN FILL
  !
  !     WTLS,WTLSB=total,branch PFT leaf+petiole C mass
  !
  CanopyLeafSheathC_pft(NZ)=0._r8
  D309: DO NB=1,NumOfBranches_pft(NZ)
    CanopyLeafSheathC_pft(NZ)=CanopyLeafSheathC_pft(NZ)+CanopyLeafSheathC_brch(NB,NZ)
  ENDDO D309
    
  !
  !     SINK STRENGTH OF BRANCHES IN EACH CANOPY AS A FRACTION
  !     OF TOTAL SINK STRENGTH OF THE CANOPY
  !
  !     iPlantBranchState_brch=branch living flag: 0=alive,1=dead
  !     iPlantPhenolPattern_pft=growth habit:0=annual,1=perennial from PFT file
  !     iPlantCalendar_brch(ipltcal_SetSeedNumber,=end date for setting final seed number
  !     BranchSinkWeight_pft=branch sink weighting factor
  !     GrothPART2LeafPetole=allocation to leaf+petiole used to modify ShootRootNonstElmConduts_pft,in annuals
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

  IF(iPlantPhenolPattern_pft(NZ).EQ.iplt_annual)THEN
    PTSHTR=ShootRootNonstElmConduts_pft(NZ)*GrothPART2LeafPetole**0.167_r8    !0.167=(1/6)~sqrt(length), GrothPART2LeafPetole is the main branch growth allocation to leaf+petiole
  ELSE
    PTSHTR=ShootRootNonstElmConduts_pft(NZ)
  ENDIF
  PTSHTR=AMIN1(PTSHTR,1._r8)
  
  D310: DO NB=1,NumOfBranches_pft(NZ)
    !exchange between branch NB and all root layers
    IF(iPlantBranchState_brch(NB,NZ).EQ.iLive)THEN
      IF(CanopyLeafSheathC_pft(NZ).GT.ZERO4Groth_pft(NZ))THEN
        BranchSinkWeight_pft(NB)=AZMAX1(CanopyLeafSheathC_brch(NB,NZ)/CanopyLeafSheathC_pft(NZ))
      ELSE
        BranchSinkWeight_pft(NB)=1.0_r8
      ENDIF

      !Roots at different depths are generally "wired" to the shoot (the source) like spokes 
      !on a wheel or branches on a river. They do not typically exchange carbon directly with each other deep underground.

      D415: DO L=NU,MaxSoiL4Root_pft(NZ)
        WTLSBX       = CanopyLeafSheathC_brch(NB,NZ)*FracShootElmAllocm(ielmc,k_fine_comp)*RootSinkWeight_pvr(L,NZ)*FWTC
        WTRTLX       = RootMycoActiveBiomC_pvr(ipltroot,L,NZ)*FracRootElmAllocm(ielmc,k_fine_comp)*BranchSinkWeight_pft(NB)*FWTS
        WTLSBB       = AZMAX1(WTLSBX,FSNK*WTRTLX)
        WTRTLR       = AZMAX1(WTRTLX,FSNK*WTLSBX)
        TwoCompMassC = WTLSBB+WTRTLR
        IF(TwoCompMassC.GT.ZERO4Groth_pft(NZ))THEN
          CPOOLB        = AZMAX1(CanopyNonstElms_brch(ielmc,NB,NZ)*RootSinkWeight_pvr(L,NZ))
          CPOOLS        = AZMAX1(RootMycoNonstElms_rpvr(ielmc,ipltroot,L,NZ)*BranchSinkWeight_pft(NB))
          CPOOLT        = CPOOLS+CPOOLB
          NonstElmGradt = (CPOOLB*WTRTLR-CPOOLS*WTLSBB)/TwoCompMassC
          XFRE          = PTSHTR*NonstElmGradt
          call ExchFluxLimiter(CanopyNonstElms_brch(ielmc,NB,NZ),RootMycoNonstElms_rpvr(ielmc,ipltroot,L,NZ),XFRE)

          ShootRootXferElm_pft(ielmc,NZ) = ShootRootXferElm_pft(ielmc,NZ)+XFRE          
          CanopyNonstElms_brch(ielmc,NB,NZ)           = CanopyNonstElms_brch(ielmc,NB,NZ)-XFRE
          RootMycoNonstElms_rpvr(ielmc,ipltroot,L,NZ) = RootMycoNonstElms_rpvr(ielmc,ipltroot,L,NZ)+XFRE

          !N & P tranfer based on stoichiometry ratio
          IF(CPOOLT.GT.ZERO4Groth_pft(NZ))THEN
            DO NE=2,NumPlantChemElms
              NonstElmBrchE = AZMAX1(CanopyNonstElms_brch(NE,NB,NZ)*RootSinkWeight_pvr(L,NZ))
              NonstElmRootE = AZMAX1(RootMycoNonstElms_rpvr(NE,ipltroot,L,NZ)*BranchSinkWeight_pft(NB))
              NonstElmGradt = (NonstElmBrchE*CPOOLS-NonstElmRootE*CPOOLB)/CPOOLT
              XFRE          = PTSHTR*NonstElmGradt

              call ExchFluxLimiter(CanopyNonstElms_brch(NE,NB,NZ),RootMycoNonstElms_rpvr(NE,ipltroot,L,NZ),XFRE)
              ShootRootXferElm_pft(NE,NZ)              = ShootRootXferElm_pft(NE,NZ)+XFRE
              CanopyNonstElms_brch(NE,NB,NZ)           = CanopyNonstElms_brch(NE,NB,NZ)-XFRE
              RootMycoNonstElms_rpvr(NE,ipltroot,L,NZ) = RootMycoNonstElms_rpvr(NE,ipltroot,L,NZ)+XFRE
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
  call PrintInfo('end '//subname)
  end associate
  end subroutine ShootRootElmTransfer  

!----------------------------------------------------------------------------------------------------
  subroutine PlantNonstElmTransfer(I,J,NZ,GrothPART2LeafPetole,RootSinkC_vr,RootSinkC,BegRemoblize)
  !
  !DESCRIPTION
  !transfer of nonstructural C/N/P 
  !
  implicit none
  integer,  intent(in) :: I,J,NZ
  integer,  intent(in) :: BegRemoblize
  real(r8), intent(in):: GrothPART2LeafPetole  !rate modifier for root-shoot nonstrucal material exchange 
  real(r8), INTENT(IN) :: RootSinkC_vr(pltpar%jroots,JZ1)
  real(r8), intent(in) :: RootSinkC(pltpar%jroots)
  real(r8) :: massE_beg(NumPlantChemElms)
  real(r8) :: massE_end(NumPlantChemElms)
  character(len=*), parameter :: subname='PlantNonstElmTransfer'

!     begin_execution
  associate(                                                       &
    iPlantPhenolPattern_pft => plt_pheno%iPlantPhenolPattern_pft  ,& !input  :plant growth habit: annual or perennial,[-]
    Myco_pft                => plt_morph%Myco_pft                 ,& !input  :mycorrhizal type (no or yes),[-]
    NumOfBranches_pft       => plt_morph%NumOfBranches_pft         & !input  :number of branches,[-]
  )
  call PrintInfo('beg '//subname)

!  call SumReserveBiomass(I,J,NZ,massE_beg)

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

  call ShootRootElmTransfer(I,J,NZ,GrothPART2LeafPetole,RootSinkC_vr,RootSinkC)

!  call SumReserveBiomass(I,J,NZ,massE_end)
  
  call PrintInfo('end '//subname)
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
    SapwoodBiomassC_brch      => plt_biom%SapwoodBiomassC_brch        ,& !input  :branch live stalk C, [gC d-2]
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
    call ExchFluxLimiter(StalkRsrvElms_brch(NE,NB,NZ),SeasonalNonstElms_pft(NE,NZ),XFRE(NE))
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
    call ExchFluxLimiter(CanopyNonstElms_brch(NE,NB,NZ),SeasonalNonstElms_pft(NE,NZ),XFRE(NE))
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
    CanopyLeafSheathC_brch  => plt_biom%CanopyLeafSheathC_brch    ,& !input  :plant branch leaf + sheath C, [g d-2]
    iPlantPhenolPattern_pft => plt_pheno%iPlantPhenolPattern_pft  ,& !input  :plant growth habit: annual or perennial,[-]
    SapwoodBiomassC_brch    => plt_biom%SapwoodBiomassC_brch      ,& !input  :branch live stalk C, [gC d-2]
    StalkRsrvElms_brch      => plt_biom%StalkRsrvElms_brch        ,& !inoput :branch reserve element mass, [g d-2]
    CanopyNonstElms_brch    => plt_biom%CanopyNonstElms_brch       & !inoput :branch nonstructural element, [g d-2]
  )

  ShootBiomC_brch = CanopyLeafSheathC_brch(NB,NZ)+SapwoodBiomassC_brch(NB,NZ)
  CPOOLT          = CanopyNonstElms_brch(ielmc,NB,NZ)+StalkRsrvElms_brch(ielmc,NB,NZ)
  IF(ShootBiomC_brch.GT.ZERO4Groth_pft(NZ))THEN
    CPOOLD=(CanopyNonstElms_brch(ielmc,NB,NZ)*SapwoodBiomassC_brch(NB,NZ) &
      -StalkRsrvElms_brch(ielmc,NB,NZ)*CanopyLeafSheathC_brch(NB,NZ))/ShootBiomC_brch
    XFRE(ielmc)                       = FXFY(iPlantPhenolPattern_pft(NZ))*CPOOLD
    CanopyNonstElms_brch(ielmc,NB,NZ) = CanopyNonstElms_brch(ielmc,NB,NZ)-XFRE(ielmc)
    StalkRsrvElms_brch(ielmc,NB,NZ)   = StalkRsrvElms_brch(ielmc,NB,NZ)+XFRE(ielmc)
  ENDIF

  IF(CPOOLT.GT.ZERO4Groth_pft(NZ))THEN
    DO NE=2,NumPlantChemElms
      NonstGradt=(CanopyNonstElms_brch(NE,NB,NZ)*StalkRsrvElms_brch(ielmc,NB,NZ) &
        -StalkRsrvElms_brch(NE,NB,NZ)*CanopyNonstElms_brch(ielmc,NB,NZ))/CPOOLT
      XFRE(NE)                       = FXFZ(iPlantPhenolPattern_pft(NZ))*NonstGradt
      call ExchFluxLimiter(CanopyNonstElms_brch(NE,NB,NZ),StalkRsrvElms_brch(NE,NB,NZ),XFRE(NE))
      CanopyNonstElms_brch(NE,NB,NZ) = CanopyNonstElms_brch(NE,NB,NZ)-XFRE(NE)
      StalkRsrvElms_brch(NE,NB,NZ)   = StalkRsrvElms_brch(NE,NB,NZ)+XFRE(NE)
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
    k_woody_comp            => pltpar%k_woody_comp                ,& !input  :woody litter complex id
    SapwoodBiomassC_brch    => plt_biom%SapwoodBiomassC_brch      ,& !input  :branch live stalk C, [gC d-2]
    VLSoilPoreMicP_vr       => plt_soilchem%VLSoilPoreMicP_vr     ,& !input  :volume of soil layer, [m3 d-2]
    RootMycoActiveBiomC_pvr => plt_biom%RootMycoActiveBiomC_pvr   ,& !input  :root layer structural C, [gC d-2]
    FracRootElmAllocm       => plt_allom%FracRootElmAllocm        ,& !input  :C woody fraction in root,[-]
    ZERO4Groth_pft          => plt_biom%ZERO4Groth_pft            ,& !input  :threshold zero for plang growth calculation, [-]
    iPlantPhenolPattern_pft => plt_pheno%iPlantPhenolPattern_pft  ,& !input  :plant growth habit: annual or perennial,[-]
    MaxSoiL4Root_pft        => plt_morph%MaxSoiL4Root_pft         ,& !input  :maximum soil layer number for all root axes,[-]
    StalkRsrvElms_brch      => plt_biom%StalkRsrvElms_brch        ,& !inoput :branch reserve element mass, [g d-2]
    RootMycoNonstElms_rpvr  => plt_biom%RootMycoNonstElms_rpvr     & !inoput :root layer nonstructural element, [g d-2]
  )
  D2050: DO L=NU,MaxSoiL4Root_pft(NZ)
    IF(VLSoilPoreMicP_vr(L).GT.ZEROS2)THEN
      !nonstructural chemical transport through sapwood/phloem
      WTRTRX = AMAX1(ZERO4Groth_pft(NZ),RootMycoActiveBiomC_pvr(ipltroot,L,NZ)*FracRootElmAllocm(ielmc,k_woody_comp))
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
            call ExchFluxLimiter(RootMycoNonstElms_rpvr(NE,ipltroot,L,NZ),StalkRsrvElms_brch(NE,NB,NZ),XFRE(NE))
            RootMycoNonstElms_rpvr(NE,ipltroot,L,NZ) = RootMycoNonstElms_rpvr(NE,ipltroot,L,NZ)-XFRE(NE)
            StalkRsrvElms_brch(NE,NB,NZ)             = StalkRsrvElms_brch(NE,NB,NZ)+XFRE(NE)
          ENDDO
        ENDIF
      ENDIF
    ENDIF
  ENDDO D2050
  end associate
  end subroutine StalkRsrvRootNonstTransfer
!----------------------------------------------------------------------------------------------------
  subroutine SumReserveBiomass(I,J,NZ,massE)
  implicit none
  integer, intent(in) :: I,J,NZ
  real(r8), intent(out) :: massE(NumPlantChemElms)
  integer :: L,NB,NE,N

  associate(                                                                 &
    NU                          => plt_site%NU                              ,& !input  :current soil surface layer number, [-]  
    StalkRsrvElms_brch          => plt_biom%StalkRsrvElms_brch              ,& !inoput :branch reserve element mass, [g d-2]  
    NumOfBranches_pft           => plt_morph%NumOfBranches_pft              ,& !input  :number of branches,[-]    
    Myco_pft                    => plt_morph%Myco_pft                       ,& !input  :mycorrhizal type (no or yes),[-]  
    MaxSoiL4Root_pft            => plt_morph%MaxSoiL4Root_pft               ,& !input  :maximum soil layer number for all root axes,[-]  
    SeasonalNonstElms_pft       => plt_biom%SeasonalNonstElms_pft           ,& !inoput :plant stored nonstructural element at current step, [g d-2]    
    RootMycoNonstElms_rpvr      => plt_biom%RootMycoNonstElms_rpvr          ,& !input :root layer nonstructural element, [g d-2]
    CanopyNonstElms_brch        => plt_biom%CanopyNonstElms_brch             & !inoput :branch nonstructural element, [g d-2]        
  )
  massE=0._r8

  DO L=NU,MaxSoiL4Root_pft(NZ)
    DO N=1,Myco_pft(NZ)    
      DO NE=1,NumPlantChemElms

        massE(NE) =massE(NE) + RootMycoNonstElms_rpvr(NE,N,L,NZ)
      ENDDO  
    ENDDO
  ENDDO

  DO NB=1,NumOfBranches_pft(NZ)
    DO NE=1,NumPlantChemElms
      massE(NE) =massE(NE) + CanopyNonstElms_brch(NE,NB,NZ)+StalkRsrvElms_brch(NE,NB,NZ)
    ENDDO
  ENDDO

  DO NE=1,NumPlantChemElms
    massE(NE) =massE(NE) + SeasonalNonstElms_pft(NE,NZ)
  ENDDO
  end associate
  end subroutine SumReserveBiomass
!----------------------------------------------------------------------------------------------------

  subroutine RepleteLowSeaStorByRoot(I,J,N,L,NZ)
  implicit none  
  integer, intent(in) :: I,J,N,L,NZ
  real(r8) :: FWTRT
  real(r8) :: WTRTTX
  real(r8) :: WTRVCX
  real(r8) :: WTRTLX
  real(r8) :: WTRTTT
  real(r8) :: NonstElmGradt
  real(r8) :: POOLEX
  real(r8) :: XFRC,CPOOLT
  integer :: NE
  character(len=*), parameter :: subname='RepleteLowSeaStorByRoot'

  associate(                                                              &
    RootElms_pft              => plt_biom%RootElms_pft                   ,& !input  :plant root element mass, [g d-2]
    ZERO4Groth_pft            => plt_biom%ZERO4Groth_pft                 ,& !input  :threshold zero for plang growth calculation, [-]
    RootMycoActiveBiomC_pvr   => plt_biom%RootMycoActiveBiomC_pvr        ,& !input  :root layer structural C, [gC d-2]
    NMaxRootBotLayer_pft      => plt_morph%NMaxRootBotLayer_pft          ,& !input  :maximum soil layer number for all root axes, [-]    
    RootMycoNonstElms_rpvr    => plt_biom%RootMycoNonstElms_rpvr         ,& !inoput :root layer nonstructural element, [g d-2]
    SeasonalNonstElms_pft     => plt_biom%SeasonalNonstElms_pft           & !inoput :plant stored nonstructural element at current step, [g d-2]
  )
  !     Note, 03/28/2024, Jinyun Tang: I need to take a careful look at the following code, particularly for WTRTLX and WTRTTX,
  !     which are essentially the same. However, assuming the flux is driven by C concentration of the seasonal storage
  !     and the nonstructural biomass in root, then WTRTTX should be defined with stalk volume.
  !     Question: where is seasonal storage located?
  !

  call PrintInfo('beg '//subname)
  IF(L.LE.NMaxRootBotLayer_pft(NZ))THEN   !within the root zone
    IF(RootMycoActiveBiomC_pvr(N,L,NZ).GT.ZERO4Groth_pft(NZ) &                      !has active root or mycorrhizae biomass
      .AND. RootElms_pft(ielmc,NZ).GT.ZERO4Groth_pft(NZ)     &                      !root has sufficient biomass
      .AND. SeasonalNonstElms_pft(ielmc,NZ).LT.XFRX*RootElms_pft(ielmc,NZ))THEN     !seasonal C storage is less than remobilizable root C
      FWTRT         = RootMycoActiveBiomC_pvr(N,L,NZ)/RootElms_pft(ielmc,NZ)
      WTRTLX        = RootMycoActiveBiomC_pvr(N,L,NZ)
      WTRTTX        = RootElms_pft(ielmc,NZ)*FWTRT
      WTRTTT        = WTRTLX+WTRTTX
      POOLEX        = AZMAX1(RootMycoNonstElms_rpvr(ielmc,N,L,NZ))
      WTRVCX        = AZMAX1(SeasonalNonstElms_pft(ielmc,NZ)*FWTRT)
      NonstElmGradt = (WTRVCX*WTRTLX-POOLEX*WTRTTX)/WTRTTT
      !root storage -> seasonal storage
      XFRC                                 = XFRY*AZMIN1(NonstElmGradt)
      RootMycoNonstElms_rpvr(ielmc,N,L,NZ) = RootMycoNonstElms_rpvr(ielmc,N,L,NZ)+XFRC
      SeasonalNonstElms_pft(ielmc,NZ)      = SeasonalNonstElms_pft(ielmc,NZ)-XFRC

      CPOOLT=WTRTLX+RootElms_pft(ielmc,NZ)

      DO NE=2,NumPlantChemElms
        POOLEX                            = AZMAX1(RootMycoNonstElms_rpvr(NE,N,L,NZ))
        WTRVCX                            = AZMAX1(SeasonalNonstElms_pft(NE,NZ)*FWTRT)
        !achor for seasonal storage is root, anchor for rootmyo is rootmyco actB
        NonstElmGradt                     = (WTRVCX*WTRTLX-POOLEX*WTRTTX)/CPOOLT
        XFRC                              = XFRY*AZMIN1(NonstElmGradt)
        call ExchFluxLimiter(SeasonalNonstElms_pft(NE,NZ),RootMycoNonstElms_rpvr(NE,N,L,NZ),XFRC)
        RootMycoNonstElms_rpvr(NE,N,L,NZ) = RootMycoNonstElms_rpvr(NE,N,L,NZ)+XFRC
        SeasonalNonstElms_pft(NE,NZ)      = SeasonalNonstElms_pft(NE,NZ)-XFRC
      ENDDO
    ENDIF
  ENDIF
  call PrintInfo('end '//subname)
  end associate      
  end subroutine RepleteLowSeaStorByRoot      
!----------------------------------------------------------------------------------------------------

  subroutine RepleteSeaStoreByStalk(I,J,NB,NZ)
 
  implicit none
  integer, intent(in) :: I,J,NB,NZ
  real(r8) :: NonstElmGradt
  real(r8) :: FracCanopyCinStalk
  real(r8) :: WVSTBX
  real(r8) :: WTRTTX,WTRSBX
  real(r8) :: WTRVCX,CPOOLT
  real(r8) :: XFRE(1:NumPlantChemElms)
  real(r8) :: ShootBiomC_brch
  integer :: NE
  character(len=*), parameter :: subname='RepleteSeaStoreByStalk'
  
  !   SapwoodBiomassC_brch,WVSTK=stalk,total stalk sapwood mass
  !   WTRT=total root mass
  !   WTRSVB,WTRSBN,WTRSBP=stalk reserve C,N,P mass
  !   XFRX=maximum storage C content for remobiln from stalk,root reserves
  !   XFRE(ielmc)=C transfer
  !   Q: why are nitrogen and phosphorus not transferred?
  associate(                                                       &
    StalkRsrvElms_brch      => plt_biom%StalkRsrvElms_brch        ,& !inoput :branch reserve element mass, [g d-2]
    CanopySapwoodC_pft      => plt_biom%CanopySapwoodC_pft        ,& !input  :canopy active stalk C, [g d-2]  
    SeasonalNonstElms_pft   => plt_biom%SeasonalNonstElms_pft     ,& !inoput :plant stored nonstructural element at current step, [g d-2]    
    SapwoodBiomassC_brch    => plt_biom%SapwoodBiomassC_brch      ,& !input  :branch live stalk C, [gC d-2]
    ZERO4Groth_pft          => plt_biom%ZERO4Groth_pft            ,& !input  :threshold zero for plang growth calculation, [-]    
    RootElms_pft            => plt_biom%RootElms_pft               & !input  :plant root element mass, [g d-2]
  )
  call PrintInfo('beg '//subname)
  IF(SapwoodBiomassC_brch(NB,NZ).GT.ZERO4Groth_pft(NZ)  &
    .AND. CanopySapwoodC_pft(NZ).GT.ZERO4Groth_pft(NZ)  &
    .AND. RootElms_pft(ielmc,NZ).GT.ZERO4Groth_pft(NZ) &
    .AND. StalkRsrvElms_brch(ielmc,NB,NZ).LE.XFRX*SapwoodBiomassC_brch(NB,NZ))THEN

    FracCanopyCinStalk = SapwoodBiomassC_brch(NB,NZ)/CanopySapwoodC_pft(NZ)
    WVSTBX             = SapwoodBiomassC_brch(NB,NZ)
    WTRTTX             = RootElms_pft(ielmc,NZ)*FracCanopyCinStalk
    ShootBiomC_brch    = WVSTBX+WTRTTX
    WTRSBX             = AZMAX1(StalkRsrvElms_brch(ielmc,NB,NZ))
    WTRVCX             = AZMAX1(SeasonalNonstElms_pft(ielmc,NZ)*FracCanopyCinStalk)
    NonstElmGradt      = (WTRVCX*WVSTBX-WTRSBX*WTRTTX)/ShootBiomC_brch
    !seasonal storage -> branch non-structrual
    XFRE(ielmc)                     = XFRY*AZMAX1(NonstElmGradt)
    StalkRsrvElms_brch(ielmc,NB,NZ) = StalkRsrvElms_brch(ielmc,NB,NZ)+XFRE(ielmc)
    SeasonalNonstElms_pft(ielmc,NZ) = SeasonalNonstElms_pft(ielmc,NZ)-XFRE(ielmc)

    CPOOLT=WVSTBX+RootElms_pft(ielmc,NZ)
    DO NE=2,NumPlantChemElms
      WTRSBX                            = AZMAX1(StalkRsrvElms_brch(ielmc,NB,NZ))
      WTRVCX                            = AZMAX1(SeasonalNonstElms_pft(NE,NZ)*FracCanopyCinStalk)
      !achor for seasonal storage is root, achor for stalkrsv is sap
      NonstElmGradt                     = (WTRVCX*WVSTBX-WTRSBX*WTRTTX)/CPOOLT
      XFRE(NE)                          = XFRY*AZMAX1(NonstElmGradt)
      call ExchFluxLimiter(SeasonalNonstElms_pft(NE,NZ),StalkRsrvElms_brch(ielmc,NB,NZ),XFRE(NE))
      StalkRsrvElms_brch(ielmc,NB,NZ)   = StalkRsrvElms_brch(ielmc,NB,NZ)+XFRE(NE)
      SeasonalNonstElms_pft(NE,NZ)      = SeasonalNonstElms_pft(NE,NZ)-XFRE(NE)
    ENDDO
  ENDIF
  call PrintInfo('end '//subname)
  end associate
  end subroutine RepleteSeaStoreByStalk  

  ![tail]
end module PlantNonstElmDynMod