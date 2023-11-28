module NoduleBGCMod

  use data_kind_mod, only : r8 => DAT_KIND_R8
  use minimathmod, only : safe_adb,AZMAX1
  use EcosimConst
  use PlantAPIData
  use GrosubPars
  implicit none
  private
  character(len=*),private, parameter :: mod_filename = &
  __FILE__
  public :: CanopyNoduleBiochemistry
  public :: RootNoduleBiomchemistry
  contains

!------------------------------------------------------------------------------------------

  subroutine CanopyNoduleBiochemistry(I,J,NZ,NB,TFN5,WFNG,UPNFC)
  implicit none
  integer, intent(in) :: I,J,NZ,NB
  real(r8), intent(in) :: TFN5,WFNG
  real(r8), intent(inout) :: UPNFC(JP1)
  integer :: M
  real(r8) :: ZADDN,ZPOOLD
  real(r8) :: CCC,CNC,CPC
  REAL(R8) :: cpoolt
  real(r8) :: CPOOLD
  real(r8) :: CCPOLN,CZPOLN
  real(r8) :: CPPOLN,CGNDL
  real(r8) :: CCNDLB
  real(r8) :: FCNPF
  real(r8) :: FXRNX
  real(r8) :: GRNDG
  real(r8) :: PPOOLD
  real(r8) :: PADDN,RCO2T
  real(r8) :: RCNDL,RMNDL,RXNDL
  real(r8) :: RGNDL,RSNDL
  real(r8) :: RGN2P,RGN2F
  real(r8) :: RUPNFB
  real(r8) :: RXNDLE(NumOfPlantChemElmnts)
  real(r8) :: RDNDLE(NumOfPlantChemElmnts)
  real(r8) :: RCNDLE(NumOfPlantChemElmnts)
  real(r8) :: RGNDG
  real(r8) :: RXNSNE(NumOfPlantChemElmnts)
  real(r8) :: RDNSNE(NumOfPlantChemElmnts)
  real(r8) :: RCNSNE(NumOfPlantChemElmnts)
  real(r8) :: SPNDLI
  real(r8) :: SPNDX
  real(r8) :: WTLSB1,WTNDB1,WTLSBT
  real(r8) :: XFRC,XFRN,XFRP
  REAL(R8) :: RCCC,RCCN,RCCP
  integer :: NE
!     begin_execution
  associate(                         &
    NU       =>  plt_site%NU       , &
    ZERO     =>  plt_site%ZERO     , &
    AREA3    =>  plt_site%AREA3    , &
    k_fine_litr=> pltpar%k_fine_litr,&
    CFOPE    =>  plt_soilchem%CFOPE, &
    iPlantNfixType   =>  plt_morph%iPlantNfixType  , &
    fTgrowCanP     =>  plt_pheno%fTgrowCanP    , &
    GrossResp_pft    =>  plt_bgcr%GrossResp_pft    , &
    ECO_ER_col     =>  plt_bgcr%ECO_ER_col     , &
    TCO2A    =>  plt_bgcr%TCO2A    , &
    Eco_AutoR_col     =>  plt_bgcr%Eco_AutoR_col     , &
    CO2NetFix_pft     =>  plt_bgcr%CO2NetFix_pft     , &
    LitterFallChemElmnt_pftvr     =>  plt_bgcr%LitterFallChemElmnt_pftvr     , &
    ifoliar  =>  pltpar%ifoliar    , &
    DMND     =>  plt_allom%DMND    , &
    CNND     =>  plt_allom%CNND    , &
    CPND     =>  plt_allom%CPND    , &
    LeafPetioleBiomassC_brch    =>  plt_biom%LeafPetioleBiomassC_brch    , &
    NonstructElmnt_brch   =>  plt_biom%NonstructElmnt_brch   , &
    NoduleNonstructElmnt_brch   =>  plt_biom%NoduleNonstructElmnt_brch   , &
    ZEROP    =>  plt_biom%ZEROP    , &
    ZEROL    =>  plt_biom%ZEROL    , &
    WTNDBE   =>  plt_biom%WTNDBE     &
  )
!     iPlantNfixType=N2 fixation: 4,5,6=rapid to slow canopy symbiosis
!
  IF(iPlantNfixType(NZ).GE.4)THEN
!
!     INITIAL INFECTION
!
!     WTNDB,WTNDBN,WTNDBP=bacterial C,N,P mass
!     WTNDI=initial bacterial mass at infection
!     AREA=grid cell area
!     CNND,CPND=bacterial N:C,P:C ratio from PFT file
!
    IF(WTNDBE(ielmc,NB,NZ).LE.0.0_r8)THEN
      WTNDBE(ielmc,NB,NZ)=WTNDBE(ielmc,NB,NZ)+WTNDI*AREA3(NU)
      WTNDBE(ielmn,NB,NZ)=WTNDBE(ielmn,NB,NZ)+WTNDI*AREA3(NU)*CNND(NZ)
      WTNDBE(ielmp,NB,NZ)=WTNDBE(ielmp,NB,NZ)+WTNDI*AREA3(NU)*CPND(NZ)
    ENDIF
!
!     O2-UNCONSTRAINED RESPIRATION RATES BY HETEROTROPHIC AEROBES
!     IN NODULE FROM SPECIFIC OXIDATION RATE, ACTIVE BIOMASS,
!     NON-STRUCTURAL C CONCENTRATION, MICROBIAL C:N:P FACTOR,
!     AND TEMPERATURE
!
!     WTNDB,WTNDBN,WTNDBP=bacterial C,N,P mass
!     CPOLNB,ZPOLNB,PPOLNB=nonstructural C,N,P in bacteria
!     CCPOLN,CZPOLN,CPPOLN=nonstructural C,N,P concn in bacteria
!     CNKI,CPKI=nonstructural N,P inhibition constant on growth
!     FCNPF=N,P constraint to bacterial activity
!
    IF(WTNDBE(ielmc,NB,NZ).GT.ZEROP(NZ))THEN
      CCPOLN=AZMAX1(NoduleNonstructElmnt_brch(ielmc,NB,NZ)/WTNDBE(ielmc,NB,NZ))
      CZPOLN=AZMAX1(NoduleNonstructElmnt_brch(ielmn,NB,NZ)/WTNDBE(ielmc,NB,NZ))
      CPPOLN=AZMAX1(NoduleNonstructElmnt_brch(ielmp,NB,NZ)/WTNDBE(ielmc,NB,NZ))
    ELSE
      CCPOLN=1.0_r8
      CZPOLN=1.0_r8
      CPPOLN=1.0_r8
    ENDIF
    IF(CCPOLN.GT.ZERO)THEN
      CCC=AZMAX1(AMIN1(1.0,safe_adb(CZPOLN,CZPOLN+CCPOLN*CNKI) &
        ,safe_adb(CPPOLN,CPPOLN+CCPOLN*CPKI)))
      CNC=AZMAX1(AMIN1(1.0,safe_adb(CCPOLN,CCPOLN+CZPOLN/CNKI)))
      CPC=AZMAX1(AMIN1(1.0,safe_adb(CCPOLN,CCPOLN+CPPOLN/CPKI)))
    ELSE
      CCC=0._r8
      CNC=0._r8
      CPC=0._r8
    ENDIF
    IF(WTNDBE(ielmc,NB,NZ).GT.ZEROP(NZ))THEN
      FCNPF=AMIN1(1.0_r8 &
        ,SQRT(WTNDBE(ielmn,NB,NZ)/(WTNDBE(ielmc,NB,NZ)*CNND(NZ))) &
        ,SQRT(WTNDBE(ielmp,NB,NZ)/(WTNDBE(ielmc,NB,NZ)*CPND(NZ))))
    ELSE
      FCNPF=1.0_r8
    ENDIF
    SPNDLI=CCPOLN/(CCPOLN+SPNDLK)
!
!     RESPIRATION FROM NON-STRUCTURAL C DETERMINED BY TEMPERATURE,
!     NON-STRUCTURAL C:N:P
!
!     RCNDL=respiration from non-structural C
!     CPOLNB,ZPOLNB,PPOLNB=nonstructural C,N,P in bacteria
!     VMXO=specific respiration rate by bacterial N2 fixers
!     WTNDB=bacterial C mass
!     fTgrowCanP=temperature function for canopy growth
!     FCNPF=N,P constraint to bacterial activity
!     WFNG=growth function of canopy water potential
!
    RCNDL=AZMAX1(AMIN1(NoduleNonstructElmnt_brch(ielmc,NB,NZ),VMXO*WTNDBE(ielmc,NB,NZ))*FCNPF*fTgrowCanP(NZ)*WFNG)
!     CPOOLNX=NoduleNonstructElmnt_brch(ielmc,NB,NZ)
!     VMXOX=VMXO*WTNDBE(ielmc,NB,NZ)*FCNPF*fTgrowCanP(NZ)*WFNG
!
!     NODULE MAINTENANCE RESPIRATION FROM SOIL TEMPERATURE,
!     NODULE STRUCTURAL N
!
!     RMNDL=bacterial maintenance respiration
!     RMPLT=specific maintenance respiration rate (g C g-1 N h-1)
!     TFN5=temperature function for canopy maintenance respiration
!     WTNDBN=bacterial N mass
!
    RMNDL=AZMAX1(RMPLT*TFN5*WTNDBE(ielmn,NB,NZ))*SPNDLI
!
!     NODULE GROWTH RESPIRATION FROM TOTAL - MAINTENANCE
!     IF > 0 DRIVES GROWTH, IF < 0 DRIVES REMOBILIZATION
!
!     RXNDL=difference between non-structural C respn and mntc respn
!     RGNDL=growth respiration unlimited by N,P
!     RSNDL=excess maintenance respiration
!
    RXNDL=RCNDL-RMNDL
    RGNDL=AZMAX1(RXNDL)
    RSNDL=AZMAX1(-RXNDL)
!
!     NODULE N2 FIXATION FROM GROWTH RESPIRATION, FIXATION ENERGY
!     REQUIREMENT AND NON-STRUCTURAL C:N:P PRODUCT INHIBITION,
!     CONSTRAINED BY MICROBIAL N REQUIREMENT
!
!     RGN2P=respiration requirement to maintain bacterial N:C ratio
!     WTNDB,WTNDBN=bacterial C,N mass
!     CNND=bacterial N:C ratio from PFT file
!     EN2F=N fixation yield from C oxidation (g N g-1 C)
!     RGNDL=growth respiration unlimited by N,P
!     RGN2F=respiration for N2 fixation
!     RUPNFB,UPNFC=branch,total N2 fixation
!
    RGN2P=AZMAX1(WTNDBE(ielmc,NB,NZ)*CNND(NZ)-WTNDBE(ielmn,NB,NZ))/EN2F
    IF(RGNDL.GT.ZEROP(NZ))THEN
      RGN2F=RGNDL*RGN2P/(RGNDL+RGN2P)
    ELSE
      RGN2F=0._r8
    ENDIF
    RUPNFB=RGN2F*EN2F
    UPNFC(NZ)=UPNFC(NZ)+RUPNFB
!
!     NODULE C,N,P REMOBILIZATION AND DECOMPOSITION AND LEAKAGE
!
!     RCCC,RCCN,RCCP=remobilization coefficient for C,N,P
!     RCCZN,RCCYN=min,max fractions for bacteria C recycling
!     RCCXN,RCCQN=max fractions for bacteria N,P recycling
!     LeafPetioleBiomassC_brch=leaf+petiole mass
!     CCNDLB=bacteria:leaf+petiole ratio
!     RDNDBX=effect of CCNDLB on bacteria decomposition rate
!     SPNDX=specific bacterial decomposition rate at current CCNDLB
!     SPNDL=specific decomposition rate by bacterial N2 fixers
!     WTNDB,WTNDBN,WTNDBP=bacterial C,N,P mass
!     RXNDLE(ielmc),RXNDLE(ielmn),RXNDLE(ielmp)=bacterial C,N,P loss from decomposition
!     RDNDLE(ielmc),RDNDLE(ielmn),RDNDLE(ielmp)=bacterial C,N,P decomposition to litterfall
!     RCNDLE(ielmc),RCNDLE(ielmn),RCNDLE(ielmp)=bacterial C,N,P decomposition to recycling
!
    RCCC=RCCZN+CCC*RCCYN
    RCCN=CNC*RCCXN
    RCCP=CPC*RCCQN
    SPNDX=SPNDL*SQRT(fTgrowCanP(NZ)*WFNG)
    DO NE=1,NumOfPlantChemElmnts
      RXNDLE(NE)=SPNDX*WTNDBE(NE,NB,NZ)
    ENDDO

    RDNDLE(ielmc)=RXNDLE(ielmc)*(1.0_r8-RCCC)
    RDNDLE(ielmn)=RXNDLE(ielmn)*(1.0_r8-RCCC)*(1.0_r8-RCCN)
    RDNDLE(ielmp)=RXNDLE(ielmp)*(1.0_r8-RCCC)*(1.0_r8-RCCP)

    DO NE=1,NumOfPlantChemElmnts
      RCNDLE(NE)=RXNDLE(NE)-RDNDLE(NE)
    ENDDO
!
!     TOTAL NON-STRUCTURAL C,N,P USED IN NODULE GROWTH
!     AND GROWTH RESPIRATION DEPENDS ON GROWTH YIELD
!     ENTERED IN 'READQ'
!
!     CGNDL=total non-structural C used in bacterial growth and growth respiration
!     CPOLNB,ZPOLNB,PPOLNB=nonstructural C,N,P in bacteria
!     RMNDL=bacterial maintenance respiration
!     RCNDL=respiration from non-structural C
!     RCNDLE(ielmc)=bacterial C decomposition to recycling
!     RGNDL=growth respiration ltd by O2
!     RGN2F=respiration for N2 fixation
!     GRNDG=bacterial growth
!     DMND=bacterial growth yield
!     RGNDG=bacterial respiration for growth and N2 fixation
!     ZADDN,PADDN=nonstructural N,P used in growth
!     CNND,CPND=bacterial N:C,P:C ratio from PFT file
!     CCPOLN,CZPOLN,CPPOLN=nonstructural C,N,P concn in bacteria
!     CZKM,CPKM=Km for nonstructural N,P uptake by bacteria
!
    CGNDL=AMIN1(NoduleNonstructElmnt_brch(ielmc,NB,NZ)-AMIN1(RMNDL,RCNDL) &
      -RGN2F+RCNDLE(ielmc),(RGNDL-RGN2F)/(1.0_r8-DMND(NZ)))
    GRNDG=CGNDL*DMND(NZ)
    RGNDG=RGN2F+CGNDL*(1.0_r8-DMND(NZ))
    ZADDN=AZMAX1(AMIN1(NoduleNonstructElmnt_brch(ielmn,NB,NZ),GRNDG*CNND(NZ)))*CZPOLN/(CZPOLN+CZKM)
    PADDN=AZMAX1(AMIN1(NoduleNonstructElmnt_brch(ielmp,NB,NZ),GRNDG*CPND(NZ)))*CPPOLN/(CPPOLN+CPKM)
!
!     NODULE SENESCENCE
!
!     RSNDL=excess maintenance respiration
!     WTNDB,WTNDBN,WTNDBP=bacterial C,N,P mass
!     RCCC,RCCN,RCCP=remobilization coefficient for C,N,P
!     RXNSNE(ielmc),RXNSNE(ielmc),RXNSNE(ielmp)=bacterial C,N,P loss from senescence
!     RDNSNE(ielmc),RDNSNE(ielmc),RDNSNE(ielmp)=bacterial C,N,P senescence to litterfall
!     RCNSNE(ielmc),RCNSNE(ielmc),RCNSNE(ielmp)=bacterial C,N,P senescence to recycling
!
    IF(RSNDL.GT.0.0.AND.WTNDBE(ielmc,NB,NZ).GT.ZEROP(NZ).AND.RCCC.GT.ZERO)THEN
      RXNSNE(ielmc)=RSNDL/RCCC
      RXNSNE(ielmn)=RXNSNE(ielmc)*WTNDBE(ielmn,NB,NZ)/WTNDBE(ielmc,NB,NZ)
      RXNSNE(ielmp)=RXNSNE(ielmc)*WTNDBE(ielmp,NB,NZ)/WTNDBE(ielmc,NB,NZ)

      RDNSNE(ielmc)=RXNSNE(ielmc)*(1.0_r8-RCCC)
      RDNSNE(ielmn)=RXNSNE(ielmn)*(1.0_r8-RCCC)*(1.0_r8-RCCN)
      RDNSNE(ielmp)=RXNSNE(ielmp)*(1.0_r8-RCCC)*(1.0_r8-RCCP)
      DO NE=1,NumOfPlantChemElmnts
        RCNSNE(NE)=RXNSNE(NE)-RDNSNE(NE)
      ENDDO
    ELSE
      RXNSNE(1:NumOfPlantChemElmnts)=0._r8
      RDNSNE(1:NumOfPlantChemElmnts)=0._r8
      RCNSNE(1:NumOfPlantChemElmnts)=0._r8
    ENDIF
!
!     TOTAL NODULE RESPIRATION
!
!     RCO2T=total C respiration
!     RMNDL=bacterial maintenance respiration
!     RCNDL=respiration from non-structural C
!     RGNDG=bacterial respiration for growth and N2 fixation
!     RCNSNE(ielmc)=bacterial C senescence to recycling
!     GrossResp_pft,TCO2A=total,above-ground PFT respiration
!     CO2NetFix_pft=PFT net CO2 fixation
!     ECO_ER_col=ecosystem respiration
!     Eco_AutoR_col=total autotrophic respiration
!
    RCO2T=AMIN1(RMNDL,RCNDL)+RGNDG+RCNSNE(ielmc)
    GrossResp_pft(NZ)=GrossResp_pft(NZ)-RCO2T
    TCO2A(NZ)=TCO2A(NZ)-RCO2T
    CO2NetFix_pft(NZ)=CO2NetFix_pft(NZ)-RCO2T
    ECO_ER_col=ECO_ER_col-RCO2T
    Eco_AutoR_col=Eco_AutoR_col-RCO2T
!
!     NODULE LITTERFALL CAUSED BY REMOBILIZATION
!
!     CSNC,ZSNC,PSNC=C,N,P litterfall from decomposition and senescence
!     CFOPC,CFOPN,CFOPC=fraction of litterfall C,N,P allocated to litter components
!     RDNDLE(ielmc),RDNDLE(ielmn),RDNDLE(ielmp)=bacterial C,N,P decomposition to litterfall
!     RDNSNE(ielmc),RDNSNE(ielmc),RDNSNE(ielmp)=bacterial C,N,P senescence to litterfall
!
    D6470: DO M=1,jsken
      DO NE=1,NumOfPlantChemElmnts
        LitterFallChemElmnt_pftvr(NE,M,k_fine_litr,0,NZ)=LitterFallChemElmnt_pftvr(NE,M,k_fine_litr,0,NZ) &
          +CFOPE(NE,ifoliar,M,NZ)*(RDNDLE(NE)+RDNSNE(NE))
      ENDDO
    ENDDO D6470
!
!     CONSUMPTION OF NON-STRUCTURAL C,N,P BY NODULE
!
!     CPOLNB,ZPOLNB,PPOLNB=nonstructural C,N,P in bacteria
!     RMNDL=bacterial maintenance respiration
!     RCNDL=respiration from non-structural C
!     RGN2F=respiration for N2 fixation
!     CGNDL=total non-structural C used in bacterial growth and growth respiration
!     RCNDLE(ielmc),RCNDLE(ielmn),RCNDLE(ielmp)=bacterial C,N,P decomposition to recycling
!     RCNSNE(ielmc),RCNSNE(ielmc),RCNSNE(ielmp)=bacterial C,N,P senescence to recycling
!     ZADDN,PADDN=nonstructural N,P used in growth
!     RUPNFB=branch N2 fixation
!
    NoduleNonstructElmnt_brch(ielmc,NB,NZ)=NoduleNonstructElmnt_brch(ielmc,NB,NZ)-AMIN1(RMNDL,RCNDL)-RGN2F-CGNDL+RCNDLE(ielmc)
    NoduleNonstructElmnt_brch(ielmn,NB,NZ)=NoduleNonstructElmnt_brch(ielmn,NB,NZ)-ZADDN+RCNDLE(ielmn)+RCNSNE(ielmn)+RUPNFB
    NoduleNonstructElmnt_brch(ielmp,NB,NZ)=NoduleNonstructElmnt_brch(ielmp,NB,NZ)-PADDN+RCNDLE(ielmp)+RCNSNE(ielmp)
!
!     UPDATE STATE VARIABLES FOR NODULE C, N, P
!
!     WTNDB,WTNDBN,WTNDBP=bacterial C,N,P mass
!     GRNDG=bacterial growth
!     RXNDLE(ielmc),RXNDLE(ielmn),RXNDLE(ielmp)=bacterial C,N,P loss from decomposition
!     RXNSNE(ielmc),RXNSNE(ielmc),RXNSNE(ielmp)=bacterial C,N,P loss from senescence
!     ZADDN,PADDN=nonstructural N,P used in growth
!
    WTNDBE(ielmc,NB,NZ)=WTNDBE(ielmc,NB,NZ)+GRNDG-RXNDLE(ielmc)-RXNSNE(ielmc)
    WTNDBE(ielmn,NB,NZ)=WTNDBE(ielmn,NB,NZ)+ZADDN-RXNDLE(ielmn)-RXNSNE(ielmn)
    WTNDBE(ielmp,NB,NZ)=WTNDBE(ielmp,NB,NZ)+PADDN-RXNDLE(ielmp)-RXNSNE(ielmp)
!
!     TRANSFER NON-STRUCTURAL C,N,P BETWEEN BRANCH AND NODULES
!     FROM NON-STRUCTURAL C,N,P CONCENTRATION DIFFERENCES
!
!     CPOOL,ZPOOL,PPOOL=branch non-structural C,N,P mass
!     LeafPetioleBiomassC_brch=leaf+petiole C mass
!     WTNDB=bacterial C mass
!     WTNDI=initial bacterial mass at infection
!     FXRN=rate constant for plant-bacteria nonstructural C,N,P exchange
!     CCNGB=parameter to calculate nonstructural C,N,P exchange
!     CCNDLB=bacteria:leaf+petiole ratio
!     XFRC,XFRN,XFRC=nonstructural C,N,P transfer
!     CPOLNB,ZPOLNB,PPOLNB=nonstructural C,N,P in bacteria
!
    IF(NonstructElmnt_brch(ielmc,NB,NZ).GT.ZEROP(NZ).AND.LeafPetioleBiomassC_brch(NB,NZ).GT.ZEROL(NZ))THEN
      CCNDLB=WTNDBE(ielmc,NB,NZ)/LeafPetioleBiomassC_brch(NB,NZ)
      WTLSB1=LeafPetioleBiomassC_brch(NB,NZ)
      WTNDB1=AMIN1(LeafPetioleBiomassC_brch(NB,NZ),AMAX1(WTNDI*AREA3(NU),WTNDBE(ielmc,NB,NZ)))
      WTLSBT=WTLSB1+WTNDB1
      IF(WTLSBT.GT.ZEROP(NZ))THEN
        FXRNX=FXRN(iPlantNfixType(NZ))/(1.0+CCNDLB/CCNGB)
!    2/(1.0+CCNDLB/(CCNGB*FXRN(iPlantNfixType(NZ))))
        CPOOLD=(NonstructElmnt_brch(ielmc,NB,NZ)*WTNDB1-NoduleNonstructElmnt_brch(ielmc,NB,NZ)*WTLSB1)/WTLSBT
        XFRC=FXRNX*CPOOLD
        NonstructElmnt_brch(ielmc,NB,NZ)=NonstructElmnt_brch(ielmc,NB,NZ)-XFRC
        NoduleNonstructElmnt_brch(ielmc,NB,NZ)=NoduleNonstructElmnt_brch(ielmc,NB,NZ)+XFRC
        CPOOLT=NonstructElmnt_brch(ielmc,NB,NZ)+NoduleNonstructElmnt_brch(ielmc,NB,NZ)
        IF(CPOOLT.GT.ZEROP(NZ))THEN
          ZPOOLD=(NonstructElmnt_brch(ielmn,NB,NZ)*NoduleNonstructElmnt_brch(ielmc,NB,NZ)-NoduleNonstructElmnt_brch(ielmn,NB,NZ)*NonstructElmnt_brch(ielmc,NB,NZ))/CPOOLT
          XFRN=FXRNX*ZPOOLD
          PPOOLD=(NonstructElmnt_brch(ielmp,NB,NZ)*NoduleNonstructElmnt_brch(ielmc,NB,NZ)-NoduleNonstructElmnt_brch(ielmp,NB,NZ)*NonstructElmnt_brch(ielmc,NB,NZ))/CPOOLT
          XFRP=FXRNX*PPOOLD
          NonstructElmnt_brch(ielmn,NB,NZ)=NonstructElmnt_brch(ielmn,NB,NZ)-XFRN
          NonstructElmnt_brch(ielmp,NB,NZ)=NonstructElmnt_brch(ielmp,NB,NZ)-XFRP
          NoduleNonstructElmnt_brch(ielmn,NB,NZ)=NoduleNonstructElmnt_brch(ielmn,NB,NZ)+XFRN
          NoduleNonstructElmnt_brch(ielmp,NB,NZ)=NoduleNonstructElmnt_brch(ielmp,NB,NZ)+XFRP
        ENDIF
      ENDIF
    ENDIF
  ENDIF
  end associate
  end subroutine CanopyNoduleBiochemistry

!------------------------------------------------------------------------------------------

  subroutine RootNoduleBiomchemistry(I,J,NZ,TFN6,WFNGR)
  implicit none
  integer , intent(in) :: I,J,NZ
  real(r8), intent(in) :: TFN6(JZ1)
  real(r8), intent(in) :: WFNGR(2,JZ1)
  integer :: L,M
  real(r8) :: ZADDN
  real(r8) :: ZPOOLD
  real(r8) :: CCC,CNC,CPC
  real(r8) :: CPOOLT
  real(r8) :: CPOOLD
  real(r8) :: CCPOLN,CZPOLN
  real(r8) :: CPPOLN
  real(r8) :: CGNDL,CPOOLNX
  real(r8) :: CCNDLR
  real(r8) :: FCNPF
  real(r8) :: FXRNX
  real(r8) :: GRNDG
  real(r8) :: PPOOLD
  real(r8) :: PADDN
  real(r8) :: RCO2T
  real(r8) :: RCO2TM
  real(r8) :: RCNDL,RMNDL,RXNDL
  real(r8) :: RGNDL,RSNDL
  real(r8) :: RGN2P,RGN2F
  real(r8) :: RXNDLE(NumOfPlantChemElmnts)
  real(r8) :: RDNDLE(NumOfPlantChemElmnts)
  real(r8) :: RCNDLE(NumOfPlantChemElmnts)
  real(r8) :: RGNDG
  real(r8) :: RXNSNE(NumOfPlantChemElmnts)
  real(r8) :: RDNSNE(NumOfPlantChemElmnts)
  real(r8) :: RCNSNE(NumOfPlantChemElmnts)
  real(r8) :: RCNDLM,RXNDLM,RGNDLM
  real(r8) :: RSNDLM
  real(r8) :: SPNDLI
  real(r8) :: SPNDX
  real(r8) :: WTRTD1,WTNDL1,WTRTDT
  real(r8) :: XFRC,XFRN,XFRP
  real(r8) :: RCCC,RCCN,RCCP
  integer  :: NE
!     begin_execution
  associate(                           &
    NU       =>   plt_site%NU        , &
    AREA3    =>   plt_site%AREA3     , &
    ZERO     =>   plt_site%ZERO      , &
    fTgrowRootP     =>   plt_pheno%fTgrowRootP     , &
    DMND     =>   plt_allom%DMND     , &
    CNND     =>   plt_allom%CNND     , &
    CPND     =>   plt_allom%CPND     , &
    k_fine_litr=> pltpar%k_fine_litr , &
    iroot    =>   pltpar%iroot       , &
    RCO2M    =>   plt_rbgc%RCO2M     , &
    RCO2N    =>   plt_rbgc%RCO2N     , &
    WFR      =>   plt_rbgc%WFR       , &
    RCO2A    =>   plt_rbgc%RCO2A     , &
    LitterFallChemElmnt_pftvr     =>   plt_bgcr%LitterFallChemElmnt_pftvr      , &
    UPNF     =>   plt_rbgc%UPNF      , &
    RUPNF    =>   plt_bgcr%RUPNF     , &
     PopuPlantRootC_vr   =>   plt_biom% PopuPlantRootC_vr    , &
    WTNDLE   =>   plt_biom%WTNDLE    , &
    ZEROP    =>   plt_biom%ZEROP     , &
    RootNoduleNonstructElmnt_vr  =>   plt_biom%RootNoduleNonstructElmnt_vr   , &
    ZEROL    =>   plt_biom%ZEROL     , &
     RootMycoNonstructElmnt_vr   =>   plt_biom% RootMycoNonstructElmnt_vr    , &
    CFOPE    =>   plt_soilchem%CFOPE , &
    iPlantNfixType   =>   plt_morph%iPlantNfixType   , &
    NIXBotRootLayer     =>   plt_morph%NIXBotRootLayer       &
  )
!     iPlantNfixType=N2 fixation: 1,2,3=rapid to slow root symbiosis
!     WTNDL,WTNDLN,WTNDLP=bacterial C,N,P mass
!     WTNDI=initial bacterial mass at infection
!     AREA=grid cell area
!     CNND,CPND=bacterial N:C,P:C ratio from PFT file
!
  IF(iPlantNfixType(NZ).GE.1.AND.iPlantNfixType(NZ).LE.3)THEN
    D5400: DO L=NU,NIXBotRootLayer(NZ)
      IF( PopuPlantRootC_vr(ipltroot,L,NZ).GT.ZEROL(NZ))THEN
!
!     INITIAL INFECTION
!
        IF(WTNDLE(ielmc,L,NZ).LE.0.0)THEN
          WTNDLE(ielmc,L,NZ)=WTNDLE(ielmc,L,NZ)+WTNDI*AREA3(NU)
          WTNDLE(ielmn,L,NZ)=WTNDLE(ielmn,L,NZ)+WTNDI*AREA3(NU)*CNND(NZ)
          WTNDLE(ielmp,L,NZ)=WTNDLE(ielmp,L,NZ)+WTNDI*AREA3(NU)*CPND(NZ)
        ENDIF
!
!     O2-UNCONSTRAINED RESPIRATION RATES BY HETEROTROPHIC AEROBES
!     IN NODULE FROM SPECIFIC OXIDATION RATE, ACTIVE BIOMASS,
!     NON-STRUCTURAL C CONCENTRATION, MICROBIAL C:N:P FACTOR,
!     AND TEMPERATURE
!
!     WTNDL,WTNDLN,WTNDLP=bacterial C,N,P mass
!     CPOOLN,ZPOOLN,PPOOLN=nonstructural C,N,P in bacteria
!     CCPOLN,CZPOLN,CPPOLN=nonstructural C,N,P concn in bacteria
!     CNKI,CPKI=nonstructural N,P inhibition constant on growth
!     FCNPF=N,P constraint to bacterial activity
!
        IF(WTNDLE(ielmc,L,NZ).GT.ZEROP(NZ))THEN
          CCPOLN=AZMAX1(RootNoduleNonstructElmnt_vr(ielmc,L,NZ)/WTNDLE(ielmc,L,NZ))
          CZPOLN=AZMAX1(RootNoduleNonstructElmnt_vr(ielmn,L,NZ)/WTNDLE(ielmc,L,NZ))
          CPPOLN=AZMAX1(RootNoduleNonstructElmnt_vr(ielmp,L,NZ)/WTNDLE(ielmc,L,NZ))
        ELSE
          CCPOLN=1.0_r8
          CZPOLN=1.0_r8
          CPPOLN=1.0_r8
        ENDIF
        IF(CCPOLN.GT.ZERO)THEN
          CCC=AZMAX1(AMIN1(1.0,safe_adb(CZPOLN,CZPOLN+CCPOLN*CNKI) &
            ,safe_adb(CPPOLN,CPPOLN+CCPOLN*CPKI)))
!          if(curday==73)write(*,*)CCPOLN,CCPOLN,CZPOLN,CNKI
          CNC=AZMAX1(AMIN1(1.0,safe_adb(CCPOLN,CCPOLN+CZPOLN/CNKI)))
          CPC=AZMAX1(AMIN1(1.0,safe_adb(CCPOLN,CCPOLN+CPPOLN/CPKI)))
        ELSE
          CCC=0._r8
          CNC=0._r8
          CPC=0._r8
        ENDIF
        IF(WTNDLE(ielmc,L,NZ).GT.ZEROP(NZ))THEN
          FCNPF=AMIN1(1.0 &
            ,SQRT(WTNDLE(ielmn,L,NZ)/(WTNDLE(ielmc,L,NZ)*CNND(NZ))) &
            ,SQRT(WTNDLE(ielmp,L,NZ)/(WTNDLE(ielmc,L,NZ)*CPND(NZ))))
        ELSE
          FCNPF=1.0_r8
        ENDIF
        SPNDLI=CCPOLN/(CCPOLN+SPNDLK)
!
!     RESPIRATION FROM NON-STRUCTURAL C DETERMINED BY TEMPERATURE,
!     NON-STRUCTURAL C:N:P
!
!     RCNDLM=respiration from non-structural C unltd by O2
!     CPOOLN,ZPOOLN,PPOOLN=nonstructural C,N,P in bacteria
!     VMXO=specific respiration rate by bacterial N2 fixers
!     WTNDL=bacterial C mass
!     fTgrowRootP=temperature function for root growth
!     FCNPF=N,P constraint to bacterial activity
!     WFNGR=growth function of root water potential
!
        RCNDLM=AZMAX1(AMIN1(RootNoduleNonstructElmnt_vr(ielmc,L,NZ) &
          ,VMXO*WTNDLE(ielmc,L,NZ))*FCNPF*fTgrowRootP(L,NZ)*WFNGR(1,L))
        CPOOLNX=RootNoduleNonstructElmnt_vr(ielmc,L,NZ)
!
!     O2-LIMITED NODULE RESPIRATION FROM 'WFR' IN 'UPTAKE'
!
!     RCNDL=respiration from non-structural C ltd by O2
!     WFR=constraint by O2 consumption on all root processes
!
        RCNDL=RCNDLM*WFR(1,L,NZ)
!
!     NODULE MAINTENANCE RESPIRATION FROM SOIL TEMPERATURE,
!     NODULE STRUCTURAL N
!
!     RMNDL=bacterial maintenance respiration
!     RMPLT=specific maintenance respiration rate (g C g-1 N h-1)
!     TFN6=temperature function for root maintenance respiration
!     WTNDLN=bacterial N mass
!
        RMNDL=AZMAX1(RMPLT*TFN6(L)*WTNDLE(ielmn,L,NZ))*SPNDLI
!
!     NODULE GROWTH RESPIRATION FROM TOTAL - MAINTENANCE
!     IF > 0 DRIVES GROWTH, IF < 0 DRIVES REMOBILIZATION
!
!     RXNDLM,RXNDL=difference between non-structural C respn and mntc respn unltd,ltd by O2
!     RGNDLM,RGNDL=growth respiration unlimited by N,P and unltd,ltd by O2
!     RSNDLM,RSNDL=excess maintenance respiration unltd,ltd by O2
!
        RXNDLM=RCNDLM-RMNDL
        RXNDL=RCNDL-RMNDL
        RGNDLM=AZMAX1(RXNDLM)
        RGNDL=AZMAX1(RXNDL)
        RSNDLM=AZMAX1(-RXNDLM)
        RSNDL=AZMAX1(-RXNDL)
!
!     NODULE N2 FIXATION FROM GROWTH RESPIRATION, FIXATION ENERGY
!     REQUIREMENT AND NON-STRUCTURAL C:N:P PRODUCT INHIBITION,
!     CONSTRAINED BY MICROBIAL N REQUIREMENT
!
!     RGN2P=respiration requirement to maintain bacterial N:C ratio
!     WTNDL,WTNDLN=bacterial C,N mass
!     CNND=bacterial N:C ratio from PFT file
!     EN2F=N fixation yield from C oxidation (g N g-1 C)
!     RGNDL=growth respiration unlimited by N,P
!     RGN2F=respiration for N2 fixation
!     RUPNF,UPNF=layer,total root N2 fixation
!
        RGN2P=AZMAX1(WTNDLE(ielmc,L,NZ)*CNND(NZ)-WTNDLE(ielmn,L,NZ))/EN2F
        IF(RGNDL.GT.ZEROP(NZ))THEN
          RGN2F=RGNDL*RGN2P/(RGNDL+RGN2P)
        ELSE
          RGN2F=0._r8
        ENDIF
        RUPNF(L,NZ)=RGN2F*EN2F
        UPNF(NZ)=UPNF(NZ)+RUPNF(L,NZ)
!
!     NODULE C,N,P REMOBILIZATION AND DECOMPOSITION
!
!     RCCC,RCCN,RCCP=remobilization coefficient for C,N,P
!     RCCZN,RCCYN=min,max fractions for bacteria C recycling
!     RCCXN,RCCQN=max fractions for bacteria N,P recycling
!     WTRTD=root C mass
!     CCNDLR=bacteria:root ratio
!     RDNDLX=effect of CCNDLR on bacteria decomposition rate
!     CCNKR=Km for bacterial vs root mass in decomposition
!     SPNDX=specific bacterial decomposition rate at current CCNDLR
!     WTNDL,WTNDLN,WTNDLP=bacterial C,N,P mass
!     RXNDLE(ielmc),RXNDLE(ielmn),RXNDLE(ielmp)=bacterial C,N,P loss from decomposition
!     RDNDLE(ielmc),RDNDLE(ielmn),RDNDLE(ielmp)=bacterial C,N,P decomposition to litterfall
!     RCNDLE(ielmc),RCNDLE(ielmn),RCNDLE(ielmp)=bacterial C,N,P decomposition to recycling
!
        RCCC=RCCZN+CCC*RCCYN
        RCCN=CNC*RCCXN
        RCCP=CPC*RCCQN
        SPNDX=SPNDL*SQRT(fTgrowRootP(L,NZ)*WFNGR(1,L))
        DO NE=1,NumOfPlantChemElmnts
          RXNDLE(NE)=SPNDX*WTNDLE(NE,L,NZ)
        ENDDO

        RDNDLE(ielmc)=RXNDLE(ielmc)*(1.0_r8-RCCC)
        RDNDLE(ielmn)=RXNDLE(ielmn)*(1.0_r8-RCCC)*(1.0_r8-RCCN)
        RDNDLE(ielmp)=RXNDLE(ielmp)*(1.0_r8-RCCC)*(1.0_r8-RCCP)
        RCNDLE(ielmc)=RXNDLE(ielmc)-RDNDLE(ielmc)
        RCNDLE(ielmn)=RXNDLE(ielmn)-RDNDLE(ielmn)
        RCNDLE(ielmp)=RXNDLE(ielmp)-RDNDLE(ielmp)
!
!     TOTAL NON-STRUCTURAL C,N,P USED IN NODULE GROWTH
!     AND GROWTH RESPIRATION DEPENDS ON GROWTH YIELD
!     ENTERED IN 'READQ'
!
!     CGNDL=total non-structural C used in bacterial growth and growth respiration
!     CPOOLN,ZPOOLN,PPOOLN=nonstructural C,N,P in bacteria
!     RMNDL=bacterial maintenance respiration
!     RCNDL=respiration from non-structural C
!     RCNDLE(ielmc)=bacterial C decomposition to recycling
!     RGNDL=growth respiration ltd by O2
!     RGN2F=respiration for N2 fixation
!     GRNDG=bacterial growth
!     DMND=bacterial growth yield
!     RGNDG=bacterial respiration for growth and N2 fixation
!     ZADDN,PADDN=nonstructural N,P used in growth
!     CNND,CPND=bacterial N:C,P:C ratio from PFT file
!     CCPOLN,CZPOLN,CPPOLN=nonstructural C,N,P concn in bacteria
!     CZKM,CPKM=Km for nonstructural N,P uptake by bacteria
!
        CGNDL=AMIN1(RootNoduleNonstructElmnt_vr(ielmc,L,NZ)-AMIN1(RMNDL,RCNDL) &
          -RGN2F+RCNDLE(ielmc),(RGNDL-RGN2F)/(1.0_r8-DMND(NZ)))
        GRNDG=CGNDL*DMND(NZ)
        RGNDG=RGN2F+CGNDL*(1.0_r8-DMND(NZ))
        ZADDN=AZMAX1(AMIN1(RootNoduleNonstructElmnt_vr(ielmn,L,NZ),GRNDG*CNND(NZ)))*CZPOLN/(CZPOLN+CZKM)
        PADDN=AZMAX1(AMIN1(RootNoduleNonstructElmnt_vr(ielmp,L,NZ),GRNDG*CPND(NZ)))*CPPOLN/(CPPOLN+CPKM)
!
!     NODULE SENESCENCE
!
!     RSNDL=excess maintenance respiration
!     WTNDL,WTNDLN,WTNDLP=bacterial C,N,P mass
!     RCCC,RCCN,RCCP=remobilization coefficient for C,N,P
!     RXNSNE(ielmc),RXNSNE(ielmc),RXNSNE(ielmp)=bacterial C,N,P loss from senescence
!     RDNSNE(ielmc),RDNSNE(ielmc),RDNSNE(ielmp)=bacterial C,N,P senescence to litterfall
!     RCNSNE(ielmc),RCNSNE(ielmc),RCNSNE(ielmp)=bacterial C,N,P senescence to recycling
!
        IF(RSNDL.GT.0.0.AND.WTNDLE(ielmc,L,NZ).GT.ZEROP(NZ).AND.RCCC.GT.ZERO)THEN
          RXNSNE(ielmc)=RSNDL/RCCC
          RXNSNE(ielmn)=RXNSNE(ielmc)*WTNDLE(ielmn,L,NZ)/WTNDLE(ielmc,L,NZ)
          RXNSNE(ielmp)=RXNSNE(ielmc)*WTNDLE(ielmp,L,NZ)/WTNDLE(ielmc,L,NZ)
          RDNSNE(ielmc)=RXNSNE(ielmc)*(1.0_r8-RCCC)
          RDNSNE(ielmn)=RXNSNE(ielmn)*(1.0_r8-RCCC)*(1.0_r8-RCCN)
          RDNSNE(ielmp)=RXNSNE(ielmp)*(1.0_r8-RCCC)*(1.0_r8-RCCP)
          DO NE=1,NumOfPlantChemElmnts
            RCNSNE(NE)=RXNSNE(NE)-RDNSNE(NE)
          ENDDO
        ELSE
          RXNSNE(1:NumOfPlantChemElmnts)=0._r8
          RDNSNE(1:NumOfPlantChemElmnts)=0._r8
          RCNSNE(1:NumOfPlantChemElmnts)=0._r8
        ENDIF
!
!     TOTAL NODULE RESPIRATION
!
!     RCO2TM,RCO2T=total C respiration unlimited,limited by O2
!     GrossResp_pft,TCO2A=total,above-ground PFT respiration
!     RMNDL=bacterial maintenance respiration
!     RCNDL=respiration from non-structural C
!     RGNDG=bacterial respiration for growth and N2 fixation
!     RCNSNE(ielmc)=bacterial C senescence to recycling
!     RCO2A=total root respiration
!     RCO2M,RCO2N,RCO2A unlimited by O2,nonstructural C
!
        RCO2TM=AMIN1(RMNDL,RCNDLM)+RGNDLM+RCNSNE(ielmc)
        RCO2T=AMIN1(RMNDL,RCNDL)+RGNDG+RCNSNE(ielmc)
        RCO2M(ipltroot,L,NZ)=RCO2M(ipltroot,L,NZ)+RCO2TM
        RCO2N(ipltroot,L,NZ)=RCO2N(ipltroot,L,NZ)+RCO2T
        RCO2A(ipltroot,L,NZ)=RCO2A(ipltroot,L,NZ)-RCO2T
!
!     NODULE LITTERFALL CAUSED BY REMOBILIZATION
!
!     CSNC,ZSNC,PSNC=C,N,P litterfall from decomposition and senescence
!     CFOPC,CFOPN,CFOPC=fraction of litterfall C,N,P allocated to litter components
!     RDNDLE(ielmc),RDNDLE(ielmn),RDNDLE(ielmp)=bacterial C,N,P decomposition to litterfall
!     RDNSNE(ielmc),RDNSNE(ielmc),RDNSNE(ielmp)=bacterial C,N,P senescence to litterfall
!
        D6370: DO M=1,jsken
          DO NE=1,NumOfPlantChemElmnts
            LitterFallChemElmnt_pftvr(NE,M,k_fine_litr,L,NZ)=LitterFallChemElmnt_pftvr(NE,M,k_fine_litr,L,NZ)&
              +CFOPE(NE,iroot,M,NZ)*(RDNDLE(NE)+RDNSNE(NE))
          ENDDO
        ENDDO D6370
!
!     CONSUMPTION OF NON-STRUCTURAL C,N,P BY NODULE
!
!     CPOOLN,ZPOOLN,PPOOLN=nonstructural C,N,P in bacteria
!     RMNDL=bacterial maintenance respiration
!     RCNDL=respiration from non-structural C
!     RGN2F=respiration for N2 fixation
!     CGNDL=total non-structural C used in bacterial growth and growth respiration
!     RCNDLE(ielmc),RCNDLE(ielmn),RCNDLE(ielmp)=bacterial C,N,P decomposition to recycling
!     RCNSNE(ielmc),RCNSNE(ielmc),RCNSNE(ielmp)=bacterial C,N,P senescence to recycling
!     ZADDN,PADDN=nonstructural N,P used in growth
!     RUPNF=root N2 fixation
!
        RootNoduleNonstructElmnt_vr(ielmc,L,NZ)=RootNoduleNonstructElmnt_vr(ielmc,L,NZ)-AMIN1(RMNDL,RCNDL)-RGN2F-CGNDL+RCNDLE(ielmc)
        RootNoduleNonstructElmnt_vr(ielmn,L,NZ)=RootNoduleNonstructElmnt_vr(ielmn,L,NZ)-ZADDN+RCNDLE(ielmn)+RCNSNE(ielmn)+RUPNF(L,NZ)
        RootNoduleNonstructElmnt_vr(ielmp,L,NZ)=RootNoduleNonstructElmnt_vr(ielmp,L,NZ)-PADDN+RCNDLE(ielmp)+RCNSNE(ielmp)
!
!     UPDATE STATE VARIABLES FOR NODULE C, N, P
!
!     WTNDL,WTNDLN,WTNDLP=bacterial C,N,P mass
!     GRNDG=bacterial growth
!     RXNDLE(ielmc),RXNDLE(ielmn),RXNDLE(ielmp)=bacterial C,N,P loss from decomposition
!     RXNSNE(ielmc),RXNSNE(ielmc),RXNSNE(ielmp)=bacterial C,N,P loss from senescence
!     ZADDN,PADDN=nonstructural N,P used in growth
!
        WTNDLE(ielmc,L,NZ)=WTNDLE(ielmc,L,NZ)+GRNDG-RXNDLE(ielmc)-RXNSNE(ielmc)
        WTNDLE(ielmn,L,NZ)=WTNDLE(ielmn,L,NZ)+ZADDN-RXNDLE(ielmn)-RXNSNE(ielmn)
        WTNDLE(ielmp,L,NZ)=WTNDLE(ielmp,L,NZ)+PADDN-RXNDLE(ielmp)-RXNSNE(ielmp)
!
!     TRANSFER NON-STRUCTURAL C,N,P BETWEEN ROOT AND NODULES
!     FROM NON-STRUCTURAL C,N,P CONCENTRATION DIFFERENCES
!
!     CPOOLR,ZPOOLR,PPOOLR=root non-structural C,N,P mass
!     WTRTD=root C mass
!     WTNDL=bacterial C mass
!     WTNDI=initial bacterial mass at infection
!     FXRN=rate constant for plant-bacteria nonstructural C,N,P exchange
!     CCNGR=parameter to calculate nonstructural C,N,P exchange
!     CCNDLR=bacteria:root ratio
!     XFRC,XFRN,XFRC=nonstructural C,N,P transfer
!     CPOOLN,ZPOOLN,PPOOLN=nonstructural C,N,P in bacteria
!
        IF( RootMycoNonstructElmnt_vr(ielmc,ipltroot,L,NZ).GT.ZEROP(NZ) &
          .AND. PopuPlantRootC_vr(ipltroot,L,NZ).GT.ZEROL(NZ))THEN
          CCNDLR=WTNDLE(ielmc,L,NZ)/ PopuPlantRootC_vr(ipltroot,L,NZ)
          WTRTD1= PopuPlantRootC_vr(ipltroot,L,NZ)
          WTNDL1=AMIN1( PopuPlantRootC_vr(ipltroot,L,NZ),AMAX1(WTNDI*AREA3(NU),WTNDLE(ielmc,L,NZ)))
          WTRTDT=WTRTD1+WTNDL1
          IF(WTRTDT.GT.ZEROP(NZ))THEN
            FXRNX=FXRN(iPlantNfixType(NZ))/(1.0_r8+CCNDLR/CCNGR)
!    2/(1.0+CCNDLR/(CCNGR*FXRN(iPlantNfixType(NZ))))
            CPOOLD=( RootMycoNonstructElmnt_vr(ielmc,ipltroot,L,NZ)*WTNDL1-RootNoduleNonstructElmnt_vr(ielmc,L,NZ)*WTRTD1)/WTRTDT
            XFRC=FXRNX*CPOOLD
             RootMycoNonstructElmnt_vr(ielmc,ipltroot,L,NZ)= RootMycoNonstructElmnt_vr(ielmc,ipltroot,L,NZ)-XFRC
            RootNoduleNonstructElmnt_vr(ielmc,L,NZ)=RootNoduleNonstructElmnt_vr(ielmc,L,NZ)+XFRC
            CPOOLT= RootMycoNonstructElmnt_vr(ielmc,ipltroot,L,NZ)+RootNoduleNonstructElmnt_vr(ielmc,L,NZ)
            IF(CPOOLT.GT.ZEROP(NZ))THEN
              ZPOOLD=( RootMycoNonstructElmnt_vr(ielmn,ipltroot,L,NZ)*RootNoduleNonstructElmnt_vr(ielmc,L,NZ) &
                -RootNoduleNonstructElmnt_vr(ielmn,L,NZ)* RootMycoNonstructElmnt_vr(ielmc,ipltroot,L,NZ))/CPOOLT
              XFRN=FXRNX*ZPOOLD
              PPOOLD=( RootMycoNonstructElmnt_vr(ielmp,ipltroot,L,NZ)*RootNoduleNonstructElmnt_vr(ielmc,L,NZ) &
                -RootNoduleNonstructElmnt_vr(ielmp,L,NZ)* RootMycoNonstructElmnt_vr(ielmc,ipltroot,L,NZ))/CPOOLT
              XFRP=FXRNX*PPOOLD
               RootMycoNonstructElmnt_vr(ielmn,ipltroot,L,NZ)= RootMycoNonstructElmnt_vr(ielmn,ipltroot,L,NZ)-XFRN
               RootMycoNonstructElmnt_vr(ielmp,ipltroot,L,NZ)= RootMycoNonstructElmnt_vr(ielmp,ipltroot,L,NZ)-XFRP
              RootNoduleNonstructElmnt_vr(ielmn,L,NZ)=RootNoduleNonstructElmnt_vr(ielmn,L,NZ)+XFRN
              RootNoduleNonstructElmnt_vr(ielmp,L,NZ)=RootNoduleNonstructElmnt_vr(ielmp,L,NZ)+XFRP
            ENDIF
          ENDIF
        ENDIF
      ENDIF
    ENDDO D5400
  ENDIF
  end associate
  end subroutine RootNoduleBiomchemistry
end module NoduleBGCMod
