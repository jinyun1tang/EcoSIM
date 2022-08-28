module NoduleBGCMod

  use data_kind_mod, only : r8 => SHR_KIND_R8
  use minimathmod, only : test_aeqb,safe_adb
  use EcosimConst
  use PlantAPIData
  use GrosubPars
  implicit none
  private

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
  real(r8) :: RXNDLC,RXNDLN,RXNDLP
  real(r8) :: RDNDLC,RDNDLN,RDNDLP
  real(r8) :: RCNDLC,RCNDLN,RCNDLP
  real(r8) :: RGNDG
  real(r8) :: RXNSNC,RXNSNN,RXNSNP
  real(r8) :: RDNSNC,RDNSNN,RDNSNP
  real(r8) :: RCNSNC,RCNSNN,RCNSNP
  real(r8) :: SPNDLI
  real(r8) :: SPNDX
  real(r8) :: WTLSB1,WTNDB1,WTLSBT
  real(r8) :: XFRC,XFRN,XFRP
  REAL(R8) :: RCCC,RCCN,RCCP
!     begin_execution
  associate(                             &
    NUs1       =>  plt_site%NUs1       , &
    ZEROs1     =>  plt_site%ZEROs1     , &
    AREA3s1    =>  plt_site%AREA3s1    , &
    CFOPCs1    =>  plt_soilchem%CFOPCs1, &
    CFOPNs1    =>  plt_soilchem%CFOPNs1, &
    CFOPPs1    =>  plt_soilchem%CFOPPs1, &
    INTYPs1    =>  plt_morph%INTYPs1   , &
    TFN3s1     =>  plt_pheno%TFN3s1    , &
    TCO2Ts1    =>  plt_bgcr%TCO2Ts1    , &
    RECOs1     =>  plt_bgcr%RECOs1     , &
    TCO2As1    =>  plt_bgcr%TCO2As1    , &
    TRAUs1     =>  plt_bgcr%TRAUs1     , &
    CNETs1     =>  plt_bgcr%CNETs1     , &
    CSNCs1     =>  plt_bgcr%CSNCs1     , &
    ZSNCs1     =>  plt_bgcr%ZSNCs1     , &
    PSNCs1     =>  plt_bgcr%PSNCs1     , &
    DMNDs1     =>  plt_allom%DMNDs1    , &
    CNNDs1     =>  plt_allom%CNNDs1    , &
    CPNDs1     =>  plt_allom%CPNDs1    , &
    WTLSBs1    =>  plt_biom%WTLSBs1    , &
    WTNDBPs1   =>  plt_biom%WTNDBPs1   , &
    CPOOLs1    =>  plt_biom%CPOOLs1    , &
    PPOOLs1    =>  plt_biom%PPOOLs1    , &
    ZPOOLs1    =>  plt_biom%ZPOOLs1    , &
    CPOLNBs1   =>  plt_biom%CPOLNBs1   , &
    ZPOLNBs1   =>  plt_biom%ZPOLNBs1   , &
    PPOLNBs1   =>  plt_biom%PPOLNBs1   , &
    ZEROPs1    =>  plt_biom%ZEROPs1    , &
    ZEROLs1    =>  plt_biom%ZEROLs1    , &
    WTNDBNs1   =>  plt_biom%WTNDBNs1   , &
    WTNDBs1    =>  plt_biom%WTNDBs1      &
  )
!     INTYP=N2 fixation: 4,5,6=rapid to slow canopy symbiosis
!
  IF(INTYPs1(NZ).GE.4)THEN
!
!     INITIAL INFECTION
!
!     WTNDB,WTNDBN,WTNDBP=bacterial C,N,P mass
!     WTNDI=initial bacterial mass at infection
!     AREA=grid cell area
!     CNND,CPND=bacterial N:C,P:C ratio from PFT file
!
    IF(WTNDBs1(NB,NZ).LE.0.0)THEN
      WTNDBs1(NB,NZ)=WTNDBs1(NB,NZ) &
        +WTNDI*AREA3s1(NUs1)
      WTNDBNs1(NB,NZ)=WTNDBNs1(NB,NZ) &
        +WTNDI*AREA3s1(NUs1)*CNNDs1(NZ)
      WTNDBPs1(NB,NZ)=WTNDBPs1(NB,NZ) &
        +WTNDI*AREA3s1(NUs1)*CPNDs1(NZ)
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
    IF(WTNDBs1(NB,NZ).GT.ZEROPs1(NZ))THEN
      CCPOLN=AMAX1(0.0,CPOLNBs1(NB,NZ)/WTNDBs1(NB,NZ))
      CZPOLN=AMAX1(0.0,ZPOLNBs1(NB,NZ)/WTNDBs1(NB,NZ))
      CPPOLN=AMAX1(0.0,PPOLNBs1(NB,NZ)/WTNDBs1(NB,NZ))
    ELSE
      CCPOLN=1.0_r8
      CZPOLN=1.0_r8
      CPPOLN=1.0_r8
    ENDIF
    IF(CCPOLN.GT.ZEROs1)THEN
      CCC=AMAX1(0.0,AMIN1(1.0,safe_adb(CZPOLN,CZPOLN+CCPOLN*CNKI) &
        ,safe_adb(CPPOLN,CPPOLN+CCPOLN*CPKI)))
      CNC=AMAX1(0.0,AMIN1(1.0,safe_adb(CCPOLN,CCPOLN+CZPOLN/CNKI)))
      CPC=AMAX1(0.0,AMIN1(1.0,safe_adb(CCPOLN,CCPOLN+CPPOLN/CPKI)))
    ELSE
      CCC=0._r8
      CNC=0._r8
      CPC=0._r8
    ENDIF
    IF(WTNDBs1(NB,NZ).GT.ZEROPs1(NZ))THEN
      FCNPF=AMIN1(1.0 &
        ,SQRT(WTNDBNs1(NB,NZ)/(WTNDBs1(NB,NZ)*CNNDs1(NZ))) &
        ,SQRT(WTNDBPs1(NB,NZ)/(WTNDBs1(NB,NZ)*CPNDs1(NZ))))
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
!     TFN3=temperature function for canopy growth
!     FCNPF=N,P constraint to bacterial activity
!     WFNG=growth function of canopy water potential
!
    RCNDL=AMAX1(0.0,AMIN1(CPOLNBs1(NB,NZ) &
      ,VMXO*WTNDBs1(NB,NZ))*FCNPF*TFN3s1(NZ)*WFNG)
!     CPOOLNX=CPOLNBs1(NB,NZ)
!     VMXOX=VMXO*WTNDBs1(NB,NZ)*FCNPF*TFN3s1(NZ)*WFNG
!
!     NODULE MAINTENANCE RESPIRATION FROM SOIL TEMPERATURE,
!     NODULE STRUCTURAL N
!
!     RMNDL=bacterial maintenance respiration
!     RMPLT=specific maintenance respiration rate (g C g-1 N h-1)
!     TFN5=temperature function for canopy maintenance respiration
!     WTNDBN=bacterial N mass
!
    RMNDL=AMAX1(0.0,RMPLT*TFN5*WTNDBNs1(NB,NZ))*SPNDLI
!
!     NODULE GROWTH RESPIRATION FROM TOTAL - MAINTENANCE
!     IF > 0 DRIVES GROWTH, IF < 0 DRIVES REMOBILIZATION
!
!     RXNDL=difference between non-structural C respn and mntc respn
!     RGNDL=growth respiration unlimited by N,P
!     RSNDL=excess maintenance respiration
!
    RXNDL=RCNDL-RMNDL
    RGNDL=AMAX1(0.0,RXNDL)
    RSNDL=AMAX1(0.0,-RXNDL)
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
    RGN2P=AMAX1(0.0,WTNDBs1(NB,NZ)*CNNDs1(NZ)-WTNDBNs1(NB,NZ))/EN2F
    IF(RGNDL.GT.ZEROPs1(NZ))THEN
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
!     WTLSB=leaf+petiole mass
!     CCNDLB=bacteria:leaf+petiole ratio
!     RDNDBX=effect of CCNDLB on bacteria decomposition rate
!     SPNDX=specific bacterial decomposition rate at current CCNDLB
!     SPNDL=specific decomposition rate by bacterial N2 fixers
!     WTNDB,WTNDBN,WTNDBP=bacterial C,N,P mass
!     RXNDLC,RXNDLN,RXNDLP=bacterial C,N,P loss from decomposition
!     RDNDLC,RDNDLN,RDNDLP=bacterial C,N,P decomposition to litterfall
!     RCNDLC,RCNDLN,RCNDLP=bacterial C,N,P decomposition to recycling
!
    RCCC=RCCZN+CCC*RCCYN
    RCCN=CNC*RCCXN
    RCCP=CPC*RCCQN
    SPNDX=SPNDL*SQRT(TFN3s1(NZ)*WFNG)
    RXNDLC=SPNDX*WTNDBs1(NB,NZ)
    RXNDLN=SPNDX*WTNDBNs1(NB,NZ)
    RXNDLP=SPNDX*WTNDBPs1(NB,NZ)
    RDNDLC=RXNDLC*(1.0_r8-RCCC)
    RDNDLN=RXNDLN*(1.0_r8-RCCC)*(1.0_r8-RCCN)
    RDNDLP=RXNDLP*(1.0_r8-RCCC)*(1.0_r8-RCCP)
    RCNDLC=RXNDLC-RDNDLC
    RCNDLN=RXNDLN-RDNDLN
    RCNDLP=RXNDLP-RDNDLP
!
!     TOTAL NON-STRUCTURAL C,N,P USED IN NODULE GROWTH
!     AND GROWTH RESPIRATION DEPENDS ON GROWTH YIELD
!     ENTERED IN 'READQ'
!
!     CGNDL=total non-structural C used in bacterial growth and growth respiration
!     CPOLNB,ZPOLNB,PPOLNB=nonstructural C,N,P in bacteria
!     RMNDL=bacterial maintenance respiration
!     RCNDL=respiration from non-structural C
!     RCNDLC=bacterial C decomposition to recycling
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
    CGNDL=AMIN1(CPOLNBs1(NB,NZ)-AMIN1(RMNDL,RCNDL) &
      -RGN2F+RCNDLC,(RGNDL-RGN2F)/(1.0_r8-DMNDs1(NZ)))
    GRNDG=CGNDL*DMNDs1(NZ)
    RGNDG=RGN2F+CGNDL*(1.0_r8-DMNDs1(NZ))
    ZADDN=AMAX1(0.0,AMIN1(ZPOLNBs1(NB,NZ) &
      ,GRNDG*CNNDs1(NZ)))*CZPOLN/(CZPOLN+CZKM)
    PADDN=AMAX1(0.0,AMIN1(PPOLNBs1(NB,NZ) &
      ,GRNDG*CPNDs1(NZ)))*CPPOLN/(CPPOLN+CPKM)
!
!     NODULE SENESCENCE
!
!     RSNDL=excess maintenance respiration
!     WTNDB,WTNDBN,WTNDBP=bacterial C,N,P mass
!     RCCC,RCCN,RCCP=remobilization coefficient for C,N,P
!     RXNSNC,RXNSNC,RXNSNP=bacterial C,N,P loss from senescence
!     RDNSNC,RDNSNC,RDNSNP=bacterial C,N,P senescence to litterfall
!     RCNSNC,RCNSNC,RCNSNP=bacterial C,N,P senescence to recycling
!
    IF(RSNDL.GT.0.0.AND.WTNDBs1(NB,NZ).GT.ZEROPs1(NZ).AND.RCCC.GT.ZEROs1)THEN
      RXNSNC=RSNDL/RCCC
      RXNSNN=RXNSNC*WTNDBNs1(NB,NZ)/WTNDBs1(NB,NZ)
      RXNSNP=RXNSNC*WTNDBPs1(NB,NZ)/WTNDBs1(NB,NZ)
      RDNSNC=RXNSNC*(1.0_r8-RCCC)
      RDNSNN=RXNSNN*(1.0_r8-RCCC)*(1.0_r8-RCCN)
      RDNSNP=RXNSNP*(1.0_r8-RCCC)*(1.0_r8-RCCP)
      RCNSNC=RXNSNC-RDNSNC
      RCNSNN=RXNSNN-RDNSNN
      RCNSNP=RXNSNP-RDNSNP
    ELSE
      RXNSNC=0._r8
      RXNSNN=0._r8
      RXNSNP=0._r8
      RDNSNC=0._r8
      RDNSNN=0._r8
      RDNSNP=0._r8
      RCNSNC=0._r8
      RCNSNN=0._r8
      RCNSNP=0._r8
    ENDIF
!
!     TOTAL NODULE RESPIRATION
!
!     RCO2T=total C respiration
!     RMNDL=bacterial maintenance respiration
!     RCNDL=respiration from non-structural C
!     RGNDG=bacterial respiration for growth and N2 fixation
!     RCNSNC=bacterial C senescence to recycling
!     TCO2T,TCO2A=total,above-ground PFT respiration
!     CNET=PFT net CO2 fixation
!     RECO=ecosystem respiration
!     TRAU=total autotrophic respiration
!
    RCO2T=AMIN1(RMNDL,RCNDL)+RGNDG+RCNSNC
    TCO2Ts1(NZ)=TCO2Ts1(NZ)-RCO2T
    TCO2As1(NZ)=TCO2As1(NZ)-RCO2T
    CNETs1(NZ)=CNETs1(NZ)-RCO2T
    RECOs1=RECOs1-RCO2T
    TRAUs1=TRAUs1-RCO2T
!
!     NODULE LITTERFALL CAUSED BY REMOBILIZATION
!
!     CSNC,ZSNC,PSNC=C,N,P litterfall from decomposition and senescence
!     CFOPC,CFOPN,CFOPC=fraction of litterfall C,N,P allocated to litter components
!     RDNDLC,RDNDLN,RDNDLP=bacterial C,N,P decomposition to litterfall
!     RDNSNC,RDNSNC,RDNSNP=bacterial C,N,P senescence to litterfall
!
    DO 6470 M=1,jsken
      CSNCs1(M,1,0,NZ)=CSNCs1(M,1,0,NZ)+CFOPCs1(1,M,NZ)*(RDNDLC+RDNSNC)
      ZSNCs1(M,1,0,NZ)=ZSNCs1(M,1,0,NZ)+CFOPNs1(1,M,NZ)*(RDNDLN+RDNSNN)
      PSNCs1(M,1,0,NZ)=PSNCs1(M,1,0,NZ)+CFOPPs1(1,M,NZ)*(RDNDLP+RDNSNP)
6470  CONTINUE
!
!     CONSUMPTION OF NON-STRUCTURAL C,N,P BY NODULE
!
!     CPOLNB,ZPOLNB,PPOLNB=nonstructural C,N,P in bacteria
!     RMNDL=bacterial maintenance respiration
!     RCNDL=respiration from non-structural C
!     RGN2F=respiration for N2 fixation
!     CGNDL=total non-structural C used in bacterial growth and growth respiration
!     RCNDLC,RCNDLN,RCNDLP=bacterial C,N,P decomposition to recycling
!     RCNSNC,RCNSNC,RCNSNP=bacterial C,N,P senescence to recycling
!     ZADDN,PADDN=nonstructural N,P used in growth
!     RUPNFB=branch N2 fixation
!
    CPOLNBs1(NB,NZ)=CPOLNBs1(NB,NZ)-AMIN1(RMNDL,RCNDL)-RGN2F-CGNDL+RCNDLC
    ZPOLNBs1(NB,NZ)=ZPOLNBs1(NB,NZ)-ZADDN+RCNDLN+RCNSNN+RUPNFB
    PPOLNBs1(NB,NZ)=PPOLNBs1(NB,NZ)-PADDN+RCNDLP+RCNSNP
!
!     UPDATE STATE VARIABLES FOR NODULE C, N, P
!
!     WTNDB,WTNDBN,WTNDBP=bacterial C,N,P mass
!     GRNDG=bacterial growth
!     RXNDLC,RXNDLN,RXNDLP=bacterial C,N,P loss from decomposition
!     RXNSNC,RXNSNC,RXNSNP=bacterial C,N,P loss from senescence
!     ZADDN,PADDN=nonstructural N,P used in growth
!
    WTNDBs1(NB,NZ)=WTNDBs1(NB,NZ)+GRNDG-RXNDLC-RXNSNC
    WTNDBNs1(NB,NZ)=WTNDBNs1(NB,NZ)+ZADDN-RXNDLN-RXNSNN
    WTNDBPs1(NB,NZ)=WTNDBPs1(NB,NZ)+PADDN-RXNDLP-RXNSNP
!
!     TRANSFER NON-STRUCTURAL C,N,P BETWEEN BRANCH AND NODULES
!     FROM NON-STRUCTURAL C,N,P CONCENTRATION DIFFERENCES
!
!     CPOOL,ZPOOL,PPOOL=branch non-structural C,N,P mass
!     WTLSB=leaf+petiole C mass
!     WTNDB=bacterial C mass
!     WTNDI=initial bacterial mass at infection
!     FXRN=rate constant for plant-bacteria nonstructural C,N,P exchange
!     CCNGB=parameter to calculate nonstructural C,N,P exchange
!     CCNDLB=bacteria:leaf+petiole ratio
!     XFRC,XFRN,XFRC=nonstructural C,N,P transfer
!     CPOLNB,ZPOLNB,PPOLNB=nonstructural C,N,P in bacteria
!
    IF(CPOOLs1(NB,NZ).GT.ZEROPs1(NZ) &
      .AND.WTLSBs1(NB,NZ).GT.ZEROLs1(NZ))THEN
      CCNDLB=WTNDBs1(NB,NZ)/WTLSBs1(NB,NZ)
      WTLSB1=WTLSBs1(NB,NZ)
      WTNDB1=AMIN1(WTLSBs1(NB,NZ),AMAX1(WTNDI*AREA3s1(NUs1),WTNDBs1(NB,NZ)))
      WTLSBT=WTLSB1+WTNDB1
      IF(WTLSBT.GT.ZEROPs1(NZ))THEN
        FXRNX=FXRN(INTYPs1(NZ))/(1.0+CCNDLB/CCNGB)
!    2/(1.0+CCNDLB/(CCNGB*FXRN(INTYPs1(NZ))))
        CPOOLD=(CPOOLs1(NB,NZ)*WTNDB1-CPOLNBs1(NB,NZ)*WTLSB1)/WTLSBT
        XFRC=FXRNX*CPOOLD
        CPOOLs1(NB,NZ)=CPOOLs1(NB,NZ)-XFRC
        CPOLNBs1(NB,NZ)=CPOLNBs1(NB,NZ)+XFRC
        CPOOLT=CPOOLs1(NB,NZ)+CPOLNBs1(NB,NZ)
        IF(CPOOLT.GT.ZEROPs1(NZ))THEN
          ZPOOLD=(ZPOOLs1(NB,NZ)*CPOLNBs1(NB,NZ)-ZPOLNBs1(NB,NZ)*CPOOLs1(NB,NZ))/CPOOLT
          XFRN=FXRNX*ZPOOLD
          PPOOLD=(PPOOLs1(NB,NZ)*CPOLNBs1(NB,NZ) &
            -PPOLNBs1(NB,NZ)*CPOOLs1(NB,NZ))/CPOOLT
          XFRP=FXRNX*PPOOLD
          ZPOOLs1(NB,NZ)=ZPOOLs1(NB,NZ)-XFRN
          PPOOLs1(NB,NZ)=PPOOLs1(NB,NZ)-XFRP
          ZPOLNBs1(NB,NZ)=ZPOLNBs1(NB,NZ)+XFRN
          PPOLNBs1(NB,NZ)=PPOLNBs1(NB,NZ)+XFRP
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
  real(r8) :: RXNDLC,RXNDLN,RXNDLP
  real(r8) :: RDNDLC,RDNDLN,RDNDLP
  real(r8) :: RCNDLC,RCNDLN,RCNDLP
  real(r8) :: RGNDG
  real(r8) :: RXNSNC,RXNSNN,RXNSNP
  real(r8) :: RDNSNC,RDNSNN,RDNSNP
  real(r8) :: RCNSNC,RCNSNN,RCNSNP
  real(r8) :: RCNDLM,RXNDLM,RGNDLM
  real(r8) :: RSNDLM
  real(r8) :: SPNDLI
  real(r8) :: SPNDX
  real(r8) :: WTRTD1,WTNDL1,WTRTDT
  real(r8) :: XFRC,XFRN,XFRP
  real(r8) :: RCCC,RCCN,RCCP
!     begin_execution
  associate(                               &
    NUs1       =>   plt_site%NUs1        , &
    AREA3s1    =>   plt_site%AREA3s1     , &
    ZEROs1     =>   plt_site%ZEROs1      , &
    TFN4s1     =>   plt_pheno%TFN4s1     , &
    DMNDs1     =>   plt_allom%DMNDs1     , &
    CNNDs1     =>   plt_allom%CNNDs1     , &
    CPNDs1     =>   plt_allom%CPNDs1     , &
    RCO2Ms1    =>   plt_rbgc%RCO2Ms1     , &
    RCO2Ns1    =>   plt_rbgc%RCO2Ns1     , &
    WFRs1      =>   plt_rbgc%WFRs1       , &
    RCO2As1    =>   plt_rbgc%RCO2As1     , &
    CSNCs1     =>   plt_bgcr%CSNCs1      , &
    ZSNCs1     =>   plt_bgcr%ZSNCs1      , &
    PSNCs1     =>   plt_bgcr%PSNCs1      , &
    UPNFs1     =>   plt_rbgc%UPNFs1      , &
    RUPNFs1    =>   plt_bgcr%RUPNFs1     , &
    ZPOOLNs1   =>   plt_biom%ZPOOLNs1    , &
    PPOOLNs1   =>   plt_biom%PPOOLNs1    , &
    WTRTDs1    =>   plt_biom%WTRTDs1     , &
    WTNDLNs1   =>   plt_biom%WTNDLNs1    , &
    WTNDLs1    =>   plt_biom%WTNDLs1     , &
    WTNDLPs1   =>   plt_biom%WTNDLPs1    , &
    ZEROPs1    =>   plt_biom%ZEROPs1     , &
    CPOOLNs1   =>   plt_biom%CPOOLNs1    , &
    ZEROLs1    =>   plt_biom%ZEROLs1     , &
    CPOOLRs1   =>   plt_biom%CPOOLRs1    , &
    ZPOOLRs1   =>   plt_biom%ZPOOLRs1    , &
    PPOOLRs1   =>   plt_biom%PPOOLRs1    , &
    CFOPCs1    =>   plt_soilchem%CFOPCs1 , &
    CFOPNs1    =>   plt_soilchem%CFOPNs1 , &
    CFOPPs1    =>   plt_soilchem%CFOPPs1 , &
    INTYPs1    =>   plt_morph%INTYPs1    , &
    NIXs1      =>   plt_morph%NIXs1        &
  )
!     INTYP=N2 fixation: 1,2,3=rapid to slow root symbiosis
!     WTNDL,WTNDLN,WTNDLP=bacterial C,N,P mass
!     WTNDI=initial bacterial mass at infection
!     AREA=grid cell area
!     CNND,CPND=bacterial N:C,P:C ratio from PFT file
!
  IF(INTYPs1(NZ).GE.1.AND.INTYPs1(NZ).LE.3)THEN
    DO 5400 L=NUs1,NIXs1(NZ)
      IF(WTRTDs1(1,L,NZ).GT.ZEROLs1(NZ))THEN
!
!     INITIAL INFECTION
!
        IF(WTNDLs1(L,NZ).LE.0.0)THEN
          WTNDLs1(L,NZ)=WTNDLs1(L,NZ)+WTNDI*AREA3s1(NUs1)
          WTNDLNs1(L,NZ)=WTNDLNs1(L,NZ)+WTNDI*AREA3s1(NUs1)*CNNDs1(NZ)
          WTNDLPs1(L,NZ)=WTNDLPs1(L,NZ)+WTNDI*AREA3s1(NUs1)*CPNDs1(NZ)
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
        IF(WTNDLs1(L,NZ).GT.ZEROPs1(NZ))THEN
          CCPOLN=AMAX1(0.0,CPOOLNs1(L,NZ)/WTNDLs1(L,NZ))
          CZPOLN=AMAX1(0.0,ZPOOLNs1(L,NZ)/WTNDLs1(L,NZ))
          CPPOLN=AMAX1(0.0,PPOOLNs1(L,NZ)/WTNDLs1(L,NZ))
        ELSE
          CCPOLN=1.0_r8
          CZPOLN=1.0_r8
          CPPOLN=1.0_r8
        ENDIF
        IF(CCPOLN.GT.ZEROs1)THEN
          CCC=AMAX1(0.0,AMIN1(1.0,safe_adb(CZPOLN,CZPOLN+CCPOLN*CNKI) &
            ,safe_adb(CPPOLN,CPPOLN+CCPOLN*CPKI)))
!          if(curday==73)write(*,*)CCPOLN,CCPOLN,CZPOLN,CNKI
          CNC=AMAX1(0.0,AMIN1(1.0,safe_adb(CCPOLN,CCPOLN+CZPOLN/CNKI)))
          CPC=AMAX1(0.0,AMIN1(1.0,safe_adb(CCPOLN,CCPOLN+CPPOLN/CPKI)))
        ELSE
          CCC=0._r8
          CNC=0._r8
          CPC=0._r8
        ENDIF
        IF(WTNDLs1(L,NZ).GT.ZEROPs1(NZ))THEN
          FCNPF=AMIN1(1.0 &
            ,SQRT(WTNDLNs1(L,NZ)/(WTNDLs1(L,NZ)*CNNDs1(NZ))) &
            ,SQRT(WTNDLPs1(L,NZ)/(WTNDLs1(L,NZ)*CPNDs1(NZ))))
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
!     TFN4=temperature function for root growth
!     FCNPF=N,P constraint to bacterial activity
!     WFNGR=growth function of root water potential
!
        RCNDLM=AMAX1(0.0,AMIN1(CPOOLNs1(L,NZ) &
          ,VMXO*WTNDLs1(L,NZ))*FCNPF*TFN4s1(L,NZ)*WFNGR(1,L))
        CPOOLNX=CPOOLNs1(L,NZ)
!
!     O2-LIMITED NODULE RESPIRATION FROM 'WFR' IN 'UPTAKE'
!
!     RCNDL=respiration from non-structural C ltd by O2
!     WFR=constraint by O2 consumption on all root processes
!
        RCNDL=RCNDLM*WFRs1(1,L,NZ)
!
!     NODULE MAINTENANCE RESPIRATION FROM SOIL TEMPERATURE,
!     NODULE STRUCTURAL N
!
!     RMNDL=bacterial maintenance respiration
!     RMPLT=specific maintenance respiration rate (g C g-1 N h-1)
!     TFN6=temperature function for root maintenance respiration
!     WTNDLN=bacterial N mass
!
        RMNDL=AMAX1(0.0,RMPLT*TFN6(L)*WTNDLNs1(L,NZ))*SPNDLI
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
        RGNDLM=AMAX1(0.0,RXNDLM)
        RGNDL=AMAX1(0.0,RXNDL)
        RSNDLM=AMAX1(0.0,-RXNDLM)
        RSNDL=AMAX1(0.0,-RXNDL)
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
        RGN2P=AMAX1(0.0,WTNDLs1(L,NZ)*CNNDs1(NZ)-WTNDLNs1(L,NZ))/EN2F
        IF(RGNDL.GT.ZEROPs1(NZ))THEN
          RGN2F=RGNDL*RGN2P/(RGNDL+RGN2P)
        ELSE
          RGN2F=0._r8
        ENDIF
        RUPNFs1(L,NZ)=RGN2F*EN2F
        UPNFs1(NZ)=UPNFs1(NZ)+RUPNFs1(L,NZ)
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
!     RXNDLC,RXNDLN,RXNDLP=bacterial C,N,P loss from decomposition
!     RDNDLC,RDNDLN,RDNDLP=bacterial C,N,P decomposition to litterfall
!     RCNDLC,RCNDLN,RCNDLP=bacterial C,N,P decomposition to recycling
!
        RCCC=RCCZN+CCC*RCCYN
        RCCN=CNC*RCCXN
        RCCP=CPC*RCCQN
        SPNDX=SPNDL*SQRT(TFN4s1(L,NZ)*WFNGR(1,L))
        RXNDLC=SPNDX*WTNDLs1(L,NZ)
        RXNDLN=SPNDX*WTNDLNs1(L,NZ)
        RXNDLP=SPNDX*WTNDLPs1(L,NZ)
        RDNDLC=RXNDLC*(1.0_r8-RCCC)
        RDNDLN=RXNDLN*(1.0_r8-RCCC)*(1.0_r8-RCCN)
        RDNDLP=RXNDLP*(1.0_r8-RCCC)*(1.0_r8-RCCP)
        RCNDLC=RXNDLC-RDNDLC
        RCNDLN=RXNDLN-RDNDLN
        RCNDLP=RXNDLP-RDNDLP
!
!     TOTAL NON-STRUCTURAL C,N,P USED IN NODULE GROWTH
!     AND GROWTH RESPIRATION DEPENDS ON GROWTH YIELD
!     ENTERED IN 'READQ'
!
!     CGNDL=total non-structural C used in bacterial growth and growth respiration
!     CPOOLN,ZPOOLN,PPOOLN=nonstructural C,N,P in bacteria
!     RMNDL=bacterial maintenance respiration
!     RCNDL=respiration from non-structural C
!     RCNDLC=bacterial C decomposition to recycling
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
        CGNDL=AMIN1(CPOOLNs1(L,NZ)-AMIN1(RMNDL,RCNDL) &
          -RGN2F+RCNDLC,(RGNDL-RGN2F)/(1.0_r8-DMNDs1(NZ)))
        GRNDG=CGNDL*DMNDs1(NZ)
        RGNDG=RGN2F+CGNDL*(1.0_r8-DMNDs1(NZ))
        ZADDN=AMAX1(0.0,AMIN1(ZPOOLNs1(L,NZ),GRNDG*CNNDs1(NZ)))*CZPOLN/(CZPOLN+CZKM)
        PADDN=AMAX1(0.0,AMIN1(PPOOLNs1(L,NZ),GRNDG*CPNDs1(NZ)))*CPPOLN/(CPPOLN+CPKM)
!
!     NODULE SENESCENCE
!
!     RSNDL=excess maintenance respiration
!     WTNDL,WTNDLN,WTNDLP=bacterial C,N,P mass
!     RCCC,RCCN,RCCP=remobilization coefficient for C,N,P
!     RXNSNC,RXNSNC,RXNSNP=bacterial C,N,P loss from senescence
!     RDNSNC,RDNSNC,RDNSNP=bacterial C,N,P senescence to litterfall
!     RCNSNC,RCNSNC,RCNSNP=bacterial C,N,P senescence to recycling
!
        IF(RSNDL.GT.0.0.AND.WTNDLs1(L,NZ).GT.ZEROPs1(NZ).AND.RCCC.GT.ZEROs1)THEN
          RXNSNC=RSNDL/RCCC
          RXNSNN=RXNSNC*WTNDLNs1(L,NZ)/WTNDLs1(L,NZ)
          RXNSNP=RXNSNC*WTNDLPs1(L,NZ)/WTNDLs1(L,NZ)
          RDNSNC=RXNSNC*(1.0_r8-RCCC)
          RDNSNN=RXNSNN*(1.0_r8-RCCC)*(1.0_r8-RCCN)
          RDNSNP=RXNSNP*(1.0_r8-RCCC)*(1.0_r8-RCCP)
          RCNSNC=RXNSNC-RDNSNC
          RCNSNN=RXNSNN-RDNSNN
          RCNSNP=RXNSNP-RDNSNP
        ELSE
          RXNSNC=0._r8
          RXNSNN=0._r8
          RXNSNP=0._r8
          RDNSNC=0._r8
          RDNSNN=0._r8
          RDNSNP=0._r8
          RCNSNC=0._r8
          RCNSNN=0._r8
          RCNSNP=0._r8
        ENDIF
!
!     TOTAL NODULE RESPIRATION
!
!     RCO2TM,RCO2T=total C respiration unlimited,limited by O2
!     TCO2T,TCO2A=total,above-ground PFT respiration
!     RMNDL=bacterial maintenance respiration
!     RCNDL=respiration from non-structural C
!     RGNDG=bacterial respiration for growth and N2 fixation
!     RCNSNC=bacterial C senescence to recycling
!     RCO2A=total root respiration
!     RCO2M,RCO2N,RCO2A unlimited by O2,nonstructural C
!
        RCO2TM=AMIN1(RMNDL,RCNDLM)+RGNDLM+RCNSNC
        RCO2T=AMIN1(RMNDL,RCNDL)+RGNDG+RCNSNC
        RCO2Ms1(1,L,NZ)=RCO2Ms1(1,L,NZ)+RCO2TM
        RCO2Ns1(1,L,NZ)=RCO2Ns1(1,L,NZ)+RCO2T
        RCO2As1(1,L,NZ)=RCO2As1(1,L,NZ)-RCO2T
!
!     NODULE LITTERFALL CAUSED BY REMOBILIZATION
!
!     CSNC,ZSNC,PSNC=C,N,P litterfall from decomposition and senescence
!     CFOPC,CFOPN,CFOPC=fraction of litterfall C,N,P allocated to litter components
!     RDNDLC,RDNDLN,RDNDLP=bacterial C,N,P decomposition to litterfall
!     RDNSNC,RDNSNC,RDNSNP=bacterial C,N,P senescence to litterfall
!
        DO 6370 M=1,jsken
          CSNCs1(M,1,L,NZ)=CSNCs1(M,1,L,NZ)+CFOPCs1(4,M,NZ)*(RDNDLC+RDNSNC)
          ZSNCs1(M,1,L,NZ)=ZSNCs1(M,1,L,NZ)+CFOPNs1(4,M,NZ)*(RDNDLN+RDNSNN)
          PSNCs1(M,1,L,NZ)=PSNCs1(M,1,L,NZ)+CFOPPs1(4,M,NZ)*(RDNDLP+RDNSNP)
6370    CONTINUE
!
!     CONSUMPTION OF NON-STRUCTURAL C,N,P BY NODULE
!
!     CPOOLN,ZPOOLN,PPOOLN=nonstructural C,N,P in bacteria
!     RMNDL=bacterial maintenance respiration
!     RCNDL=respiration from non-structural C
!     RGN2F=respiration for N2 fixation
!     CGNDL=total non-structural C used in bacterial growth and growth respiration
!     RCNDLC,RCNDLN,RCNDLP=bacterial C,N,P decomposition to recycling
!     RCNSNC,RCNSNC,RCNSNP=bacterial C,N,P senescence to recycling
!     ZADDN,PADDN=nonstructural N,P used in growth
!     RUPNF=root N2 fixation
!
        CPOOLNs1(L,NZ)=CPOOLNs1(L,NZ)-AMIN1(RMNDL,RCNDL)-RGN2F-CGNDL+RCNDLC
        ZPOOLNs1(L,NZ)=ZPOOLNs1(L,NZ)-ZADDN+RCNDLN+RCNSNN+RUPNFs1(L,NZ)
        PPOOLNs1(L,NZ)=PPOOLNs1(L,NZ)-PADDN+RCNDLP+RCNSNP
!
!     UPDATE STATE VARIABLES FOR NODULE C, N, P
!
!     WTNDL,WTNDLN,WTNDLP=bacterial C,N,P mass
!     GRNDG=bacterial growth
!     RXNDLC,RXNDLN,RXNDLP=bacterial C,N,P loss from decomposition
!     RXNSNC,RXNSNC,RXNSNP=bacterial C,N,P loss from senescence
!     ZADDN,PADDN=nonstructural N,P used in growth
!
        WTNDLs1(L,NZ)=WTNDLs1(L,NZ)+GRNDG-RXNDLC-RXNSNC
        WTNDLNs1(L,NZ)=WTNDLNs1(L,NZ)+ZADDN-RXNDLN-RXNSNN
        WTNDLPs1(L,NZ)=WTNDLPs1(L,NZ)+PADDN-RXNDLP-RXNSNP
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
        IF(CPOOLRs1(1,L,NZ).GT.ZEROPs1(NZ) &
          .AND.WTRTDs1(1,L,NZ).GT.ZEROLs1(NZ))THEN
          CCNDLR=WTNDLs1(L,NZ)/WTRTDs1(1,L,NZ)
          WTRTD1=WTRTDs1(1,L,NZ)
          WTNDL1=AMIN1(WTRTDs1(1,L,NZ) &
          ,AMAX1(WTNDI*AREA3s1(NUs1),WTNDLs1(L,NZ)))
          WTRTDT=WTRTD1+WTNDL1
          IF(WTRTDT.GT.ZEROPs1(NZ))THEN
            FXRNX=FXRN(INTYPs1(NZ))/(1.0+CCNDLR/CCNGR)
!    2/(1.0+CCNDLR/(CCNGR*FXRN(INTYPs1(NZ))))
            CPOOLD=(CPOOLRs1(1,L,NZ)*WTNDL1-CPOOLNs1(L,NZ)*WTRTD1)/WTRTDT
            XFRC=FXRNX*CPOOLD
            CPOOLRs1(1,L,NZ)=CPOOLRs1(1,L,NZ)-XFRC
            CPOOLNs1(L,NZ)=CPOOLNs1(L,NZ)+XFRC
            CPOOLT=CPOOLRs1(1,L,NZ)+CPOOLNs1(L,NZ)
            IF(CPOOLT.GT.ZEROPs1(NZ))THEN
              ZPOOLD=(ZPOOLRs1(1,L,NZ)*CPOOLNs1(L,NZ) &
                -ZPOOLNs1(L,NZ)*CPOOLRs1(1,L,NZ))/CPOOLT
              XFRN=FXRNX*ZPOOLD
              PPOOLD=(PPOOLRs1(1,L,NZ)*CPOOLNs1(L,NZ) &
                -PPOOLNs1(L,NZ)*CPOOLRs1(1,L,NZ))/CPOOLT
              XFRP=FXRNX*PPOOLD
              ZPOOLRs1(1,L,NZ)=ZPOOLRs1(1,L,NZ)-XFRN
              PPOOLRs1(1,L,NZ)=PPOOLRs1(1,L,NZ)-XFRP
              ZPOOLNs1(L,NZ)=ZPOOLNs1(L,NZ)+XFRN
              PPOOLNs1(L,NZ)=PPOOLNs1(L,NZ)+XFRP
            ENDIF
          ENDIF
        ENDIF
      ENDIF
5400  CONTINUE
  ENDIF
  end associate
  end subroutine RootNoduleBiomchemistry
end module NoduleBGCMod
