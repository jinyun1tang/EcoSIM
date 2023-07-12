module NoduleBGCMod

  use data_kind_mod, only : r8 => DAT_KIND_R8
  use minimathmod, only : safe_adb,AZMAX1
  use EcosimConst
  use PlantAPIData
  use GrosubPars
  implicit none
  private
  character(len=*),private, parameter :: mod_filename = __FILE__
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
  real(r8) :: RXNDLE(npelms)
  real(r8) :: RDNDLE(npelms)
  real(r8) :: RCNDLE(npelms)
  real(r8) :: RGNDG
  real(r8) :: RXNSNE(npelms)
  real(r8) :: RDNSNE(npelms)
  real(r8) :: RCNSNE(npelms)
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
    INTYP    =>  plt_morph%INTYP   , &
    TFN3     =>  plt_pheno%TFN3    , &
    TCO2T    =>  plt_bgcr%TCO2T    , &
    RECO     =>  plt_bgcr%RECO     , &
    TCO2A    =>  plt_bgcr%TCO2A    , &
    TRAU     =>  plt_bgcr%TRAU     , &
    CNET     =>  plt_bgcr%CNET     , &
    ESNC     =>  plt_bgcr%ESNC     , &
    ifoliar  =>  pltpar%ifoliar    , &
    DMND     =>  plt_allom%DMND    , &
    CNND     =>  plt_allom%CNND    , &
    CPND     =>  plt_allom%CPND    , &
    WTLSB    =>  plt_biom%WTLSB    , &
    EPOOL    =>  plt_biom%EPOOL    , &
    EPOLNB   =>  plt_biom%EPOLNB   , &
    ZEROP    =>  plt_biom%ZEROP    , &
    ZEROL    =>  plt_biom%ZEROL    , &
    WTNDBE   =>  plt_biom%WTNDBE     &
  )
!     INTYP=N2 fixation: 4,5,6=rapid to slow canopy symbiosis
!
  IF(INTYP(NZ).GE.4)THEN
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
      CCPOLN=AZMAX1(EPOLNB(ielmc,NB,NZ)/WTNDBE(ielmc,NB,NZ))
      CZPOLN=AZMAX1(EPOLNB(ielmn,NB,NZ)/WTNDBE(ielmc,NB,NZ))
      CPPOLN=AZMAX1(EPOLNB(ielmp,NB,NZ)/WTNDBE(ielmc,NB,NZ))
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
!     TFN3=temperature function for canopy growth
!     FCNPF=N,P constraint to bacterial activity
!     WFNG=growth function of canopy water potential
!
    RCNDL=AZMAX1(AMIN1(EPOLNB(ielmc,NB,NZ),VMXO*WTNDBE(ielmc,NB,NZ))*FCNPF*TFN3(NZ)*WFNG)
!     CPOOLNX=EPOLNB(ielmc,NB,NZ)
!     VMXOX=VMXO*WTNDBE(ielmc,NB,NZ)*FCNPF*TFN3(NZ)*WFNG
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
!     WTLSB=leaf+petiole mass
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
    SPNDX=SPNDL*SQRT(TFN3(NZ)*WFNG)
    DO NE=1,npelms
      RXNDLE(NE)=SPNDX*WTNDBE(NE,NB,NZ)
    ENDDO

    RDNDLE(ielmc)=RXNDLE(ielmc)*(1.0_r8-RCCC)
    RDNDLE(ielmn)=RXNDLE(ielmn)*(1.0_r8-RCCC)*(1.0_r8-RCCN)
    RDNDLE(ielmp)=RXNDLE(ielmp)*(1.0_r8-RCCC)*(1.0_r8-RCCP)

    DO NE=1,npelms
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
    CGNDL=AMIN1(EPOLNB(ielmc,NB,NZ)-AMIN1(RMNDL,RCNDL) &
      -RGN2F+RCNDLE(ielmc),(RGNDL-RGN2F)/(1.0_r8-DMND(NZ)))
    GRNDG=CGNDL*DMND(NZ)
    RGNDG=RGN2F+CGNDL*(1.0_r8-DMND(NZ))
    ZADDN=AZMAX1(AMIN1(EPOLNB(ielmn,NB,NZ),GRNDG*CNND(NZ)))*CZPOLN/(CZPOLN+CZKM)
    PADDN=AZMAX1(AMIN1(EPOLNB(ielmp,NB,NZ),GRNDG*CPND(NZ)))*CPPOLN/(CPPOLN+CPKM)
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
      DO NE=1,npelms
        RCNSNE(NE)=RXNSNE(NE)-RDNSNE(NE)
      ENDDO
    ELSE
      RXNSNE(1:npelms)=0._r8
      RDNSNE(1:npelms)=0._r8
      RCNSNE(1:npelms)=0._r8
    ENDIF
!
!     TOTAL NODULE RESPIRATION
!
!     RCO2T=total C respiration
!     RMNDL=bacterial maintenance respiration
!     RCNDL=respiration from non-structural C
!     RGNDG=bacterial respiration for growth and N2 fixation
!     RCNSNE(ielmc)=bacterial C senescence to recycling
!     TCO2T,TCO2A=total,above-ground PFT respiration
!     CNET=PFT net CO2 fixation
!     RECO=ecosystem respiration
!     TRAU=total autotrophic respiration
!
    RCO2T=AMIN1(RMNDL,RCNDL)+RGNDG+RCNSNE(ielmc)
    TCO2T(NZ)=TCO2T(NZ)-RCO2T
    TCO2A(NZ)=TCO2A(NZ)-RCO2T
    CNET(NZ)=CNET(NZ)-RCO2T
    RECO=RECO-RCO2T
    TRAU=TRAU-RCO2T
!
!     NODULE LITTERFALL CAUSED BY REMOBILIZATION
!
!     CSNC,ZSNC,PSNC=C,N,P litterfall from decomposition and senescence
!     CFOPC,CFOPN,CFOPC=fraction of litterfall C,N,P allocated to litter components
!     RDNDLE(ielmc),RDNDLE(ielmn),RDNDLE(ielmp)=bacterial C,N,P decomposition to litterfall
!     RDNSNE(ielmc),RDNSNE(ielmc),RDNSNE(ielmp)=bacterial C,N,P senescence to litterfall
!
    D6470: DO M=1,jsken
      DO NE=1,npelms
        ESNC(NE,M,k_fine_litr,0,NZ)=ESNC(NE,M,k_fine_litr,0,NZ) &
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
    EPOLNB(ielmc,NB,NZ)=EPOLNB(ielmc,NB,NZ)-AMIN1(RMNDL,RCNDL)-RGN2F-CGNDL+RCNDLE(ielmc)
    EPOLNB(ielmn,NB,NZ)=EPOLNB(ielmn,NB,NZ)-ZADDN+RCNDLE(ielmn)+RCNSNE(ielmn)+RUPNFB
    EPOLNB(ielmp,NB,NZ)=EPOLNB(ielmp,NB,NZ)-PADDN+RCNDLE(ielmp)+RCNSNE(ielmp)
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
!     WTLSB=leaf+petiole C mass
!     WTNDB=bacterial C mass
!     WTNDI=initial bacterial mass at infection
!     FXRN=rate constant for plant-bacteria nonstructural C,N,P exchange
!     CCNGB=parameter to calculate nonstructural C,N,P exchange
!     CCNDLB=bacteria:leaf+petiole ratio
!     XFRC,XFRN,XFRC=nonstructural C,N,P transfer
!     CPOLNB,ZPOLNB,PPOLNB=nonstructural C,N,P in bacteria
!
    IF(EPOOL(ielmc,NB,NZ).GT.ZEROP(NZ).AND.WTLSB(NB,NZ).GT.ZEROL(NZ))THEN
      CCNDLB=WTNDBE(ielmc,NB,NZ)/WTLSB(NB,NZ)
      WTLSB1=WTLSB(NB,NZ)
      WTNDB1=AMIN1(WTLSB(NB,NZ),AMAX1(WTNDI*AREA3(NU),WTNDBE(ielmc,NB,NZ)))
      WTLSBT=WTLSB1+WTNDB1
      IF(WTLSBT.GT.ZEROP(NZ))THEN
        FXRNX=FXRN(INTYP(NZ))/(1.0+CCNDLB/CCNGB)
!    2/(1.0+CCNDLB/(CCNGB*FXRN(INTYP(NZ))))
        CPOOLD=(EPOOL(ielmc,NB,NZ)*WTNDB1-EPOLNB(ielmc,NB,NZ)*WTLSB1)/WTLSBT
        XFRC=FXRNX*CPOOLD
        EPOOL(ielmc,NB,NZ)=EPOOL(ielmc,NB,NZ)-XFRC
        EPOLNB(ielmc,NB,NZ)=EPOLNB(ielmc,NB,NZ)+XFRC
        CPOOLT=EPOOL(ielmc,NB,NZ)+EPOLNB(ielmc,NB,NZ)
        IF(CPOOLT.GT.ZEROP(NZ))THEN
          ZPOOLD=(EPOOL(ielmn,NB,NZ)*EPOLNB(ielmc,NB,NZ)-EPOLNB(ielmn,NB,NZ)*EPOOL(ielmc,NB,NZ))/CPOOLT
          XFRN=FXRNX*ZPOOLD
          PPOOLD=(EPOOL(ielmp,NB,NZ)*EPOLNB(ielmc,NB,NZ)-EPOLNB(ielmp,NB,NZ)*EPOOL(ielmc,NB,NZ))/CPOOLT
          XFRP=FXRNX*PPOOLD
          EPOOL(ielmn,NB,NZ)=EPOOL(ielmn,NB,NZ)-XFRN
          EPOOL(ielmp,NB,NZ)=EPOOL(ielmp,NB,NZ)-XFRP
          EPOLNB(ielmn,NB,NZ)=EPOLNB(ielmn,NB,NZ)+XFRN
          EPOLNB(ielmp,NB,NZ)=EPOLNB(ielmp,NB,NZ)+XFRP
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
  real(r8) :: RXNDLE(npelms)
  real(r8) :: RDNDLE(npelms)
  real(r8) :: RCNDLE(npelms)
  real(r8) :: RGNDG
  real(r8) :: RXNSNE(npelms)
  real(r8) :: RDNSNE(npelms)
  real(r8) :: RCNSNE(npelms)
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
    TFN4     =>   plt_pheno%TFN4     , &
    DMND     =>   plt_allom%DMND     , &
    CNND     =>   plt_allom%CNND     , &
    CPND     =>   plt_allom%CPND     , &
    k_fine_litr=> pltpar%k_fine_litr , &
    iroot    =>   pltpar%iroot       , &
    RCO2M    =>   plt_rbgc%RCO2M     , &
    RCO2N    =>   plt_rbgc%RCO2N     , &
    WFR      =>   plt_rbgc%WFR       , &
    RCO2A    =>   plt_rbgc%RCO2A     , &
    ESNC     =>   plt_bgcr%ESNC      , &
    UPNF     =>   plt_rbgc%UPNF      , &
    RUPNF    =>   plt_bgcr%RUPNF     , &
    WTRTD    =>   plt_biom%WTRTD     , &
    WTNDLE   =>   plt_biom%WTNDLE    , &
    ZEROP    =>   plt_biom%ZEROP     , &
    EPOOLN   =>   plt_biom%EPOOLN    , &
    ZEROL    =>   plt_biom%ZEROL     , &
    EPOOLR   =>   plt_biom%EPOOLR    , &
    CFOPE    =>   plt_soilchem%CFOPE , &
    INTYP    =>   plt_morph%INTYP    , &
    NIX      =>   plt_morph%NIX        &
  )
!     INTYP=N2 fixation: 1,2,3=rapid to slow root symbiosis
!     WTNDL,WTNDLN,WTNDLP=bacterial C,N,P mass
!     WTNDI=initial bacterial mass at infection
!     AREA=grid cell area
!     CNND,CPND=bacterial N:C,P:C ratio from PFT file
!
  IF(INTYP(NZ).GE.1.AND.INTYP(NZ).LE.3)THEN
    D5400: DO L=NU,NIX(NZ)
      IF(WTRTD(ipltroot,L,NZ).GT.ZEROL(NZ))THEN
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
          CCPOLN=AZMAX1(EPOOLN(ielmc,L,NZ)/WTNDLE(ielmc,L,NZ))
          CZPOLN=AZMAX1(EPOOLN(ielmn,L,NZ)/WTNDLE(ielmc,L,NZ))
          CPPOLN=AZMAX1(EPOOLN(ielmp,L,NZ)/WTNDLE(ielmc,L,NZ))
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
!     TFN4=temperature function for root growth
!     FCNPF=N,P constraint to bacterial activity
!     WFNGR=growth function of root water potential
!
        RCNDLM=AZMAX1(AMIN1(EPOOLN(ielmc,L,NZ) &
          ,VMXO*WTNDLE(ielmc,L,NZ))*FCNPF*TFN4(L,NZ)*WFNGR(1,L))
        CPOOLNX=EPOOLN(ielmc,L,NZ)
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
        SPNDX=SPNDL*SQRT(TFN4(L,NZ)*WFNGR(1,L))
        DO NE=1,npelms
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
        CGNDL=AMIN1(EPOOLN(ielmc,L,NZ)-AMIN1(RMNDL,RCNDL) &
          -RGN2F+RCNDLE(ielmc),(RGNDL-RGN2F)/(1.0_r8-DMND(NZ)))
        GRNDG=CGNDL*DMND(NZ)
        RGNDG=RGN2F+CGNDL*(1.0_r8-DMND(NZ))
        ZADDN=AZMAX1(AMIN1(EPOOLN(ielmn,L,NZ),GRNDG*CNND(NZ)))*CZPOLN/(CZPOLN+CZKM)
        PADDN=AZMAX1(AMIN1(EPOOLN(ielmp,L,NZ),GRNDG*CPND(NZ)))*CPPOLN/(CPPOLN+CPKM)
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
          DO NE=1,npelms
            RCNSNE(NE)=RXNSNE(NE)-RDNSNE(NE)
          ENDDO
        ELSE
          RXNSNE(1:npelms)=0._r8
          RDNSNE(1:npelms)=0._r8
          RCNSNE(1:npelms)=0._r8
        ENDIF
!
!     TOTAL NODULE RESPIRATION
!
!     RCO2TM,RCO2T=total C respiration unlimited,limited by O2
!     TCO2T,TCO2A=total,above-ground PFT respiration
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
          DO NE=1,npelms
            ESNC(NE,M,k_fine_litr,L,NZ)=ESNC(NE,M,k_fine_litr,L,NZ)&
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
        EPOOLN(ielmc,L,NZ)=EPOOLN(ielmc,L,NZ)-AMIN1(RMNDL,RCNDL)-RGN2F-CGNDL+RCNDLE(ielmc)
        EPOOLN(ielmn,L,NZ)=EPOOLN(ielmn,L,NZ)-ZADDN+RCNDLE(ielmn)+RCNSNE(ielmn)+RUPNF(L,NZ)
        EPOOLN(ielmp,L,NZ)=EPOOLN(ielmp,L,NZ)-PADDN+RCNDLE(ielmp)+RCNSNE(ielmp)
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
        IF(EPOOLR(ielmc,ipltroot,L,NZ).GT.ZEROP(NZ) &
          .AND.WTRTD(ipltroot,L,NZ).GT.ZEROL(NZ))THEN
          CCNDLR=WTNDLE(ielmc,L,NZ)/WTRTD(ipltroot,L,NZ)
          WTRTD1=WTRTD(ipltroot,L,NZ)
          WTNDL1=AMIN1(WTRTD(ipltroot,L,NZ),AMAX1(WTNDI*AREA3(NU),WTNDLE(ielmc,L,NZ)))
          WTRTDT=WTRTD1+WTNDL1
          IF(WTRTDT.GT.ZEROP(NZ))THEN
            FXRNX=FXRN(INTYP(NZ))/(1.0_r8+CCNDLR/CCNGR)
!    2/(1.0+CCNDLR/(CCNGR*FXRN(INTYP(NZ))))
            CPOOLD=(EPOOLR(ielmc,ipltroot,L,NZ)*WTNDL1-EPOOLN(ielmc,L,NZ)*WTRTD1)/WTRTDT
            XFRC=FXRNX*CPOOLD
            EPOOLR(ielmc,ipltroot,L,NZ)=EPOOLR(ielmc,ipltroot,L,NZ)-XFRC
            EPOOLN(ielmc,L,NZ)=EPOOLN(ielmc,L,NZ)+XFRC
            CPOOLT=EPOOLR(ielmc,ipltroot,L,NZ)+EPOOLN(ielmc,L,NZ)
            IF(CPOOLT.GT.ZEROP(NZ))THEN
              ZPOOLD=(EPOOLR(ielmn,ipltroot,L,NZ)*EPOOLN(ielmc,L,NZ) &
                -EPOOLN(ielmn,L,NZ)*EPOOLR(ielmc,ipltroot,L,NZ))/CPOOLT
              XFRN=FXRNX*ZPOOLD
              PPOOLD=(EPOOLR(ielmp,ipltroot,L,NZ)*EPOOLN(ielmc,L,NZ) &
                -EPOOLN(ielmp,L,NZ)*EPOOLR(ielmc,ipltroot,L,NZ))/CPOOLT
              XFRP=FXRNX*PPOOLD
              EPOOLR(ielmn,ipltroot,L,NZ)=EPOOLR(ielmn,ipltroot,L,NZ)-XFRN
              EPOOLR(ielmp,ipltroot,L,NZ)=EPOOLR(ielmp,ipltroot,L,NZ)-XFRP
              EPOOLN(ielmn,L,NZ)=EPOOLN(ielmn,L,NZ)+XFRN
              EPOOLN(ielmp,L,NZ)=EPOOLN(ielmp,L,NZ)+XFRP
            ENDIF
          ENDIF
        ENDIF
      ENDIF
    ENDDO D5400
  ENDIF
  end associate
  end subroutine RootNoduleBiomchemistry
end module NoduleBGCMod
