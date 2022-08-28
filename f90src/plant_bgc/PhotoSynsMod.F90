module PhotoSynsMod
  use data_kind_mod, only : r8 => SHR_KIND_R8
  use EcosimConst
  use GrosubPars
  use PlantAPIData
implicit none
  private
  public :: ComputeGPP

  contains

!------------------------------------------------------------------------------------------

  subroutine ComputeGPP_C3(K,NB,NZ,WFNG,WFNC,CH2O3,CO2F,CH2O)

  implicit none
  integer, intent(in) :: K,NB,NZ
  real(r8), intent(in) :: WFNG,WFNC
  real(r8), intent(inout) :: CH2O3(25),CO2F,CH2O
  integer :: L,NN,M,N
  real(r8) :: WFNB
  real(r8) :: CO2X,CO2C,CO2Y
  real(r8) :: CBXNX
  real(r8) :: DIFF
  real(r8) :: EGROX
  real(r8) :: ETLF
  real(r8) :: EGRO
  real(r8) :: GSL
  real(r8) :: PARX,PARJ
  real(r8) :: RS,RSL
  real(r8) :: VL,VGROX
  real(r8) :: VA,VG
!begin_execution
  associate(                          &
  ZEROs1       => plt_site%ZEROs1   , &
  IGTYPs1      => plt_pheno%IGTYPs1 , &
  FDBKs1       => plt_photo%FDBKs1  , &
  CO2Qs1       => plt_photo%CO2Qs1  , &
  SURFXs1      => plt_photo%SURFXs1 , &
  XKCO2Os1     => plt_photo%XKCO2Os1, &
  SCO2s1       => plt_photo%SCO2s1  , &
  CBXNs1       => plt_photo%CBXNs1  , &
  ETGROs1      => plt_photo%ETGROs1 , &
  VGROs1       => plt_photo%VGROs1  , &
  CO2Is1       => plt_photo%CO2Is1  , &
  RCMXs1       => plt_photo%RCMXs1  , &
  XKCO24s1     => plt_photo%XKCO24s1, &
  VCGROs1      => plt_photo%VCGROs1 , &
  COMPLs1      => plt_photo%COMPLs1 , &
  FMOLs1       => plt_photo%FMOLs1  , &
  DCO2s1       => plt_photo%DCO2s1  , &
  ZEROPs1      => plt_biom%ZEROPs1  , &
  ARLFLs1      => plt_morph%ARLFLs1 , &
  PARDIFs1     => plt_rad%PARDIFs1  , &
  PARs1        => plt_rad%PARs1     , &
  TAU0s1       => plt_rad%TAU0s1    , &
  TAUSs1       => plt_rad%TAUSs1      &
  )

! FOR EACH CANOPY LAYER
!
  DO 210 L=JC1,1,-1
    IF(ARLFLs1(L,K,NB,NZ).GT.ZEROPs1(NZ))THEN
!
!     FOR EACH LEAF AZIMUTH AND INCLINATION
!
      DO 215 N=1,JLI1
        DO 220 M=1,JSA1
!
!         CO2 FIXATION BY SUNLIT LEAVES
!
!         SURFX=unself-shaded leaf surface area
!
          IF(SURFXs1(N,L,K,NB,NZ).GT.ZEROPs1(NZ))THEN
            IF(PARs1(N,M,L,NZ).GT.0.0)THEN
!
!             C3 CARBOXYLATION REACTIONS IN MESOPHYLL
!
!             QNTM=quantum efficiency
!             PAR=direct PAR flux
!             ETGRO=light saturated e- transport rate from stomate.f
!             ETLF=light-limited e- transport rate
!             CURV=shape parameter for e- transport response to PAR
!             EGRO=light-limited rubisco carboxylation rate
!             CBXN=rubisco caboxylation efficiency
!             VL=rubisco carboxylation rate limited by light,CO2,N,P
!             VGRO=rubisco carboxylation rate limited by CO2 from stomate.f
!             FDBK=N,P feedback inhibition on C4 CO2 fixation
!
              PARX=QNTM*PARs1(N,M,L,NZ)
              PARJ=PARX+ETGROs1(K,NB,NZ)
              ETLF=(PARJ-SQRT(PARJ**2-CURV4*PARX*ETGROs1(K,NB,NZ)))/CURV2
              EGRO=ETLF*CBXNs1(K,NB,NZ)
              VL=AMIN1(VGROs1(K,NB,NZ),EGRO)*FDBKs1(NB,NZ)
!
!             STOMATAL EFFECT OF WATER DEFICIT IN MESOPHYLL
!
!             RS,RSL=leaf stomatal resistance to CO2 at zero,current water potential
!             RCMN=minimum stomatal resistance to CO2 (s m-1)
!             RCMX=cuticular resistance to CO2 from startq.f (s m-1)
!             DCO2=difference between atmosph and intercellular CO2 concn (umol m-3)
!             GSL=leaf stomatal conductance (mol m-2 s-1)
!             WFNC=stomatal resistance function of canopy turgor
!             FMOL=number of moles of air per m3
!
              IF(VL.GT.ZEROs1)THEN
                RS=AMIN1(RCMXs1(NZ),AMAX1(RCMN,DCO2s1(NZ)/VL))
                RSL=RS+(RCMXs1(NZ)-RS)*WFNC
                GSL=1.0_r8/RSL*FMOLs1(NZ)
!
!               EFFECT OF WATER DEFICIT IN MESOPHYLL
!
!               IGTYP=growth type:0=bryophyte,1=graminoid,2=shrub,tree
!               WFNB=non-stomatal effects of water stress on CO2 fixation
!
                IF(IGTYPs1(NZ).NE.0)THEN
                  WFNB=SQRT(RS/RSL)
                ELSE
                  WFNB=WFNG
                ENDIF
!
!               CONVERGENCE SOLUTION FOR CO2I AT WHICH CARBOXYLATION
!               EQUALS DIFFUSION IN MESOPHYLL
!
!               CO2I=intercellular,mesophyll CO2 concentration at zero water potential
!               CO2X,CO2C=intercellular,mesophyll CO2 concentration during convergence
!               SCO2=solubility of CO2 (uM/(umol mol-1))
!               COMPL=C3 CO2 compensation point (uM)
!               CBXNX=rubisco carboxylation efficiency
!               ELEC3=e- requirement for CO2 fixn by rubisco
!               VCGRO,VGROX=rubisco carboxylation rate unlimited,limited by CO2
!               XKCO2O=Km for rubisco carboxylation
!               EGROX=light-limited rubisco carboxylation rate
!               ETLF=light-limited e- transport rate
!               VL=rubisco carboxylation rate limited by light,CO2,N,P,water stress
!               VG=CO2 diffusion rate limited by water stress
!               GSL=leaf stomatal conductance (mol m-2 s-1)
!
                CO2X=CO2Is1(NZ)
                DO 225 NN=1,100
                  CO2C=CO2X*SCO2s1(NZ)
                  CO2Y=AMAX1(0.0,CO2C-COMPLs1(K,NB,NZ))
                  CBXNX=CO2Y/(ELEC3*CO2C+10.5*COMPLs1(K,NB,NZ))
                  VGROX=VCGROs1(K,NB,NZ)*CO2Y/(CO2C+XKCO2Os1(NZ))
                  EGROX=ETLF*CBXNX
                  VL=AMIN1(VGROX,EGROX)*WFNB*FDBKs1(NB,NZ)
                  VG=(CO2Qs1(NZ)-CO2X)*GSL
                  IF(VL+VG.GT.ZEROs1)THEN
                    DIFF=(VL-VG)/(VL+VG)
                    IF(ABS(DIFF).LT.0.005)exit
                    VA=0.95*VG+0.05*VL
                    CO2X=CO2Qs1(NZ)-VA/GSL
                  ELSE
                    VL=0._r8
                    exit
                  ENDIF
225             CONTINUE
!
!               ACCUMULATE C3 FIXATION PRODUCT IN MESOPHYLL
!
!               CH2O3=total C4 CO2 fixation
!               SURFX=unself-shaded leaf surface area
!               TAUS=fraction of direct radiation transmitted from layer above
!
                CH2O3(K)=CH2O3(K)+VL*SURFXs1(N,L,K,NB,NZ) &
                  *TAUSs1(L+1)
!               ICO2I=MAX(1,MIN(400,INT(CO2X)))
!               VCO2(ICO2I,I,NZ)=VCO2(ICO2I,I,NZ)
!              2+(VL*SURFXs1(N,L,K,NB,NZ)*TAUSs1(L+1))*0.0432

              ENDIF
            ENDIF
!
!           CO2 FIXATION IN MESOPHYLL BY SHADED LEAVES
!
            IF(PARDIFs1(N,M,L,NZ).GT.0.0)THEN
!
!             C3 CARBOXYLATION REACTIONS USING VARIABLES FROM 'STOMATE'
!
!             QNTM=quantum efficiency
!             PARDIFs1=diffuse PAR flux
!             ETGRO=light saturated e- transport rate from stomate.f
!             ETLF=light-limited e- transport rate
!             CURV=shape parameter for e- transport response to PAR
!             EGRO=light-limited rubisco carboxylation rate
!             CBXN=rubisco caboxylation efficiency
!             VL=rubisco carboxylation rate limited by light,CO2,N,P
!             VGRO=rubisco carboxylation rate limited by CO2 from stomate.f
!             FDBK=N,P feedback inhibition on C3 CO2 fixation
!
              PARX=QNTM*PARDIFs1(N,M,L,NZ)
              PARJ=PARX+ETGROs1(K,NB,NZ)
              ETLF=(PARJ-SQRT(PARJ**2-CURV4*PARX*ETGROs1(K,NB,NZ)))/CURV2
              EGRO=ETLF*CBXNs1(K,NB,NZ)
              VL=AMIN1(VGROs1(K,NB,NZ),EGRO)*FDBKs1(NB,NZ)
!
!             STOMATAL EFFECT OF WATER DEFICIT IN MESOPHYLL
!
!             RS,RSL=leaf stomatal resistance to CO2 at zero,current water potential
!             RCMN=minimum stomatal resistance to CO2 (s m-1)
!             RCMX=cuticular resistance to CO2 from startq.f (s m-1)
!             DCO2=difference between atmosph and intercellular CO2 concn (umol m-3)
!             GSL=leaf stomatal conductance (mol m-2 s-1)
!             WFNC=stomatal resistance function of canopy turgor
!             FMOL=number of moles of air per m3
!
              IF(VL.GT.ZEROs1)THEN
                RS=AMIN1(RCMXs1(NZ),AMAX1(RCMN,DCO2s1(NZ)/VL))
                RSL=RS+(RCMXs1(NZ)-RS)*WFNC
                GSL=1.0_r8/RSL*FMOLs1(NZ)
!
!               EFFECT OF WATER DEFICIT IN MESOPHYLL
!
!               IGTYP=growth type:0=bryophyte,1=graminoid,2=shrub,tree
!               WFNB=non-stomatal effects of water stress on C3 CO2 fixation
!
                IF(IGTYPs1(NZ).NE.0)THEN
                  WFNB=SQRT(RS/RSL)
                ELSE
                  WFNB=WFNG
                ENDIF
!
!               CONVERGENCE SOLUTION FOR CO2I AT WHICH CARBOXYLATION
!               EQUALS DIFFUSION IN MESOPHYLL
!
!               CO2I=intercellular,mesophyll CO2 concentration at zero water potential
!               CO2X,CO2C=intercellular,mesophyll CO2 concentration during convergence
!               SCO2=solubility of CO2 (uM/(umol mol-1))
!               COMPL=C3 CO2 compensation point (uM)
!               CBXNX=rubisco caboxylation efficiency
!               ELEC3=e- requirement for CO2 fixn by rubisco carboxylase
!               VCGRO,VGROX=rubisco carboxylation rate unlimited,limited by CO2
!               XKCO2O=Km for rubisco carboxylation from stomate.f (uM)
!               EGROX=light-limited rubisco carboxylation rate
!               ETLF=light-limited e- transport rate
!               VL=rubisco carboxylation rate limited by light,CO2,N,P,water stress
!               VG=CO2 diffusion rate limited by water stress
!               GSL=leaf stomatal conductance (mol m-2 s-1)
!
                CO2X=CO2Is1(NZ)
                DO 235 NN=1,100
                  CO2C=CO2X*SCO2s1(NZ)
                  CO2Y=AMAX1(0.0,CO2C-COMPLs1(K,NB,NZ))
                  CBXNX=CO2Y/(ELEC3*CO2C+10.5*COMPLs1(K,NB,NZ))
                  VGROX=VCGROs1(K,NB,NZ)*CO2Y/(CO2C+XKCO2Os1(NZ))
                  EGROX=ETLF*CBXNX
                  VL=AMIN1(VGROX,EGROX)*WFNB*FDBKs1(NB,NZ)
                  VG=(CO2Qs1(NZ)-CO2X)*GSL
                  IF(VL+VG.GT.ZEROs1)THEN
                    DIFF=(VL-VG)/(VL+VG)
                    IF(ABS(DIFF).LT.0.005)exit
                    VA=0.95*VG+0.05*VL
                    CO2X=CO2Qs1(NZ)-VA/GSL
                  ELSE
                    VL=0._r8
                    exit
                  ENDIF
235             CONTINUE
!
!               ACCUMULATE C3 FIXATION PRODUCT IN MESOPHYLL
!
!               CH2O3=total C3 CO2 fixation
!               SURFX=unself-shaded leaf surface area
!               TAU0s1=fraction of diffuse radiation transmitted from layer above
!
                CH2O3(K)=CH2O3(K)+VL*SURFXs1(N,L,K,NB,NZ) &
                  *TAU0s1(L+1)
!               ICO2I=MAX(1,MIN(400,INT(CO2X)))
!               VCO2(ICO2I,I,NZ)=VCO2(ICO2I,I,NZ)
!              2+(VL*SURFXs1(N,L,K,NB,NZ)*TAU0s1(L+1))*0.0432
              ENDIF
            ENDIF
          ENDIF
220     CONTINUE
215   CONTINUE
    ENDIF
210   CONTINUE
  CO2F=CO2F+CH2O3(K)
  CH2O=CH2O+CH2O3(K)
  end associate
  end subroutine ComputeGPP_C3


!------------------------------------------------------------------------------------------

  subroutine ComputeGPP_C4(K,NB,NZ,WFNG,WFNC,CH2O3,CH2O4,CO2F,CH2O)
  implicit none
  integer, intent(in) :: K,NB,NZ
  real(r8), intent(in):: WFNG,WFNC
  real(r8), intent(inout) :: CH2O3(25),CH2O4(25),CO2F,CH2O
  integer :: L,NN,M,N
  real(r8) :: WFN4
  real(r8) :: WFNB
  real(r8) :: CO2X,CO2C,CO2Y
  real(r8) :: CBXNX
  real(r8) :: DIFF
  real(r8) :: ETLF4
  real(r8) :: EGRO4
  real(r8) :: EGROX
  real(r8) :: ETLF
  real(r8) :: EGRO
  real(r8) :: GSL
  real(r8) :: PARX,PARJ
  real(r8) :: RS,RSL,VL
  real(r8) :: VGROX
  real(r8) :: VA,VG
! begin_execution
  associate(                          &
  IGTYPs1      => plt_pheno%IGTYPs1 , &
  ZEROPs1      => plt_biom%ZEROPs1  , &
  XKCO24s1     => plt_photo%XKCO24s1, &
  FDBK4s1      => plt_photo%FDBK4s1 , &
  SCO2s1       => plt_photo%SCO2s1  , &
  CBXN4s1      => plt_photo%CBXN4s1 , &
  SURFXs1      => plt_photo%SURFXs1 , &
  VGRO4s1      => plt_photo%VGRO4s1 , &
  RCMXs1       => plt_photo%RCMXs1  , &
  DCO2s1       => plt_photo%DCO2s1  , &
  FMOLs1       => plt_photo%FMOLs1  , &
  ETGROs1      => plt_photo%ETGROs1 , &
  VGROs1       => plt_photo%VGROs1  , &
  ETGR4s1      => plt_photo%ETGR4s1 , &
  CO2Qs1       => plt_photo%CO2Qs1  , &
  FDBKs1       => plt_photo%FDBKs1  , &
  CO2Is1       => plt_photo%CO2Is1  , &
  CBXNs1       => plt_photo%CBXNs1  , &
  VCGR4s1      => plt_photo%VCGR4s1 , &
  ARLFLs1      => plt_morph%ARLFLs1 , &
  ZEROs1       => plt_site%ZEROs1   , &
  PARDIFs1     => plt_rad%PARDIFs1  , &
  PARs1        => plt_rad%PARs1     , &
  TAU0s1       => plt_rad%TAU0s1    , &
  TAUSs1       => plt_rad%TAUSs1      &
  )
! FOR EACH CANOPY LAYER
!
  DO 110 L=JC1,1,-1
    IF(ARLFLs1(L,K,NB,NZ).GT.ZEROPs1(NZ))THEN
!
!     FOR EACH LEAF AZIMUTH AND INCLINATION
!
      DO 115 N =1,JLI1
        DO 120 M =1,JSA1
!
!         CO2 FIXATION IN MESOPHYLL BY SUNLIT LEAVES
!
!         SURFX=unself-shaded leaf surface area
!
          IF(SURFXs1(N,L,K,NB,NZ).GT.ZEROPs1(NZ))THEN
            IF(PARs1(N,M,L,NZ).GT.0.0)THEN
!
!             C4 CARBOXYLATION REACTIONS IN MESOPHYLL
!
!             QNTM=quantum efficiency
!             PAR=direct PAR flux
!             ETGR4=light saturated e- transport rate from stomate.f
!             ETLF4=light-limited e- transport rate
!             CURV=shape parameter for e- transport response to PAR
!             EGRO4=light-limited PEP carboxylation rate
!             CBXN4=PEP caboxylation efficiency
!             VL=PEP carboxylation rate limited by light,CO2,N,P
!             VGRO4=PEP carboxylation rate limited by CO2 from stomate.f
!             FDBK4=N,P feedback inhibition on C4 CO2 fixation
!
              PARX=QNTM*PARs1(N,M,L,NZ)
              PARJ=PARX+ETGR4s1(K,NB,NZ)
              ETLF4=(PARJ-SQRT(PARJ**2-CURV4*PARX*ETGR4s1(K,NB,NZ)))/CURV2
              EGRO4=ETLF4*CBXN4s1(K,NB,NZ)
              VL=AMIN1(VGRO4s1(K,NB,NZ),EGRO4)*FDBK4s1(K,NB,NZ)
!
!             STOMATAL EFFECT OF WATER DEFICIT IN MESOPHYLL
!
!             RS,RSL=leaf stomatal resistance to CO2 at zero,current water potential
!             RCMN=minimum stomatal resistance to CO2 (s m-1)
!             RCMX=cuticular resistance to CO2 from startq.f (s m-1)
!             DCO2=difference between atmosph and intercellular CO2 concn (umol m-3)
!             GSL=leaf stomatal conductance (mol m-2 s-1)
!             WFNC=stomatal resistance function of canopy turgor
!             FMOL=number of moles of air per m3
!
              IF(VL.GT.ZEROs1)THEN
                RS=AMIN1(RCMXs1(NZ),AMAX1(RCMN,DCO2s1(NZ)/VL))
                RSL=RS+(RCMXs1(NZ)-RS)*WFNC
                GSL=1.0_r8/RSL*FMOLs1(NZ)
!
!               EFFECT OF WATER DEFICIT IN MESOPHYLL
!
!               IGTYP=growth type:0=bryophyte,1=graminoid,2=shrub,tree
!               WFNB=non-stomatal effects of water stress on C4,C3 CO2 fixation
!
                IF(IGTYPs1(NZ).NE.0)THEN
                  WFN4=RS/RSL
                  WFNB=SQRT(RS/RSL)
                ELSE
                  WFN4=WFNG
                  WFNB=WFNG
                ENDIF
!
!               CONVERGENCE SOLUTION FOR CO2I AT WHICH CARBOXYLATION
!               EQUALS DIFFUSION IN MESOPHYLL
!
!               CO2I=intercellular,mesophyll CO2 concentration at zero water potential
!               CO2X,CO2C=intercellular,mesophyll CO2 concentration during convergence
!               SCO2=solubility of CO2 (uM/(umol mol-1))
!               COMP4=C4 CO2 compensation point (uM)
!               CBXNX=PEP carboxylation efficiency
!               ELEC4=e- requirement for CO2 fixn by PEP carboxylase
!               VCGR4,VGROX=PEP carboxylation rate unlimited,limited by CO2
!               XKCO24=Km for VCMX4 from PFT file (uM)
!               EGROX=light-limited PEP carboxylation rate
!               ETLF4=light-limited e- transport rate
!               VL=PEP carboxylation rate limited by light,CO2,N,P,water stress
!               VG=CO2 diffusion rate limited by water stress
!               GSL=leaf stomatal conductance (mol m-2 s-1)
!
                CO2X=CO2Is1(NZ)
                DO 125 NN=1,100
                  CO2C=CO2X*SCO2s1(NZ)
                  CO2Y=AMAX1(0.0,CO2C-COMP4)
                  CBXNX=CO2Y/(ELEC4*CO2C+10.5*COMP4)
                  VGROX=VCGR4s1(K,NB,NZ)*CO2Y/(CO2C+XKCO24s1(NZ))
                  EGROX=ETLF4*CBXNX
                  VL=AMIN1(VGROX,EGROX)*WFN4*FDBK4s1(K,NB,NZ)
                  VG=(CO2Qs1(NZ)-CO2X)*GSL
                  IF(VL+VG.GT.ZEROs1)THEN
                    DIFF=(VL-VG)/(VL+VG)
                    IF(ABS(DIFF).LT.0.005)exit
                    VA=0.95*VG+0.05*VL
                    CO2X=CO2Qs1(NZ)-VA/GSL
                  ELSE
                    VL=0._r8
                    exit
                  ENDIF
125             CONTINUE
!
!               ACCUMULATE C4 FIXATION PRODUCT IN MESOPHYLL
!
!               CH2O4=total C4 CO2 fixation
!               SURFX=unself-shaded leaf surface area
!               TAUS=fraction of direct radiation transmitted from layer above
!
                CH2O4(K)=CH2O4(K)+VL*SURFXs1(N,L,K,NB,NZ) &
                  *TAUSs1(L+1)
!               ICO2I=MAX(1,MIN(400,INT(CO2X)))
!               VCO2(ICO2I,I,NZ)=VCO2(ICO2I,I,NZ)
!              2+(VL*SURFXs1(N,L,K,NB,NZ)*TAUSs1(L+1))*0.0432
!
!               C3 CARBOXYLATION REACTIONS IN BUNDLE SHEATH OF C4 PLANTS
!
!               ETGRO=light saturated e- transport rate from stomate.f
!               ETLF=light-limited e- transport rate
!               CURV=shape parameter for e- transport response to PAR
!               EGRO=light-limited rubisco carboxylation rate
!               CBXN=rubisco caboxylation efficiency
!               VL=rubisco carboxylation rate limited by light,CO2,N,P
!               VGRO=rubisco carboxylation rate limited by CO2 from stomate.f
!               FDBK=N,P feedback inhibition on C3 CO2 fixation
!
                PARJ=PARX+ETGROs1(K,NB,NZ)
                ETLF=(PARJ-SQRT(PARJ**2-CURV4*PARX*ETGROs1(K,NB,NZ)))/CURV2
                EGRO=ETLF*CBXNs1(K,NB,NZ)
                VL=AMIN1(VGROs1(K,NB,NZ),EGRO)*WFNB*FDBKs1(NB,NZ)
!
!               ACCUMULATE C3 FIXATION PRODUCT IN BUNDLE SHEATH
!
!               CH2O3=total C3 CO2 fixation
!               SURFX=unself-shaded leaf surface area
!               TAUS=fraction of direct radiation transmitted from layer above
!
                CH2O3(K)=CH2O3(K)+VL*SURFXs1(N,L,K,NB,NZ)*TAUSs1(L+1)

              ENDIF
            ENDIF
!
!           CO2 FIXATION IN MESOPHYLL BY SHADED LEAVES
!
            IF(PARDIFs1(N,M,L,NZ).GT.0.0)THEN
!
!           C4 CARBOXYLATION REACTIONS IN MESOPHYLL
!
!           QNTM=quantum efficiency
!           PARDIFs1=diffuse PAR flux
!           ETGR4=light saturated e- transport rate from stomate.f
!           ETLF4=light-limited e- transport rate
!           CURV=shape parameter for e- transport response to PAR
!           EGRO4=light-limited PEP carboxylation rate
!           CBXN4=PEP caboxylation efficiency
!           VL=PEP carboxylation rate limited by light,CO2,N,P
!           VGRO4=PEP carboxylation rate limited by CO2 from stomate.f
!           FDBK4=N,P feedback inhibition on C4 CO2 fixation
!
              PARX=QNTM*PARDIFs1(N,M,L,NZ)
              PARJ=PARX+ETGR4s1(K,NB,NZ)
              ETLF4=(PARJ-SQRT(PARJ**2-CURV4*PARX*ETGR4s1(K,NB,NZ)))/CURV2
              EGRO4=ETLF4*CBXN4s1(K,NB,NZ)
              VL=AMIN1(VGRO4s1(K,NB,NZ),EGRO4)*FDBK4s1(K,NB,NZ)
!
!             STOMATAL EFFECT OF WATER DEFICIT IN MESOPHYLL
!
!             RS,RSL=leaf stomatal resistance to CO2 at zero,current water potential
!             RCMN=minimum stomatal resistance to CO2 (s m-1)
!             RCMX=cuticular resistance to CO2 from startq.f (s m-1)
!             DCO2=difference between atmosph and intercellular CO2 concn (umol m-3)
!             GSL=leaf stomatal conductance (mol m-2 s-1)
!             WFNC=stomatal resistance function of canopy turgor
!             FMOL=number of moles of air per m3
!
              IF(VL.GT.ZEROs1)THEN
                RS=AMIN1(RCMXs1(NZ),AMAX1(RCMN,DCO2s1(NZ)/VL))
                RSL=RS+(RCMXs1(NZ)-RS)*WFNC
                GSL=1.0_r8/RSL*FMOLs1(NZ)
!
!               EFFECT OF WATER DEFICIT IN MESOPHYLL
!
!               IGTYP=growth type:0=bryophyte,1=graminoid,2=shrub,tree
!               WFN4,WFNB=non-stomatal effects of water stress on C4,C3 CO2 fixation
!
                IF(IGTYPs1(NZ).NE.0)THEN
                  WFN4=(RS/RSL)**1.00
                  WFNB=SQRT(RS/RSL)
                ELSE
                  WFN4=WFNG
                  WFNB=WFNG
                ENDIF
!
!               CONVERGENCE SOLUTION FOR CO2I AT WHICH CARBOXYLATION
!               EQUALS DIFFUSION IN MESOPHYLL
!
!               CO2I=intercellular,mesophyll CO2 concentration at zero water potential
!               CO2X,CO2C=intercellular,mesophyll CO2 concentration during convergence
!               SCO2=solubility of CO2 (uM/(umol mol-1))
!               COMP4=C4 CO2 compensation point (uM)
!               CBXNX=PEP caboxylation efficiency
!               ELEC4=e- requirement for CO2 fixn by PEP carboxylase
!               VCGR4,VGRO4=PEP carboxylation rate unlimited,limited by CO2
!               XKCO24=Km for VCMX4 from PFT file (uM)
!               EGROX=light-limited PEP carboxylation rate
!               ETLF4=light-limited e- transport rate
!               VL=PEP carboxylation rate limited by light,CO2,N,P,water stress
!               VG=CO2 diffusion rate limited by water stress
!               GSL=leaf stomatal conductance (mol m-2 s-1)
!
                CO2X=CO2Is1(NZ)
                DO 135 NN=1,100
                  CO2C=CO2X*SCO2s1(NZ)
                  CO2Y=AMAX1(0.0,CO2C-COMP4)
                  CBXNX=CO2Y/(ELEC4*CO2C+10.5*COMP4)
                  VGROX=VCGR4s1(K,NB,NZ)*CO2Y/(CO2C+XKCO24s1(NZ))
                  EGROX=ETLF4*CBXNX
                  VL=AMIN1(VGROX,EGROX)*WFN4*FDBK4s1(K,NB,NZ)
                  VG=(CO2Qs1(NZ)-CO2X)*GSL
                  IF(VL+VG.GT.ZEROs1)THEN
                    DIFF=(VL-VG)/(VL+VG)
                    IF(ABS(DIFF).LT.0.005)exit
                    VA=0.95*VG+0.05*VL
                    CO2X=CO2Qs1(NZ)-VA/GSL
                  ELSE
                    VL=0._r8
                    exit
                  ENDIF
135             CONTINUE
!
!               ACCUMULATE C4 FIXATION PRODUCT IN MESOPHYLL
!
!               CH2O4=total C4 CO2 fixation
!               SURFX=unself-shaded leaf surface area
!               TAU0s1=fraction of diffuse radiation transmitted from layer above
!
                CH2O4(K)=CH2O4(K)+VL*SURFXs1(N,L,K,NB,NZ) &
                  *TAU0s1(L+1)
!               ICO2I=MAX(1,MIN(400,INT(CO2X)))
!               VCO2(ICO2I,I,NZ)=VCO2(ICO2I,I,NZ)
!              2+(VL*SURFXs1(N,L,K,NB,NZ)*TAU0s1(L+1))*0.0432
!
!               C3 CARBOXYLATION REACTIONS IN IN BUNDLE SHEATH OF C4 PLANTS
!
!               ETGRO=light saturated e- transport rate from stomate.f
!               ETLF=light-limited e- transport rate
!               CURV=shape parameter for e- transport response to PAR
!               EGRO=light-limited rubisco carboxylation rate
!               CBXN=rubisco caboxylation efficiency
!               VL=rubisco carboxylation rate limited by light,CO2,N,P
!               VGRO=rubisco carboxylation rate limited by CO2 from stomate.f
!               FDBK=N,P feedback inhibition on C3 CO2 fixation
!
                PARJ=PARX+ETGROs1(K,NB,NZ)
                ETLF=(PARJ-SQRT(PARJ**2-CURV4*PARX*ETGROs1(K,NB,NZ)))/CURV2
                EGRO=ETLF*CBXNs1(K,NB,NZ)
                VL=AMIN1(VGROs1(K,NB,NZ),EGRO)*WFNB*FDBKs1(NB,NZ)
!
!               ACCUMULATE C3 FIXATION PRODUCT IN BUNDLE SHEATH
!
!               CH2O3=total C3 CO2 fixation
!               SURFX=unself-shaded leaf surface area
!               TAU0s1=fraction of diffuse radiation transmitted from layer above
!
                CH2O3(K)=CH2O3(K)+VL*SURFXs1(N,L,K,NB,NZ) &
                  *TAU0s1(L+1)
              ENDIF
            ENDIF
          ENDIF
120     CONTINUE
115   CONTINUE
    ENDIF
110   CONTINUE
  CO2F=CO2F+CH2O4(K)
  CH2O=CH2O+CH2O3(K)
  end associate
  end subroutine ComputeGPP_C4


!------------------------------------------------------------------------------------------

  subroutine ComputeGPP(NB,NZ,WFNG,WFNC,CH2O3,CH2O4,CH2O)
  implicit none
  integer, intent(in) :: NB,NZ
  real(r8), intent(in) :: WFNG
  real(r8), intent(in) :: WFNC
  real(r8), intent(out) :: CH2O3(JNODS1),CH2O4(JNODS1)
  real(r8), intent(out) :: CH2O
  real(r8) :: CO2F,ZADDB,PADDB

  integer :: K

! begin_execution
  associate(                           &
    CO2Qs1    =>  plt_photo%CO2Qs1   , &
    ICTYPs1   =>  plt_photo%ICTYPs1  , &
    VCGR4s1   =>  plt_photo%VCGR4s1  , &
    VCGROs1   =>  plt_photo%VCGROs1  , &
    IGTYPs1   =>  plt_pheno%IGTYPs1  , &
    SSINs1    =>  plt_rad%SSINs1     , &
    RADPs1    =>  plt_rad%RADPs1     , &
    ZEROPs1   =>  plt_biom%ZEROPs1   , &
    ARLF1s1   =>  plt_morph%ARLF1s1  , &
    FDBKs1    =>  plt_photo%FDBKs1     &
  )

  IF(abs(FDBKs1(NB,NZ)).GT.0._r8)THEN
    IF(SSINs1.GT.0.0_r8.AND.RADPs1(NZ).GT.0.0_r8.AND.CO2Qs1(NZ).GT.0.0_r8)THEN
      CO2F=0._r8
      CH2O=0._r8
      IF(IGTYPs1(NZ).NE.0.OR.WFNC.GT.0.0_r8)THEN
!
!         FOR EACH NODE
!
        DO 100 K=1,JNODS1
          CH2O3(K)=0._r8
          CH2O4(K)=0._r8
          IF(ARLF1s1(K,NB,NZ).GT.ZEROPs1(NZ))THEN
!
!             C4 PHOTOSYNTHESIS
!
!             ARLF,ARLFL=leaf area
!             ICTYP=photosynthesis type:3=C3,4=C4 from PFT file
!             VCGR4=PEP carboxylation rate unlimited by CO2
!
            IF(ICTYPs1(NZ).EQ.4.AND.VCGR4s1(K,NB,NZ).GT.0.0)THEN
!
              CALL ComputeGPP_C4(K,NB,NZ,WFNG,WFNC,CH2O3,CH2O4,CO2F,CH2O)
!
!               C3 PHOTOSYNTHESIS
!
            ELSEIF(ICTYPs1(NZ).NE.4.AND.VCGROs1(K,NB,NZ).GT.0.0)THEN
              call ComputeGPP_C3(K,NB,NZ,WFNG,WFNC,CH2O3,CO2F,CH2O)

            ENDIF
          ENDIF
100     CONTINUE
!
!         CO2F,CH2O=total CO2 fixation,CH2O production
!
        CO2F=CO2F*0.0432_r8
        CH2O=CH2O*0.0432_r8
!
!         CONVERT UMOL M-2 S-1 TO G C M-2 H-1
!
        DO 150 K=1,JNODS1
          CH2O3(K)=CH2O3(K)*0.0432
          CH2O4(K)=CH2O4(K)*0.0432
150     CONTINUE
      ELSE
        CO2F=0._r8
        CH2O=0._r8
        IF(ICTYPs1(NZ).EQ.4)THEN
          DO 155 K=1,JNODS1
            CH2O3(K)=0._r8
            CH2O4(K)=0._r8
155       CONTINUE
        ENDIF
      ENDIF
    ELSE
      CO2F=0._r8
      CH2O=0._r8
      IF(ICTYPs1(NZ).EQ.4)THEN
        DO 160 K=1,JNODS1
          CH2O3(K)=0._r8
          CH2O4(K)=0._r8
160     CONTINUE
      ENDIF
    ENDIF
  ELSE
    CO2F=0._r8
    CH2O=0._r8
    IF(ICTYPs1(NZ).EQ.4)THEN
      DO 165 K=1,JNODS1
        CH2O3(K)=0._r8
        CH2O4(K)=0._r8
165   CONTINUE
    ENDIF
  ENDIF
  end associate
  end subroutine ComputeGPP
end module PhotoSynsMod
