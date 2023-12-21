module PhotoSynsMod
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use EcosimConst
  use GrosubPars
  use minimathmod, only : AZMAX1
  use PlantMathFuncMod
  use PlantAPIData
implicit none
  private
  character(len=*),private, parameter :: mod_filename = &
  __FILE__
  public :: ComputeGPP

  contains

!------------------------------------------------------------------------------------------

  subroutine ComputeGPP_C3(K,NB,NZ,WFNG,Stomata_Activity,CH2O3K)

  implicit none
  integer, intent(in) :: K,NB,NZ
  real(r8), intent(in) :: WFNG,Stomata_Activity
  real(r8), intent(out) :: CH2O3K
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
  ZERO                         => plt_site%ZERO   , &
  iPlantMorphologyType_pft     => plt_pheno%iPlantMorphologyType_pft, &
  RubiscoActivity_brch         => plt_photo%RubiscoActivity_brch  , &
  CanopyGasCO2_pft             => plt_photo%CanopyGasCO2_pft  , &
  LeafAUnshaded_zsec           => plt_photo%LeafAUnshaded_zsec , &
  Km4RubiscoCarboxy_pft        => plt_photo%Km4RubiscoCarboxy_pft, &
  CO2Solubility_pft            => plt_photo%CO2Solubility_pft , &
  RubiscoCarboxyEff_node       => plt_photo%RubiscoCarboxyEff_node , &
  LigthSatCarboxyRate_node     => plt_photo%LigthSatCarboxyRate_node , &
  CO2lmtRubiscoCarboxyRate_node=> plt_photo%CO2lmtRubiscoCarboxyRate_node , &
  LeafIntracellularCO2_pft     => plt_photo%LeafIntracellularCO2_pft , &
  CO2CuticleResist_pft         => plt_photo%CO2CuticleResist_pft  , &
  Km4PEPCarboxy_pft            => plt_photo%Km4PEPCarboxy_pft, &
  Vmax4RubiscoCarboxy_pft      => plt_photo%Vmax4RubiscoCarboxy_pft , &
  CO2CompenPoint_node          => plt_photo%CO2CompenPoint_node , &
  AirConc_pft                  => plt_photo%AirConc_pft  , &
  DiffCO2Atmos2Intracel_pft    => plt_photo%DiffCO2Atmos2Intracel_pft , &
  ZEROP                        => plt_biom%ZEROP  , &
  CanopyLeafAreaByLayer_pft    => plt_morph%CanopyLeafAreaByLayer_pft , &
  PARDiffus_zsec               => plt_rad%PARDiffus_zsec  , &
  PARDirect_zsec               => plt_rad%PARDirect_zsec    , &
  TAU_RadThru                  => plt_rad%TAU_RadThru    , &
  TAU_RadCapt                  => plt_rad%TAU_RadCapt      &
  )

  CH2O3K=0._r8
! FOR EACH CANOPY LAYER
!
  D210: DO L=NumOfCanopyLayers1,1,-1
    IF(CanopyLeafAreaByLayer_pft(L,K,NB,NZ).GT.ZEROP(NZ))THEN
!
!     FOR EACH LEAF AZIMUTH AND INCLINATION
!
      D215: DO N=1,NumOfLeafZenithSectors1
        D220: DO M=1,NumOfSkyAzimuSects1
!
!         CO2 FIXATION BY SUNLIT LEAVES
!
!         LeafAUnshaded_zsec=unself-shaded leaf surface area
!
          IF(LeafAUnshaded_zsec(N,L,K,NB,NZ).GT.ZEROP(NZ))THEN
            IF(PARDirect_zsec(N,M,L,NZ).GT.0.0_r8)THEN
!
!             C3 CARBOXYLATION REACTIONS IN MESOPHYLL
!
!             QNTM=quantum efficiency
!             PAR=direct PAR flux
!             LigthSatCarboxyRate_node=light saturated e- transport rate from stomate.f
!             ETLF=light-limited e- transport rate
!             CURV=shape parameter for e- transport response to PAR
!             EGRO=light-limited rubisco carboxylation rate
!             RubiscoCarboxyEff_node=rubisco caboxylation efficiency
!             VL=rubisco carboxylation rate limited by light,CO2,N,P
!             CO2lmtRubiscoCarboxyRate_node=rubisco carboxylation rate limited by CO2 from stomate.f
!             RubiscoActivity_brch=N,P feedback inhibition on C4 CO2 fixation
!
              PARX=QNTM*PARDirect_zsec(N,M,L,NZ)
              PARJ=PARX+LigthSatCarboxyRate_node(K,NB,NZ)
              ETLF=(PARJ-SQRT(PARJ**2-CURV4*PARX*LigthSatCarboxyRate_node(K,NB,NZ)))/CURV2
              EGRO=ETLF*RubiscoCarboxyEff_node(K,NB,NZ)
              VL=AMIN1(CO2lmtRubiscoCarboxyRate_node(K,NB,NZ),EGRO)*RubiscoActivity_brch(NB,NZ)
!
!             STOMATAL EFFECT OF WATER DEFICIT IN MESOPHYLL
!
!             RS,RSL=leaf stomatal resistance to CO2 at zero,current water potential
!             RCMN=minimum stomatal resistance to CO2 (s m-1)
!             CO2CuticleResist_pft=cuticular resistance to CO2 from startq.f (s m-1)
!             DiffCO2Atmos2Intracel_pft=difference between atmosph and intercellular CO2 concn (umol m-3)
!             GSL=leaf stomatal conductance (mol m-2 s-1)
!             Stomata_Activity=stomatal resistance function of canopy turgor
!             AirConc_pft=number of moles of air per m3
!
              IF(VL.GT.ZERO)THEN
                RS=AMIN1(CO2CuticleResist_pft(NZ),AMAX1(RCMN,DiffCO2Atmos2Intracel_pft(NZ)/VL))
                RSL=RS+(CO2CuticleResist_pft(NZ)-RS)*Stomata_Activity
                GSL=1.0_r8/RSL*AirConc_pft(NZ)
!
!               EFFECT OF WATER DEFICIT IN MESOPHYLL
!
!               iPlantMorphologyType_pft=growth type:0=bryophyte,1=graminoid,2=shrub,tree
!               WFNB=non-stomatal effects of water stress on CO2 fixation
!
                IF(.not.is_plant_bryophyte(iPlantMorphologyType_pft(NZ)))THEN
                  WFNB=SQRT(RS/RSL)
                ELSE
                  WFNB=WFNG
                ENDIF
!
!               CONVERGENCE SOLUTION FOR LeafIntracellularCO2_pftAT WHICH CARBOXYLATION
!               EQUALS DIFFUSION IN MESOPHYLL
!
!               LeafIntracellularCO2_pft=intercellular,mesophyll CO2 concentration at zero water potential
!               CO2X,CO2C=intercellular,mesophyll CO2 concentration during convergence
!               CO2Solubility_pft=solubility of CO2 (uM/(umol mol-1))
!               CO2CompenPoint_node=C3 CO2 compensation point (uM)
!               CBXNX=rubisco carboxylation efficiency
!               ELEC3=e- requirement for CO2 fixn by rubisco
!               Vmax4RubiscoCarboxy_pft,VGROX=rubisco carboxylation rate unlimited,limited by CO2
!               Km4RubiscoCarboxy_pft=Km for rubisco carboxylation
!               EGROX=light-limited rubisco carboxylation rate
!               ETLF=light-limited e- transport rate
!               VL=rubisco carboxylation rate limited by light,CO2,N,P,water stress
!               VG=CO2 diffusion rate limited by water stress
!               GSL=leaf stomatal conductance (mol m-2 s-1)
!
                CO2X=LeafIntracellularCO2_pft(NZ)
                D225: DO NN=1,100
                  CO2C=CO2X*CO2Solubility_pft(NZ)
                  CO2Y=AZMAX1(CO2C-CO2CompenPoint_node(K,NB,NZ))
                  CBXNX=CO2Y/(ELEC3*CO2C+10.5_r8*CO2CompenPoint_node(K,NB,NZ))
                  VGROX=Vmax4RubiscoCarboxy_pft(K,NB,NZ)*CO2Y/(CO2C+Km4RubiscoCarboxy_pft(NZ))
                  EGROX=ETLF*CBXNX
                  VL=AMIN1(VGROX,EGROX)*WFNB*RubiscoActivity_brch(NB,NZ)
                  VG=(CanopyGasCO2_pft(NZ)-CO2X)*GSL
                  
                  IF(VL+VG.GT.ZERO)THEN
                    DIFF=(VL-VG)/(VL+VG)
                    IF(ABS(DIFF).LT.0.005_r8)exit
                    VA=0.95_r8*VG+0.05_r8*VL
                    CO2X=CanopyGasCO2_pft(NZ)-VA/GSL
                  ELSE
                    VL=0._r8
                    exit
                  ENDIF
                ENDDO D225
!
!               ACCUMULATE C3 FIXATION PRODUCT IN MESOPHYLL
!
!               CH2O3=total C3 CO2 fixation
!               LeafAUnshaded_zsec=unself-shaded leaf surface area
!               TAU_RadCapt=fraction of direct radiation transmitted from layer above
!
                CH2O3K=CH2O3K+VL*LeafAUnshaded_zsec(N,L,K,NB,NZ)*TAU_RadCapt(L+1)
!               ICO2I=MAX(1,MIN(400,INT(CO2X)))
!               VCO2(ICO2I,I,NZ)=VCO2(ICO2I,I,NZ)
!              2+(VL*LeafAUnshaded_zsec(N,L,K,NB,NZ)*TAU_RadCapt(L+1))*0.0432

              ENDIF
            ENDIF
!
!           CO2 FIXATION IN MESOPHYLL BY SHADED LEAVES
!
            IF(PARDiffus_zsec(N,M,L,NZ).GT.0.0_r8)THEN
!
!             C3 CARBOXYLATION REACTIONS USING VARIABLES FROM 'STOMATE'
!
!             QNTM=quantum efficiency
!             PARDiffus_zsec=diffuse PAR flux
!             LigthSatCarboxyRate_node=light saturated e- transport rate from stomate.f
!             ETLF=light-limited e- transport rate
!             CURV=shape parameter for e- transport response to PAR
!             EGRO=light-limited rubisco carboxylation rate
!             RubiscoCarboxyEff_node=rubisco caboxylation efficiency
!             VL=rubisco carboxylation rate limited by light,CO2,N,P
!             CO2lmtRubiscoCarboxyRate_node=rubisco carboxylation rate limited by CO2 from stomate.f
!             RubiscoActivity_brch=N,P feedback inhibition on C3 CO2 fixation
!
              PARX=QNTM*PARDiffus_zsec(N,M,L,NZ)
              PARJ=PARX+LigthSatCarboxyRate_node(K,NB,NZ)
              ETLF=(PARJ-SQRT(PARJ**2-CURV4*PARX*LigthSatCarboxyRate_node(K,NB,NZ)))/CURV2
              EGRO=ETLF*RubiscoCarboxyEff_node(K,NB,NZ)
              VL=AMIN1(CO2lmtRubiscoCarboxyRate_node(K,NB,NZ),EGRO)*RubiscoActivity_brch(NB,NZ)
!
!             STOMATAL EFFECT OF WATER DEFICIT IN MESOPHYLL
!
!             RS,RSL=leaf stomatal resistance to CO2 at zero,current water potential
!             RCMN=minimum stomatal resistance to CO2 (s m-1)
!             CO2CuticleResist_pft=cuticular resistance to CO2 from startq.f (s m-1)
!             DiffCO2Atmos2Intracel_pft=difference between atmosph and intercellular CO2 concn (umol m-3)
!             GSL=leaf stomatal conductance (mol m-2 s-1)
!             Stomata_Activity=stomatal resistance function of canopy turgor
!             AirConc_pft=number of moles of air per m3
!
              IF(VL.GT.ZERO)THEN
                RS=AMIN1(CO2CuticleResist_pft(NZ),AMAX1(RCMN,DiffCO2Atmos2Intracel_pft(NZ)/VL))
                RSL=RS+(CO2CuticleResist_pft(NZ)-RS)*Stomata_Activity
                GSL=1.0_r8/RSL*AirConc_pft(NZ)
!
!               EFFECT OF WATER DEFICIT IN MESOPHYLL
!
!               iPlantMorphologyType_pft=growth type:0=bryophyte,1=graminoid,2=shrub,tree
!               WFNB=non-stomatal effects of water stress on C3 CO2 fixation
!
                IF(.not.is_plant_bryophyte(iPlantMorphologyType_pft(NZ)))THEN
                  WFNB=SQRT(RS/RSL)
                ELSE
                  WFNB=WFNG
                ENDIF
!
!               CONVERGENCE SOLUTION FOR LeafIntracellularCO2_pftAT WHICH CARBOXYLATION
!               EQUALS DIFFUSION IN MESOPHYLL
!
!               LeafIntracellularCO2_pft=intercellular,mesophyll CO2 concentration at zero water potential
!               CO2X,CO2C=intercellular,mesophyll CO2 concentration during convergence
!               CO2Solubility_pft=solubility of CO2 (uM/(umol mol-1))
!               CO2CompenPoint_node=C3 CO2 compensation point (uM)
!               CBXNX=rubisco caboxylation efficiency
!               ELEC3=e- requirement for CO2 fixn by rubisco carboxylase
!               Vmax4RubiscoCarboxy_pft,VGROX=rubisco carboxylation rate unlimited,limited by CO2
!               Km4RubiscoCarboxy_pft=Km for rubisco carboxylation from stomate.f (uM)
!               EGROX=light-limited rubisco carboxylation rate
!               ETLF=light-limited e- transport rate
!               VL=rubisco carboxylation rate limited by light,CO2,N,P,water stress
!               VG=CO2 diffusion rate limited by water stress
!               GSL=leaf stomatal conductance (mol m-2 s-1)
!
                CO2X=LeafIntracellularCO2_pft(NZ)
                D235: DO NN=1,100
                  CO2C=CO2X*CO2Solubility_pft(NZ)
                  CO2Y=AZMAX1(CO2C-CO2CompenPoint_node(K,NB,NZ))
                  CBXNX=CO2Y/(ELEC3*CO2C+10.5_r8*CO2CompenPoint_node(K,NB,NZ))
                  VGROX=Vmax4RubiscoCarboxy_pft(K,NB,NZ)*CO2Y/(CO2C+Km4RubiscoCarboxy_pft(NZ))
                  EGROX=ETLF*CBXNX
                  VL=AMIN1(VGROX,EGROX)*WFNB*RubiscoActivity_brch(NB,NZ)
                  VG=(CanopyGasCO2_pft(NZ)-CO2X)*GSL
                  IF(VL+VG.GT.ZERO)THEN
                    DIFF=(VL-VG)/(VL+VG)
                    IF(ABS(DIFF).LT.0.005_r8)exit
                    VA=0.95_r8*VG+0.05_r8*VL
                    CO2X=CanopyGasCO2_pft(NZ)-VA/GSL
                  ELSE
                    VL=0._r8
                    exit
                  ENDIF
                ENDDO D235
!
!               ACCUMULATE C3 FIXATION PRODUCT IN MESOPHYLL
!
!               CH2O3=total C3 CO2 fixation
!               LeafAUnshaded_zsec=unself-shaded leaf surface area
!               TAU_RadThru=fraction of diffuse radiation transmitted from layer above
!
                CH2O3K=CH2O3K+VL*LeafAUnshaded_zsec(N,L,K,NB,NZ)*TAU_RadThru(L+1)
!               ICO2I=MAX(1,MIN(400,INT(CO2X)))
!               VCO2(ICO2I,I,NZ)=VCO2(ICO2I,I,NZ)
!              2+(VL*LeafAUnshaded_zsec(N,L,K,NB,NZ)*TAU_RadThru(L+1))*0.0432
              ENDIF
            ENDIF
          ENDIF
        ENDDO D220
      ENDDO D215
    ENDIF
  ENDDO D210
  end associate
  end subroutine ComputeGPP_C3


!------------------------------------------------------------------------------------------

  subroutine ComputeGPP_C4(K,NB,NZ,WFNG,Stomata_Activity,CH2O3K,CH2O4K)
  implicit none
  integer, intent(in) :: K,NB,NZ
  real(r8), intent(in):: WFNG,Stomata_Activity
  real(r8), intent(out) :: CH2O3K,CH2O4K
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
  iPlantMorphologyType_pft          => plt_pheno%iPlantMorphologyType_pft, &
  ZEROP                             => plt_biom%ZEROP  , &
  Km4PEPCarboxy_pft                 => plt_photo%Km4PEPCarboxy_pft, &
  NutrientCtrlonC4Carboxy_node      => plt_photo%NutrientCtrlonC4Carboxy_node , &
  CO2Solubility_pft                 => plt_photo%CO2Solubility_pft , &
  C4CarboxyEff_node                 => plt_photo%C4CarboxyEff_node , &
  LeafAUnshaded_zsec                => plt_photo%LeafAUnshaded_zsec , &
  CO2lmtPEPCarboxyRate_node         => plt_photo%CO2lmtPEPCarboxyRate_node, &
  CO2CuticleResist_pft              => plt_photo%CO2CuticleResist_pft  , &
  DiffCO2Atmos2Intracel_pft         => plt_photo%DiffCO2Atmos2Intracel_pft , &
  AirConc_pft                       => plt_photo%AirConc_pft  , &
  LigthSatCarboxyRate_node          => plt_photo%LigthSatCarboxyRate_node , &
  CO2lmtRubiscoCarboxyRate_node     => plt_photo%CO2lmtRubiscoCarboxyRate_node , &
  LigthSatC4CarboxyRate_node        => plt_photo%LigthSatC4CarboxyRate_node , &
  CanopyGasCO2_pft                  => plt_photo%CanopyGasCO2_pft  , &
  RubiscoActivity_brch              => plt_photo%RubiscoActivity_brch  , &
  LeafIntracellularCO2_pft          => plt_photo%LeafIntracellularCO2_pft , &
  RubiscoCarboxyEff_node            => plt_photo%RubiscoCarboxyEff_node , &
  Vmax4PEPCarboxy_pft               => plt_photo%Vmax4PEPCarboxy_pft , &
  CanopyLeafAreaByLayer_pft         => plt_morph%CanopyLeafAreaByLayer_pft , &
  ZERO                              => plt_site%ZERO   , &
  PARDiffus_zsec                    => plt_rad%PARDiffus_zsec  , &
  PARDirect_zsec                    => plt_rad%PARDirect_zsec    , &
  TAU_RadThru                       => plt_rad%TAU_RadThru    , &
  TAU_RadCapt                       => plt_rad%TAU_RadCapt      &
  )

  CH2O3K=0._r8;CH2O4K=0._r8
! FOR EACH CANOPY LAYER
!
  D110: DO L=NumOfCanopyLayers1,1,-1
    IF(CanopyLeafAreaByLayer_pft(L,K,NB,NZ).GT.ZEROP(NZ))THEN
!
!     FOR EACH LEAF AZIMUTH AND INCLINATION
!
      D115: DO N =1,NumOfLeafZenithSectors1
        D120: DO M =1,NumOfSkyAzimuSects1
!
!         CO2 FIXATION IN MESOPHYLL BY SUNLIT LEAVES
!
!         LeafAUnshaded_zsec=unself-shaded leaf surface area
!
          IF(LeafAUnshaded_zsec(N,L,K,NB,NZ).GT.ZEROP(NZ))THEN
            IF(PARDirect_zsec(N,M,L,NZ).GT.0.0)THEN
!
!             C4 CARBOXYLATION REACTIONS IN MESOPHYLL
!
!             QNTM=quantum efficiency
!             PAR=direct PAR flux
!             LigthSatC4CarboxyRate_node=light saturated e- transport rate from stomate.f
!             ETLF4=light-limited e- transport rate
!             CURV=shape parameter for e- transport response to PAR
!             EGRO4=light-limited PEP carboxylation rate
!             C4CarboxyEff_node=PEP caboxylation efficiency
!             VL=PEP carboxylation rate limited by light,CO2,N,P
!             CO2lmtPEPCarboxyRate_node=PEP carboxylation rate limited by CO2 from stomate.f
!             NutrientCtrlonC4Carboxy_node=N,P feedback inhibition on C4 CO2 fixation
!
              PARX=QNTM*PARDirect_zsec(N,M,L,NZ)
              PARJ=PARX+LigthSatC4CarboxyRate_node(K,NB,NZ)
              ETLF4=(PARJ-SQRT(PARJ**2-CURV4*PARX*LigthSatC4CarboxyRate_node(K,NB,NZ)))/CURV2
              EGRO4=ETLF4*C4CarboxyEff_node(K,NB,NZ)
              VL=AMIN1(CO2lmtPEPCarboxyRate_node(K,NB,NZ),EGRO4)*NutrientCtrlonC4Carboxy_node(K,NB,NZ)
!
!             STOMATAL EFFECT OF WATER DEFICIT IN MESOPHYLL
!
!             RS,RSL=leaf stomatal resistance to CO2 at zero,current water potential
!             RCMN=minimum stomatal resistance to CO2 (s m-1)
!             CO2CuticleResist_pft=cuticular resistance to CO2 from startq.f (s m-1)
!             DiffCO2Atmos2Intracel_pft=difference between atmosph and intercellular CO2 concn (umol m-3)
!             GSL=leaf stomatal conductance (mol m-2 s-1)
!             Stomata_Activity=stomatal resistance function of canopy turgor
!             AirConc_pft=number of moles of air per m3
!
              IF(VL.GT.ZERO)THEN
                RS=AMIN1(CO2CuticleResist_pft(NZ),AMAX1(RCMN,DiffCO2Atmos2Intracel_pft(NZ)/VL))
                RSL=RS+(CO2CuticleResist_pft(NZ)-RS)*Stomata_Activity
                GSL=1.0_r8/RSL*AirConc_pft(NZ)
!
!               EFFECT OF WATER DEFICIT IN MESOPHYLL
!
!               iPlantMorphologyType_pft=growth type:0=bryophyte,1=graminoid,2=shrub,tree
!               WFNB=non-stomatal effects of water stress on C4,C3 CO2 fixation
!
                IF(.not.is_plant_bryophyte(iPlantMorphologyType_pft(NZ)))THEN
                  WFN4=RS/RSL
                  WFNB=SQRT(RS/RSL)
                ELSE
                  WFN4=WFNG
                  WFNB=WFNG
                ENDIF
!
!               CONVERGENCE SOLUTION FOR LeafIntracellularCO2_pftAT WHICH CARBOXYLATION
!               EQUALS DIFFUSION IN MESOPHYLL
!
!               LeafIntracellularCO2_pft=intercellular,mesophyll CO2 concentration at zero water potential
!               CO2X,CO2C=intercellular,mesophyll CO2 concentration during convergence
!               CO2Solubility_pft=solubility of CO2 (uM/(umol mol-1))
!               COMP4=C4 CO2 compensation point (uM)
!               CBXNX=PEP carboxylation efficiency
!               ELEC4=e- requirement for CO2 fixn by PEP carboxylase
!               Vmax4PEPCarboxy_pft,VGROX=PEP carboxylation rate unlimited,limited by CO2
!               Km4PEPCarboxy_pft=Km for VmaxPEPCarboxyRef_pft from PFT file (uM)
!               EGROX=light-limited PEP carboxylation rate
!               ETLF4=light-limited e- transport rate
!               VL=PEP carboxylation rate limited by light,CO2,N,P,water stress
!               VG=CO2 diffusion rate limited by water stress
!               GSL=leaf stomatal conductance (mol m-2 s-1)
!
                CO2X=LeafIntracellularCO2_pft(NZ)
                D125: DO NN=1,100
                  CO2C=CO2X*CO2Solubility_pft(NZ)
                  CO2Y=AZMAX1(CO2C-COMP4)
                  CBXNX=CO2Y/(ELEC4*CO2C+10.5_r8*COMP4)
                  VGROX=Vmax4PEPCarboxy_pft(K,NB,NZ)*CO2Y/(CO2C+Km4PEPCarboxy_pft(NZ))
                  EGROX=ETLF4*CBXNX
                  VL=AMIN1(VGROX,EGROX)*WFN4*NutrientCtrlonC4Carboxy_node(K,NB,NZ)
                  VG=(CanopyGasCO2_pft(NZ)-CO2X)*GSL
                  IF(VL+VG.GT.ZERO)THEN
                    DIFF=(VL-VG)/(VL+VG)
                    IF(ABS(DIFF).LT.0.005_r8)exit
                    VA=0.95_r8*VG+0.05_r8*VL
                    CO2X=CanopyGasCO2_pft(NZ)-VA/GSL
                  ELSE
                    VL=0._r8
                    exit
                  ENDIF
                ENDDO D125
!
!               ACCUMULATE C4 FIXATION PRODUCT IN MESOPHYLL
!
!               CH2O4=total C4 CO2 fixation
!               LeafAUnshaded_zsec=unself-shaded leaf surface area
!               TAU_RadCapt=fraction of direct radiation transmitted from layer above
!
                CH2O4K=CH2O4K+VL*LeafAUnshaded_zsec(N,L,K,NB,NZ) &
                  *TAU_RadCapt(L+1)
!               ICO2I=MAX(1,MIN(400,INT(CO2X)))
!               VCO2(ICO2I,I,NZ)=VCO2(ICO2I,I,NZ)
!              2+(VL*LeafAUnshaded_zsec(N,L,K,NB,NZ)*TAU_RadCapt(L+1))*0.0432
!
!               C3 CARBOXYLATION REACTIONS IN BUNDLE SHEATH OF C4 PLANTS
!
!               LigthSatCarboxyRate_node=light saturated e- transport rate from stomate.f
!               ETLF=light-limited e- transport rate
!               CURV=shape parameter for e- transport response to PAR
!               EGRO=light-limited rubisco carboxylation rate
!               RubiscoCarboxyEff_node=rubisco caboxylation efficiency
!               VL=rubisco carboxylation rate limited by light,CO2,N,P
!               CO2lmtRubiscoCarboxyRate_node=rubisco carboxylation rate limited by CO2 from stomate.f
!               RubiscoActivity_brch=N,P feedback inhibition on C3 CO2 fixation
!
                PARJ=PARX+LigthSatCarboxyRate_node(K,NB,NZ)
                ETLF=(PARJ-SQRT(PARJ**2-CURV4*PARX*LigthSatCarboxyRate_node(K,NB,NZ)))/CURV2
                EGRO=ETLF*RubiscoCarboxyEff_node(K,NB,NZ)
                VL=AMIN1(CO2lmtRubiscoCarboxyRate_node(K,NB,NZ),EGRO)*WFNB*RubiscoActivity_brch(NB,NZ)
!
!               ACCUMULATE C3 FIXATION PRODUCT IN BUNDLE SHEATH
!
!               CH2O3=total C3 CO2 fixation
!               LeafAUnshaded_zsec=unself-shaded leaf surface area
!               TAU_RadCapt=fraction of direct radiation transmitted from layer above
!
                CH2O3K=CH2O3K+VL*LeafAUnshaded_zsec(N,L,K,NB,NZ)*TAU_RadCapt(L+1)

              ENDIF
            ENDIF
!
!           CO2 FIXATION IN MESOPHYLL BY SHADED LEAVES
!
            IF(PARDiffus_zsec(N,M,L,NZ).GT.0.0)THEN
!
!           C4 CARBOXYLATION REACTIONS IN MESOPHYLL
!
!           QNTM=quantum efficiency
!           PARDiffus_zsec=diffuse PAR flux
!           LigthSatC4CarboxyRate_node=light saturated e- transport rate from stomate.f
!           ETLF4=light-limited e- transport rate
!           CURV=shape parameter for e- transport response to PAR
!           EGRO4=light-limited PEP carboxylation rate
!           C4CarboxyEff_node=PEP caboxylation efficiency
!           VL=PEP carboxylation rate limited by light,CO2,N,P
!           CO2lmtPEPCarboxyRate_node=PEP carboxylation rate limited by CO2 from stomate.f
!           NutrientCtrlonC4Carboxy_node=N,P feedback inhibition on C4 CO2 fixation
!
              PARX=QNTM*PARDiffus_zsec(N,M,L,NZ)
              PARJ=PARX+LigthSatC4CarboxyRate_node(K,NB,NZ)
              ETLF4=(PARJ-SQRT(PARJ**2-CURV4*PARX*LigthSatC4CarboxyRate_node(K,NB,NZ)))/CURV2
              EGRO4=ETLF4*C4CarboxyEff_node(K,NB,NZ)
              VL=AMIN1(CO2lmtPEPCarboxyRate_node(K,NB,NZ),EGRO4)*NutrientCtrlonC4Carboxy_node(K,NB,NZ)
!
!             STOMATAL EFFECT OF WATER DEFICIT IN MESOPHYLL
!
!             RS,RSL=leaf stomatal resistance to CO2 at zero,current water potential
!             RCMN=minimum stomatal resistance to CO2 (s m-1)
!             CO2CuticleResist_pft=cuticular resistance to CO2 from startq.f (s m-1)
!             DiffCO2Atmos2Intracel_pft=difference between atmosph and intercellular CO2 concn (umol m-3)
!             GSL=leaf stomatal conductance (mol m-2 s-1)
!             Stomata_Activity=stomatal resistance function of canopy turgor
!             AirConc_pft=number of moles of air per m3
!
              IF(VL.GT.ZERO)THEN
                RS=AMIN1(CO2CuticleResist_pft(NZ),AMAX1(RCMN,DiffCO2Atmos2Intracel_pft(NZ)/VL))
                RSL=RS+(CO2CuticleResist_pft(NZ)-RS)*Stomata_Activity
                GSL=1.0_r8/RSL*AirConc_pft(NZ)
!
!               EFFECT OF WATER DEFICIT IN MESOPHYLL
!
!               iPlantMorphologyType_pft=growth type:0=bryophyte,1=graminoid,2=shrub,tree
!               WFN4,WFNB=non-stomatal effects of water stress on C4,C3 CO2 fixation
!
                IF(.not.is_plant_bryophyte(iPlantMorphologyType_pft(NZ)))THEN
                  WFN4=(RS/RSL)**1.00_r8
                  WFNB=SQRT(RS/RSL)
                ELSE
                  WFN4=WFNG
                  WFNB=WFNG
                ENDIF
!
!               CONVERGENCE SOLUTION FOR LeafIntracellularCO2_pftAT WHICH CARBOXYLATION
!               EQUALS DIFFUSION IN MESOPHYLL
!
!               LeafIntracellularCO2_pft=intercellular,mesophyll CO2 concentration at zero water potential
!               CO2X,CO2C=intercellular,mesophyll CO2 concentration during convergence
!               CO2Solubility_pft=solubility of CO2 (uM/(umol mol-1))
!               COMP4=C4 CO2 compensation point (uM)
!               CBXNX=PEP caboxylation efficiency
!               ELEC4=e- requirement for CO2 fixn by PEP carboxylase
!               Vmax4PEPCarboxy_pft,CO2lmtPEPCarboxyRate_node=PEP carboxylation rate unlimited,limited by CO2
!               Km4PEPCarboxy_pft=Km for VmaxPEPCarboxyRef_pft from PFT file (uM)
!               EGROX=light-limited PEP carboxylation rate
!               ETLF4=light-limited e- transport rate
!               VL=PEP carboxylation rate limited by light,CO2,N,P,water stress
!               VG=CO2 diffusion rate limited by water stress
!               GSL=leaf stomatal conductance (mol m-2 s-1)
!
                CO2X=LeafIntracellularCO2_pft(NZ)
                D135: DO NN=1,100
                  CO2C=CO2X*CO2Solubility_pft(NZ)
                  CO2Y=AZMAX1(CO2C-COMP4)
                  CBXNX=CO2Y/(ELEC4*CO2C+10.5_r8*COMP4)
                  VGROX=Vmax4PEPCarboxy_pft(K,NB,NZ)*CO2Y/(CO2C+Km4PEPCarboxy_pft(NZ))
                  EGROX=ETLF4*CBXNX
                  VL=AMIN1(VGROX,EGROX)*WFN4*NutrientCtrlonC4Carboxy_node(K,NB,NZ)
                  VG=(CanopyGasCO2_pft(NZ)-CO2X)*GSL
                  IF(VL+VG.GT.ZERO)THEN
                    DIFF=(VL-VG)/(VL+VG)
                    !exit evaluation
                    IF(ABS(DIFF).LT.0.005_r8)exit
                    VA=0.95_r8*VG+0.05_r8*VL
                    CO2X=CanopyGasCO2_pft(NZ)-VA/GSL
                  ELSE
                    VL=0._r8
                    exit
                  ENDIF
                ENDDO D135
!
!               ACCUMULATE C4 FIXATION PRODUCT IN MESOPHYLL
!
!               CH2O4=total C4 CO2 fixation
!               LeafAUnshaded_zsec=unself-shaded leaf surface area
!               TAU_RadThru=fraction of diffuse radiation transmitted from layer above
!
                CH2O4K=CH2O4K+VL*LeafAUnshaded_zsec(N,L,K,NB,NZ)*TAU_RadThru(L+1)
!               ICO2I=MAX(1,MIN(400,INT(CO2X)))
!               VCO2(ICO2I,I,NZ)=VCO2(ICO2I,I,NZ)
!              2+(VL*LeafAUnshaded_zsec(N,L,K,NB,NZ)*TAU_RadThru(L+1))*0.0432
!
!               C3 CARBOXYLATION REACTIONS IN IN BUNDLE SHEATH OF C4 PLANTS
!
!               LigthSatCarboxyRate_node=light saturated e- transport rate from stomate.f
!               ETLF=light-limited e- transport rate
!               CURV=shape parameter for e- transport response to PAR
!               EGRO=light-limited rubisco carboxylation rate
!               RubiscoCarboxyEff_node=rubisco caboxylation efficiency
!               VL=rubisco carboxylation rate limited by light,CO2,N,P
!               CO2lmtRubiscoCarboxyRate_node=rubisco carboxylation rate limited by CO2 from stomate.f
!               RubiscoActivity_brch=N,P feedback inhibition on C3 CO2 fixation
!
                PARJ=PARX+LigthSatCarboxyRate_node(K,NB,NZ)
                ETLF=(PARJ-SQRT(PARJ**2-CURV4*PARX*LigthSatCarboxyRate_node(K,NB,NZ)))/CURV2
                EGRO=ETLF*RubiscoCarboxyEff_node(K,NB,NZ)
                VL=AMIN1(CO2lmtRubiscoCarboxyRate_node(K,NB,NZ),EGRO)*WFNB*RubiscoActivity_brch(NB,NZ)
!
!               ACCUMULATE C3 FIXATION PRODUCT IN BUNDLE SHEATH
!
!               CH2O3=total C3 CO2 fixation
!               LeafAUnshaded_zsec=unself-shaded leaf surface area
!               TAU_RadThru=fraction of diffuse radiation transmitted from layer above
!
                CH2O3K=CH2O3K+VL*LeafAUnshaded_zsec(N,L,K,NB,NZ)*TAU_RadThru(L+1)
              ENDIF
            ENDIF
          ENDIF
        ENDDO D120
      ENDDO D115
    ENDIF
  ENDDO D110
  end associate
  end subroutine ComputeGPP_C4


!------------------------------------------------------------------------------------------

  subroutine ComputeGPP(NB,NZ,WFNG,Stomata_Activity,CH2O3,CH2O4,CH2O,CO2F)
  implicit none
  integer, intent(in) :: NB,NZ
  real(r8), intent(in) :: WFNG
  real(r8), intent(in) :: Stomata_Activity    !between 0. and 1., a function of canopy turgor
  real(r8), intent(out) :: CH2O3(MaxNodesPerBranch1),CH2O4(MaxNodesPerBranch1)
  real(r8), intent(out) :: CO2F,CH2O
  real(r8) :: ZADDB,PADDB

  integer :: K

! begin_execution
  associate(                                                     &
    CanopyGasCO2_pft          =>  plt_photo%CanopyGasCO2_pft   , &
    iPlantPhotosynthesisType  =>  plt_photo%iPlantPhotosynthesisType , &
    Vmax4PEPCarboxy_pft       =>  plt_photo%Vmax4PEPCarboxy_pft  , &
    Vmax4RubiscoCarboxy_pft   =>  plt_photo%Vmax4RubiscoCarboxy_pft  , &
    iPlantMorphologyType_pft  =>  plt_pheno%iPlantMorphologyType_pft , &
    SineSolarIncliAngle            =>  plt_rad%SineSolarIncliAngle     , &
    RadPARbyCanopy_pft        =>  plt_rad%RadPARbyCanopy_pft     , &
    ZEROP                     =>  plt_biom%ZEROP   , &
    LeafAreaNode_brch         =>  plt_morph%LeafAreaNode_brch  , &
    RubiscoActivity_brch     =>  plt_photo%RubiscoActivity_brch     &
  )
!  write(193,*)'rubisco',NB,NZ,RubiscoActivity_brch(NB,NZ)
  IF(abs(RubiscoActivity_brch(NB,NZ)).GT.0._r8)THEN
    IF(SineSolarIncliAngle.GT.0.0_r8.AND.RadPARbyCanopy_pft(NZ).GT.0.0_r8 &
      .AND.CanopyGasCO2_pft(NZ).GT.0.0_r8)THEN
      CO2F=0._r8
      CH2O=0._r8
      IF(.not.is_plant_bryophyte(iPlantMorphologyType_pft(NZ)).OR.Stomata_Activity.GT.0.0_r8)THEN
!
!         FOR EACH NODE
!
        D100: DO K=1,MaxNodesPerBranch1
          CH2O3(K)=0._r8
          CH2O4(K)=0._r8
!          write(192,*)'leafarea',K,NB,NZ,LeafAreaNode_brch(K,NB,NZ)
          IF(LeafAreaNode_brch(K,NB,NZ).GT.ZEROP(NZ))THEN
!
!             C4 PHOTOSYNTHESIS
!
!             LeafAreaNode_brch,CanopyLeafAreaByLayer_pft=leaf area
!             iPlantPhotosynthesisType=photosynthesis type:3=C3,4=C4 from PFT file
!             Vmax4PEPCarboxy_pft=PEP carboxylation rate unlimited by CO2
!
            IF(iPlantPhotosynthesisType(NZ).EQ.ic4_photo.AND.Vmax4PEPCarboxy_pft(K,NB,NZ).GT.0.0_r8)THEN
!
              CALL ComputeGPP_C4(K,NB,NZ,WFNG,Stomata_Activity,CH2O3(K),CH2O4(K))
              CO2F=CO2F+CH2O4(K)
              CH2O=CH2O+CH2O3(K)
!              write(192,*)'ComputeGPP_C4',K,NB,NZ,CH2O3(K),CH2O4(K)
              
!
!               C3 PHOTOSYNTHESIS
!
            ELSEIF(iPlantPhotosynthesisType(NZ).NE.ic4_photo &
              .AND.Vmax4RubiscoCarboxy_pft(K,NB,NZ).GT.0.0_r8)THEN
              call ComputeGPP_C3(K,NB,NZ,WFNG,Stomata_Activity,CH2O3(K))
              CO2F=CO2F+CH2O3(K)
              CH2O=CH2O+CH2O3(K)
!              write(191,*)'ComputeGPP_C3',K,NB,NZ,CH2O3(K)
            ENDIF
          ENDIF
        ENDDO D100
!
!         CO2F,CH2O=total CO2 fixation,CH2O production
!
        CO2F=CO2F*0.0432_r8
        CH2O=CH2O*0.0432_r8
!
!         CONVERT UMOL M-2 S-1 TO G C M-2 H-1
!
        D150: DO K=1,MaxNodesPerBranch1
          CH2O3(K)=CH2O3(K)*0.0432_r8
          CH2O4(K)=CH2O4(K)*0.0432_r8
        ENDDO D150
      ELSE
        CO2F=0._r8
        CH2O=0._r8
        IF(iPlantPhotosynthesisType(NZ).EQ.ic4_photo)THEN
          D155: DO K=1,MaxNodesPerBranch1
            CH2O3(K)=0._r8
            CH2O4(K)=0._r8
          ENDDO D155
        ENDIF
      ENDIF
    ELSE
      CO2F=0._r8
      CH2O=0._r8
      !C4
      IF(iPlantPhotosynthesisType(NZ).EQ.ic4_photo)THEN
        D160: DO K=1,MaxNodesPerBranch1
          CH2O3(K)=0._r8
          CH2O4(K)=0._r8
        ENDDO D160
      ENDIF
    ENDIF
  ELSE
    CO2F=0._r8
    CH2O=0._r8
    IF(iPlantPhotosynthesisType(NZ).EQ.ic4_photo)THEN
      D165: DO K=1,MaxNodesPerBranch1
        CH2O3(K)=0._r8
        CH2O4(K)=0._r8
      ENDDO D165
    ENDIF
  ENDIF
  end associate
  end subroutine ComputeGPP
end module PhotoSynsMod
