module ExtractsMod
!!
!Description:
!     THIS SUBROUTINE AGGREGATES ALL SOIL-PLANT C,N,P EXCHANGES
!     FROM 'UPTAKE' AMD 'GROSUB' AND SENDS RESULTS TO 'REDIST'
!
  use data_kind_mod, only : r8 => SHR_KIND_R8
  use EcosimConst
  use GrosubPars
  use PlantAPIData
  implicit none

  private

  public :: extracts
  contains

  SUBROUTINE extracts(I,J)
!     execution begins here
  implicit none

  integer, intent(in) :: I, J

  integer :: NZ

  call TotalLitterfall()

  DO NZ=1,NPs1
    IF(plt_pheno%IFLGCs1(NZ).EQ.1)THEN

      call TotalLeafArea(NZ)

      call TotalGasandSoluteUptake(NZ)

      call CanopyFluxesandFixation(NZ)

    ENDIF
  ENDDO
  RETURN
  END subroutine extracts
!------------------------------------------------------------------------------------------

  subroutine TotalLitterfall()

  implicit none
  integer :: NZ,L,K,M

  associate(                             &
   WGLFTs1    => plt_biom%WGLFTs1      , &
   WTSTGTs1   => plt_biom%WTSTGTs1     , &
   WTSTGs1    => plt_biom%WTSTGs1      , &
   NIs1       => plt_morph%NIs1        , &
   ARSTTs1    => plt_morph%ARSTTs1     , &
   ARLFTs1    =>  plt_morph%ARLFTs1    , &
   ARSTCs1    =>  plt_morph%ARSTCs1    , &
   ARLFCs1    =>  plt_morph%ARLFCs1      &
  )
  DO NZ=1,NP0s1
!
!   TOTAL LITTERFALL OF ALL PLANT SPECIES
!
!   ZCSNC,ZZSNC,ZPSNC=total C,N,P litterfall
!   HCSNC,HZSNC,HPSNC=hourly PFT C,N,P litterfall from grosub.f
!   WTSTGT=total standing dead C,N,P mass
!   WTSTG=PFT standing dead C,N,P mass
!   CSNC,ZSNC,PSNC=cumulative PFT C,N,P litterfall from grosub.f
!   CSNT,ZSNT,PSNT=cumulative total C,N,P litterfall
!
    ZCSNCs1=ZCSNCs1+HCSNCs1(NZ)
    ZZSNCs1=ZZSNCs1+HZSNCs1(NZ)
    ZPSNCs1=ZPSNCs1+HPSNCs1(NZ)
    WTSTGTs1=WTSTGTs1+WTSTGs1(NZ)
    DO  L=0,NIs1(NZ)
      DO K=0,1
        DO  M=1,jcplx11
          CSNTs1(M,K,L)=CSNTs1(M,K,L)+CSNCs1(M,K,L,NZ)
          ZSNTs1(M,K,L)=ZSNTs1(M,K,L)+ZSNCs1(M,K,L,NZ)
          PSNTs1(M,K,L)=PSNTs1(M,K,L)+PSNCs1(M,K,L,NZ)
        enddo
      ENDDO
    ENDDO
  ENDDO
  ARLFCs1=0._r8
  ARSTCs1=0._r8
  DO  L=1,JC1
    ARLFTs1(L)=0._r8
    WGLFTs1(L)=0._r8
    ARSTTs1(L)=0._r8
  ENDDO
  end associate
  end subroutine TotalLitterfall
!------------------------------------------------------------------------------------------

  subroutine TotalLeafArea(NZ)
!
!     TOTAL LEAF AREA OF ALL PLANT SPECIES
!
!     ARLFT,ARSTT=total leaf,stalk area of combined canopy layer
!     ARLFV,ARSTV=PFT leaf,stalk area in canopy layer
!     WGLFT=total leaf C of combined canopy layer
!     WGLFV=PFT leaf C in canopy layer
!
  implicit none
  integer, intent(in) :: NZ
  integer :: L
  associate(                              &
    WGLFTs1    => plt_biom%WGLFTs1      , &
    ARLFTs1    =>  plt_morph%ARLFTs1    , &
    ARSTVs1    =>  plt_morph%ARSTVs1    , &
    ARSTTs1    => plt_morph%ARSTTs1     , &
    ARLFVs1    => plt_morph%ARLFVs1       &
  )
  DO L=1,JC1
    ARLFTs1(L)=ARLFTs1(L)+ARLFVs1(L,NZ)
    WGLFTs1(L)=WGLFTs1(L)+WGLFVs1(L,NZ)
    ARSTTs1(L)=ARSTTs1(L)+ARSTVs1(L,NZ)
  ENDDO
  end associate
  end subroutine TotalLeafArea
!------------------------------------------------------------------------------------------

  subroutine TotalGasandSoluteUptake(NZ)
!
!     TOTAL GAS AND SOLUTE UPTAKE BY ALL PLANT SPECIES
!

  implicit none
  integer, intent(in) :: NZ

  integer :: N,L,K

  associate(                       &
    TKSs1   => plt_ew%TKSs1      , &
    MYs1    => plt_morph%MYs1    , &
    NIs1    => plt_morph%NIs1      &
  )

  DO N=1,MYs1(NZ)
    DO L=NUs1,NIs1(NZ)
!
!     TOTAL ROOT DENSITY
!
!     RTDNT=total root length density
!     RTDNP=PFT root length density per plant
!     UPWTR=total water uptake
!     UPWTR=PFT root water uptake
!     TUPHT=total convective heat in root water uptake
!     TKS=soil temperature
!     PP=PFT population
!
      IF(N.EQ.1)THEN
        RTDNTs1(L)=RTDNTs1(L)+RTDNPs1(N,L,NZ)*PPs1(NZ)/AREA3s1(L)
      ENDIF
!
!     TOTAL WATER UPTAKE
!
      TUPWTRs1(L)=TUPWTRs1(L)+UPWTRs1(N,L,NZ)
      TUPHTs1(L)=TUPHTs1(L)+UPWTRs1(N,L,NZ)*cpw*TKSs1(L)
!
!     ROOT GAS CONTENTS FROM FLUXES IN 'UPTAKE'
!
!     *A,*P=PFT root gaseous, aqueous gas content
!     gas code:CO=CO2,OX=O2,CH=CH4,N2=N2O,NH=NH3,H2=H2
!     R*FLA=root gaseous-atmosphere CO2 exchange
!     R*DFA=root aqueous-gaseous CO2 exchange
!
      CO2As1(N,L,NZ)=CO2As1(N,L,NZ)+RCOFLAs1(N,L,NZ)-RCODFAs1(N,L,NZ)
      OXYAs1(N,L,NZ)=OXYAs1(N,L,NZ)+ROXFLAs1(N,L,NZ)-ROXDFAs1(N,L,NZ)
      CH4As1(N,L,NZ)=CH4As1(N,L,NZ)+RCHFLAs1(N,L,NZ)-RCHDFAs1(N,L,NZ)
      Z2OAs1(N,L,NZ)=Z2OAs1(N,L,NZ)+RN2FLAs1(N,L,NZ)-RN2DFAs1(N,L,NZ)
      ZH3As1(N,L,NZ)=ZH3As1(N,L,NZ)+RNHFLAs1(N,L,NZ)-RNHDFAs1(N,L,NZ)
      H2GAs1(N,L,NZ)=H2GAs1(N,L,NZ)+RHGFLAs1(N,L,NZ)-RHGDFAs1(N,L,NZ)
      CO2Ps1(N,L,NZ)=CO2Ps1(N,L,NZ)+RCODFAs1(N,L,NZ)+RCO2Ps1(N,L,NZ)
      OXYPs1(N,L,NZ)=OXYPs1(N,L,NZ)+ROXDFAs1(N,L,NZ)-RUPOXPs1(N,L,NZ)
      CH4Ps1(N,L,NZ)=CH4Ps1(N,L,NZ)+RCHDFAs1(N,L,NZ)+RUPCHSs1(N,L,NZ)
      Z2OPs1(N,L,NZ)=Z2OPs1(N,L,NZ)+RN2DFAs1(N,L,NZ)+RUPN2Ss1(N,L,NZ)
      ZH3Ps1(N,L,NZ)=ZH3Ps1(N,L,NZ)+RNHDFAs1(N,L,NZ)+RUPN3Ss1(N,L,NZ)+RUPN3Bs1(N,L,NZ)
      H2GPs1(N,L,NZ)=H2GPs1(N,L,NZ)+RHGDFAs1(N,L,NZ)+RUPHGSs1(N,L,NZ)
!
!     TOTAL ROOT GAS CONTENTS
!
!     TL*P=total root gas content
!     *A,*P=PFT root gaseous, aqueous gas content
!
      TLCO2Ps1(L)=TLCO2Ps1(L)+CO2Ps1(N,L,NZ)+CO2As1(N,L,NZ)
      TLOXYPs1(L)=TLOXYPs1(L)+OXYPs1(N,L,NZ)+OXYAs1(N,L,NZ)
      TLCH4Ps1(L)=TLCH4Ps1(L)+CH4Ps1(N,L,NZ)+CH4As1(N,L,NZ)
      TLN2OPs1(L)=TLN2OPs1(L)+Z2OPs1(N,L,NZ)+Z2OAs1(N,L,NZ)
      TLNH3Ps1(L)=TLNH3Ps1(L)+ZH3Ps1(N,L,NZ)+ZH3As1(N,L,NZ)
      TLH2GPs1(L)=TLH2GPs1(L)+H2GPs1(N,L,NZ)+H2GAs1(N,L,NZ)
!
!     TOTAL ROOT BOUNDARY GAS FLUXES
!
!     T*FLA=total root gaseous-atmosphere CO2 exchange
!     R*FLA=PFT root gaseous-atmosphere CO2 exchange
!     gas code:CO=CO2,OX=O2,CH=CH4,N2=N2O,NH=NH3,H2=H2
!     TUP*S,TUP*B=total root-soil gas, solute exchange in non-band,band
!     RUP*S,RUP*B*=PFT root-soil gas, solute exchange in non-band,band
!     gas code:CO=CO2,OX=O2,CH=CH4,N2=N2O,NH=NH3,H2=H2
!     solute code:NH4=NH4,NO3=NO3,H2P=H2PO4,H1P=H1PO4 in non-band
!                :NHB=NH4,NOB=NO3,H2B=H2PO4,H1B=H1PO4 in band
!
      TCOFLAs1(L)=TCOFLAs1(L)+RCOFLAs1(N,L,NZ)
      TOXFLAs1(L)=TOXFLAs1(L)+ROXFLAs1(N,L,NZ)
      TCHFLAs1(L)=TCHFLAs1(L)+RCHFLAs1(N,L,NZ)
      TN2FLAs1(L)=TN2FLAs1(L)+RN2FLAs1(N,L,NZ)
      TNHFLAs1(L)=TNHFLAs1(L)+RNHFLAs1(N,L,NZ)
      THGFLAs1(L)=THGFLAs1(L)+RHGFLAs1(N,L,NZ)
      TCO2Ps1(L)=TCO2Ps1(L)-RCO2Ps1(N,L,NZ)
      TUPOXPs1(L)=TUPOXPs1(L)+RUPOXPs1(N,L,NZ)
      TCO2Ss1(L)=TCO2Ss1(L)+RCO2Ss1(N,L,NZ)
      TUPOXSs1(L)=TUPOXSs1(L)+RUPOXSs1(N,L,NZ)
      TUPCHSs1(L)=TUPCHSs1(L)+RUPCHSs1(N,L,NZ)
      TUPN2Ss1(L)=TUPN2Ss1(L)+RUPN2Ss1(N,L,NZ)
      TUPN3Ss1(L)=TUPN3Ss1(L)+RUPN3Ss1(N,L,NZ)
      TUPN3Bs1(L)=TUPN3Bs1(L)+RUPN3Bs1(N,L,NZ)
      TUPHGSs1(L)=TUPHGSs1(L)+RUPHGSs1(N,L,NZ)
      TUPNH4s1(L)=TUPNH4s1(L)+RUPNH4s1(N,L,NZ)
      TUPNO3s1(L)=TUPNO3s1(L)+RUPNO3s1(N,L,NZ)
      TUPH2Ps1(L)=TUPH2Ps1(L)+RUPH2Ps1(N,L,NZ)
      TUPH1Ps1(L)=TUPH1Ps1(L)+RUPH1Ps1(N,L,NZ)
      TUPNHBs1(L)=TUPNHBs1(L)+RUPNHBs1(N,L,NZ)
      TUPNOBs1(L)=TUPNOBs1(L)+RUPNOBs1(N,L,NZ)
      TUPH2Bs1(L)=TUPH2Bs1(L)+RUPH2Bs1(N,L,NZ)
      TUPH1Bs1(L)=TUPH1Bs1(L)+RUPH1Bs1(N,L,NZ)
!
!     TOTAL ROOT C,N,P EXUDATION
!
!     TDFOMC,TDFOMN,TDFOMP=total nonstructl C,N,P exchange
!     RDFOMC,RDFOMN,RDFOMP=PFT nonstructl C,N,P exchange
!
      DO K=0,jcplx11
        TDFOMCs1(K,L)=TDFOMCs1(K,L)-RDFOMCs1(N,K,L,NZ)
        TDFOMNs1(K,L)=TDFOMNs1(K,L)-RDFOMNs1(N,K,L,NZ)
        TDFOMPs1(K,L)=TDFOMPs1(K,L)-RDFOMPs1(N,K,L,NZ)
      ENDDO
!
!     TOTAL ROOT O2, NH4, NO3, PO4 UPTAKE CONTRIBUTES TO
!     TOTAL ROOT + MICROBIAL UPTAKE USED TO CALCULATE
!     COMPETITION CONSTRAINTS
!
!     ROXYX=O2 demand by all microbial,root,myco populations
!     RNH4X=NH4 demand in non-band by all microbial,root,myco populations
!     RNO3X=NO3 demand in non-band by all microbial,root,myco populations
!     RPO4X=H2PO4 demand in non-band by all microbial,root,myco populations
!     RP14X=HPO4 demand in non-band by all microbial,root,myco populations
!     RNHBX=NH4 demand in band by all microbial,root,myco populations
!     RN3BX=NO3 demand in band by all microbial,root,myco populations
!     RPOBX=H2PO4 demand in band by all microbial,root,myco populations
!     RP1BX=HPO4 demand in band by all microbial,root,myco populations
!     ROXYP=O2 demand by each root,myco population
!     RUNNHP=NH4 demand in non-band by each root population
!     RUNNOP=NO3 demand in non-band by each root population
!     RUPP2P=H2PO4 demand in non-band by each root population
!     RUPP1P=HPO4 demand in non-band by each root population
!     RUNNBP=NH4 demand in band by each root population
!     RUNNXB=NO3 demand in band by each root population
!     RUPP2B=H2PO4 demand in band by each root population
!     RUPP1B=HPO4 demand in band by each root population
!
      ROXYXs1(L)=ROXYXs1(L)+ROXYPs1(N,L,NZ)
      RNH4Xs1(L)=RNH4Xs1(L)+RUNNHPs1(N,L,NZ)
      RNO3Xs1(L)=RNO3Xs1(L)+RUNNOPs1(N,L,NZ)
      RPO4Xs1(L)=RPO4Xs1(L)+RUPP2Ps1(N,L,NZ)
      RP14Xs1(L)=RP14Xs1(L)+RUPP1Ps1(N,L,NZ)
      RNHBXs1(L)=RNHBXs1(L)+RUNNBPs1(N,L,NZ)
      RN3BXs1(L)=RN3BXs1(L)+RUNNXPs1(N,L,NZ)
      RPOBXs1(L)=RPOBXs1(L)+RUPP2Bs1(N,L,NZ)
      RP1BXs1(L)=RP1BXs1(L)+RUPP1Bs1(N,L,NZ)
    ENDDO
  ENDDO
  end associate
  end subroutine TotalGasandSoluteUptake
!------------------------------------------------------------------------------------------

  subroutine CanopyFluxesandFixation(NZ)
!
!     TOTAL ROOT N2 FIXATION BY ALL PLANT SPECIES
!
!     TUPNF=total root N2 fixation
!     RUPNF=PFT root N2 fixation
!

  implicit none
  integer, intent(in) :: NZ
  integer :: L, NB
  real(r8) :: ENGYC
  associate(                       &
    TVOLWPs1=> plt_ew%TVOLWPs1   , &
    TKCs1   => plt_ew%TKCs1      , &
    TKSs1   => plt_ew%TKSs1      , &
    ENGYXs1 => plt_ew%ENGYXs1    , &
    TSHs1   => plt_ew%TSHs1      , &
    TEVAPCs1=> plt_ew%TEVAPCs1   , &
    TENGYCs1=> plt_ew%TENGYCs1   , &
    TEVAPPs1=> plt_ew%TEVAPPs1   , &
    THFLXCs1=> plt_ew%THFLXCs1   , &
    THRMCs1 => plt_ew%THRMCs1    , &
    TKAs1   => plt_ew%TKAs1      , &
    TLEs1   => plt_ew%TLEs1      , &
    TVOLWCs1=> plt_ew%TVOLWCs1   , &
    ARSTCs1 => plt_morph%ARSTCs1 , &
    ARLFCs1 => plt_morph%ARLFCs1 , &
    NIs1    => plt_morph%NIs1    , &
    NBRs1   => plt_morph%NBRs1   , &
    ARSTPs1 => plt_morph%ARSTPs1 , &
    ARLFPs1 => plt_morph%ARLFPs1 , &
    RAD1s1  => plt_rad%RAD1s1    , &
    THRM1s1 => plt_rad%THRM1s1   , &
    TRNs1   => plt_rad%TRNs1       &
  )
  DO L=NUs1,NIs1(NZ)
    TUPNFs1(L)=TUPNFs1(L)+RUPNFs1(L,NZ)
  ENDDO
!
!     TOTAL ENERGY, WATER, CO2 FLUXES
!
!     TRN=total net SW+LW absorbed by canopy
!     RAD1=PFT net SW+LW absorbed by canopy
!     TLE=total canopy latent heat flux
!     EFLXC=PFT canopy latent heat flux
!     TSH=total canopy sensible heat flux
!     SFLXC=PFT canopy sensible heat flux
!     TGH=total canopy storage heat flux
!     HFLXC=PFT canopy storage heat flux
!     TCCAN=total net CO2 fixation
!     CNET=PFT net CO2 fixation
!     TVOLWP,TVOLWC=total water volume in canopy,on canopy surfaces
!     VOLWP,VOLWC=PFT water volume in canopy,on canopy surfaces
!     TEVAPP,TEVAPC=total water flux to,from canopy,canopy surfaces
!     EVAPC,EP=water flux to,from canopy surfaces, inside canopy
!     TENGYC=total canopy water heat content
!     ENGYC=PFT canopy water heat content
!     ARLFC,ARSTC=total leaf,stalk area
!     ARLFP,ARSTP=PFT leaf,stalk area
!     ZCSNC,ZZSNC,ZPSNC=total net root-soil C,N,P exchange
!     HCUPTK,HZUPTK,HPUPTK=PFT net root-soil C,N,P exchange
!     TBALC,TBALN,TBALP=total C,N,P balance
!     BALC,BALN,BALP=PFT C,N,P balance
!     TCO2Z,TOXYZ,TCH4Z,TN2OZ,TNH3Z,TH2GZ=total loss of root CO2, O2, CH4, N2O, NH3, H2
!     RCO2Z,ROXYZ,RCH4Z,RN2OZ,RNH3Z,RH2GZ=PFT loss of root CO2, O2, CH4, N2O, NH3, H2
!
  TRNs1=TRNs1+RAD1s1(NZ)
  TLEs1=TLEs1+EFLXCs1(NZ)
  TSHs1=TSHs1+SFLXCs1(NZ)
  TGHs1=TGHs1+HFLXCs1(NZ)
  TCCANs1=TCCANs1+CNETs1(NZ)
  CTRANs1(NZ)=CTRANs1(NZ)+EPs1(NZ)+EVAPCs1(NZ)
  TVOLWPs1=TVOLWPs1+VOLWPs1(NZ)
  TVOLWCs1=TVOLWCs1+VOLWCs1(NZ)
  TEVAPPs1=TEVAPPs1+EPs1(NZ)+EVAPCs1(NZ)
  TEVAPCs1=TEVAPCs1+EVAPCs1(NZ)
  ENGYC=cpw*(VOLWCs1(NZ)+FLWCs1(NZ)+EVAPCs1(NZ))*TKCs1(NZ)
  TENGYCs1=TENGYCs1+ENGYC
  THFLXCs1=THFLXCs1+ENGYC-ENGYXs1(NZ)-(FLWCs1(NZ)*cpw*TKAs1)
  ENGYXs1(NZ)=ENGYC
  THRMCs1=THRMCs1+THRM1s1(NZ)
  ARLFCs1=ARLFCs1+ARLFPs1(NZ)
  ARSTCs1=ARSTCs1+ARSTPs1(NZ)
  ZCSNCs1=ZCSNCs1-HCUPTKs1(NZ)
  ZZSNCs1=ZZSNCs1-HZUPTKs1(NZ)
  ZPSNCs1=ZPSNCs1-HPUPTKs1(NZ)
  TBALCs1=TBALCs1+BALCs1(NZ)
  TBALNs1=TBALNs1+BALNs1(NZ)
  TBALPs1=TBALPs1+BALPs1(NZ)
  TCO2Zs1=TCO2Zs1+RCO2Zs1(NZ)
  TOXYZs1=TOXYZs1+ROXYZs1(NZ)
  TCH4Zs1=TCH4Zs1+RCH4Zs1(NZ)
  TN2OZs1=TN2OZs1+RN2OZs1(NZ)
  TNH3Zs1=TNH3Zs1+RNH3Zs1(NZ)
  TH2GZs1=TH2GZs1+RH2GZs1(NZ)
!
!     TOTAL CANOPY NH3 EXCHANGE AND EXUDATION
!
!     RNH3B,RNH3C=PFT NH3 flux between atmosphere and branch,canopy
!     TNH3C=total NH3 flux between atmosphere and canopy
!
  RNH3Cs1(NZ)=0._r8
  DO NB=1,NBRs1(NZ)
    RNH3Cs1(NZ)=RNH3Cs1(NZ)+RNH3Bs1(NB,NZ)
    TNH3Cs1(NZ)=TNH3Cs1(NZ)+RNH3Bs1(NB,NZ)
  ENDDO
  end associate
  end subroutine CanopyFluxesandFixation

  end module ExtractsMod
