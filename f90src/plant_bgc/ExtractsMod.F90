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
  character(len=*),private, parameter :: mod_filename = __FILE__

  public :: extracts
  contains

  SUBROUTINE extracts(I,J)
!     execution begins here
  implicit none

  integer, intent(in) :: I, J

  integer :: NZ

  call TotalLitterfall()

  DO NZ=1,plt_site%NP
    IF(plt_pheno%IFLGC(NZ).EQ.1)THEN

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
  integer :: NE
  associate(                             &
   NP0      => plt_site%NP0        , &
   WGLFT    => plt_biom%WGLFT      , &
   WTSTGT   => plt_biom%WTSTGT     , &
   WTSTGE   => plt_biom%WTSTGE     , &
   HESNC    => plt_bgcr%HESNC      , &
   ZESNC    => plt_bgcr%ZESNC      , &
   ESNT     => plt_bgcr%ESNT       , &
   ESNC     => plt_bgcr%ESNC       , &
   NI       => plt_morph%NI        , &
   ARSTT    => plt_morph%ARSTT     , &
   ARLFT    =>  plt_morph%ARLFT    , &
   ARSTC    =>  plt_morph%ARSTC    , &
   ARLFC    =>  plt_morph%ARLFC      &
  )
  DO NZ=1,NP0
!
!   TOTAL LITTERFALL OF ALL PLANT SPECIES
!
!   ZCSNC,ZZSNC,ZPSNC=total C,N,P litterfall
!   HCSNC,HZSNC,HPSNC=hourly PFT C,N,P litterfall from grosub.f
!   WTSTGT=total standing dead C,N,P mass
!   WTSTG=PFT standing dead C,N,P mass
!   ESNC,=cumulative PFT C,N,P litterfall from grosub.f
!   ESNT,=cumulative total C,N,P litterfall
!
    DO NE=1,npelms
      ZESNC(NE)=ZESNC(NE)+HESNC(NE,NZ)
    ENDDO
    WTSTGT=WTSTGT+WTSTGE(ielmc,NZ)
    DO  L=0,NI(NZ)
      DO K=1,pltpar%n_pltlitrk
        DO NE=1,npelms
          DO  M=1,pltpar%jsken
            ESNT(M,NE,K,L)=ESNT(M,NE,K,L)+ESNC(M,NE,K,L,NZ)
          enddo
        ENDDO
      ENDDO
    ENDDO
  ENDDO
  ARLFC=0._r8
  ARSTC=0._r8
  DO  L=1,JC1
    ARLFT(L)=0._r8
    WGLFT(L)=0._r8
    ARSTT(L)=0._r8
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
    WGLFV    => plt_biom%WGLFV      , &
    WGLFT    => plt_biom%WGLFT      , &
    ARLFT    =>  plt_morph%ARLFT    , &
    ARSTV    =>  plt_morph%ARSTV    , &
    ARSTT    => plt_morph%ARSTT     , &
    ARLFV    => plt_morph%ARLFV       &
  )
  DO L=1,JC1
    ARLFT(L)=ARLFT(L)+ARLFV(L,NZ)
    WGLFT(L)=WGLFT(L)+WGLFV(L,NZ)
    ARSTT(L)=ARSTT(L)+ARSTV(L,NZ)
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

  integer :: N,L,K,NTG

  associate(                       &
    NU    => plt_site%NU     , &
    AREA3 => plt_site%AREA3  , &
    PP    => plt_site%PP     , &
    RUPP1B=> plt_rbgc%RUPP1B , &
    RUPP2B=> plt_rbgc%RUPP2B , &
    RUNNXP=> plt_rbgc%RUNNXP , &
    trcg_RDFA=> plt_rbgc%trcg_RDFA , &
    trcg_RFLA=> plt_rbgc%trcg_RFLA , &
    RCO2P => plt_rbgc%RCO2P  , &
    RUPCHS=> plt_rbgc%RUPCHS , &
    RUPOXP=> plt_rbgc%RUPOXP , &
    RUPN2S=> plt_rbgc%RUPN2S , &
    RUPN3B=> plt_rbgc%RUPN3B , &
    RUPHGS=> plt_rbgc%RUPHGS , &
    RUPN3S=> plt_rbgc%RUPN3S , &
    RUPNOB=> plt_rbgc%RUPNOB , &
    RCO2S => plt_rbgc%RCO2S  , &
    RUPOXS=> plt_rbgc%RUPOXS , &
    RUPNH4=> plt_rbgc%RUPNH4 , &
    RUPNO3=> plt_rbgc%RUPNO3 , &
    RUPH1P=> plt_rbgc%RUPH1P , &
    RUPH2P=> plt_rbgc%RUPH2P , &
    RUPNHB=> plt_rbgc%RUPNHB , &
    RUPH2B=> plt_rbgc%RUPH2B , &
    RUPH1B=> plt_rbgc%RUPH1B , &
    trcg_TFLA=> plt_rbgc%trcg_TFLA , &
    trcg_TLP=> plt_rbgc%trcg_TLP , &
    ROXYP => plt_rbgc%ROXYP  , &
    RDFOME=> plt_rbgc%RDFOME , &
    RUNNHP=> plt_rbgc%RUNNHP , &
    RUNNOP=> plt_rbgc%RUNNOP , &
    RUPP2P=> plt_rbgc%RUPP2P , &
    RUNNBP=> plt_rbgc%RUNNBP , &
    RUPP1P=> plt_rbgc%RUPP1P , &
    RNO3X => plt_bgcr%RNO3X  , &
    RNH4X => plt_bgcr%RNH4X  , &
    RPO4X => plt_bgcr%RPO4X  , &
    RN3BX => plt_bgcr%RN3BX  , &
    RP14X => plt_bgcr%RP14X  , &
    RNHBX => plt_bgcr%RNHBX  , &
    ROXYX => plt_bgcr%ROXYX  , &
    TDFOMP=> plt_bgcr%TDFOMP , &
    TDFOMN=> plt_bgcr%TDFOMN , &
    TDFOMC=> plt_bgcr%TDFOMC , &
    TUPH1B=> plt_bgcr%TUPH1B , &
    TUPH2B=> plt_bgcr%TUPH2B , &
    TUPNOB=> plt_bgcr%TUPNOB , &
    TUPNHB=> plt_bgcr%TUPNHB , &
    TUPH1P=> plt_bgcr%TUPH1P , &
    TUPNO3=> plt_bgcr%TUPNO3 , &
    TUPH2P=> plt_bgcr%TUPH2P , &
    TUPN3S=> plt_bgcr%TUPN3S , &
    TUPN3B=> plt_bgcr%TUPN3B , &
    TUPNH4=> plt_bgcr%TUPNH4 , &
    TUPHGS=> plt_bgcr%TUPHGS , &
    TUPOXS=> plt_bgcr%TUPOXS , &
    TUPCHS=> plt_bgcr%TUPCHS , &
    TUPN2S=> plt_bgcr%TUPN2S , &
    TUPOXP=> plt_bgcr%TUPOXP , &
    TCO2S => plt_bgcr%TCO2S  , &
    TCO2P => plt_bgcr%TCO2P  , &
    RPOBX => plt_bgcr%RPOBX  , &
    RP1BX => plt_bgcr%RP1BX  , &
    TKS   => plt_ew%TKS      , &
    TUPHT => plt_ew%TUPHT    , &
    TUPWTR=> plt_ew%TUPWTR   , &
    UPWTR => plt_ew%UPWTR    , &
    trcg_rootml  => plt_rbgc%trcg_rootml,&
    trcs_rootml => plt_rbgc%trcs_rootml, &
    RTDNP => plt_morph%RTDNP , &
    RTDNT => plt_morph%RTDNT , &
    MY    => plt_morph%MY    , &
    NI    => plt_morph%NI      &
  )

  DO N=1,MY(NZ)
    DO L=NU,NI(NZ)
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
        RTDNT(L)=RTDNT(L)+RTDNP(N,L,NZ)*PP(NZ)/AREA3(L)
      ENDIF
!
!     TOTAL WATER UPTAKE
!
      TUPWTR(L)=TUPWTR(L)+UPWTR(N,L,NZ)
      TUPHT(L)=TUPHT(L)+UPWTR(N,L,NZ)*cpw*TKS(L)
!
!     ROOT GAS CONTENTS FROM FLUXES IN 'UPTAKE'
!
!     *A,*P=PFT root gaseous, aqueous gas content
!     gas code:CO=CO2,OX=O2,CH=CH4,N2=N2O,NH=NH3,H2=H2
!     R*FLA=root gaseous-atmosphere CO2 exchange
!     R*DFA=root aqueous-gaseous CO2 exchange
!
      DO NTG=idg_beg,idg_end-1
        trcg_rootml(NTG,N,L,NZ)=trcg_rootml(NTG,N,L,NZ)+trcg_RFLA(NTG,N,L,NZ)-trcg_RDFA(NTG,N,L,NZ)
      ENDDO

      trcs_rootml(idg_CO2,N,L,NZ)=trcs_rootml(idg_CO2,N,L,NZ)+trcg_RDFA(idg_CO2,N,L,NZ)+RCO2P(N,L,NZ)
      trcs_rootml(idg_O2,N,L,NZ)=trcs_rootml(idg_O2,N,L,NZ)+trcg_RDFA(idg_O2,N,L,NZ)-RUPOXP(N,L,NZ)
      trcs_rootml(idg_CH4,N,L,NZ)=trcs_rootml(idg_CH4,N,L,NZ)+trcg_RDFA(idg_CH4,N,L,NZ)+RUPCHS(N,L,NZ)
      trcs_rootml(idg_N2O,N,L,NZ)=trcs_rootml(idg_N2O,N,L,NZ)+trcg_RDFA(idg_N2O,N,L,NZ)+RUPN2S(N,L,NZ)
      trcs_rootml(idg_NH3,N,L,NZ)=trcs_rootml(idg_NH3,N,L,NZ)+trcg_RDFA(idg_NH3,N,L,NZ)+RUPN3S(N,L,NZ)+RUPN3B(N,L,NZ)
      trcs_rootml(idg_H2,N,L,NZ)=trcs_rootml(idg_H2,N,L,NZ)+trcg_RDFA(idg_H2,N,L,NZ)+RUPHGS(N,L,NZ)
!
!     TOTAL ROOT GAS CONTENTS
!
!     TL*P=total root gas content
!     *A,*P=PFT root gaseous, aqueous gas content
!
    DO NTG=idg_beg,idg_end-1
      trcg_TLP(NTG,L)=trcg_TLP(NTG,L)+trcs_rootml(NTG,N,L,NZ)+trcg_rootml(NTG,N,L,NZ)
    ENDDO
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
      DO NTG=idg_beg,idg_end-1
        trcg_TFLA(NTG,L)=trcg_TFLA(NTG,L)+trcg_RFLA(NTG,N,L,NZ)
      ENDDO

      TCO2P(L)=TCO2P(L)-RCO2P(N,L,NZ)
      TUPOXP(L)=TUPOXP(L)+RUPOXP(N,L,NZ)
      TCO2S(L)=TCO2S(L)+RCO2S(N,L,NZ)
      TUPOXS(L)=TUPOXS(L)+RUPOXS(N,L,NZ)
      TUPCHS(L)=TUPCHS(L)+RUPCHS(N,L,NZ)
      TUPN2S(L)=TUPN2S(L)+RUPN2S(N,L,NZ)
      TUPN3S(L)=TUPN3S(L)+RUPN3S(N,L,NZ)
      TUPN3B(L)=TUPN3B(L)+RUPN3B(N,L,NZ)
      TUPHGS(L)=TUPHGS(L)+RUPHGS(N,L,NZ)
      TUPNH4(L)=TUPNH4(L)+RUPNH4(N,L,NZ)
      TUPNO3(L)=TUPNO3(L)+RUPNO3(N,L,NZ)
      TUPH2P(L)=TUPH2P(L)+RUPH2P(N,L,NZ)
      TUPH1P(L)=TUPH1P(L)+RUPH1P(N,L,NZ)
      TUPNHB(L)=TUPNHB(L)+RUPNHB(N,L,NZ)
      TUPNOB(L)=TUPNOB(L)+RUPNOB(N,L,NZ)
      TUPH2B(L)=TUPH2B(L)+RUPH2B(N,L,NZ)
      TUPH1B(L)=TUPH1B(L)+RUPH1B(N,L,NZ)
!
!     TOTAL ROOT C,N,P EXUDATION
!
!     TDFOMC,TDFOMN,TDFOMP=total nonstructl C,N,P exchange
!     RDFOMC,RDFOMN,RDFOMP=PFT nonstructl C,N,P exchange
!
      DO K=1,jcplx
        TDFOMC(K,L)=TDFOMC(K,L)-RDFOME(ielmc,N,K,L,NZ)
        TDFOMN(K,L)=TDFOMN(K,L)-RDFOME(ielmn,N,K,L,NZ)
        TDFOMP(K,L)=TDFOMP(K,L)-RDFOME(ielmp,N,K,L,NZ)
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
      ROXYX(L)=ROXYX(L)+ROXYP(N,L,NZ)
      RNH4X(L)=RNH4X(L)+RUNNHP(N,L,NZ)
      RNO3X(L)=RNO3X(L)+RUNNOP(N,L,NZ)
      RPO4X(L)=RPO4X(L)+RUPP2P(N,L,NZ)
      RP14X(L)=RP14X(L)+RUPP1P(N,L,NZ)
      RNHBX(L)=RNHBX(L)+RUNNBP(N,L,NZ)
      RN3BX(L)=RN3BX(L)+RUNNXP(N,L,NZ)
      RPOBX(L)=RPOBX(L)+RUPP2B(N,L,NZ)
      RP1BX(L)=RP1BX(L)+RUPP1B(N,L,NZ)
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
  integer :: L, NB,NE,NTG
  real(r8) :: ENGYC

  associate(                       &
    TBALE => plt_site%TBALE  , &
    BALE  => plt_site%BALE   , &
    TNH3C => plt_bgcr%TNH3C  , &
    RNH3C => plt_bgcr%RNH3C  , &
    TCCAN => plt_bgcr%TCCAN  , &
    ZESNC => plt_bgcr%ZESNC  , &
    RFGas_root => plt_bgcr%RFGas_root, &
    RUPNF => plt_bgcr%RUPNF  , &
    CNET  => plt_bgcr%CNET   , &
    CTRAN => plt_ew%CTRAN    , &
    TH2GZ => plt_bgcr%TH2GZ  , &
    RNH3B => plt_rbgc%RNH3B  , &
    HEUPTK=> plt_rbgc%HEUPTK , &
    TUPNF => plt_rbgc%TUPNF  , &
    TRFGas_root => plt_rbgc%TRFGas_root  , &
    EP    => plt_ew%EP       , &
    FLWC  => plt_ew%FLWC     , &
    EVAPC => plt_ew%EVAPC    , &
    VOLWC => plt_ew%VOLWC    , &
    VOLWP => plt_ew%VOLWP    , &
    TGH   => plt_ew%TGH      , &
    SFLXC => plt_ew%SFLXC    , &
    EFLXC => plt_ew%EFLXC    , &
    TVOLWP=> plt_ew%TVOLWP   , &
    TKC   => plt_ew%TKC      , &
    TKS   => plt_ew%TKS      , &
    ENGYX => plt_ew%ENGYX    , &
    TSH   => plt_ew%TSH      , &
    TEVAPC=> plt_ew%TEVAPC   , &
    TENGYC=> plt_ew%TENGYC   , &
    TEVAPP=> plt_ew%TEVAPP   , &
    THFLXC=> plt_ew%THFLXC   , &
    THRMC => plt_ew%THRMC    , &
    TKA   => plt_ew%TKA      , &
    HFLXC => plt_ew%HFLXC    , &
    TLE   => plt_ew%TLE      , &
    TVOLWC=> plt_ew%TVOLWC   , &
    NU    => plt_site%NU     , &
    ARSTC => plt_morph%ARSTC , &
    ARLFC => plt_morph%ARLFC , &
    NI    => plt_morph%NI    , &
    NBR   => plt_morph%NBR   , &
    ARSTP => plt_morph%ARSTP , &
    ARLFP => plt_morph%ARLFP , &
    RAD1  => plt_rad%RAD1    , &
    THRM1 => plt_rad%THRM1   , &
    TRN   => plt_rad%TRN       &
  )
  DO L=NU,NI(NZ)
    TUPNF(L)=TUPNF(L)+RUPNF(L,NZ)
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
!     TRFGas_root=total loss of root CO2, O2, CH4, N2O, NH3, H2
!     RFGas_root=PFT loss of root CO2, O2, CH4, N2O, NH3, H2
!
  TRN=TRN+RAD1(NZ)
  TLE=TLE+EFLXC(NZ)
  TSH=TSH+SFLXC(NZ)
  TGH=TGH+HFLXC(NZ)
  TCCAN=TCCAN+CNET(NZ)
  CTRAN(NZ)=CTRAN(NZ)+EP(NZ)+EVAPC(NZ)
  TVOLWP=TVOLWP+VOLWP(NZ)
  TVOLWC=TVOLWC+VOLWC(NZ)
  TEVAPP=TEVAPP+EP(NZ)+EVAPC(NZ)
  TEVAPC=TEVAPC+EVAPC(NZ)
  ENGYC=cpw*(VOLWC(NZ)+FLWC(NZ)+EVAPC(NZ))*TKC(NZ)
  TENGYC=TENGYC+ENGYC
  THFLXC=THFLXC+ENGYC-ENGYX(NZ)-(FLWC(NZ)*cpw*TKA)
  ENGYX(NZ)=ENGYC
  THRMC=THRMC+THRM1(NZ)
  ARLFC=ARLFC+ARLFP(NZ)
  ARSTC=ARSTC+ARSTP(NZ)
  DO NE=1,npelms
    ZESNC(NE)=ZESNC(NE)-HEUPTK(NE,NZ)
    TBALE(NE)=TBALE(NE)+BALE(NE,NZ)
  ENDDO

  DO NTG=idg_beg,idg_end-1
    TRFGas_root(NTG)=TRFGas_root(NTG)+RFGas_root(NTG,NZ)
  ENDDO
!
!     TOTAL CANOPY NH3 EXCHANGE AND EXUDATION
!
!     RNH3B,RNH3C=PFT NH3 flux between atmosphere and branch,canopy
!     TNH3C=total NH3 flux between atmosphere and canopy
!
  RNH3C(NZ)=0._r8
  DO NB=1,NBR(NZ)
    RNH3C(NZ)=RNH3C(NZ)+RNH3B(NB,NZ)
    TNH3C(NZ)=TNH3C(NZ)+RNH3B(NB,NZ)
  ENDDO
  end associate
  end subroutine CanopyFluxesandFixation

  end module ExtractsMod
