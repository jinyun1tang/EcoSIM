module ExtractMod
!!
!Description:
!     THIS SUBROUTINE AGGREGATES ALL SOIL-PLANT C,N,P EXCHANGES
!     FROM 'UPTAKE' AMD 'GROSUB' AND SENDS RESULTS TO 'REDIST'
!
  use data_kind_mod, only : r8 => SHR_KIND_R8
  use EcosimConst
  use GridDataType
  use FlagDataType
  use SoilHeatDatatype
  use SoilWaterDataType
  use PlantDataStateType
  use PlantDataRateType
  use PlantDataCharType
  use ClimForcDataType
  use PhenologyDataType
  use CanopyDataType
  USE RootDataType
  use EcoSimSumDataType
  use SOMDataType
  use SoilBGCDataType
  use EcosysBGCFluxType
  implicit none

  private


  real(r8) :: ENGYC

  integer :: K,L,M,NX,NY,NZ,N,NB

  public :: extract
  contains

  SUBROUTINE extract(I,J,NHW,NHE,NVN,NVS)
!     execution begins here
  implicit none

  integer, intent(in) :: I, J
  integer, intent(in) :: NHW,NHE,NVN,NVS

  DO 9995 NX=NHW,NHE
    DO 9990 NY=NVN,NVS

      call TotalLitterfall(NY,NX)

      DO 9975 NZ=1,NP(NY,NX)
        IF(IFLGC(NZ,NY,NX).EQ.1)THEN

          call TotalLeafArea(NZ,NY,NX)

          call TotalGasandSoluteUptake(NZ,NY,NX)

          call CanopyFluxesandFixation(NZ,NY,NX)

        ENDIF
9975  CONTINUE
9990  CONTINUE
9995  CONTINUE
  RETURN
  END subroutine extract
!------------------------------------------------------------------------------------------

  subroutine TotalLitterfall(NY,NX)

  implicit none
  integer, intent(in) :: NY,NX

  DO 9985 NZ=1,NP0(NY,NX)
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
    ZCSNC(NY,NX)=ZCSNC(NY,NX)+HCSNC(NZ,NY,NX)
    ZZSNC(NY,NX)=ZZSNC(NY,NX)+HZSNC(NZ,NY,NX)
    ZPSNC(NY,NX)=ZPSNC(NY,NX)+HPSNC(NZ,NY,NX)
    WTSTGT(NY,NX)=WTSTGT(NY,NX)+WTSTG(NZ,NY,NX)
    DO 90 L=0,NI(NZ,NY,NX)
      DO 95 K=0,1
        DO  M=1,4
          CSNT(M,K,L,NY,NX)=CSNT(M,K,L,NY,NX)+CSNC(M,K,L,NZ,NY,NX)
          ZSNT(M,K,L,NY,NX)=ZSNT(M,K,L,NY,NX)+ZSNC(M,K,L,NZ,NY,NX)
          PSNT(M,K,L,NY,NX)=PSNT(M,K,L,NY,NX)+PSNC(M,K,L,NZ,NY,NX)
        enddo
95    CONTINUE
90  CONTINUE
9985  CONTINUE
  ARLFC(NY,NX)=0._r8
  ARSTC(NY,NX)=0._r8
  DO 915 L=1,JC
    ARLFT(L,NY,NX)=0._r8
    WGLFT(L,NY,NX)=0._r8
    ARSTT(L,NY,NX)=0._r8
915   CONTINUE
  end subroutine TotalLitterfall
!------------------------------------------------------------------------------------------

  subroutine TotalLeafArea(NZ,NY,NX)
!
!     TOTAL LEAF AREA OF ALL PLANT SPECIES
!
!     ARLFT,ARSTT=total leaf,stalk area of combined canopy layer
!     ARLFV,ARSTV=PFT leaf,stalk area in canopy layer
!     WGLFT=total leaf C of combined canopy layer
!     WGLFV=PFT leaf C in canopy layer
!
  implicit none
  integer, intent(in) :: NZ,NY,NX

  DO 910 L=1,JC
    ARLFT(L,NY,NX)=ARLFT(L,NY,NX)+ARLFV(L,NZ,NY,NX)
    WGLFT(L,NY,NX)=WGLFT(L,NY,NX)+WGLFV(L,NZ,NY,NX)
    ARSTT(L,NY,NX)=ARSTT(L,NY,NX)+ARSTV(L,NZ,NY,NX)
910   CONTINUE
  end subroutine TotalLeafArea
!------------------------------------------------------------------------------------------

  subroutine TotalGasandSoluteUptake(NZ,NY,NX)
!
!     TOTAL GAS AND SOLUTE UPTAKE BY ALL PLANT SPECIES
!

  implicit none
  integer, intent(in) :: NZ, NY, NX

  DO 100 N=1,MY(NZ,NY,NX)
    DO L=NU(NY,NX),NI(NZ,NY,NX)
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
        RTDNT(L,NY,NX)=RTDNT(L,NY,NX)+RTDNP(N,L,NZ,NY,NX)*PP(NZ,NY,NX)/AREA(3,L,NY,NX)
      ENDIF
!
!     TOTAL WATER UPTAKE
!
      TUPWTR(L,NY,NX)=TUPWTR(L,NY,NX)+UPWTR(N,L,NZ,NY,NX)
      TUPHT(L,NY,NX)=TUPHT(L,NY,NX)+UPWTR(N,L,NZ,NY,NX)*cpw*TKS(L,NY,NX)
!
!     ROOT GAS CONTENTS FROM FLUXES IN 'UPTAKE'
!
!     *A,*P=PFT root gaseous, aqueous gas content
!     gas code:CO=CO2,OX=O2,CH=CH4,N2=N2O,NH=NH3,H2=H2
!     R*FLA=root gaseous-atmosphere CO2 exchange
!     R*DFA=root aqueous-gaseous CO2 exchange
!
      CO2A(N,L,NZ,NY,NX)=CO2A(N,L,NZ,NY,NX)+RCOFLA(N,L,NZ,NY,NX)-RCODFA(N,L,NZ,NY,NX)
      OXYA(N,L,NZ,NY,NX)=OXYA(N,L,NZ,NY,NX)+ROXFLA(N,L,NZ,NY,NX)-ROXDFA(N,L,NZ,NY,NX)
      CH4A(N,L,NZ,NY,NX)=CH4A(N,L,NZ,NY,NX)+RCHFLA(N,L,NZ,NY,NX)-RCHDFA(N,L,NZ,NY,NX)
      Z2OA(N,L,NZ,NY,NX)=Z2OA(N,L,NZ,NY,NX)+RN2FLA(N,L,NZ,NY,NX)-RN2DFA(N,L,NZ,NY,NX)
      ZH3A(N,L,NZ,NY,NX)=ZH3A(N,L,NZ,NY,NX)+RNHFLA(N,L,NZ,NY,NX)-RNHDFA(N,L,NZ,NY,NX)
      H2GA(N,L,NZ,NY,NX)=H2GA(N,L,NZ,NY,NX)+RHGFLA(N,L,NZ,NY,NX)-RHGDFA(N,L,NZ,NY,NX)
      CO2P(N,L,NZ,NY,NX)=CO2P(N,L,NZ,NY,NX)+RCODFA(N,L,NZ,NY,NX)+RCO2P(N,L,NZ,NY,NX)
      OXYP(N,L,NZ,NY,NX)=OXYP(N,L,NZ,NY,NX)+ROXDFA(N,L,NZ,NY,NX)-RUPOXP(N,L,NZ,NY,NX)
      CH4P(N,L,NZ,NY,NX)=CH4P(N,L,NZ,NY,NX)+RCHDFA(N,L,NZ,NY,NX)+RUPCHS(N,L,NZ,NY,NX)
      Z2OP(N,L,NZ,NY,NX)=Z2OP(N,L,NZ,NY,NX)+RN2DFA(N,L,NZ,NY,NX)+RUPN2S(N,L,NZ,NY,NX)
      ZH3P(N,L,NZ,NY,NX)=ZH3P(N,L,NZ,NY,NX)+RNHDFA(N,L,NZ,NY,NX) &
        +RUPN3S(N,L,NZ,NY,NX)+RUPN3B(N,L,NZ,NY,NX)
      H2GP(N,L,NZ,NY,NX)=H2GP(N,L,NZ,NY,NX)+RHGDFA(N,L,NZ,NY,NX)+RUPHGS(N,L,NZ,NY,NX)
!
!     TOTAL ROOT GAS CONTENTS
!
!     TL*P=total root gas content
!     *A,*P=PFT root gaseous, aqueous gas content
!
      TLCO2P(L,NY,NX)=TLCO2P(L,NY,NX)+CO2P(N,L,NZ,NY,NX)+CO2A(N,L,NZ,NY,NX)
      TLOXYP(L,NY,NX)=TLOXYP(L,NY,NX)+OXYP(N,L,NZ,NY,NX)+OXYA(N,L,NZ,NY,NX)
      TLCH4P(L,NY,NX)=TLCH4P(L,NY,NX)+CH4P(N,L,NZ,NY,NX)+CH4A(N,L,NZ,NY,NX)
      TLN2OP(L,NY,NX)=TLN2OP(L,NY,NX)+Z2OP(N,L,NZ,NY,NX)+Z2OA(N,L,NZ,NY,NX)
      TLNH3P(L,NY,NX)=TLNH3P(L,NY,NX)+ZH3P(N,L,NZ,NY,NX)+ZH3A(N,L,NZ,NY,NX)
      TLH2GP(L,NY,NX)=TLH2GP(L,NY,NX)+H2GP(N,L,NZ,NY,NX)+H2GA(N,L,NZ,NY,NX)
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
      TCOFLA(L,NY,NX)=TCOFLA(L,NY,NX)+RCOFLA(N,L,NZ,NY,NX)
      TOXFLA(L,NY,NX)=TOXFLA(L,NY,NX)+ROXFLA(N,L,NZ,NY,NX)
      TCHFLA(L,NY,NX)=TCHFLA(L,NY,NX)+RCHFLA(N,L,NZ,NY,NX)
      TN2FLA(L,NY,NX)=TN2FLA(L,NY,NX)+RN2FLA(N,L,NZ,NY,NX)
      TNHFLA(L,NY,NX)=TNHFLA(L,NY,NX)+RNHFLA(N,L,NZ,NY,NX)
      THGFLA(L,NY,NX)=THGFLA(L,NY,NX)+RHGFLA(N,L,NZ,NY,NX)
      TCO2P(L,NY,NX)=TCO2P(L,NY,NX)-RCO2P(N,L,NZ,NY,NX)
      TUPOXP(L,NY,NX)=TUPOXP(L,NY,NX)+RUPOXP(N,L,NZ,NY,NX)
      TCO2S(L,NY,NX)=TCO2S(L,NY,NX)+RCO2S(N,L,NZ,NY,NX)
      TUPOXS(L,NY,NX)=TUPOXS(L,NY,NX)+RUPOXS(N,L,NZ,NY,NX)
      TUPCHS(L,NY,NX)=TUPCHS(L,NY,NX)+RUPCHS(N,L,NZ,NY,NX)
      TUPN2S(L,NY,NX)=TUPN2S(L,NY,NX)+RUPN2S(N,L,NZ,NY,NX)
      TUPN3S(L,NY,NX)=TUPN3S(L,NY,NX)+RUPN3S(N,L,NZ,NY,NX)
      TUPN3B(L,NY,NX)=TUPN3B(L,NY,NX)+RUPN3B(N,L,NZ,NY,NX)
      TUPHGS(L,NY,NX)=TUPHGS(L,NY,NX)+RUPHGS(N,L,NZ,NY,NX)
      TUPNH4(L,NY,NX)=TUPNH4(L,NY,NX)+RUPNH4(N,L,NZ,NY,NX)
      TUPNO3(L,NY,NX)=TUPNO3(L,NY,NX)+RUPNO3(N,L,NZ,NY,NX)
      TUPH2P(L,NY,NX)=TUPH2P(L,NY,NX)+RUPH2P(N,L,NZ,NY,NX)
      TUPH1P(L,NY,NX)=TUPH1P(L,NY,NX)+RUPH1P(N,L,NZ,NY,NX)
      TUPNHB(L,NY,NX)=TUPNHB(L,NY,NX)+RUPNHB(N,L,NZ,NY,NX)
      TUPNOB(L,NY,NX)=TUPNOB(L,NY,NX)+RUPNOB(N,L,NZ,NY,NX)
      TUPH2B(L,NY,NX)=TUPH2B(L,NY,NX)+RUPH2B(N,L,NZ,NY,NX)
      TUPH1B(L,NY,NX)=TUPH1B(L,NY,NX)+RUPH1B(N,L,NZ,NY,NX)
!
!     TOTAL ROOT C,N,P EXUDATION
!
!     TDFOMC,TDFOMN,TDFOMP=total nonstructl C,N,P exchange
!     RDFOMC,RDFOMN,RDFOMP=PFT nonstructl C,N,P exchange
!
      DO 195 K=0,4
        TDFOMC(K,L,NY,NX)=TDFOMC(K,L,NY,NX)-RDFOMC(N,K,L,NZ,NY,NX)
        TDFOMN(K,L,NY,NX)=TDFOMN(K,L,NY,NX)-RDFOMN(N,K,L,NZ,NY,NX)
        TDFOMP(K,L,NY,NX)=TDFOMP(K,L,NY,NX)-RDFOMP(N,K,L,NZ,NY,NX)
195   CONTINUE
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
      ROXYX(L,NY,NX)=ROXYX(L,NY,NX)+ROXYP(N,L,NZ,NY,NX)
      RNH4X(L,NY,NX)=RNH4X(L,NY,NX)+RUNNHP(N,L,NZ,NY,NX)
      RNO3X(L,NY,NX)=RNO3X(L,NY,NX)+RUNNOP(N,L,NZ,NY,NX)
      RPO4X(L,NY,NX)=RPO4X(L,NY,NX)+RUPP2P(N,L,NZ,NY,NX)
      RP14X(L,NY,NX)=RP14X(L,NY,NX)+RUPP1P(N,L,NZ,NY,NX)
      RNHBX(L,NY,NX)=RNHBX(L,NY,NX)+RUNNBP(N,L,NZ,NY,NX)
      RN3BX(L,NY,NX)=RN3BX(L,NY,NX)+RUNNXP(N,L,NZ,NY,NX)
      RPOBX(L,NY,NX)=RPOBX(L,NY,NX)+RUPP2B(N,L,NZ,NY,NX)
      RP1BX(L,NY,NX)=RP1BX(L,NY,NX)+RUPP1B(N,L,NZ,NY,NX)
    ENDDO
100   CONTINUE
  end subroutine TotalGasandSoluteUptake
!------------------------------------------------------------------------------------------

  subroutine CanopyFluxesandFixation(NZ,NY,NX)
!
!     TOTAL ROOT N2 FIXATION BY ALL PLANT SPECIES
!
!     TUPNF=total root N2 fixation
!     RUPNF=PFT root N2 fixation
!

  implicit none
  integer, intent(in) :: NZ, NY, NX

  DO 85 L=NU(NY,NX),NI(NZ,NY,NX)
    TUPNF(L,NY,NX)=TUPNF(L,NY,NX)+RUPNF(L,NZ,NY,NX)
85  CONTINUE
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
  TRN(NY,NX)=TRN(NY,NX)+RAD1(NZ,NY,NX)
  TLE(NY,NX)=TLE(NY,NX)+EFLXC(NZ,NY,NX)
  TSH(NY,NX)=TSH(NY,NX)+SFLXC(NZ,NY,NX)
  TGH(NY,NX)=TGH(NY,NX)+HFLXC(NZ,NY,NX)
  TCCAN(NY,NX)=TCCAN(NY,NX)+CNET(NZ,NY,NX)
  CTRAN(NZ,NY,NX)=CTRAN(NZ,NY,NX)+EP(NZ,NY,NX)+EVAPC(NZ,NY,NX)
  TVOLWP(NY,NX)=TVOLWP(NY,NX)+VOLWP(NZ,NY,NX)
  TVOLWC(NY,NX)=TVOLWC(NY,NX)+VOLWC(NZ,NY,NX)
  TEVAPP(NY,NX)=TEVAPP(NY,NX)+EP(NZ,NY,NX)+EVAPC(NZ,NY,NX)
  TEVAPC(NY,NX)=TEVAPC(NY,NX)+EVAPC(NZ,NY,NX)
  ENGYC=cpw*(VOLWC(NZ,NY,NX)+FLWC(NZ,NY,NX)+EVAPC(NZ,NY,NX))*TKC(NZ,NY,NX)
  TENGYC(NY,NX)=TENGYC(NY,NX)+ENGYC
  THFLXC(NY,NX)=THFLXC(NY,NX)+ENGYC-ENGYX(NZ,NY,NX)-(FLWC(NZ,NY,NX)*cpw*TKA(NY,NX))
  ENGYX(NZ,NY,NX)=ENGYC
  THRMC(NY,NX)=THRMC(NY,NX)+THRM1(NZ,NY,NX)
  ARLFC(NY,NX)=ARLFC(NY,NX)+ARLFP(NZ,NY,NX)
  ARSTC(NY,NX)=ARSTC(NY,NX)+ARSTP(NZ,NY,NX)
  ZCSNC(NY,NX)=ZCSNC(NY,NX)-HCUPTK(NZ,NY,NX)
  ZZSNC(NY,NX)=ZZSNC(NY,NX)-HZUPTK(NZ,NY,NX)
  ZPSNC(NY,NX)=ZPSNC(NY,NX)-HPUPTK(NZ,NY,NX)
  TBALC=TBALC+BALC(NZ,NY,NX)
  TBALN=TBALN+BALN(NZ,NY,NX)
  TBALP=TBALP+BALP(NZ,NY,NX)
  TCO2Z(NY,NX)=TCO2Z(NY,NX)+RCO2Z(NZ,NY,NX)
  TOXYZ(NY,NX)=TOXYZ(NY,NX)+ROXYZ(NZ,NY,NX)
  TCH4Z(NY,NX)=TCH4Z(NY,NX)+RCH4Z(NZ,NY,NX)
  TN2OZ(NY,NX)=TN2OZ(NY,NX)+RN2OZ(NZ,NY,NX)
  TNH3Z(NY,NX)=TNH3Z(NY,NX)+RNH3Z(NZ,NY,NX)
  TH2GZ(NY,NX)=TH2GZ(NY,NX)+RH2GZ(NZ,NY,NX)
!
!     TOTAL CANOPY NH3 EXCHANGE AND EXUDATION
!
!     RNH3B,RNH3C=PFT NH3 flux between atmosphere and branch,canopy
!     TNH3C=total NH3 flux between atmosphere and canopy
!
  RNH3C(NZ,NY,NX)=0._r8
  DO 80 NB=1,NBR(NZ,NY,NX)
    RNH3C(NZ,NY,NX)=RNH3C(NZ,NY,NX)+RNH3B(NB,NZ,NY,NX)
    TNH3C(NZ,NY,NX)=TNH3C(NZ,NY,NX)+RNH3B(NB,NZ,NY,NX)
80    CONTINUE
  end subroutine CanopyFluxesandFixation

  end module ExtractMod
