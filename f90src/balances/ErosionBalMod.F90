module ErosionBalMod
  use data_kind_mod, only : r8 => SHR_KIND_R8
  use SoilPropertyDataType
  use RootDataType
  use EcoSiMParDataMod, only : micpar
  USE EcoSIMCtrlDataType
  use MicrobialDataType
  USE SOMDataType
  use SedimentDataType
  use GridDataType
  USE AqueChemDatatype
  use FertilizerDataType
implicit none
  private
  character(len=*), parameter :: mod_filename = __FILE__
  public :: SinkChemicals
  contains


!------------------------------------------------------------------------------------------

  subroutine SinkChemicals(NY,NX)
  implicit none
  integer, intent(in) :: NY,NX

  integer :: L,LL,M,N,K,NGL

  real(r8) :: FSINK,FSAN,FSIL,FCLA,FCEC,FAEC
  real(r8) :: FNH4, FNH3,FNHU,FNO3
  real(r8) :: FN4,FHY,FAL,FFE,FCA,FMG,FNA,FKA
  real(r8) :: FHC,FAL2,FFE2,FOH0,FOH1,FOH2
  real(r8) :: FCAC,FCAS,FALP,FFEP
  real(r8) :: FCPD,FCPH,FCPM
  real(r8) :: FOMC,FOMN,FOMP,FH1P,FH2P,FALO,FFEO
  real(r8) :: FORC,FORN,FORP,FOHC,FOHN,FOHP,FOHA
  real(r8) :: FOSC,FOSA,FOSN,FOSP

! begin_execution
! SINK ALL SOLID C,N,P IN POND
!
! BKDS=bulk density
! DLYR=layer depth
! VLS=hourly sinking rate from hour1.f
! FSINK=hourly rate for sediment sinking
!
  D9885: DO L=NL(NY,NX)-1,0,-1
    IF(BKDS(L,NY,NX).LE.ZERO.AND.DLYR(3,L,NY,NX).GT.ZERO)THEN
      !sinking from water to sediment layer
      D9880: DO LL=L+1,NL(NY,NX)
        IF(DLYR(3,LL,NY,NX).GT.ZEROS(NY,NX))exit
      ENDDO D9880

      FSINK=AMIN1(1.0_r8,VLS(NY,NX)/DLYR(3,L,NY,NX))
!
!     SOIL MINERALS
!
!     *ER=sediment flux from erosion.f
!     sediment code:XSED=total,XSAN=sand,XSIL=silt,XCLA=clay
!
      FSAN=FSINK*SAND(L,NY,NX)
      FSIL=FSINK*SILT(L,NY,NX)
      FCLA=FSINK*CLAY(L,NY,NX)
      FCEC=FSINK*XCEC(L,NY,NX)
      FAEC=FSINK*XAEC(L,NY,NX)
      SAND(L,NY,NX)=SAND(L,NY,NX)-FSAN
      SILT(L,NY,NX)=SILT(L,NY,NX)-FSIL
      CLAY(L,NY,NX)=CLAY(L,NY,NX)-FCLA
      XCEC(L,NY,NX)=XCEC(L,NY,NX)-FCEC
      XAEC(L,NY,NX)=XAEC(L,NY,NX)-FAEC
      SAND(LL,NY,NX)=SAND(LL,NY,NX)+FSAN
      SILT(LL,NY,NX)=SILT(LL,NY,NX)+FSIL
      CLAY(LL,NY,NX)=CLAY(LL,NY,NX)+FCLA
      XCEC(LL,NY,NX)=XCEC(LL,NY,NX)+FCEC
      XAEC(LL,NY,NX)=XAEC(LL,NY,NX)+FAEC
!
!     FERTILIZER POOLS
!
!     *ER=sediment flux from erosion.f
!     sediment code:NH4,NH3,NHU,NO3=NH4,NH3,urea,NO3
!
      FNH4=FSINK*ZNH4FA(L,NY,NX)
      FNH3=FSINK*ZNH3FA(L,NY,NX)
      FNHU=FSINK*ZNHUFA(L,NY,NX)
      FNO3=FSINK*ZNO3FA(L,NY,NX)
      ZNH4FA(L,NY,NX)=ZNH4FA(L,NY,NX)-FNH4
      ZNH3FA(L,NY,NX)=ZNH3FA(L,NY,NX)-FNH3
      ZNHUFA(L,NY,NX)=ZNHUFA(L,NY,NX)-FNHU
      ZNO3FA(L,NY,NX)=ZNO3FA(L,NY,NX)-FNO3
      ZNH4FA(LL,NY,NX)=ZNH4FA(LL,NY,NX)+FNH4
      ZNH3FA(LL,NY,NX)=ZNH3FA(LL,NY,NX)+FNH3
      ZNHUFA(LL,NY,NX)=ZNHUFA(LL,NY,NX)+FNHU
      ZNO3FA(LL,NY,NX)=ZNO3FA(LL,NY,NX)+FNO3
!
!     EXCHANGEABLE CATIONS AND ANIONS
!
!     sediment code
!     :XN4=adsorbed NH4
!     :XHY,XAL,XFE,XCA,XMG,XNA,XKA,XHC,AL2,FE2
!      =adsorbed H,Al,Fe,Ca,Mg,Na,K,HCO3,AlOH2,FeOH2
!     :XOH0,XOH1,XOH2=adsorbed R-,R-OH,R-OH2
!     :XH1P,XH2P=adsorbed HPO4,H2PO4
!
      FN4=FSINK*XN4(L,NY,NX)
      FHY=FSINK*XHY(L,NY,NX)
      FAL=FSINK*XAL(L,NY,NX)
      FFE=FSINK*XFE(L,NY,NX)
      FCA=FSINK*XCA(L,NY,NX)
      FMG=FSINK*XMG(L,NY,NX)
      FNA=FSINK*XNA(L,NY,NX)
      FKA=FSINK*XKA(L,NY,NX)
      FHC=FSINK*XHC(L,NY,NX)
      FAL2=FSINK*XALO2(L,NY,NX)
      FFE2=FSINK*XFEO2(L,NY,NX)
      FOH0=FSINK*XOH0(L,NY,NX)
      FOH1=FSINK*XOH1(L,NY,NX)
      FOH2=FSINK*XOH2(L,NY,NX)
      FH1P=FSINK*XH1P(L,NY,NX)
      FH2P=FSINK*XH2P(L,NY,NX)
      XN4(L,NY,NX)=XN4(L,NY,NX)-FN4
      XHY(L,NY,NX)=XHY(L,NY,NX)-FHY
      XAL(L,NY,NX)=XAL(L,NY,NX)-FAL
      XFE(L,NY,NX)=XFE(L,NY,NX)-FFE
      XCA(L,NY,NX)=XCA(L,NY,NX)-FCA
      XMG(L,NY,NX)=XMG(L,NY,NX)-FMG
      XNA(L,NY,NX)=XNA(L,NY,NX)-FNA
      XKA(L,NY,NX)=XKA(L,NY,NX)-FKA
      XHC(L,NY,NX)=XHC(L,NY,NX)-FHC
      XALO2(L,NY,NX)=XALO2(L,NY,NX)-FAL2
      XFEO2(L,NY,NX)=XFEO2(L,NY,NX)-FFE2
      XOH0(L,NY,NX)=XOH0(L,NY,NX)-FOH0
      XOH1(L,NY,NX)=XOH1(L,NY,NX)-FOH1
      XOH2(L,NY,NX)=XOH2(L,NY,NX)-FOH2
      XH1P(L,NY,NX)=XH1P(L,NY,NX)-FH1P
      XH2P(L,NY,NX)=XH2P(L,NY,NX)-FH2P
      XN4(LL,NY,NX)=XN4(LL,NY,NX)+FN4
      XHY(LL,NY,NX)=XHY(LL,NY,NX)+FHY
      XAL(LL,NY,NX)=XAL(LL,NY,NX)+FAL
      XFE(LL,NY,NX)=XFE(LL,NY,NX)+FFE
      XCA(LL,NY,NX)=XCA(LL,NY,NX)+FCA
      XMG(LL,NY,NX)=XMG(LL,NY,NX)+FMG
      XNA(LL,NY,NX)=XNA(LL,NY,NX)+FNA
      XKA(LL,NY,NX)=XKA(LL,NY,NX)+FKA
      XHC(LL,NY,NX)=XHC(LL,NY,NX)+FHC
      XALO2(LL,NY,NX)=XALO2(LL,NY,NX)+FAL2
      XFEO2(LL,NY,NX)=XFEO2(LL,NY,NX)+FFE2
      XOH0(LL,NY,NX)=XOH0(LL,NY,NX)+FOH0
      XOH1(LL,NY,NX)=XOH1(LL,NY,NX)+FOH1
      XOH2(LL,NY,NX)=XOH2(LL,NY,NX)+FOH2
      XH1P(LL,NY,NX)=XH1P(LL,NY,NX)+FH1P
      XH2P(LL,NY,NX)=XH2P(LL,NY,NX)+FH2P
!
!     PRECIPITATES
!
!     sediment code
!     :PALO,PFEO=precip AlOH,FeOH
!     :PCAC,PCAS=precip CaCO3,CaSO4
!     :PALP,PFEP=precip AlPO4,FEPO4
!     :PCPM,PCPD,PCPH=precip CaH2PO4,CaHPO4,apatite
!
      FALO=FSINK*PALOH(L,NY,NX)
      FFEO=FSINK*PFEOH(L,NY,NX)
      FCAC=FSINK*PCACO(L,NY,NX)
      FCAS=FSINK*PCASO(L,NY,NX)
      FALP=FSINK*PALPO(L,NY,NX)
      FFEP=FSINK*PFEPO(L,NY,NX)
      FCPD=FSINK*PCAPD(L,NY,NX)
      FCPH=FSINK*PCAPH(L,NY,NX)
      FCPM=FSINK*PCAPM(L,NY,NX)
      PALOH(L,NY,NX)=PALOH(L,NY,NX)-FALO
      PFEOH(L,NY,NX)=PFEOH(L,NY,NX)-FFEO
      PCACO(L,NY,NX)=PCACO(L,NY,NX)-FCAC
      PCASO(L,NY,NX)=PCASO(L,NY,NX)-FCAS
      PALPO(L,NY,NX)=PALPO(L,NY,NX)-FALP
      PFEPO(L,NY,NX)=PFEPO(L,NY,NX)-FFEP
      PCAPD(L,NY,NX)=PCAPD(L,NY,NX)-FCPD
      PCAPH(L,NY,NX)=PCAPH(L,NY,NX)-FCPH
      PCAPM(L,NY,NX)=PCAPM(L,NY,NX)-FCPM
      PALOH(LL,NY,NX)=PALOH(LL,NY,NX)+FALO
      PFEOH(LL,NY,NX)=PFEOH(LL,NY,NX)+FFEO
      PCACO(LL,NY,NX)=PCACO(LL,NY,NX)+FCAC
      PCASO(LL,NY,NX)=PCASO(LL,NY,NX)+FCAS
      PALPO(LL,NY,NX)=PALPO(LL,NY,NX)+FALP
      PFEPO(LL,NY,NX)=PFEPO(LL,NY,NX)+FFEP
      PCAPD(LL,NY,NX)=PCAPD(LL,NY,NX)+FCPD
      PCAPH(LL,NY,NX)=PCAPH(LL,NY,NX)+FCPH
      PCAPM(LL,NY,NX)=PCAPM(LL,NY,NX)+FCPM
!
!     MICROBIAL C,N,P
!
      D1970: DO K=1,micpar%n_litrsfk

!         OMC,OMN,OMP=microbial C,N,P
!         ORC,ORN,ORP=microbial residue C,N,P
!         OHC,OHN,OHP,OHA=adsorbed C,N,P,acetate
!         OSC,OAA,OSN,OSP=SOC,colonized SOC,SON,SOP
!
        D1960: DO N=1,NFGs
          DO NGL=JGnio(N),JGnfo(N)
            DO M=1,nlbiomcp
              FOMC=FSINK*OMC(M,NGL,K,L,NY,NX)
              FOMN=FSINK*OMN(M,NGL,K,L,NY,NX)
              FOMP=FSINK*OMP(M,NGL,K,L,NY,NX)
              OMC(M,NGL,K,LL,NY,NX)=OMC(M,NGL,K,LL,NY,NX)+FOMC
              OMN(M,NGL,K,LL,NY,NX)=OMN(M,NGL,K,LL,NY,NX)+FOMN
              OMP(M,NGL,K,LL,NY,NX)=OMP(M,NGL,K,LL,NY,NX)+FOMP
              OMC(M,NGL,K,L,NY,NX)=OMC(M,NGL,K,L,NY,NX)-FOMC
              OMN(M,NGL,K,L,NY,NX)=OMN(M,NGL,K,L,NY,NX)-FOMN
              OMP(M,NGL,K,L,NY,NX)=OMP(M,NGL,K,L,NY,NX)-FOMP
            enddo
          enddo
        ENDDO D1960
      ENDDO D1970

!         OMC,OMN,OMP=microbial C,N,P
!         ORC,ORN,ORP=microbial residue C,N,P
!         OHC,OHN,OHP,OHA=adsorbed C,N,P,acetate
!         OSC,OAA,OSN,OSP=SOC,colonized SOC,SON,SOP
!
      DO N=1,NFGs
        DO NGL=JGniA(N),JGnfA(N)
          DO M=1,nlbiomcp
            FOMC=FSINK*OMCff(M,NGL,L,NY,NX)
            FOMN=FSINK*OMNff(M,NGL,L,NY,NX)
            FOMP=FSINK*OMPff(M,NGL,L,NY,NX)
            OMCff(M,NGL,LL,NY,NX)=OMCff(M,NGL,LL,NY,NX)+FOMC
            OMNff(M,NGL,LL,NY,NX)=OMNff(M,NGL,LL,NY,NX)+FOMN
            OMPff(M,NGL,LL,NY,NX)=OMPff(M,NGL,LL,NY,NX)+FOMP
            OMCff(M,NGL,L,NY,NX)=OMCff(M,NGL,L,NY,NX)-FOMC
            OMNff(M,NGL,L,NY,NX)=OMNff(M,NGL,L,NY,NX)-FOMN
            OMPff(M,NGL,L,NY,NX)=OMPff(M,NGL,L,NY,NX)-FOMP
          enddo
        enddo
      ENDDO

!
!     MICROBIAL RESIDUE C,N,P
!
      D1900: DO K=1,micpar%n_litrsfk
        D1940: DO M=1,ndbiomcp
          FORC=FSINK*ORC(M,K,L,NY,NX)
          FORN=FSINK*ORN(M,K,L,NY,NX)
          FORP=FSINK*ORP(M,K,L,NY,NX)
          ORC(M,K,LL,NY,NX)=ORC(M,K,LL,NY,NX)+FORC
          ORN(M,K,LL,NY,NX)=ORN(M,K,LL,NY,NX)+FORN
          ORP(M,K,LL,NY,NX)=ORP(M,K,LL,NY,NX)+FORP
          ORC(M,K,L,NY,NX)=ORC(M,K,L,NY,NX)-FORC
          ORN(M,K,L,NY,NX)=ORN(M,K,L,NY,NX)-FORN
          ORP(M,K,L,NY,NX)=ORP(M,K,L,NY,NX)-FORP
        ENDDO D1940
!
!       ADSORBED C,N,P
!
        FOHC=FSINK*OHC(K,L,NY,NX)
        FOHN=FSINK*OHN(K,L,NY,NX)
        FOHP=FSINK*OHP(K,L,NY,NX)
        FOHA=FSINK*OHA(K,L,NY,NX)
        OHC(K,LL,NY,NX)=OHC(K,LL,NY,NX)+FOHC
        OHN(K,LL,NY,NX)=OHN(K,LL,NY,NX)+FOHN
        OHP(K,LL,NY,NX)=OHP(K,LL,NY,NX)+FOHP
        OHA(K,LL,NY,NX)=OHA(K,LL,NY,NX)+FOHA
        OHC(K,L,NY,NX)=OHC(K,L,NY,NX)-FOHC
        OHN(K,L,NY,NX)=OHN(K,L,NY,NX)-FOHN
        OHP(K,L,NY,NX)=OHP(K,L,NY,NX)-FOHP
        OHA(K,L,NY,NX)=OHA(K,L,NY,NX)-FOHA
!
!       SOC,N,P
!
        D1930: DO M=1,jsken
          FOSC=FSINK*OSC(M,K,L,NY,NX)
          FOSA=FSINK*OSA(M,K,L,NY,NX)
          FOSN=FSINK*OSN(M,K,L,NY,NX)
          FOSP=FSINK*OSP(M,K,L,NY,NX)
          OSC(M,K,LL,NY,NX)=OSC(M,K,LL,NY,NX)+FOSC
          OSA(M,K,LL,NY,NX)=OSA(M,K,LL,NY,NX)+FOSA
          OSN(M,K,LL,NY,NX)=OSN(M,K,LL,NY,NX)+FOSN
          OSP(M,K,LL,NY,NX)=OSP(M,K,LL,NY,NX)+FOSP
          OSC(M,K,L,NY,NX)=OSC(M,K,L,NY,NX)-FOSC
          OSA(M,K,L,NY,NX)=OSA(M,K,L,NY,NX)-FOSA
          OSN(M,K,L,NY,NX)=OSN(M,K,L,NY,NX)-FOSN
          OSP(M,K,L,NY,NX)=OSP(M,K,L,NY,NX)-FOSP
        ENDDO D1930
      ENDDO D1900
    ENDIF
  ENDDO D9885
  end subroutine SinkChemicals

end module ErosionBalMod
