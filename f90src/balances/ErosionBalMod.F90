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
      FCEC=FSINK*trcx_solml(idx_CEC,L,NY,NX)
      FAEC=FSINK*trcx_solml(idx_AEC,L,NY,NX)
      SAND(L,NY,NX)=SAND(L,NY,NX)-FSAN
      SILT(L,NY,NX)=SILT(L,NY,NX)-FSIL
      CLAY(L,NY,NX)=CLAY(L,NY,NX)-FCLA
      trcx_solml(idx_CEC,L,NY,NX)=trcx_solml(idx_CEC,L,NY,NX)-FCEC
      trcx_solml(idx_AEC,L,NY,NX)=trcx_solml(idx_AEC,L,NY,NX)-FAEC
      SAND(LL,NY,NX)=SAND(LL,NY,NX)+FSAN
      SILT(LL,NY,NX)=SILT(LL,NY,NX)+FSIL
      CLAY(LL,NY,NX)=CLAY(LL,NY,NX)+FCLA
      trcx_solml(idx_CEC,LL,NY,NX)=trcx_solml(idx_CEC,LL,NY,NX)+FCEC
      trcx_solml(idx_AEC,LL,NY,NX)=trcx_solml(idx_AEC,LL,NY,NX)+FAEC
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
      FN4=FSINK*trcx_solml(idx_NH4,L,NY,NX)
      FHY=FSINK*trcx_solml(idx_Hp,L,NY,NX)
      FAL=FSINK*trcx_solml(idx_Al,L,NY,NX)
      FFE=FSINK*trcx_solml(idx_Fe,L,NY,NX)
      FCA=FSINK*trcx_solml(idx_Ca,L,NY,NX)
      FMG=FSINK*trcx_solml(idx_Mg,L,NY,NX)
      FNA=FSINK*trcx_solml(idx_Na,L,NY,NX)
      FKA=FSINK*trcx_solml(idx_K,L,NY,NX)
      FHC=FSINK*trcx_solml(idx_COOH,L,NY,NX)
      FAL2=FSINK*trcx_solml(idx_AlOH2,L,NY,NX)
      FFE2=FSINK*trcx_solml(idx_FeOH2,L,NY,NX)
      FOH0=FSINK*trcx_solml(idx_OHe,L,NY,NX)
      FOH1=FSINK*trcx_solml(idx_OH,L,NY,NX)
      FOH2=FSINK*trcx_solml(idx_OHp,L,NY,NX)
      FH1P=FSINK*trcx_solml(idx_HPO4,L,NY,NX)
      FH2P=FSINK*trcx_solml(idx_H2PO4,L,NY,NX)
      trcx_solml(idx_NH4,L,NY,NX)=trcx_solml(idx_NH4,L,NY,NX)-FN4
      trcx_solml(idx_Hp,L,NY,NX)=trcx_solml(idx_Hp,L,NY,NX)-FHY
      trcx_solml(idx_Al,L,NY,NX)=trcx_solml(idx_Al,L,NY,NX)-FAL
      trcx_solml(idx_Fe,L,NY,NX)=trcx_solml(idx_Fe,L,NY,NX)-FFE
      trcx_solml(idx_Ca,L,NY,NX)=trcx_solml(idx_Ca,L,NY,NX)-FCA
      trcx_solml(idx_Mg,L,NY,NX)=trcx_solml(idx_Mg,L,NY,NX)-FMG
      trcx_solml(idx_Na,L,NY,NX)=trcx_solml(idx_Na,L,NY,NX)-FNA
      trcx_solml(idx_K,L,NY,NX)=trcx_solml(idx_K,L,NY,NX)-FKA
      trcx_solml(idx_COOH,L,NY,NX)=trcx_solml(idx_COOH,L,NY,NX)-FHC
      trcx_solml(idx_AlOH2,L,NY,NX)=trcx_solml(idx_AlOH2,L,NY,NX)-FAL2
      trcx_solml(idx_FeOH2,L,NY,NX)=trcx_solml(idx_FeOH2,L,NY,NX)-FFE2
      trcx_solml(idx_OHe,L,NY,NX)=trcx_solml(idx_OHe,L,NY,NX)-FOH0
      trcx_solml(idx_OH,L,NY,NX)=trcx_solml(idx_OH,L,NY,NX)-FOH1
      trcx_solml(idx_OHp,L,NY,NX)=trcx_solml(idx_OHp,L,NY,NX)-FOH2
      trcx_solml(idx_HPO4,L,NY,NX)=trcx_solml(idx_HPO4,L,NY,NX)-FH1P
      trcx_solml(idx_H2PO4,L,NY,NX)=trcx_solml(idx_H2PO4,L,NY,NX)-FH2P
      trcx_solml(idx_NH4,LL,NY,NX)=trcx_solml(idx_NH4,LL,NY,NX)+FN4
      trcx_solml(idx_Hp,LL,NY,NX)=trcx_solml(idx_Hp,LL,NY,NX)+FHY
      trcx_solml(idx_Al,LL,NY,NX)=trcx_solml(idx_Al,LL,NY,NX)+FAL
      trcx_solml(idx_Fe,LL,NY,NX)=trcx_solml(idx_Fe,LL,NY,NX)+FFE
      trcx_solml(idx_Ca,LL,NY,NX)=trcx_solml(idx_Ca,LL,NY,NX)+FCA
      trcx_solml(idx_Mg,LL,NY,NX)=trcx_solml(idx_Mg,LL,NY,NX)+FMG
      trcx_solml(idx_Na,LL,NY,NX)=trcx_solml(idx_Na,LL,NY,NX)+FNA
      trcx_solml(idx_K,LL,NY,NX)=trcx_solml(idx_K,LL,NY,NX)+FKA
      trcx_solml(idx_COOH,LL,NY,NX)=trcx_solml(idx_COOH,LL,NY,NX)+FHC
      trcx_solml(idx_AlOH2,LL,NY,NX)=trcx_solml(idx_AlOH2,LL,NY,NX)+FAL2
      trcx_solml(idx_FeOH2,LL,NY,NX)=trcx_solml(idx_FeOH2,LL,NY,NX)+FFE2
      trcx_solml(idx_OHe,LL,NY,NX)=trcx_solml(idx_OHe,LL,NY,NX)+FOH0
      trcx_solml(idx_OH,LL,NY,NX)=trcx_solml(idx_OH,LL,NY,NX)+FOH1
      trcx_solml(idx_OHp,LL,NY,NX)=trcx_solml(idx_OHp,LL,NY,NX)+FOH2
      trcx_solml(idx_HPO4,LL,NY,NX)=trcx_solml(idx_HPO4,LL,NY,NX)+FH1P
      trcx_solml(idx_H2PO4,LL,NY,NX)=trcx_solml(idx_H2PO4,LL,NY,NX)+FH2P
!
!     PRECIPITATES
!
!     sediment code
!     :PALO,PFEO=precip AlOH,FeOH
!     :PCAC,PCAS=precip CaCO3,CaSO4
!     :PALP,PFEP=precip AlPO4,FEPO4
!     :PCPM,PCPD,PCPH=precip CaH2PO4,CaHPO4,apatite
!
      FALO=FSINK*trcp_salml(idsp_AlOH3,L,NY,NX)
      FFEO=FSINK*trcp_salml(idsp_FeOH3,L,NY,NX)
      FCAC=FSINK*trcp_salml(idsp_CaCO3,L,NY,NX)
      FCAS=FSINK*trcp_salml(idsp_CaSO4,L,NY,NX)
      FALP=FSINK*trcp_salml(idsp_AlPO4,L,NY,NX)
      FFEP=FSINK*trcp_salml(idsp_FePO4,L,NY,NX)
      FCPD=FSINK*trcp_salml(idsp_CaHPO4,L,NY,NX)
      FCPH=FSINK*trcp_salml(idsp_HA,L,NY,NX)
      FCPM=FSINK*trcp_salml(idsp_CaH2PO4,L,NY,NX)
      trcp_salml(idsp_AlOH3,L,NY,NX)=trcp_salml(idsp_AlOH3,L,NY,NX)-FALO
      trcp_salml(idsp_FeOH3,L,NY,NX)=trcp_salml(idsp_FeOH3,L,NY,NX)-FFEO
      trcp_salml(idsp_CaCO3,L,NY,NX)=trcp_salml(idsp_CaCO3,L,NY,NX)-FCAC
      trcp_salml(idsp_CaSO4,L,NY,NX)=trcp_salml(idsp_CaSO4,L,NY,NX)-FCAS
      trcp_salml(idsp_AlPO4,L,NY,NX)=trcp_salml(idsp_AlPO4,L,NY,NX)-FALP
      trcp_salml(idsp_FePO4,L,NY,NX)=trcp_salml(idsp_FePO4,L,NY,NX)-FFEP
      trcp_salml(idsp_CaHPO4,L,NY,NX)=trcp_salml(idsp_CaHPO4,L,NY,NX)-FCPD
      trcp_salml(idsp_HA,L,NY,NX)=trcp_salml(idsp_HA,L,NY,NX)-FCPH
      trcp_salml(idsp_CaH2PO4,L,NY,NX)=trcp_salml(idsp_CaH2PO4,L,NY,NX)-FCPM
      trcp_salml(idsp_AlOH3,LL,NY,NX)=trcp_salml(idsp_AlOH3,LL,NY,NX)+FALO
      trcp_salml(idsp_FeOH3,LL,NY,NX)=trcp_salml(idsp_FeOH3,LL,NY,NX)+FFEO
      trcp_salml(idsp_CaCO3,LL,NY,NX)=trcp_salml(idsp_CaCO3,LL,NY,NX)+FCAC
      trcp_salml(idsp_CaSO4,LL,NY,NX)=trcp_salml(idsp_CaSO4,LL,NY,NX)+FCAS
      trcp_salml(idsp_AlPO4,LL,NY,NX)=trcp_salml(idsp_AlPO4,LL,NY,NX)+FALP
      trcp_salml(idsp_FePO4,LL,NY,NX)=trcp_salml(idsp_FePO4,LL,NY,NX)+FFEP
      trcp_salml(idsp_CaHPO4,LL,NY,NX)=trcp_salml(idsp_CaHPO4,LL,NY,NX)+FCPD
      trcp_salml(idsp_HA,LL,NY,NX)=trcp_salml(idsp_HA,LL,NY,NX)+FCPH
      trcp_salml(idsp_CaH2PO4,LL,NY,NX)=trcp_salml(idsp_CaH2PO4,LL,NY,NX)+FCPM
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
