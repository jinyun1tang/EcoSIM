module ErosionBalMod
  use data_kind_mod, only : r8 => DAT_KIND_R8
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
  real(r8) :: FNX
  real(r8) :: FOMC,FOMN,FOMP
  real(r8) :: FORC,FORN,FORP,FOHC,FOHN,FOHP,FOHA
  real(r8) :: FOSC,FOSA,FOSN,FOSP
  integer  :: NTF,NTP,NTX
  real(r8) :: FPX,FSNX
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
      SAND(L,NY,NX)=SAND(L,NY,NX)-FSAN
      SILT(L,NY,NX)=SILT(L,NY,NX)-FSIL
      CLAY(L,NY,NX)=CLAY(L,NY,NX)-FCLA
      SAND(LL,NY,NX)=SAND(LL,NY,NX)+FSAN
      SILT(LL,NY,NX)=SILT(LL,NY,NX)+FSIL
      CLAY(LL,NY,NX)=CLAY(LL,NY,NX)+FCLA
!
!     FERTILIZER POOLS
!
!     *ER=sediment flux from erosion.f
!     sediment code:NH4,NH3,NHU,NO3=NH4,NH3,urea,NO3
!
      DO NTF=ifertn_beg,ifertn_end
        FNX=FSINK*FertN_soil(NTF,L,NY,NX)
        FertN_soil(NTF,L,NY,NX)=FertN_soil(NTF,L,NY,NX)-FNX
        FertN_soil(NTF,LL,NY,NX)=FertN_soil(NTF,LL,NY,NX)+FNX
      ENDDO
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
! exclude banded nutrients
      DO NTX=idx_CEC,idx_cation_soil_end
        FSNX=FSINK*trcx_solml(NTX,L,NY,NX)
        trcx_solml(NTX,L,NY,NX)=trcx_solml(NTX,L,NY,NX)-FSNX
        trcx_solml(NTX,LL,NY,NX)=trcx_solml(NTX,LL,NY,NX)+FSNX
      ENDDO

      DO NTX=idx_AEC,idx_anion_soil_end
        FSNX=FSINK*trcx_solml(NTX,L,NY,NX)
        trcx_solml(NTX,L,NY,NX)=trcx_solml(NTX,L,NY,NX)-FSNX
        trcx_solml(NTX,LL,NY,NX)=trcx_solml(NTX,LL,NY,NX)+FSNX
      ENDDO

!
!     PRECIPITATES
!
!     sediment code
!     :PALO,PFEO=precip AlOH,FeOH
!     :PCAC,PCAS=precip CaCO3,CaSO4
!     :PALP,PFEP=precip AlPO4,FEPO4
!     :PCPM,PCPD,PCPH=precip CaH2PO4,CaHPO4,apatite
!
! only for non-banded precipitates
      DO NTP=idsp_beg,idsp_beg_band-1
        FPX=FSINK*trcp_salml(NTP,L,NY,NX)
        trcp_salml(NTP,L,NY,NX)=trcp_salml(NTP,L,NY,NX)-FPX
        trcp_salml(NTP,LL,NY,NX)=trcp_salml(NTP,LL,NY,NX)+FPX
      ENDDO

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
