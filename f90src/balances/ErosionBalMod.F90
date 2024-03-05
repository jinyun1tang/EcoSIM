module ErosionBalMod
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use SoilPropertyDataType
  use RootDataType
  use EcoSiMParDataMod, only : micpar
  USE EcoSIMCtrlDataType
  use EcoSIMCtrlMod
  use MicrobialDataType
  USE SOMDataType
  use SedimentDataType
  use GridDataType
  USE AqueChemDatatype
  use FlagDataType
  use FertilizerDataType
  use EcoSIMConfig, only : nlbiomcp=>NumLiveMicrbCompts
  use TFlxTypeMod
implicit none
  private
  character(len=*), parameter :: mod_filename = &
  __FILE__
  public :: SinkChemicals
  public :: ZeroErosionArray
  contains

!------------------------------------------------------------------------------------------

  subroutine SinkChemicals(NY,NX)
  implicit none
  integer, intent(in) :: NY,NX

  integer :: L,LL,M,N,K,NGL,MID,idom,NE
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
    IF(SoiBulkDensity(L,NY,NX).LE.ZERO.AND.DLYR(3,L,NY,NX).GT.ZERO)THEN
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
!     :PCPM,PCPD,PCPH=precip CaH4P2O8,CaHPO4,apatite
!
! only for non-banded precipitates
      DO NTP=idsp_beg,idsp_beg_band-1
        FPX=FSINK*trcp_salml(NTP,L,NY,NX)
        trcp_salml(NTP,L,NY,NX)=trcp_salml(NTP,L,NY,NX)-FPX
        trcp_salml(NTP,LL,NY,NX)=trcp_salml(NTP,LL,NY,NX)+FPX
      ENDDO

!     MICROBIAL C,N,P
!
      D1970: DO K=1,micpar%NumOfLitrCmplxs

!         OMC,OMN,OMP=microbial C,N,P
!         ORC,ORN,ORP=microbial residue C,N,P
!         OHC,OHN,OHP,OHA=adsorbed C,N,P,acetate
!         OSC,OAA,OSN,OSP=SOC,colonized SOC,SON,SOP
!
        D1960: DO N=1,NumMicbFunGroups
          DO NGL=JGnio(N),JGnfo(N)
            DO M=1,nlbiomcp
              MID=micpar%get_micb_id(M,NGL)
              DO NE=1,NumPlantChemElms
                FOMC=FSINK*OMEhetr(NE,MID,K,L,NY,NX)
                OMEhetr(NE,MID,K,LL,NY,NX)=OMEhetr(NE,MID,K,LL,NY,NX)+FOMC                
                OMEhetr(NE,MID,K,L,NY,NX)=OMEhetr(NE,MID,K,L,NY,NX)-FOMC
              ENDDO
            enddo
          enddo
        ENDDO D1960
      ENDDO D1970

!         OMC,OMN,OMP=microbial C,N,P
!         ORC,ORN,ORP=microbial residue C,N,P
!         OHC,OHN,OHP,OHA=adsorbed C,N,P,acetate
!         OSC,OAA,OSN,OSP=SOC,colonized SOC,SON,SOP
!
      DO N=1,NumMicbFunGroups
        DO NGL=JGniA(N),JGnfA(N)
          DO M=1,nlbiomcp
            MID=micpar%get_micb_id(M,NGL)        
            DO NE=1,NumPlantChemElms              
              FOMC=FSINK*OMEauto(NE,MID,L,NY,NX)
              OMEauto(NE,MID,LL,NY,NX)=OMEauto(NE,MID,LL,NY,NX)+FOMC
              OMEauto(NE,MID,L,NY,NX)=OMEauto(NE,MID,L,NY,NX)-FOMC
            ENDDO
          enddo
        enddo
      ENDDO

!
!     MICROBIAL RESIDUE C,N,P
!
      D1900: DO K=1,micpar%NumOfLitrCmplxs
        D1940: DO M=1,ndbiomcp
          DO NE=1,NumPlantChemElms                   
            FORC=FSINK*ORM(NE,M,K,L,NY,NX)
            ORM(NE,M,K,LL,NY,NX)=ORM(NE,M,K,LL,NY,NX)+FORC
            ORM(NE,M,K,L,NY,NX)=ORM(NE,M,K,L,NY,NX)-FORC
          ENDDO
        ENDDO D1940
!
!       ADSORBED C,N,P
!
        DO idom=idom_beg,idom_end
          FOHC=FSINK*OHM(idom,K,L,NY,NX)
          OHM(idom,K,LL,NY,NX)=OHM(idom,K,LL,NY,NX)+FOHC
          OHM(idom,K,L,NY,NX)=OHM(idom,K,L,NY,NX)-FOHC
        ENDDO
!
!       SOC,N,P
!
        D1930: DO M=1,jsken
          FOSA=FSINK*OSA(M,K,L,NY,NX)
          OSA(M,K,LL,NY,NX)=OSA(M,K,LL,NY,NX)+FOSA
          OSA(M,K,L,NY,NX)=OSA(M,K,L,NY,NX)-FOSA
          DO NE=1,NumPlantChemElms
            FOSC=FSINK*OSM(NE,M,K,L,NY,NX)
            OSM(NE,M,K,LL,NY,NX)=OSM(NE,M,K,LL,NY,NX)+FOSC
            OSM(NE,M,K,L,NY,NX)=OSM(NE,M,K,L,NY,NX)-FOSC
          ENDDO
        ENDDO D1930
      ENDDO D1900
    ENDIF
  ENDDO D9885
  end subroutine SinkChemicals
!------------------------------------------------------------------------------------------

  subroutine ZeroErosionArray(NY,NX)
  implicit none
  integer, intent(in) :: NY,NX

  IF(iErosionMode.EQ.ieros_frzthaweros.OR.iErosionMode.EQ.ieros_frzthawsomeros)THEN
    tErosionSedmLoss(NY,NX)=0.0_r8
    TSANER(NY,NX)=0.0_r8
    TSILER(NY,NX)=0.0_r8
    TCLAER(NY,NX)=0.0_r8
    TNH4ER(NY,NX)=0.0_r8
    TNH3ER(NY,NX)=0.0_r8
    TNHUER(NY,NX)=0.0_r8
    TNO3ER(NY,NX)=0.0_r8
    TNH4EB(NY,NX)=0.0_r8
    TNH3EB(NY,NX)=0.0_r8
    TNHUEB(NY,NX)=0.0_r8
    TNO3EB(NY,NX)=0.0_r8

    trcx_TER(idx_beg:idx_end,NY,NX)=0.0_r8
    trcp_TER(idsp_beg:idsp_end,NY,NX)=0.0_r8

    TOMEERhetr(1:NumPlantChemElms,1:NumLiveHeterBioms,1:jcplx,NY,NX)=0.0_r8
    TOMEERauto(1:NumPlantChemElms,1:NumLiveAutoBioms,NY,NX)=0.0_r8

    TORMER(1:NumPlantChemElms,1:ndbiomcp,1:jcplx,NY,NX)=0.0_r8
    TOHMER(idom_beg:idom_end,1:jcplx,NY,NX)=0.0_r8
    TOSMER(1:NumPlantChemElms,1:jsken,1:jcplx,NY,NX)=0.0_r8
    TOSAER(1:jsken,1:jcplx,NY,NX)=0.0_r8
  ENDIF
  end subroutine ZeroErosionArray

end module ErosionBalMod
