module TillageMixMod
  use data_kind_mod, only : r8 => SHR_KIND_R8
  use EcoSiMParDataMod, only : micpar
  use EcosimConst
  use SOMDataType
  USE GridDataType
  USE GridConsts
  USE SoilBGCDataType
  use FertilizerDataType
  use AqueChemDatatype
  USE MicrobialDataType
  use SurfLitterDataType
  use SoilHeatDataType
  use SoilWaterDataType
  use SoilPhysDataType
  USE SoilPropertyDataType
  USE ClimForcDataType
  USE FlagDataType
  use EcoSIMCtrlDataType
  USE EcoSimSumDataType
  use EcoSIMConfig , only : ndbiomcp => ndbiomcpc
  implicit none
  character(len=*),private, parameter :: mod_filename = __FILE__
  public :: ApplyTillageMixing
  contains

!------------------------------------------------------------------------------------------

  subroutine ApplyTillageMixing(I,J,NY,NX)
  implicit none
  integer, intent(in) :: I,J,NY,NX

  integer :: K,N,M,L,LL,NGL,NTS

  REAL(R8) :: TOSGC(jsken,0:micpar%n_litrsfk),TOSGA(jsken,0:micpar%n_litrsfk)
  real(r8) :: TOSGN(jsken,0:micpar%n_litrsfk),TOSGP(jsken,0:micpar%n_litrsfk)
  real(r8) :: TORXC(ndbiomcp,0:micpar%n_litrsfk)
  real(r8) :: TORXN(ndbiomcp,0:micpar%n_litrsfk),TORXP(ndbiomcp,0:micpar%n_litrsfk)
  real(r8) :: TOQGC(0:micpar%n_litrsfk),TOQGN(0:micpar%n_litrsfk)
  real(r8) :: TOQGP(0:micpar%n_litrsfk),TOQHC(0:micpar%n_litrsfk)
  real(r8) :: TOQHN(0:micpar%n_litrsfk),TOQHP(0:micpar%n_litrsfk)
  real(r8) :: TOHGC(0:micpar%n_litrsfk),TOHGN(0:micpar%n_litrsfk)
  real(r8) :: TOHGP(0:micpar%n_litrsfk),TOHGA(0:micpar%n_litrsfk)
  real(r8) :: TOQGA(0:micpar%n_litrsfk),TOQHA(0:micpar%n_litrsfk)

  real(r8) :: TOMGCff(nlbiomcp,NMICBSA)
  real(r8) :: TOMGNff(nlbiomcp,NMICBSA)
  real(r8) :: TOMGPff(nlbiomcp,NMICBSA)
  real(r8) :: TOMGC(nlbiomcp,NMICBSO,1:jcplx)
  real(r8) :: TOMGN(nlbiomcp,NMICBSO,1:jcplx)
  real(r8) :: TOMGP(nlbiomcp,NMICBSO,1:jcplx)
  REAL(R8) :: TOMC(nlbiomcp,NMICBSO,1:jcplx)
  REAL(R8) :: TOMN(nlbiomcp,NMICBSO,1:jcplx)
  REAL(R8) :: TOMP(nlbiomcp,NMICBSO,1:jcplx)
  REAL(R8) :: TOMCff(nlbiomcp,NMICBSA)
  REAL(R8) :: TOMNff(nlbiomcp,NMICBSA)
  REAL(R8) :: TOMPff(nlbiomcp,NMICBSA)
  real(r8) :: TORC(ndbiomcp,1:jcplx)
  real(r8) :: TORN(ndbiomcp,1:jcplx)
  real(r8) :: TORP(ndbiomcp,1:jcplx)
  real(r8) :: TOQC(1:jcplx)
  real(r8) :: TOQN(1:jcplx)
  real(r8) :: TOQP(1:jcplx)
  real(r8) :: TOQA(1:jcplx)
  real(r8) :: TOHC(1:jcplx)
  real(r8) :: TOHN(1:jcplx)
  real(r8) :: TOHP(1:jcplx)
  real(r8) :: TOHA(1:jcplx)
  real(r8) :: TOSC(jsken,1:jcplx)
  real(r8) :: TOSA(jsken,1:jcplx)
  real(r8) :: TOSN(jsken,1:jcplx)
  real(r8) :: TOSP(jsken,1:jcplx)
  real(r8) :: CORP,CORP0,XCORP0,DCORPZ
  real(r8) :: TBKDX,TFC,TWP,TSCNV,TSCNH,TSAND
  REAL(R8) :: ZNHUX0,ZNHUXI,ZNFNX0
  real(r8) :: TSILT,TCLAY,TGKC4,TGKCA,TGKCM,TGKCN,TGKCK
  real(r8) :: TP_soil(idsp_psoi_beg:idsp_psoi_end)
  real(r8) :: TFertNG_soil(ifertn_beg:ifertn_end),TZNFNG,TVOLWR
  real(r8) :: TZ2OG,TZ2OS,TZNH3G,TH2GG,TH2GS,TZNFN2,TZNFNI,TCO2GS
  real(r8) :: TCH4GS,TOXYGS,TZ2GSG,TZ2OGS,TH2GGS,TNH4GS,TNH3GS
  real(r8) :: TX_anion(idx_AEC+1:idx_anion_soil_end)
  real(r8) :: TNO3GS,TNO2GS,TP14GS,TPO4GS,TXN4G
  real(r8) :: TX_solml(idx_beg:idx_end)
  real(r8) :: TP_salml(idsp_beg:idsp_end)
  real(r8) :: TG_gasml(idg_beg:idg_end-1)
  real(r8) :: TSA_solml(idsa_beg:idsa_end)
  real(r8) :: TS_solml(ids_beg:ids_end)
  real(r8) :: TS0_solml(ids_beg:ids_end)   !surface mass for incoporation
  real(r8) :: TfertN_band(ifertn_beg:ifertnb_end)
  real(r8) :: TNH4S,TNH4B,TNH3S,TNH3B,TNO3S
  real(r8) :: TNO3B,TNO2S,TNO2B,TZAL,TZFE,TZHY,TZCA,TZMG,TZNA
  real(r8) :: TNFNIH,TfertN_soil(ifertn_beg:ifertn_end)
  real(r8) :: TENGYR,TVOLW,TENGY
  real(r8) :: ENGYV,ENGYL,ENGYM
  real(r8) :: DC,DN,DP,OC,ON,OP
  real(r8) :: TVOLI,HFLXD
  real(r8) :: FI,TI,TX,TL
  integer  :: NTX,NTP,NTG,NTSA,NTN
!     begin_execution
!
  IF(J.EQ.INT(ZNOON(NY,NX)).AND.XCORP(NY,NX).LT.1.0.AND.DCORP(I,NY,NX).GT.0.0)THEN
!
!     EXTENT OF MIXING
!
    IFLGS(NY,NX)=1
    CORP=1.0_r8-XCORP(NY,NX)
    ENGYP(NY,NX)=0.0_r8
!
!     TEMPORARY ACCUMULATORS
!
    TBKDX=0.0_r8
    TFC=0.0_r8
    TWP=0.0_r8
    TSCNV=0.0_r8
    TSCNH=0.0_r8
    TSAND=0.0_r8
    TSILT=0.0_r8
    TCLAY=0.0_r8
    TGKC4=0.0_r8
    TGKCA=0.0_r8
    TGKCM=0.0_r8
    TGKCN=0.0_r8
    TGKCK=0.0_r8
    TVOLW=0.0_r8
    TVOLI=0.0_r8
!     TVOLP=0.0_r8
!     TVOLA=0.0_r8
    TENGY=0.0_r8
    TNFNIH=0.0_r8

    TfertN_soil(ifertn_beg:ifertn_end)=0.0_r8
    TfertN_band(ifertnb_beg:ifertnb_end)=0.0_r8

    TS_solml(ids_beg:ids_end)=0._r8
    TS0_solml(ids_beg:ids_end)=0._r8
    TX_solml(idx_beg:idx_end)=0._r8
    TP_salml(idsp_beg:idsp_end)=0._r8
    TG_gasml(idg_beg:idg_end-1)=0._r8

    TOMC=0.0_r8
    TOMN=0.0_r8
    TOMP=0.0_r8
    TOMCff=0.0_r8
    TOMNff=0.0_r8
    TOMPff=0.0_r8
    TORC=0.0_r8
    TORN=0.0_r8
    TORP=0.0_r8
    TOQC=0.0_r8
    TOQN=0.0_r8
    TOQP=0.0_r8
    TOQA=0.0_r8
    TOHC=0.0_r8
    TOHN=0.0_r8
    TOHP=0.0_r8
    TOHA=0.0_r8
    TOSC=0.0_r8
    TOSA=0.0_r8
    TOSN=0.0_r8
    TOSP=0.0_r8
    TZNFN2=0.0_r8
    TZNFNI=0.0_r8
    ZNHUX0=0.0_r8
    ZNHUXI=0.0_r8
    ZNFNX0=0.0_r8
!
!     ACCUMULATE STATE VARIABLES IN SURFACE RESIDUE FOR ADDITION
!     TO SOIL IN TILLAGE MIXING ZONE
!
!     IF(ORGC(0,NY,NX).GT.ZEROS(NY,NX))THEN
!     XCORP0=AMAX1(XCORP(NY,NX),AMIN1(1.0,
!    2(VHCPRX(NY,NX)/cpo)/ORGC(0,NY,NX)))
    XCORP0=AMAX1(0.001_r8,XCORP(NY,NX))
!     ELSE
!     XCORP0=1.0
!     ENDIF
    CORP0=1.0_r8-XCORP0   !fraction mixed into next layer

    DC=0.0_r8
    DN=0.0_r8
    DP=0.0_r8

    DO  K=1,micpar%n_litrsfk
      DO  N=1,NFGs
        DO NGL=JGnio(N),JGnfo(N)
          DO  M=1,nlbiomcp
            TOMGC(M,NGL,K)=OMC(M,NGL,K,0,NY,NX)*CORP0
            TOMGN(M,NGL,K)=OMN(M,NGL,K,0,NY,NX)*CORP0
            TOMGP(M,NGL,K)=OMP(M,NGL,K,0,NY,NX)*CORP0
            OMC(M,NGL,K,0,NY,NX)=OMC(M,NGL,K,0,NY,NX)*XCORP0
            OMN(M,NGL,K,0,NY,NX)=OMN(M,NGL,K,0,NY,NX)*XCORP0
            OMP(M,NGL,K,0,NY,NX)=OMP(M,NGL,K,0,NY,NX)*XCORP0
            DC=DC+OMC(M,NGL,K,0,NY,NX)
            DN=DN+OMN(M,NGL,K,0,NY,NX)
            DP=DP+OMP(M,NGL,K,0,NY,NX)
          enddo
        enddo
      ENDDO
    ENDDO

    DO  N=1,NFGs
      DO NGL=JGniA(N),JGnfA(N)
        DO  M=1,nlbiomcp
          TOMGCff(M,NGL)=OMCff(M,NGL,0,NY,NX)*CORP0
          TOMGNff(M,NGL)=OMNff(M,NGL,0,NY,NX)*CORP0
          TOMGPff(M,NGL)=OMPff(M,NGL,0,NY,NX)*CORP0
          OMCff(M,NGL,0,NY,NX)=OMCff(M,NGL,0,NY,NX)*XCORP0
          OMNff(M,NGL,0,NY,NX)=OMNff(M,NGL,0,NY,NX)*XCORP0
          OMPff(M,NGL,0,NY,NX)=OMPff(M,NGL,0,NY,NX)*XCORP0
          DC=DC+OMCff(M,NGL,0,NY,NX)
          DN=DN+OMNff(M,NGL,0,NY,NX)
          DP=DP+OMPff(M,NGL,0,NY,NX)
        enddo
      enddo
    ENDDO

    DO K=1,micpar%n_litrsfk
      DO M=1,ndbiomcp
        TORXC(M,K)=ORC(M,K,0,NY,NX)*CORP0
        TORXN(M,K)=ORN(M,K,0,NY,NX)*CORP0
        TORXP(M,K)=ORP(M,K,0,NY,NX)*CORP0
        ORC(M,K,0,NY,NX)=ORC(M,K,0,NY,NX)*XCORP0
        ORN(M,K,0,NY,NX)=ORN(M,K,0,NY,NX)*XCORP0
        ORP(M,K,0,NY,NX)=ORP(M,K,0,NY,NX)*XCORP0
        DC=DC+ORC(M,K,0,NY,NX)
        DN=DN+ORN(M,K,0,NY,NX)
        DP=DP+ORP(M,K,0,NY,NX)
      ENDDO
      TOQGC(K)=OQC(K,0,NY,NX)*CORP0
      TOQGN(K)=OQN(K,0,NY,NX)*CORP0
      TOQGP(K)=OQP(K,0,NY,NX)*CORP0
      TOQGA(K)=OQA(K,0,NY,NX)*CORP0
      TOQHC(K)=OQCH(K,0,NY,NX)*CORP0
      TOQHN(K)=OQNH(K,0,NY,NX)*CORP0
      TOQHP(K)=OQPH(K,0,NY,NX)*CORP0
      TOQHA(K)=OQAH(K,0,NY,NX)*CORP0
      TOHGC(K)=OHC(K,0,NY,NX)*CORP0
      TOHGN(K)=OHN(K,0,NY,NX)*CORP0
      TOHGP(K)=OHP(K,0,NY,NX)*CORP0
      TOHGA(K)=OHA(K,0,NY,NX)*CORP0
  !
  !     REDUCE SURFACE RESIDUE STATE VARIABLES FOR INCORPORATION
  !
      OQC(K,0,NY,NX)=OQC(K,0,NY,NX)*XCORP0
      OQN(K,0,NY,NX)=OQN(K,0,NY,NX)*XCORP0
      OQP(K,0,NY,NX)=OQP(K,0,NY,NX)*XCORP0
      OQA(K,0,NY,NX)=OQA(K,0,NY,NX)*XCORP0
      OQCH(K,0,NY,NX)=OQCH(K,0,NY,NX)*XCORP0
      OQNH(K,0,NY,NX)=OQNH(K,0,NY,NX)*XCORP0
      OQPH(K,0,NY,NX)=OQPH(K,0,NY,NX)*XCORP0
      OQAH(K,0,NY,NX)=OQAH(K,0,NY,NX)*XCORP0
      OHC(K,0,NY,NX)=OHC(K,0,NY,NX)*XCORP0
      OHN(K,0,NY,NX)=OHN(K,0,NY,NX)*XCORP0
      OHP(K,0,NY,NX)=OHP(K,0,NY,NX)*XCORP0
      OHA(K,0,NY,NX)=OHA(K,0,NY,NX)*XCORP0
      DC=DC+OQC(K,0,NY,NX)+OQCH(K,0,NY,NX)+OHC(K,0,NY,NX) &
        +OQA(K,0,NY,NX)+OQAH(K,0,NY,NX)+OHA(K,0,NY,NX)
      DN=DN+OQN(K,0,NY,NX)+OQNH(K,0,NY,NX)+OHN(K,0,NY,NX)
      DP=DP+OQP(K,0,NY,NX)+OQPH(K,0,NY,NX)+OHP(K,0,NY,NX)
      DO  M=1,jsken
        TOSGC(M,K)=OSC(M,K,0,NY,NX)*CORP0
        TOSGA(M,K)=OSA(M,K,0,NY,NX)*CORP0
        TOSGN(M,K)=OSN(M,K,0,NY,NX)*CORP0
        TOSGP(M,K)=OSP(M,K,0,NY,NX)*CORP0
        OSC(M,K,0,NY,NX)=OSC(M,K,0,NY,NX)*XCORP0
        OSA(M,K,0,NY,NX)=OSA(M,K,0,NY,NX)*XCORP0
        OSN(M,K,0,NY,NX)=OSN(M,K,0,NY,NX)*XCORP0
        OSP(M,K,0,NY,NX)=OSP(M,K,0,NY,NX)*XCORP0
        DC=DC+OSC(M,K,0,NY,NX)
        DN=DN+OSN(M,K,0,NY,NX)
        DP=DP+OSP(M,K,0,NY,NX)
      ENDDO
    ENDDO

    DO NTS=ids_beg,idg_NH3
      TS0_solml(NTS)=trc_solml(NTS,0,NY,NX)*CORP0
    ENDDO
    DO NTS=ids_nut_beg,ids_nuts_end
      TS0_solml(NTS)=trc_solml(NTS,0,NY,NX)*CORP0
    ENDDO

    TXN4G=trcx_solml(idx_NH4,0,NY,NX)*CORP0

    DO NTX= idx_AEC+1,idx_anion_soil_end
      TX_anion(NTX)=trcx_solml(NTX,0,NY,NX)*CORP0
    ENDDO

    DO NTP=idsp_psoi_beg,idsp_psoi_end
      TP_soil(NTP)=trcp_salml(NTP,0,NY,NX)*CORP0
    ENDDO

    DO NTN=ifertn_beg,ifertn_end
      TFertNG_soil(NTN)=FertN_soil(NTN,0,NY,NX)*CORP0
    ENDDO

    TZNFNG=ZNFNI(0,NY,NX)*CORP0
    TVOLWR=VOLW(0,NY,NX)*CORP0
    HFLXD=cpo*ORGC(0,NY,NX)*CORP0*TKS(0,NY,NX)
    HEATIN=HEATIN-HFLXD
    HEATSO=HEATSO-HFLXD
    TENGYR=cpw*TVOLWR*TKS(0,NY,NX)
    ORGC(0,NY,NX)=DC
    ORGN(0,NY,NX)=DN
    ORGR(0,NY,NX)=DC

    DO NTS=ids_beg,idg_NH3
      trc_solml(NTS,0,NY,NX)=trc_solml(NTS,0,NY,NX)*XCORP0
    ENDDO

    DO NTS=ids_nut_beg,ids_nuts_end
      trc_solml(NTS,0,NY,NX)=trc_solml(NTS,0,NY,NX)*XCORP0
    ENDDO

    trcx_solml(idx_NH4,0,NY,NX)=trcx_solml(idx_NH4,0,NY,NX)*XCORP0

    DO NTX=idx_AEC+1,idx_anion_soil_end
      trcx_solml(idx_OHe,0,NY,NX)=trcx_solml(idx_OHe,0,NY,NX)*XCORP0
    ENDDO

    trcp_salml(idsp_AlPO4,0,NY,NX)=trcp_salml(idsp_AlPO4,0,NY,NX)*XCORP0
    trcp_salml(idsp_FePO4,0,NY,NX)=trcp_salml(idsp_FePO4,0,NY,NX)*XCORP0
    trcp_salml(idsp_CaHPO4,0,NY,NX)=trcp_salml(idsp_CaHPO4,0,NY,NX)*XCORP0
    trcp_salml(idsp_HA,0,NY,NX)=trcp_salml(idsp_HA,0,NY,NX)*XCORP0
    trcp_salml(idsp_CaH2PO4,0,NY,NX)=trcp_salml(idsp_CaH2PO4,0,NY,NX)*XCORP0

    DO NTN=ifertn_beg,ifertn_end
      FertN_soil(NTN,0,NY,NX)=FertN_soil(NTN,0,NY,NX)*XCORP0
    ENDDO

    VOLW(0,NY,NX)=VOLW(0,NY,NX)*XCORP0
    VHCP(0,NY,NX)=cpo*ORGC(0,NY,NX)+cpw*VOLW(0,NY,NX)+cpi*VOLI(0,NY,NX)
    VOLR(NY,NX)=VOLR(NY,NX)*XCORP0
    VOLT(0,NY,NX)=VOLT(0,NY,NX)*XCORP0
    ZNHUX0=AMAX1(ZNHUX0,ZNHU0(0,NY,NX))
    ZNHUXI=AMAX1(ZNHUXI,ZNHUI(0,NY,NX))
    ZNFNX0=AMAX1(ZNFNX0,ZNFN0(0,NY,NX))
    LL=NU(NY,NX)
!
!     REDISTRIBUTE SOIL STATE VARIABLES DURING TILLAGE
!
    DCORPZ=AMIN1(DCORP(I,NY,NX),CDPTHZ(NL(NY,NX),NY,NX))
!
!     ACCUMULATE SOIL STATE VARIABLES WITHIN TILLAGE MIXING ZONE
!
    DO L=NU(NY,NX),NL(NY,NX)
      IF(CDPTHZ(L,NY,NX)-DLYR(3,L,NY,NX).LT.DCORPZ.AND.DLYR(3,L,NY,NX).GT.ZERO)THEN
        TL=AMIN1(DLYR(3,L,NY,NX),DCORPZ-(CDPTHZ(L,NY,NX)-DLYR(3,L,NY,NX)))
        FI=TL/DCORPZ
        TI=TL/DLYR(3,L,NY,NX)
        TBKDX=TBKDX+FI*BKDSI(L,NY,NX)
        TFC=TFC+FI*FC(L,NY,NX)
        TWP=TWP+FI*WP(L,NY,NX)
        TSCNV=TSCNV+FI*SCNV(L,NY,NX)
        TSCNH=TSCNH+FI*SCNH(L,NY,NX)
        TSAND=TSAND+TI*SAND(L,NY,NX)
        TSILT=TSILT+TI*SILT(L,NY,NX)
        TCLAY=TCLAY+TI*CLAY(L,NY,NX)
        TGKC4=TGKC4+FI*GKC4(L,NY,NX)
        TGKCA=TGKCA+FI*GKCA(L,NY,NX)
        TGKCM=TGKCM+FI*GKCM(L,NY,NX)
        TGKCN=TGKCN+FI*GKCN(L,NY,NX)
        TGKCK=TGKCK+FI*GKCK(L,NY,NX)
        TVOLW=TVOLW+TI*VOLW(L,NY,NX)
        TVOLI=TVOLI+TI*VOLI(L,NY,NX)
!     TVOLP=TVOLP+TI*VOLP(L,NY,NX)
!     TVOLA=TVOLA+TI*VOLA(L,NY,NX)
        TENGY=TENGY+TI*(cpw*(VOLW(L,NY,NX)+VOLWH(L,NY,NX)) &
          +cpi*(VOLI(L,NY,NX)+VOLIH(L,NY,NX)))*TKS(L,NY,NX)
        DO NTN=ifertn_beg,ifertn_end
          TfertN_soil(NTN)=TfertN_soil(NTN)+TI*FertN_soil(NTN,L,NY,NX)
        ENDDO

        DO NTN=ifertnb_beg,ifertnb_end
          TfertN_band(NTN)=TfertN_band(NTN)+TI*FertN_band(NTN,L,NY,NX)
        ENDDO

        DO NTS=ids_beg,ids_end
          TS_solml(NTS)=TS_solml(NTS)+TI*trc_solml(NTS,L,NY,NX)
        ENDDO

        DO NTSA=idsa_beg,idsa_end
          TSA_solml(NTSA)=TSA_solml(NTSA)+TI*trcsa_solml(NTSA,L,NY,NX)
        ENDDO
!cation
        DO NTX=idx_CEC,idx_cation_end
          TX_solml(NTX)=TX_solml(NTX)+TI*trcx_solml(NTX,L,NY,NX)
        ENDDO
!anion
        DO NTX=idx_AEC,idx_end
          TX_solml(NTX)=TX_solml(NTX)+TI*trcx_solml(NTX,L,NY,NX)
        ENDDO

        DO NTP=idsp_beg,idsp_end
          TP_salml(NTP)=TP_salml(NTP)+TI*trcp_salml(NTP,L,NY,NX)
        ENDDO

        DO NTG=idg_beg,idg_end-1
          TG_gasml(NTG)=TG_gasml(NTG)+TI*trc_gasml(NTG,L,NY,NX)
        ENDDO

        DO  K=1,jcplx
          DO  N=1,NFGs
            DO NGL=JGnio(N),JGnfo(N)
              DO  M=1,nlbiomcp
                TOMC(M,NGL,K)=TOMC(M,NGL,K)+TI*OMC(M,NGL,K,L,NY,NX)
                TOMN(M,NGL,K)=TOMN(M,NGL,K)+TI*OMN(M,NGL,K,L,NY,NX)
                TOMP(M,NGL,K)=TOMP(M,NGL,K)+TI*OMP(M,NGL,K,L,NY,NX)
              enddo
            enddo
          enddo
        ENDDO

        DO  N=1,NFGs
          DO NGL=JGniA(N),JGnfA(N)
            DO  M=1,nlbiomcp
              TOMCff(M,NGL)=TOMCff(M,NGL)+TI*OMCff(M,NGL,L,NY,NX)
              TOMNff(M,NGL)=TOMNff(M,NGL)+TI*OMNff(M,NGL,L,NY,NX)
              TOMPff(M,NGL)=TOMPff(M,NGL)+TI*OMPff(M,NGL,L,NY,NX)
            enddo
          enddo
        enddo

      DO K=1,jcplx
        DO  M=1,ndbiomcp
          TORC(M,K)=TORC(M,K)+TI*ORC(M,K,L,NY,NX)
          TORN(M,K)=TORN(M,K)+TI*ORN(M,K,L,NY,NX)
          TORP(M,K)=TORP(M,K)+TI*ORP(M,K,L,NY,NX)
        ENDDO
        TOQC(K)=TOQC(K)+TI*OQC(K,L,NY,NX)
        TOQN(K)=TOQN(K)+TI*OQN(K,L,NY,NX)
        TOQP(K)=TOQP(K)+TI*OQP(K,L,NY,NX)
        TOQA(K)=TOQA(K)+TI*OQA(K,L,NY,NX)
        TOHC(K)=TOHC(K)+TI*OHC(K,L,NY,NX)
        TOHN(K)=TOHN(K)+TI*OHN(K,L,NY,NX)
        TOHP(K)=TOHP(K)+TI*OHP(K,L,NY,NX)
        TOHA(K)=TOHA(K)+TI*OHA(K,L,NY,NX)
        DO  M=1,jsken
          TOSC(M,K)=TOSC(M,K)+TI*OSC(M,K,L,NY,NX)
          TOSA(M,K)=TOSA(M,K)+TI*OSA(M,K,L,NY,NX)
          TOSN(M,K)=TOSN(M,K)+TI*OSN(M,K,L,NY,NX)
          TOSP(M,K)=TOSP(M,K)+TI*OSP(M,K,L,NY,NX)
        ENDDO
      ENDDO
      ZNHUX0=AMAX1(ZNHUX0,ZNHU0(L,NY,NX))
      ZNHUXI=AMAX1(ZNHUXI,ZNHUI(L,NY,NX))
      ZNFNX0=AMAX1(ZNFNX0,ZNFN0(L,NY,NX))
      TZNFNI=TZNFNI+ZNFNI(L,NY,NX)
      LL=L
    ENDIF
  ENDDO
!
!     CHANGE SOIL STATE VARIABLES IN TILLAGE MIXING ZONE
!     TO ACCOUNT FOR REDISTRIBUTION FROM MIXING
!   but non-P related precipitated species are not mixed, Jinyun Tang, Nov 30,2022

    D2000: DO  L=NU(NY,NX),LL
      IF(DLYR(3,L,NY,NX).GT.ZERO)THEN
        TL=AMIN1(DLYR(3,L,NY,NX),DCORPZ-(CDPTHZ(L,NY,NX)-DLYR(3,L,NY,NX)))
        FI=TL/DCORPZ
        TI=TL/DLYR(3,L,NY,NX)
        TX=1.0_r8-TI
        BKDSI(L,NY,NX)=TI*(BKDSI(L,NY,NX)+CORP*(TBKDX-BKDSI(L,NY,NX)))+TX*BKDSI(L,NY,NX)
        FC(L,NY,NX)=TI*(FC(L,NY,NX)+CORP*(TFC-FC(L,NY,NX)))+TX*FC(L,NY,NX)
        WP(L,NY,NX)=TI*(WP(L,NY,NX)+CORP*(TWP-WP(L,NY,NX)))+TX*WP(L,NY,NX)
        SCNV(L,NY,NX)=TI*(SCNV(L,NY,NX)+CORP*(TSCNV-SCNV(L,NY,NX)))+TX*SCNV(L,NY,NX)
        SCNH(L,NY,NX)=TI*(SCNH(L,NY,NX)+CORP*(TSCNH-SCNH(L,NY,NX)))+TX*SCNH(L,NY,NX)
        SAND(L,NY,NX)=TI*SAND(L,NY,NX)+CORP*(FI*TSAND-TI*SAND(L,NY,NX))+TX*SAND(L,NY,NX)
        SILT(L,NY,NX)=TI*SILT(L,NY,NX)+CORP*(FI*TSILT-TI*SILT(L,NY,NX))+TX*SILT(L,NY,NX)
        CLAY(L,NY,NX)=TI*CLAY(L,NY,NX)+CORP*(FI*TCLAY-TI*CLAY(L,NY,NX))+TX*CLAY(L,NY,NX)

        GKC4(L,NY,NX)=TI*(GKC4(L,NY,NX)+CORP*(TGKC4-GKC4(L,NY,NX)))+TX*GKC4(L,NY,NX)
        GKCA(L,NY,NX)=TI*(GKCA(L,NY,NX)+CORP*(TGKCA-GKCA(L,NY,NX)))+TX*GKCA(L,NY,NX)
        GKCM(L,NY,NX)=TI*(GKCM(L,NY,NX)+CORP*(TGKCM-GKCM(L,NY,NX)))+TX*GKCM(L,NY,NX)
        GKCN(L,NY,NX)=TI*(GKCN(L,NY,NX)+CORP*(TGKCN-GKCN(L,NY,NX)))+TX*GKCN(L,NY,NX)
        GKCK(L,NY,NX)=TI*(GKCK(L,NY,NX)+CORP*(TGKCK-GKCK(L,NY,NX)))+TX*GKCK(L,NY,NX)
        ENGYM=VHCM(L,NY,NX)*TKS(L,NY,NX)
        ENGYV=(cpw*(VOLW(L,NY,NX)+VOLWH(L,NY,NX))+cpi*(VOLI(L,NY,NX)+VOLIH(L,NY,NX)))*TKS(L,NY,NX)
        VOLW(L,NY,NX)=TI*VOLW(L,NY,NX)+CORP*(FI*TVOLW-TI*VOLW(L,NY,NX))+TX*VOLW(L,NY,NX)+FI*TVOLWR
        VOLI(L,NY,NX)=TI*VOLI(L,NY,NX)+CORP*(FI*TVOLI-TI*VOLI(L,NY,NX))+TX*VOLI(L,NY,NX)
      VOLWX(L,NY,NX)=VOLW(L,NY,NX)
!     VOLW(L,NY,NX)=VOLW(L,NY,NX)+CORP*VOLWH(L,NY,NX)
!     VOLI(L,NY,NX)=VOLI(L,NY,NX)+CORP*VOLIH(L,NY,NX)
!     VOLA(L,NY,NX)=VOLA(L,NY,NX)+CORP*VOLAH(L,NY,NX)
!     VOLWH(L,NY,NX)=XCORP(NY,NX)*VOLWH(L,NY,NX)
!     VOLIH(L,NY,NX)=XCORP(NY,NX)*VOLIH(L,NY,NX)
!     VOLAH(L,NY,NX)=XCORP(NY,NX)*VOLAH(L,NY,NX)
!     FHOL(L,NY,NX)=XCORP(NY,NX)*FHOL(L,NY,NX)
        ENGYL=TI*ENGYV+CORP*(FI*TENGY-TI*ENGYV)+TX*ENGYV+FI*TENGYR
        VHCP(L,NY,NX)=VHCM(L,NY,NX)+cpw*(VOLW(L,NY,NX)+VOLWH(L,NY,NX)) &
          +cpi*(VOLI(L,NY,NX)+VOLIH(L,NY,NX))
        TKS(L,NY,NX)=(ENGYM+ENGYL)/VHCP(L,NY,NX)
        TCS(L,NY,NX)=TKS(L,NY,NX)-TC2K
        DO NTN=ifertn_beg,ifertn_end
          FertN_soil(NTN,L,NY,NX)=TI*FertN_soil(NTN,L,NY,NX) &
            +CORP*(FI*TfertN_soil(NTN)-TI*FertN_soil(NTN,L,NY,NX))&
            +TX*FertN_soil(NTN,L,NY,NX)
        ENDDO

        DO NTN=ifertnb_beg,ifertnb_end
          FertN_band(NTN,L,NY,NX)=TI*FertN_band(NTN,L,NY,NX) &
            +CORP*(FI*TfertN_band(NTN)-TI*FertN_band(NTN,L,NY,NX)) &
            +TX*FertN_band(NTN,L,NY,NX)
        ENDDO

        !SALT
        DO NTSA=idsa_beg,idsa_end
          trcsa_solml(NTSA,L,NY,NX)=TI*trcsa_solml(NTSA,L,NY,NX)  &
            +CORP*(FI*TSA_solml(NTSA)-TI*trcsa_solml(NTSA,L,NY,NX)) &
            +TX*trcsa_solml(NTSA,L,NY,NX)+CORP*trcsa_soHml(NTSA,L,NY,NX)
        ENDDO

        ! solute
        DO NTS=ids_beg,ids_end
          trc_solml(NTS,L,NY,NX)=TI*trc_solml(NTS,L,NY,NX) &
            +CORP*(FI*TS_solml(NTS)-TI*trc_solml(NTS,L,NY,NX)) &
            +TX*trc_solml(NTS,L,NY,NX)+CORP*trc_soHml(NTS,L,NY,NX)
        ENDDO

        DO NTX=idx_CEC,idx_cation_end
          trcx_solml(NTX,L,NY,NX)=TI*trcx_solml(NTX,L,NY,NX)+ &
            CORP*(FI*TX_solml(NTX)-TI*trcx_solml(NTX,L,NY,NX))+TX*trcx_solml(NTX,L,NY,NX)
        ENDDO

        DO NTX=idx_AEC,idx_end
          trcx_solml(NTX,L,NY,NX)=TI*trcx_solml(NTX,L,NY,NX)+ &
            CORP*(FI*TX_solml(NTX)-TI*trcx_solml(NTX,L,NY,NX))+TX*trcx_solml(NTX,L,NY,NX)
        ENDDO

        DO NTP=idsp_beg,idsp_end
          trcp_salml(NTP,L,NY,NX)=TI*trcp_salml(NTP,L,NY,NX)+CORP*(FI*TP_salml(NTP) &
            -TI*trcp_salml(NTP,L,NY,NX))+TX*trcp_salml(NTP,L,NY,NX)
        ENDDO

        DO NTG=idg_beg,idg_end-1
          trc_gasml(NTG,L,NY,NX)=TI*trc_gasml(NTG,L,NY,NX)+CORP*(FI*TG_gasml(NTG) &
            -TI*trc_gasml(NTG,L,NY,NX))+TX*trc_gasml(NTG,L,NY,NX)
        ENDDO

        call MixSoluteH(L,NY,NX)

        DO  K=1,jcplx
          DO  N=1,NFGs
            DO NGL=JGnio(N),JGnfo(N)
              DO  M=1,nlbiomcp
                OMC(M,NGL,K,L,NY,NX)=TI*OMC(M,NGL,K,L,NY,NX)+CORP*(FI*TOMC(M,NGL,K) &
                  -TI*OMC(M,NGL,K,L,NY,NX))+TX*OMC(M,NGL,K,L,NY,NX)
                OMN(M,NGL,K,L,NY,NX)=TI*OMN(M,NGL,K,L,NY,NX)+CORP*(FI*TOMN(M,NGL,K) &
                  -TI*OMN(M,NGL,K,L,NY,NX))+TX*OMN(M,NGL,K,L,NY,NX)
                OMP(M,NGL,K,L,NY,NX)=TI*OMP(M,NGL,K,L,NY,NX)+CORP*(FI*TOMP(M,NGL,K) &
                  -TI*OMP(M,NGL,K,L,NY,NX))+TX*OMP(M,NGL,K,L,NY,NX)
              enddo
            enddo
          enddo
        ENDDO
        DO  N=1,NFGs
          DO NGL=JGniA(N),JGnfA(N)
            DO  M=1,nlbiomcp
              OMCff(M,NGL,L,NY,NX)=TI*OMCff(M,NGL,L,NY,NX)+CORP*(FI*TOMCff(M,NGL) &
                -TI*OMCff(M,NGL,L,NY,NX))+TX*OMCff(M,NGL,L,NY,NX)
              OMNff(M,NGL,L,NY,NX)=TI*OMNff(M,NGL,L,NY,NX)+CORP*(FI*TOMNff(M,NGL) &
                -TI*OMNff(M,NGL,L,NY,NX))+TX*OMNff(M,NGL,L,NY,NX)
              OMPff(M,NGL,L,NY,NX)=TI*OMPff(M,NGL,L,NY,NX)+CORP*(FI*TOMPff(M,NGL) &
                -TI*OMPff(M,NGL,L,NY,NX))+TX*OMPff(M,NGL,L,NY,NX)
            enddo
          enddo
        enddo

        DO  K=1,jcplx
          DO  M=1,ndbiomcp
            ORC(M,K,L,NY,NX)=TI*ORC(M,K,L,NY,NX)+CORP*(FI*TORC(M,K) &
              -TI*ORC(M,K,L,NY,NX))+TX*ORC(M,K,L,NY,NX)
            ORN(M,K,L,NY,NX)=TI*ORN(M,K,L,NY,NX)+CORP*(FI*TORN(M,K) &
              -TI*ORN(M,K,L,NY,NX))+TX*ORN(M,K,L,NY,NX)
            ORP(M,K,L,NY,NX)=TI*ORP(M,K,L,NY,NX)+CORP*(FI*TORP(M,K) &
              -TI*ORP(M,K,L,NY,NX))+TX*ORP(M,K,L,NY,NX)
          ENDDO
          OQC(K,L,NY,NX)=TI*OQC(K,L,NY,NX)+CORP*(FI*TOQC(K) &
            -TI*OQC(K,L,NY,NX))+TX*OQC(K,L,NY,NX)+CORP*OQCH(K,L,NY,NX)
          OQN(K,L,NY,NX)=TI*OQN(K,L,NY,NX)+CORP*(FI*TOQN(K) &
            -TI*OQN(K,L,NY,NX))+TX*OQN(K,L,NY,NX)+CORP*OQNH(K,L,NY,NX)
          OQP(K,L,NY,NX)=TI*OQP(K,L,NY,NX)+CORP*(FI*TOQP(K) &
            -TI*OQP(K,L,NY,NX))+TX*OQP(K,L,NY,NX)+CORP*OQPH(K,L,NY,NX)
          OQA(K,L,NY,NX)=TI*OQA(K,L,NY,NX)+CORP*(FI*TOQA(K) &
            -TI*OQA(K,L,NY,NX))+TX*OQA(K,L,NY,NX)+CORP*OQAH(K,L,NY,NX)

          OQCH(K,L,NY,NX)=XCORP(NY,NX)*OQCH(K,L,NY,NX)
          OQNH(K,L,NY,NX)=XCORP(NY,NX)*OQNH(K,L,NY,NX)
          OQPH(K,L,NY,NX)=XCORP(NY,NX)*OQPH(K,L,NY,NX)
          OQAH(K,L,NY,NX)=XCORP(NY,NX)*OQAH(K,L,NY,NX)

          OHC(K,L,NY,NX)=TI*OHC(K,L,NY,NX)+CORP*(FI*TOHC(K) &
            -TI*OHC(K,L,NY,NX))+TX*OHC(K,L,NY,NX)
          OHN(K,L,NY,NX)=TI*OHN(K,L,NY,NX)+CORP*(FI*TOHN(K) &
            -TI*OHN(K,L,NY,NX))+TX*OHN(K,L,NY,NX)
          OHP(K,L,NY,NX)=TI*OHP(K,L,NY,NX)+CORP*(FI*TOHP(K) &
            -TI*OHP(K,L,NY,NX))+TX*OHP(K,L,NY,NX)
          OHA(K,L,NY,NX)=TI*OHA(K,L,NY,NX)+CORP*(FI*TOHA(K) &
            -TI*OHA(K,L,NY,NX))+TX*OHA(K,L,NY,NX)
          DO  M=1,jsken
            OSC(M,K,L,NY,NX)=TI*OSC(M,K,L,NY,NX)+CORP*(FI*TOSC(M,K) &
              -TI*OSC(M,K,L,NY,NX))+TX*OSC(M,K,L,NY,NX)
            OSA(M,K,L,NY,NX)=TI*OSA(M,K,L,NY,NX)+CORP*(FI*TOSA(M,K) &
              -TI*OSA(M,K,L,NY,NX))+TX*OSA(M,K,L,NY,NX)
            OSN(M,K,L,NY,NX)=TI*OSN(M,K,L,NY,NX)+CORP*(FI*TOSN(M,K) &
              -TI*OSN(M,K,L,NY,NX))+TX*OSN(M,K,L,NY,NX)
            OSP(M,K,L,NY,NX)=TI*OSP(M,K,L,NY,NX)+CORP*(FI*TOSP(M,K) &
              -TI*OSP(M,K,L,NY,NX))+TX*OSP(M,K,L,NY,NX)
          ENDDO
        ENDDO
!
!     ADD STATE VARIABLES IN SURFACE RESIDUE INCORPORATED
!     WITHIN TILLAGE MIXING ZONE
!
        DO  K=1,micpar%n_litrsfk
          DO  N=1,NFGs
            DO NGL=JGnio(N),JGnfo(N)
              DO M=1,nlbiomcp
                OMC(M,NGL,K,L,NY,NX)=OMC(M,NGL,K,L,NY,NX)+FI*TOMGC(M,NGL,K)
                OMN(M,NGL,K,L,NY,NX)=OMN(M,NGL,K,L,NY,NX)+FI*TOMGN(M,NGL,K)
                OMP(M,NGL,K,L,NY,NX)=OMP(M,NGL,K,L,NY,NX)+FI*TOMGP(M,NGL,K)
              enddo
            enddo
          ENDDO
        ENDDO

        DO  N=1,NFGs
          DO NGL=JGniA(N),JGnfA(N)
            DO M=1,nlbiomcp
              OMCff(M,NGL,L,NY,NX)=OMCff(M,NGL,L,NY,NX)+FI*TOMGCff(M,NGL)
              OMNff(M,NGL,L,NY,NX)=OMNff(M,NGL,L,NY,NX)+FI*TOMGNff(M,NGL)
              OMPff(M,NGL,L,NY,NX)=OMPff(M,NGL,L,NY,NX)+FI*TOMGPff(M,NGL)
            enddo
          enddo
        ENDDO

        DO K=1,micpar%n_litrsfk
          DO  M=1,ndbiomcp
            ORC(M,K,L,NY,NX)=ORC(M,K,L,NY,NX)+FI*TORXC(M,K)
            ORN(M,K,L,NY,NX)=ORN(M,K,L,NY,NX)+FI*TORXN(M,K)
            ORP(M,K,L,NY,NX)=ORP(M,K,L,NY,NX)+FI*TORXP(M,K)
          ENDDO
          OQC(K,L,NY,NX)=OQC(K,L,NY,NX)+FI*TOQGC(K)
          OQN(K,L,NY,NX)=OQN(K,L,NY,NX)+FI*TOQGN(K)
          OQP(K,L,NY,NX)=OQP(K,L,NY,NX)+FI*TOQGP(K)
          OQA(K,L,NY,NX)=OQA(K,L,NY,NX)+FI*TOQGA(K)
          OQCH(K,L,NY,NX)=OQCH(K,L,NY,NX)+FI*TOQHC(K)
          OQNH(K,L,NY,NX)=OQNH(K,L,NY,NX)+FI*TOQHN(K)
          OQPH(K,L,NY,NX)=OQPH(K,L,NY,NX)+FI*TOQHP(K)
          OQAH(K,L,NY,NX)=OQAH(K,L,NY,NX)+FI*TOQHA(K)
          OHC(K,L,NY,NX)=OHC(K,L,NY,NX)+FI*TOHGC(K)
          OHN(K,L,NY,NX)=OHN(K,L,NY,NX)+FI*TOHGN(K)
          OHP(K,L,NY,NX)=OHP(K,L,NY,NX)+FI*TOHGP(K)
          OHA(K,L,NY,NX)=OHA(K,L,NY,NX)+FI*TOHGA(K)
          DO  M=1,jsken
            OSC(M,K,L,NY,NX)=OSC(M,K,L,NY,NX)+FI*TOSGC(M,K)
            OSA(M,K,L,NY,NX)=OSA(M,K,L,NY,NX)+FI*TOSGA(M,K)
            OSN(M,K,L,NY,NX)=OSN(M,K,L,NY,NX)+FI*TOSGN(M,K)
            OSP(M,K,L,NY,NX)=OSP(M,K,L,NY,NX)+FI*TOSGP(M,K)
          ENDDO
        ENDDO
        OC=0.0_r8
        ON=0.0_r8
        OP=0.0_r8
        DC=0.0_r8
        DN=0.0_r8
        DP=0.0_r8

        DO  K=1,jcplx
          DO  N=1,NFGs
            DO NGL=JGnio(N),JGnfo(N)
              DO  M=1,nlbiomcp
                OC=OC+OMC(M,NGL,K,L,NY,NX)
                ON=ON+OMN(M,NGL,K,L,NY,NX)
                OP=OP+OMP(M,NGL,K,L,NY,NX)
              enddo
            enddo
          enddo
        ENDDO

        DO  K=1,micpar%n_litrsfk
          DO  N=1,NFGs
            DO NGL=JGnio(N),JGnfo(N)
              DO  M=1,nlbiomcp
                DC=DC+OMC(M,NGL,K,L,NY,NX)
                DN=DN+OMN(M,NGL,K,L,NY,NX)
                DP=DP+OMP(M,NGL,K,L,NY,NX)
              enddo
            enddo
          enddo
        ENDDO
        DO  N=1,NFGs
          DO NGL=JGniA(N),JGnfA(N)
            DO  M=1,nlbiomcp
              OC=OC+OMCff(M,NGL,L,NY,NX)
              ON=ON+OMNff(M,NGL,L,NY,NX)
              OP=OP+OMPff(M,NGL,L,NY,NX)
            enddo
          enddo
        enddo

        DO K=1,jcplx
          DO  M=1,ndbiomcp
            OC=OC+ORC(M,K,L,NY,NX)
            ON=ON+ORN(M,K,L,NY,NX)
            OP=OP+ORP(M,K,L,NY,NX)
            IF(micpar%is_litter(K))THEN
              DC=DC+ORC(M,K,L,NY,NX)
              DN=DN+ORN(M,K,L,NY,NX)
              DP=DP+ORP(M,K,L,NY,NX)
            ENDIF
          ENDDO

          OC=OC+OQC(K,L,NY,NX)+OQCH(K,L,NY,NX)+OHC(K,L,NY,NX) &
            +OQA(K,L,NY,NX)+OQAH(K,L,NY,NX)+OHA(K,L,NY,NX)
          ON=ON+OQN(K,L,NY,NX)+OQNH(K,L,NY,NX)+OHN(K,L,NY,NX)
          OP=OP+OQP(K,L,NY,NX)+OQPH(K,L,NY,NX)+OHP(K,L,NY,NX)
          IF(micpar%is_litter(K))THEN
            DC=DC+OQC(K,L,NY,NX)+OQCH(K,L,NY,NX)+OHC(K,L,NY,NX) &
                +OQA(K,L,NY,NX)+OQAH(K,L,NY,NX)+OHA(K,L,NY,NX)
            DN=DN+OQN(K,L,NY,NX)+OQNH(K,L,NY,NX)+OHN(K,L,NY,NX)
            DC=DC+OQC(K,L,NY,NX)+OQCH(K,L,NY,NX)+OHC(K,L,NY,NX)
          ENDIF
          DO  M=1,jsken
            OC=OC+OSC(M,K,L,NY,NX)
            ON=ON+OSN(M,K,L,NY,NX)
            OP=OP+OSP(M,K,L,NY,NX)
            IF(micpar%is_litter(K))THEN
              DC=DC+OSC(M,K,L,NY,NX)
              DN=DN+OSN(M,K,L,NY,NX)
              DP=DP+OSP(M,K,L,NY,NX)
            ENDIF
          ENDDO
        ENDDO
        ORGC(L,NY,NX)=OC
        ORGN(L,NY,NX)=ON
        ORGR(L,NY,NX)=DC

        DO NTS=ids_beg,idg_NH3
          trc_solml(NTS,L,NY,NX)=trc_solml(NTS,L,NY,NX)+FI*TS0_solml(NTS)
        ENDDO

        DO NTS=ids_nut_beg,ids_nuts_end
          trc_solml(NTS,L,NY,NX)=trc_solml(NTS,L,NY,NX)+FI*TS0_solml(NTS)
        ENDDO

        trcx_solml(idx_NH4,L,NY,NX)=trcx_solml(idx_NH4,L,NY,NX)+FI*TXN4G


        DO NTX=idx_AEC+1,idx_anion_soil_end
          trcx_solml(NTX,L,NY,NX)=trcx_solml(NTX,L,NY,NX)+FI*TX_anion(NTX)
        ENDDO

        DO NTP=idsp_psoi_beg,idsp_psoi_end
          trcp_salml(NTP,L,NY,NX)=trcp_salml(NTP,L,NY,NX)+FI*TP_soil(NTP)
        ENDDO

        DO NTN=ifertn_beg,ifertn_end
          FertN_soil(NTN,L,NY,NX)=FertN_soil(NTN,L,NY,NX)+FI*TFertNG_soil(NTN)
        ENDDO

        ZNHU0(L,NY,NX)=ZNHUX0
        ZNHUI(L,NY,NX)=ZNHUXI
        ZNFN0(L,NY,NX)=ZNFNX0
        ZNFNI(L,NY,NX)=(TI*ZNFNI(L,NY,NX)+CORP*(FI*TZNFNI &
          -TI*ZNFNI(L,NY,NX))+TX*ZNFNI(L,NY,NX)+FI*TZNFNG)/FI
        TZNFN2=TZNFN2+ZNFNI(L,NY,NX)
      ENDIF
    ENDDO D2000

! nitrogen inhibitor
    ZNFN0(0,NY,NX)=ZNFNX0
    ZNFNI(0,NY,NX)=ZNFNI(0,NY,NX)*XCORP0
    TZNFN2=TZNFN2+TZNFNG
    TZNFNI=TZNFNI+TZNFNG
    DO  L=NU(NY,NX),LL
      IF(TZNFN2.GT.ZERO)THEN
      ZNFNI(L,NY,NX)=ZNFNI(L,NY,NX)*TZNFNI/TZNFN2
      ZNFNI(L,NY,NX)=ZNFNI(L,NY,NX)+0.5_r8*(ZNFN0(L,NY,NX)-ZNFNI(L,NY,NX))
      ENDIF
    ENDDO
  ENDIF
  end subroutine ApplyTillageMixing


!--------------------------------------------------------------------------------

  subroutine MixSoluteH(L,NY,NX)
  implicit none
  integer, intent(in) :: L,NY,NX
  integer :: NTS,NTSA

  DO NTS=ids_beg,ids_end
    trc_soHml(NTS,L,NY,NX)=XCORP(NY,NX)*trc_soHml(NTS,L,NY,NX)
  ENDDO

  DO NTSA=idsa_beg,idsa_end
    trcsa_soHml(ntsa,L,NY,NX)=XCORP(NY,NX)*trcsa_soHml(ntsa,L,NY,NX)
  ENDDO
  end subroutine MixSoluteH


end module TillageMixMod
