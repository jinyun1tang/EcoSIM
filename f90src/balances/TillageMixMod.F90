module TillageMixMod
  use data_kind_mod, only : r8 => DAT_KIND_R8
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
  use EcoSIMConfig , only : ndbiomcp => NumDeadMicrbCompts
  use UnitMod, only : units
  implicit none
  character(len=*),private, parameter :: mod_filename =&
     __FILE__
  public :: ApplyTillageMixing
  contains

!------------------------------------------------------------------------------------------

  subroutine ApplyTillageMixing(I,J,NY,NX)
  implicit none
  integer, intent(in) :: I,J,NY,NX

  integer :: K,N,M,L,LL,NGL,NTS

  REAL(R8) :: TOSGC(jsken,0:micpar%NumOfLitrCmplxs),TOSGA(jsken,0:micpar%NumOfLitrCmplxs)
  real(r8) :: TOSGN(jsken,0:micpar%NumOfLitrCmplxs),TOSGP(jsken,0:micpar%NumOfLitrCmplxs)
  real(r8) :: TORXC(ndbiomcp,0:micpar%NumOfLitrCmplxs)
  real(r8) :: TORXN(ndbiomcp,0:micpar%NumOfLitrCmplxs),TORXP(ndbiomcp,0:micpar%NumOfLitrCmplxs)
  real(r8) :: TOQGC(0:micpar%NumOfLitrCmplxs),TOQGN(0:micpar%NumOfLitrCmplxs)
  real(r8) :: TOQGP(0:micpar%NumOfLitrCmplxs),TOQHC(0:micpar%NumOfLitrCmplxs)
  real(r8) :: TOQHN(0:micpar%NumOfLitrCmplxs),TOQHP(0:micpar%NumOfLitrCmplxs)
  real(r8) :: TOHGC(0:micpar%NumOfLitrCmplxs),TOHGN(0:micpar%NumOfLitrCmplxs)
  real(r8) :: TOHGP(0:micpar%NumOfLitrCmplxs),TOHGA(0:micpar%NumOfLitrCmplxs)
  real(r8) :: TOQGA(0:micpar%NumOfLitrCmplxs),TOQHA(0:micpar%NumOfLitrCmplxs)

  real(r8) :: TOMGCff(1:NumLiveAutoBioms)
  real(r8) :: TOMGNff(1:NumLiveAutoBioms)
  real(r8) :: TOMGPff(1:NumLiveAutoBioms)
  real(r8) :: TOMGC(1:NumLiveHeterBioms,1:jcplx)
  real(r8) :: TOMGN(1:NumLiveHeterBioms,1:jcplx)
  real(r8) :: TOMGP(1:NumLiveHeterBioms,1:jcplx)
  REAL(R8) :: TOMEheter(NumPlantChemElms,NumLiveHeterBioms,1:jcplx)
  REAL(R8) :: TOMEauto(NumPlantChemElms,NumLiveAutoBioms)
  real(r8) :: TORM(NumPlantChemElms,ndbiomcp,1:jcplx)
  real(r8) :: TDOM(idom_beg:idom_end,1:jcplx)
  real(r8) :: TOHM(idom_beg:idom_end,1:jcplx)
  real(r8) :: TOSM(1:NumPlantChemElms,jsken,1:jcplx)
  real(r8) :: TOSA(jsken,1:jcplx)
  real(r8) :: CORP,CORP0,XCORP0,DCORPZ
  real(r8) :: TBKDX,TFC,TWP,TSatHydroCondVert,TSCNH,TSAND
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
  real(r8) :: TSA_solml(idsalt_beg:idsalt_end)
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
  integer  :: NTX,NTP,NTG,NTSA,NTN,idom,MID,NE

!     begin_execution
!
  !DCORP=soil mixing fraction with tillage
  !XCORP=factor for surface litter incorporation and soil mixing
  IF(J.EQ.INT(SolarNoonHour_col(NY,NX)).AND.XCORP(NY,NX).LT.1.0.AND.DCORP(I,NY,NX).GT.0.0_r8)THEN
!
!     EXTENT OF MIXING
!
    IFLGS(NY,NX)=1
    CORP=1.0_r8-XCORP(NY,NX)
    EnergyImpact4Erosion(NY,NX)=0.0_r8
!
!     TEMPORARY ACCUMULATORS
!
    TBKDX=0.0_r8
    TFC=0.0_r8
    TWP=0.0_r8
    TSatHydroCondVert=0.0_r8
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

    TOMEheter=0.0_r8
    TOMEauto=0.0_r8
    TORM=0.0_r8
    TDOM=0.0_r8
    TOHM=0.0_r8
    TOSM=0.0_r8
    TOSA=0.0_r8
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
!    2(VHeatCapLitR(NY,NX)/cpo)/ORGC(0,NY,NX)))
    XCORP0=AMAX1(0.001_r8,XCORP(NY,NX))
!     ELSE
!     XCORP0=1.0
!     ENDIF
    CORP0=1.0_r8-XCORP0   !fraction mixed into next layer

    DC=0.0_r8
    DN=0.0_r8
    DP=0.0_r8

    DO  K=1,micpar%NumOfLitrCmplxs
      DO  N=1,NumMicbFunGroups
        DO NGL=JGnio(N),JGnfo(N)
          DO  M=1,nlbiomcp
            MID=micpar%get_micb_id(M,NGL)
            TOMGC(MID,K)=OMEheter(ielmc,MID,K,0,NY,NX)*CORP0
            TOMGN(MID,K)=OMEheter(ielmn,MID,K,0,NY,NX)*CORP0
            TOMGP(MID,K)=OMEheter(ielmp,MID,K,0,NY,NX)*CORP0
            OMEheter(ielmc,MID,K,0,NY,NX)=OMEheter(ielmc,MID,K,0,NY,NX)*XCORP0
            OMEheter(ielmn,MID,K,0,NY,NX)=OMEheter(ielmn,MID,K,0,NY,NX)*XCORP0
            OMEheter(ielmp,MID,K,0,NY,NX)=OMEheter(ielmp,MID,K,0,NY,NX)*XCORP0
            DC=DC+OMEheter(ielmc,MID,K,0,NY,NX)
            DN=DN+OMEheter(ielmn,MID,K,0,NY,NX)
            DP=DP+OMEheter(ielmp,MID,K,0,NY,NX)
          enddo
        enddo
      ENDDO
    ENDDO

    DO  N=1,NumMicbFunGroups
      DO NGL=JGniA(N),JGnfA(N)
        DO  M=1,nlbiomcp
          MID=micpar%get_micb_id(M,NGL)
          TOMGCff(MID)=OMEauto(ielmc,MID,0,NY,NX)*CORP0
          TOMGNff(MID)=OMEauto(ielmn,MID,0,NY,NX)*CORP0
          TOMGPff(MID)=OMEauto(ielmp,MID,0,NY,NX)*CORP0
          OMEauto(ielmc,MID,0,NY,NX)=OMEauto(ielmc,MID,0,NY,NX)*XCORP0
          OMEauto(ielmn,MID,0,NY,NX)=OMEauto(ielmn,MID,0,NY,NX)*XCORP0
          OMEauto(ielmp,MID,0,NY,NX)=OMEauto(ielmp,MID,0,NY,NX)*XCORP0
          DC=DC+OMEauto(ielmc,MID,0,NY,NX)
          DN=DN+OMEauto(ielmn,MID,0,NY,NX)
          DP=DP+OMEauto(ielmp,MID,0,NY,NX)
        enddo
      enddo
    ENDDO

    DO K=1,micpar%NumOfLitrCmplxs
      DO M=1,ndbiomcp
        TORXC(M,K)=ORM(ielmc,M,K,0,NY,NX)*CORP0
        TORXN(M,K)=ORM(ielmn,M,K,0,NY,NX)*CORP0
        TORXP(M,K)=ORM(ielmp,M,K,0,NY,NX)*CORP0
        ORM(ielmc,M,K,0,NY,NX)=ORM(ielmc,M,K,0,NY,NX)*XCORP0
        ORM(ielmn,M,K,0,NY,NX)=ORM(ielmn,M,K,0,NY,NX)*XCORP0
        ORM(ielmp,M,K,0,NY,NX)=ORM(ielmp,M,K,0,NY,NX)*XCORP0
        DC=DC+ORM(ielmc,M,K,0,NY,NX)
        DN=DN+ORM(ielmn,M,K,0,NY,NX)
        DP=DP+ORM(ielmp,M,K,0,NY,NX)
      ENDDO
      TOQGC(K)=DOM(idom_doc,K,0,NY,NX)*CORP0
      TOQGN(K)=DOM(idom_don,K,0,NY,NX)*CORP0
      TOQGP(K)=DOM(idom_dop,K,0,NY,NX)*CORP0
      TOQGA(K)=DOM(idom_acetate,K,0,NY,NX)*CORP0
      TOQHC(K)=DOM_Macp(idom_doc,K,0,NY,NX)*CORP0
      TOQHN(K)=DOM_Macp(idom_don,K,0,NY,NX)*CORP0
      TOQHP(K)=DOM_Macp(idom_dop,K,0,NY,NX)*CORP0
      TOQHA(K)=DOM_Macp(idom_acetate,K,0,NY,NX)*CORP0
      TOHGC(K)=OHM(ielmc,K,0,NY,NX)*CORP0
      TOHGN(K)=OHM(ielmn,K,0,NY,NX)*CORP0
      TOHGP(K)=OHM(ielmp,K,0,NY,NX)*CORP0
      TOHGA(K)=OHM(idom_acetate,K,0,NY,NX)*CORP0
  !
  !     REDUCE SURFACE RESIDUE STATE VARIABLES FOR INCORPORATION
  !
      DOM(idom_doc,K,0,NY,NX)=DOM(idom_doc,K,0,NY,NX)*XCORP0
      DOM(idom_don,K,0,NY,NX)=DOM(idom_don,K,0,NY,NX)*XCORP0
      DOM(idom_dop,K,0,NY,NX)=DOM(idom_dop,K,0,NY,NX)*XCORP0
      DOM(idom_acetate,K,0,NY,NX)=DOM(idom_acetate,K,0,NY,NX)*XCORP0
      DOM_Macp(idom_doc,K,0,NY,NX)=DOM_Macp(idom_doc,K,0,NY,NX)*XCORP0
      DOM_Macp(idom_don,K,0,NY,NX)=DOM_Macp(idom_don,K,0,NY,NX)*XCORP0
      DOM_Macp(idom_dop,K,0,NY,NX)=DOM_Macp(idom_dop,K,0,NY,NX)*XCORP0
      DOM_Macp(idom_acetate,K,0,NY,NX)=DOM_Macp(idom_acetate,K,0,NY,NX)*XCORP0
      OHM(ielmc,K,0,NY,NX)=OHM(ielmc,K,0,NY,NX)*XCORP0
      OHM(ielmn,K,0,NY,NX)=OHM(ielmn,K,0,NY,NX)*XCORP0
      OHM(ielmp,K,0,NY,NX)=OHM(ielmp,K,0,NY,NX)*XCORP0
      OHM(idom_acetate,K,0,NY,NX)=OHM(idom_acetate,K,0,NY,NX)*XCORP0
      DC=DC+DOM(idom_doc,K,0,NY,NX)+DOM_Macp(idom_doc,K,0,NY,NX)+OHM(ielmc,K,0,NY,NX) &
        +DOM(idom_acetate,K,0,NY,NX)+DOM_Macp(idom_acetate,K,0,NY,NX)+OHM(idom_acetate,K,0,NY,NX)
      DN=DN+DOM(idom_don,K,0,NY,NX)+DOM_Macp(idom_don,K,0,NY,NX)+OHM(ielmn,K,0,NY,NX)
      DP=DP+DOM(idom_dop,K,0,NY,NX)+DOM_Macp(idom_dop,K,0,NY,NX)+OHM(ielmp,K,0,NY,NX)
      DO  M=1,jsken
        TOSGC(M,K)=OSM(ielmc,M,K,0,NY,NX)*CORP0
        TOSGA(M,K)=OSA(M,K,0,NY,NX)*CORP0
        TOSGN(M,K)=OSM(ielmn,M,K,0,NY,NX)*CORP0
        TOSGP(M,K)=OSM(ielmp,M,K,0,NY,NX)*CORP0
        OSM(ielmc,M,K,0,NY,NX)=OSM(ielmc,M,K,0,NY,NX)*XCORP0
        OSA(M,K,0,NY,NX)=OSA(M,K,0,NY,NX)*XCORP0
        OSM(ielmn,M,K,0,NY,NX)=OSM(ielmn,M,K,0,NY,NX)*XCORP0
        OSM(ielmp,M,K,0,NY,NX)=OSM(ielmp,M,K,0,NY,NX)*XCORP0
        DC=DC+OSM(ielmc,M,K,0,NY,NX)
        DN=DN+OSM(ielmn,M,K,0,NY,NX)
        DP=DP+OSM(ielmp,M,K,0,NY,NX)
      ENDDO
    ENDDO

    DO NTS=ids_beg,idg_NH3
      TS0_solml(NTS)=trc_solml_vr(NTS,0,NY,NX)*CORP0
    ENDDO
    DO NTS=ids_nut_beg,ids_nuts_end
      TS0_solml(NTS)=trc_solml_vr(NTS,0,NY,NX)*CORP0
    ENDDO

    TXN4G=trcx_solml(idx_NH4,0,NY,NX)*CORP0

    DO NTX= idx_AEC+1,idx_anion_soil_end
      TX_anion(NTX)=trcx_solml(NTX,0,NY,NX)*CORP0
    ENDDO

    DO NTP=idsp_psoi_beg,idsp_psoi_end
      TP_soil(NTP)=trcp_salml(NTP,0,NY,NX)*CORP0
    ENDDO

    DO NTN=ifertn_beg,ifertn_end
      TFertNG_soil(NTN)=FertN_soil_vr(NTN,0,NY,NX)*CORP0
    ENDDO

    TZNFNG=ZNFNI(0,NY,NX)*CORP0
    TVOLWR=VLWatMicP(0,NY,NX)*CORP0
    HFLXD=cpo*ORGC(0,NY,NX)*CORP0*TKS(0,NY,NX)
    HEATIN=HEATIN-HFLXD
    HeatStoreLandscape=HeatStoreLandscape-HFLXD
    TENGYR=cpw*TVOLWR*TKS(0,NY,NX)
    ORGC(0,NY,NX)=DC
    ORGN(0,NY,NX)=DN
    ORGR(0,NY,NX)=DC

    DO NTS=ids_beg,idg_NH3
      trc_solml_vr(NTS,0,NY,NX)=trc_solml_vr(NTS,0,NY,NX)*XCORP0
    ENDDO

    DO NTS=ids_nut_beg,ids_nuts_end
      trc_solml_vr(NTS,0,NY,NX)=trc_solml_vr(NTS,0,NY,NX)*XCORP0
    ENDDO

    trcx_solml(idx_NH4,0,NY,NX)=trcx_solml(idx_NH4,0,NY,NX)*XCORP0

    DO NTX=idx_AEC+1,idx_anion_soil_end
      trcx_solml(idx_OHe,0,NY,NX)=trcx_solml(idx_OHe,0,NY,NX)*XCORP0
    ENDDO

    trcp_salml(idsp_AlPO4,0,NY,NX)=trcp_salml(idsp_AlPO4,0,NY,NX)*XCORP0
    trcp_salml(idsp_FePO4,0,NY,NX)=trcp_salml(idsp_FePO4,0,NY,NX)*XCORP0
    trcp_salml(idsp_CaHPO4,0,NY,NX)=trcp_salml(idsp_CaHPO4,0,NY,NX)*XCORP0
    trcp_salml(idsp_HA,0,NY,NX)=trcp_salml(idsp_HA,0,NY,NX)*XCORP0
    trcp_salml(idsp_CaH4P2O8,0,NY,NX)=trcp_salml(idsp_CaH4P2O8,0,NY,NX)*XCORP0

    DO NTN=ifertn_beg,ifertn_end
      FertN_soil_vr(NTN,0,NY,NX)=FertN_soil_vr(NTN,0,NY,NX)*XCORP0
    ENDDO

    VLWatMicP(0,NY,NX)=VLWatMicP(0,NY,NX)*XCORP0
    VHeatCapacity(0,NY,NX)=cpo*ORGC(0,NY,NX)+cpw*VLWatMicP(0,NY,NX)+cpi*VLiceMicP(0,NY,NX)
    VLitR(NY,NX)=VLitR(NY,NX)*XCORP0
    VGeomLayer(0,NY,NX)=VGeomLayer(0,NY,NX)*XCORP0
    ZNHUX0=AMAX1(ZNHUX0,ZNHU0(0,NY,NX))
    ZNHUXI=AMAX1(ZNHUXI,ZNHUI(0,NY,NX))
    ZNFNX0=AMAX1(ZNFNX0,ZNFN0(0,NY,NX))
    LL=NU(NY,NX)
!
!     REDISTRIBUTE SOIL STATE VARIABLES DURING TILLAGE
!
    DCORPZ=AMIN1(DCORP(I,NY,NX),CumSoilThickness(NL(NY,NX),NY,NX))
!
!     ACCUMULATE SOIL STATE VARIABLES WITHIN TILLAGE MIXING ZONE
!
    DO L=NU(NY,NX),NL(NY,NX)
      IF(CumSoilThickness(L,NY,NX)-DLYR(3,L,NY,NX).LT.DCORPZ.AND.DLYR(3,L,NY,NX).GT.ZERO)THEN
        TL=AMIN1(DLYR(3,L,NY,NX),DCORPZ-(CumSoilThickness(L,NY,NX)-DLYR(3,L,NY,NX)))
        FI=TL/DCORPZ
        TI=TL/DLYR(3,L,NY,NX)
        TBKDX=TBKDX+FI*SoiBulkDensityt0(L,NY,NX)
        TFC=TFC+FI*FieldCapacity(L,NY,NX)
        TWP=TWP+FI*WiltPoint(L,NY,NX)
        TSatHydroCondVert=TSatHydroCondVert+FI*SatHydroCondVert(L,NY,NX)
        TSCNH=TSCNH+FI*SatHydroCondHrzn(L,NY,NX)
        TSAND=TSAND+TI*SAND(L,NY,NX)
        TSILT=TSILT+TI*SILT(L,NY,NX)
        TCLAY=TCLAY+TI*CLAY(L,NY,NX)
        TGKC4=TGKC4+FI*GKC4(L,NY,NX)
        TGKCA=TGKCA+FI*GKCA(L,NY,NX)
        TGKCM=TGKCM+FI*GKCM(L,NY,NX)
        TGKCN=TGKCN+FI*GKCN(L,NY,NX)
        TGKCK=TGKCK+FI*GKCK(L,NY,NX)
        TVOLW=TVOLW+TI*VLWatMicP(L,NY,NX)
        TVOLI=TVOLI+TI*VLiceMicP(L,NY,NX)
!     TVOLP=TVOLP+TI*VLsoiAirP(L,NY,NX)
!     TVOLA=TVOLA+TI*VLMicP(L,NY,NX)
        TENGY=TENGY+TI*(cpw*(VLWatMicP(L,NY,NX)+VLWatMacP(L,NY,NX)) &
          +cpi*(VLiceMicP(L,NY,NX)+VLiceMacP(L,NY,NX)))*TKS(L,NY,NX)
        DO NTN=ifertn_beg,ifertn_end
          TfertN_soil(NTN)=TfertN_soil(NTN)+TI*FertN_soil_vr(NTN,L,NY,NX)
        ENDDO

        DO NTN=ifertnb_beg,ifertnb_end
          TfertN_band(NTN)=TfertN_band(NTN)+TI*FertN_Band_vr(NTN,L,NY,NX)
        ENDDO

        DO NTS=ids_beg,ids_end
          TS_solml(NTS)=TS_solml(NTS)+TI*trc_solml_vr(NTS,L,NY,NX)
        ENDDO

        DO NTSA=idsalt_beg,idsalt_end
          TSA_solml(NTSA)=TSA_solml(NTSA)+TI*trcSalt_solml(NTSA,L,NY,NX)
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
          TG_gasml(NTG)=TG_gasml(NTG)+TI*trc_gasml_vr(NTG,L,NY,NX)
        ENDDO

        DO  K=1,jcplx
          DO  N=1,NumMicbFunGroups
            DO NGL=JGnio(N),JGnfo(N)
              DO  M=1,nlbiomcp
                MID=micpar%get_micb_id(M,NGL)
                TOMEheter(ielmc,MID,K)=TOMEheter(ielmc,MID,K)+TI*OMEheter(ielmc,MID,K,L,NY,NX)
                TOMEheter(ielmn,MID,K)=TOMEheter(ielmn,MID,K)+TI*OMEheter(ielmn,MID,K,L,NY,NX)
                TOMEheter(ielmp,MID,K)=TOMEheter(ielmp,MID,K)+TI*OMEheter(ielmp,MID,K,L,NY,NX)
              enddo
            enddo
          enddo
        ENDDO

        DO  N=1,NumMicbFunGroups
          DO NGL=JGniA(N),JGnfA(N)
            DO  M=1,nlbiomcp
              MID=micpar%get_micb_id(M,NGL)
              TOMEauto(ielmc,MID)=TOMEauto(ielmc,MID)+TI*OMEauto(ielmc,MID,L,NY,NX)
              TOMEauto(ielmn,MID)=TOMEauto(ielmn,MID)+TI*OMEauto(ielmn,MID,L,NY,NX)
              TOMEauto(ielmp,MID)=TOMEauto(ielmp,MID)+TI*OMEauto(ielmp,MID,L,NY,NX)
            enddo
          enddo
        enddo

      DO K=1,jcplx
        DO  M=1,ndbiomcp
          DO NE=1,NumPlantChemElms
            TORM(NE,M,K)=TORM(NE,M,K)+TI*ORM(NE,M,K,L,NY,NX)
          ENDDO
        ENDDO
        do idom=idom_beg,idom_end
          TDOM(idom,K)=TDOM(idom,K)+TI*DOM(idom,K,L,NY,NX)
          TOHM(idom,K)=TOHM(idom,K)+TI*OHM(idom,K,L,NY,NX)
        enddo
        DO  M=1,jsken
          TOSM(ielmc,M,K)=TOSM(ielmc,M,K)+TI*OSM(ielmc,M,K,L,NY,NX)
          TOSA(M,K)=TOSA(M,K)+TI*OSA(M,K,L,NY,NX)
          TOSM(ielmn,M,K)=TOSM(ielmn,M,K)+TI*OSM(ielmn,M,K,L,NY,NX)
          TOSM(ielmp,M,K)=TOSM(ielmp,M,K)+TI*OSM(ielmp,M,K,L,NY,NX)
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
        TL=AMIN1(DLYR(3,L,NY,NX),DCORPZ-(CumSoilThickness(L,NY,NX)-DLYR(3,L,NY,NX)))
        FI=TL/DCORPZ
        TI=TL/DLYR(3,L,NY,NX)
        TX=1.0_r8-TI
        SoiBulkDensityt0(L,NY,NX)=TI*(SoiBulkDensityt0(L,NY,NX)+CORP*(TBKDX-SoiBulkDensityt0(L,NY,NX))) &
          +TX*SoiBulkDensityt0(L,NY,NX)
        FieldCapacity(L,NY,NX)=TI*(FieldCapacity(L,NY,NX)+CORP*(TFC-FieldCapacity(L,NY,NX)))+TX*FieldCapacity(L,NY,NX)
        WiltPoint(L,NY,NX)=TI*(WiltPoint(L,NY,NX)+CORP*(TWP-WiltPoint(L,NY,NX)))+TX*WiltPoint(L,NY,NX)
        SatHydroCondVert(L,NY,NX)=TI*(SatHydroCondVert(L,NY,NX)+CORP*(TSatHydroCondVert-SatHydroCondVert(L,NY,NX)))+TX*SatHydroCondVert(L,NY,NX)
        SatHydroCondHrzn(L,NY,NX)=TI*(SatHydroCondHrzn(L,NY,NX)+CORP*(TSCNH-SatHydroCondHrzn(L,NY,NX)))+TX*SatHydroCondHrzn(L,NY,NX)
        SAND(L,NY,NX)=TI*SAND(L,NY,NX)+CORP*(FI*TSAND-TI*SAND(L,NY,NX))+TX*SAND(L,NY,NX)
        SILT(L,NY,NX)=TI*SILT(L,NY,NX)+CORP*(FI*TSILT-TI*SILT(L,NY,NX))+TX*SILT(L,NY,NX)
        CLAY(L,NY,NX)=TI*CLAY(L,NY,NX)+CORP*(FI*TCLAY-TI*CLAY(L,NY,NX))+TX*CLAY(L,NY,NX)

        GKC4(L,NY,NX)=TI*(GKC4(L,NY,NX)+CORP*(TGKC4-GKC4(L,NY,NX)))+TX*GKC4(L,NY,NX)
        GKCA(L,NY,NX)=TI*(GKCA(L,NY,NX)+CORP*(TGKCA-GKCA(L,NY,NX)))+TX*GKCA(L,NY,NX)
        GKCM(L,NY,NX)=TI*(GKCM(L,NY,NX)+CORP*(TGKCM-GKCM(L,NY,NX)))+TX*GKCM(L,NY,NX)
        GKCN(L,NY,NX)=TI*(GKCN(L,NY,NX)+CORP*(TGKCN-GKCN(L,NY,NX)))+TX*GKCN(L,NY,NX)
        GKCK(L,NY,NX)=TI*(GKCK(L,NY,NX)+CORP*(TGKCK-GKCK(L,NY,NX)))+TX*GKCK(L,NY,NX)
        
        ENGYM=VHeatCapacitySoilM(L,NY,NX)*TKS(L,NY,NX)
        ENGYV=(cpw*(VLWatMicP(L,NY,NX)+VLWatMacP(L,NY,NX))+cpi*(VLiceMicP(L,NY,NX)+VLiceMacP(L,NY,NX)))*TKS(L,NY,NX)
        VLWatMicP(L,NY,NX)=TI*VLWatMicP(L,NY,NX)+CORP*(FI*TVOLW-TI*VLWatMicP(L,NY,NX))+TX*VLWatMicP(L,NY,NX)+FI*TVOLWR
        VLiceMicP(L,NY,NX)=TI*VLiceMicP(L,NY,NX)+CORP*(FI*TVOLI-TI*VLiceMicP(L,NY,NX))+TX*VLiceMicP(L,NY,NX)
      VLWatMicPX(L,NY,NX)=VLWatMicP(L,NY,NX)
!     VLWatMicP(L,NY,NX)=VLWatMicP(L,NY,NX)+CORP*VLWatMacP(L,NY,NX)
!     VLiceMicP(L,NY,NX)=VLiceMicP(L,NY,NX)+CORP*VLiceMacP(L,NY,NX)
!     VLMicP(L,NY,NX)=VLMicP(L,NY,NX)+CORP*VLMacP(L,NY,NX)
!     VLWatMacP(L,NY,NX)=XCORP(NY,NX)*VLWatMacP(L,NY,NX)
!     VLiceMacP(L,NY,NX)=XCORP(NY,NX)*VLiceMacP(L,NY,NX)
!     VLMacP(L,NY,NX)=XCORP(NY,NX)*VLMacP(L,NY,NX)
!     SoilFracAsMacP(L,NY,NX)=XCORP(NY,NX)*SoilFracAsMacP(L,NY,NX)
        ENGYL=TI*ENGYV+CORP*(FI*TENGY-TI*ENGYV)+TX*ENGYV+FI*TENGYR
        VHeatCapacity(L,NY,NX)=VHeatCapacitySoilM(L,NY,NX)+cpw*(VLWatMicP(L,NY,NX)+VLWatMacP(L,NY,NX)) &
          +cpi*(VLiceMicP(L,NY,NX)+VLiceMacP(L,NY,NX))
        TKS(L,NY,NX)=(ENGYM+ENGYL)/VHeatCapacity(L,NY,NX)
        TCS(L,NY,NX)=units%Kelvin2Celcius(TKS(L,NY,NX))
        DO NTN=ifertn_beg,ifertn_end
          FertN_soil_vr(NTN,L,NY,NX)=TI*FertN_soil_vr(NTN,L,NY,NX) &
            +CORP*(FI*TfertN_soil(NTN)-TI*FertN_soil_vr(NTN,L,NY,NX))&
            +TX*FertN_soil_vr(NTN,L,NY,NX)
        ENDDO

        DO NTN=ifertnb_beg,ifertnb_end
          FertN_Band_vr(NTN,L,NY,NX)=TI*FertN_Band_vr(NTN,L,NY,NX) &
            +CORP*(FI*TFertN_Band(NTN)-TI*FertN_Band_vr(NTN,L,NY,NX)) &
            +TX*FertN_Band_vr(NTN,L,NY,NX)
        ENDDO

        !SALT
        DO NTSA=idsalt_beg,idsalt_end
          trcSalt_solml(NTSA,L,NY,NX)=TI*trcSalt_solml(NTSA,L,NY,NX)  &
            +CORP*(FI*TSA_solml(NTSA)-TI*trcSalt_solml(NTSA,L,NY,NX)) &
            +TX*trcSalt_solml(NTSA,L,NY,NX)+CORP*trcSalt_soHml(NTSA,L,NY,NX)
        ENDDO

        ! solute
        DO NTS=ids_beg,ids_end
          trc_solml_vr(NTS,L,NY,NX)=TI*trc_solml_vr(NTS,L,NY,NX) &
            +CORP*(FI*TS_solml(NTS)-TI*trc_solml_vr(NTS,L,NY,NX)) &
            +TX*trc_solml_vr(NTS,L,NY,NX)+CORP*trc_soHml(NTS,L,NY,NX)
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
          trc_gasml_vr(NTG,L,NY,NX)=TI*trc_gasml_vr(NTG,L,NY,NX)+CORP*(FI*TG_gasml(NTG) &
            -TI*trc_gasml_vr(NTG,L,NY,NX))+TX*trc_gasml_vr(NTG,L,NY,NX)
        ENDDO

        call MixSoluteH(L,NY,NX)

        DO  K=1,jcplx
          DO  N=1,NumMicbFunGroups
            DO NGL=JGnio(N),JGnfo(N)
              DO  M=1,nlbiomcp
                MID=micpar%get_micb_id(M,NGL)
                OMEheter(ielmc,MID,K,L,NY,NX)=TI*OMEheter(ielmc,MID,K,L,NY,NX)+CORP*(FI*TOMEheter(ielmc,MID,K) &
                  -TI*OMEheter(ielmc,MID,K,L,NY,NX))+TX*OMEheter(ielmc,MID,K,L,NY,NX)
                OMEheter(ielmn,MID,K,L,NY,NX)=TI*OMEheter(ielmn,MID,K,L,NY,NX)+CORP*(FI*TOMEheter(ielmn,MID,K) &
                  -TI*OMEheter(ielmn,MID,K,L,NY,NX))+TX*OMEheter(ielmn,MID,K,L,NY,NX)
                OMEheter(ielmp,MID,K,L,NY,NX)=TI*OMEheter(ielmp,MID,K,L,NY,NX)+CORP*(FI*TOMEheter(ielmp,MID,K) &
                  -TI*OMEheter(ielmp,MID,K,L,NY,NX))+TX*OMEheter(ielmp,MID,K,L,NY,NX)
              enddo
            enddo
          enddo
        ENDDO
        DO  N=1,NumMicbFunGroups
          DO NGL=JGniA(N),JGnfA(N)
            DO  M=1,nlbiomcp
              MID=micpar%get_micb_id(M,NGL)
              OMEauto(ielmc,MID,L,NY,NX)=TI*OMEauto(ielmc,MID,L,NY,NX)+CORP*(FI*TOMEauto(ielmc,MID) &
                -TI*OMEauto(ielmc,MID,L,NY,NX))+TX*OMEauto(ielmc,MID,L,NY,NX)
              OMEauto(ielmn,MID,L,NY,NX)=TI*OMEauto(ielmn,MID,L,NY,NX)+CORP*(FI*TOMEauto(ielmn,MID) &
                -TI*OMEauto(ielmn,MID,L,NY,NX))+TX*OMEauto(ielmn,MID,L,NY,NX)
              OMEauto(ielmp,MID,L,NY,NX)=TI*OMEauto(ielmp,MID,L,NY,NX)+CORP*(FI*TOMEauto(ielmp,MID) &
                -TI*OMEauto(ielmp,MID,L,NY,NX))+TX*OMEauto(ielmp,MID,L,NY,NX)
            enddo
          enddo
        enddo

        DO  K=1,jcplx
          DO  M=1,ndbiomcp
            ORM(ielmc,M,K,L,NY,NX)=TI*ORM(ielmc,M,K,L,NY,NX)+CORP*(FI*TORM(ielmc,M,K) &
              -TI*ORM(ielmc,M,K,L,NY,NX))+TX*ORM(ielmc,M,K,L,NY,NX)
            ORM(ielmn,M,K,L,NY,NX)=TI*ORM(ielmn,M,K,L,NY,NX)+CORP*(FI*TORM(ielmn,M,K) &
              -TI*ORM(ielmn,M,K,L,NY,NX))+TX*ORM(ielmn,M,K,L,NY,NX)
            ORM(ielmp,M,K,L,NY,NX)=TI*ORM(ielmp,M,K,L,NY,NX)+CORP*(FI*TORM(ielmp,M,K) &
              -TI*ORM(ielmp,M,K,L,NY,NX))+TX*ORM(ielmp,M,K,L,NY,NX)
          ENDDO
          DOM(idom_doc,K,L,NY,NX)=TI*DOM(idom_doc,K,L,NY,NX)+CORP*(FI*TDOM(idom_doc,K) &
            -TI*DOM(idom_doc,K,L,NY,NX))+TX*DOM(idom_doc,K,L,NY,NX)+CORP*DOM_Macp(idom_doc,K,L,NY,NX)
          DOM(idom_don,K,L,NY,NX)=TI*DOM(idom_don,K,L,NY,NX)+CORP*(FI*TDOM(idom_don,K) &
            -TI*DOM(idom_don,K,L,NY,NX))+TX*DOM(idom_don,K,L,NY,NX)+CORP*DOM_Macp(idom_don,K,L,NY,NX)
          DOM(idom_dop,K,L,NY,NX)=TI*DOM(idom_dop,K,L,NY,NX)+CORP*(FI*TDOM(idom_dop,K) &
            -TI*DOM(idom_dop,K,L,NY,NX))+TX*DOM(idom_dop,K,L,NY,NX)+CORP*DOM_Macp(idom_dop,K,L,NY,NX)
          DOM(idom_acetate,K,L,NY,NX)=TI*DOM(idom_acetate,K,L,NY,NX)+CORP*(FI*TDOM(idom_acetate,K) &
            -TI*DOM(idom_acetate,K,L,NY,NX))+TX*DOM(idom_acetate,K,L,NY,NX)+CORP*DOM_Macp(idom_acetate,K,L,NY,NX)

          DOM_Macp(idom_doc,K,L,NY,NX)=XCORP(NY,NX)*DOM_Macp(idom_doc,K,L,NY,NX)
          DOM_Macp(idom_don,K,L,NY,NX)=XCORP(NY,NX)*DOM_Macp(idom_don,K,L,NY,NX)
          DOM_Macp(idom_dop,K,L,NY,NX)=XCORP(NY,NX)*DOM_Macp(idom_dop,K,L,NY,NX)
          DOM_Macp(idom_acetate,K,L,NY,NX)=XCORP(NY,NX)*DOM_Macp(idom_acetate,K,L,NY,NX)

          OHM(ielmc,K,L,NY,NX)=TI*OHM(ielmc,K,L,NY,NX)+CORP*(FI*TOHM(ielmc,K) &
            -TI*OHM(ielmc,K,L,NY,NX))+TX*OHM(ielmc,K,L,NY,NX)
          OHM(ielmn,K,L,NY,NX)=TI*OHM(ielmn,K,L,NY,NX)+CORP*(FI*TOHM(ielmn,K) &
            -TI*OHM(ielmn,K,L,NY,NX))+TX*OHM(ielmn,K,L,NY,NX)
          OHM(ielmp,K,L,NY,NX)=TI*OHM(ielmp,K,L,NY,NX)+CORP*(FI*TOHM(ielmp,K) &
            -TI*OHM(ielmp,K,L,NY,NX))+TX*OHM(ielmp,K,L,NY,NX)
          OHM(idom_acetate,K,L,NY,NX)=TI*OHM(idom_acetate,K,L,NY,NX)+CORP*(FI*TOHM(idom_acetate,K) &
            -TI*OHM(idom_acetate,K,L,NY,NX))+TX*OHM(idom_acetate,K,L,NY,NX)
          DO  M=1,jsken
            OSM(ielmc,M,K,L,NY,NX)=TI*OSM(ielmc,M,K,L,NY,NX)+CORP*(FI*TOSM(ielmc,M,K) &
              -TI*OSM(ielmc,M,K,L,NY,NX))+TX*OSM(ielmc,M,K,L,NY,NX)
            OSA(M,K,L,NY,NX)=TI*OSA(M,K,L,NY,NX)+CORP*(FI*TOSA(M,K) &
              -TI*OSA(M,K,L,NY,NX))+TX*OSA(M,K,L,NY,NX)
            OSM(ielmn,M,K,L,NY,NX)=TI*OSM(ielmn,M,K,L,NY,NX)+CORP*(FI*TOSM(ielmn,M,K) &
              -TI*OSM(ielmn,M,K,L,NY,NX))+TX*OSM(ielmn,M,K,L,NY,NX)
            OSM(ielmp,M,K,L,NY,NX)=TI*OSM(ielmp,M,K,L,NY,NX)+CORP*(FI*TOSM(ielmp,M,K) &
              -TI*OSM(ielmp,M,K,L,NY,NX))+TX*OSM(ielmp,M,K,L,NY,NX)
          ENDDO
        ENDDO
!
!     ADD STATE VARIABLES IN SURFACE RESIDUE INCORPORATED
!     WITHIN TILLAGE MIXING ZONE
!
        DO  K=1,micpar%NumOfLitrCmplxs
          DO  N=1,NumMicbFunGroups
            DO NGL=JGnio(N),JGnfo(N)
              DO M=1,nlbiomcp
                MID=micpar%get_micb_id(M,NGL)
                OMEheter(ielmc,MID,K,L,NY,NX)=OMEheter(ielmc,MID,K,L,NY,NX)+FI*TOMGC(MID,K)
                OMEheter(ielmn,MID,K,L,NY,NX)=OMEheter(ielmn,MID,K,L,NY,NX)+FI*TOMGN(MID,K)
                OMEheter(ielmp,MID,K,L,NY,NX)=OMEheter(ielmp,MID,K,L,NY,NX)+FI*TOMGP(MID,K)
              enddo
            enddo
          ENDDO
        ENDDO

        DO  N=1,NumMicbFunGroups
          DO NGL=JGniA(N),JGnfA(N)
            DO M=1,nlbiomcp
              MID=micpar%get_micb_id(M,NGL)
              OMEauto(ielmc,MID,L,NY,NX)=OMEauto(ielmc,MID,L,NY,NX)+FI*TOMGCff(MID)
              OMEauto(ielmn,MID,L,NY,NX)=OMEauto(ielmn,MID,L,NY,NX)+FI*TOMGNff(MID)
              OMEauto(ielmp,MID,L,NY,NX)=OMEauto(ielmp,MID,L,NY,NX)+FI*TOMGPff(MID)
            enddo
          enddo
        ENDDO

        DO K=1,micpar%NumOfLitrCmplxs
          DO  M=1,ndbiomcp
            ORM(ielmc,M,K,L,NY,NX)=ORM(ielmc,M,K,L,NY,NX)+FI*TORXC(M,K)
            ORM(ielmn,M,K,L,NY,NX)=ORM(ielmn,M,K,L,NY,NX)+FI*TORXN(M,K)
            ORM(ielmp,M,K,L,NY,NX)=ORM(ielmp,M,K,L,NY,NX)+FI*TORXP(M,K)
          ENDDO
          DOM(idom_doc,K,L,NY,NX)=DOM(idom_doc,K,L,NY,NX)+FI*TOQGC(K)
          DOM(idom_don,K,L,NY,NX)=DOM(idom_don,K,L,NY,NX)+FI*TOQGN(K)
          DOM(idom_dop,K,L,NY,NX)=DOM(idom_dop,K,L,NY,NX)+FI*TOQGP(K)
          DOM(idom_acetate,K,L,NY,NX)=DOM(idom_acetate,K,L,NY,NX)+FI*TOQGA(K)
          DOM_Macp(idom_doc,K,L,NY,NX)=DOM_Macp(idom_doc,K,L,NY,NX)+FI*TOQHC(K)
          DOM_Macp(idom_don,K,L,NY,NX)=DOM_Macp(idom_don,K,L,NY,NX)+FI*TOQHN(K)
          DOM_Macp(idom_dop,K,L,NY,NX)=DOM_Macp(idom_dop,K,L,NY,NX)+FI*TOQHP(K)
          DOM_Macp(idom_acetate,K,L,NY,NX)=DOM_Macp(idom_acetate,K,L,NY,NX)+FI*TOQHA(K)
          OHM(ielmc,K,L,NY,NX)=OHM(ielmc,K,L,NY,NX)+FI*TOHGC(K)
          OHM(ielmn,K,L,NY,NX)=OHM(ielmn,K,L,NY,NX)+FI*TOHGN(K)
          OHM(ielmp,K,L,NY,NX)=OHM(ielmp,K,L,NY,NX)+FI*TOHGP(K)
          OHM(idom_acetate,K,L,NY,NX)=OHM(idom_acetate,K,L,NY,NX)+FI*TOHGA(K)
          DO  M=1,jsken
            OSM(ielmc,M,K,L,NY,NX)=OSM(ielmc,M,K,L,NY,NX)+FI*TOSGC(M,K)
            OSA(M,K,L,NY,NX)=OSA(M,K,L,NY,NX)+FI*TOSGA(M,K)
            OSM(ielmn,M,K,L,NY,NX)=OSM(ielmn,M,K,L,NY,NX)+FI*TOSGN(M,K)
            OSM(ielmp,M,K,L,NY,NX)=OSM(ielmp,M,K,L,NY,NX)+FI*TOSGP(M,K)
          ENDDO
        ENDDO
        OC=0.0_r8
        ON=0.0_r8
        OP=0.0_r8
        DC=0.0_r8
        DN=0.0_r8
        DP=0.0_r8

        DO  K=1,jcplx
          DO  N=1,NumMicbFunGroups
            DO NGL=JGnio(N),JGnfo(N)
              DO  M=1,nlbiomcp
                MID=micpar%get_micb_id(M,NGL)
                OC=OC+OMEheter(ielmc,MID,K,L,NY,NX)
                ON=ON+OMEheter(ielmn,MID,K,L,NY,NX)
                OP=OP+OMEheter(ielmp,MID,K,L,NY,NX)
              enddo
            enddo
          enddo
        ENDDO

        DO  K=1,micpar%NumOfLitrCmplxs
          DO  N=1,NumMicbFunGroups
            DO NGL=JGnio(N),JGnfo(N)
              DO  M=1,nlbiomcp
                MID=micpar%get_micb_id(M,NGL)
                DC=DC+OMEheter(ielmc,MID,K,L,NY,NX)
                DN=DN+OMEheter(ielmn,MID,K,L,NY,NX)
                DP=DP+OMEheter(ielmp,MID,K,L,NY,NX)
              enddo
            enddo
          enddo
        ENDDO
        DO  N=1,NumMicbFunGroups
          DO NGL=JGniA(N),JGnfA(N)
            DO  M=1,nlbiomcp
              MID=micpar%get_micb_id(M,NGL)
              OC=OC+OMEauto(ielmc,MID,L,NY,NX)
              ON=ON+OMEauto(ielmn,MID,L,NY,NX)
              OP=OP+OMEauto(ielmp,MID,L,NY,NX)
            enddo
          enddo
        enddo

        DO K=1,jcplx
          DO  M=1,ndbiomcp
            OC=OC+ORM(ielmc,M,K,L,NY,NX)
            ON=ON+ORM(ielmn,M,K,L,NY,NX)
            OP=OP+ORM(ielmp,M,K,L,NY,NX)
            IF(micpar%is_litter(K))THEN
              DC=DC+ORM(ielmc,M,K,L,NY,NX)
              DN=DN+ORM(ielmn,M,K,L,NY,NX)
              DP=DP+ORM(ielmp,M,K,L,NY,NX)
            ENDIF
          ENDDO

          OC=OC+DOM(idom_doc,K,L,NY,NX)+DOM_Macp(idom_doc,K,L,NY,NX)+OHM(ielmc,K,L,NY,NX) &
            +DOM(idom_acetate,K,L,NY,NX)+DOM_Macp(idom_acetate,K,L,NY,NX)+OHM(idom_acetate,K,L,NY,NX)
          ON=ON+DOM(idom_don,K,L,NY,NX)+DOM_Macp(idom_don,K,L,NY,NX)+OHM(ielmn,K,L,NY,NX)
          OP=OP+DOM(idom_dop,K,L,NY,NX)+DOM_Macp(idom_dop,K,L,NY,NX)+OHM(ielmp,K,L,NY,NX)
          IF(micpar%is_litter(K))THEN
            DC=DC+DOM(idom_doc,K,L,NY,NX)+DOM_Macp(idom_doc,K,L,NY,NX)+OHM(ielmc,K,L,NY,NX) &
                +DOM(idom_acetate,K,L,NY,NX)+DOM_Macp(idom_acetate,K,L,NY,NX)+OHM(idom_acetate,K,L,NY,NX)
            DN=DN+DOM(idom_don,K,L,NY,NX)+DOM_Macp(idom_don,K,L,NY,NX)+OHM(ielmn,K,L,NY,NX)
            DC=DC+DOM(idom_doc,K,L,NY,NX)+DOM_Macp(idom_doc,K,L,NY,NX)+OHM(ielmc,K,L,NY,NX)
          ENDIF
          DO  M=1,jsken
            OC=OC+OSM(ielmc,M,K,L,NY,NX)
            ON=ON+OSM(ielmn,M,K,L,NY,NX)
            OP=OP+OSM(ielmp,M,K,L,NY,NX)
            IF(micpar%is_litter(K))THEN
              DC=DC+OSM(ielmc,M,K,L,NY,NX)
              DN=DN+OSM(ielmn,M,K,L,NY,NX)
              DP=DP+OSM(ielmp,M,K,L,NY,NX)
            ENDIF
          ENDDO
        ENDDO
        ORGC(L,NY,NX)=OC
        ORGN(L,NY,NX)=ON
        ORGR(L,NY,NX)=DC

        DO NTS=ids_beg,idg_NH3
          trc_solml_vr(NTS,L,NY,NX)=trc_solml_vr(NTS,L,NY,NX)+FI*TS0_solml(NTS)
        ENDDO

        DO NTS=ids_nut_beg,ids_nuts_end
          trc_solml_vr(NTS,L,NY,NX)=trc_solml_vr(NTS,L,NY,NX)+FI*TS0_solml(NTS)
        ENDDO

        trcx_solml(idx_NH4,L,NY,NX)=trcx_solml(idx_NH4,L,NY,NX)+FI*TXN4G


        DO NTX=idx_AEC+1,idx_anion_soil_end
          trcx_solml(NTX,L,NY,NX)=trcx_solml(NTX,L,NY,NX)+FI*TX_anion(NTX)
        ENDDO

        DO NTP=idsp_psoi_beg,idsp_psoi_end
          trcp_salml(NTP,L,NY,NX)=trcp_salml(NTP,L,NY,NX)+FI*TP_soil(NTP)
        ENDDO

        DO NTN=ifertn_beg,ifertn_end
          FertN_soil_vr(NTN,L,NY,NX)=FertN_soil_vr(NTN,L,NY,NX)+FI*TFertNG_soil(NTN)
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

  DO NTSA=idsalt_beg,idsalt_end
    trcSalt_soHml(ntsa,L,NY,NX)=XCORP(NY,NX)*trcSalt_soHml(ntsa,L,NY,NX)
  ENDDO
  end subroutine MixSoluteH


end module TillageMixMod
