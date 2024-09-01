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
  use SoilBGCNLayMod, only : sumLitrOMLayL,sumORGMLayL
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

  REAL(R8) :: TOSGE(1:NumPlantChemElms,jsken,0:micpar%NumOfLitrCmplxs),TOSGA(jsken,0:micpar%NumOfLitrCmplxs)
  real(r8) :: TORXE(1:NumPlantChemElms,ndbiomcp,0:micpar%NumOfLitrCmplxs)
  real(r8) :: TDOMG(idom_beg:idom_end,0:micpar%NumOfLitrCmplxs),TDOMH(idom_beg:idom_end,0:micpar%NumOfLitrCmplxs)
  real(r8) :: TDOMHG(idom_beg:idom_end,0:micpar%NumOfLitrCmplxs)

  real(r8) :: TOMGAutor(1:NumPlantChemElms,1:NumLiveAutoBioms)
  real(r8) :: TBOMGE(1:NumPlantChemElms,1:NumLiveHeterBioms,1:jcplx)
  REAL(R8) :: TmBiomeHeter(NumPlantChemElms,NumLiveHeterBioms,1:jcplx)
  REAL(R8) :: TmBiomeAutor(NumPlantChemElms,NumLiveAutoBioms)
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
  real(r8) :: ORGM(1:NumPlantChemElms),litrOM(1:NumPlantChemElms)
  real(r8) :: TVOLI,HFLXD
  real(r8) :: FI,TI,TX,TL
  integer  :: NTX,NTP,NTG,NTSA,NTN,idom,MID,NE

!     begin_execution
!
  !DCORP=soil mixing fraction with tillage
  !XCORP=factor for surface litter incorporation and soil mixing
  
  IF(J.EQ.INT(SolarNoonHour_col(NY,NX)) .AND. XCORP(NY,NX).LT.1.0 .AND. DCORP(I,NY,NX).GT.0.0_r8)THEN
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

    TmBiomeHeter=0.0_r8
    TmBiomeAutor=0.0_r8
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
!     IF(SoilOrgM_vr(ielmc,0,NY,NX).GT.ZEROS(NY,NX))THEN
!     XCORP0=AMAX1(XCORP(NY,NX),AMIN1(1.0,
!    2(VHeatCapLitR(NY,NX)/cpo)/SoilOrgM_vr(ielmc,0,NY,NX)))
    XCORP0=AMAX1(0.001_r8,XCORP(NY,NX))
!     ELSE
!     XCORP0=1.0
!     ENDIF
    CORP0=1.0_r8-XCORP0   !fraction mixed into next layer

    DO  K=1,micpar%NumOfLitrCmplxs
      DO  N=1,NumMicbFunGrupsPerCmplx
        DO NGL=JGnio(N),JGnfo(N)
          DO  M=1,nlbiomcp
            MID=micpar%get_micb_id(M,NGL)
            DO NE=1,NumPlantChemElms
              TBOMGE(NE,MID,K)=mBiomeHeter_vr(NE,MID,K,0,NY,NX)*CORP0
              mBiomeHeter_vr(NE,MID,K,0,NY,NX)=mBiomeHeter_vr(NE,MID,K,0,NY,NX)*XCORP0
            ENDDO
          enddo
        enddo
      ENDDO
    ENDDO

    DO  N=1,NumMicbFunGrupsPerCmplx
      DO NGL=JGniA(N),JGnfA(N)
        DO  M=1,nlbiomcp
          MID=micpar%get_micb_id(M,NGL)
          DO NE=1,NumPlantChemElms
            TOMGAutor(NE,MID)=mBiomeAutor_vr(NE,MID,0,NY,NX)*CORP0
            mBiomeAutor_vr(NE,MID,0,NY,NX)=mBiomeAutor_vr(NE,MID,0,NY,NX)*XCORP0
          ENDDO
        enddo
      enddo
    ENDDO

    DO K=1,micpar%NumOfLitrCmplxs
      DO M=1,ndbiomcp
        DO NE=1,NumPlantChemElms
          TORXE(NE,M,K)=OMBioResdu_vr(NE,M,K,0,NY,NX)*CORP0
          OMBioResdu_vr(NE,M,K,0,NY,NX)=OMBioResdu_vr(NE,M,K,0,NY,NX)*XCORP0
        ENDDO

      ENDDO

      DO idom=idom_beg,idom_end
        TDOMG(idom,K)=DOM_vr(idom,K,0,NY,NX)*CORP0
        TDOMH(idom,K)=DOM_MacP_vr(idom,K,0,NY,NX)*CORP0
        TDOMHG(idom,K)=SorbedOM_vr(idom,K,0,NY,NX)*CORP0
  !
  !     REDUCE SURFACE RESIDUE STATE VARIABLES FOR INCORPORATION
  !
        DOM_vr(idom,K,0,NY,NX)=DOM_vr(idom,K,0,NY,NX)*XCORP0
        
        DOM_MacP_vr(idom,K,0,NY,NX)=DOM_MacP_vr(idom,K,0,NY,NX)*XCORP0
        SorbedOM_vr(idom,K,0,NY,NX)=SorbedOM_vr(idom,K,0,NY,NX)*XCORP0
      enddo 

      DO  M=1,jsken
        TOSGA(M,K)=SolidOMAct_vr(M,K,0,NY,NX)*CORP0      
        SolidOMAct_vr(M,K,0,NY,NX)=SolidOMAct_vr(M,K,0,NY,NX)*XCORP0

        DO NE=1,NumPlantChemElms
          TOSGE(NE,M,K)=SolidOM_vr(NE,M,K,0,NY,NX)*CORP0
          SolidOM_vr(NE,M,K,0,NY,NX)=SolidOM_vr(NE,M,K,0,NY,NX)*XCORP0
        ENDDO

      ENDDO
    ENDDO

    DO NTS=ids_beg,idg_NH3
      TS0_solml(NTS)=trc_solml_vr(NTS,0,NY,NX)*CORP0
    ENDDO
    DO NTS=ids_nut_beg,ids_nuts_end
      TS0_solml(NTS)=trc_solml_vr(NTS,0,NY,NX)*CORP0
    ENDDO

    TXN4G=trcx_solml_vr(idx_NH4,0,NY,NX)*CORP0

    DO NTX= idx_AEC+1,idx_anion_soil_end
      TX_anion(NTX)=trcx_solml_vr(NTX,0,NY,NX)*CORP0
    ENDDO

    DO NTP=idsp_psoi_beg,idsp_psoi_end
      TP_soil(NTP)=trcp_saltpml_vr(NTP,0,NY,NX)*CORP0
    ENDDO

    DO NTN=ifertn_beg,ifertn_end
      TFertNG_soil(NTN)=FertN_soil_vr(NTN,0,NY,NX)*CORP0
    ENDDO

    TZNFNG=ZNFNI(0,NY,NX)*CORP0
    TVOLWR=VLWatMicP_vr(0,NY,NX)*CORP0
    HFLXD=cpo*SoilOrgM_vr(ielmc,0,NY,NX)*CORP0*TKS_vr(0,NY,NX)
    HEATIN_lnd=HEATIN_lnd-HFLXD
    HeatStore_lnd=HeatStore_lnd-HFLXD
    TENGYR=cpw*TVOLWR*TKS_vr(0,NY,NX)

    call sumLitrOMLayL(0,NY,NX,litrOM)

    SoilOrgM_vr(1:NumPlantChemElms,0,NY,NX)=litrOM(1:NumPlantChemElms)

    OMLitrC_vr(0,NY,NX)=litrOM(ielmc)

    DO NTS=ids_beg,idg_NH3
      trc_solml_vr(NTS,0,NY,NX)=trc_solml_vr(NTS,0,NY,NX)*XCORP0
    ENDDO

    DO NTS=ids_nut_beg,ids_nuts_end
      trc_solml_vr(NTS,0,NY,NX)=trc_solml_vr(NTS,0,NY,NX)*XCORP0
    ENDDO

    trcx_solml_vr(idx_NH4,0,NY,NX)=trcx_solml_vr(idx_NH4,0,NY,NX)*XCORP0

    DO NTX=idx_AEC+1,idx_anion_soil_end
      trcx_solml_vr(idx_OHe,0,NY,NX)=trcx_solml_vr(idx_OHe,0,NY,NX)*XCORP0
    ENDDO

    trcp_saltpml_vr(idsp_AlPO4,0,NY,NX)    = trcp_saltpml_vr(idsp_AlPO4,0,NY,NX)*XCORP0
    trcp_saltpml_vr(idsp_FePO4,0,NY,NX)    = trcp_saltpml_vr(idsp_FePO4,0,NY,NX)*XCORP0
    trcp_saltpml_vr(idsp_CaHPO4,0,NY,NX)   = trcp_saltpml_vr(idsp_CaHPO4,0,NY,NX)*XCORP0
    trcp_saltpml_vr(idsp_HA,0,NY,NX)       = trcp_saltpml_vr(idsp_HA,0,NY,NX)*XCORP0
    trcp_saltpml_vr(idsp_CaH4P2O8,0,NY,NX) = trcp_saltpml_vr(idsp_CaH4P2O8,0,NY,NX)*XCORP0

    DO NTN=ifertn_beg,ifertn_end
      FertN_soil_vr(NTN,0,NY,NX)=FertN_soil_vr(NTN,0,NY,NX)*XCORP0
    ENDDO

    VLWatMicP_vr(0,NY,NX)=VLWatMicP_vr(0,NY,NX)*XCORP0
    VHeatCapacity_vr(0,NY,NX)=cpo*SoilOrgM_vr(ielmc,0,NY,NX)+cpw*VLWatMicP_vr(0,NY,NX)+cpi*VLiceMicP_vr(0,NY,NX)
    VLitR(NY,NX)=VLitR(NY,NX)*XCORP0
    VGeomLayer_vr(0,NY,NX)=VGeomLayer_vr(0,NY,NX)*XCORP0
    ZNHUX0=AMAX1(ZNHUX0,ZNHU0(0,NY,NX))
    ZNHUXI=AMAX1(ZNHUXI,ZNHUI(0,NY,NX))
    ZNFNX0=AMAX1(ZNFNX0,ZNFN0(0,NY,NX))
    LL=NU(NY,NX)
!
!     REDISTRIBUTE SOIL STATE VARIABLES DURING TILLAGE
!
    DCORPZ=AMIN1(DCORP(I,NY,NX),CumSoilThickness_vr(NL(NY,NX),NY,NX))
!
!     ACCUMULATE SOIL STATE VARIABLES WITHIN TILLAGE MIXING ZONE
!
    DO L=NU(NY,NX),NL(NY,NX)
      IF(CumSoilThickness_vr(L,NY,NX)-DLYR(3,L,NY,NX).LT.DCORPZ.AND.DLYR(3,L,NY,NX).GT.ZERO)THEN
        TL=AMIN1(DLYR(3,L,NY,NX),DCORPZ-(CumSoilThickness_vr(L,NY,NX)-DLYR(3,L,NY,NX)))
        FI=TL/DCORPZ
        TI=TL/DLYR(3,L,NY,NX)
        TBKDX=TBKDX+FI*SoiBulkDensityt0_vr(L,NY,NX)
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
        TVOLW=TVOLW+TI*VLWatMicP_vr(L,NY,NX)
        TVOLI=TVOLI+TI*VLiceMicP_vr(L,NY,NX)
!     TVOLP=TVOLP+TI*VLsoiAirP_vr(L,NY,NX)
!     TVOLA=TVOLA+TI*VLMicP_vr(L,NY,NX)
        TENGY=TENGY+TI*(cpw*(VLWatMicP_vr(L,NY,NX)+VLWatMacP_vr(L,NY,NX)) &
          +cpi*(VLiceMicP_vr(L,NY,NX)+VLiceMacP_vr(L,NY,NX)))*TKS_vr(L,NY,NX)
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
          TSA_solml(NTSA)=TSA_solml(NTSA)+TI*trcSalt_solml_vr(NTSA,L,NY,NX)
        ENDDO
!cation
        DO NTX=idx_CEC,idx_cation_end
          TX_solml(NTX)=TX_solml(NTX)+TI*trcx_solml_vr(NTX,L,NY,NX)
        ENDDO
!anion
        DO NTX=idx_AEC,idx_end
          TX_solml(NTX)=TX_solml(NTX)+TI*trcx_solml_vr(NTX,L,NY,NX)
        ENDDO

        DO NTP=idsp_beg,idsp_end
          TP_salml(NTP)=TP_salml(NTP)+TI*trcp_saltpml_vr(NTP,L,NY,NX)
        ENDDO

        DO NTG=idg_beg,idg_end-1
          TG_gasml(NTG)=TG_gasml(NTG)+TI*trc_gasml_vr(NTG,L,NY,NX)
        ENDDO

        DO  K=1,jcplx
          DO  N=1,NumMicbFunGrupsPerCmplx
            DO NGL=JGnio(N),JGnfo(N)
              DO  M=1,nlbiomcp
                MID=micpar%get_micb_id(M,NGL)
                DO NE=1,NumPlantChemElms
                  TmBiomeHeter(NE,MID,K)=TmBiomeHeter(NE,MID,K)+TI*mBiomeHeter_vr(NE,MID,K,L,NY,NX)
                ENDDO
              enddo
            enddo
          enddo
        ENDDO

        DO  N=1,NumMicbFunGrupsPerCmplx
          DO NGL=JGniA(N),JGnfA(N)
            DO  M=1,nlbiomcp
              MID=micpar%get_micb_id(M,NGL)
              DO NE=1,NumPlantChemElms
              TmBiomeAutor(NE,MID)=TmBiomeAutor(NE,MID)+TI*mBiomeAutor_vr(NE,MID,L,NY,NX)
              ENDDO
            enddo
          enddo
        enddo

      DO K=1,jcplx
        DO  M=1,ndbiomcp
          DO NE=1,NumPlantChemElms
            TORM(NE,M,K)=TORM(NE,M,K)+TI*OMBioResdu_vr(NE,M,K,L,NY,NX)
          ENDDO
        ENDDO
        do idom=idom_beg,idom_end
          TDOM(idom,K)=TDOM(idom,K)+TI*DOM_vr(idom,K,L,NY,NX)
          TOHM(idom,K)=TOHM(idom,K)+TI*SorbedOM_vr(idom,K,L,NY,NX)
        enddo
        DO  M=1,jsken
          TOSA(M,K)=TOSA(M,K)+TI*SolidOMAct_vr(M,K,L,NY,NX)
          DO NE=1,NumPlantChemElms
            TOSM(NE,M,K)=TOSM(NE,M,K)+TI*SolidOM_vr(NE,M,K,L,NY,NX)
          ENDDO
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
        TL=AMIN1(DLYR(3,L,NY,NX),DCORPZ-(CumSoilThickness_vr(L,NY,NX)-DLYR(3,L,NY,NX)))
        FI=TL/DCORPZ
        TI=TL/DLYR(3,L,NY,NX)
        TX=1.0_r8-TI
        SoiBulkDensityt0_vr(L,NY,NX) = TI*(SoiBulkDensityt0_vr(L,NY,NX)+CORP*(TBKDX-SoiBulkDensityt0_vr(L,NY,NX))) &
          +TX*SoiBulkDensityt0_vr(L,NY,NX)
        FieldCapacity(L,NY,NX)    = TI*(FieldCapacity(L,NY,NX)+CORP*(TFC-FieldCapacity(L,NY,NX)))+TX*FieldCapacity(L,NY,NX)
        WiltPoint(L,NY,NX)        = TI*(WiltPoint(L,NY,NX)+CORP*(TWP-WiltPoint(L,NY,NX)))+TX*WiltPoint(L,NY,NX)
        SatHydroCondVert(L,NY,NX) = TI*(SatHydroCondVert(L,NY,NX)+CORP*(TSatHydroCondVert-SatHydroCondVert(L,NY,NX)))&
          +TX*SatHydroCondVert(L,NY,NX)
        SatHydroCondHrzn(L,NY,NX) = TI*(SatHydroCondHrzn(L,NY,NX)+CORP*(TSCNH-SatHydroCondHrzn(L,NY,NX))) &
          +TX*SatHydroCondHrzn(L,NY,NX)
        SAND(L,NY,NX)=TI*SAND(L,NY,NX)+CORP*(FI*TSAND-TI*SAND(L,NY,NX))+TX*SAND(L,NY,NX)
        SILT(L,NY,NX)=TI*SILT(L,NY,NX)+CORP*(FI*TSILT-TI*SILT(L,NY,NX))+TX*SILT(L,NY,NX)
        CLAY(L,NY,NX)=TI*CLAY(L,NY,NX)+CORP*(FI*TCLAY-TI*CLAY(L,NY,NX))+TX*CLAY(L,NY,NX)

        GKC4(L,NY,NX)=TI*(GKC4(L,NY,NX)+CORP*(TGKC4-GKC4(L,NY,NX)))+TX*GKC4(L,NY,NX)
        GKCA(L,NY,NX)=TI*(GKCA(L,NY,NX)+CORP*(TGKCA-GKCA(L,NY,NX)))+TX*GKCA(L,NY,NX)
        GKCM(L,NY,NX)=TI*(GKCM(L,NY,NX)+CORP*(TGKCM-GKCM(L,NY,NX)))+TX*GKCM(L,NY,NX)
        GKCN(L,NY,NX)=TI*(GKCN(L,NY,NX)+CORP*(TGKCN-GKCN(L,NY,NX)))+TX*GKCN(L,NY,NX)
        GKCK(L,NY,NX)=TI*(GKCK(L,NY,NX)+CORP*(TGKCK-GKCK(L,NY,NX)))+TX*GKCK(L,NY,NX)
        
        ENGYM=VHeatCapacitySoilM(L,NY,NX)*TKS_vr(L,NY,NX)
        ENGYV=(cpw*(VLWatMicP_vr(L,NY,NX)+VLWatMacP_vr(L,NY,NX))+cpi*(VLiceMicP_vr(L,NY,NX)+VLiceMacP_vr(L,NY,NX)))*TKS_vr(L,NY,NX)
        VLWatMicP_vr(L,NY,NX)=TI*VLWatMicP_vr(L,NY,NX)+CORP*(FI*TVOLW-TI*VLWatMicP_vr(L,NY,NX))+TX*VLWatMicP_vr(L,NY,NX)+FI*TVOLWR
        VLiceMicP_vr(L,NY,NX)=TI*VLiceMicP_vr(L,NY,NX)+CORP*(FI*TVOLI-TI*VLiceMicP_vr(L,NY,NX))+TX*VLiceMicP_vr(L,NY,NX)
        VLWatMicPX_vr(L,NY,NX)=VLWatMicP_vr(L,NY,NX)
!     VLWatMicP_vr(L,NY,NX)=VLWatMicP_vr(L,NY,NX)+CORP*VLWatMacP_vr(L,NY,NX)
!     VLiceMicP_vr(L,NY,NX)=VLiceMicP_vr(L,NY,NX)+CORP*VLiceMacP_vr(L,NY,NX)
!     VLMicP_vr(L,NY,NX)=VLMicP_vr(L,NY,NX)+CORP*VLMacP_vr(L,NY,NX)
!     VLWatMacP_vr(L,NY,NX)=XCORP(NY,NX)*VLWatMacP_vr(L,NY,NX)
!     VLiceMacP_vr(L,NY,NX)=XCORP(NY,NX)*VLiceMacP_vr(L,NY,NX)
!     VLMacP_vr(L,NY,NX)=XCORP(NY,NX)*VLMacP_vr(L,NY,NX)
!     SoilFracAsMacP_vr(L,NY,NX)=XCORP(NY,NX)*SoilFracAsMacP_vr(L,NY,NX)
        ENGYL=TI*ENGYV+CORP*(FI*TENGY-TI*ENGYV)+TX*ENGYV+FI*TENGYR
        VHeatCapacity_vr(L,NY,NX)=VHeatCapacitySoilM(L,NY,NX)+cpw*(VLWatMicP_vr(L,NY,NX)+VLWatMacP_vr(L,NY,NX)) &
          +cpi*(VLiceMicP_vr(L,NY,NX)+VLiceMacP_vr(L,NY,NX))
        TKS_vr(L,NY,NX)=(ENGYM+ENGYL)/VHeatCapacity_vr(L,NY,NX)
        TCS(L,NY,NX)=units%Kelvin2Celcius(TKS_vr(L,NY,NX))

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
          trcSalt_solml_vr(NTSA,L,NY,NX)=TI*trcSalt_solml_vr(NTSA,L,NY,NX)  &
            +CORP*(FI*TSA_solml(NTSA)-TI*trcSalt_solml_vr(NTSA,L,NY,NX)) &
            +TX*trcSalt_solml_vr(NTSA,L,NY,NX)+CORP*trcSalt_soHml_vr(NTSA,L,NY,NX)
        ENDDO

        ! solute
        DO NTS=ids_beg,ids_end
          trc_solml_vr(NTS,L,NY,NX)=TI*trc_solml_vr(NTS,L,NY,NX) &
            +CORP*(FI*TS_solml(NTS)-TI*trc_solml_vr(NTS,L,NY,NX)) &
            +TX*trc_solml_vr(NTS,L,NY,NX)+CORP*trc_soHml_vr(NTS,L,NY,NX)
        ENDDO

        DO NTX=idx_CEC,idx_cation_end
          trcx_solml_vr(NTX,L,NY,NX)=TI*trcx_solml_vr(NTX,L,NY,NX)+ &
            CORP*(FI*TX_solml(NTX)-TI*trcx_solml_vr(NTX,L,NY,NX))+TX*trcx_solml_vr(NTX,L,NY,NX)
        ENDDO

        DO NTX=idx_AEC,idx_end
          trcx_solml_vr(NTX,L,NY,NX)=TI*trcx_solml_vr(NTX,L,NY,NX)+ &
            CORP*(FI*TX_solml(NTX)-TI*trcx_solml_vr(NTX,L,NY,NX))+TX*trcx_solml_vr(NTX,L,NY,NX)
        ENDDO

        DO NTP=idsp_beg,idsp_end
          trcp_saltpml_vr(NTP,L,NY,NX)=TI*trcp_saltpml_vr(NTP,L,NY,NX)+CORP*(FI*TP_salml(NTP) &
            -TI*trcp_saltpml_vr(NTP,L,NY,NX))+TX*trcp_saltpml_vr(NTP,L,NY,NX)
        ENDDO

        DO NTG=idg_beg,idg_end-1
          trc_gasml_vr(NTG,L,NY,NX)=TI*trc_gasml_vr(NTG,L,NY,NX)+CORP*(FI*TG_gasml(NTG) &
            -TI*trc_gasml_vr(NTG,L,NY,NX))+TX*trc_gasml_vr(NTG,L,NY,NX)
        ENDDO

        call MixSoluteMacpore(L,NY,NX)

        DO  K=1,jcplx
          DO  N=1,NumMicbFunGrupsPerCmplx
            DO NGL=JGnio(N),JGnfo(N)
              DO  M=1,nlbiomcp
                MID=micpar%get_micb_id(M,NGL)
                DO NE=1,NumPlantChemElms
                  mBiomeHeter_vr(NE,MID,K,L,NY,NX)=TI*mBiomeHeter_vr(NE,MID,K,L,NY,NX)+CORP*(FI*TmBiomeHeter(NE,MID,K) &
                    -TI*mBiomeHeter_vr(NE,MID,K,L,NY,NX))+TX*mBiomeHeter_vr(NE,MID,K,L,NY,NX)
                  ENDDO
              enddo
            enddo
          enddo
        ENDDO
        DO  N=1,NumMicbFunGrupsPerCmplx
          DO NGL=JGniA(N),JGnfA(N)
            DO  M=1,nlbiomcp
              MID=micpar%get_micb_id(M,NGL)
              DO NE=1,NumPlantChemElms
                mBiomeAutor_vr(NE,MID,L,NY,NX)=TI*mBiomeAutor_vr(NE,MID,L,NY,NX)+CORP*(FI*TmBiomeAutor(NE,MID) &
                  -TI*mBiomeAutor_vr(NE,MID,L,NY,NX))+TX*mBiomeAutor_vr(NE,MID,L,NY,NX)
              ENDDO  
            enddo
          enddo
        enddo

        DO  K=1,jcplx
          DO  M=1,ndbiomcp
            DO NE=1,NumPlantChemElms
              OMBioResdu_vr(NE,M,K,L,NY,NX)=TI*OMBioResdu_vr(NE,M,K,L,NY,NX)+CORP*(FI*TORM(NE,M,K) &
                -TI*OMBioResdu_vr(NE,M,K,L,NY,NX))+TX*OMBioResdu_vr(NE,M,K,L,NY,NX)
            ENDDO  
          ENDDO

          DO idom=idom_beg,idom_end
            DOM_vr(idom,K,L,NY,NX)=TI*DOM_vr(idom,K,L,NY,NX)+CORP*(FI*TDOM(idom,K) &
              -TI*DOM_vr(idom,K,L,NY,NX))+TX*DOM_vr(idom,K,L,NY,NX)+CORP*DOM_MacP_vr(idom,K,L,NY,NX)
            DOM_MacP_vr(idom,K,L,NY,NX)=XCORP(NY,NX)*DOM_MacP_vr(idom,K,L,NY,NX)
            SorbedOM_vr(idom,K,L,NY,NX)=TI*SorbedOM_vr(idom,K,L,NY,NX)+CORP*(FI*TOHM(idom,K) &
              -TI*SorbedOM_vr(idom,K,L,NY,NX))+TX*SorbedOM_vr(idom,K,L,NY,NX)
          ENDDO  

          DO  M=1,jsken
            SolidOMAct_vr(M,K,L,NY,NX)=TI*SolidOMAct_vr(M,K,L,NY,NX)+CORP*(FI*TOSA(M,K) &
              -TI*SolidOMAct_vr(M,K,L,NY,NX))+TX*SolidOMAct_vr(M,K,L,NY,NX)
            DO NE=1,NumPlantChemElms  
              SolidOM_vr(NE,M,K,L,NY,NX)=TI*SolidOM_vr(NE,M,K,L,NY,NX)+CORP*(FI*TOSM(NE,M,K) &
                -TI*SolidOM_vr(NE,M,K,L,NY,NX))+TX*SolidOM_vr(NE,M,K,L,NY,NX)
            ENDDO

          ENDDO
        ENDDO
!
!     ADD STATE VARIABLES IN SURFACE RESIDUE INCORPORATED
!     WITHIN TILLAGE MIXING ZONE
!
        DO  K=1,micpar%NumOfLitrCmplxs
          DO  N=1,NumMicbFunGrupsPerCmplx
            DO NGL=JGnio(N),JGnfo(N)
              DO M=1,nlbiomcp
                MID=micpar%get_micb_id(M,NGL)         
                DO NE=1,NumPlantChemElms       
                  mBiomeHeter_vr(NE,MID,K,L,NY,NX)=mBiomeHeter_vr(NE,MID,K,L,NY,NX)+FI*TBOMGE(NE,MID,K)
                ENDDO
              enddo
            enddo
          ENDDO
        ENDDO

        DO  N=1,NumMicbFunGrupsPerCmplx
          DO NGL=JGniA(N),JGnfA(N)
            DO M=1,nlbiomcp
              MID=micpar%get_micb_id(M,NGL)
              DO NE=1,NumPlantChemElms
                mBiomeAutor_vr(NE,MID,L,NY,NX)=mBiomeAutor_vr(NE,MID,L,NY,NX)+FI*TOMGAutor(NE,MID)
              ENDDO
            enddo
          enddo
        ENDDO

        DO K=1,micpar%NumOfLitrCmplxs
          DO  M=1,ndbiomcp
            DO NE=1,NumPlantChemElms
              OMBioResdu_vr(NE,M,K,L,NY,NX)=OMBioResdu_vr(NE,M,K,L,NY,NX)+FI*TORXE(NE,M,K)
            ENDDO
          ENDDO
          DO idom=idom_beg,idom_end
            DOM_vr(idom,K,L,NY,NX)=DOM_vr(idom,K,L,NY,NX)+FI*TDOMG(idom,K)
            DOM_MacP_vr(idom,K,L,NY,NX)=DOM_MacP_vr(idom,K,L,NY,NX)+FI*TDOMH(idom,K)
            SorbedOM_vr(idom,K,L,NY,NX)=SorbedOM_vr(idom,K,L,NY,NX)+FI*TDOMHG(idom,K)
          ENDDO
          DO  M=1,jsken
            SolidOMAct_vr(M,K,L,NY,NX)=SolidOMAct_vr(M,K,L,NY,NX)+FI*TOSGA(M,K)
            DO NE=1,NumPlantChemElms
              SolidOM_vr(NE,M,K,L,NY,NX)=SolidOM_vr(NE,M,K,L,NY,NX)+FI*TOSGE(NE,M,K)
            ENDDO
          ENDDO
        ENDDO

        call sumORGMLayL(L,NY,NX,ORGM)

        call sumLitrOMLayL(L,NY,NX,litrOM)

        SoilOrgM_vr(ielmc,L,NY,NX)=ORGM(ielmc)
        SoilOrgM_vr(ielmn,L,NY,NX)=ORGM(ielmn)
        OMLitrC_vr(L,NY,NX)=litrOM(ielmc)

        DO NTS=ids_beg,idg_NH3
          trc_solml_vr(NTS,L,NY,NX)=trc_solml_vr(NTS,L,NY,NX)+FI*TS0_solml(NTS)
        ENDDO

        DO NTS=ids_nut_beg,ids_nuts_end
          trc_solml_vr(NTS,L,NY,NX)=trc_solml_vr(NTS,L,NY,NX)+FI*TS0_solml(NTS)
        ENDDO

        trcx_solml_vr(idx_NH4,L,NY,NX)=trcx_solml_vr(idx_NH4,L,NY,NX)+FI*TXN4G


        DO NTX=idx_AEC+1,idx_anion_soil_end
          trcx_solml_vr(NTX,L,NY,NX)=trcx_solml_vr(NTX,L,NY,NX)+FI*TX_anion(NTX)
        ENDDO

        DO NTP=idsp_psoi_beg,idsp_psoi_end
          trcp_saltpml_vr(NTP,L,NY,NX)=trcp_saltpml_vr(NTP,L,NY,NX)+FI*TP_soil(NTP)
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

  subroutine MixSoluteMacpore(L,NY,NX)
  implicit none
  integer, intent(in) :: L,NY,NX
  integer :: NTS,NTSA

  DO NTS=ids_beg,ids_end
    trc_soHml_vr(NTS,L,NY,NX)=XCORP(NY,NX)*trc_soHml_vr(NTS,L,NY,NX)
  ENDDO

  DO NTSA=idsalt_beg,idsalt_end
    trcSalt_soHml_vr(ntsa,L,NY,NX)=XCORP(NY,NX)*trcSalt_soHml_vr(ntsa,L,NY,NX)
  ENDDO
  end subroutine MixSoluteMacpore


end module TillageMixMod
