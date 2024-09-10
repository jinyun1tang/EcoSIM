module TillageMixMod
  use data_kind_mod,    only : r8 => DAT_KIND_R8
  use EcoSiMParDataMod, only : micpar
  use EcoSIMConfig,     only : ndbiomcp => NumDeadMicrbCompts
  use UnitMod,          only : units
  use minimathmod,      only : AZMAX1
  use SoilBGCNLayMod,   only : sumLitrOMLayL, sumORGMLayL
  use EcoSIMCtrlMod,    only : salt_model
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
  real(r8) :: FracTillCorp,CORP0,XCORP0,DCORPZ
  real(r8) :: TWP,TSatHydroCondVert
  REAL(R8) :: ZNHUX0,ZNHUXI,ZNFNX0
  real(r8) :: TGKCK
  real(r8) :: TP_soil(idsp_psoi_beg:idsp_psoi_end)
  real(r8) :: TFertNG_soil(ifertn_beg:ifertn_end),TZNFNG,TVOLWR
  real(r8) :: TZ2OG,TZ2OS,TZNH3G,TH2GG,TH2GS,TZNFN2,TZNFNI,TCO2GS
  real(r8) :: TCH4GS,TOXYGS,TZ2GSG,TZ2OGS,TH2GGS,TNH4GS,TNH3GS
  real(r8) :: TX_anion(idx_AEC+1:idx_anion_soil_end)
  real(r8) :: TNO3GS,TNO2GS,TP14GS,TPO4GS,TXN4G
  real(r8) :: TX_solml(idx_beg:idx_end)
  real(r8) :: TP_salml(idsp_beg:idsp_end)
  real(r8) :: TG_gasml(idg_beg:idg_end-1)
  real(r8) :: TS0_solml(ids_beg:ids_end)   !surface mass for incoporation
  real(r8) :: TfertN_band(ifertn_beg:ifertnb_end)
  real(r8) :: TNH4S,TNH4B,TNH3S,TNH3B,TNO3S
  real(r8) :: TNO3B,TNO2S,TNO2B,TZAL,TZFE,TZHY,TZCA,TZMG,TZNA
  real(r8) :: TNFNIH
  real(r8) :: TENGYR,TVOLW,TENGY
  real(r8) :: ENGYV,ENGYL,ENGYM
  real(r8) :: ORGM(1:NumPlantChemElms),litrOM(1:NumPlantChemElms)
  real(r8) :: TVOLI,HFLXD
  real(r8) :: TX,TL
  integer  :: NTX,NTP,NTG,NTSA,NTN,idom,MID,NE
  real(r8) :: FI(JZ),TI(JZ),TI1,FI1
  logical :: TillTest,tillTestL

!     begin_execution
!
  !DCORP=soil mixing fraction with tillage
  !XCORP=factor for surface litter incorporation and soil mixing

  TillTest=J.EQ.INT(SolarNoonHour_col(NY,NX)) .AND. XTillCorp_col(NY,NX).LT.1.0_r8 .AND. DepzCorp_col(I,NY,NX).GT.0.0_r8
  if(.not. TillTest)return
!
  iResetSoilProf_col(NY,NX) = itrue
!   tillage depth cannot exceed whole soil column thickness
  DCORPZ = AMIN1(DepzCorp_col(I,NY,NX),CumSoilThickness_vr(NL(NY,NX),NY,NX))

!
!     ACCUMULATE STATE VARIABLES IN SURFACE RESIDUE FOR ADDITION
!     TO SOIL IN TILLAGE MIXING ZONE
!
!     IF(SoilOrgM_vr(ielmc,0,NY,NX).GT.ZEROS(NY,NX))THEN
!     XCORP0=AMAX1(XTillCorp_col(NY,NX),AMIN1(1.0,
!    2(VHeatCapLitRMin_col(NY,NX)/cpo)/SoilOrgM_vr(ielmc,0,NY,NX)))
  XCORP0=AMAX1(0.001_r8,XTillCorp_col(NY,NX))
!     ELSE
!     XCORP0=1.0
!     ENDIF
!   XCORP0: Fraction of material not being mixed

  call DeriveTillageProf(NY,NX,DCORPZ,LL,FI,TI)

  call Mix1D(XTillCorp_col(NY,NX),TI,FI,NU(NY,NX),LL,SoiBulkDensityt0_vr(1:JZ,NY,NX))

  call Mix1D(XTillCorp_col(NY,NX),TI,FI,NU(NY,NX),LL,FieldCapacity_vr(1:JZ,NY,NX))

  call Mix1D(XTillCorp_col(NY,NX),TI,FI,NU(NY,NX),LL,WiltPoint_vr(1:JZ,NY,NX))

  call Mix1D(XTillCorp_col(NY,NX),TI,FI,NU(NY,NX),LL,SatHydroCondVert_vr(1:JZ,NY,NX))

  call Mix1D(XTillCorp_col(NY,NX),TI,FI,NU(NY,NX),LL,SatHydroCondHrzn_vr(1:JZ,NY,NX))

  call Mix1D(XTillCorp_col(NY,NX),TI,FI,NU(NY,NX),LL,SAND(1:JZ,NY,NX))

  call Mix1D(XTillCorp_col(NY,NX),TI,FI,NU(NY,NX),LL,SILT(1:JZ,NY,NX))

  call Mix1D(XTillCorp_col(NY,NX),TI,FI,NU(NY,NX),LL,CLAY(1:JZ,NY,NX))

  call Mix1D(XTillCorp_col(NY,NX),TI,FI,NU(NY,NX),LL,GKC4(1:JZ,NY,NX))

  call Mix1D(XTillCorp_col(NY,NX),TI,FI,NU(NY,NX),LL,GKCA(1:JZ,NY,NX))

  call Mix1D(XTillCorp_col(NY,NX),TI,FI,NU(NY,NX),LL,GKCM(1:JZ,NY,NX))

  call Mix1D(XTillCorp_col(NY,NX),TI,FI,NU(NY,NX),LL,GKCN(1:JZ,NY,NX))

  call Mix1D(XTillCorp_col(NY,NX),TI,FI,NU(NY,NX),LL,GKCK(1:JZ,NY,NX))

  CALL Mix2D(XTillCorp_col(NY,NX),TI,FI,NU(NY,NX),LL,FertN_soil_vr(ifertn_beg:ifertn_end,0:JZ,NY,NX),XCORP0)

  CALL Mix2D(XTillCorp_col(NY,NX),TI,FI,NU(NY,NX),LL,trc_solml_vr(ids_beg:idg_NH3,0:JZ,NY,NX),XCORP0)

  CALL Mix2D(XTillCorp_col(NY,NX),TI,FI,NU(NY,NX),LL,trc_solml_vr(ids_nut_beg:ids_nuts_end,0:JZ,NY,NX),XCORP0)

  CALL Mix2D(XTillCorp_col(NY,NX),TI,FI,NU(NY,NX),LL,trc_solml_vr(idg_NH3B:ids_nutb_end,1:JZ,NY,NX))

  if(salt_model)then
    CALL Mix2D(XTillCorp_col(NY,NX),TI,FI,NU(NY,NX),LL, trcSalt_solml_vr(idsalt_beg:idsalt_end,0:JZ,NY,NX),XCORP0)

    CALL Mix2D(XTillCorp_col(NY,NX),TI,FI,NU(NY,NX),LL, trcSalt_solml_vr(idsalt_pband_beg:idsalt_pband_end,1:JZ,NY,NX))
  endif

  CALL Mix2D(XTillCorp_col(NY,NX),TI,FI,NU(NY,NX),LL, trcp_saltpml_vr(idsp_beg:idsp_CaSO4,1:JZ,NY,NX))

  CALL Mix2D(XTillCorp_col(NY,NX),TI,FI,NU(NY,NX),LL, trcp_saltpml_vr(idsp_beg_band:idsp_end,1:JZ,NY,NX))

  CALL Mix2D(XTillCorp_col(NY,NX),TI,FI,NU(NY,NX),LL, trcp_saltpml_vr(idsp_psoi_beg:idsp_psoi_end,0:JZ,NY,NX),XCORP0)
!
!     EXTENT OF MIXING

    !
    FracTillCorp                = 1.0_r8-XTillCorp_col(NY,NX)
    EnergyImpact4Erosion(NY,NX) = 0.0_r8
!
!     TEMPORARY ACCUMULATORS
!
    TVOLW             = 0.0_r8
    TVOLI             = 0.0_r8
  ! TVOLP             = 0.0_r8
  ! TVOLA             = 0.0_r8
    TENGY             = 0.0_r8
    TNFNIH            = 0.0_r8

    TfertN_band(ifertnb_beg:ifertnb_end)=0.0_r8

    TS0_solml(ids_beg:ids_end)  = 0._r8
    TX_solml(idx_beg:idx_end)   = 0._r8
    TP_salml(idsp_beg:idsp_end) = 0._r8
    TG_gasml(idg_beg:idg_end-1) = 0._r8

    TmBiomeHeter = 0.0_r8
    TmBiomeAutor = 0.0_r8
    TORM         = 0.0_r8
    TDOM         = 0.0_r8
    TOHM         = 0.0_r8
    TOSM         = 0.0_r8
    TOSA         = 0.0_r8
    TZNFN2       = 0.0_r8
    TZNFNI       = 0.0_r8
    ZNHUX0       = 0.0_r8
    ZNHUXI       = 0.0_r8
    ZNFNX0       = 0.0_r8

!   CORP0: Fraction of material being mixed
    CORP0=1.0_r8-XCORP0   

    DO  K=1,micpar%NumOfLitrCmplxs
      DO  N=1,NumMicbFunGrupsPerCmplx
        DO NGL=JGnio(N),JGnfo(N)
          DO  M=1,nlbiomcp
            MID=micpar%get_micb_id(M,NGL)
            DO NE=1,NumPlantChemElms
              TBOMGE(NE,MID,K)                 = mBiomeHeter_vr(NE,MID,K,0,NY,NX)*CORP0
              mBiomeHeter_vr(NE,MID,K,0,NY,NX) = mBiomeHeter_vr(NE,MID,K,0,NY,NX)*XCORP0
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
            TOMGAutor(NE,MID)              = mBiomeAutor_vr(NE,MID,0,NY,NX)*CORP0
            mBiomeAutor_vr(NE,MID,0,NY,NX) = mBiomeAutor_vr(NE,MID,0,NY,NX)*XCORP0
          ENDDO
        enddo
      enddo
    ENDDO

    DO K=1,micpar%NumOfLitrCmplxs
      DO M=1,ndbiomcp
        DO NE=1,NumPlantChemElms
          TORXE(NE,M,K)                 = OMBioResdu_vr(NE,M,K,0,NY,NX)*CORP0
          OMBioResdu_vr(NE,M,K,0,NY,NX) = OMBioResdu_vr(NE,M,K,0,NY,NX)*XCORP0
        ENDDO
      ENDDO

      DO idom=idom_beg,idom_end
        TDOMG(idom,K)  = DOM_vr(idom,K,0,NY,NX)*CORP0
        TDOMH(idom,K)  = DOM_MacP_vr(idom,K,0,NY,NX)*CORP0
        TDOMHG(idom,K) = SorbedOM_vr(idom,K,0,NY,NX)*CORP0
  !
  !     REDUCE SURFACE RESIDUE STATE VARIABLES FOR INCORPORATION
  !
        DOM_vr(idom,K,0,NY,NX)      = DOM_vr(idom,K,0,NY,NX)*XCORP0
        DOM_MacP_vr(idom,K,0,NY,NX) = DOM_MacP_vr(idom,K,0,NY,NX)*XCORP0
        SorbedOM_vr(idom,K,0,NY,NX) = SorbedOM_vr(idom,K,0,NY,NX)*XCORP0
      enddo 

      DO  M=1,jsken
        TOSGA(M,K)                 = SolidOMAct_vr(M,K,0,NY,NX)*CORP0
        SolidOMAct_vr(M,K,0,NY,NX) = SolidOMAct_vr(M,K,0,NY,NX)*XCORP0

        DO NE=1,NumPlantChemElms
          TOSGE(NE,M,K)              = SolidOM_vr(NE,M,K,0,NY,NX)*CORP0
          SolidOM_vr(NE,M,K,0,NY,NX) = SolidOM_vr(NE,M,K,0,NY,NX)*XCORP0
        ENDDO

      ENDDO
    ENDDO

    TXN4G=trcx_solml_vr(idx_NH4,0,NY,NX)*CORP0

    DO NTX= idx_AEC+1,idx_anion_soil_end
      TX_anion(NTX)=trcx_solml_vr(NTX,0,NY,NX)*CORP0
    ENDDO


    TZNFNG        = ZNFNI(0,NY,NX)*CORP0
    TVOLWR        = VLWatMicP_vr(0,NY,NX)*CORP0
    HFLXD         = cpo*SoilOrgM_vr(ielmc,0,NY,NX)*CORP0*TKS_vr(0,NY,NX)
    HEATIN_lnd    = HEATIN_lnd-HFLXD
    HeatStore_lnd = HeatStore_lnd-HFLXD
    TENGYR        = cpw*TVOLWR*TKS_vr(0,NY,NX)

    call sumLitrOMLayL(0,NY,NX,litrOM)

    SoilOrgM_vr(1:NumPlantChemElms,0,NY,NX)=litrOM(1:NumPlantChemElms)

    OMLitrC_vr(0,NY,NX)=litrOM(ielmc)

    trcx_solml_vr(idx_NH4,0,NY,NX)=trcx_solml_vr(idx_NH4,0,NY,NX)*XCORP0

    DO NTX=idx_AEC+1,idx_anion_soil_end
      trcx_solml_vr(NTX,0,NY,NX)=trcx_solml_vr(NTX,0,NY,NX)*XCORP0
    ENDDO

    DO NTN=ifertn_beg,ifertn_end
      FertN_soil_vr(NTN,0,NY,NX)=FertN_soil_vr(NTN,0,NY,NX)*XCORP0
    ENDDO

    VLWatMicP_vr(0,NY,NX)     = VLWatMicP_vr(0,NY,NX)*XCORP0
    VHeatCapacity_vr(0,NY,NX) = cpo*SoilOrgM_vr(ielmc,0,NY,NX)+cpw*VLWatMicP_vr(0,NY,NX)+cpi*VLiceMicP_vr(0,NY,NX)
    VLitR_col(NY,NX)          = VLitR_col(NY,NX)*XCORP0
    VGeomLayer_vr(0,NY,NX)    = VGeomLayer_vr(0,NY,NX)*XCORP0
    ZNHUX0                    = AMAX1(ZNHUX0,ZNHU0(0,NY,NX))
    ZNHUXI                    = AMAX1(ZNHUXI,ZNHUI(0,NY,NX))
    ZNFNX0                    = AMAX1(ZNFNX0,ZNFN0(0,NY,NX))
!
!     REDISTRIBUTE SOIL STATE VARIABLES DURING TILLAGE
!
!     ACCUMULATE SOIL STATE VARIABLES WITHIN TILLAGE MIXING ZONE
!
    D100: DO L=NU(NY,NX),NL(NY,NX)
      !if current layer is within the tillage depth
      tillTestL=CumSoilThickness_vr(L,NY,NX)-DLYR(3,L,NY,NX).LT.DCORPZ .AND. DLYR(3,L,NY,NX).GT.ZERO
      if(.not.tillTestL)exit
!     TL: fraction of layer L over the tillage zone
!     TI: fraction of layer L being tilled
      TL = AMIN1(DLYR(3,L,NY,NX),DCORPZ-(CumSoilThickness_vr(L,NY,NX)-DLYR(3,L,NY,NX)))
      FI1 = TL/DCORPZ
      TI1 = AMIN1(TL/DLYR(3,L,NY,NX),1._r8)

      TVOLW = TVOLW+TI1*VLWatMicP_vr(L,NY,NX)
      TVOLI = TVOLI+TI1*VLiceMicP_vr(L,NY,NX)
! TVOLP = TVOLP+TI1*VLsoiAirP_vr(L,NY,NX)
! TVOLA = TVOLA+TI1*VLMicP_vr(L,NY,NX)
      TENGY = TENGY+TI1*(cpw*(VLWatMicP_vr(L,NY,NX)+VLWatMacP_vr(L,NY,NX)) &
        +cpi*(VLiceMicP_vr(L,NY,NX)+VLiceMacP_vr(L,NY,NX)))*TKS_vr(L,NY,NX)


      DO NTN=ifertnb_beg,ifertnb_end
        TfertN_band(NTN)=TfertN_band(NTN)+TI1*FertN_Band_vr(NTN,L,NY,NX)
      ENDDO


!cation
      DO NTX=idx_CEC,idx_cation_end
        TX_solml(NTX)=TX_solml(NTX)+TI1*trcx_solml_vr(NTX,L,NY,NX)
      ENDDO
!anion
      DO NTX=idx_AEC,idx_end
        TX_solml(NTX)=TX_solml(NTX)+TI1*trcx_solml_vr(NTX,L,NY,NX)
      ENDDO

      DO NTG=idg_beg,idg_end-1
        TG_gasml(NTG)=TG_gasml(NTG)+TI1*trc_gasml_vr(NTG,L,NY,NX)
      ENDDO

      DO  K=1,jcplx
        DO  N=1,NumMicbFunGrupsPerCmplx
          DO NGL=JGnio(N),JGnfo(N)
            DO  M=1,nlbiomcp
              MID=micpar%get_micb_id(M,NGL)
              DO NE=1,NumPlantChemElms
                TmBiomeHeter(NE,MID,K)=TmBiomeHeter(NE,MID,K)+TI1*mBiomeHeter_vr(NE,MID,K,L,NY,NX)
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
            TmBiomeAutor(NE,MID)=TmBiomeAutor(NE,MID)+TI1*mBiomeAutor_vr(NE,MID,L,NY,NX)
            ENDDO
          enddo
        enddo
      enddo

      DO K=1,jcplx
        DO  M=1,ndbiomcp
          DO NE=1,NumPlantChemElms
            TORM(NE,M,K)=TORM(NE,M,K)+TI1*OMBioResdu_vr(NE,M,K,L,NY,NX)
          ENDDO
        ENDDO
        do idom=idom_beg,idom_end
          TDOM(idom,K)=TDOM(idom,K)+TI1*DOM_vr(idom,K,L,NY,NX)
          TOHM(idom,K)=TOHM(idom,K)+TI1*SorbedOM_vr(idom,K,L,NY,NX)
        enddo
        DO  M=1,jsken
          TOSA(M,K)=TOSA(M,K)+TI1*SolidOMAct_vr(M,K,L,NY,NX)
          DO NE=1,NumPlantChemElms
            TOSM(NE,M,K)=TOSM(NE,M,K)+TI1*SolidOM_vr(NE,M,K,L,NY,NX)
          ENDDO
        ENDDO
      ENDDO
      ZNHUX0=AMAX1(ZNHUX0,ZNHU0(L,NY,NX))
      ZNHUXI=AMAX1(ZNHUXI,ZNHUI(L,NY,NX))
      ZNFNX0=AMAX1(ZNFNX0,ZNFN0(L,NY,NX))
      TZNFNI=TZNFNI+ZNFNI(L,NY,NX)
      LL=L      
    ENDDO D100
!
!     CHANGE SOIL STATE VARIABLES IN TILLAGE MIXING ZONE
!     TO ACCOUNT FOR REDISTRIBUTION FROM MIXING
!   but non-P related precipitated species are not mixed, Jinyun Tang, Nov 30,2022
!   LL is the last layer tilled. 

    D2000: DO  L=NU(NY,NX),LL
      IF(DLYR(3,L,NY,NX).GT.ZERO)THEN
        !TL: soil thickness in tillage zone
        !TI: fraction of current layer being tilled.
        !FI: fraction of current layer over the overall mixing zone
        !TX: fraction of current layer not tilled
        TL=AMIN1(DLYR(3,L,NY,NX),DCORPZ-(CumSoilThickness_vr(L,NY,NX)-DLYR(3,L,NY,NX)))
        FI1=TL/DCORPZ
        TI1=AMIN1(TL/DLYR(3,L,NY,NX),1._r8)
        TX=AZMAX1(1.0_r8-TI1)
        
        ENGYM                  = VHeatCapacitySoilM_vr(L,NY,NX)*TKS_vr(L,NY,NX)
        ENGYV                  = (cpw*(VLWatMicP_vr(L,NY,NX)+VLWatMacP_vr(L,NY,NX)) &
          +cpi*(VLiceMicP_vr(L,NY,NX)+VLiceMacP_vr(L,NY,NX)))*TKS_vr(L,NY,NX)
        VLWatMicP_vr(L,NY,NX)  = TI1*VLWatMicP_vr(L,NY,NX)+FracTillCorp*(FI1*TVOLW-TI1*VLWatMicP_vr(L,NY,NX)) &
          +TX*VLWatMicP_vr(L,NY,NX)+FI1*TVOLWR
        VLiceMicP_vr(L,NY,NX)  = TI1*VLiceMicP_vr(L,NY,NX)+FracTillCorp*(FI1*TVOLI-TI1*VLiceMicP_vr(L,NY,NX)) &
          +TX*VLiceMicP_vr(L,NY,NX)
        VLWatMicPX_vr(L,NY,NX) = VLWatMicP_vr(L,NY,NX)

!     VLWatMicP_vr(L,NY,NX)=VLWatMicP_vr(L,NY,NX)+FracTillCorp*VLWatMacP_vr(L,NY,NX)
!     VLiceMicP_vr(L,NY,NX)=VLiceMicP_vr(L,NY,NX)+FracTillCorp*VLiceMacP_vr(L,NY,NX)
!     VLMicP_vr(L,NY,NX)=VLMicP_vr(L,NY,NX)+FracTillCorp*VLMacP_vr(L,NY,NX)
!     VLWatMacP_vr(L,NY,NX)=XTillCorp_col(NY,NX)*VLWatMacP_vr(L,NY,NX)
!     VLiceMacP_vr(L,NY,NX)=XTillCorp_col(NY,NX)*VLiceMacP_vr(L,NY,NX)
!     VLMacP_vr(L,NY,NX)=XTillCorp_col(NY,NX)*VLMacP_vr(L,NY,NX)
!     SoilFracAsMacP_vr(L,NY,NX)=XTillCorp_col(NY,NX)*SoilFracAsMacP_vr(L,NY,NX)
        ENGYL=TI1*ENGYV+FracTillCorp*(FI1*TENGY-TI1*ENGYV)+TX*ENGYV+FI1*TENGYR

        VHeatCapacity_vr(L,NY,NX)=VHeatCapacitySoilM_vr(L,NY,NX)+cpw*(VLWatMicP_vr(L,NY,NX)+VLWatMacP_vr(L,NY,NX)) &
          +cpi*(VLiceMicP_vr(L,NY,NX)+VLiceMacP_vr(L,NY,NX))
        TKS_vr(L,NY,NX)=(ENGYM+ENGYL)/VHeatCapacity_vr(L,NY,NX)        
        TCS(L,NY,NX)=units%Kelvin2Celcius(TKS_vr(L,NY,NX))


        DO NTN=ifertnb_beg,ifertnb_end
          FertN_Band_vr(NTN,L,NY,NX)=TI1*FertN_Band_vr(NTN,L,NY,NX) &
            +FracTillCorp*(FI1*TFertN_Band(NTN)-TI1*FertN_Band_vr(NTN,L,NY,NX)) &
            +TX*FertN_Band_vr(NTN,L,NY,NX)
        ENDDO

        !SALT
        DO NTSA=idsalt_beg,idsalt_end
          trcSalt_solml_vr(NTSA,L,NY,NX)=trcSalt_solml_vr(NTSA,L,NY,NX)  &
            +FracTillCorp*trcSalt_soHml_vr(NTSA,L,NY,NX)
        ENDDO

        ! solute
        DO NTS=ids_beg,ids_end
          trc_solml_vr(NTS,L,NY,NX)=trc_solml_vr(NTS,L,NY,NX)+FracTillCorp*trc_soHml_vr(NTS,L,NY,NX)
        ENDDO

        DO NTX=idx_CEC,idx_cation_end
          trcx_solml_vr(NTX,L,NY,NX)=TI1*trcx_solml_vr(NTX,L,NY,NX)+ &
            FracTillCorp*(FI1*TX_solml(NTX)-TI1*trcx_solml_vr(NTX,L,NY,NX))+TX*trcx_solml_vr(NTX,L,NY,NX)
        ENDDO

        DO NTX=idx_AEC,idx_end
          trcx_solml_vr(NTX,L,NY,NX)=TI1*trcx_solml_vr(NTX,L,NY,NX)+ &
            FracTillCorp*(FI1*TX_solml(NTX)-TI1*trcx_solml_vr(NTX,L,NY,NX))+TX*trcx_solml_vr(NTX,L,NY,NX)
        ENDDO


        DO NTG=idg_beg,idg_end-1
          trc_gasml_vr(NTG,L,NY,NX)=TI1*trc_gasml_vr(NTG,L,NY,NX)+FracTillCorp*(FI1*TG_gasml(NTG) &
            -TI1*trc_gasml_vr(NTG,L,NY,NX))+TX*trc_gasml_vr(NTG,L,NY,NX)
        ENDDO

        call MixSoluteMacpore(L,NY,NX)

        DO  K=1,jcplx
          DO  N=1,NumMicbFunGrupsPerCmplx
            DO NGL=JGnio(N),JGnfo(N)
              DO  M=1,nlbiomcp
                MID=micpar%get_micb_id(M,NGL)
                DO NE=1,NumPlantChemElms
                  mBiomeHeter_vr(NE,MID,K,L,NY,NX)=TI1*mBiomeHeter_vr(NE,MID,K,L,NY,NX)+FracTillCorp*(FI1*TmBiomeHeter(NE,MID,K) &
                    -TI1*mBiomeHeter_vr(NE,MID,K,L,NY,NX))+TX*mBiomeHeter_vr(NE,MID,K,L,NY,NX)
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
                mBiomeAutor_vr(NE,MID,L,NY,NX)=TI1*mBiomeAutor_vr(NE,MID,L,NY,NX)+FracTillCorp*(FI1*TmBiomeAutor(NE,MID) &
                  -TI1*mBiomeAutor_vr(NE,MID,L,NY,NX))+TX*mBiomeAutor_vr(NE,MID,L,NY,NX)
              ENDDO  
            enddo
          enddo
        enddo

        DO  K=1,jcplx
          DO  M=1,ndbiomcp
            DO NE=1,NumPlantChemElms
              OMBioResdu_vr(NE,M,K,L,NY,NX)=TI1*OMBioResdu_vr(NE,M,K,L,NY,NX)+FracTillCorp*(FI1*TORM(NE,M,K) &
                -TI1*OMBioResdu_vr(NE,M,K,L,NY,NX))+TX*OMBioResdu_vr(NE,M,K,L,NY,NX)
            ENDDO  
          ENDDO

          DO idom=idom_beg,idom_end
            DOM_vr(idom,K,L,NY,NX)=TI1*DOM_vr(idom,K,L,NY,NX)+FracTillCorp*(FI1*TDOM(idom,K) &
              -TI1*DOM_vr(idom,K,L,NY,NX))+TX*DOM_vr(idom,K,L,NY,NX)+FracTillCorp*DOM_MacP_vr(idom,K,L,NY,NX)              
            DOM_MacP_vr(idom,K,L,NY,NX)=XTillCorp_col(NY,NX)*DOM_MacP_vr(idom,K,L,NY,NX)
            SorbedOM_vr(idom,K,L,NY,NX)=TI1*SorbedOM_vr(idom,K,L,NY,NX)+FracTillCorp*(FI1*TOHM(idom,K) &
              -TI1*SorbedOM_vr(idom,K,L,NY,NX))+TX*SorbedOM_vr(idom,K,L,NY,NX)
          ENDDO  

          DO  M=1,jsken
            SolidOMAct_vr(M,K,L,NY,NX)=TI1*SolidOMAct_vr(M,K,L,NY,NX)+FracTillCorp*(FI1*TOSA(M,K) &
              -TI1*SolidOMAct_vr(M,K,L,NY,NX))+TX*SolidOMAct_vr(M,K,L,NY,NX)
            DO NE=1,NumPlantChemElms  
              SolidOM_vr(NE,M,K,L,NY,NX)=TI1*SolidOM_vr(NE,M,K,L,NY,NX)+FracTillCorp*(FI1*TOSM(NE,M,K) &
                -TI1*SolidOM_vr(NE,M,K,L,NY,NX))+TX*SolidOM_vr(NE,M,K,L,NY,NX)
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
                  mBiomeHeter_vr(NE,MID,K,L,NY,NX)=mBiomeHeter_vr(NE,MID,K,L,NY,NX)+FI1*TBOMGE(NE,MID,K)
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
                mBiomeAutor_vr(NE,MID,L,NY,NX)=mBiomeAutor_vr(NE,MID,L,NY,NX)+FI1*TOMGAutor(NE,MID)
              ENDDO
            enddo
          enddo
        ENDDO

        DO K=1,micpar%NumOfLitrCmplxs
          DO  M=1,ndbiomcp
            DO NE=1,NumPlantChemElms
              OMBioResdu_vr(NE,M,K,L,NY,NX)=OMBioResdu_vr(NE,M,K,L,NY,NX)+FI1*TORXE(NE,M,K)
            ENDDO
          ENDDO
          DO idom=idom_beg,idom_end
            DOM_vr(idom,K,L,NY,NX)=DOM_vr(idom,K,L,NY,NX)+FI1*TDOMG(idom,K)
            DOM_MacP_vr(idom,K,L,NY,NX)=DOM_MacP_vr(idom,K,L,NY,NX)+FI1*TDOMH(idom,K)
            SorbedOM_vr(idom,K,L,NY,NX)=SorbedOM_vr(idom,K,L,NY,NX)+FI1*TDOMHG(idom,K)
          ENDDO
          DO  M=1,jsken
            SolidOMAct_vr(M,K,L,NY,NX)=SolidOMAct_vr(M,K,L,NY,NX)+FI1*TOSGA(M,K)
            DO NE=1,NumPlantChemElms
              SolidOM_vr(NE,M,K,L,NY,NX)=SolidOM_vr(NE,M,K,L,NY,NX)+FI1*TOSGE(NE,M,K)
            ENDDO
          ENDDO
        ENDDO

        call sumORGMLayL(L,NY,NX,ORGM)

        call sumLitrOMLayL(L,NY,NX,litrOM)

        SoilOrgM_vr(:,L,NY,NX)=ORGM(:)

        OMLitrC_vr(L,NY,NX)=litrOM(ielmc)


        if(trc_solml_vr(idg_O2,1,1,1)<0._r8)then
          print*,'tillage problem',trc_solml_vr(idg_O2,1,1,1)
          stop
        endif

        trcx_solml_vr(idx_NH4,L,NY,NX)=trcx_solml_vr(idx_NH4,L,NY,NX)+FI1*TXN4G

        DO NTX=idx_AEC+1,idx_anion_soil_end
          trcx_solml_vr(NTX,L,NY,NX)=trcx_solml_vr(NTX,L,NY,NX)+FI1*TX_anion(NTX)
        ENDDO

        DO NTN=ifertn_beg,ifertn_end
          FertN_soil_vr(NTN,L,NY,NX)=FertN_soil_vr(NTN,L,NY,NX)+FI1*TFertNG_soil(NTN)
        ENDDO

        ZNHU0(L,NY,NX)=ZNHUX0
        ZNHUI(L,NY,NX)=ZNHUXI
        ZNFN0(L,NY,NX)=ZNFNX0
        ZNFNI(L,NY,NX)=(TI1*ZNFNI(L,NY,NX)+FracTillCorp*(FI1*TZNFNI &
          -TI1*ZNFNI(L,NY,NX))+TX*ZNFNI(L,NY,NX)+FI1*TZNFNG)/FI1
        TZNFN2=TZNFN2+ZNFNI(L,NY,NX)
      ENDIF
    ENDDO D2000

! nitrogen inhibitor
    ZNFN0(0,NY,NX) = ZNFNX0
    ZNFNI(0,NY,NX) = ZNFNI(0,NY,NX)*XCORP0
    TZNFN2         = TZNFN2+TZNFNG
    TZNFNI         = TZNFNI+TZNFNG
    DO  L=NU(NY,NX),LL
      IF(TZNFN2.GT.ZERO)THEN
        ZNFNI(L,NY,NX)=ZNFNI(L,NY,NX)*TZNFNI/TZNFN2
        ZNFNI(L,NY,NX)=ZNFNI(L,NY,NX)+0.5_r8*(ZNFN0(L,NY,NX)-ZNFNI(L,NY,NX))
      ENDIF
    ENDDO

  end subroutine ApplyTillageMixing


!--------------------------------------------------------------------------------

  subroutine MixSoluteMacpore(L,NY,NX)
  !
  !mix solutes in macropores
  !assuming macropores do not change
  implicit none
  integer, intent(in) :: L,NY,NX
  integer :: NTS,NTSA

  DO NTS=ids_beg,ids_end
    trc_soHml_vr(NTS,L,NY,NX)=XTillCorp_col(NY,NX)*trc_soHml_vr(NTS,L,NY,NX)
  ENDDO

  DO NTSA=idsalt_beg,idsalt_end
    trcSalt_soHml_vr(ntsa,L,NY,NX)=XTillCorp_col(NY,NX)*trcSalt_soHml_vr(ntsa,L,NY,NX)
  ENDDO
  end subroutine MixSoluteMacpore
!--------------------------------------------------------------------------------

  subroutine Mix1D(xdCorp,TI,FI,NU,LL,mass1d,XCORP0)

  implicit none
  real(r8), intent(in) :: xdCorp    !fraction of layer not mixed  
  real(r8), intent(in) :: TI(1:JZ)  !fraction of layer L over the tillage zone
  real(r8), intent(in) :: FI(1:JZ)  !fraction of layer L being tilled
  integer , intent(in) :: NU,LL
  real(r8), dimension(:), intent(inout) :: mass1D
  real(r8), optional, intent(in) :: XCORP0
  real(r8)  :: tmass1D,tmass1DR,TX,CORP0,dCorp
  integer :: L,L1,L0

  L0=0  
  if(present(XCORP0))then
    L0         = 1
    CORP0      = 1.0_r8-XCORP0
    tmass1DR   = mass1D(L0)*CORP0
    mass1D(L0) = mass1D(L0)*XCORP0
  endif

  tmass1D = 0._r8
  dCorp   = 1.0_r8-xdCorp
  DO L=NU,LL    
    L1=L+L0
    tmass1D = tmass1D+TI(L)*mass1D(L1)*dCorp
  ENDDO

  D100: DO L=NU,LL
    L1=L+L0
    if(TI(L)<1._r8)then
      TX=1._r8-TI(L)
      mass1D(L1)=TI(L)*mass1D(L1)*xdCorp+FI(L)*tmass1D+TX*mass1D(L1)
    else
      mass1D(L1)=mass1D(L1)*xdCorp+FI(L)*tmass1D
    endif
    if(present(XCORP0))then
      mass1D(L1)=mass1D(L1)+FI(L)*tmass1DR
    endif
  ENDDO D100

  end subroutine Mix1D

!--------------------------------------------------------------------------------

  subroutine Mix2D(xdCorp,TI,FI,NU,LL,mass2d,XCORP0)

  implicit none
  real(r8), intent(in) :: xdCorp    !fraction of layer not mixed  
  real(r8), intent(in) :: TI(1:JZ)  !fraction of layer L over the tillage zone
  real(r8), intent(in) :: FI(1:JZ)  !fraction of layer L being tilled
  integer , intent(in) :: NU,LL

  real(r8), dimension(:,:), intent(inout) :: mass2D
  real(r8), optional, intent(in) :: XCORP0
  real(r8)  :: tmass2DR(size(mass2D,1))
  real(r8)  :: tmass2D(size(mass2D,1))
  real(r8)  :: dCorp

  integer :: sz,ii,L,L0,L1
  real(r8) :: TX,CORP0   !fraction of layer L not tilled  

  sz=size(mass2D,1)
  L0=0
  if(present(XCORP0))then
    L0           = 1
    CORP0        = 1.0_r8-XCORP0
    tmass2DR     = mass2D(:,L0)*CORP0
    mass2D(:,L0) = mass2D(:,L0)*XCORP0
  endif

  tmass2D = 0._r8
  dCorp   = 1.0_r8-xdCorp
  DO L=NU,LL   
    L1=L+L0    
    DO ii=1,sz 
      tmass2D(ii) = tmass2D(ii)+TI(L)*mass2D(ii,L1)*dCorp
    enddo
  ENDDO

  D100: DO L=NU,LL
    L1=L+L0  
    do ii=1,sz
      if(TI(L)<1._r8)then
        TX=1._r8-TI(L)
        mass2D(ii,L1)=TI(L)*mass2D(ii,L1)*xdCorp+FI(L)*tmass2D(ii)+TX*mass2D(ii,L1)
      else
        mass2D(ii,L1)=mass2D(ii,L1)*xdCorp+FI(L)*tmass2D(ii)
      endif          
      if(present(XCORP0))then
        mass2D(ii,L1)=mass2D(ii,L1)+FI(L)*tmass2DR(ii)
      endif
    enddo
  ENDDO D100

  end subroutine Mix2D

!--------------------------------------------------------------------------------

  subroutine Mix3D(xdCorp,TI,FI,NU,LL,mass3d,XCORP0)

  implicit none
  real(r8), intent(in) :: xdCorp    !fraction of layer not mixed
  real(r8), intent(in) :: TI(1:JZ)  !fraction of layer L over the tillage zone
  real(r8), intent(in) :: FI(1:JZ)  !fraction of layer L being tilled
  integer , intent(in) :: NU,LL
  real(r8), dimension(:,:,:), intent(inout) :: mass3D
  real(r8), optional, intent(in) :: XCORP0
  real(r8)  :: tmass3D(size(mass3D,1),size(mass3D,2))
  real(r8)  :: tmass3DR(size(mass3D,1),size(mass3D,2))
  integer :: sz(3), ii, jj, L0,L1,L
  real(r8) :: TX,CORP0   !fraction of layer L not tilled  
  real(r8)  :: dCorp  

  sz=shape(mass3D)
  L0=0
  if(present(XCORP0))then
    L0             = 1
    CORP0          = 1.0_r8-XCORP0
    tmass3DR       = mass3D(:,:,L0)*CORP0
    mass3D(:,:,L0) = mass3D(:,:,L0)*XCORP0
  endif

  tmass3D = 0._r8
  dCorp   = 1.0_r8-xdCorp
  DO L=NU,LL  
    L1=L+L0  
    DO jj=1,sz(2) 
      DO ii=1,sz(1) 
        tmass3D(ii,jj) = tmass3D(ii,jj)+TI(L)*mass3D(ii,jj,L1)*dCorp
      enddo
    ENDDO
  ENDDO

  D100: DO L=NU,LL
    L1=L+L0  
    do jj=1,sz(2)
      do ii=1,sz(1)
          if(TI(L)<1._r8)then
            TX=1._r8-TI(L)
            mass3D(ii,jj,L1)=TI(L)*(mass3D(ii,jj,L1)*xdCorp+FI(L)*tmass3D(ii,jj))+TX*mass3D(ii,jj,L1)
          else
            mass3D(ii,jj,L1)=mass3D(ii,jj,L1)*xdCorp+FI(L)*tmass3D(ii,jj)
          endif          
          if(present(XCORP0))then
            mass3D(ii,jj,L1)=mass3D(ii,jj,L1)+FI(L1)*tmass3DR(ii,jj)
          endif
      enddo
    enddo
  ENDDO D100

  end subroutine Mix3D

!--------------------------------------------------------------------------------

  subroutine Mix4D(xdCorp,TI,FI,NU,LL,mass4d,XCORP0)

  implicit none
  real(r8), intent(in) :: xdCorp    !fraction of layer not mixed
  real(r8), intent(in) :: TI(1:JZ)  !fraction of layer L over the tillage zone
  real(r8), intent(in) :: FI(1:JZ)  !fraction of layer L being tilled
  integer,  intent(in) :: NU,LL
  real(r8), dimension(:,:,:,:), intent(inout) :: mass4D
  real(r8), optional, intent(in) :: XCORP0
  real(r8)  :: tmass4DR(size(mass4D,1),size(mass4D,2),size(mass4D,3))  
  real(r8)  :: tmass4D(size(mass4D,1),size(mass4D,2),size(mass4D,3))  
  real(r8) :: TX,CORP0   !fraction of layer L not tilled  
  integer :: sz(4), ii, jj,kk, L,L0,L1
  real(r8)  :: dCorp  
  
  sz=shape(mass4D)
  L0=0
  if(present(XCORP0))then
    L0               = 1
    CORP0            = 1.0_r8-XCORP0
    tmass4DR         = mass4D(:,:,:,L0)*CORP0
    mass4D(:,:,:,L0) = mass4D(:,:,:,L0)*XCORP0
  endif

  tmass4D = 0._r8
  dCorp   = 1.0_r8-xdCorp
  DO L=NU,LL  
    L1=L+L0  
    do kk=1,sz(3)
      DO jj=1,sz(2) 
        DO ii=1,sz(1) 
          tmass4D(ii,jj,kk) = tmass4D(ii,jj,kk)+TI(L)*mass4D(ii,jj,kk,L1)*dCorp
        enddo
      ENDDO
    enddo
  ENDDO

  D100: DO L=NU,LL
    L1=L+L0  
    do kk=1,sz(3)
      do jj=1,sz(2)
        do ii=1,sz(1)
          if(TI(L)<1._r8)then
            TX=1._r8-TI(L)
            mass4D(ii,jj,kk,L1)=TI(L)*mass4D(ii,jj,kk,L1)*xdCorp+FI(L)*tmass4D(ii,jj,kk)+TX*mass4D(ii,jj,kk,L1)
          else
            mass4D(ii,jj,kk,L1)=mass4D(ii,jj,kk,L1)*xdCorp+FI(L)*tmass4D(ii,jj,kk)
          endif     
          if(present(XCORP0))then
            mass4D(ii,jj,kk,L1)=mass4D(ii,jj,kk,L1)+FI(L)*tmass4DR(ii,jj,kk)
          endif
        enddo
      enddo
    enddo
  ENDDO D100

  end subroutine Mix4D


 !--------------------------------------------------------------------------------
  subroutine DeriveTillageProf(NY,NX,DCORPZ,LL,FI,TI)
  !derive the mass allocation profile for tillage
  implicit none
  integer , intent(in) :: NY,NX
  real(r8), intent(in) :: DCORPZ
  integer , intent(out) :: LL
  real(r8), intent(out) :: FI(1:JZ)   !fraction of current layer over the overall mixing zone
  real(r8), intent(out) :: TI(1:JZ)   !fraction of layer L being tilled
  integer :: L
  logical :: tillTestL
  real(r8) :: TL  !fraction of layer L over the tillage zone

  LL=NU(NY,NX)
  D100: DO L=NU(NY,NX),NL(NY,NX)
    !if current layer is within the tillage depth
    tillTestL=CumSoilThickness_vr(L,NY,NX)-DLYR(3,L,NY,NX).LT.DCORPZ .AND. DLYR(3,L,NY,NX).GT.ZERO
    if(.not.tillTestL)exit
    TL = AMIN1(DLYR(3,L,NY,NX),DCORPZ-(CumSoilThickness_vr(L,NY,NX)-DLYR(3,L,NY,NX)))
    FI(L) = TL/DCORPZ
    TI(L) = AMIN1(TL/DLYR(3,L,NY,NX),1._r8)
    LL=L  
  ENDDO D100
  end subroutine DeriveTillageProf

end module TillageMixMod
