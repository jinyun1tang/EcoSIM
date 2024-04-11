module NitroDisturbMod

! USES:
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use abortutils  , only : endrun
  use minimathmod, only : safe_adb
  use MicrobialDataType
  use EcoSiMParDataMod, only : micpar
  use NitroPars
  use SOMDataType
  use ChemTranspDataType
  use FertilizerDataType
  use GridConsts
  use SoilPhysDataType
  use FlagDataType
  use SoilHeatDatatype
  use SoilWaterDataType
  use ClimForcDataType
  use EcosimConst
  use PlantMgmtDataType
  use EcoSimSumDataType
  use EcosimBGCFluxType
  use SoilBGCDataType
  use AqueChemDatatype
  use GridDataType
  implicit none

  private

  character(len=*), parameter :: mod_filename = &
  __FILE__

  public :: SOMRemovalByDisturbance
  contains
!------------------------------------------------------------------------------------------

  subroutine SOMRemovalByDisturbance(I,J,NY,NX)
!
!     Description:
!
  implicit none
  integer, intent(in) :: I,J,NY,NX

  integer :: L,K,M,N,IFLGJ,NLL,NGL,NTF,MID
  real(r8) :: DC,DN,DP
  real(r8) :: DCORPC
  real(r8) :: FORGCX
  real(r8) :: HFLXD
  real(r8) :: OC,ON,OP,OCH,ONH,OPH,ONX
  REAL(R8) :: OPX,OCA,OAH
  real(r8) :: ONL(4,1:jcplx),OPL(4,1:jcplx)
  real(r8) :: DCORPC1
!     begin_execution

  IF(J.EQ.INT(SolarNoonHour_col(NY,NX)).AND.(ITILL(I,NY,NX).EQ.21.OR.ITILL(I,NY,NX).EQ.22))THEN
    IF(ITILL(I,NY,NX).EQ.22)THEN
      IFLGS(NY,NX)=1
      IFLGJ=0
      NLL=-1
      D2945: DO L=0,NL(NY,NX)

        IF(L.EQ.0.OR.L.GE.NUM(NY,NX))THEN
          IF(IFLGJ.EQ.1)THEN
            exit
          ELSEIF(THETW(L,NY,NX).GT.FVLWB.OR.CORGC(L,NY,NX).LE.FORGC)THEN
            IFLGJ=1
          ELSE
            NLL=L
          ENDIF
        ENDIF
      ENDDO D2945
    ELSE
      NLL=0
    ENDIF

    D2950: DO L=0,NLL
      IF(NLL.GE.0)THEN
        IF(ITILL(I,NY,NX).EQ.22)THEN
          IF(L.EQ.0)THEN
            FORGCX=0.0_r8
          ELSE
            FORGCX=FORGC
          ENDIF
          DCORPC=AMIN1(0.999_r8,DCORP(I,NY,NX))*(CORGC(L,NY,NX)-FORGCX) &
            /(AMAX1(CORGC(L,NY,NX),orgcden)-FORGCX)
        ELSE
          DCORPC=AMIN1(0.999_r8,DCORP(I,NY,NX))
        ENDIF
!     VOLWOU=VOLWOU+DCORPC*VLWatMicP(L,NY,NX)
!     HEATOU=HEATOU+DCORPC*4.19*TKS(L,NY,NX)*VLWatMicP(L,NY,NX)
!     VLWatMicP(L,NY,NX)=VLWatMicP(L,NY,NX)-DCORPC*VLWatMicP(L,NY,NX)
        OC=0.0_r8
        ON=0.0_r8
        OP=0.0_r8
        DC=0.0_r8
        DN=0.0_r8
        DP=0.0_r8
        ONL(1:4,1:jcplx)=0._r8
        OPL(1:4,1:jcplx)=0._r8

        D2970: DO K=1,jcplx
          IF(L.NE.0.OR.(micpar%is_litter(K)))THEN
!
!     REMOVE MICROBIAL BIOMASS
!
            D2960: DO N=1,NumMicbFunGroups
              DO NGL=JGnio(N),JGnfo(N)
                DO M=1,nlbiomcp
                  MID=micpar%get_micb_id(M,NGL)
                  OCH=DCORPC*OMEhetr(ielmc,MID,K,L,NY,NX)
                  ONH=DCORPC*OMEhetr(ielmn,MID,K,L,NY,NX)
                  OPH=DCORPC*OMEhetr(ielmp,MID,K,L,NY,NX)
                  ONX=EFIRE(1,ITILL(I,NY,NX))*ONH
                  OPX=EFIRE(2,ITILL(I,NY,NX))*OPH
                  IF(micpar%is_litter(K))THEN
                    ONL(4,K)=ONL(4,K)+ONH-ONX
                    OPL(4,K)=OPL(4,K)+OPH-OPX
                  ELSE
                    ONL(1,K)=ONL(1,K)+ONH-ONX
                    OPL(1,K)=OPL(1,K)+OPH-OPX
                  ENDIF
                  OMEhetr(ielmc,MID,K,L,NY,NX)=OMEhetr(ielmc,MID,K,L,NY,NX)-OCH
                  OMEhetr(ielmn,MID,K,L,NY,NX)=OMEhetr(ielmn,MID,K,L,NY,NX)-ONH
                  OMEhetr(ielmp,MID,K,L,NY,NX)=OMEhetr(ielmp,MID,K,L,NY,NX)-OPH
                  DC=DC+OMEhetr(ielmc,MID,K,L,NY,NX)
                  DN=DN+OMEhetr(ielmn,MID,K,L,NY,NX)
                  DP=DP+OMEhetr(ielmp,MID,K,L,NY,NX)
                  OC=OC+OCH
                  ON=ON+ONX
                  OP=OP+OPX
                enddo
              enddo
            ENDDO D2960
          ENDIF
        ENDDO D2970


!          IF(L.NE.0.OR.(K.NE.3.AND.K.NE.4))THEN
!
!     REMOVE MICROBIAL BIOMASS
!
        DO  N=1,NumMicbFunGroups
          DO NGL=JGniA(N),JGnfA(N)
            DO M=1,nlbiomcp
              MID=micpar%get_micb_id(M,NGL)
              OCH=DCORPC*OMEauto(ielmc,MID,L,NY,NX)
              ONH=DCORPC*OMEauto(ielmn,MID,L,NY,NX)
              OPH=DCORPC*OMEauto(ielmp,MID,L,NY,NX)
              ONX=EFIRE(1,ITILL(I,NY,NX))*ONH
              OPX=EFIRE(2,ITILL(I,NY,NX))*OPH
              ONL(4,1)=ONL(4,1)+ONH-ONX
              OPL(4,1)=OPL(4,1)+OPH-OPX
              OMEauto(ielmc,MID,L,NY,NX)=OMEauto(ielmc,MID,L,NY,NX)-OCH
              OMEauto(ielmn,MID,L,NY,NX)=OMEauto(ielmn,MID,L,NY,NX)-ONH
              OMEauto(ielmp,MID,L,NY,NX)=OMEauto(ielmp,MID,L,NY,NX)-OPH
              DC=DC+OMEauto(ielmc,MID,L,NY,NX)
              DN=DN+OMEauto(ielmn,MID,L,NY,NX)
              DP=DP+OMEauto(ielmp,MID,L,NY,NX)
              OC=OC+OCH
              ON=ON+ONX
              OP=OP+OPX
            enddo
          enddo
        ENDDO
!          ENDIF

!
!     REMOVE MICROBIAL RESIDUE
!
        D2900: DO K=1,jcplx
          IF(L.NE.0.OR.(micpar%is_litter(K)))THEN
            D2940: DO M=1,ndbiomcp
              OCH=DCORPC*ORM(ielmc,M,K,L,NY,NX)
              ONH=DCORPC*ORM(ielmn,M,K,L,NY,NX)
              OPH=DCORPC*ORM(ielmp,M,K,L,NY,NX)
              ONX=EFIRE(1,ITILL(I,NY,NX))*ONH
              OPX=EFIRE(2,ITILL(I,NY,NX))*OPH
              IF(micpar%is_litter(K))THEN
                ONL(4,K)=ONL(4,K)+ONH-ONX
                OPL(4,K)=OPL(4,K)+OPH-OPX
              ELSE
                ONL(1,K)=ONL(1,K)+ONH-ONX
                OPL(1,K)=OPL(1,K)+OPH-OPX
              ENDIF
              ORM(ielmc,M,K,L,NY,NX)=ORM(ielmc,M,K,L,NY,NX)-OCH
              ORM(ielmn,M,K,L,NY,NX)=ORM(ielmn,M,K,L,NY,NX)-ONH
              ORM(ielmp,M,K,L,NY,NX)=ORM(ielmp,M,K,L,NY,NX)-OPH
              DC=DC+ORM(ielmc,M,K,L,NY,NX)
              DN=DN+ORM(ielmn,M,K,L,NY,NX)
              DP=DP+ORM(ielmp,M,K,L,NY,NX)
              OC=OC+OCH
              ON=ON+ONX
              OP=OP+OPX
            ENDDO D2940
!
!     REMOVE DOC, DON, DOP
!
            OCH=DCORPC*DOM(idom_doc,K,L,NY,NX)
            OCA=DCORPC*DOM(idom_acetate,K,L,NY,NX)
            ONH=DCORPC*DOM(idom_don,K,L,NY,NX)
            OPH=DCORPC*DOM(idom_dop,K,L,NY,NX)
            ONX=EFIRE(1,ITILL(I,NY,NX))*ONH
            OPX=EFIRE(2,ITILL(I,NY,NX))*OPH
            IF(micpar%is_litter(K))THEN
              ONL(4,K)=ONL(4,K)+ONH-ONX
              OPL(4,K)=OPL(4,K)+OPH-OPX
            ELSE
              ONL(1,K)=ONL(1,K)+ONH-ONX
              OPL(1,K)=OPL(1,K)+OPH-OPX
            ENDIF
            DOM(idom_doc,K,L,NY,NX)=DOM(idom_doc,K,L,NY,NX)-OCH
            DOM(idom_acetate,K,L,NY,NX)=DOM(idom_acetate,K,L,NY,NX)-OCA
            DOM(idom_don,K,L,NY,NX)=DOM(idom_don,K,L,NY,NX)-ONH
            DOM(idom_dop,K,L,NY,NX)=DOM(idom_dop,K,L,NY,NX)-OPH
            OC=OC+OCH+OCA
            ON=ON+ONX
            OP=OP+OPX
            OCH=DCORPC*DOM_Macp(idom_doc,K,L,NY,NX)
            ONH=DCORPC*DOM_Macp(idom_don,K,L,NY,NX)
            OPH=DCORPC*DOM_Macp(idom_dop,K,L,NY,NX)
            OAH=DCORPC*DOM_Macp(idom_acetate,K,L,NY,NX)
            ONX=EFIRE(1,ITILL(I,NY,NX))*ONH
            OPX=EFIRE(2,ITILL(I,NY,NX))*OPH
            IF(micpar%is_litter(K))THEN
              ONL(4,K)=ONL(4,K)+ONH-ONX
              OPL(4,K)=OPL(4,K)+OPH-OPX
            ELSE
              ONL(1,K)=ONL(1,K)+ONH-ONX
              OPL(1,K)=OPL(1,K)+OPH-OPX
            ENDIF
            DOM_Macp(idom_doc,K,L,NY,NX)=DOM_Macp(idom_doc,K,L,NY,NX)-OCH
            DOM_Macp(idom_don,K,L,NY,NX)=DOM_Macp(idom_don,K,L,NY,NX)-ONH
            DOM_Macp(idom_dop,K,L,NY,NX)=DOM_Macp(idom_dop,K,L,NY,NX)-OPH
            DOM_Macp(idom_acetate,K,L,NY,NX)=DOM_Macp(idom_acetate,K,L,NY,NX)-OAH
            OC=OC+OCH+OAH
            ON=ON+ONX
            OP=OP+OPX
!
!     REMOVE ADSORBED OM
!
            OCH=DCORPC*OHM(ielmc,K,L,NY,NX)
            ONH=DCORPC*OHM(ielmn,K,L,NY,NX)
            OPH=DCORPC*OHM(ielmp,K,L,NY,NX)
            OAH=DCORPC*OHM(idom_acetate,K,L,NY,NX)
            ONX=EFIRE(1,ITILL(I,NY,NX))*ONH
            OPX=EFIRE(2,ITILL(I,NY,NX))*OPH
            IF(micpar%is_litter(K))THEN
              ONL(4,K)=ONL(4,K)+ONH-ONX
              OPL(4,K)=OPL(4,K)+OPH-OPX
            ELSE
              ONL(1,K)=ONL(1,K)+ONH-ONX
              OPL(1,K)=OPL(1,K)+OPH-OPX
            ENDIF
            OHM(ielmc,K,L,NY,NX)=OHM(ielmc,K,L,NY,NX)-OCH
            OHM(ielmn,K,L,NY,NX)=OHM(ielmn,K,L,NY,NX)-ONH
            OHM(ielmp,K,L,NY,NX)=OHM(ielmp,K,L,NY,NX)-OPH
            OHM(idom_acetate,K,L,NY,NX)=OHM(idom_acetate,K,L,NY,NX)-OAH
            DC=DC+DOM(idom_doc,K,L,NY,NX)+DOM_Macp(idom_doc,K,L,NY,NX)+OHM(ielmc,K,L,NY,NX) &
              +DOM(idom_acetate,K,L,NY,NX)+DOM_Macp(idom_acetate,K,L,NY,NX)+OHM(idom_acetate,K,L,NY,NX)
            DN=DN+DOM(idom_don,K,L,NY,NX)+DOM_Macp(idom_don,K,L,NY,NX)+OHM(ielmn,K,L,NY,NX)
            DP=DP+DOM(idom_dop,K,L,NY,NX)+DOM_Macp(idom_dop,K,L,NY,NX)+OHM(ielmp,K,L,NY,NX)
            OC=OC+OCH
            ON=ON+ONX
            OP=OP+OPX
!
!     REMOVE RESIDUE
!
            D2930: DO M=1,jsken
              OCH=DCORPC*OSM(ielmc,M,K,L,NY,NX)
              ONH=DCORPC*OSM(ielmn,M,K,L,NY,NX)
              OPH=DCORPC*OSM(ielmp,M,K,L,NY,NX)
              OCA=DCORPC*OSA(M,K,L,NY,NX)
              ONX=EFIRE(1,ITILL(I,NY,NX))*ONH
              OPX=EFIRE(2,ITILL(I,NY,NX))*OPH
              ONL(M,K)=ONL(M,K)+ONH-ONX
              OPL(M,K)=OPL(M,K)+OPH-OPX
              OSM(ielmc,M,K,L,NY,NX)=OSM(ielmc,M,K,L,NY,NX)-OCH
              OSM(ielmn,M,K,L,NY,NX)=OSM(ielmn,M,K,L,NY,NX)-ONH
              OSM(ielmp,M,K,L,NY,NX)=OSM(ielmp,M,K,L,NY,NX)-OPH
              OSA(M,K,L,NY,NX)=OSA(M,K,L,NY,NX)-OCA              
              DC=DC+OSM(ielmc,M,K,L,NY,NX)
              DN=DN+OSM(ielmn,M,K,L,NY,NX)
              DP=DP+OSM(ielmp,M,K,L,NY,NX)
              OC=OC+OCH
              ON=ON+ONX
              OP=OP+OPX
            ENDDO D2930
          ENDIF
        ENDDO D2900
!
!     ADD UNBURNED N,P TO ORG N, ORG P
!
        D2905: DO K=1,jcplx
          DO  M=1,jsken
            OSM(ielmn,M,K,L,NY,NX)=OSM(ielmn,M,K,L,NY,NX)+ONL(M,K)
            OSM(ielmp,M,K,L,NY,NX)=OSM(ielmp,M,K,L,NY,NX)+OPL(M,K)
            DN=DN+ONL(M,K)
            DP=DP+OPL(M,K)
          enddo
        ENDDO D2905
!
!     REMOVE FERTILIZER IN RESIDUE
!
        IF(ITILL(I,NY,NX).EQ.21)THEN
          ON=ON+DCORPC*(trc_solml_vr(ids_NH4,L,NY,NX)+trc_solml_vr(idg_NH3,L,NY,NX) &
            +trc_solml_vr(ids_NO3,L,NY,NX)+trc_solml_vr(ids_NO2,L,NY,NX))
          OP=OP+DCORPC*(trc_solml_vr(ids_H1PO4,L,NY,NX)+trc_solml_vr(ids_H2PO4,L,NY,NX))
          DCORPC1=1.0_r8-DCORPC
          trc_solml_vr(ids_NH4,L,NY,NX)=DCORPC1*trc_solml_vr(ids_NH4,L,NY,NX)
          trc_solml_vr(idg_NH3,L,NY,NX)=DCORPC1*trc_solml_vr(idg_NH3,L,NY,NX)
          trc_solml_vr(ids_NO3,L,NY,NX)=DCORPC1*trc_solml_vr(ids_NO3,L,NY,NX)
          trc_solml_vr(ids_NO2,L,NY,NX)=DCORPC1*trc_solml_vr(ids_NO2,L,NY,NX)
          trc_solml_vr(ids_H1PO4,L,NY,NX)=DCORPC1*trc_solml_vr(ids_H1PO4,L,NY,NX)
          trc_solml_vr(ids_H2PO4,L,NY,NX)=DCORPC1*trc_solml_vr(ids_H2PO4,L,NY,NX)
          trcx_solml(idx_NH4,L,NY,NX)  =DCORPC1*trcx_solml(idx_NH4,L,NY,NX)
          trcp_salml(idsp_AlPO4,L,NY,NX)=DCORPC1*trcp_salml(idsp_AlPO4,L,NY,NX)
          trcp_salml(idsp_FePO4,L,NY,NX)=DCORPC1*trcp_salml(idsp_FePO4,L,NY,NX)
          trcp_salml(idsp_CaHPO4,L,NY,NX)=DCORPC1*trcp_salml(idsp_CaHPO4,L,NY,NX)
          trcp_salml(idsp_HA,L,NY,NX)=DCORPC1*trcp_salml(idsp_HA,L,NY,NX)
          trcp_salml(idsp_CaH4P2O8,L,NY,NX)=DCORPC1*trcp_salml(idsp_CaH4P2O8,L,NY,NX)

          DO NTF=ifertn_beg,ifertn_end
            FertN_soil(NTF,L,NY,NX)=DCORPC1*FertN_soil(NTF,L,NY,NX)
          ENDDO
        ENDIF
        ORGC(L,NY,NX)=DC
        ORGN(L,NY,NX)=DN
        IF(L.EQ.0)THEN
          HFLXD=4.19E-06_r8*(ORGCX(L,NY,NX)-ORGC(L,NY,NX))*TKS(L,NY,NX)
          HEATOU=HEATOU+HFLXD
        ENDIF
!     IF(L.EQ.0)THEN
!     VHeatCapacity(0,NY,NX)=2.496E-06*ORGC(0,NY,NX)+4.19*VLWatMicP(0,NY,NX)
!    2+1.9274*VLiceMicP(0,NY,NX)
!     ELSE
!     VHeatCapacity(L,NY,NX)=VHeatCapacitySoilM(L,NY,NX)+4.19*(VLWatMicP(L,NY,NX)+VLWatMacP(L,NY,NX))
!    2+1.9274*(VLiceMicP(L,NY,NX)+VLiceMacP(L,NY,NX))
!     ENDIF
        IF(ITILL(I,NY,NX).EQ.21)THEN
          TOMOU(ielmc)=TOMOU(ielmc)+OC
          TOMOU(ielmn)=TOMOU(ielmn)+ON
          TOMOU(ielmp)=TOMOU(ielmp)+OP
          HydroSufDOCFlx_col(NY,NX)=HydroSufDOCFlx_col(NY,NX)+OC
          HydroSufDONFlx_col(NY,NX)=HydroSufDONFlx_col(NY,NX)+ON
          HydroSufDOPFlx_col(NY,NX)=HydroSufDOPFlx_col(NY,NX)+OP
          Eco_NBP_col(NY,NX)=Eco_NBP_col(NY,NX)-OC
        ELSEIF(ITILL(I,NY,NX).EQ.22)THEN
          CO2GIN=CO2GIN-OC
          OXYGIN=OXYGIN+2.667_r8*OC
          OXYGOU=OXYGOU+2.667_r8*OC
          TOMOU(ielmn)=TOMOU(ielmn)+ON
          TOMOU(ielmp)=TOMOU(ielmp)+OP
          CO2byFire_col(NY,NX)=CO2byFire_col(NY,NX)-(1.0_r8-FCH4F)*OC
          CH4byFire_col(NY,NX)=CH4byFire_col(NY,NX)-FCH4F*OC
          O2byFire_col(NY,NX)=O2byFire_col(NY,NX)+(1.0_r8-FCH4F)*2.667_r8*OC
          NH3byFire_col(NY,NX)=NH3byFire_col(NY,NX)-ON
          N2ObyFire_col(NY,NX)=N2ObyFire_col(NY,NX)-0.0_r8
          PO4byFire_col(NY,NX)=PO4byFire_col(NY,NX)-OP
          Eco_NBP_col(NY,NX)=Eco_NBP_col(NY,NX)-OC
        ENDIF
      ENDIF
    ENDDO D2950
  ENDIF
  end subroutine SOMRemovalByDisturbance
!------------------------------------------------------------------------------------------


end module NitroDisturbMod
