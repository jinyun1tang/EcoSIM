module SoilDisturbMod

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
!ITILL,DCORP=disturbance type,intensity
!ITILL=soil disturbance type 1-20:tillage,21=litter removal,22=fire,23-24=drainage
  implicit none
  integer, intent(in) :: I,J,NY,NX

  integer :: L,K,M,N,IFLGJ,NLL,NGL,NTF,MID,ids,NE,idom
  real(r8) :: DC,DN,DP
  real(r8) :: DCORPC
  real(r8) :: FORGCX
  real(r8) :: HFLXD
  real(r8) :: ONX
  REAL(R8) :: OPX,OCA,OAH
  real(r8) :: ONL(4,1:jcplx),OPL(4,1:jcplx)
  real(r8) :: rmDOM(idom_beg:idom_end)
  real(r8) :: rmBiom(1:NumPlantChemElms )
  real(r8) :: OMelm(1:NumPlantChemElms)
  real(r8) :: dOMelm(1:NumPlantChemElms)
  real(r8) :: DCORPC1

!     begin_execution

  IF(J.EQ.INT(SolarNoonHour_col(NY,NX)) .AND. (iSoilDisturbType_col(I,NY,NX).EQ.21 .OR. iSoilDisturbType_col(I,NY,NX).EQ.22))THEN
    IF(iSoilDisturbType_col(I,NY,NX).EQ.22)THEN
      iResetSoilProf_col(NY,NX)=itrue
      IFLGJ=0
      NLL=-1

      D2945: DO L=0,NL(NY,NX)
        IF(L.EQ.0 .OR. L.GE.NUM(NY,NX))THEN
          IF(IFLGJ.EQ.1)THEN
            exit
          ELSEIF(THETW_vr(L,NY,NX).GT.VolMaxSoilMoist4Fire .OR. CSoilOrgM_vr(ielmc,L,NY,NX).LE.FORGC)THEN
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
        IF(iSoilDisturbType_col(I,NY,NX).EQ.22)THEN
          IF(L.EQ.0)THEN
            FORGCX=0.0_r8
          ELSE
            FORGCX=FORGC
          ENDIF
          DCORPC=AMIN1(0.999_r8,DepzCorp_col(I,NY,NX))*(CSoilOrgM_vr(ielmc,L,NY,NX)-FORGCX) &
            /(AMAX1(CSoilOrgM_vr(ielmc,L,NY,NX),orgcden)-FORGCX)
        ELSE
          DCORPC=AMIN1(0.999_r8,DepzCorp_col(I,NY,NX))
        ENDIF
!     QH2OLoss_lnds=QH2OLoss_lnds+DCORPC*VLWatMicP_vr(L,NY,NX)
!     HeatOut_lnds=HeatOut_lnds+DCORPC*4.19*TKS_vr(L,NY,NX)*VLWatMicP_vr(L,NY,NX)
!     VLWatMicP_vr(L,NY,NX)=VLWatMicP_vr(L,NY,NX)-DCORPC*VLWatMicP_vr(L,NY,NX)
        OMelm=0._r8
        dOMelm=0._r8

        ONL(1:4,1:jcplx)=0._r8
        OPL(1:4,1:jcplx)=0._r8

        D2970: DO K=1,jcplx
          IF(L.NE.0.OR.(micpar%is_litter(K)))THEN
!
!     REMOVE MICROBIAL BIOMASS
!
            D2960: DO N=1,NumMicbFunGrupsPerCmplx
              DO NGL=JGnio(N),JGnfo(N)
                DO M=1,nlbiomcp
                  MID=micpar%get_micb_id(M,NGL)
                  DO NE=1,NumPlantChemElms 
                    rmBiom(NE)=DCORPC*mBiomeHeter_vr(NE,MID,K,L,NY,NX)
                  enddo

                  ONX=EFIRE(1,iSoilDisturbType_col(I,NY,NX))*rmBiom(ielmn)
                  OPX=EFIRE(2,iSoilDisturbType_col(I,NY,NX))*rmBiom(ielmp)
                  IF(micpar%is_litter(K))THEN
                    ONL(4,K)=ONL(4,K)+rmBiom(ielmn)-ONX
                    OPL(4,K)=OPL(4,K)+rmBiom(ielmp)-OPX
                  ELSE
                    ONL(1,K)=ONL(1,K)+rmBiom(ielmn)-ONX
                    OPL(1,K)=OPL(1,K)+rmBiom(ielmp)-OPX
                  ENDIF
                  DO  NE=1,NumPlantChemElms 
                    mBiomeHeter_vr(NE,MID,K,L,NY,NX)=mBiomeHeter_vr(NE,MID,K,L,NY,NX)-rmBiom(NE)
                    dOMelm(NE)=dOMelm(NE)+mBiomeHeter_vr(NE,MID,K,L,NY,NX)                    
                  ENDDO
                  OMelm(ielmc)=OMelm(ielmc)+rmBiom(ielmc)
                  OMelm(ielmn)=OMelm(ielmn)+ONX
                  OMelm(ielmp)=OMelm(ielmp)+OPX
                enddo
              enddo
            ENDDO D2960
          ENDIF
        ENDDO D2970


!          IF(L.NE.0.OR.(K.NE.3.AND.K.NE.4))THEN
!
!     REMOVE MICROBIAL BIOMASS
!
        DO  N=1,NumMicbFunGrupsPerCmplx
          DO NGL=JGniA(N),JGnfA(N)
            DO M=1,nlbiomcp
              MID=micpar%get_micb_id(M,NGL)
              DO NE=1,NumPlantChemElms
                rmBiom(NE)=DCORPC*mBiomeAutor_vr(NE,MID,L,NY,NX)
              enddo

              ONX=EFIRE(1,iSoilDisturbType_col(I,NY,NX))*rmBiom(ielmn)
              OPX=EFIRE(2,iSoilDisturbType_col(I,NY,NX))*rmBiom(ielmp)

              ONL(4,1)=ONL(4,1)+rmBiom(ielmn)-ONX
              OPL(4,1)=OPL(4,1)+rmBiom(ielmp)-OPX

              DO NE=1,NumPlantChemElms
                mBiomeAutor_vr(NE,MID,L,NY,NX)=mBiomeAutor_vr(NE,MID,L,NY,NX)-rmBiom(NE)
                dOMelm(NE)=dOMelm(NE)+mBiomeAutor_vr(NE,MID,L,NY,NX)
              enddo
              OMelm(ielmc)=OMelm(ielmc)+rmBiom(ielmc)
              OMelm(ielmn)=OMelm(ielmn)+ONX
              OMelm(ielmp)=OMelm(ielmp)+OPX              
            enddo
          enddo
        ENDDO
!          ENDIF

!
!     REMOVE MICROBIAL RESIDUE
!
        D2900: DO K=1,jcplx
          IF(L.NE.0 .OR. (micpar%is_litter(K)))THEN
            D2940: DO M=1,ndbiomcp
              DO NE=1,NumPlantChemElms
                rmBiom(NE)=DCORPC*OMBioResdu_vr(NE,M,K,L,NY,NX)
              enddo

              ONX=EFIRE(1,iSoilDisturbType_col(I,NY,NX))*rmBiom(ielmn)
              OPX=EFIRE(2,iSoilDisturbType_col(I,NY,NX))*rmBiom(ielmp)
              IF(micpar%is_litter(K))THEN
                ONL(4,K)=ONL(4,K)+rmBiom(ielmn)-ONX
                OPL(4,K)=OPL(4,K)+rmBiom(ielmp)-OPX
              ELSE
                ONL(1,K)=ONL(1,K)+rmBiom(ielmn)-ONX
                OPL(1,K)=OPL(1,K)+rmBiom(ielmp)-OPX
              ENDIF
              DO NE=1,NumPlantChemElms
                OMBioResdu_vr(NE,M,K,L,NY,NX)=OMBioResdu_vr(NE,M,K,L,NY,NX)-rmBiom(NE)
                dOMelm(NE)=dOMelm(NE)+OMBioResdu_vr(NE,M,K,L,NY,NX)
              enddo
              OMelm(ielmc)=OMelm(ielmc)+rmBiom(ielmc)
              OMelm(ielmn)=OMelm(ielmn)+ONX
              OMelm(ielmp)=OMelm(ielmp)+OPX                            
            ENDDO D2940
!
!     REMOVE DOC, DON, DOP
!
            DO idom=idom_beg,idom_end
              rmDOM(idom)=DCORPC*DOM_vr(idom,K,L,NY,NX)
            enddo

            ONX=EFIRE(1,iSoilDisturbType_col(I,NY,NX))*rmDOM(ielmn)
            OPX=EFIRE(2,iSoilDisturbType_col(I,NY,NX))*rmDOM(ielmp)
            IF(micpar%is_litter(K))THEN
              ONL(4,K)=ONL(4,K)+rmDOM(ielmn)-ONX
              OPL(4,K)=OPL(4,K)+rmDOM(ielmp)-OPX
            ELSE
              ONL(1,K)=ONL(1,K)+rmDOM(ielmn)-ONX
              OPL(1,K)=OPL(1,K)+rmDOM(ielmp)-OPX
            ENDIF
            DO idom=idom_beg,idom_end
              DOM_vr(idom,K,L,NY,NX)=DOM_vr(idom,K,L,NY,NX)-rmDOM(idom)
            enddo

            OMelm(ielmc)=OMelm(ielmc)+rmDOM(idom_doc)+rmDOM(idom_acetate)            
            OMelm(ielmn)=OMelm(ielmn)+ONX
            OMelm(ielmp)=OMelm(ielmp)+OPX

            DO idom=idom_beg,idom_end
              rmDOM(idom)=DCORPC*DOM_MacP_vr(idom,K,L,NY,NX)
            enddo

            ONX=EFIRE(1,iSoilDisturbType_col(I,NY,NX))*rmDOM(ielmn)
            OPX=EFIRE(2,iSoilDisturbType_col(I,NY,NX))*rmDOM(ielmp)
            IF(micpar%is_litter(K))THEN
              ONL(4,K)=ONL(4,K)+rmDOM(ielmn)-ONX
              OPL(4,K)=OPL(4,K)+rmDOM(ielmp)-OPX
            ELSE
              ONL(1,K)=ONL(1,K)+rmDOM(ielmn)-ONX
              OPL(1,K)=OPL(1,K)+rmDOM(ielmp)-OPX
            ENDIF
            DO idom=idom_beg,idom_end
              DOM_MacP_vr(idom,K,L,NY,NX)=DOM_MacP_vr(idom,K,L,NY,NX)-rmDOM(idom)
            enddo

            OMelm(ielmc)=OMelm(ielmc)+rmDOM(idom_doc)+rmDOM(idom_acetate)            
            OMelm(ielmn)=OMelm(ielmn)+ONX
            OMelm(ielmp)=OMelm(ielmp)+OPX
!
!     REMOVE ADSORBED OM
!
            DO idom=idom_beg,idom_end
              rmDOM(idom)=DCORPC*SorbedOM_vr(idom,K,L,NY,NX)
            enddo

            ONX=EFIRE(1,iSoilDisturbType_col(I,NY,NX))*rmDOM(ielmn)
            OPX=EFIRE(2,iSoilDisturbType_col(I,NY,NX))*rmDOM(ielmp)
            IF(micpar%is_litter(K))THEN
              ONL(4,K)=ONL(4,K)+rmDOM(ielmn)-ONX
              OPL(4,K)=OPL(4,K)+rmDOM(ielmp)-OPX
            ELSE
              ONL(1,K)=ONL(1,K)+rmDOM(ielmn)-ONX
              OPL(1,K)=OPL(1,K)+rmDOM(ielmp)-OPX
            ENDIF
            DO idom=idom_beg,idom_end
              SorbedOM_vr(idom,K,L,NY,NX)=SorbedOM_vr(idom,K,L,NY,NX)-rmDOM(idom)
            enddo

            dOMelm(ielmc)=dOMelm(ielmc)+DOM_vr(idom_doc,K,L,NY,NX)+DOM_MacP_vr(idom_doc,K,L,NY,NX)+SorbedOM_vr(ielmc,K,L,NY,NX) &
              +DOM_vr(idom_acetate,K,L,NY,NX)+DOM_MacP_vr(idom_acetate,K,L,NY,NX)+SorbedOM_vr(idom_acetate,K,L,NY,NX)
            dOMelm(ielmn)=dOMelm(ielmn)+DOM_vr(idom_don,K,L,NY,NX)+DOM_MacP_vr(idom_don,K,L,NY,NX)+SorbedOM_vr(ielmn,K,L,NY,NX)
            dOMelm(ielmp)=dOMelm(ielmp)+DOM_vr(idom_dop,K,L,NY,NX)+DOM_MacP_vr(idom_dop,K,L,NY,NX)+SorbedOM_vr(ielmp,K,L,NY,NX)

            OMelm(ielmc)=OMelm(ielmc)+rmDOM(idom_doc)+rmDOM(idom_acetate)
            OMelm(ielmn)=OMelm(ielmn)+ONX
            OMelm(ielmp)=OMelm(ielmp)+OPX!
!     REMOVE RESIDUE
!
            D2930: DO M=1,jsken
              do NE=1,NumPlantChemElms
                rmBiom(NE)=DCORPC*SolidOM_vr(NE,M,K,L,NY,NX)
              enddo

              OCA=DCORPC*SolidOMAct_vr(M,K,L,NY,NX)
              ONX=EFIRE(1,iSoilDisturbType_col(I,NY,NX))*rmBiom(ielmn)
              OPX=EFIRE(2,iSoilDisturbType_col(I,NY,NX))*rmBiom(ielmp)
              ONL(M,K)=ONL(M,K)+rmBiom(ielmn)-ONX
              OPL(M,K)=OPL(M,K)+rmBiom(ielmp)-OPX

              do NE=1,NumPlantChemElms
                SolidOM_vr(NE,M,K,L,NY,NX)=SolidOM_vr(NE,M,K,L,NY,NX)-rmBiom(NE)
                dOMelm(NE)=dOMelm(NE)+SolidOM_vr(NE,M,K,L,NY,NX)
              enddo
                            
              SolidOMAct_vr(M,K,L,NY,NX)=SolidOMAct_vr(M,K,L,NY,NX)-OCA              
              OMelm(ielmc)=OMelm(ielmc)+rmBiom(ielmc)
              OMelm(ielmn)=OMelm(ielmn)+ONX
              OMelm(ielmp)=OMelm(ielmp)+OPX      
            ENDDO D2930
          ENDIF
        ENDDO D2900
!
!     ADD UNBURNED N,P TO ORG N, ORG P
!
        D2905: DO K=1,jcplx
          DO  M=1,jsken
            SolidOM_vr(ielmn,M,K,L,NY,NX)=SolidOM_vr(ielmn,M,K,L,NY,NX)+ONL(M,K)
            SolidOM_vr(ielmp,M,K,L,NY,NX)=SolidOM_vr(ielmp,M,K,L,NY,NX)+OPL(M,K)
            DN=DN+ONL(M,K)
            DP=DP+OPL(M,K)
          enddo
        ENDDO D2905
!
!     REMOVE FERTILIZER IN RESIDUE
!
        IF(iSoilDisturbType_col(I,NY,NX).EQ.21)THEN
          OMelm(ielmn)=OMelm(ielmn)+DCORPC*(trcs_solml_vr(ids_NH4,L,NY,NX)+trcs_solml_vr(idg_NH3,L,NY,NX) &
            +trcs_solml_vr(ids_NO3,L,NY,NX)+trcs_solml_vr(ids_NO2,L,NY,NX))
          OMelm(ielmp)=OMelm(ielmp)+DCORPC*(trcs_solml_vr(ids_H1PO4,L,NY,NX)+trcs_solml_vr(ids_H2PO4,L,NY,NX))
          
          DCORPC1=1.0_r8-DCORPC

          trcs_solml_vr(idg_NH3,L,NY,NX)=DCORPC1*trcs_solml_vr(idg_NH3,L,NY,NX)
          do ids=ids_nut_beg,ids_nuts_end
            trcs_solml_vr(ids,L,NY,NX)=DCORPC1*trcs_solml_vr(ids,L,NY,NX)
          enddo
          trcx_solml_vr(idx_NH4,L,NY,NX)  =DCORPC1*trcx_solml_vr(idx_NH4,L,NY,NX)
          trcp_saltpml_vr(idsp_AlPO4,L,NY,NX)    = DCORPC1*trcp_saltpml_vr(idsp_AlPO4,L,NY,NX)
          trcp_saltpml_vr(idsp_FePO4,L,NY,NX)    = DCORPC1*trcp_saltpml_vr(idsp_FePO4,L,NY,NX)
          trcp_saltpml_vr(idsp_CaHPO4,L,NY,NX)   = DCORPC1*trcp_saltpml_vr(idsp_CaHPO4,L,NY,NX)
          trcp_saltpml_vr(idsp_HA,L,NY,NX)       = DCORPC1*trcp_saltpml_vr(idsp_HA,L,NY,NX)
          trcp_saltpml_vr(idsp_CaH4P2O8,L,NY,NX) = DCORPC1*trcp_saltpml_vr(idsp_CaH4P2O8,L,NY,NX)

          DO NTF=ifertn_beg,ifertn_end
            FertN_mole_soil_vr(NTF,L,NY,NX)=DCORPC1*FertN_mole_soil_vr(NTF,L,NY,NX)
          ENDDO
        ENDIF
        SoilOrgM_vr(1:NumPlantChemElms,L,NY,NX)=dOMelm(1:NumPlantChemElms)
        IF(L.EQ.0)THEN
          HFLXD=4.19E-06_r8*(ORGCX_vr(L,NY,NX)-SoilOrgM_vr(ielmc,L,NY,NX))*TKS_vr(L,NY,NX)
          HeatOut_lnds=HeatOut_lnds+HFLXD
        ENDIF
!     IF(L.EQ.0)THEN
!     VHeatCapacity_vr(0,NY,NX)=2.496E-06*SoilOrgM_vr(ielmc,0,NY,NX)+4.19*VLWatMicP_vr(0,NY,NX)
!    2+1.9274*VLiceMicP_vr(0,NY,NX)
!     ELSE
!     VHeatCapacity_vr(L,NY,NX)=VHeatCapacitySoilM_vr(L,NY,NX)+4.19*(VLWatMicP_vr(L,NY,NX)+VLWatMacP_vr(L,NY,NX))
!    2+1.9274*(VLiceMicP_vr(L,NY,NX)+VLiceMacP_vr(L,NY,NX))
!     ENDIF
        IF(iSoilDisturbType_col(I,NY,NX).EQ.21)THEN
          DO NE=1,NumPlantChemElms
            TOMOU_lnds(NE)=TOMOU_lnds(NE)+OMelm(NE)
          ENDDO

          HydroSufDOCFlx_col(NY,NX)=HydroSufDOCFlx_col(NY,NX)+OMelm(ielmc)
          HydroSufDONFlx_CumYr_col(NY,NX)=HydroSufDONFlx_CumYr_col(NY,NX)+OMelm(ielmn)
          HydroSufDOPFlx_CumYr_col(NY,NX)=HydroSufDOPFlx_CumYr_col(NY,NX)+OMelm(ielmp)
          
          Eco_NBP_CumYr_col(NY,NX)=Eco_NBP_CumYr_col(NY,NX)-OMelm(ielmc)
        ELSEIF(iSoilDisturbType_col(I,NY,NX).EQ.22)THEN
          SurfGas_lnd(idg_CO2)=SurfGas_lnd(idg_CO2)-OMelm(ielmc)
          SurfGas_lnd(idg_O2)=SurfGas_lnd(idg_O2)+2.667_r8*OMelm(ielmc)
          OXYGOU=OXYGOU+2.667_r8*OMelm(ielmc)

          TOMOU_lnds(ielmn)=TOMOU_lnds(ielmn)+OMelm(ielmn)
          TOMOU_lnds(ielmp)=TOMOU_lnds(ielmp)+OMelm(ielmp)
          CO2byFire_CumYr_col(NY,NX)=CO2byFire_CumYr_col(NY,NX)-(1.0_r8-FrcAsCH4byFire)*OMelm(ielmc)
          CH4byFire_CumYr_col(NY,NX)=CH4byFire_CumYr_col(NY,NX)-FrcAsCH4byFire*OMelm(ielmc)
          O2byFire_CumYr_col(NY,NX)=O2byFire_CumYr_col(NY,NX)+(1.0_r8-FrcAsCH4byFire)*2.667_r8*OMelm(ielmc)
          NH3byFire_CumYr_col(NY,NX)=NH3byFire_CumYr_col(NY,NX)-OMelm(ielmn)
          N2ObyFire_CumYr_col(NY,NX)=N2ObyFire_CumYr_col(NY,NX)-0.0_r8
          PO4byFire_CumYr_col(NY,NX)=PO4byFire_CumYr_col(NY,NX)-OMelm(ielmp)
          Eco_NBP_CumYr_col(NY,NX)=Eco_NBP_CumYr_col(NY,NX)-OMelm(ielmc)
        ENDIF
      ENDIF
    ENDDO D2950
  ENDIF
  end subroutine SOMRemovalByDisturbance
!------------------------------------------------------------------------------------------


end module SoilDisturbMod
