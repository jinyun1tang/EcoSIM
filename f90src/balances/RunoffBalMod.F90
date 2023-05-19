module RunoffBalMod
  use data_kind_mod, only : r8 => DAT_KIND_R8
  USE EcoSimSumDataType
  use minimathmod, only : AZMAX1
  USE ChemTranspDataType
  use GridConsts
  use SnowDataType
  use FlagDataType
  use SoilBGCDataType
  USE SedimentDataType
  use MicrobialDataType
  USE SoilWaterDataType
  USE AqueChemDatatype
  USE EcoSIMCtrlDataType
  USE GridDataType
  use EcoSIMConfig, only : jcplx => jcplxc
  use EcosimConst, only : patomw,natomw
implicit none
  private
  character(len=*), parameter :: mod_filename = __FILE__
  public :: RunoffBal

  real(r8) :: SG

  contains

  subroutine RunoffBal(I,J,NY,NX,NHW,NHE,NVN,NVS)
  implicit none
  integer, intent(in) :: I,J, NY,NX,NHW,NHE,NVN,NVS
  integer :: N,N1,N2,N3,N4,N5,N6
  integer :: L,NN
  real(r8) :: XN
  real(r8) :: CXR,ZXR,PXR
  real(r8) :: ZGR

  SG=0._r8
  D9985: DO L=NU(NY,NX),NL(NY,NX)
!
!     LOCATE EXTERNAL BOUNDARIES
!
!     N2,N1=NY,NX of source grid cell
!     N5,N4=NY,NX of destination grid cell
!
    N1=NX
    N2=NY
    D9980: DO N=NCN(NY,NX),3
      D9975: DO NN=1,2
        IF(N.EQ.1)THEN
          IF(NN.EQ.1)THEN
            IF(NX.EQ.NHE)THEN
              N4=NX+1
              N5=NY
              N6=L
              XN=-1.0_r8
            ELSE
              cycle
            ENDIF
          ELSEIF(NN.EQ.2)THEN
            IF(NX.EQ.NHW)THEN
              N4=NX
              N5=NY
              N6=L
              XN=1.0_r8
            ELSE
              cycle
            ENDIF
          ENDIF
        ELSEIF(N.EQ.2)THEN
          IF(NN.EQ.1)THEN
            IF(NY.EQ.NVS)THEN
              N4=NX
              N5=NY+1
              N6=L
              XN=-1.0_r8
            ELSE
              cycle
            ENDIF
          ELSEIF(NN.EQ.2)THEN
            IF(NY.EQ.NVN)THEN
              N4=NX
              N5=NY
              N6=L
              XN=1.0_r8
            ELSE
              cycle
            ENDIF
          ENDIF
        ELSEIF(N.EQ.3)THEN
          IF(NN.EQ.1)THEN
            IF(L.EQ.NL(NY,NX))THEN
              N4=NX
              N5=NY
              N6=L+1
              XN=-1.0_r8
            ELSE
              cycle
            ENDIF
          ELSEIF(NN.EQ.2)THEN
            cycle
          ENDIF
        ENDIF
!
        call RunoffXBoundaryFluxes(L,N,NY,NX,N1,N2,N4,N5,NN,XN,CXR,ZXR,PXR,ZGR)
    !
        call SubsurfaceBoundaryFluxes(I,J,N,NY,NX,N1,N2,N3,N4,N5,N6,XN)
!
    !     WATER, HEAT, SOLUTES IN SNOW DRIFT
        call WaterHeatSoluteBySnowDrift(N,N4,N5,L,NY,NX,CXR,ZXR,PXR,ZGR,XN)
      ENDDO D9975
    ENDDO D9980
  ENDDO D9985
  end subroutine RunoffBal

!------------------------------------------------------------------------------------------

  subroutine RunoffXBoundaryFluxes(L,N,NY,NX,N1,N2,N4,N5,NN,XN,CXR,ZXR,PXR,ZGR)
  implicit none
  integer, intent(in) :: L,N,NY,NX
  integer, intent(in) :: N1,N2,N4,N5,NN
  real(r8), intent(in) :: XN
  real(r8), intent(out) :: CXR,ZXR,PXR,ZGR

  real(r8) :: SS1,SS2,SS3,SS4,SSR
  real(r8) :: COR,ZOR,POR
  real(r8) :: CXE,COE,COD
  real(r8) :: ECHY,ECOH,ECAL,ECFE,ECCA,ECMG,ECNA,ECKA
  real(r8) :: ECCO,ECHC,ECSO,ECCL,ECNO
  real(r8) :: ECNDQ
  real(r8) :: SEF,SEX,SEP,SET
  real(r8) :: PXE,PPE,POE
  real(r8) :: ZXE,ZPE,ZOE
  real(r8) :: WX,HGR,PSS
  real(r8) :: WQRN,HQRN
  real(r8) :: OXR,ER
  integer :: K,M,NO,NGL

!     begin_execution
!     RUNOFF BOUNDARY FLUXES OF WATER AND HEAT
!
!     QR,QS,QW,QI=runoff from surface water, snowpack snow,water,ice from watsub.f
!     CRUN,URUN=cumulative water and snow runoff
!     HEATOU=cumulative heat loss through lateral and lower boundaries
!
  IF(N.NE.3.AND.L.EQ.NU(NY,NX))THEN
    WQRN=XN*QR(N,NN,N5,N4)
    WQRH(N2,N1)=WQRH(N2,N1)+WQRN
    IF(ABS(WQRN).GT.ZEROS(N5,N4))THEN
      CRUN=CRUN-WQRN
      URUN(NY,NX)=URUN(NY,NX)-WQRN
      HQRN=XN*HQR(N,NN,N5,N4)
      HEATOU=HEATOU-HQRN
!
!     RUNOFF BOUNDARY FLUXES OF C, N AND P
!
!     X*QRS,X*QSS=solute in runoff, snow drift from trnsfr.f
!     solute code:CO=CO2,CH=CH4,OX=O2,NG=N2,N2=N2O,HG=H2
!             :OC=DOC,OA=acetate,ON=DON,OP=DOP
!             :N4=NH4,N3=NH3,NO=NO3,NX=NO2,PI=HPO4,P4=H2PO4 in non-band
!     XN=direction indicator
!     TCOU,OXYGOU,H2GOU,TZOU,TPOU=cumulative C,O2,H2,N,P loss through lateral and lower boundaries
!     UDOCQ,UDICQ=dissolved organic,inorganic C loss through runoff
!     UDONQ,UDINQ=dissolved organic,inorganic N loss through runoff
!     UDOPQ,UDIPQ=dissolved organic,inorganic P loss through runoff
!
      CXR=XN*(trcg_XRS(idg_CO2,N,NN,N5,N4)+trcg_XRS(idg_CH4,N,NN,N5,N4))
      ZXR=XN*(trcn_XRS(ids_NH4,N,NN,N5,N4)+trcg_XRS(idg_NH3,N,NN,N5,N4) &
        +trcn_XRS(ids_NO3,N,NN,N5,N4)+trcn_XRS(ids_NO2,N,NN,N5,N4))
      ZGR=XN*(trcg_XRS(idg_N2O,N,NN,N5,N4)+trcg_XRS(idg_N2,N,NN,N5,N4))
      PXR=XN*(trcn_XRS(ids_H2PO4,N,NN,N5,N4)+trcn_XRS(ids_H1PO4,N,NN,N5,N4))
      COR=0.0_r8
      ZOR=0.0_r8
      POR=0.0_r8
      D2575: DO K=1,jcplx
        COR=COR+XN*(XOCQRS(K,N,NN,N5,N4)+XOAQRS(K,N,NN,N5,N4))
        ZOR=ZOR+XN*XONQRS(K,N,NN,N5,N4)
        POR=POR+XN*XOPQRS(K,N,NN,N5,N4)
      ENDDO D2575
      TCOU=TCOU-CXR-COR
      TZOU=TZOU-ZXR-ZOR-ZGR
      TPOU=TPOU-PXR-POR
      UDOCQ(NY,NX)=UDOCQ(NY,NX)-COR
      UDICQ(NY,NX)=UDICQ(NY,NX)-CXR
      UDONQ(NY,NX)=UDONQ(NY,NX)-ZOR
      UDINQ(NY,NX)=UDINQ(NY,NX)-ZXR-ZGR
      UDOPQ(NY,NX)=UDOPQ(NY,NX)-POR
      UDIPQ(NY,NX)=UDIPQ(NY,NX)-PXR
      OXR=XN*trcg_XRS(idg_O2,N,NN,N5,N4)
      OXYGOU=OXYGOU-OXR
      HGR=XN*trcg_XRS(idg_H2,N,NN,N5,N4)
      H2GOU=H2GOU+HGR
!
!     RUNOFF BOUNDARY FLUXES OF SOLUTES
!
!     ISALTG=salt flag
!     XQR*,XQS*=solute loss in runoff,snow drift from trnsfrs.f
!     salt code: *HY*=H+,*OH*=OH-,*AL*=Al3+,*FE*=Fe3+,*CA*=Ca2+,*MG*=Mg2+
!          :*NA*=Na+,*KA*=K+,*SO4*=SO42-,*CL*=Cl-,*CO3*=CO32-,*HCO3*=HCO3-
!          :*CO2*=CO2,*ALO1*=AlOH2-,*ALOH2=AlOH2-,*ALOH3*=AlOH3
!          :*ALOH4*=AlOH4+,*ALS*=AlSO4+,*FEO1*=FeOH2-,*FEOH2=F3OH2-
!          :*FEOH3*=FeOH3,*FEOH4*=FeOH4+,*FES*=FeSO4+,*CAO*=CaOH
!          :*CAC*=CaCO3,*CAH*=CaHCO3-,*CAS*=CaSO4,*MGO*=MgOH,*MGC*=MgCO3
!          :*MHG*=MgHCO3-,*MGS*=MgSO4,*NAC*=NaCO3-,*NAS*=NaSO4-,*KAS*=KSO4-
!     phosphorus code: *H0P*=PO43-,*H3P*=H3PO4,*F1P*=FeHPO42-,*F2P*=F1H2PO4-
!          :*C0P*=CaPO4-,*C1P*=CaHPO4,*C2P*=CaH2PO4+,*M1P*=MgHPO4,*COO*=COOH-
!          :*1=non-band,*B=band
!     XN=direction indicator
!     TPOU,TIONOU=total P,salt loss through lateral and lower boundaries
!
      IF(salt_model)THEN
        PSS=XN*patomw*(trcsa_XQR(idsa_H0PO4,N,NN,N5,N4) &
          +trcsa_XQR(idsa_CaPO4,N,NN,N5,N4)+trcsa_XQR(idsa_FeHPO4,N,NN,N5,N4) &
          +trcsa_XQR(idsa_CaHPO4,N,NN,N5,N4) &
          +trcsa_XQR(idsa_MgHPO4,N,NN,N5,N4)+trcsa_XQR(idsa_H3PO4,N,NN,N5,N4) &
          +trcsa_XQR(idsa_FeH2PO4,N,NN,N5,N4) &
          +trcsa_XQR(idsa_CaH2PO4,N,NN,N5,N4))
        SS1=XN*(trcsa_XQR(idsa_Al,N,NN,N5,N4)+trcsa_XQR(idsa_Fe,N,NN,N5,N4) &
          +trcsa_XQR(idsa_Hp,N,NN,N5,N4) &
          +trcsa_XQR(idsa_Ca,N,NN,N5,N4)+trcsa_XQR(idsa_Mg,N,NN,N5,N4) &
          +trcsa_XQR(idsa_Na,N,NN,N5,N4) &
          +trcsa_XQR(idsa_K,N,NN,N5,N4)+trcsa_XQR(idsa_OH,N,NN,N5,N4) &
          +trcsa_XQR(idsa_SO4,N,NN,N5,N4) &
          +trcsa_XQR(idsa_Cl,N,NN,N5,N4)+trcsa_XQR(idsa_CO3,N,NN,N5,N4) &
          +trcsa_XQR(idsa_H0PO4,N,NN,N5,N4))
        SS2=XN*2.0*(trcsa_XQR(idsa_HCO3,N,NN,N5,N4)+trcsa_XQR(idsa_AlOH,N,NN,N5,N4) &
          +trcsa_XQR(idsa_AlSO4,N,NN,N5,N4)+trcsa_XQR(idsa_FeOH,N,NN,N5,N4) &
          +trcsa_XQR(idsa_FeSO4,N,NN,N5,N4) &
          +trcsa_XQR(idsa_CaOH2,N,NN,N5,N4)+trcsa_XQR(idsa_CaCO3,N,NN,N5,N4) &
          +trcsa_XQR(idsa_CaSO4,N,NN,N5,N4) &
          +trcsa_XQR(idsa_MgOH2,N,NN,N5,N4)+trcsa_XQR(idsa_MgCO3,N,NN,N5,N4) &
          +trcsa_XQR(idsa_MgSO4,N,NN,N5,N4) &
          +trcsa_XQR(idsa_NaCO3,N,NN,N5,N4)+trcsa_XQR(idsa_NaSO4,N,NN,N5,N4) &
          +trcsa_XQR(idsa_KSO4,N,NN,N5,N4)+trcsa_XQR(idsa_CaPO4,N,NN,N5,N4))
        SS3=XN*3.0*(trcsa_XQR(idsa_AlOH2,N,NN,N5,N4)+trcsa_XQR(idsa_FeOH2,N,NN,N5,N4) &
          +trcsa_XQR(idsa_CaHCO3,N,NN,N5,N4)+trcsa_XQR(idsa_MgHCO3,N,NN,N5,N4)+trcsa_XQR(idsa_FeHPO4,N,NN,N5,N4) &
          +trcsa_XQR(idsa_CaHPO4,N,NN,N5,N4)+trcsa_XQR(idsa_MgHPO4,N,NN,N5,N4))
        SS4=XN*4.0*(trcsa_XQR(idsa_AlOH3,N,NN,N5,N4)+trcsa_XQR(idsa_FeOH3,N,NN,N5,N4) &
          +trcsa_XQR(idsa_H3PO4,N,NN,N5,N4)+trcsa_XQR(idsa_FeH2PO4,N,NN,N5,N4)+trcsa_XQR(idsa_CaH2PO4,N,NN,N5,N4)) &
          +XN*5.0*(trcsa_XQR(idsa_AlOH4,N,NN,N5,N4)+trcsa_XQR(idsa_FeOH4,N,NN,N5,N4))
        PSS=PSS+XN*patomw*trcsa_XQS(idsa_H0PO4,N,N5,N4)
        TPOU=TPOU-PSS
        SSR=SS1+SS2+SS3+SS4
        TIONOU=TIONOU-SSR
        UIONOU(NY,NX)=UIONOU(NY,NX)-SSR
!     WRITE(20,3336)'SSR',I,J,N,N5,N4,SSR,SS1,SS2,SS3,SS4,TIONOU
!3336  FORMAT(A8,5I6,20F16.9)
!
!       SURFACE RUNOFF ELECTRICAL CONDUCTIVITY
!
!       QR=surface runoff from watsub.f
!       XQR*=solute in runoff from trnsfrs.f
!       ECNDQ=electrical conductivity of runoff
!
        WX=QR(N,NN,N5,N4)
        IF(ABS(WX).GT.ZEROS(N5,N4))THEN
          ECHY=0.337*AZMAX1(trcsa_XQR(idsa_Hp,N,NN,N5,N4)/WX)
          ECOH=0.192*AZMAX1(trcsa_XQR(idsa_OH,N,NN,N5,N4)/WX)
          ECAL=0.056*AZMAX1(trcsa_XQR(idsa_Al,N,NN,N5,N4)*3.0/WX)
          ECFE=0.051*AZMAX1(trcsa_XQR(idsa_Fe,N,NN,N5,N4)*3.0/WX)
          ECCA=0.060*AZMAX1(trcsa_XQR(idsa_Ca,N,NN,N5,N4)*2.0/WX)
          ECMG=0.053*AZMAX1(trcsa_XQR(idsa_Mg,N,NN,N5,N4)*2.0/WX)
          ECNA=0.050*AZMAX1(trcsa_XQR(idsa_Na,N,NN,N5,N4)/WX)
          ECKA=0.070*AZMAX1(trcsa_XQR(idsa_K,N,NN,N5,N4)/WX)
          ECCO=0.072*AZMAX1(trcsa_XQR(idsa_CO3,N,NN,N5,N4)*2.0/WX)
          ECHC=0.044*AZMAX1(trcsa_XQR(idsa_HCO3,N,NN,N5,N4)/WX)
          ECSO=0.080*AZMAX1(trcsa_XQR(idsa_SO4,N,NN,N5,N4)*2.0/WX)
          ECCL=0.076*AZMAX1(trcsa_XQR(idsa_Cl,N,NN,N5,N4)/WX)
          ECNO=0.071*AZMAX1(trcn_XRS(ids_NO3,N,NN,N5,N4)/(WX*natomw))
          ECNDQ=ECHY+ECOH+ECAL+ECFE+ECCA+ECMG+ECNA+ECKA &
            +ECCO+ECHC+ECSO+ECCL+ECNO
!     WRITE(*,9991)'ECNDQ',IYRC,I,J,N4,N5,N,NN,WX,ECNDQ
!9991  FORMAT(A8,7I4,2E12.4)
        ELSE
          ECNDQ=0.0_r8
        ENDIF
      ENDIF
!
!     RUNOFF BOUNDARY FLUXES OF SEDIMENT FROM EROSION
!
!     IERSNG=erosion flag
!     *ER=sediment flux from erosion.f
!     sediment code:XSED=total,XSAN=sand,XSIL=silt,XCLA=clay
!     TSEDOU,USEDOU=cumulative sediment loss through lateral and lower boundaries
!
      IF(N.NE.3.AND.IERSNG.EQ.1.OR.IERSNG.EQ.3)THEN
        IF(ABS(XSEDER(N,NN,N5,N4)).GT.ZEROS(N5,N4))THEN
          ER=XN*XSEDER(N,NN,N5,N4)
          TSEDOU=TSEDOU-ER
          USEDOU(NY,NX)=USEDOU(NY,NX)-ER
!
!
!         RUNOFF BOUNDARY FLUXES OF ORGANIC MATTER FROM EROSION

!         *ER=sediment flux from erosion.f
!         sediment code:OMC,OMN,OMP=microbial C,N,P; ORC=microbial residue C,N,P
!                      :OHC,OHN,OHP=adsorbed C,N,P; OSC,OSN,OSP=humus C,N,P
!         TSEDOU,USEDOU=cumulative sediment loss through lateral and lower boundaries
!         UDOCQ,UDICQ=dissolved organic,inorganic C loss through lateral and lower boundaries
!         UDONQ,UDINQ=dissolved organic,inorganic N loss through lateral and lower boundaries
!         UDOPQ,UDIPQ=dissolved organic,inorganic P loss through lateral and lower boundaries
!         TCOU,TZOU,TPOU=total C,N,P loss through lateral and lower boundaries

!         MICROBIAL C IN RUNOFF SEDIMENT
!
          CXE=0.0_r8
          ZXE=XN*natomw*(trcx_XER(idx_NH4,N,NN,N5,N4)+trcx_XER(idx_NH4B,N,NN,N5,N4))
          ZPE=XN*natomw*(XNH4ER(N,NN,N5,N4)+XNH3ER(N,NN,N5,N4) &
            +XNHUER(N,NN,N5,N4)+XNO3ER(N,NN,N5,N4)+XNH4EB(N,NN,N5,N4) &
            +XNH3EB(N,NN,N5,N4)+XNHUEB(N,NN,N5,N4)+XNO3EB(N,NN,N5,N4))
          PXE=XN*patomw*(trcx_XER(idx_HPO4,N,NN,N5,N4)+trcx_XER(idx_H2PO4,N,NN,N5,N4) &
            +trcx_XER(idx_HPO4B,N,NN,N5,N4)+trcx_XER(idx_H2PO4B,N,NN,N5,N4))
          PPE=XN*patomw*(1._r8*(trcp_ER(idsp_AlPO4,N,NN,N5,N4)+trcp_ER(idsp_FePO4,N,NN,N5,N4) &
            +trcp_ER(idsp_CaHPO4,N,NN,N5,N4)+trcp_ER(idsp_AlPO4B,N,NN,N5,N4) &
            +trcp_ER(idsp_FePO4B,N,NN,N5,N4)+trcp_ER(idsp_CaHPO4B,N,NN,N5,N4)) &
            +2.0_r8*(trcp_ER(idsp_CaH2PO4,N,NN,N5,N4)+trcp_ER(idsp_CaH2PO4B,N,NN,N5,N4)) &
            +3.0_r8*(trcp_ER(idsp_HA,N,NN,N5,N4)+trcp_ER(idsp_HAB,N,NN,N5,N4)))
          COE=0.0_r8
          ZOE=0.0_r8
          POE=0.0_r8
          D3580: DO K=1,jcplx
            DO NO=1,NFGs
              DO M=1,nlbiomcp
                DO NGL=JGnio(NO),JGnfo(NO)
                  COE=COE+XN*OMCER(M+(NGL-1)*nlbiomcp,K,N,NN,N5,N4)
                  ZOE=ZOE+XN*OMNER(M+(NGL-1)*nlbiomcp,K,N,NN,N5,N4)
                  POE=POE+XN*OMPER(M+(NGL-1)*nlbiomcp,K,N,NN,N5,N4)
                enddo
              enddo
            enddo
          ENDDO D3580
          DO NO=1,NFGs
            DO M=1,nlbiomcp
              DO NGL=JGniA(NO),JGnfA(NO)
                COE=COE+XN*OMCERff(M+(NGL-1)*nlbiomcp,N,NN,N5,N4)
                ZOE=ZOE+XN*OMNERff(M+(NGL-1)*nlbiomcp,N,NN,N5,N4)
                POE=POE+XN*OMPERff(M+(NGL-1)*nlbiomcp,N,NN,N5,N4)
              enddo
            enddo
          enddo

!
    !     MICROBIAL RESIDUE C IN RUNOFF SEDIMENT
!
          D3575: DO K=1,jcplx
            D3570: DO M=1,ndbiomcp
              COE=COE+XN*ORCER(M,K,N,NN,N5,N4)
              ZOE=ZOE+XN*ORNER(M,K,N,NN,N5,N4)
              POE=POE+XN*ORPER(M,K,N,NN,N5,N4)
            ENDDO D3570
!
        !   DOC, ADSORBED AND HUMUS C IN RUNOFF SEDIMENT
!
            COE=COE+XN*(OHCER(K,N,NN,N5,N4)+OHAER(K,N,NN,N5,N4))
            ZOE=ZOE+XN*OHNER(K,N,NN,N5,N4)
            POE=POE+XN*OHPER(K,N,NN,N5,N4)
            D3565: DO M=1,jsken
              COE=COE+XN*OSCER(M,K,N,NN,N5,N4)
              ZOE=ZOE+XN*OSNER(M,K,N,NN,N5,N4)
              POE=POE+XN*OSPER(M,K,N,NN,N5,N4)
            ENDDO D3565
          ENDDO D3575
          TCOU=TCOU-COE-CXE
          TZOU=TZOU-ZOE-ZXE-ZPE
          TPOU=TPOU-POE-PXE-PPE
          UDOCQ(NY,NX)=UDOCQ(NY,NX)-COE
          UDICQ(NY,NX)=UDICQ(NY,NX)-CXE
          UDONQ(NY,NX)=UDONQ(NY,NX)-ZOE
          UDINQ(NY,NX)=UDINQ(NY,NX)-ZXE-ZPE
          UDOPQ(NY,NX)=UDOPQ(NY,NX)-POE
          UDIPQ(NY,NX)=UDIPQ(NY,NX)-PXE-PPE
!     WRITE(*,6635)'POE',I,J,N4,N5,N,NN
!    2,COE,CXE,ZOE,ZXE,ZPE
!    3,POE,PXE,PPE,TPOU,XSEDER(N,NN,N5,N4)
!    3,XN,TCOU,TZOU,TPOU
!6635  FORMAT(A8,6I4,20F17.8)
!
!         ADSORBED AND PRECIPITATED SALTS IN RUNOFF SEDIMENTS

!         ISALTG=salt flag
!         *ER=sediment flux from erosion.f
!         sediment code
!           :NH4,NH3,NHU,NO3=fertilizer NH4,NH3,urea,NO3 in non-band
!           :NH4B,NH3B,NHUB,NO3B=fertilizer NH4,NH3,urea,NO3 in band
!           :XN4,XNB=adsorbed NH4 in non-band,band
!           :XHY,XAL,XFE,XCA,XMG,XNA,XKA,XHC,AL2,FE2
!            =adsorbed H,Al,Fe,Ca,Mg,Na,K,HCO3,AlOH2,FeOH2
!           :XOH0,XOH1,XOH2=adsorbed R-,R-OH,R-OH2 in non-band
!           :XOH0B,XOH1B,XOH2B=adsorption sites R-,R-OH,R-OH2 in band
!           :XH1P,XH2P=adsorbed HPO4,H2PO4 in non-band
!           :XH1PB,XP2PB=adsorbed HPO4,H2PO4 in band
!           :PALO,PFEO=precip AlOH,FeOH
!           :PCAC,PCAS=precip CaCO3,CaSO4
!           :PALP,PFEP=precip AlPO4,FEPO4 in non-band
!           :PALPB,PFEPB=precip AlPO4,FEPO4 in band
!           :PCPM,PCPD,PCPH=precip CaH2PO4,CaHPO4,apatite in non-band
!           :PCPMB,PCPDB,PCPHB=precip CaH2PO4,CaHPO4,apatite in band
!         TIONOU,UIONOU=total salt loss through lateral and lower boundaries
!
          IF(salt_model)THEN
            SEF=XN*(XNH3ER(N,NN,N5,N4)+XNHUER(N,NN,N5,N4)+XNO3ER(N,NN,N5,N4) &
              +XNH3EB(N,NN,N5,N4)+XNHUEB(N,NN,N5,N4)+XNO3EB(N,NN,N5,N4)) &
              +2.0*(XNH4ER(N,NN,N5,N4)+XNH4EB(N,NN,N5,N4))
            SEX=XN*(trcx_XER(idx_Hp,N,NN,N5,N4)+trcx_XER(idx_Al,N,NN,N5,N4) &
              +trcx_XER(idx_Fe,N,NN,N5,N4)+trcx_XER(idx_Ca,N,NN,N5,N4)+trcx_XER(idx_Mg,N,NN,N5,N4) &
              +trcx_XER(idx_Na,N,NN,N5,N4)+trcx_XER(idx_K,N,NN,N5,N4)+trcx_XER(idx_COOH,N,NN,N5,N4) &
              +trcx_XER(idx_OHe,N,NN,N5,N4)+trcx_XER(idx_OHeB,N,NN,N5,N4)) &
              +XN*2.0*(trcx_XER(idx_NH4,N,NN,N5,N4)+trcx_XER(idx_NH4B,N,NN,N5,N4) &
              +trcx_XER(idx_OH,N,NN,N5,N4)+trcx_XER(idx_OHB,N,NN,N5,N4)) &
              +XN*3.0*(trcx_XER(idx_AlOH2,N,NN,N5,N4)+trcx_XER(idx_FeOH2,N,NN,N5,N4) &
              +trcx_XER(idx_OHp,N,NN,N5,N4)+trcx_XER(idx_OHpB,N,NN,N5,N4) &
              +trcx_XER(idx_HPO4,N,NN,N5,N4)+trcx_XER(idx_HPO4B,N,NN,N5,N4)) &
              +XN*4.0*(trcx_XER(idx_H2PO4,N,NN,N5,N4)+trcx_XER(idx_H2PO4B,N,NN,N5,N4))
            SEP=XN*2.0*(trcp_ER(idsp_CaCO3,N,NN,N5,N4)+trcp_ER(idsp_CaSO4,N,NN,N5,N4) &
              +trcp_ER(idsp_AlPO4,N,NN,N5,N4)+trcp_ER(idsp_FePO4,N,NN,N5,N4) &
              +trcp_ER(idsp_AlPO4B,N,NN,N5,N4)+trcp_ER(idsp_FePO4B,N,NN,N5,N4)) &
              +XN*3.0*(trcp_ER(idsp_CaHPO4,N,NN,N5,N4)+trcp_ER(idsp_CaHPO4B,N,NN,N5,N4)) &
              +XN*4.0*(trcp_ER(idsp_AlOH3,N,NN,N5,N4)+trcp_ER(idsp_FeOH3,N,NN,N5,N4)) &
              +XN*7.0*(trcp_ER(idsp_CaH2PO4,N,NN,N5,N4)+trcp_ER(idsp_CaH2PO4B,N,NN,N5,N4)) &
              +XN*9.0*(trcp_ER(idsp_HA,N,NN,N5,N4)+trcp_ER(idsp_HAB,N,NN,N5,N4))
            SET=SEF+SEX+SEP
            TIONOU=TIONOU-SET
            UIONOU(NY,NX)=UIONOU(NY,NX)-SET
!     WRITE(*,3342)'SET',I,J,N4,N5,NN,N,SET,SEF,SEX,SEP,TIONOU
!3342  FORMAT(A8,6I4,12F18.6)
          ENDIF
        ENDIF
      ENDIF
    ENDIF
  ENDIF
  end subroutine RunoffXBoundaryFluxes
!------------------------------------------------------------------------------------------

  subroutine SubsurfaceBoundaryFluxes(I,J,N,NY,NX,N1,N2,N3,N4,N5,N6,XN)
  implicit none
  integer, intent(in) :: I,J,N,NY,NX
  integer, intent(in) :: N1,N2,N3,N4,N5,N6
  real(r8), intent(in) :: XN
  real(r8) :: ECHY,ECOH,ECAL,ECFE,ECCA,ECMG,ECNA,ECKA
  real(r8) :: ECCO,ECHC,ECSO,ECCL,ECNO
  real(r8) :: ECNDX
  real(r8) :: WO,SO
  real(r8) :: COD,ZOD,POD,HOD,OOD
  real(r8) :: CXD,ZXD,PXD,ZGD
  real(r8) :: SSD,SHD,PHD,PQD
  real(r8) :: WX
  integer :: K
!     begin_execution
!     SUBSURFACE BOUNDARY FLUXES OF WATER AND HEAT
!
!     FLW,FLWH,HFLW=micropore,macropore,heat flux through lateral and lower boundaries from watsub.f
!     VOLWOU,HEATOU=cumulative water, heat loss through lateral and lower boundaries
!     UVOLO,HVOLO=cumulative,hourly water loss through lateral and lower boundaries
!
  IF(NCN(NY,NX).NE.3.OR.N.EQ.3)THEN
    HEATOU=HEATOU-XN*HFLW(N,N6,N5,N4)
    WO=XN*(FLW(N,N6,N5,N4)+FLWH(N,N6,N5,N4))
    IF(abs(WO)>0._r8)THEN
      VOLWOU=VOLWOU-WO
      HVOLO(N2,N1)=HVOLO(N2,N1)-WO
      UVOLO(N2,N1)=UVOLO(N2,N1)-WO
!     IF((I/10)*10.EQ.I.AND.J.EQ.15)THEN
!     WRITE(*,3488)'UVOLO',I,J,N6,N5,N4,N,XN,WO
!    2,UVOLO(NY,NX),FLW(N,N6,N5,N4),FLWH(N,N6,N5,N4)
!3488  FORMAT(A8,6I4,12E12.4)
!     ENDIF
!
!     SUBSURFACE BOUNDARY FLUXES OF CO2 AND DOC
!
!     X*FLS,X*FHS=solute flux in macropores,micropores from trnsfr.f
!     X*FLG=convective+diffusive gas flux from trnsfr.f
!     TCOU=cumulative C loss through lateral and lower boundaries
!     UDOCD,UDICD=dissolved organic,inorganic C loss through subsurface boundaries
!
!     SUBSURFACE BOUNDARY FLUXES OF N2O, N2, NH4, NH3, NO3, NO2 AND DON
!
!     X*FLS,X*FHS,X*FLB,X*FHB=solute flux in macropores,micropores in non-band,band from trnsfr.f
!     X*FLG=convective+diffusive gas flux from trnsfr.f
!     TZOU=cumulative N loss through lateral and lower boundaries
!     UDOND,UDIND=dissolved organic,inorganic N loss through subsurface boundaries
!
!     SUBSURFACE BOUNDARY FLUXES OF PO4 AND DOP
!
!     X*FLS,X*FHS,X*FLB,X*FHB=solute flux in macropores,micropores in non-band,band from trnsfr.f
!     TPOU=cumulative P loss through lateral and lower boundaries
!     UDOPD,UDIPD=dissolved organic,inorganic P loss through subsurface boundaries
!
      COD=0.0_r8
      ZOD=0.0_r8
      POD=0.0_r8
      D450: DO K=1,jcplx
        COD=COD+XN*(XOCFLS(K,N,N6,N5,N4)+XOAFLS(K,N,N6,N5,N4) &
          +XOCFHS(K,N,N6,N5,N4)+XOAFHS(K,N,N6,N5,N4))
        ZOD=ZOD+XN*(XONFLS(K,N,N6,N5,N4)+XONFHS(K,N,N6,N5,N4))
        POD=POD+XN*(XOPFLS(K,N,N6,N5,N4)+XOPFHS(K,N,N6,N5,N4))
      ENDDO D450
      CXD=XN*(trcs_XFLS(idg_CO2,N,N6,N5,N4)+trcs_XFHS(idg_CO2,N,N6,N5,N4) &
        +R3GasADTFlx(idg_CO2,N,N6,N5,N4)+trcs_XFLS(idg_CH4,N,N6,N5,N4) &
        +trcs_XFHS(idg_CH4,N,N6,N5,N4)+R3GasADTFlx(idg_CH4,N,N6,N5,N4))
      ZXD=XN*(trcs_XFLS(ids_NH4,N,N6,N5,N4)+trcs_XFLS(idg_NH3,N,N6,N5,N4) &
        +trcs_XFLS(ids_NO3,N,N6,N5,N4) &
        +trcs_XFLS(ids_NH4B,N,N6,N5,N4)+trcs_XFLS(idg_NH3B,N,N6,N5,N4)&
        +trcs_XFLS(ids_NO3B,N,N6,N5,N4) &
        +trcs_XFLS(ids_NO2,N,N6,N5,N4)+trcs_XFLS(ids_NO2B,N,N6,N5,N4) &
        +trcs_XFHS(ids_NH4,N,N6,N5,N4)+trcs_XFHS(idg_NH3,N,N6,N5,N4) &
        +trcs_XFHS(ids_NO3,N,N6,N5,N4) &
        +trcs_XFHS(ids_NH4B,N,N6,N5,N4)+trcs_XFHS(idg_NH3B,N,N6,N5,N4) &
        +trcs_XFHS(ids_NO3B,N,N6,N5,N4) &
        +trcs_XFHS(ids_NO2,N,N6,N5,N4)+trcs_XFHS(ids_NO2B,N,N6,N5,N4))
      ZGD=XN*(trcs_XFLS(idg_N2,N,N6,N5,N4)+R3GasADTFlx(idg_N2,N,N6,N5,N4) &
        +trcs_XFHS(idg_N2,N,N6,N5,N4) &
        +trcs_XFLS(idg_N2O,N,N6,N5,N4)+R3GasADTFlx(idg_N2O,N,N6,N5,N4) &
        +trcs_XFHS(idg_N2O,N,N6,N5,N4) &
        +R3GasADTFlx(idg_NH3,N,N6,N5,N4))
      PXD=XN*(trcs_XFLS(ids_H2PO4,N,N6,N5,N4)+trcs_XFLS(ids_H2PO4B,N,N6,N5,N4) &
        +trcs_XFHS(ids_H2PO4,N,N6,N5,N4)+trcs_XFHS(ids_H2PO4B,N,N6,N5,N4)&
        +trcs_XFLS(ids_H1PO4,N,N6,N5,N4) &
        +trcs_XFLS(ids_H1PO4B,N,N6,N5,N4)+trcs_XFHS(ids_H1PO4,N,N6,N5,N4) &
        +trcs_XFHS(ids_H1PO4B,N,N6,N5,N4))

      TCOU=TCOU-COD-CXD
      TZOU=TZOU-ZOD-ZXD-ZGD
      TPOU=TPOU-POD-PXD
      UDOCD(N2,N1)=UDOCD(N2,N1)-COD
      UDICD(N2,N1)=UDICD(N2,N1)-CXD
      UDOND(N2,N1)=UDOND(N2,N1)-ZOD
      UDIND(N2,N1)=UDIND(N2,N1)-ZXD
      UDOPD(N2,N1)=UDOPD(N2,N1)-POD
      UDIPD(N2,N1)=UDIPD(N2,N1)-PXD
!
!     SUBSURFACE BOUNDARY FLUXES OF O2
!
!     X*FLS,X*FHS=solute flux in macropores,micropores from trnsfr.f
!     X*FLG=convective+diffusive gas flux from trnsfr.f
!     OXYGOU,H2GOU=cumulative O2,H2 loss through lateral and lower boundaries
!
      OOD=XN*(trcs_XFLS(idg_O2,N,N6,N5,N4)+trcs_XFHS(idg_O2,N,N6,N5,N4)+R3GasADTFlx(idg_O2,N,N6,N5,N4))
      OXYGOU=OXYGOU-OOD
      HOD=XN*(trcs_XFLS(idg_H2,N,N6,N5,N4)+trcs_XFHS(idg_H2,N,N6,N5,N4)+R3GasADTFlx(idg_H2,N,N6,N5,N4))
      H2GOU=H2GOU-HOD
!
!     SUBSURFACE BOUNDARY FLUXES OF SOLUTES
!
!     X*FLS=hourly convective + diffusive solute flux through micropores from trnsfrs.f
!     X*FLW,X*FLB= hourly convective + diffusive solute flux through micropores in non-band,band from trnsfrs.f
!     X*FHS=hourly convective + diffusive solute flux through macropores from trnsfrs.f
!     X*FHW,X*FHB= hourly convective + diffusive solute flux through macropores in non-band,band from trnsfrs.f
!     TIONOU,UIONOU=total salt loss through lateral and lower boundaries
!
      IF(salt_model)THEN
        PQD=XN*patomw*(trcsa_XFLS(idsa_H0PO4,N,N6,N5,N4)+trcsa_XFLS(idsa_H0PO4B,N,N6,N5,N4) &
          +trcsa_XFLS(idsa_CaPO4,N,N6,N5,N4)+trcsa_XFLS(idsa_CaPO4B,N,N6,N5,N4) &
          +trcsa_XFLS(idsa_FeHPO4,N,N6,N5,N4) &
          +trcsa_XFLS(idsa_CaHPO4,N,N6,N5,N4)+trcsa_XFLS(idsa_MgHPO4,N,N6,N5,N4) &
          +trcsa_XFLS(idsa_FeHPO4B,N,N6,N5,N4) &
          +trcsa_XFLS(idsa_CaHPO4B,N,N6,N5,N4)+trcsa_XFLS(idsa_MgHPO4B,N,N6,N5,N4) &
          +trcsa_XFLS(idsa_H3PO4,N,N6,N5,N4) &
          +trcsa_XFLS(idsa_FeH2PO4,N,N6,N5,N4)+trcsa_XFLS(idsa_CaH2PO4,N,N6,N5,N4) &
          +trcsa_XFLS(idsa_H3PO4B,N,N6,N5,N4) &
          +trcsa_XFLS(idsa_FeH2PO4B,N,N6,N5,N4)+trcsa_XFLS(idsa_CaH2PO4B,N,N6,N5,N4))

        PHD=XN*patomw*(trcsa_XFHS(idsa_H0PO4,N,N6,N5,N4)+trcsa_XFHS(idsa_H0PO4B,N,N6,N5,N4) &
          +trcsa_XFHS(idsa_CaPO4,N,N6,N5,N4)+trcsa_XFHS(idsa_CaPO4B,N,N6,N5,N4) &
          +trcsa_XFHS(idsa_FeHPO4,N,N6,N5,N4) &
          +trcsa_XFHS(idsa_CaHPO4,N,N6,N5,N4)+trcsa_XFHS(idsa_MgHPO4,N,N6,N5,N4) &
          +trcsa_XFHS(idsa_FeHPO4B,N,N6,N5,N4) &
          +trcsa_XFHS(idsa_CaHPO4B,N,N6,N5,N4)+trcsa_XFHS(idsa_MgHPO4B,N,N6,N5,N4) &
          +trcsa_XFHS(idsa_H3PO4,N,N6,N5,N4) &
          +trcsa_XFHS(idsa_FeH2PO4,N,N6,N5,N4)+trcsa_XFHS(idsa_CaH2PO4,N,N6,N5,N4) &
          +trcsa_XFHS(idsa_H3PO4B,N,N6,N5,N4) &
          +trcsa_XFHS(idsa_FeH2PO4B,N,N6,N5,N4)+trcsa_XFHS(idsa_CaH2PO4B,N,N6,N5,N4))

        TPOU=TPOU-PQD-PHD
        SSD=XN*(trcsa_XFLS(idsa_Al,N,N6,N5,N4)+trcsa_XFLS(idsa_Fe,N,N6,N5,N4) &
          +trcsa_XFLS(idsa_Hp,N,N6,N5,N4) &
          +trcsa_XFLS(idsa_Ca,N,N6,N5,N4)+trcsa_XFLS(idsa_Mg,N,N6,N5,N4) &
          +trcsa_XFLS(idsa_Na,N,N6,N5,N4) &
          +trcsa_XFLS(idsa_K,N,N6,N5,N4)+trcsa_XFLS(idsa_OH,N,N6,N5,N4) &
          +trcsa_XFLS(idsa_SO4,N,N6,N5,N4) &
          +trcsa_XFLS(idsa_Cl,N,N6,N5,N4)+trcsa_XFLS(idsa_CO3,N,N6,N5,N4) &
          +trcsa_XFLS(idsa_H0PO4,N,N6,N5,N4) &
          +trcsa_XFLS(idsa_H0PO4B,N,N6,N5,N4)+2.0*(trcsa_XFLS(idsa_HCO3,N,N6,N5,N4) &
          +trcsa_XFLS(idsa_AlOH,N,N6,N5,N4) &
          +trcsa_XFLS(idsa_AlSO4,N,N6,N5,N4)+trcsa_XFLS(idsa_FeOH,N,N6,N5,N4) &
          +trcsa_XFLS(idsa_FeSO4,N,N6,N5,N4) &
          +trcsa_XFLS(idsa_CaOH2,N,N6,N5,N4)+trcsa_XFLS(idsa_CaCO3,N,N6,N5,N4) &
          +trcsa_XFLS(idsa_CaSO4,N,N6,N5,N4)+trcsa_XFLS(idsa_MgOH2,N,N6,N5,N4) &
          +trcsa_XFLS(idsa_MgCO3,N,N6,N5,N4) &
          +trcsa_XFLS(idsa_MgSO4,N,N6,N5,N4)+trcsa_XFLS(idsa_NaCO3,N,N6,N5,N4) &
          +trcsa_XFLS(idsa_NaSO4,N,N6,N5,N4) &
          +trcsa_XFLS(idsa_KSO4,N,N6,N5,N4)+trcsa_XFLS(idsa_CaPO4,N,N6,N5,N4) &
          +trcsa_XFLS(idsa_CaPO4B,N,N6,N5,N4)) &
          +3.0*(trcsa_XFLS(idsa_AlOH2,N,N6,N5,N4) &
          +trcsa_XFLS(idsa_FeOH2,N,N6,N5,N4)+trcsa_XFLS(idsa_CaHCO3,N,N6,N5,N4) &
          +trcsa_XFLS(idsa_MgHCO3,N,N6,N5,N4) &
          +trcsa_XFLS(idsa_FeHPO4,N,N6,N5,N4)+trcsa_XFLS(idsa_CaHPO4,N,N6,N5,N4) &
          +trcsa_XFLS(idsa_MgHPO4,N,N6,N5,N4) &
          +trcsa_XFLS(idsa_FeHPO4B,N,N6,N5,N4)+trcsa_XFLS(idsa_CaHPO4B,N,N6,N5,N4) &
          +trcsa_XFLS(idsa_MgHPO4B,N,N6,N5,N4)) &
          +4.0*(trcsa_XFLS(idsa_AlOH3,N,N6,N5,N4)+trcsa_XFLS(idsa_FeOH3,N,N6,N5,N4) &
          +trcsa_XFLS(idsa_H3PO4,N,N6,N5,N4) &
          +trcsa_XFLS(idsa_FeH2PO4,N,N6,N5,N4)+trcsa_XFLS(idsa_CaH2PO4,N,N6,N5,N4) &
          +trcsa_XFLS(idsa_H3PO4B,N,N6,N5,N4) &
          +trcsa_XFLS(idsa_FeH2PO4B,N,N6,N5,N4)+trcsa_XFLS(idsa_CaH2PO4B,N,N6,N5,N4)) &
          +5.0*(trcsa_XFLS(idsa_AlOH4,N,N6,N5,N4)+trcsa_XFLS(idsa_FeOH4,N,N6,N5,N4)))

        SHD=XN*(trcsa_XFHS(idsa_Al,N,N6,N5,N4)+trcsa_XFHS(idsa_Fe,N,N6,N5,N4) &
          +trcsa_XFHS(idsa_Hp,N,N6,N5,N4) &
          +trcsa_XFHS(idsa_Ca,N,N6,N5,N4)+trcsa_XFHS(idsa_Mg,N,N6,N5,N4) &
          +trcsa_XFHS(idsa_Na,N,N6,N5,N4) &
          +trcsa_XFHS(idsa_K,N,N6,N5,N4)+trcsa_XFHS(idsa_OH,N,N6,N5,N4) &
          +trcsa_XFHS(idsa_SO4,N,N6,N5,N4) &
          +trcsa_XFHS(idsa_Cl,N,N6,N5,N4)+trcsa_XFHS(idsa_CO3,N,N6,N5,N4) &
          +trcsa_XFHS(idsa_H0PO4,N,N6,N5,N4) &
          +trcsa_XFHS(idsa_H0PO4B,N,N6,N5,N4) &
          +2.0*(trcsa_XFHS(idsa_HCO3,N,N6,N5,N4)+trcsa_XFHS(idsa_AlOH,N,N6,N5,N4) &
          +trcsa_XFHS(idsa_AlSO4,N,N6,N5,N4)+trcsa_XFHS(idsa_FeOH,N,N6,N5,N4) &
          +trcsa_XFHS(idsa_FeSO4,N,N6,N5,N4) &
          +trcsa_XFHS(idsa_CaOH2,N,N6,N5,N4)+trcsa_XFHS(idsa_CaCO3,N,N6,N5,N4) &
          +trcsa_XFHS(idsa_CaSO4,N,N6,N5,N4)+trcsa_XFHS(idsa_MgOH2,N,N6,N5,N4) &
          +trcsa_XFHS(idsa_MgCO3,N,N6,N5,N4) &
          +trcsa_XFHS(idsa_MgSO4,N,N6,N5,N4)+trcsa_XFHS(idsa_NaCO3,N,N6,N5,N4) &
          +trcsa_XFHS(idsa_NaSO4,N,N6,N5,N4) &
          +trcsa_XFHS(idsa_KSO4,N,N6,N5,N4)+trcsa_XFHS(idsa_CaPO4,N,N6,N5,N4) &
          +trcsa_XFHS(idsa_CaPO4B,N,N6,N5,N4)) &
          +3.0*(trcsa_XFHS(idsa_AlOH2,N,N6,N5,N4) &
          +trcsa_XFHS(idsa_FeOH2,N,N6,N5,N4)+trcsa_XFHS(idsa_CaHCO3,N,N6,N5,N4) &
          +trcsa_XFHS(idsa_MgHCO3,N,N6,N5,N4) &
          +trcsa_XFHS(idsa_FeHPO4,N,N6,N5,N4)+trcsa_XFHS(idsa_CaHPO4,N,N6,N5,N4) &
          +trcsa_XFHS(idsa_MgHPO4,N,N6,N5,N4) &
          +trcsa_XFHS(idsa_FeHPO4B,N,N6,N5,N4)+trcsa_XFHS(idsa_CaHPO4B,N,N6,N5,N4) &
          +trcsa_XFHS(idsa_MgHPO4B,N,N6,N5,N4)) &
          +4.0*(trcsa_XFHS(idsa_AlOH3,N,N6,N5,N4)+trcsa_XFHS(idsa_FeOH3,N,N6,N5,N4) &
          +trcsa_XFHS(idsa_H3PO4,N,N6,N5,N4) &
          +trcsa_XFHS(idsa_FeH2PO4,N,N6,N5,N4)+trcsa_XFHS(idsa_CaH2PO4,N,N6,N5,N4) &
          +trcsa_XFHS(idsa_H3PO4B,N,N6,N5,N4) &
          +trcsa_XFHS(idsa_FeH2PO4B,N,N6,N5,N4)+trcsa_XFHS(idsa_CaH2PO4B,N,N6,N5,N4)) &
          +5.0*(trcsa_XFHS(idsa_AlOH4,N,N6,N5,N4)+trcsa_XFHS(idsa_AlOH4,N,N6,N5,N4)))

        SO=SSD+SHD
        TIONOU=TIONOU-SO
        UIONOU(N2,N1)=UIONOU(N2,N1)-SO

!
!     SUBSURFACE FLUX ELECTRICAL CONDUCTIVITY
!
!     FLW,FLWH=micropore,macropore flux through lateral and lower boundaries from watsub.f
!     X*FLS=hourly convective + diffusive solute flux through micropores from trnsfrs.f
!     X*FLW,X*FLB= hourly convective + diffusive solute flux through micropores in non-band,band from trnsfrs.f
!     X*FHS=hourly convective + diffusive solute flux through macropores from trnsfrs.f
!     X*FHW,X*FHB= hourly convective + diffusive solute flux through macropores in non-band,band from trnsfrs.f
!     ECNDQ=electrical conductivity of water flux
!
        WX=FLW(N,N6,N5,N4)+FLWH(N,N6,N5,N4)
        IF(ABS(WX).GT.ZEROS(N2,N1))THEN
          ECHY=0.337*AZMAX1((trcsa_XFLS(idsa_Hp,N,N6,N5,N4)+trcsa_XFHS(idsa_Hp,N,N6,N5,N4))/WX)
          ECOH=0.192*AZMAX1((trcsa_XFLS(idsa_OH,N,N6,N5,N4)+trcsa_XFHS(idsa_OH,N,N6,N5,N4))/WX)
          ECAL=0.056*AZMAX1((trcsa_XFLS(idsa_Al,N,N6,N5,N4)+trcsa_XFHS(idsa_Ca,N,N6,N5,N4))*3.0/WX)
          ECFE=0.051*AZMAX1((trcsa_XFLS(idsa_Fe,N,N6,N5,N4)+trcsa_XFHS(idsa_Fe,N,N6,N5,N4))*3.0/WX)
          ECCA=0.060*AZMAX1((trcsa_XFLS(idsa_Ca,N,N6,N5,N4)+trcsa_XFHS(idsa_Ca,N,N6,N5,N4))*2.0/WX)
          ECMG=0.053*AZMAX1((trcsa_XFLS(idsa_Mg,N,N6,N5,N4)+trcsa_XFHS(idsa_Mg,N,N6,N5,N4))*2.0/WX)
          ECNA=0.050*AZMAX1((trcsa_XFLS(idsa_Na,N,N6,N5,N4)+trcsa_XFHS(idsa_Na,N,N6,N5,N4))/WX)
          ECKA=0.070*AZMAX1((trcsa_XFLS(idsa_K,N,N6,N5,N4)+trcsa_XFHS(idsa_K,N,N6,N5,N4))/WX)
          ECCO=0.072*AZMAX1((trcsa_XFLS(idsa_CO3,N,N6,N5,N4)+trcsa_XFHS(idsa_CO3,N,N6,N5,N4))*2.0/WX)
          ECHC=0.044*AZMAX1((trcsa_XFLS(idsa_HCO3,N,N6,N5,N4)+trcsa_XFHS(idsa_HCO3,N,N6,N5,N4))/WX)
          ECSO=0.080*AZMAX1((trcsa_XFLS(idsa_SO4,N,N6,N5,N4)+trcsa_XFHS(idsa_SO4,N,N6,N5,N4))*2.0/WX)
          ECCL=0.076*AZMAX1((trcsa_XFLS(idsa_Cl,N,N6,N5,N4)+trcsa_XFHS(idsa_Cl,N,N6,N5,N4))/WX)
          ECNO=0.071*AZMAX1((trcs_XFLS(ids_NO3,N,N6,N5,N4)+trcs_XFHS(ids_NO3,N,N6,N5,N4))/(WX*natomw))
          ECNDX=ECHY+ECOH+ECAL+ECFE+ECCA+ECMG+ECNA+ECKA+ECCO+ECHC+ECSO+ECCL+ECNO
!     IF((I/10)*10.EQ.I.AND.J.EQ.15)THEN
!     WRITE(*,9992)'ECNDX',IYRC,I,J,N4,N5,N6,N,WX,ECNDX
!    2,FLW(N,N6,N5,N4),FLWH(N,N6,N5,N4)
!9992  FORMAT(A8,7I4,4E12.4)
!     ENDIF
        ELSE
          ECNDX=0.0_r8
        ENDIF
      ENDIF
      SG=SG+trcs_XFLS(idg_H2,N,N6,N5,N4)+R3GasADTFlx(idg_H2,N,N6,N5,N4)
    ENDIF
  ENDIF
  end subroutine SubsurfaceBoundaryFluxes

!------------------------------------------------------------------------------------------

  subroutine WaterHeatSoluteBySnowDrift(N,N4,N5,L,NY,NX,CXR,ZXR,PXR,ZGR,XN)
  implicit none
  integer, intent(in) :: N,L,NY,NX,N4,N5
  REAL(R8), INTENT(IN) :: CXR,ZXR,PXR,ZGR
  real(r8), intent(in) :: XN
  real(r8) :: CXS,ZXS,ZGS,PXS
  real(r8) :: SS1,SS2,SS3,SS4,SSR
  real(r8) :: WQRS,HQRS,OXS,PSS
!     begin_execution
  IF(N.NE.3.AND.L.EQ.NU(NY,NX))THEN
    WQRS=XN*(QS(N,N5,N4)+QW(N,N5,N4)+QI(N,N5,N4))
    IF(ABS(WQRS).GT.ZEROS(N5,N4))THEN
      CRUN=CRUN-WQRS
      URUN(NY,NX)=URUN(NY,NX)-WQRS
      HQRS=XN*HQS(N,N5,N4)
      HEATOU=HEATOU-HQRS
      CXS=XN*(XCOQSS(N,N5,N4)+XCHQSS(N,N5,N4))
      ZXS=XN*(XN4QSS(N,N5,N4)+XN3QSS(N,N5,N4)+XNOQSS(N,N5,N4))
      ZGS=XN*(XN2QSS(N,N5,N4)+XNGQSS(N,N5,N4))
      PXS=XN*(XP4QSS(N,N5,N4)+XP1QSS(N,N5,N4))
      TCOU=TCOU-CXS
      TZOU=TZOU-ZXS-ZGS
      TPOU=TPOU-PXS
      UDICQ(NY,NX)=UDICQ(NY,NX)-CXR
      UDINQ(NY,NX)=UDINQ(NY,NX)-ZXR-ZGR
      UDIPQ(NY,NX)=UDIPQ(NY,NX)-PXR
      OXS=XN*XOXQSS(N,N5,N4)
      OXYGOU=OXYGOU-OXS
      IF(salt_model)THEN
        PSS=XN*patomw*(trcsa_XQS(idsa_CaPO4,N,N5,N4) &
          +trcsa_XQS(idsa_FeHPO4,N,N5,N4)+trcsa_XQS(idsa_CaHPO4,N,N5,N4) &
          +trcsa_XQS(idsa_MgHPO4,N,N5,N4)+trcsa_XQS(idsa_H3PO4,N,N5,N4) &
          +trcsa_XQS(idsa_FeH2PO4,N,N5,N4)+trcsa_XQS(idsa_CaH2PO4,N,N5,N4))
        TPOU=TPOU-PSS
        SS1=XN*(trcsa_XQS(idsa_Al,N,N5,N4)+trcsa_XQS(idsa_Fe,N,N5,N4) &
          +trcsa_XQS(idsa_Hp,N,N5,N4)+trcsa_XQS(idsa_Ca,N,N5,N4) &
          +trcsa_XQS(idsa_Mg,N,N5,N4)+trcsa_XQS(idsa_Na,N,N5,N4) &
          +trcsa_XQS(idsa_K,N,N5,N4)+trcsa_XQS(idsa_OH,N,N5,N4) &
          +trcsa_XQS(idsa_SO4,N,N5,N4)+trcsa_XQS(idsa_Cl,N,N5,N4) &
          +trcsa_XQS(idsa_CO3,N,N5,N4)+trcsa_XQS(idsa_H0PO4,N,N5,N4))
        SS2=XN*2.0*(trcsa_XQS(idsa_HCO3,N,N5,N4)+trcsa_XQS(idsa_AlOH,N,N5,N4) &
          +trcsa_XQS(idsa_AlSO4,N,N5,N4) &
          +trcsa_XQS(idsa_FeOH,N,N5,N4)+trcsa_XQS(idsa_FeSO4,N,N5,N4) &
          +trcsa_XQS(idsa_CaOH2,N,N5,N4)+trcsa_XQS(idsa_CaCO3,N,N5,N4) &
          +trcsa_XQS(idsa_CaSO4,N,N5,N4)+trcsa_XQS(idsa_MgOH2,N,N5,N4)&
          +trcsa_XQS(idsa_MgCO3,N,N5,N4)+trcsa_XQS(idsa_MgSO4,N,N5,N4) &
          +trcsa_XQS(idsa_NaCO3,N,N5,N4)+trcsa_XQS(idsa_NaSO4,N,N5,N4) &
          +trcsa_XQS(idsa_KSO4,N,N5,N4)+trcsa_XQS(idsa_CaPO4,N,N5,N4))
        SS3=XN*3.0*(trcsa_XQS(idsa_AlOH2,N,N5,N4)+trcsa_XQS(idsa_FeOH2,N,N5,N4)&
          +trcsa_XQS(idsa_CaHCO3,N,N5,N4) &
          +trcsa_XQS(idsa_MgHCO3,N,N5,N4)+trcsa_XQS(idsa_FeHPO4,N,N5,N4) &
          +trcsa_XQS(idsa_CaHPO4,N,N5,N4)+trcsa_XQS(idsa_MgHPO4,N,N5,N4))
        SS4=XN*4.0*(trcsa_XQS(idsa_AlOH3,N,N5,N4)+trcsa_XQS(idsa_FeOH3,N,N5,N4) &
          +trcsa_XQS(idsa_H3PO4,N,N5,N4)+trcsa_XQS(idsa_FeH2PO4,N,N5,N4)&
          +trcsa_XQS(idsa_CaH2PO4,N,N5,N4)) &
          +XN*5.0*(trcsa_XQS(idsa_AlOH4,N,N5,N4)+trcsa_XQS(idsa_FeOH4,N,N5,N4))
        SSR=SS1+SS2+SS3+SS4
        TIONOU=TIONOU-SSR
        UIONOU(NY,NX)=UIONOU(NY,NX)-SSR
      ENDIF
    ENDIF
  ENDIF
  end subroutine WaterHeatSoluteBySnowDrift

end module RunoffBalMod
