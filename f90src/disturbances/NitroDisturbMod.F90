module NitroDisturbMod

! USES:
  use data_kind_mod, only : r8 => SHR_KIND_R8
  use abortutils  , only : endrun
  use minimathmod, only : safe_adb
  use MicrobialDataType
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
  use PlantMngmtDataType
  use EcoSimSumDataType
  use EcosimBGCFluxType
  use SoilBGCDataType
  use AqueChemDatatype
  use GridDataType
  implicit none

  private

  character(len=*), parameter :: mod_filename = __FILE__

  public :: SOMRemovalByDisturbance
  contains
!------------------------------------------------------------------------------------------

  subroutine SOMRemovalByDisturbance(I,J,NY,NX)
!
!     Description:
!
  implicit none
  integer, intent(in) :: I,J,NY,NX

  integer :: L,K,M,N,IFLGJ,NLL,NGL
  real(r8) :: DC,DN,DP
  real(r8) :: DCORPC
  real(r8) :: FORGCX
  real(r8) :: HFLXD
  real(r8) :: OC,ON,OP,OCH,ONH,OPH,ONX
  REAL(R8) :: OPX,OCA,OAH
  real(r8) :: ONL(4,0:jcplx1),OPL(4,0:jcplx1)
!     begin_execution

  IF(J.EQ.INT(ZNOON(NY,NX)).AND.(ITILL(I,NY,NX).EQ.21 &
    .OR.ITILL(I,NY,NX).EQ.22))THEN
    IF(ITILL(I,NY,NX).EQ.22)THEN
      IFLGS(NY,NX)=1
      IFLGJ=0
      NLL=-1
      DO 2945 L=0,NL(NY,NX)

        IF(L.EQ.0.OR.L.GE.NUM(NY,NX))THEN
          IF(IFLGJ.EQ.1)THEN
            exit
          ELSEIF(THETW(L,NY,NX).GT.FVLWB.OR.CORGC(L,NY,NX).LE.FORGC)THEN
            IFLGJ=1
          ELSE
            NLL=L
          ENDIF
        ENDIF
2945  CONTINUE
      ELSE
        NLL=0
      ENDIF
      DO 2950 L=0,NLL
        IF(NLL.GE.0)THEN
          IF(ITILL(I,NY,NX).EQ.22)THEN
            IF(L.EQ.0)THEN
              FORGCX=0.0_r8
            ELSE
              FORGCX=FORGC
            ENDIF
            DCORPC=AMIN1(0.999,DCORP(I,NY,NX))*(CORGC(L,NY,NX)-FORGCX) &
              /(AMAX1(CORGC(L,NY,NX),0.55E+06)-FORGCX)
          ELSE
            DCORPC=AMIN1(0.999,DCORP(I,NY,NX))
          ENDIF
!     VOLWOU=VOLWOU+DCORPC*VOLW(L,NY,NX)
!     HEATOU=HEATOU+DCORPC*4.19*TKS(L,NY,NX)*VOLW(L,NY,NX)
!     VOLW(L,NY,NX)=VOLW(L,NY,NX)-DCORPC*VOLW(L,NY,NX)
        OC=0.0_r8
        ON=0.0_r8
        OP=0.0_r8
        DC=0.0_r8
        DN=0.0_r8
        DP=0.0_r8
        ONL(1:4,0:jcplx1)=0._r8
        OPL(1:4,0:jcplx1)=0._r8

        DO 2970 K=0,4
          IF(L.NE.0.OR.(K.NE.3.AND.K.NE.4))THEN
!
!     REMOVE MICROBIAL BIOMASS
!
            DO 2960 N=1,7
              DO NGL=1,JG
                DO M=1,3
                  OCH=DCORPC*OMC(M,NGL,N,K,L,NY,NX)
                  ONH=DCORPC*OMN(M,NGL,N,K,L,NY,NX)
                  OPH=DCORPC*OMP(M,NGL,N,K,L,NY,NX)
                  ONX=EFIRE(1,ITILL(I,NY,NX))*ONH
                  OPX=EFIRE(2,ITILL(I,NY,NX))*OPH
                  IF(K.LE.2)THEN
                    ONL(4,K)=ONL(4,K)+ONH-ONX
                    OPL(4,K)=OPL(4,K)+OPH-OPX
                  ELSEIF(K.LE.4)THEN
                    ONL(1,K)=ONL(1,K)+ONH-ONX
                    OPL(1,K)=OPL(1,K)+OPH-OPX
                  ENDIF
                  OMC(M,NGL,N,K,L,NY,NX)=OMC(M,NGL,N,K,L,NY,NX)-OCH
                  OMN(M,NGL,N,K,L,NY,NX)=OMN(M,NGL,N,K,L,NY,NX)-ONH
                  OMP(M,NGL,N,K,L,NY,NX)=OMP(M,NGL,N,K,L,NY,NX)-OPH
                  DC=DC+OMC(M,NGL,N,K,L,NY,NX)
                  DN=DN+OMN(M,NGL,N,K,L,NY,NX)
                  DP=DP+OMP(M,NGL,N,K,L,NY,NX)
                  OC=OC+OCH
                  ON=ON+ONX
                  OP=OP+OPX
                enddo
              enddo
2960        CONTINUE
          ENDIF
2970    CONTINUE


!          IF(L.NE.0.OR.(K.NE.3.AND.K.NE.4))THEN
!
!     REMOVE MICROBIAL BIOMASS
!
            DO  N=1,7
              DO NGL=1,JG
                DO M=1,3
                  OCH=DCORPC*OMCff(M,NGL,N,L,NY,NX)
                  ONH=DCORPC*OMNff(M,NGL,N,L,NY,NX)
                  OPH=DCORPC*OMPff(M,NGL,N,L,NY,NX)
                  ONX=EFIRE(1,ITILL(I,NY,NX))*ONH
                  OPX=EFIRE(2,ITILL(I,NY,NX))*OPH
                  ONL(4,1)=ONL(4,1)+ONH-ONX
                  OPL(4,1)=OPL(4,1)+OPH-OPX
                  OMCff(M,NGL,N,L,NY,NX)=OMCff(M,NGL,N,L,NY,NX)-OCH
                  OMNff(M,NGL,N,L,NY,NX)=OMNff(M,NGL,N,L,NY,NX)-ONH
                  OMPff(M,NGL,N,L,NY,NX)=OMPff(M,NGL,N,L,NY,NX)-OPH
                  DC=DC+OMCff(M,NGL,N,L,NY,NX)
                  DN=DN+OMNff(M,NGL,N,L,NY,NX)
                  DP=DP+OMPff(M,NGL,N,L,NY,NX)
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
        DO 2900 K=0,4
          IF(L.NE.0.OR.(K.NE.3.AND.K.NE.4))THEN
            DO 2940 M=1,2
              OCH=DCORPC*ORC(M,K,L,NY,NX)
              ONH=DCORPC*ORN(M,K,L,NY,NX)
              OPH=DCORPC*ORP(M,K,L,NY,NX)
              ONX=EFIRE(1,ITILL(I,NY,NX))*ONH
              OPX=EFIRE(2,ITILL(I,NY,NX))*OPH
              IF(K.LE.2)THEN
                ONL(4,K)=ONL(4,K)+ONH-ONX
                OPL(4,K)=OPL(4,K)+OPH-OPX
              ELSE
                ONL(1,K)=ONL(1,K)+ONH-ONX
                OPL(1,K)=OPL(1,K)+OPH-OPX
              ENDIF
              ORC(M,K,L,NY,NX)=ORC(M,K,L,NY,NX)-OCH
              ORN(M,K,L,NY,NX)=ORN(M,K,L,NY,NX)-ONH
              ORP(M,K,L,NY,NX)=ORP(M,K,L,NY,NX)-OPH
              DC=DC+ORC(M,K,L,NY,NX)
              DN=DN+ORN(M,K,L,NY,NX)
              DP=DP+ORP(M,K,L,NY,NX)
              OC=OC+OCH
              ON=ON+ONX
              OP=OP+OPX
2940        CONTINUE
!
!     REMOVE DOC, DON, DOP
!
            OCH=DCORPC*OQC(K,L,NY,NX)
            OCA=DCORPC*OQA(K,L,NY,NX)
            ONH=DCORPC*OQN(K,L,NY,NX)
            OPH=DCORPC*OQP(K,L,NY,NX)
            ONX=EFIRE(1,ITILL(I,NY,NX))*ONH
            OPX=EFIRE(2,ITILL(I,NY,NX))*OPH
            IF(K.LE.2)THEN
              ONL(4,K)=ONL(4,K)+ONH-ONX
              OPL(4,K)=OPL(4,K)+OPH-OPX
            ELSE
              ONL(1,K)=ONL(1,K)+ONH-ONX
              OPL(1,K)=OPL(1,K)+OPH-OPX
            ENDIF
            OQC(K,L,NY,NX)=OQC(K,L,NY,NX)-OCH
            OQA(K,L,NY,NX)=OQA(K,L,NY,NX)-OCA
            OQN(K,L,NY,NX)=OQN(K,L,NY,NX)-ONH
            OQP(K,L,NY,NX)=OQP(K,L,NY,NX)-OPH
            OC=OC+OCH+OCA
            ON=ON+ONX
            OP=OP+OPX
            OCH=DCORPC*OQCH(K,L,NY,NX)
            ONH=DCORPC*OQNH(K,L,NY,NX)
            OPH=DCORPC*OQPH(K,L,NY,NX)
            OAH=DCORPC*OQAH(K,L,NY,NX)
            ONX=EFIRE(1,ITILL(I,NY,NX))*ONH
            OPX=EFIRE(2,ITILL(I,NY,NX))*OPH
            IF(K.LE.2)THEN
              ONL(4,K)=ONL(4,K)+ONH-ONX
              OPL(4,K)=OPL(4,K)+OPH-OPX
            ELSE
              ONL(1,K)=ONL(1,K)+ONH-ONX
              OPL(1,K)=OPL(1,K)+OPH-OPX
            ENDIF
            OQCH(K,L,NY,NX)=OQCH(K,L,NY,NX)-OCH
            OQNH(K,L,NY,NX)=OQNH(K,L,NY,NX)-ONH
            OQPH(K,L,NY,NX)=OQPH(K,L,NY,NX)-OPH
            OQAH(K,L,NY,NX)=OQAH(K,L,NY,NX)-OAH
            OC=OC+OCH+OAH
            ON=ON+ONX
            OP=OP+OPX
!
!     REMOVE ADSORBED OM
!
            OCH=DCORPC*OHC(K,L,NY,NX)
            ONH=DCORPC*OHN(K,L,NY,NX)
            OPH=DCORPC*OHP(K,L,NY,NX)
            OAH=DCORPC*OHA(K,L,NY,NX)
            ONX=EFIRE(1,ITILL(I,NY,NX))*ONH
            OPX=EFIRE(2,ITILL(I,NY,NX))*OPH
            IF(K.LE.2)THEN
              ONL(4,K)=ONL(4,K)+ONH-ONX
              OPL(4,K)=OPL(4,K)+OPH-OPX
            ELSE
              ONL(1,K)=ONL(1,K)+ONH-ONX
              OPL(1,K)=OPL(1,K)+OPH-OPX
            ENDIF
            OHC(K,L,NY,NX)=OHC(K,L,NY,NX)-OCH
            OHN(K,L,NY,NX)=OHN(K,L,NY,NX)-ONH
            OHP(K,L,NY,NX)=OHP(K,L,NY,NX)-OPH
            OHA(K,L,NY,NX)=OHA(K,L,NY,NX)-OAH
            DC=DC+OQC(K,L,NY,NX)+OQCH(K,L,NY,NX)+OHC(K,L,NY,NX) &
              +OQA(K,L,NY,NX)+OQAH(K,L,NY,NX)+OHA(K,L,NY,NX)
            DN=DN+OQN(K,L,NY,NX)+OQNH(K,L,NY,NX)+OHN(K,L,NY,NX)
            DP=DP+OQP(K,L,NY,NX)+OQPH(K,L,NY,NX)+OHP(K,L,NY,NX)
            OC=OC+OCH
            ON=ON+ONX
            OP=OP+OPX
!
!     REMOVE RESIDUE
!
            DO 2930 M=1,jsken
              OCH=DCORPC*OSC(M,K,L,NY,NX)
              OCA=DCORPC*OSA(M,K,L,NY,NX)
              ONH=DCORPC*OSN(M,K,L,NY,NX)
              OPH=DCORPC*OSP(M,K,L,NY,NX)
              ONX=EFIRE(1,ITILL(I,NY,NX))*ONH
              OPX=EFIRE(2,ITILL(I,NY,NX))*OPH
              ONL(M,K)=ONL(M,K)+ONH-ONX
              OPL(M,K)=OPL(M,K)+OPH-OPX
              OSC(M,K,L,NY,NX)=OSC(M,K,L,NY,NX)-OCH
              OSA(M,K,L,NY,NX)=OSA(M,K,L,NY,NX)-OCA
              OSN(M,K,L,NY,NX)=OSN(M,K,L,NY,NX)-ONH
              OSP(M,K,L,NY,NX)=OSP(M,K,L,NY,NX)-OPH
              DC=DC+OSC(M,K,L,NY,NX)
              DN=DN+OSN(M,K,L,NY,NX)
              DP=DP+OSP(M,K,L,NY,NX)
              OC=OC+OCH
              ON=ON+ONX
              OP=OP+OPX
2930        CONTINUE
          ENDIF
2900    CONTINUE
!
!     ADD UNBURNED N,P TO ORG N, ORG P
!
        DO 2905 K=0,4
          DO  M=1,jsken
            OSN(M,K,L,NY,NX)=OSN(M,K,L,NY,NX)+ONL(M,K)
            OSP(M,K,L,NY,NX)=OSP(M,K,L,NY,NX)+OPL(M,K)
            DN=DN+ONL(M,K)
            DP=DP+OPL(M,K)
          enddo
2905  CONTINUE
!
!     REMOVE FERTILIZER IN RESIDUE
!
      IF(ITILL(I,NY,NX).EQ.21)THEN
        ON=ON+DCORPC*(ZNH4S(L,NY,NX)+ZNH3S(L,NY,NX) &
          +ZNO3S(L,NY,NX)+ZNO2S(L,NY,NX))
        OP=OP+DCORPC*(H1PO4(L,NY,NX)+H2PO4(L,NY,NX))
        ZNH4S(L,NY,NX)=(1.0-DCORPC)*ZNH4S(L,NY,NX)
        ZNH3S(L,NY,NX)=(1.0-DCORPC)*ZNH3S(L,NY,NX)
        ZNO3S(L,NY,NX)=(1.0-DCORPC)*ZNO3S(L,NY,NX)
        ZNO2S(L,NY,NX)=(1.0-DCORPC)*ZNO2S(L,NY,NX)
        H1PO4(L,NY,NX)=(1.0-DCORPC)*H1PO4(L,NY,NX)
        H2PO4(L,NY,NX)=(1.0-DCORPC)*H2PO4(L,NY,NX)
        XN4(L,NY,NX)=(1.0-DCORPC)*XN4(L,NY,NX)
        PALPO(L,NY,NX)=(1.0-DCORPC)*PALPO(L,NY,NX)
        PFEPO(L,NY,NX)=(1.0-DCORPC)*PFEPO(L,NY,NX)
        PCAPD(L,NY,NX)=(1.0-DCORPC)*PCAPD(L,NY,NX)
        PCAPH(L,NY,NX)=(1.0-DCORPC)*PCAPH(L,NY,NX)
        PCAPM(L,NY,NX)=(1.0-DCORPC)*PCAPM(L,NY,NX)
        ZNH4FA(L,NY,NX)=(1.0-DCORPC)*ZNH4FA(L,NY,NX)
        ZNH3FA(L,NY,NX)=(1.0-DCORPC)*ZNH3FA(L,NY,NX)
        ZNHUFA(L,NY,NX)=(1.0-DCORPC)*ZNHUFA(L,NY,NX)
        ZNO3FA(L,NY,NX)=(1.0-DCORPC)*ZNO3FA(L,NY,NX)
      ENDIF
      ORGC(L,NY,NX)=DC
      ORGN(L,NY,NX)=DN
      IF(L.EQ.0)THEN
        HFLXD=4.19E-06*(ORGCX(L,NY,NX)-ORGC(L,NY,NX))*TKS(L,NY,NX)
        HEATOU=HEATOU+HFLXD
      ENDIF
!     IF(L.EQ.0)THEN
!     VHCP(0,NY,NX)=2.496E-06*ORGC(0,NY,NX)+4.19*VOLW(0,NY,NX)
!    2+1.9274*VOLI(0,NY,NX)
!     ELSE
!     VHCP(L,NY,NX)=VHCM(L,NY,NX)+4.19*(VOLW(L,NY,NX)+VOLWH(L,NY,NX))
!    2+1.9274*(VOLI(L,NY,NX)+VOLIH(L,NY,NX))
!     ENDIF
      IF(ITILL(I,NY,NX).EQ.21)THEN
        TCOU=TCOU+OC
        TZOU=TZOU+ON
        TPOU=TPOU+OP
        UDOCQ(NY,NX)=UDOCQ(NY,NX)+OC
        UDONQ(NY,NX)=UDONQ(NY,NX)+ON
        UDOPQ(NY,NX)=UDOPQ(NY,NX)+OP
        TNBP(NY,NX)=TNBP(NY,NX)-OC
      ELSEIF(ITILL(I,NY,NX).EQ.22)THEN
        CO2GIN=CO2GIN-OC
        OXYGIN=OXYGIN+2.667*OC
        OXYGOU=OXYGOU+2.667*OC
        TZOU=TZOU+ON
        TPOU=TPOU+OP
        UCO2F(NY,NX)=UCO2F(NY,NX)-(1.0-FCH4F)*OC
        UCH4F(NY,NX)=UCH4F(NY,NX)-FCH4F*OC
        UOXYF(NY,NX)=UOXYF(NY,NX)+(1.0-FCH4F)*2.667*OC
        UNH3F(NY,NX)=UNH3F(NY,NX)-ON
        UN2OF(NY,NX)=UN2OF(NY,NX)-0.0
        UPO4F(NY,NX)=UPO4F(NY,NX)-OP
        TNBP(NY,NX)=TNBP(NY,NX)-OC
      ENDIF
    ENDIF
2950  CONTINUE
  ENDIF
  end subroutine SOMRemovalByDisturbance
!------------------------------------------------------------------------------------------


end module NitroDisturbMod
