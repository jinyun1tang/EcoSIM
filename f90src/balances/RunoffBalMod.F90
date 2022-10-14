module RunoffBalMod
  use data_kind_mod, only : r8 => SHR_KIND_R8
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
      ENDDO D9975
!
    !     WATER, HEAT, SOLUTES IN SNOW DRIFT
      call WaterHeatSoluteBySnowDrift(N,N4,N5,L,NY,NX,CXR,ZXR,PXR,ZGR,XN)
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
      CXR=XN*(XCOQRS(N,NN,N5,N4)+XCHQRS(N,NN,N5,N4))
      ZXR=XN*(XN4QRW(N,NN,N5,N4)+XN3QRW(N,NN,N5,N4) &
        +XNOQRW(N,NN,N5,N4)+XNXQRS(N,NN,N5,N4))
      ZGR=XN*(XN2QRS(N,NN,N5,N4)+XNGQRS(N,NN,N5,N4))
      PXR=XN*(XP4QRW(N,NN,N5,N4)+XP1QRW(N,NN,N5,N4))
      COR=0.0_r8
      ZOR=0.0_r8
      POR=0.0_r8
      DO 2575 K=0,jcplx1
        COR=COR+XN*(XOCQRS(K,N,NN,N5,N4)+XOAQRS(K,N,NN,N5,N4))
        ZOR=ZOR+XN*XONQRS(K,N,NN,N5,N4)
        POR=POR+XN*XOPQRS(K,N,NN,N5,N4)
2575  CONTINUE
      TCOU=TCOU-CXR-COR
      TZOU=TZOU-ZXR-ZOR-ZGR
      TPOU=TPOU-PXR-POR
      UDOCQ(NY,NX)=UDOCQ(NY,NX)-COR
      UDICQ(NY,NX)=UDICQ(NY,NX)-CXR
      UDONQ(NY,NX)=UDONQ(NY,NX)-ZOR
      UDINQ(NY,NX)=UDINQ(NY,NX)-ZXR-ZGR
      UDOPQ(NY,NX)=UDOPQ(NY,NX)-POR
      UDIPQ(NY,NX)=UDIPQ(NY,NX)-PXR
      OXR=XN*XOXQRS(N,NN,N5,N4)
      OXYGOU=OXYGOU-OXR
      HGR=XN*XHGQRS(N,NN,N5,N4)
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
      IF(ISALTG.NE.0)THEN
        PSS=XN*31.0*(XQRH0P(N,NN,N5,N4) &
          +XQRC0P(N,NN,N5,N4)+XQRF1P(N,NN,N5,N4)+XQRC1P(N,NN,N5,N4) &
          +XQRM1P(N,NN,N5,N4)+XQRH3P(N,NN,N5,N4)+XQRF2P(N,NN,N5,N4) &
          +XQRC2P(N,NN,N5,N4))
        SS1=XN*(XQRAL(N,NN,N5,N4)+XQRFE(N,NN,N5,N4)+XQRHY(N,NN,N5,N4) &
          +XQRCA(N,NN,N5,N4)+XQRMG(N,NN,N5,N4)+XQRNA(N,NN,N5,N4) &
          +XQRKA(N,NN,N5,N4)+XQROH(N,NN,N5,N4)+XQRSO(N,NN,N5,N4) &
          +XQRCL(N,NN,N5,N4)+XQRC3(N,NN,N5,N4)+XQRH0P(N,NN,N5,N4))
        SS2=XN*2.0*(XQRHC(N,NN,N5,N4)+XQRAL1(N,NN,N5,N4) &
          +XQRALS(N,NN,N5,N4)+XQRFE1(N,NN,N5,N4)+XQRFES(N,NN,N5,N4) &
          +XQRCAO(N,NN,N5,N4)+XQRCAC(N,NN,N5,N4)+XQRCAS(N,NN,N5,N4) &
          +XQRMGO(N,NN,N5,N4)+XQRMGC(N,NN,N5,N4)+XQRMGS(N,NN,N5,N4) &
          +XQRNAC(N,NN,N5,N4)+XQRNAS(N,NN,N5,N4)+XQRKAS(N,NN,N5,N4) &
          +XQRC0P(N,NN,N5,N4))
        SS3=XN*3.0*(XQRAL2(N,NN,N5,N4)+XQRFE2(N,NN,N5,N4) &
          +XQRCAH(N,NN,N5,N4)+XQRMGH(N,NN,N5,N4)+XQRF1P(N,NN,N5,N4) &
          +XQRC1P(N,NN,N5,N4)+XQRM1P(N,NN,N5,N4))
        SS4=XN*4.0*(XQRAL3(N,NN,N5,N4)+XQRFE3(N,NN,N5,N4) &
          +XQRH3P(N,NN,N5,N4)+XQRF2P(N,NN,N5,N4)+XQRC2P(N,NN,N5,N4)) &
          +XN*5.0*(XQRAL4(N,NN,N5,N4)+XQRFE4(N,NN,N5,N4))
        PSS=PSS+XN*31.0*XQSH0P(N,N5,N4)
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
          ECHY=0.337*AZMAX1(XQRHY(N,NN,N5,N4)/WX)
          ECOH=0.192*AZMAX1(XQROH(N,NN,N5,N4)/WX)
          ECAL=0.056*AZMAX1(XQRAL(N,NN,N5,N4)*3.0/WX)
          ECFE=0.051*AZMAX1(XQRFE(N,NN,N5,N4)*3.0/WX)
          ECCA=0.060*AZMAX1(XQRCA(N,NN,N5,N4)*2.0/WX)
          ECMG=0.053*AZMAX1(XQRMG(N,NN,N5,N4)*2.0/WX)
          ECNA=0.050*AZMAX1(XQRNA(N,NN,N5,N4)/WX)
          ECKA=0.070*AZMAX1(XQRKA(N,NN,N5,N4)/WX)
          ECCO=0.072*AZMAX1(XQRC3(N,NN,N5,N4)*2.0/WX)
          ECHC=0.044*AZMAX1(XQRHC(N,NN,N5,N4)/WX)
          ECSO=0.080*AZMAX1(XQRSO(N,NN,N5,N4)*2.0/WX)
          ECCL=0.076*AZMAX1(XQRCL(N,NN,N5,N4)/WX)
          ECNO=0.071*AZMAX1(XNOQRW(N,NN,N5,N4)/(WX*14.0))
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
          ZXE=XN*14.0*(XN4ER(N,NN,N5,N4)+XNBER(N,NN,N5,N4))
          ZPE=XN*14.0*(XNH4ER(N,NN,N5,N4)+XNH3ER(N,NN,N5,N4) &
            +XNHUER(N,NN,N5,N4)+XNO3ER(N,NN,N5,N4)+XNH4EB(N,NN,N5,N4) &
            +XNH3EB(N,NN,N5,N4)+XNHUEB(N,NN,N5,N4)+XNO3EB(N,NN,N5,N4))
          PXE=XN*31.0*(XH1PER(N,NN,N5,N4)+XH2PER(N,NN,N5,N4) &
            +XH1PEB(N,NN,N5,N4)+XH2PEB(N,NN,N5,N4))
          PPE=XN*(31.0*(PALPER(N,NN,N5,N4)+PFEPER(N,NN,N5,N4) &
            +PCPDER(N,NN,N5,N4)+PALPEB(N,NN,N5,N4) &
            +PFEPEB(N,NN,N5,N4)+PCPDEB(N,NN,N5,N4)) &
            +62.0*(PCPMER(N,NN,N5,N4)+PCPMEB(N,NN,N5,N4)) &
            +93.0*(PCPHER(N,NN,N5,N4)+PCPHEB(N,NN,N5,N4)))
          COE=0.0_r8
          ZOE=0.0_r8
          POE=0.0_r8
          D3580: DO K=0,jcplx1
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
          D3575: DO K=0,jcplx1
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
          IF(ISALTG.NE.0)THEN
            SEF=XN*(XNH3ER(N,NN,N5,N4)+XNHUER(N,NN,N5,N4)+XNO3ER(N,NN,N5,N4) &
              +XNH3EB(N,NN,N5,N4)+XNHUEB(N,NN,N5,N4)+XNO3EB(N,NN,N5,N4)) &
              +2.0*(XNH4ER(N,NN,N5,N4)+XNH4EB(N,NN,N5,N4))
            SEX=XN*(XHYER(N,NN,N5,N4)+XALER(N,NN,N5,N4) &
              +XFEER(N,NN,N5,N4)+XCAER(N,NN,N5,N4)+XMGER(N,NN,N5,N4) &
              +XNAER(N,NN,N5,N4)+XKAER(N,NN,N5,N4)+XHCER(N,NN,N5,N4) &
              +XOH0ER(N,NN,N5,N4)+XOH0EB(N,NN,N5,N4)) &
              +XN*2.0*(XN4ER(N,NN,N5,N4)+XNBER(N,NN,N5,N4) &
              +XOH1ER(N,NN,N5,N4)+XOH1EB(N,NN,N5,N4)) &
              +XN*3.0*(XAL2ER(N,NN,N5,N4)+XFE2ER(N,NN,N5,N4) &
              +XOH2ER(N,NN,N5,N4)+XOH2EB(N,NN,N5,N4) &
              +XH1PER(N,NN,N5,N4)+XH1PEB(N,NN,N5,N4)) &
              +XN*4.0*(XH2PER(N,NN,N5,N4)+XH2PEB(N,NN,N5,N4))
            SEP=XN*2.0*(PCACER(N,NN,N5,N4)+PCASER(N,NN,N5,N4) &
              +PALPER(N,NN,N5,N4)+PFEPER(N,NN,N5,N4) &
              +PALPEB(N,NN,N5,N4)+PFEPEB(N,NN,N5,N4)) &
              +XN*3.0*(PCPDER(N,NN,N5,N4)+PCPDEB(N,NN,N5,N4)) &
              +XN*4.0*(PALOER(N,NN,N5,N4)+PFEOER(N,NN,N5,N4)) &
              +XN*7.0*(PCPMER(N,NN,N5,N4)+PCPMEB(N,NN,N5,N4)) &
              +XN*9.0*(PCPHER(N,NN,N5,N4)+PCPHEB(N,NN,N5,N4))
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
      DO 450 K=0,jcplx1
        COD=COD+XN*(XOCFLS(K,N,N6,N5,N4)+XOAFLS(K,N,N6,N5,N4) &
          +XOCFHS(K,N,N6,N5,N4)+XOAFHS(K,N,N6,N5,N4))
        ZOD=ZOD+XN*(XONFLS(K,N,N6,N5,N4)+XONFHS(K,N,N6,N5,N4))
        POD=POD+XN*(XOPFLS(K,N,N6,N5,N4)+XOPFHS(K,N,N6,N5,N4))
450   CONTINUE
      CXD=XN*(XCOFLS(N,N6,N5,N4)+XCOFHS(N,N6,N5,N4) &
        +XCOFLG(N,N6,N5,N4)+XCHFLS(N,N6,N5,N4) &
        +XCHFHS(N,N6,N5,N4)+XCHFLG(N,N6,N5,N4))
      ZXD=XN*(XN4FLW(N,N6,N5,N4)+XN3FLW(N,N6,N5,N4)+XNOFLW(N,N6,N5,N4) &
        +XN4FLB(N,N6,N5,N4)+XN3FLB(N,N6,N5,N4)+XNOFLB(N,N6,N5,N4) &
        +XNXFLS(N,N6,N5,N4)+XNXFLB(N,N6,N5,N4) &
        +XN4FHW(N,N6,N5,N4)+XN3FHW(N,N6,N5,N4)+XNOFHW(N,N6,N5,N4) &
        +XN4FHB(N,N6,N5,N4)+XN3FHB(N,N6,N5,N4)+XNOFHB(N,N6,N5,N4) &
        +XNXFHS(N,N6,N5,N4)+XNXFHB(N,N6,N5,N4))
      ZGD=XN*(XNGFLS(N,N6,N5,N4)+XNGFLG(N,N6,N5,N4)+XNGFHS(N,N6,N5,N4) &
        +XN2FLS(N,N6,N5,N4)+XN2FLG(N,N6,N5,N4)+XN2FHS(N,N6,N5,N4) &
        +XN3FLG(N,N6,N5,N4))
      PXD=XN*(XH2PFS(N,N6,N5,N4)+XH2BFB(N,N6,N5,N4) &
        +XH2PHS(N,N6,N5,N4)+XH2BHB(N,N6,N5,N4)+XH1PFS(N,N6,N5,N4) &
        +XH1BFB(N,N6,N5,N4)+XH1PHS(N,N6,N5,N4)+XH1BHB(N,N6,N5,N4))
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
      OOD=XN*(XOXFLS(N,N6,N5,N4)+XOXFHS(N,N6,N5,N4)+XOXFLG(N,N6,N5,N4))
      OXYGOU=OXYGOU-OOD
      HOD=XN*(XHGFLS(N,N6,N5,N4)+XHGFHS(N,N6,N5,N4)+XHGFLG(N,N6,N5,N4))
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
      IF(ISALTG.NE.0)THEN
        PQD=XN*31.0*(XH0PFS(N,N6,N5,N4)+XH0BFB(N,N6,N5,N4) &
          +XC0PFS(N,N6,N5,N4)+XC0BFB(N,N6,N5,N4)+XF1PFS(N,N6,N5,N4) &
          +XC1PFS(N,N6,N5,N4)+XM1PFS(N,N6,N5,N4)+XF1BFB(N,N6,N5,N4) &
          +XC1BFB(N,N6,N5,N4)+XM1BFB(N,N6,N5,N4)+XH3PFS(N,N6,N5,N4) &
          +XF2PFS(N,N6,N5,N4)+XC2PFS(N,N6,N5,N4)+XH3BFB(N,N6,N5,N4) &
          +XF2BFB(N,N6,N5,N4)+XC2BFB(N,N6,N5,N4))
        PHD=XN*31.0*(XH0PHS(N,N6,N5,N4)+XH0BHB(N,N6,N5,N4) &
          +XC0PHS(N,N6,N5,N4)+XC0BHB(N,N6,N5,N4)+XF1PHS(N,N6,N5,N4) &
          +XC1PHS(N,N6,N5,N4)+XM1PHS(N,N6,N5,N4)+XF1BHB(N,N6,N5,N4) &
          +XC1BHB(N,N6,N5,N4)+XM1BHB(N,N6,N5,N4)+XH3PHS(N,N6,N5,N4) &
          +XF2PHS(N,N6,N5,N4)+XC2PHS(N,N6,N5,N4)+XH3BHB(N,N6,N5,N4) &
          +XF2BHB(N,N6,N5,N4)+XC2BHB(N,N6,N5,N4))
        TPOU=TPOU-PQD-PHD
        SSD=XN*(XALFLS(N,N6,N5,N4)+XFEFLS(N,N6,N5,N4)+XHYFLS(N,N6,N5,N4) &
          +XCAFLS(N,N6,N5,N4)+XMGFLS(N,N6,N5,N4)+XNAFLS(N,N6,N5,N4) &
          +XKAFLS(N,N6,N5,N4)+XOHFLS(N,N6,N5,N4)+XSOFLS(N,N6,N5,N4) &
          +XCLFLS(N,N6,N5,N4)+XC3FLS(N,N6,N5,N4)+XH0PFS(N,N6,N5,N4) &
          +XH0BFB(N,N6,N5,N4)+2.0*(XHCFLS(N,N6,N5,N4)+XAL1FS(N,N6,N5,N4) &
          +XALSFS(N,N6,N5,N4)+XFE1FS(N,N6,N5,N4)+XFESFS(N,N6,N5,N4) &
          +XCAOFS(N,N6,N5,N4)+XCACFS(N,N6,N5,N4) &
          +XCASFS(N,N6,N5,N4)+XMGOFS(N,N6,N5,N4)+XMGCFS(N,N6,N5,N4) &
          +XMGSFS(N,N6,N5,N4)+XNACFS(N,N6,N5,N4)+XNASFS(N,N6,N5,N4) &
          +XKASFS(N,N6,N5,N4)+XC0PFS(N,N6,N5,N4)+XC0BFB(N,N6,N5,N4)) &
          +3.0*(XAL2FS(N,N6,N5,N4) &
          +XFE2FS(N,N6,N5,N4)+XCAHFS(N,N6,N5,N4)+XMGHFS(N,N6,N5,N4) &
          +XF1PFS(N,N6,N5,N4)+XC1PFS(N,N6,N5,N4)+XM1PFS(N,N6,N5,N4) &
          +XF1BFB(N,N6,N5,N4)+XC1BFB(N,N6,N5,N4)+XM1BFB(N,N6,N5,N4)) &
          +4.0*(XAL3FS(N,N6,N5,N4)+XFE3FS(N,N6,N5,N4)+XH3PFS(N,N6,N5,N4) &
          +XF2PFS(N,N6,N5,N4)+XC2PFS(N,N6,N5,N4)+XH3BFB(N,N6,N5,N4) &
          +XF2BFB(N,N6,N5,N4)+XC2BFB(N,N6,N5,N4)) &
          +5.0*(XAL4FS(N,N6,N5,N4)+XFE4FS(N,N6,N5,N4)))
        SHD=XN*(XALFHS(N,N6,N5,N4)+XFEFHS(N,N6,N5,N4)+XHYFHS(N,N6,N5,N4) &
          +XCAFHS(N,N6,N5,N4)+XMGFHS(N,N6,N5,N4)+XNAFHS(N,N6,N5,N4) &
          +XKAFHS(N,N6,N5,N4)+XOHFHS(N,N6,N5,N4)+XSOFHS(N,N6,N5,N4) &
          +XCLFHS(N,N6,N5,N4)+XC3FHS(N,N6,N5,N4)+XH0PHS(N,N6,N5,N4) &
          +XH0BHB(N,N6,N5,N4) &
          +2.0*(XHCFHS(N,N6,N5,N4)+XAL1HS(N,N6,N5,N4) &
          +XALSHS(N,N6,N5,N4)+XFE1HS(N,N6,N5,N4)+XFESHS(N,N6,N5,N4) &
          +XCAOHS(N,N6,N5,N4)+XCACHS(N,N6,N5,N4) &
          +XCASHS(N,N6,N5,N4)+XMGOHS(N,N6,N5,N4)+XMGCHS(N,N6,N5,N4) &
          +XMGSHS(N,N6,N5,N4)+XNACHS(N,N6,N5,N4)+XNASHS(N,N6,N5,N4) &
          +XKASHS(N,N6,N5,N4)+XC0PHS(N,N6,N5,N4)+XC0BHB(N,N6,N5,N4)) &
          +3.0*(XAL2HS(N,N6,N5,N4) &
          +XFE2HS(N,N6,N5,N4)+XCAHHS(N,N6,N5,N4)+XMGHHS(N,N6,N5,N4) &
          +XF1PHS(N,N6,N5,N4)+XC1PHS(N,N6,N5,N4)+XM1PHS(N,N6,N5,N4) &
          +XF1BHB(N,N6,N5,N4)+XC1BHB(N,N6,N5,N4)+XM1BHB(N,N6,N5,N4)) &
          +4.0*(XAL3HS(N,N6,N5,N4)+XFE3HS(N,N6,N5,N4)+XH3PHS(N,N6,N5,N4) &
          +XF2PHS(N,N6,N5,N4)+XC2PHS(N,N6,N5,N4)+XH3BHB(N,N6,N5,N4) &
          +XF2BHB(N,N6,N5,N4)+XC2BHB(N,N6,N5,N4)) &
          +5.0*(XAL4HS(N,N6,N5,N4)+XAL4HS(N,N6,N5,N4)))
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
          ECHY=0.337*AZMAX1((XHYFLS(N,N6,N5,N4)+XHYFHS(N,N6,N5,N4))/WX)
          ECOH=0.192*AZMAX1((XOHFLS(N,N6,N5,N4)+XOHFHS(N,N6,N5,N4))/WX)
          ECAL=0.056*AZMAX1((XALFLS(N,N6,N5,N4)+XCAFHS(N,N6,N5,N4))*3.0/WX)
          ECFE=0.051*AZMAX1((XFEFLS(N,N6,N5,N4)+XFEFHS(N,N6,N5,N4))*3.0/WX)
          ECCA=0.060*AZMAX1((XCAFLS(N,N6,N5,N4)+XCAFHS(N,N6,N5,N4))*2.0/WX)
          ECMG=0.053*AZMAX1((XMGFLS(N,N6,N5,N4)+XMGFHS(N,N6,N5,N4))*2.0/WX)
          ECNA=0.050*AZMAX1((XNAFLS(N,N6,N5,N4)+XNAFHS(N,N6,N5,N4))/WX)
          ECKA=0.070*AZMAX1((XKAFLS(N,N6,N5,N4)+XKAFHS(N,N6,N5,N4))/WX)
          ECCO=0.072*AZMAX1((XC3FLS(N,N6,N5,N4)+XC3FHS(N,N6,N5,N4))*2.0/WX)
          ECHC=0.044*AZMAX1((XHCFLS(N,N6,N5,N4)+XHCFHS(N,N6,N5,N4))/WX)
          ECSO=0.080*AZMAX1((XSOFLS(N,N6,N5,N4)+XSOFHS(N,N6,N5,N4))*2.0/WX)
          ECCL=0.076*AZMAX1((XCLFLS(N,N6,N5,N4)+XCLFHS(N,N6,N5,N4))/WX)
          ECNO=0.071*AZMAX1((XNOFLW(N,N6,N5,N4)+XNOFHW(N,N6,N5,N4))/(WX*14.0))
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
      SG=SG+XHGFLS(N,N6,N5,N4)+XHGFLG(N,N6,N5,N4)
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
      IF(ISALTG.NE.0)THEN
        PSS=XN*31.0*(XQSC0P(N,N5,N4)+XQSF1P(N,N5,N4)+XQSC1P(N,N5,N4) &
          +XQSM1P(N,N5,N4)+XQSH3P(N,N5,N4)+XQSF2P(N,N5,N4)+XQSC2P(N,N5,N4))
        TPOU=TPOU-PSS
        SS1=XN*(XQSAL(N,N5,N4)+XQSFE(N,N5,N4) &
          +XQSHY(N,N5,N4)+XQSCA(N,N5,N4)+XQSMG(N,N5,N4)+XQSNA(N,N5,N4) &
          +XQSKA(N,N5,N4)+XQSOH(N,N5,N4)+XQSSO(N,N5,N4)+XQSCL(N,N5,N4) &
          +XQSC3(N,N5,N4)+XQSH0P(N,N5,N4))
        SS2=XN*2.0*(XQSHC(N,N5,N4)+XQSAL1(N,N5,N4)+XQSALS(N,N5,N4) &
          +XQSFE1(N,N5,N4)+XQSFES(N,N5,N4)+XQSCAO(N,N5,N4)+XQSCAC(N,N5,N4) &
          +XQSCAS(N,N5,N4)+XQSMGO(N,N5,N4)+XQSMGC(N,N5,N4)+XQSMGS(N,N5,N4) &
          +XQSNAC(N,N5,N4)+XQSNAS(N,N5,N4)+XQSKAS(N,N5,N4)+XQSC0P(N,N5,N4))
        SS3=XN*3.0*(XQSAL2(N,N5,N4)+XQSFE2(N,N5,N4)+XQSCAH(N,N5,N4) &
          +XQSMGH(N,N5,N4)+XQSF1P(N,N5,N4)+XQSC1P(N,N5,N4)+XQSM1P(N,N5,N4))
        SS4=XN*4.0*(XQSAL3(N,N5,N4)+XQSFE3(N,N5,N4) &
          +XQSH3P(N,N5,N4)+XQSF2P(N,N5,N4)+XQSC2P(N,N5,N4)) &
          +XN*5.0*(XQSAL4(N,N5,N4)+XQSFE4(N,N5,N4))
        SSR=SS1+SS2+SS3+SS4
        TIONOU=TIONOU-SSR
        UIONOU(NY,NX)=UIONOU(NY,NX)-SSR
      ENDIF
    ENDIF
  ENDIF
  end subroutine WaterHeatSoluteBySnowDrift

end module RunoffBalMod
