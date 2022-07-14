module GeochemAPI

  use data_kind_mod, only : r8 => SHR_KIND_R8
  use SoluteChemDataType, only : solute_flx_type, chem_var_type
  use SoluteMod
  use AqueChemDatatype
  use SoilPropertyDataType
  use SoilBGCDataType
  use EcoSIMCtrlDataType
  use SOMDataType
  use EcoSIMSolverPar
  use SoilWaterDataType
  use GridDataType
  implicit none

  public :: soluteModel
  contains


  subroutine soluteModel(I,J,NHW,NHE,NVN,NVS)
!
!     THIS SUBROUTINE CALCULATES ALL SOLUTE TRANSFORMATIONS
!     FROM THERMODYNAMIC EQUILIBRIA
!
  implicit none

  integer, intent(in) :: I, J
  integer, intent(in) :: NHW, NHE, NVN, NVS

  ! declaration of local variables
  integer :: L,NY,NX,NPI
  real(r8) :: RHP1,RHP2,RN3S,RN4S
  type(solute_flx_type) :: solflx
  type(chem_var_type) :: chemvar

!     begin_execution
  NPI=INT(NPH/2)
  DO   NX=NHW,NHE
    DO   NY=NVN,NVS
      DO   L=NU(NY,NX),NL(NY,NX)
        IF(VOLX(L,NY,NX).GT.ZEROS2(NY,NX).AND.VOLWM(NPH,L,NY,NX).GT.ZEROS2(NY,NX))THEN
!
!     WATER VOLUME IN NON-BAND AND BAND SOIL ZONES
!
!     VOLWM=soil water volume
!     VLNH4,VLNHB=fractions of soil volume in NH4 non-band,band
!     VLNO3,VLNOB=fractions of soil volume in N03 non-band,band
!     VLPO4,VLPOB=fractions of soil volume in H2PO4 non-band,band
!     BKVL=soil mass
!
          chemvar%VOLWNH=VOLWM(NPH,L,NY,NX)*VLNH4(L,NY,NX)
          chemvar%VOLWNB=VOLWM(NPH,L,NY,NX)*VLNHB(L,NY,NX)
          chemvar%VOLWNO=VOLWM(NPH,L,NY,NX)*VLNO3(L,NY,NX)
          chemvar%VOLWNZ=VOLWM(NPH,L,NY,NX)*VLNOB(L,NY,NX)
          chemvar%VOLWPO=VOLWM(NPH,L,NY,NX)*VLPO4(L,NY,NX)
          chemvar%VOLWPB=VOLWM(NPH,L,NY,NX)*VLPOB(L,NY,NX)
          IF(BKVL(L,NY,NX).GT.ZEROS(NY,NX))THEN
            chemvar%BKVLX=BKVL(L,NY,NX)
            chemvar%BKVLNH=BKVL(L,NY,NX)*VLNH4(L,NY,NX)
            chemvar%BKVLNB=BKVL(L,NY,NX)*VLNHB(L,NY,NX)
            chemvar%BKVLNO=BKVL(L,NY,NX)*VLNO3(L,NY,NX)
            chemvar%BKVLNZ=BKVL(L,NY,NX)*VLNOB(L,NY,NX)
            chemvar%BKVLPO=BKVL(L,NY,NX)*VLPO4(L,NY,NX)
            chemvar%BKVLPB=BKVL(L,NY,NX)*VLPOB(L,NY,NX)
          ELSE
            chemvar%BKVLX=VOLA(L,NY,NX)
            chemvar%BKVLNH=chemvar%VOLWNH
            chemvar%BKVLNB=chemvar%VOLWNB
            chemvar%BKVLNO=chemvar%VOLWNO
            chemvar%BKVLNZ=chemvar%VOLWNZ
            chemvar%BKVLPO=chemvar%VOLWPO
            chemvar%BKVLPB=chemvar%VOLWPB
          ENDIF

          call UpdateSoilFertlizer(L,NY,NX,chemvar)

          call GeochemAPISend(L,NY,NX,chemvar,solflx)

          call GeoChemEquilibria(chemvar, solflx)

          call GeochemAPIRecv(L,NY,NX,solflx)

          call UpdateFertilizerBand(L,NY,NX)
        ENDIF
      ENDDO
!
!     SURFACE RESIDUE
!
      call UpdateSoluteInSurfaceResidue(NX,NY)
!
    ENDDO
  ENDDO
  end subroutine soluteModel

!------------------------------------------------------------------------

  subroutine GeochemAPISend(L,NY,NX,chemvar,solflx)

  implicit none
  integer, intent(in) :: L,NY,NX
  type(solute_flx_type), intent(inout) :: solflx
  type(chem_var_type), intent(inout) :: chemvar

  chemvar%ZEROS=ZEROS(NY,NX)
  chemvar%XCEC=XCEC(L,NY,NX)
  chemvar%PH=PH(L,NY,NX)
  chemvar%CAL=CAL(L,NY,NX)
  chemvar%CFE=CFE(L,NY,NX)
  chemvar%VOLWM=VOLWM(NPH,L,NY,NX)
  chemvar%ZMG=ZMG(L,NY,NX)
  chemvar%ZNA=ZNA(L,NY,NX)
  chemvar%ZKA=ZKA(L,NY,NX)
  chemvar%CCO2S=CCO2S(L,NY,NX)
  chemvar%CCA=CCA(L,NY,NX)
  chemvar%ZEROS2=ZEROS2(NY,NX)
  chemvar%BKVL=BKVL(L,NY,NX)
  chemvar%XAEC=XAEC(L,NY,NX)
  chemvar%VLNH4=VLNH4(L,NY,NX)
  chemvar%GKC4=GKC4(L,NY,NX)
  chemvar%VLNHB=VLNHB(L,NY,NX)
  chemvar%GKCA=GKCA(L,NY,NX)
  chemvar%GKCH=GKCH(L,NY,NX)
  chemvar%GKCM=GKCM(L,NY,NX)
  chemvar%GKCK=GKCK(L,NY,NX)
  chemvar%GKCN=GKCN(L,NY,NX)
  chemvar%ZNO3S=ZNO3S(L,NY,NX)
  chemvar%ZNO3B=ZNO3B(L,NY,NX)
  chemvar%XZHYS=XZHYS(L,NY,NX)
  chemvar%ZHY=ZHY(L,NY,NX)
  chemvar%ZOH=ZOH(L,NY,NX)
  chemvar%ZAL=ZAL(L,NY,NX)
  chemvar%ZFE=ZFE(L,NY,NX)
  chemvar%ZCA=ZCA(L,NY,NX)
  chemvar%ZSO4=ZSO4(L,NY,NX)
  chemvar%ZCL=ZCL(L,NY,NX)
  chemvar%ZCO3=ZCO3(L,NY,NX)
  chemvar%ZHCO3=ZHCO3(L,NY,NX)
  chemvar%CO2S=CO2S(L,NY,NX)
  chemvar%ZALOH1=ZALOH1(L,NY,NX)
  chemvar%ZALOH2=ZALOH2(L,NY,NX)
  chemvar%ZALOH3=ZALOH3(L,NY,NX)
  chemvar%ZALOH4=ZALOH4(L,NY,NX)
  chemvar%ZALS=ZALS(L,NY,NX)
  chemvar%ZFEOH1=ZFEOH1(L,NY,NX)
  chemvar%ZFEOH2=ZFEOH2(L,NY,NX)
  chemvar%ZFEOH3=ZFEOH3(L,NY,NX)
  chemvar%ZFEOH4=ZFEOH4(L,NY,NX)
  chemvar%ZFES=ZFES(L,NY,NX)
  chemvar%ZCAO=ZCAO(L,NY,NX)
  chemvar%ZCAC=ZCAC(L,NY,NX)
  chemvar%ZCAH=ZCAH(L,NY,NX)
  chemvar%ZCAS=ZCAS(L,NY,NX)
  chemvar%ZMGO=ZMGO(L,NY,NX)
  chemvar%ZMGC=ZMGC(L,NY,NX)
  chemvar%ZMGH=ZMGH(L,NY,NX)
  chemvar%ZMGS=ZMGS(L,NY,NX)
  chemvar%ZNAC=ZNAC(L,NY,NX)
  chemvar%ZNAS=ZNAS(L,NY,NX)
  chemvar%ZKAS=ZKAS(L,NY,NX)
  chemvar%H0PO4=H0PO4(L,NY,NX)
  chemvar%H3PO4=H3PO4(L,NY,NX)
  chemvar%ZFE1P=ZFE1P(L,NY,NX)
  chemvar%ZFE2P=ZFE2P(L,NY,NX)
  chemvar%ZCA0P=ZCA0P(L,NY,NX)
  chemvar%ZCA1P=ZCA1P(L,NY,NX)
  chemvar%ZCA2P=ZCA2P(L,NY,NX)
  chemvar%ZMG1P=ZMG1P(L,NY,NX)
  chemvar%H0POB=H0POB(L,NY,NX)
  chemvar%H3POB=H3POB(L,NY,NX)
  chemvar%ZFE1PB=ZFE1PB(L,NY,NX)
  chemvar%ZFE2PB=ZFE2PB(L,NY,NX)
  chemvar%ZCA0PB=ZCA0PB(L,NY,NX)
  chemvar%ZCA1PB=ZCA1PB(L,NY,NX)
  chemvar%ZCA2PB=ZCA2PB(L,NY,NX)
  chemvar%ZMG1PB=ZMG1PB(L,NY,NX)
  chemvar%XHY=XHY(L,NY,NX)
  chemvar%XAL=XAL(L,NY,NX)
  chemvar%XFE=XFE(L,NY,NX)
  chemvar%XCA=XCA(L,NY,NX)
  chemvar%XMG=XMG(L,NY,NX)
  chemvar%XNA=XNA(L,NY,NX)
  chemvar%XKA=XKA(L,NY,NX)
  chemvar%XHC=XHC(L,NY,NX)
  chemvar%XALO2=XALO2(L,NY,NX)
  chemvar%XFEO2=XFEO2(L,NY,NX)
  chemvar%ORGC=ORGC(L,NY,NX)
  chemvar%PALOH=PALOH(L,NY,NX)
  chemvar%PFEOH=PFEOH(L,NY,NX)
  chemvar%PCACO=PCACO(L,NY,NX)
  chemvar%PCASO=PCASO(L,NY,NX)
  chemvar%VLPOB=VLPOB(L,NY,NX)
  chemvar%VLPO4=VLPO4(L,NY,NX)
  chemvar%VLNOB=VLNOB(L,NY,NX)
  chemvar%VLNO3=VLNO3(L,NY,NX)

  solflx%TRCACO=TRCACO(L,NY,NX)
  solflx%TRNAC=TRNAC(L,NY,NX)
  solflx%TRMGC=TRMGC(L,NY,NX)
  solflx%TRCAC=TRCAC(L,NY,NX)
  solflx%TRMGH=TRMGH(L,NY,NX)
  solflx%TRCAH=TRCAH(L,NY,NX)
  solflx%TRHCO=TRHCO(L,NY,NX)
  solflx%TRXHC=TRXHC(L,NY,NX)
  solflx%TRC2P=TRC2P(L,NY,NX)
  solflx%TRF2P=TRF2P(L,NY,NX)
  solflx%TRC2B=TRC2B(L,NY,NX)
  solflx%TRF2B=TRF2B(L,NY,NX)
  solflx%TRM1P=TRM1P(L,NY,NX)
  solflx%TRC1P=TRC1P(L,NY,NX)
  solflx%TRF1P=TRF1P(L,NY,NX)
  solflx%TRM1B=TRM1B(L,NY,NX)
  solflx%TRC1B=TRC1B(L,NY,NX)
  solflx%TRF1B=TRF1B(L,NY,NX)
  solflx%TRC0B=TRC0B(L,NY,NX)
  solflx%TRH0B=TRH0B(L,NY,NX)
  solflx%TRC0P=TRC0P(L,NY,NX)
  solflx%TRH0P=TRH0P(L,NY,NX)
  solflx%TRFEPO=TRFEPO(L,NY,NX)
  solflx%TRALPO=TRALPO(L,NY,NX)
  solflx%TRFEPB=TRFEPB(L,NY,NX)
  solflx%TRALPB=TRALPB(L,NY,NX)
  solflx%TRCAPM=TRCAPM(L,NY,NX)
  solflx%TRCAPD=TRCAPD(L,NY,NX)
  solflx%TRCPMB=TRCPMB(L,NY,NX)
  solflx%TRCPDB=TRCPDB(L,NY,NX)
  solflx%TRCPHB=TRCPHB(L,NY,NX)
  solflx%TRCAPH=TRCAPH(L,NY,NX)
  solflx%TRFE1=TRFE1(L,NY,NX)
  solflx%TRAL=TRAL(L,NY,NX)
  solflx%TRAL1=TRAL1(L,NY,NX)
  solflx%TRFE2=TRFE2(L,NY,NX)
  solflx%TRAL2=TRAL2(L,NY,NX)
  solflx%TRFE3=TRFE3(L,NY,NX)
  solflx%TRAL3=TRAL3(L,NY,NX)
  solflx%TRFE4=TRFE4(L,NY,NX)
  solflx%TRAL4=TRAL4(L,NY,NX)
  solflx%TRMGO=TRMGO(L,NY,NX)
  solflx%TRCAO=TRCAO(L,NY,NX)
  solflx%TRFEOH=TRFEOH(L,NY,NX)
  solflx%TRALOH=TRALOH(L,NY,NX)
  solflx%TRXFE2=TRXFE2(L,NY,NX)
  solflx%TRALS =TRALS(L,NY,NX)
  solflx%TRB1P =TRB1P(L,NY,NX)
  solflx%TRB2P =TRB2P(L,NY,NX)
  solflx%TRBH0 =TRBH0(L,NY,NX)
  solflx%TRBH1 =TRBH1(L,NY,NX)
  solflx%TRBH2 =TRBH2(L,NY,NX)
  solflx%TRCA  =TRCA(L,NY,NX)
  solflx%TRCAS =TRCAS(L,NY,NX)
  solflx%TRCASO=TRCASO(L,NY,NX)
  solflx%TRCO2 =TRCO2(L,NY,NX)
  solflx%TRCO3 =TRCO3(L,NY,NX)
  solflx%TRFE  =TRFE(L,NY,NX)
  solflx%TRFES =TRFES(L,NY,NX)
  solflx%TRH1B =TRH1B(L,NY,NX)
  solflx%TRH2B =TRH2B(L,NY,NX)
  solflx%TRH1P =TRH1P(L,NY,NX)
  solflx%TRH2P =TRH2P(L,NY,NX)
  solflx%TRH3B =TRH3B(L,NY,NX)
  solflx%TRH3P =TRH3P(L,NY,NX)
  solflx%TRHY  =TRHY(L,NY,NX)
  solflx%TRKA  =TRKA(L,NY,NX)
  solflx%TRKAS =TRKAS(L,NY,NX)
  solflx%TRMG  =TRMG(L,NY,NX)
  solflx%TRMGS =TRMGS(L,NY,NX)
  solflx%TRN3B =TRN3B(L,NY,NX)
  solflx%TRN3S =TRN3S(L,NY,NX)
  solflx%TRN4B =TRN4B(L,NY,NX)
  solflx%TRN4S =TRN4S(L,NY,NX)
  solflx%TRNA  =TRNA(L,NY,NX)
  solflx%TRNAS =TRNAS(L,NY,NX)
  solflx%TROH  =TROH(L,NY,NX)
  solflx%TRSO4 =TRSO4(L,NY,NX)
  solflx%TRX1P =TRX1P(L,NY,NX)
  solflx%TRX2P =TRX2P(L,NY,NX)
  solflx%TRXAL =TRXAL(L,NY,NX)
  solflx%TRXAL2=TRXAL2(L,NY,NX)
  solflx%TRXCA =TRXCA(L,NY,NX)
  solflx%TRXFE =TRXFE(L,NY,NX)
  solflx%TRXH0 =TRXH0(L,NY,NX)
  solflx%TRXH1 =TRXH1(L,NY,NX)
  solflx%TRXH2 =TRXH2(L,NY,NX)
  solflx%TRXHY =TRXHY(L,NY,NX)
  solflx%TRXKA =TRXKA(L,NY,NX)
  solflx%TRXMG =TRXMG(L,NY,NX)
  solflx%TRXN4 =TRXN4(L,NY,NX)
  solflx%TRXNA =TRXNA(L,NY,NX)
  solflx%TRXNB =TRXNB(L,NY,NX)
  solflx%TBCO2 =TBCO2(L,NY,NX)
  solflx%TBION =TBION(L,NY,NX)
  solflx%TRH2O =TRH2O(L,NY,NX)
!  Integers

  end subroutine GeochemAPISend

!-----------------------------------------------------------------------

  subroutine GeochemAPIRecv(L,NY,NX,solflx)

  implicit none
  integer, intent(in) :: L,NY,NX
  type(solute_flx_type), intent(in) :: solflx

  TRN4S(L,NY,NX)=solflx%TRN4S
  TRN4B(L,NY,NX)=solflx%TRN4B
  TRN3S(L,NY,NX)=solflx%TRN3S
  TRN3B(L,NY,NX)=solflx%TRN3B
  TRH1P(L,NY,NX)=solflx%TRH1P
  TRH2P(L,NY,NX)=solflx%TRH2P
  TRH1B(L,NY,NX)=solflx%TRH1B
  TRH2B(L,NY,NX)=solflx%TRH2B
  TRXN4(L,NY,NX)=solflx%TRXN4
  TRXNB(L,NY,NX)=solflx%TRXNB
  TRXH1(L,NY,NX)=solflx%TRXH1
  TRXH2(L,NY,NX)=solflx%TRXH2
  TRX1P(L,NY,NX)=solflx%TRX1P
  TRX2P(L,NY,NX)=solflx%TRX2P
  TRBH1(L,NY,NX)=solflx%TRBH1
  TRBH2(L,NY,NX)=solflx%TRBH2
  TRB1P(L,NY,NX)=solflx%TRB1P
  TRB2P(L,NY,NX)=solflx%TRB2P
  TRALPO(L,NY,NX)=solflx%TRALPO
  TRFEPO(L,NY,NX)=solflx%TRFEPO
  TRCAPD(L,NY,NX)=solflx%TRCAPD
  TRCAPH(L,NY,NX)=solflx%TRCAPH
  TRCAPM(L,NY,NX)=solflx%TRCAPM
  TRALPB(L,NY,NX)=solflx%TRALPB
  TRFEPB(L,NY,NX)=solflx%TRFEPB
  TRCPDB(L,NY,NX)=solflx%TRCPDB
  TRCPHB(L,NY,NX)=solflx%TRCPHB
  TRCPMB(L,NY,NX)=solflx%TRCPMB
  TRAL(L,NY,NX)=solflx%TRAL
  TRFE(L,NY,NX)=solflx%TRFE
  TRHY(L,NY,NX)=solflx%TRHY
  TRCA(L,NY,NX)=solflx%TRCA
  TRMG(L,NY,NX)=solflx%TRMG
  TRNA(L,NY,NX)=solflx%TRNA
  TRKA(L,NY,NX)=solflx%TRKA
  TROH(L,NY,NX)=solflx%TROH
  TRSO4(L,NY,NX)=solflx%TRSO4
  TRCO3(L,NY,NX)=solflx%TRCO3
  TRHCO(L,NY,NX)=solflx%TRHCO
  TRCO2(L,NY,NX)=solflx%TRCO2
  TRAL1(L,NY,NX)=solflx%TRAL1
  TRAL2(L,NY,NX)=solflx%TRAL2
  TRAL3(L,NY,NX)=solflx%TRAL3
  TRAL4(L,NY,NX)=solflx%TRAL4
  TRALS(L,NY,NX)=solflx%TRALS
  TRFE1(L,NY,NX)=solflx%TRFE1
  TRFE2(L,NY,NX)=solflx%TRFE2
  TRFE3(L,NY,NX)=solflx%TRFE3
  TRFE4(L,NY,NX)=solflx%TRFE4
  TRFES(L,NY,NX)=solflx%TRFES
  TRCAO(L,NY,NX)=solflx%TRCAO
  TRCAC(L,NY,NX)=solflx%TRCAC
  TRCAH(L,NY,NX)=solflx%TRCAH
  TRCAS(L,NY,NX)=solflx%TRCAS
  TRMGO(L,NY,NX)=solflx%TRMGO
  TRMGC(L,NY,NX)=solflx%TRMGC
  TRMGH(L,NY,NX)=solflx%TRMGH
  TRMGS(L,NY,NX)=solflx%TRMGS
  TRNAC(L,NY,NX)=solflx%TRNAC
  TRNAS(L,NY,NX)=solflx%TRNAS
  TRKAS(L,NY,NX)=solflx%TRKAS
  TRH0P(L,NY,NX)=solflx%TRH0P
  TRH3P(L,NY,NX)=solflx%TRH3P
  TRF1P(L,NY,NX)=solflx%TRF1P
  TRF2P(L,NY,NX)=solflx%TRF2P
  TRC0P(L,NY,NX)=solflx%TRC0P
  TRC1P(L,NY,NX)=solflx%TRC1P
  TRC2P(L,NY,NX)=solflx%TRC2P
  TRM1P(L,NY,NX)=solflx%TRM1P
  TRH0B(L,NY,NX)=solflx%TRH0B
  TRH3B(L,NY,NX)=solflx%TRH3B
  TRF1B(L,NY,NX)=solflx%TRF1B
  TRF2B(L,NY,NX)=solflx%TRF2B
  TRC0B(L,NY,NX)=solflx%TRC0B
  TRC1B(L,NY,NX)=solflx%TRC1B
  TRC2B(L,NY,NX)=solflx%TRC2B
  TRM1B(L,NY,NX)=solflx%TRM1B
  TRXHY(L,NY,NX)=solflx%TRXHY
  TRXAL(L,NY,NX)=solflx%TRXAL
  TRXFE(L,NY,NX)=solflx%TRXFE
  TRXCA(L,NY,NX)=solflx%TRXCA
  TRXMG(L,NY,NX)=solflx%TRXMG
  TRXNA(L,NY,NX)=solflx%TRXNA
  TRXKA(L,NY,NX)=solflx%TRXKA
  TRXHC(L,NY,NX)=solflx%TRXHC
  TRXAL2(L,NY,NX)=solflx%TRXAL2
  TRXFE2(L,NY,NX)=solflx%TRXFE2
  TRXH0(L,NY,NX)=solflx%TRXH0
  TRBH0(L,NY,NX)=solflx%TRBH0
  TRALOH(L,NY,NX)=solflx%TRALOH
  TRFEOH(L,NY,NX)=solflx%TRFEOH
  TRCACO(L,NY,NX)=solflx%TRCACO
  TRCASO(L,NY,NX)=solflx%TRCASO
  TRH2O(L,NY,NX)=solflx%TRH2O
  TBION(L,NY,NX)=solflx%TBION
  TBCO2(L,NY,NX)=solflx%TBCO2
  end subroutine GeochemAPIRecv

end module GeochemAPI
