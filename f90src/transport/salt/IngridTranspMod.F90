module IngridTranspMod
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use TrnsfrsDataMod
  use EcoSIMCtrlDataType
  use SoilWaterDataType
  use GridDataType
  use SurfLitterDataType
  use SoilBGCDataType
  use SurfSoilDataType
  use SnowDataType
  use ChemTranspDataType
  use EcoSIMSolverPar
  use minimathmod, only : AZMAX1,AZMIN1
  implicit none
  private
  CHARACTER(LEN=*), PARAMETER :: MOD_FILENAME=__FILE__

  real(r8), PARAMETER :: XFRS=0.05_r8
  public :: SaltModelSoluteHydroFlux
  contains


!------------------------------------------------------------------------------------------

  subroutine SaltModelSoluteHydroFlux(I,M,NHW,NHE,NVN,NVS)
!
!     Description:
!
  implicit none
  integer, intent(in) :: I,M,NHW,NHE,NVN,NVS

  integer :: NY,NX,N1,N2
  real(r8) :: FLWRM1
  real(r8) :: trcsa_RFLS1(idsa_beg:idsab_end)
  real(r8) :: trcsa_RFL(idsa_beg:idsab_end)
  real(r8) :: trcsa_RFLS0(idsa_beg:idsa_end)
  real(r8) :: trcsa_DFV(idsa_beg:idsab_end)
!     begin_execution

  DO NX=NHW,NHE
    DO  NY=NVN,NVS
!
!     ADD SINKS FROM SOLUTE.F
      call SoluteSinksInSoil(NY,NX)
!
!     SOLUTE FLUXES FROM MELTING SNOWPACK TO
!     SOIL SURFACE FROM SNOWMELT IN 'WATSUB' AND
!     CONCENTRATIONS IN SNOWPACK
      call SoluteFluxInSnowpack(M,NY,NX,trcsa_RFLS1,trcsa_RFLS0)
!
!     CONVECTIVE SOLUTE EXCHANGE BETWEEN RESIDUE AND SOIL SURFACE
!
      FLWRM1=FLWRM(M,NY,NX)
!
!     FLWRM=litter-soil water flux from watsub.f

      IF(FLWRM1.GT.0.0_r8)THEN
        call Residue2TopsoilSoluteAdvExch(M,NY,NX,FLWRM1,trcsa_RFL)
!
      ELSE
        call Topsoil2ResidueSoluteAdvExch(M,NY,NX,FLWRM1,trcsa_RFL)
      ENDIF
!
!     DIFFUSIVE FLUXES OF GASES AND SOLUTES BETWEEN RESIDUE AND
!     SOIL SURFACE FROM AQUEOUS DIFFUSIVITIES
!     AND CONCENTRATION DIFFERENCES
!
!     VOLT,DLYR,AREA=soil surface volume, thickness, area
!     VOLWM=micropore water-filled porosity from watsub.f
!
      IF((VOLT(0,NY,NX).GT.ZEROS(NY,NX).AND.VOLWM(M,0,NY,NX).GT.ZEROS2(NY,NX)) &
        .AND.(VOLWM(M,NU(NY,NX),NY,NX).GT.ZEROS2(NY,NX)))THEN

        call TopsoilResidueSolutedifusExch(M,NY,NX,FLWRM1,trcsa_DFV)
      ELSE
        trcsa_DFV(idsa_beg:idsab_end)=0.0_r8
      ENDIF
!
!     TOTAL MICROPORE AND MACROPORE SOLUTE TRANSPORT FLUXES BETWEEN
!     ADJACENT GRID CELLS = CONVECTIVE + DIFFUSIVE FLUXES

      call TopsoilResidueFluxAdvPlusDifus(NY,NX,trcsa_RFLS1,&
        trcsa_RFL,trcsa_RFLS0,trcsa_DFV)
!
      call AccumHourlyTopsoilReisdueFlux(NY,NX,trcsa_RFL,trcsa_DFV,&
        trcsa_RFLS0,trcsa_RFLS1)
!
!     MACROPORE-MICROPORE SOLUTE EXCHANGE IN SOIL
!     SURFACE LAYER FROM WATER EXCHANGE IN 'WATSUB' AND
!     FROM MACROPORE OR MICROPORE SOLUTE CONCENTRATIONS
!
!     FINHM=macro-micropore water transfer from watsub.f
!     VOLWM,VOLWHM=micropore,macropore water volume
!     RFL*=convective macropore-micropore solute transfer
!     VLNH4,VLNO3,VLPO4=non-band NH4,NO3,PO4 volume fraction
!     VLNHB,VLNOB,VLPOB=band NH4,NO3,PO4 volume fraction
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
!     *H2,*2=macropore,micropore solute content
!
!     MACROPORE TO MICROPORE TRANSFER
!
      IF(FINHM(M,NU(NY,NX),NY,NX).GT.0.0)THEN
        call MacToMicPoreSoluteAdvExchange(M,NY,NX)
!
!     MICROPORE TO MACROPORE TRANSFER
!
      ELSEIF(FINHM(M,NU(NY,NX),NY,NX).LT.0.0)THEN
        call MicToMacPoreSoluteAdvExchange(M,NY,NX)
!
!     NO MACROPORE TO MICROPORE TRANSFER
!
      ELSE
        trcsa_RFL(idsa_beg:idsab_end)=0.0_r8
      ENDIF
!
!     DIFFUSIVE FLUXES OF SOLUTES BETWEEN MICROPORES AND
!     MACROPORES FROM AQUEOUS DIFFUSIVITIES AND CONCENTRATION DIFFERENCES
!
      IF(VOLWHM(M,NU(NY,NX),NY,NX).GT.ZEROS2(NY,NX))THEN
        call MacMicPoreSoluteDifusExchange(M,NY,NX,trcsa_DFV)
      ELSE
        trcsa_DFV(idsa_beg:idsab_end)=0.0_r8
      ENDIF
!
!     TOTAL CONVECTIVE +DIFFUSIVE TRANSFER BETWEEN MACROPOES AND MICROPORES
!      call MacMicPoreFluxAdvPlusDifus(NY,NX,trcsa_DFV,trcsa_RFL)
!
!     ACCUMULATE HOURLY FLUXES FOR USE IN REDIST.F
      call AccumHourlyMicMacPoreFlux(NY,NX)
!
!     SOLUTE TRANSPORT FROM WATER OVERLAND FLOW
!     IN 'WATSUB' AND FROM SOLUTE CONCENTRATIONS
!     IN SOIL SURFACE LAYER
!
!     QRM=runoff from watsub.f
!     RQR*0=solute in runoff
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
!     VOLWM=litter water volume from watsub.f
!     *S2=litter solute content
!     N2,N1=NY,NX of source grid cell
!     N5,N4=NY,NX of destination grid cell
!     X*QRS=accumulated hourly solute in runoff
!
!
      N1=NX
      N2=NY
      IF(QRM(M,N2,N1).GT.ZEROS(N2,N1))THEN
        call SoluteFluxBySurfaceOutflow(M,N1,N2)
      ELSE
        trcsa_RQR0(idsa_beg:idsa_end,N2,N1)=0.0_r8
      ENDIF
!
      call UpdateSoluteInSurfNeighbors(M,N1,N2,NY,NX,NHW,NHE,NVN,NVS)
!
!     SOLUTE FLUXES BETWEEN ADJACENT GRID CELLS
      call UpdateSoluteInSubsurfNeighbors(M,NY,NX,NHW,NHE,NVN,NVS)

    ENDDO
  ENDDO
  end subroutine SaltModelSoluteHydroFlux

!------------------------------------------------------------------------------------------

  subroutine SoluteSinksInSoil(NY,NX)
!
!     Description:
!
  implicit none
  integer, intent(in) :: NY,NX

  integer :: L,NTSA
!     begin_execution
!
!     RZ*2=solute flux at time step for flux calculations
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
!
  DO L=NU(NY,NX),NL(NY,NX)
    DO NTSA=idsa_beg,idsab_end
      trcsa_solml2(NTSA,L,NY,NX)=trcsa_solml2(NTSA,L,NY,NX)-trcsa_solml2R(NTSA,L,NY,NX)
    ENDDO
  ENDDO
  end subroutine SoluteSinksInSoil
!------------------------------------------------------------------------------------------

  subroutine SoluteFluxInSnowpack(M,NY,NX,trcsa_RFLS1,trcsa_RFLS0)
!
!     Description:
!
  implicit none
  integer, intent(in) :: M,NY,NX
  integer :: L,ICHKL,L2,NTSA
  real(r8) :: VFLWW,VFLWR,VFLWS,VFLWPO4,VFLWPOB
  real(r8), intent(out) :: trcsa_RFLS1(idsa_beg:idsab_end)
  real(r8), intent(out) :: trcsa_RFLS0(idsa_beg:idsa_end)
!     begin_execution
!
!     VHCPWM,VHCPWX=current,minimum volumetric heat capacity of snowpack
!     VOLWSL=snowpack water content
!     FLQWM=snowpack water flux
!     R*BLS=solute flux in snowpack
!     X*BLS=hourly solute flux in snowpack
!
  ICHKL=0
  DO L=1,JS
    IF(VHCPWM(M,L,NY,NX).GT.VHCPWX(NY,NX))THEN
      L2=MIN(JS,L+1)
      IF(L.LT.JS.AND.VHCPWM(M,L2,NY,NX).GT.VHCPWX(NY,NX))THEN
        IF(VOLWSL(L,NY,NX).GT.ZEROS2(NY,NX))THEN
          VFLWW=AZMAX1(AMIN1(1.0_r8,FLQWM(M,L2,NY,NX)/VOLWSL(L,NY,NX)))
        ELSE
          VFLWW=1.0_r8
        ENDIF

        DO NTSA=idsa_beg,idsa_end
          trcsa_RBLS(NTSA,L2,NY,NX)=trcsa_sosml2(NTSA,L,NY,NX)*VFLWW
          trcsa_XBLS(NTSA,L2,NY,NX)=trcsa_XBLS(NTSA,L2,NY,NX)+trcsa_RBLS(NTSA,L2,NY,NX)
        ENDDO

      ELSE
        IF(L.LT.JS)THEN
          trcsa_RBLS(idsa_beg:idsa_end,L2,NY,NX)=0.0_r8
        ENDIF
!
!     SNOWPACK SOLUTE DISCHARGE TO SURFACE LITTER, SOIL SURFACE
!
!     VOLWSL=snowpack water content
!     FLQRM,FLQSM,FLQHM=total water flux to litter,soil micropore,macropore
!     CVRD,BARE=litter cover fraction,1-CVRD
!     VLNH4,VLNO3,VLPO4=non-band NH4,NO3,PO4 volume fraction
!     VLNHB,VLNOB,VLPOB=band NH4,NO3,PO4 volume fraction
!     R*S0=solute flux to surface litter
!     R*S1,R*B1=solute flux to soil surface non-band,band
!
        IF(ICHKL.EQ.0)THEN
          IF(VOLWSL(L,NY,NX).GT.ZEROS2(NY,NX))THEN
            VFLWR=AZMAX1(AMIN1(1.0_r8,FLQRM(M,NY,NX)/VOLWSL(L,NY,NX)))
            VFLWS=AZMAX1(AMIN1(1.0_r8,(FLQSM(M,NY,NX)+FLQHM(M,NY,NX))/VOLWSL(L,NY,NX)))
          ELSE
            VFLWR=CVRD(NY,NX)
            VFLWS=BARE(NY,NX)
          ENDIF
          VFLWPO4=VFLWS*trcs_VLN(ids_H1PO4,NU(NY,NX),NY,NX)
          VFLWPOB=VFLWS*trcs_VLN(ids_H1PO4B,NU(NY,NX),NY,NX)

          DO NTSA=idsa_beg,idsa_end
            trcsa_RFLS0(NTSA)=trcsa_sosml2(NTSA,L,NY,NX)*VFLWR
          ENDDO

          DO NTSA=idsa_beg,idsa_psoil_beg-1
            trcsa_RFLS1(NTSA)=trcsa_sosml2(NTSA,L,NY,NX)*VFLWS
          ENDDO

          DO NTSA=idsa_psoil_beg,idsa_psoil_end
            trcsa_RFLS1(NTSA)=trcsa_sosml2(NTSA,L,NY,NX)*VFLWPO4
          ENDDO

          DO NTSA=idsa_pband_beg,idsa_pband_end
            trcsa_RFLS1(NTSA)=trcsa_sosml2(NTSA,L,NY,NX)*VFLWPOB
          ENDDO
          ICHKL=1
        ENDIF
      ENDIF
    ELSE
      trcsa_RFLS0(idsa_beg:idsa_end)=0.0_r8
      trcsa_RFLS1(idsa_beg:idsab_end)=0.0_r8
    ENDIF
  ENDDO
  end subroutine SoluteFluxInSnowpack
!------------------------------------------------------------------------------------------

  subroutine Residue2TopsoilSoluteAdvExch(M,NY,NX,FLWRM1,trcsa_RFL)
!
!     Description:
!
  implicit none
  integer, intent(in) :: M,NY,NX
  real(r8),intent(in) :: FLWRM1
  real(r8),intent(out) :: trcsa_RFL(idsa_beg:idsab_end)
  real(r8) :: VFLW
  integer :: NTSA

!     begin_execution
!
!     IF WATER FLUX FROM 'WATSUB' IS FROM RESIDUE TO
!     SOIL SURFACE THEN CONVECTIVE TRANSPORT IS THE PRODUCT
!     OF WATER FLUX AND MICROPORE GAS OR SOLUTE CONCENTRATIONS
!     IN RESIDUE
!
!     VOLWM=litter water volume
!     RFL*=soil-litter convective solute flux
!     Z*2=litter solute content
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
!     VLNH4,VLNO3,VLPO4=non-band NH4,NO3,PO4 volume fraction
!     VLNHB,VLNOB,VLPOB=band NH4,NO3,PO4 volume fraction
!
  IF(VOLWM(M,0,NY,NX).GT.ZEROS2(NY,NX))THEN
    VFLW=AZMAX1(AMIN1(VFLWX,FLWRM1/VOLWM(M,0,NY,NX)))
  ELSE
    VFLW=VFLWX
  ENDIF
  DO NTSA=idsa_beg,idsa_psoil_beg-1
    trcsa_RFL(NTSA)=VFLW*AZMAX1(trcsa_solml2(NTSA,0,NY,NX))
  ENDDO

  DO NTSA=idsa_psoil_beg,idsa_psoil_end
    trcsa_RFL(NTSA)=VFLW*AZMAX1(trcsa_solml2(NTSA,0,NY,NX)) &
      *trcs_VLN(ids_H1PO4,NU(NY,NX),NY,NX)
  ENDDO

  DO NTSA=idsa_pband_beg,idsa_pband_end
    trcsa_RFL(NTSA)=VFLW*AZMAX1(trcsa_solml2(NTSA,0,NY,NX)) &
      *trcs_VLN(ids_H1PO4B,NU(NY,NX),NY,NX)
  ENDDO
  end subroutine Residue2TopsoilSoluteAdvExch
!------------------------------------------------------------------------------------------

  subroutine Topsoil2ResidueSoluteAdvExch(M,NY,NX,FLWRM1,trcsa_RFL)
!
!     Description:
!
  implicit none
  integer, intent(in) :: M,NY,NX
  real(r8),intent(in) :: FLWRM1
  real(r8),intent(out) :: trcsa_RFL(idsa_beg:idsab_end)
  integer :: NTSA
  real(r8) :: VFLW

!     begin_execution
!
!     IF WATER FLUX FROM 'WATSUB' IS TO RESIDUE FROM
!     SOIL SURFACE THEN CONVECTIVE TRANSPORT IS THE PRODUCT
!     OF WATER FLUX AND MICROPORE GAS OR SOLUTE CONCENTRATIONS
!     IN SOIL SURFACE
!
!     VOLWM=litter water volume
!     RFL*=soil-litter convective solute flux
!     Z*2=soil solute content
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
!     VLNH4,VLNO3,VLPO4=non-band NH4,NO3,PO4 volume fraction
!     VLNHB,VLNOB,VLPOB=band NH4,NO3,PO4 volume fraction
  IF(VOLWM(M,NU(NY,NX),NY,NX).GT.ZEROS2(NY,NX))THEN
    VFLW=AZMIN1(AMAX1(-VFLWX,FLWRM1/VOLWM(M,NU(NY,NX),NY,NX)))
  ELSE
    VFLW=-VFLWX
  ENDIF

  DO NTSA=idsa_beg,idsab_end
    trcsa_RFL(NTSA)=VFLW*AZMAX1(trcsa_solml2(NTSA,NU(NY,NX),NY,NX))
  ENDDO
  end subroutine Topsoil2ResidueSoluteAdvExch
!------------------------------------------------------------------------------------------

  subroutine TopsoilResidueSolutedifusExch(M,NY,NX,FLWRM1,trcsa_DFV)
!
!     Description:
!
  implicit none
  integer, intent(in) :: M,NY,NX
  real(r8),intent(in) :: FLWRM1
  real(r8),intent(out) :: trcsa_DFV(idsa_beg:idsab_end)
  real(r8) :: DLYR1,DLYR2,TORTL,DISPN,DIFPO,DIFAL,DIFFE,DIFCA,DIFMG,DIFNA,DIFKA
  real(r8) :: trcsa_DIFC(idsa_beg:idsa_psoil_beg-1)
  real(r8) :: trcsa_solCl1(idsa_beg:idsab_end)
  real(r8) :: trcsa_solCl2(idsa_beg:idsab_end)
  real(r8) :: VOLWPB,VOLWPA,TORT0,TORT1
  real(r8) :: DLYR0
  integer :: NTSA
!     begin_execution
!
!     MICROPORE CONCENTRATIONS FROM WATER IN RESIDUE AND SOIL SURFACE
!
!     C*1,C*2=solute concentration in litter, soil
!     Z*1,Z*2=solute content in litter, soil
!
  VOLWPA=VOLWM(M,NU(NY,NX),NY,NX)*trcs_VLN(ids_H1PO4,NU(NY,NX),NY,NX)
  VOLWPB=VOLWM(M,NU(NY,NX),NY,NX)*trcs_VLN(ids_H1PO4B,NU(NY,NX),NY,NX)

  DO NTSA=idsa_beg,idsa_end
    trcsa_solCl1(NTSA)=AZMAX1(trcsa_solml2(NTSA,0,NY,NX)/VOLWM(M,0,NY,NX))
  ENDDO

  DO NTSA=idsa_beg,idsa_psoil_beg-1
    trcsa_solCl2(NTSA)=AZMAX1(trcsa_solml2(NTSA,NU(NY,NX),NY,NX)/VOLWM(M,NU(NY,NX),NY,NX))
  ENDDO

  IF(VOLWPA.GT.ZEROS2(NY,NX))THEN
    DO NTSA=idsa_psoil_beg,idsa_psoil_end
      trcsa_solCl2(NTSA)=AZMAX1(trcsa_solml2(NTSA,NU(NY,NX),NY,NX)/VOLWPA)
    ENDDO
  ELSE
    trcsa_solCl2(idsa_psoiL_beg:idsa_psoil_end)=0.0_r8
  ENDIF
  IF(VOLWPB.GT.ZEROS2(NY,NX))THEN
    DO NTSA=idsa_pband_beg,idsa_pband_end
      trcsa_solCl2(NTSA)=AZMAX1(trcsa_solml2(NTSA,NU(NY,NX),NY,NX)/VOLWPB)
    ENDDO
  ELSE
    DO NTSA=0,idsa_nuts
      trcsa_solCl2(idsa_H0PO4B+NTSA)=trcsa_solCl2(idsa_H0PO4+NTSA)
    ENDDO
  ENDIF
!
!     DIFFUSIVITIES IN RESIDUE AND SOIL SURFACE
!
!     DLYR0,DLYR1=litter, soil surface thickness
!     TORT=tortuosity from hour1.f
!     CVRD,BARE=litter cover fraction,1-CVRD
!     DISP=dispersivity parameter
!     FLWRM=litter-soil water flux from watsub.f
!     DIF*=aqueous diffusivity-dispersivity in soil surface
!     *SGL*=solute diffusivity from hour1.f
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
!     DIF*=aqueous diffusivity-dispersivity between litter and soil surface
!
  DLYR0=AMAX1(ZERO2,DLYR(3,0,NY,NX))
  TORT0=TORT(M,0,NY,NX)*CVRD(NY,NX)
  DLYR1=AMAX1(ZERO2,DLYR(3,NU(NY,NX),NY,NX))
  TORT1=TORT(M,NU(NY,NX),NY,NX)
  TORTL=AMIN1(1.0_r8,(TORT0+TORT1)/(DLYR0+DLYR1))
  DISPN=DISP(3,NU(NY,NX),NY,NX)*AMIN1(VFLWX,ABS(FLWRM1/AREA(3,NU(NY,NX),NY,NX)))

  DIFPO=(POSGL2(NU(NY,NX),NY,NX)*TORTL+DISPN)*AREA(3,NU(NY,NX),NY,NX)

  DIFAL=(ALSGL2(NU(NY,NX),NY,NX)*TORTL+DISPN)*AREA(3,NU(NY,NX),NY,NX)
  trcsa_DIFC((/idsa_Al,idsa_AlOH,idsa_AlOH2,idsa_AlOH3,idsa_AlOH4,idsa_AlSO4/))=DIFAL

  DIFFE=(FESGL2(NU(NY,NX),NY,NX)*TORTL+DISPN)*AREA(3,NU(NY,NX),NY,NX)
  trcsa_DIFC((/idsa_Fe,idsa_FeOH,idsa_FeOH2,idsa_FeOH3,idsa_FeOH4,idsa_FeSO4/))=DIFFE

  trcsa_DIFC(idsa_Hp)=(HYSGL2(NU(NY,NX),NY,NX)*TORTL+DISPN)*AREA(3,NU(NY,NX),NY,NX)

  DIFCA=(CASGL2(NU(NY,NX),NY,NX)*TORTL+DISPN)*AREA(3,NU(NY,NX),NY,NX)
  trcsa_DIFC((/idsa_Ca,idsa_CaOH2,idsa_CO3,idsa_HCO3,idsa_CaSO4/))=DIFCA

  DIFMG=(GMSGL2(NU(NY,NX),NY,NX)*TORTL+DISPN)*AREA(3,NU(NY,NX),NY,NX)
  trcsa_DIFC((/idsa_Mg,idsa_MgOH2,idsa_MgCO3,idsa_MgHCO3,idsa_SO4/))=DIFMG

  DIFNA=(ANSGL2(NU(NY,NX),NY,NX)*TORTL+DISPN)*AREA(3,NU(NY,NX),NY,NX)
  trcsa_DIFC((/idsa_Na,idsa_NaCO3,idsa_NaSO4/))=DIFNA

  DIFKA=(AKSGL2(NU(NY,NX),NY,NX)*TORTL+DISPN)*AREA(3,NU(NY,NX),NY,NX)
  trcsa_DIFC((/idsa_K,idsa_KSO4/))=DIFKA

  trcsa_DIFC(idsa_OH)=(OHSGL2(NU(NY,NX),NY,NX)*TORTL+DISPN)*AREA(3,NU(NY,NX),NY,NX)
  trcsa_DIFC(idsa_SO4)=(SOSGL2(NU(NY,NX),NY,NX)*TORTL+DISPN)*AREA(3,NU(NY,NX),NY,NX)
  trcsa_DIFC(idsa_Cl)=(CLSXL2(NU(NY,NX),NY,NX)*TORTL+DISPN)*AREA(3,NU(NY,NX),NY,NX)
  trcsa_DIFC(idsa_CO3)=(C3SGL2(NU(NY,NX),NY,NX)*TORTL+DISPN)*AREA(3,NU(NY,NX),NY,NX)
  trcsa_DIFC(idsa_HCO3)=(HCSGL2(NU(NY,NX),NY,NX)*TORTL+DISPN)*AREA(3,NU(NY,NX),NY,NX)
!
!     DIFFUSIVE FLUXES BETWEEN CURRENT AND ADJACENT GRID CELL
!     MICROPORES
!
!     DFV*S,DFV*B=diffusive solute flux between litter and soil surface in non-band,band
!     DIF*=aqueous diffusivity-dispersivity between litter and soil surface
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
!     C*1,C*2=solute concentration in litter,soil surface
!

  DO NTSA=idsa_beg,idsa_psoil_beg-1
    trcsa_DFV(NTSA)=trcsa_DIFC(NTSA)*(trcsa_solCl1(NTSA)-trcsa_solCl2(NTSA))
  ENDDO

  DO NTSA=idsa_psoil_beg,idsa_psoil_end
    trcsa_DFV(NTSA)=DIFPO*(trcsa_solCl1(NTSA)-trcsa_solCl2(NTSA)) &
      *trcs_VLN(ids_H1PO4,NU(NY,NX),NY,NX)
  ENDDO

  DO NTSA=0,idsa_nuts
    trcsa_DFV(idsa_H0PO4B+NTSA)=DIFPO*(trcsa_solCl1(idsa_H0PO4+NTSA)&
      -trcsa_solCl2(idsa_H0PO4B+NTSA))*trcs_VLN(ids_H1PO4B,NU(NY,NX),NY,NX)
  ENDDO
  end subroutine TopsoilResidueSolutedifusExch
!------------------------------------------------------------------------------------------

  subroutine TopsoilResidueFluxAdvPlusDifus(NY,NX,trcsa_RFLS1,&
    trcsa_RFL,trcsa_RFLS0,trcsa_DFV)
!
!     Description:
!
  implicit none
  integer, intent(in) :: NY,NX
  real(r8), intent(in) :: trcsa_RFLS1(idsa_beg:idsab_end)
  real(r8), intent(in) :: trcsa_RFL(idsa_beg:idsab_end)
  real(r8), intent(in) :: trcsa_RFLS0(idsa_beg:idsa_end)
  real(r8), intent(in) :: trcsa_DFV(idsa_beg:idsab_end)
  integer :: NTSA
!     begin_execution
!
!     R*FLS=convective + diffusive solute flux between litter, soil surface
!     R*FLW,R*FLB=convective + diffusive solute flux into soil in non-band,band
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
!     R*FL0,R*FL1=convective flux into surface litter, soil surface
!     RFL*=convective flux between surface litter and soil surface
!     DFV*=diffusive solute flux between litter and soil surface
!
  DO NTSA=idsa_beg,idsa_end
    trcsa_RFLS(NTSA,3,0,NY,NX)=trcsa_RFL0(NTSA,NY,NX)+trcsa_RFLS0(NTSA) &
      -trcsa_RFL(NTSA)-trcsa_DFV(NTSA)
  ENDDO

  DO NTSA=idsa_beg,idsab_end
    trcsa_RFLS(NTSA,3,NU(NY,NX),NY,NX)=trcsa_RFL1(NTSA,NY,NX) &
      +trcsa_RFLS1(NTSA)+trcsa_RFL(NTSA)+trcsa_DFV(NTSA)
  ENDDO

  end subroutine TopsoilResidueFluxAdvPlusDifus
!------------------------------------------------------------------------------------------

  subroutine AccumHourlyTopsoilReisdueFlux(NY,NX,trcsa_RFL,&
    trcsa_DFV,trcsa_RFLS0,trcsa_RFLS1)
!
!     Description:
!
  implicit none
  integer, intent(in) :: NY,NX
  REAL(R8), intent(in) :: trcsa_RFL(idsa_beg:idsab_end)
  real(r8), intent(in) :: trcsa_DFV(idsa_beg:idsab_end)
  real(r8), intent(in) :: trcsa_RFLS0(idsa_beg:idsa_end)
  real(r8), intent(in) :: trcsa_RFLS1(idsa_beg:idsab_end)
  integer :: NTSA
!     begin_execution
!     ACCUMULATE HOURLY FLUXES FOR USE IN REDIST.F
!
!     X*FLS=hourly convective + diffusive solute flux
!     X*FLW,X*FLB= hourly convective + diffusive solute flux in non-band,band
!
  DO NTSA=idsa_beg,idsa_end
    trcsa_XFLS(NTSA,3,0,NY,NX)=trcsa_XFLS(NTSA,3,0,NY,NX)+trcsa_RFLS0(NTSA) &
      -trcsa_RFL(NTSA)-trcsa_DFV(NTSA)
  ENDDO

  DO NTSA=idsa_beg,idsab_end
    trcsa_XFLS(NTSA,3,NU(NY,NX),NY,NX)=trcsa_XFLS(NTSA,3,NU(NY,NX),NY,NX) &
      +trcsa_RFLS1(NTSA)+trcsa_RFL(NTSA)+trcsa_DFV(NTSA)
  ENDDO
  end subroutine AccumHourlyTopsoilReisdueFlux
!------------------------------------------------------------------------------------------

  subroutine MacToMicPoreSoluteAdvExchange(M,NY,NX)
!
!     Description:
!
  implicit none
  integer, intent(in) :: M,NY,NX
  real(r8) :: trcsa_RFL(idsa_beg:idsab_end)
  real(r8) :: VFLW
  integer :: NTSA
!     begin_execution

  IF(VOLWHM(M,NU(NY,NX),NY,NX).GT.ZEROS2(NY,NX))THEN
    VFLW=AZMAX1(AMIN1(VFLWX,FINHM(M,NU(NY,NX),NY,NX)/VOLWHM(M,NU(NY,NX),NY,NX)))
  ELSE
    VFLW=VFLWX
  ENDIF

  DO NTSA=idsa_beg,idsa_psoil_beg-1
    trcsa_RFL(NTSA)=VFLW*AZMAX1(trcsa_soHml2(NTSA,NU(NY,NX),NY,NX))
  ENDDO

  DO NTSA=idsa_psoil_beg,idsa_psoil_end
    trcsa_RFL(NTSA)=VFLW*AZMAX1(trcsa_soHml2(NTSA,NU(NY,NX),NY,NX)) &
      *trcs_VLN(ids_H1PO4,NU(NY,NX),NY,NX)
  ENDDO

  DO NTSA=idsa_pband_beg,idsa_pband_end
    trcsa_RFL(NTSA)=VFLW*AZMAX1(trcsa_soHml2(NTSA,NU(NY,NX),NY,NX)) &
      *trcs_VLN(ids_H1PO4B,NU(NY,NX),NY,NX)
  ENDDO

  DO NTSA=idsa_beg,idsab_end
    trcsa_RFXS(NTSA,NU(NY,NX),NY,NX)=trcsa_RFL(NTSA)
  ENDDO
  end subroutine MacToMicPoreSoluteAdvExchange
!------------------------------------------------------------------------------------------

  subroutine MicToMacPoreSoluteAdvExchange(M,NY,NX)
!
!     Description:
!
  implicit none
  integer, intent(in) :: M,NY,NX
  real(r8) :: trcsa_RFL(idsa_beg:idsab_end)
  integer :: NTSA
  real(r8) :: VFLW
!     begin_execution
  IF(VOLWM(M,NU(NY,NX),NY,NX).GT.ZEROS2(NY,NX))THEN
    VFLW=AZMIN1(AMAX1(-VFLWX,FINHM(M,NU(NY,NX),NY,NX)/VOLWM(M,NU(NY,NX),NY,NX)))
  ELSE
    VFLW=-VFLWX
  ENDIF

  DO NTSA=idsa_beg,idsa_psoil_beg-1
    trcsa_RFL(idsa_Al)=VFLW*AZMAX1(trcsa_solml2(idsa_Al,NU(NY,NX),NY,NX))
  ENDDO

  DO NTSA=idsa_psoil_beg,idsa_psoil_end
    trcsa_RFL(NTSA)=VFLW*AZMAX1(trcsa_solml2(NTSA,NU(NY,NX),NY,NX)) &
      *trcs_VLN(ids_H1PO4,NU(NY,NX),NY,NX)
  ENDDO

  DO NTSA=idsa_pband_beg,idsa_pband_end
    trcsa_RFL(NTSA)=VFLW*AZMAX1(trcsa_solml2(NTSA,NU(NY,NX),NY,NX)) &
      *trcs_VLN(ids_H1PO4B,NU(NY,NX),NY,NX)
  ENDDO

  DO NTSA=idsa_beg,idsab_end
    trcsa_RFXS(NTSA,NU(NY,NX),NY,NX)=trcsa_RFL(NTSA)
  ENDDO
  end subroutine MicToMacPoreSoluteAdvExchange

!------------------------------------------------------------------------------------------

  subroutine MacMicPoreSoluteDifusExchange(M,NY,NX,trcsa_DFV)
!
!     Description:
!
  implicit none
  integer, intent(in) :: M,NY,NX
  real(r8), intent(out) :: trcsa_DFV(idsa_beg:idsab_end)
  real(r8) :: VOLWHS,VOLWT,VOLWMNU
  integer :: NTSA
!     begin_execution
!
!     VOLWM,VOLWHM=micropore,macropore water volume
!     XFRS*VOLT=maximum macropore volume for solute transfer
!     DFV*=diffusive macropore-micropore solute transfer
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
!     VLNH4,VLNO3,VLPO4=non-band NH4,NO3,PO4 volume fraction
!     VLNHB,VLNOB,VLPOB=band NH4,NO3,PO4 volume fraction
!     XNPH=1/no. of cycles h-1 for water, heat and solute flux calculations
!     *H2,*2=macropore,micropore solute content
  VOLWMNU=VOLWM(M,NU(NY,NX),NY,NX)
  VOLWHS=AMIN1(XFRS*VOLT(NU(NY,NX),NY,NX),VOLWHM(M,NU(NY,NX),NY,NX))
  VOLWT=VOLWM(M,NU(NY,NX),NY,NX)+VOLWHS

  DO NTSA=idsa_beg,idsa_psoil_beg-1
    trcsa_DFV(NTSA)=XNPH*(AZMAX1(trcsa_soHml2(NTSA,NU(NY,NX),NY,NX))*VOLWMNU &
      -AZMAX1(trcsa_solml2(NTSA,NU(NY,NX),NY,NX))*VOLWHS)/VOLWT
  ENDDO

  DO NTSA=idsa_psoil_beg,idsa_psoil_end
    trcsa_DFV(NTSA)=XNPH*(AZMAX1(trcsa_soHml2(NTSA,NU(NY,NX),NY,NX))*VOLWMNU &
      -AZMAX1(trcsa_solml2(NTSA,NU(NY,NX),NY,NX))*VOLWHS)/VOLWT*trcs_VLN(ids_H1PO4,NU(NY,NX),NY,NX)
  ENDDO

  DO NTSA=idsa_pband_beg,idsa_pband_end
    trcsa_DFV(NTSA)=XNPH*(AZMAX1(trcsa_soHml2(NTSA,NU(NY,NX),NY,NX))*VOLWMNU &
      -AZMAX1(trcsa_solml2(NTSA,NU(NY,NX),NY,NX))*VOLWHS)/VOLWT*trcs_VLN(ids_H1PO4B,NU(NY,NX),NY,NX)
  ENDDO


!
!     R*FXS=convective + diffusive solute flux between macropores and micropores
!     RFL*=convective flux between macropores and micropores
!     DFV*=diffusive solute flux between macropores and micropores
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
!
  DO NTSA=idsa_beg,idsab_end
    trcsa_RFXS(NTSA,NU(NY,NX),NY,NX)=trcsa_RFXS(NTSA,NU(NY,NX),NY,NX)+trcsa_DFV(NTSA)
  ENDDO

  end subroutine MacMicPoreSoluteDifusExchange

!------------------------------------------------------------------------------------------

  subroutine AccumHourlyMicMacPoreFlux(NY,NX)
!
!     Description:
!
  implicit none
  integer, intent(in) :: NY,NX
  INTEGER :: NTSA
!     begin_execution
!
!     X*FXS=hourly convective + diffusive solute flux between macropores and micropores
!     R*FXS=total convective + diffusive solute flux between macropores and micropores
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
!
  DO NTSA=idsa_beg,idsab_end
    trcsa_XFXS(NTSA,NU(NY,NX),NY,NX)=trcsa_XFXS(NTSA,NU(NY,NX),NY,NX)+trcsa_RFXS(NTSA,NU(NY,NX),NY,NX)
  ENDDO
  end subroutine AccumHourlyMicMacPoreFlux
!------------------------------------------------------------------------------------------

  subroutine SoluteFluxBySurfaceOutflow(M,N1,N2)
!
!     Description:
!
  implicit none
  integer, intent(in) :: M,N1,N2
  real(r8) :: VFLW
  INTEGER :: NTSA
!     begin_execution

  IF(VOLWM(M,0,N2,N1).GT.ZEROS2(N2,N1))THEN
    VFLW=AMIN1(VFLWX,QRM(M,N2,N1)/VOLWM(M,0,N2,N1))
  ELSE
    VFLW=VFLWX
  ENDIF

  DO NTSA=idsa_beg,idsa_end
    trcsa_RQR0(NTSA,N2,N1)=VFLW*AZMAX1(trcsa_solml2(NTSA,0,N2,N1))
  ENDDO
  end subroutine SoluteFluxBySurfaceOutflow

!------------------------------------------------------------------------------------------

  subroutine UpdateSoluteInSurfNeighbors(M,N1,N2,NY,NX,NHW,NHE,NVN,NVS)
!
!     Description:
!
  implicit none
  integer, intent(in) :: M,N1,N2,NY,NX,NHW,NHE,NVN,NVS
  integer :: N,NN,N4,N5,N4B,N5B,NTSA
  real(r8) :: FQRM,VFLW
!     begin_execution
!     LOCATE INTERNAL BOUNDARIES BETWEEN ADJACENT GRID CELLS
!
  DO N=1,2
    DO  NN=1,2
      IF(N.EQ.1)THEN
        IF(NX.EQ.NHE.AND.NN.EQ.1.OR.NX.EQ.NHW.AND.NN.EQ.2)THEN
          cycle
        ELSE
          N4=NX+1
          N5=NY
          N4B=NX-1
          N5B=NY
        ENDIF
      ELSEIF(N.EQ.2)THEN
        IF(NY.EQ.NVS.AND.NN.EQ.1.OR.NY.EQ.NVN.AND.NN.EQ.2)THEN
          cycle
        ELSE
          N4=NX
          N5=NY+1
          N4B=NX
          N5B=NY-1
        ENDIF
      ENDIF
!
!     QRM=runoff from watsub.f
!     RQR*=solute in runoff
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
!
!     IF OVERLAND FLOW IS FROM CURRENT TO ADJACENT GRID CELL
!
      IF(QRM(M,N2,N1).GT.ZEROS(N2,N1))THEN
        IF(NN.EQ.1)THEN
          FQRM=QRMN(M,N,2,N5,N4)/QRM(M,N2,N1)

          DO NTSA=idsa_beg,idsa_end
            trcsa_RQR(NTSA,N,2,N5,N4)=trcsa_RQR0(NTSA,N2,N1)*FQRM
          ENDDO
!
!     ACCUMULATE HOURLY FLUXES FOR USE IN REDIST.F
!
!     XQR*=hourly solute in runoff
!     RQR*=solute in runoff
!
          DO NTSA=idsa_beg,idsa_end
            trcsa_XQR(NTSA,N,2,N5,N4)=trcsa_XQR(NTSA,N,2,N5,N4)+trcsa_RQR(NTSA,N,2,N5,N4)
          ENDDO
        ELSE
          trcsa_RQR(idsa_beg:idsa_end,N,2,N5,N4)=0.0_r8
        ENDIF
!
!     IF OVERLAND FLOW IS FROM CURRENT TO ADJACENT GRID CELL
!
        IF(NN.EQ.2)THEN
          IF(N4B.GT.0.AND.N5B.GT.0)THEN
            FQRM=QRMN(M,N,1,N5B,N4B)/QRM(M,N2,N1)
            DO NTSA=idsa_beg,idsa_end
              trcsa_RQR(NTSA,N,1,N5B,N4B)=trcsa_RQR0(NTSA,N2,N1)*FQRM
            ENDDO
      !
!     ACCUMULATE HOURLY FLUXES FOR USE IN REDIST.F
!
!     XQR*=hourly solute in runoff
!     RQR*=solute in runoff
!
            DO NTSA=idsa_beg,idsa_end
              trcsa_XQR(NTSA,N,1,N5B,N4B)=trcsa_XQR(NTSA,N,1,N5B,N4B)+trcsa_RQR(NTSA,N,1,N5B,N4B)
            ENDDO
          ELSE
            trcsa_RQR(idsa_beg:idsa_end,N,1,N5B,N4B)=0.0_r8
          ENDIF
        ENDIF
      ELSE
        trcsa_RQR(idsa_beg:idsa_end,N,2,N5,N4)=0.0_r8

        IF(N4B.GT.0.AND.N5B.GT.0)THEN
          trcsa_RQR(idsa_beg:idsa_end,N,1,N5B,N4B)=0.0_r8
        ENDIF
      ENDIF
!
!     SNOW DRIFT
!
!     QSM=snow transfer from watsub.f
!     RQS*=solute flux in snow transfer
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
!     VOLS=volume of snowpack from watsub.f
!     *W2=solute content of snowpack
!
      IF(NN.EQ.1)THEN
!
!     IF NO SNOW DRIFT THEN NO TRANSPORT
!
        IF(ABS(QSM(M,N,N5,N4)).LE.ZEROS2(N2,N1))THEN
          trcsa_RQ(idsa_beg:idsab_end,N,N5,N4)=0.0_r8
!
!     IF DRIFT IS FROM CURRENT TO ADJACENT GRID CELL
!
        ELSEIF(QSM(M,N,N5,N4).GT.0.0)THEN
          IF(VOLS(N2,N1).GT.ZEROS2(NY,NX))THEN
            VFLW=AZMAX1(AMIN1(VFLWX,QSM(M,N,N5,N4)/VOLS(N2,N1)))
          ELSE
            VFLW=VFLWX
          ENDIF

          DO NTSA=idsa_beg,idsab_end
            trcsa_RQ(NTSA,N,N5,N4)=VFLW*AZMAX1(trcsa_sosml2(NTSA,1,N2,N1))
          ENDDO
!
!     IF DRIFT IS TO CURRENT FROM ADJACENT GRID CELL
!
        ELSEIF(QSM(M,N,N5,N4).LT.-ZEROS2(N2,N1))THEN
          IF(VOLS(N5,N4).GT.ZEROS2(N5,N4))THEN
            VFLW=AZMIN1(AMAX1(-VFLWX,QSM(M,N,N5,N4)/VOLS(N5,N4)))
          ELSE
            VFLW=-VFLWX
          ENDIF
          DO NTSA=idsa_beg,idsab_end
            trcsa_RQ(NTSA,N,N5,N4)=VFLW*AZMAX1(trcsa_sosml2(NTSA,1,N5,N4))
          ENDDO
        ENDIF
!
!     ACCUMULATE HOURLY FLUXES FOR USE IN REDIST.F
!
!     XQS*=hourly solute in snow transfer
!     RQS*=solute in snow transfer
!
        DO NTSA=idsa_beg,idsa_end
          trcsa_XQS(NTSA,N,N5,N4)=trcsa_XQS(NTSA,N,N5,N4)+trcsa_RQ(NTSA,N,N5,N4)
        ENDDO
      ENDIF
    enddo
  enddo
  end subroutine UpdateSoluteInSurfNeighbors
!------------------------------------------------------------------------------------------
   subroutine SoluteAdvMicropore(M,N,N1,N2,N3,N4,N5,N6)
   implicit none
   integer, intent(in) :: M,N,N1,N2,N3,N4,N5,N6
   real(r8) :: VFLW
   integer :: NTSA
   real(r8) :: trcsa_RFL(idsa_beg:idsab_end)
!
!     IF MICROPORE WATER FLUX FROM 'WATSUB' IS FROM CURRENT TO
!     ADJACENT GRID CELL THEN CONVECTIVE TRANSPORT IS THE PRODUCT
!     OF WATER FLUX AND MICROPORE GAS OR SOLUTE CONCENTRATIONS
!     IN CURRENT GRID CELL
!
!     FLWM=water flux through soil micropore from watsub.f
!     VOLWM=micropore water-filled porosity from watsub.f
!     RFL*S=solute diffusive flux through micropore
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
!     *S2,*B2=micropore solute content in non-band,band
!
  IF(FLWM(M,N,N6,N5,N4).GT.0.0)THEN
!
!     IF MICROPORE WATER FLUX FROM 'WATSUB' IS FROM CURRENT TO
!     ADJACENT GRID CELL THEN CONVECTIVE TRANSPORT IS THE PRODUCT
!     OF WATER FLUX AND MICROPORE GAS OR SOLUTE CONCENTRATIONS
!     IN CURRENT GRID CELL
!
    IF(VOLWM(M,N3,N2,N1).GT.ZEROS2(N2,N1))THEN
      VFLW=AZMAX1(AMIN1(VFLWX,FLWM(M,N,N6,N5,N4)/VOLWM(M,N3,N2,N1)))
    ELSE
      VFLW=VFLWX
    ENDIF

    DO NTSA=idsa_beg,idsa_psoil_beg-1
        trcsa_RFL(NTSA)=VFLW*AZMAX1(trcsa_solml2(NTSA,N3,N2,N1))
    ENDDO

    DO NTSA=idsa_psoil_beg,idsa_psoil_end
      trcsa_RFL(NTSA)=VFLW*AZMAX1(trcsa_solml2(NTSA,N3,N2,N1))*trcs_VLN(ids_H1PO4,N3,N2,N1)
    ENDDO

    DO NTSA=idsa_pband_beg,idsa_pband_end
      trcsa_RFL(NTSA)=VFLW*AZMAX1(trcsa_solml2(NTSA,N3,N2,N1))*trcs_VLN(ids_H1PO4B,N3,N2,N1)
    ENDDO
!
!     IF MICROPORE WATER FLUX FROM 'WATSUB' IS TO CURRENT FROM
!     ADJACENT GRID CELL THEN CONVECTIVE TRANSPORT IS THE PRODUCT
!     OF WATER FLUX AND MICROPORE GAS OR SOLUTE CONCENTRATIONS
!     IN ADJACENT GRID CELL
!
  ELSE
    IF(VOLWM(M,N6,N5,N4).GT.ZEROS2(N5,N4))THEN
      VFLW=AZMIN1(AMAX1(-VFLWX,FLWM(M,N,N6,N5,N4)/VOLWM(M,N6,N5,N4)))
    ELSE
      VFLW=-VFLWX
    ENDIF
    DO NTSA=idsa_beg,idsa_psoil_beg-1
      trcsa_RFL(NTSA)=VFLW*AZMAX1(trcsa_solml2(NTSA,N6,N5,N4))
    ENDDO

    DO NTSA=idsa_psoil_beg,idsa_psoil_end
      trcsa_RFL(NTSA)=VFLW*AZMAX1(trcsa_solml2(NTSA,N6,N5,N4))*trcs_VLN(ids_H1PO4,N6,N5,N4)
    ENDDO

    DO NTSA=idsa_pband_beg,idsa_pband_end
      trcsa_RFL(NTSA)=VFLW*AZMAX1(trcsa_solml2(NTSA,N6,N5,N4))*trcs_VLN(ids_H1PO4B,N6,N5,N4)
    ENDDO
  ENDIF

!     RFL*=convective flux through micropores
!     DFV*=diffusive solute flux through micropores
  DO NTSA=idsa_beg,idsab_end
    trcsa_RFLS(NTSA,N,N6,N5,N4)=trcsa_RFL(NTSA)
  ENDDO
  end subroutine SoluteAdvMicropore
!------------------------------------------------------------------------------------------
  subroutine SoluteDifsMicropore(M,N,N1,N2,N3,N4,N5,N6,THETW1)
  implicit none
  integer, intent(in) :: M,N,N1,N2,N3,N4,N5,N6
  REAL(R8), INTENT(IN):: THETW1(JZ,JY,JX)
  REAL(R8) :: VLWPA1,VLWPB1,VLWPA2,VLWPB2
  real(r8) :: DLYR1,DLYR2,TORTL,DISPN,DIFPO,DIFAL,DIFFE,DIFCA,DIFMG,DIFNA,DIFKA
  real(r8) :: trcsa_DIFC(idsa_beg:idsa_psoil_beg-1)
  real(r8) :: trcsa_solCl1(idsa_beg:idsab_end)
  real(r8) :: trcsa_solCl2(idsa_beg:idsab_end)
  real(r8) :: trcsa_DFV(idsa_beg:idsab_end)
  integer :: NTSA
!     DIFFUSIVE FLUXES OF SOLUTES BETWEEN CURRENT AND
!     ADJACENT GRID CELL MICROPORES FROM AQUEOUS DIFFUSIVITIES
!     AND CONCENTRATION DIFFERENCES
!
  IF(THETW1(N3,N2,N1).GT.THETY(N3,N2,N1) &
    .AND.THETW1(N6,N5,N4).GT.THETY(N6,N5,N4) &
    .AND.VOLWM(M,N3,N2,N1).GT.ZEROS2(N2,N1) &
    .AND.VOLWM(M,N6,N5,N4).GT.ZEROS2(N5,N4))THEN
!
!     VOLWM=micropore water-filled porosity from watsub.f
!     THETW=volumetric water content
!
!     MICROPORE CONCENTRATIONS FROM WATER-FILLED POROSITY
!     IN CURRENT AND ADJACENT GRID CELLS
!
!     C*1,C*2=solute concentration in source,destination layer
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
!     *S2,*B2=soil solute content in non-band,band
!
    VLWPA1=VOLWM(M,N3,N2,N1)*trcs_VLN(ids_H1PO4,N3,N2,N1)
    VLWPB1=VOLWM(M,N3,N2,N1)*trcs_VLN(ids_H1PO4B,N3,N2,N1)
    VLWPA2=VOLWM(M,N6,N5,N4)*trcs_VLN(ids_H1PO4,N6,N5,N4)
    VLWPB2=VOLWM(M,N6,N5,N4)*trcs_VLN(ids_H1PO4B,N6,N5,N4)

    DO NTSA=idsa_beg,idsa_psoil_beg-1
      trcsa_solCl1(NTSA)=AZMAX1(trcsa_solml2(NTSA,N3,N2,N1)/VOLWM(M,N3,N2,N1))
    ENDDO

    IF(VLWPA1.GT.ZEROS(N2,N1))THEN
      DO NTSA=idsa_psoil_beg,idsa_psoil_end
        trcsa_solCl1(NTSA)=AZMAX1(trcsa_solml2(NTSA,N3,N2,N1)/VLWPA1)
      ENDDO
    ELSE
      trcsa_solCl1(idsa_psoil_beg:idsa_psoil_end)=0.0_r8
    ENDIF

    IF(VLWPB1.GT.ZEROS(N2,N1))THEN
      DO NTSA=idsa_pband_beg,idsa_pband_end
        trcsa_solCl1(NTSA)=AZMAX1(trcsa_solml2(NTSA,N3,N2,N1)/VLWPB1)
      ENDDO
    ELSE
      DO NTSA=0,idsa_nuts
        trcsa_solCl1(idsa_H0PO4B+NTSA)=trcsa_solCl1(idsa_H0PO4+NTSA)
      ENDDO
    ENDIF

    DO NTSA=idsa_beg,idsa_psoil_beg-1
      trcsa_solCl2(NTSA)=AZMAX1(trcsa_solml2(NTSA,N6,N5,N4)/VOLWM(M,N6,N5,N4))
    ENDDO

    IF(VLWPA2.GT.ZEROS(N5,N4))THEN
      DO NTSA=idsa_psoil_beg,idsa_psoil_end
        trcsa_solCl2(NTSA)=AZMAX1(trcsa_solml2(NTSA,N6,N5,N4)/VOLWM(M,N6,N5,N4))
      ENDDO
    ELSE
      trcsa_solCl2(idsa_psoil_beg:idsa_psoil_end)=0.0_r8
    ENDIF
    IF(VLWPB2.GT.ZEROS(N5,N4))THEN
      DO NTSA=idsa_pband_beg,idsa_pband_end
        trcsa_solCl2(NTSA)=AZMAX1(trcsa_solml2(NTSA,N6,N5,N4)/VOLWM(M,N6,N5,N4))
      ENDDO
    ELSE
      DO NTSA=0,idsa_nuts
        trcsa_solCl2(idsa_H0PO4B+NTSA)=trcsa_solCl2(idsa_H0PO4+NTSA)
      ENDDO
    ENDIF
!
!     DIFFUSIVITIES IN CURRENT AND ADJACENT GRID CELL MICROPORES
!
!     DLYR=soil layer thickness
!     TORT=micropore tortuosity from hour1.f
!     DISP=dispersivity parameter
!     FLWM=water flux through soil micropore from watsub.f
!     DIF*=aqueous diffusivity-dispersivity through micropore
!     *SGL2=solute diffusivity from hour1.f
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
!     XDPTH=cross-sectional area/distance between layers
!     C*1,C*2=micropore solute concentration in source,destination layer
!     DFV*=diffusive solute transfer through soil micropore
!     VLNH4,VLNO3,VLPO4=non-band NH4,NO3,PO4 volume fraction
!     VLNHB,VLNOB,VLPOB=band NH4,NO3,PO4 volume fraction
!

    DLYR1=AMAX1(ZERO2,DLYR(N,N3,N2,N1))
    DLYR2=AMAX1(ZERO2,DLYR(N,N6,N5,N4))
    TORTL=(TORT(M,N3,N2,N1)*DLYR1+TORT(M,N6,N5,N4)*DLYR2)/(DLYR1+DLYR2)

    DISPN=DISP(N,N6,N5,N4)*AMIN1(VFLWX,ABS(FLWM(M,N,N6,N5,N4)/AREA(N,N6,N5,N4)))
    DIFPO=(POSGL2(N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)

    DIFAL=(ALSGL2(N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)
    trcsa_DIFC((/idsa_Al,idsa_AlOH,idsa_AlOH2,idsa_AlOH3,idsa_AlOH4,idsa_AlSO4/))=DIFAL

    DIFFE=(FESGL2(N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)
    trcsa_DIFC((/idsa_Fe,idsa_FeOH,idsa_FeOH2,idsa_FeOH3,idsa_FeOH4,idsa_FeSO4/))=DIFFE

    trcsa_DIFC(idsa_Hp)=(HYSGL2(N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)

    DIFCA=(CASGL2(N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)
    trcsa_DIFC((/idsa_Ca,idsa_CaOH2,idsa_CO3,idsa_HCO3,idsa_CaSO4/))=DIFCA

    DIFMG=(GMSGL2(N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)
    trcsa_DIFC((/idsa_Mg,idsa_MgOH2,idsa_MgCO3,idsa_MgHCO3,idsa_SO4/))=DIFMG

    DIFNA=(ANSGL2(N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)
    trcsa_DIFC((/idsa_Na,idsa_NaCO3,idsa_NaSO4/))=DIFNA

    DIFKA=(AKSGL2(N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)
    trcsa_DIFC((/idsa_K,idsa_KSO4/))=DIFKA

    trcsa_DIFC(idsa_OH)=(OHSGL2(N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)
    trcsa_DIFC(idsa_SO4)=(SOSGL2(N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)
    trcsa_DIFC(idsa_Cl)=(CLSXL2(N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)
    trcsa_DIFC(idsa_CO3)=(C3SGL2(N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)
    trcsa_DIFC(idsa_HCO3)=(HCSGL2(N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)
!
!     DIFFUSIVE FLUXES BETWEEN CURRENT AND ADJACENT GRID CELL
!     MICROPORES
!

    DO NTSA=idsa_beg,idsa_psoil_beg-1
      trcsa_DFV(NTSA)=trcsa_DIFC(NTSA)*(trcsa_solCl1(NTSA)-trcsa_solCl2(NTSA))
    ENDDO

    DO NTSA=idsa_psoil_beg,idsa_psoiL_end
      trcsa_DFV(NTSA)=DIFPO*(trcsa_solCl1(NTSA)-trcsa_solCl2(NTSA))*trcs_VLN(ids_H1PO4,N6,N5,N4)
    ENDDO

    DO NTSA=idsa_pband_beg,idsa_pband_end
      trcsa_DFV(NTSA)=DIFPO*(trcsa_solCl1(NTSA)-trcsa_solCl2(NTSA))*trcs_VLN(ids_H1PO4B,N6,N5,N4)
    ENDDO
  ELSE
    trcsa_DFV(idsa_beg:idsab_end)=0.0_r8
  ENDIF
!     RFL*=convective flux through micropores
!     DFV*=diffusive solute flux through micropores

  DO NTSA=idsa_beg,idsab_end
    trcsa_RFLS(NTSA,N,N6,N5,N4)=trcsa_RFLS(NTSA,N,N6,N5,N4)+trcsa_DFV(NTSA)
  ENDDO

  end subroutine SoluteDifsMicropore
!----------------------------------------------------------------------
  subroutine SoluteAdvMacropore(M,N,N1,N2,N3,N4,N5,N6)
  implicit none
  integer, intent(in) :: M,N,N1,N2,N3,N4,N5,N6
  integer :: NTSA
  real(r8) :: trcsa_RFH(idsa_beg:idsab_end)
  real(r8) :: VFLW
!     FLWHM=water flux through soil macropore from watsub.f
!
  IF(FLWHM(M,N,N6,N5,N4).GT.0.0)THEN
!
!     IF MACROPORE WATER FLUX FROM 'WATSUB' IS FROM CURRENT TO
!     ADJACENT GRID CELL THEN CONVECTIVE TRANSPORT IS THE PRODUCT
!     OF WATER FLUX AND MACROPORE GAS OR SOLUTE CONCENTRATIONS
!     IN CURRENT GRID CELL
!
!     VOLWHM=macropore water-filled porosity from watsub.f
!     VOLWAH=macropore porosity
!     RFH*=solute diffusive flux through macropore
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
!     *SH2,*BH2=macropore solute content in non-band,band
!     R*FXS=convective + diffusive solute flux between macropores and micropores
!     VLNH4,VLNO3,VLPO4=non-band NH4,NO3,PO4 volume fraction
!     VLNHB,VLNOB,VLPOB=band NH4,NO3,PO4 volume fraction
!
    IF(VOLWHM(M,N3,N2,N1).GT.ZEROS2(N2,N1))THEN
      VFLW=AZMAX1(AMIN1(VFLWX,FLWHM(M,N,N6,N5,N4)/VOLWHM(M,N3,N2,N1)))
    ELSE
      VFLW=VFLWX
    ENDIF
!
!     ACCOUNT FOR MACROPORE-MICROPORE EXCHANGE IN VERTICAL FLUX
!
    IF(N.EQ.3.AND.VOLAH(N6,N5,N4).GT.VOLWHM(M,N6,N5,N4))THEN
      DO NTSA=idsa_beg,idsa_psoil_beg-1
        trcsa_RFH(NTSA)=VFLW*AZMAX1((trcsa_soHml2(NTSA,N3,N2,N1) &
          -AZMIN1(trcsa_RFXS(NTSA,NU(N2,N1),N2,N1))))
      ENDDO

      DO NTSA=idsa_psoil_beg,idsa_psoil_end
        trcsa_RFH(NTSA)=VFLW*AZMAX1((trcsa_soHml2(NTSA,N3,N2,N1) &
          -AZMIN1(trcsa_RFXS(NTSA,NU(N2,N1),N2,N1))))*trcs_VLN(ids_H1PO4,N3,N2,N1)
      ENDDO

      DO NTSA=idsa_pband_beg,idsa_pband_end
        trcsa_RFH(NTSA)=VFLW*AZMAX1((trcsa_soHml2(NTSA,N3,N2,N1) &
          -AZMIN1(trcsa_RFXS(NTSA,NU(N2,N1),N2,N1))))*trcs_VLN(ids_H1PO4B,N3,N2,N1)
      ENDDO
!
!     OTHERWISE
!
    ELSE

      DO NTSA=idsa_beg,idsa_psoil_beg-1
        trcsa_RFH(NTSA)=VFLW*AZMAX1(trcsa_soHml2(NTSA,N3,N2,N1))
      ENDDO

      DO NTSA=idsa_psoil_beg,idsa_psoil_end
        trcsa_RFH(NTSA)=VFLW*AZMAX1(trcsa_soHml2(NTSA,N3,N2,N1)) &
          *trcs_VLN(ids_H1PO4,N6,N5,N4)
      ENDDO

      DO NTSA=idsa_pband_beg,idsa_pband_end
        trcsa_RFH(NTSA)=VFLW*AZMAX1(trcsa_soHml2(NTSA,N3,N2,N1)) &
          *trcs_VLN(ids_H1PO4B,N6,N5,N4)
      ENDDO
    ENDIF
  ELSEIF(FLWHM(M,N,N6,N5,N4).LT.0.0)THEN
!
!     IF MACROPORE WATER FLUX FROM 'WATSUB' IS FROM ADJACENT TO
!     CURRENT GRID CELL THEN CONVECTIVE TRANSPORT IS THE PRODUCT
!     OF WATER FLUX AND MACROPORE GAS OR SOLUTE CONCENTRATIONS
!     IN ADJACENT GRID CELL
!
    IF(VOLWHM(M,N6,N5,N4).GT.ZEROS2(N5,N4))THEN
      VFLW=AZMIN1(AMAX1(-VFLWX,FLWHM(M,N,N6,N5,N4)/VOLWHM(M,N6,N5,N4)))
    ELSE
      VFLW=-VFLWX
    ENDIF

    DO NTSA=idsa_beg,idsa_psoil_beg-1
      trcsa_RFH(NTSA)=VFLW*AZMAX1(trcsa_soHml2(NTSA,N6,N5,N4))
    ENDDO
    DO NTSA=idsa_psoil_beg,idsa_psoil_end
      trcsa_RFH(NTSA)=VFLW*AZMAX1(trcsa_soHml2(NTSA,N6,N5,N4)) &
        *trcs_VLN(ids_H1PO4,N6,N5,N4)
    ENDDO

    DO NTSA=idsa_pband_beg,idsa_pband_end
      trcsa_RFH(NTSA)=VFLW*AZMAX1(trcsa_soHml2(NTSA,N6,N5,N4)) &
        *trcs_VLN(ids_H1PO4B,N6,N5,N4)
    ENDDO
!
!     NO MACROPORE FLUX
!
  ELSE
    trcsa_RFH(idsa_beg:idsab_end)=0.0_r8
  ENDIF
!     RFH*=convective flux through macropores
!     DFH*=diffusive solute flux through macropores
!

  DO NTSA=idsa_beg,idsab_end
    trcsa_RFHS(NTSA,N,N6,N5,N4)=trcsa_RFH(NTSA)
  ENDDO

  end subroutine SoluteAdvMacropore

!----------------------------------------------------------------------
  subroutine SoluteDifsMacropore(M,N,N1,N2,N3,N4,N5,N6)
  implicit none
  integer, intent(in) :: M,N,N1,N2,N3,N4,N5,N6
  real(r8) :: DLYR1,DLYR2,TORTL,DISPN,DIFPO,DIFAL,DIFFE,DIFCA,DIFMG,DIFNA,DIFKA
  real(r8) :: trcsa_DIFC(idsa_beg:idsa_psoil_beg-1)
  real(r8) :: trcsa_DFH(idsa_beg:idsab_end)
  real(r8) :: trcsa_solCl1(idsa_beg:idsab_end)
  real(r8) :: trcsa_solCl2(idsa_beg:idsab_end)
  integer  :: NTSA

!     DIFFUSIVE FLUXES OF GASES AND SOLUTES BETWEEN CURRENT AND
!     ADJACENT GRID CELL MACROPORES FROM AQUEOUS DIFFUSIVITIES
!     AND CONCENTRATION DIFFERENCES
!
!     VOLWHM=macropore water-filled porosity from watsub.f
!     THETY=hygroscopic water content
!     VOLAH=total macropore volume

  IF(VOLWHM(M,N3,N2,N1).GT.THETY(N3,N2,N1)*VOLAH(N3,N2,N1) &
    .AND.VOLWHM(M,N6,N5,N4).GT.THETY(N6,N5,N4)*VOLAH(N6,N5,N4))THEN
!
!     MACROPORE CONCENTRATIONS IN CURRENT AND ADJACENT GRID CELLS
!
!     C*H1,C*H2=macropore solute concentration in source,destination layer
!     *H2=macropore solute content
!     VOLWHM=macropore water content
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
!
    DO NTSA=idsa_beg,idsab_end
      trcsa_solCl1(NTSA)=AZMAX1(trcsa_soHml2(NTSA,N3,N2,N1)/VOLWHM(M,N3,N2,N1))
      trcsa_solCl2(NTSA)=AZMAX1(trcsa_soHml2(NTSA,N6,N5,N4)/VOLWHM(M,N6,N5,N4))
    ENDDO
!
!     DIFFUSIVITIES IN CURRENT AND ADJACENT GRID CELL MACROPORES
!
!     DLYR=soil layer thickness
!     TORTH=macropore tortuosity from hour1.f
!     DISP=dispersivity parameter
!     FLWHM=water flux through soil macropore from watsub.f
!     DIF*=aqueous diffusivity-dispersivity through macropore
!     *SGL2=solute diffusivity from hour1.f
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
!     XDPTH=cross-sectional area/distance between layers
!     C*H1,C*H2=macropore solute concentration in source,destination layer
!     DFH*=diffusive solute transfer through soil macropore
!     VLNH4,VLNO3,VLPO4=non-band NH4,NO3,PO4 volume fraction
!     VLNHB,VLNOB,VLPOB=band NH4,NO3,PO4 volume fraction
!
    DLYR1=AMAX1(ZERO2,DLYR(N,N3,N2,N1))
    DLYR2=AMAX1(ZERO2,DLYR(N,N6,N5,N4))
    TORTL=(TORTH(M,N3,N2,N1)*DLYR1+TORTH(M,N6,N5,N4)*DLYR2)/(DLYR1+DLYR2)
    DISPN=DISP(N,N6,N5,N4)*AMIN1(VFLWX,ABS(FLWHM(M,N,N6,N5,N4)/AREA(N,N6,N5,N4)))
    DIFPO=(POSGL2(N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)

    DIFAL=(ALSGL2(N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)
    trcsa_DIFC((/idsa_Al,idsa_AlOH,idsa_AlOH2,idsa_AlOH3,idsa_AlOH4,idsa_AlSO4/))=DIFAL

    DIFFE=(FESGL2(N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)
    trcsa_DIFC((/idsa_Fe,idsa_FeOH,idsa_FeOH2,idsa_FeOH3,idsa_FeOH4,idsa_FeSO4/))=DIFFE

    trcsa_DIFC(idsa_Hp)=(HYSGL2(N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)

    DIFCA=(CASGL2(N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)
    trcsa_DIFC((/idsa_Ca,idsa_CaOH2,idsa_CO3,idsa_HCO3,idsa_CaSO4/))=DIFCA

    DIFMG=(GMSGL2(N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)
    trcsa_DIFC((/idsa_Mg,idsa_MgOH2,idsa_MgCO3,idsa_MgHCO3,idsa_SO4/))=DIFMG

    DIFNA=(ANSGL2(N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)
    trcsa_DIFC((/idsa_Na,idsa_NaCO3,idsa_NaSO4/))=DIFNA

    DIFKA=(AKSGL2(N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)
    trcsa_DIFC((/idsa_K,idsa_KSO4/))=DIFKA

    trcsa_DIFC(idsa_OH)=(OHSGL2(N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)
    trcsa_DIFC(idsa_SO4)=(SOSGL2(N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)
    trcsa_DIFC(idsa_Cl)=(CLSXL2(N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)
    trcsa_DIFC(idsa_CO3)=(C3SGL2(N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)
    trcsa_DIFC(idsa_HCO3)=(HCSGL2(N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)
!
!     DIFFUSIVE FLUXES BETWEEN CURRENT AND ADJACENT GRID CELL
!     MACROPORES
!

    DO NTSA=idsa_beg,idsa_psoil_beg-1
      trcsa_DFH(NTSA)=trcsa_DIFC(NTSA)*(trcsa_solCl1(NTSA)-trcsa_solCl2(NTSA))
    ENDDO

    DO NTSA=idsa_psoil_beg,idsa_psoil_end
      trcsa_DFH(NTSA)=DIFPO*(trcsa_solCl1(NTSA)-trcsa_solCl2(NTSA))*trcs_VLN(ids_H1PO4,N6,N5,N4)
    ENDDO

    DO NTSA=idsa_pband_beg,idsa_pband_end
      trcsa_DFH(NTSA)=DIFPO*(trcsa_solCl1(NTSA)-trcsa_solCl2(NTSA))*trcs_VLN(ids_H1PO4B,N6,N5,N4)
    ENDDO
  ELSE
    trcsa_DFH(idsa_beg:idsab_end)=0._r8
  ENDIF
!     RFH*=convective flux through macropores
!     DFH*=diffusive solute flux through macropores
!

  DO NTSA=idsa_beg,idsab_end
    trcsa_RFHS(NTSA,N,N6,N5,N4)=trcsa_RFHS(NTSA,N,N6,N5,N4)+trcsa_DFH(NTSA)
  ENDDO
  end subroutine SoluteDifsMacropore
!----------------------------------------------------------------------
  subroutine SoluteAdvExchMicMacpores(M,N,NY,NX,N1,N2,N3,N4,N5,N6)
  implicit none
  integer, intent(in) :: M,N,NY,NX,N1,N2,N3,N4,N5,N6

  real(r8) :: VFLW
  integer :: NTSA
  real(r8) :: trcsa_RFL(idsa_beg:idsab_end)
!     MACROPORE-MICROPORE CONVECTIVE SOLUTE EXCHANGE IN SOIL
!     LAYER FROM WATER EXCHANGE IN 'WATSUB' AND
!     FROM MACROPORE OR MICROPORE SOLUTE CONCENTRATIONS
!
!     FINHM=macro-micropore water transfer from watsub.f
!     VOLWM,VOLWHM=micropore,macropore water volume
!     RFL*=convective macropore-micropore solute transfer
!     VLNH4,VLNO3,VLPO4=non-band NH4,NO3,PO4 volume fraction
!     VLNHB,VLNOB,VLPOB=band NH4,NO3,PO4 volume fraction
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
!     *H2,*2=macropore,micropore solute content
!
!     MACROPORE TO MICROPORE TRANSFER
!
  IF(FINHM(M,N6,N5,N4).GT.0.0)THEN
    IF(VOLWHM(M,N6,N5,N4).GT.ZEROS2(NY,NX))THEN
      VFLW=AZMAX1(AMIN1(VFLWX,FINHM(M,N6,N5,N4)/VOLWHM(M,N6,N5,N4)))
    ELSE
      VFLW=VFLWX
    ENDIF
    DO NTSA=idsa_beg,idsa_psoil_beg-1
      trcsa_RFL(NTSA)=VFLW*AZMAX1(trcsa_soHml2(NTSA,N6,N5,N4))
    ENDDO

    DO NTSA=idsa_psoil_beg,idsa_psoil_end
      trcsa_RFL(NTSA)=VFLW*AZMAX1(trcsa_soHml2(NTSA,N6,N5,N4))*trcs_VLN(ids_H1PO4,N6,N5,N4)
    ENDDO

    DO NTSA=idsa_pband_beg,idsa_pband_end
      trcsa_RFL(NTSA)=VFLW*AZMAX1(trcsa_soHml2(NTSA,N6,N5,N4))*trcs_VLN(ids_H1PO4B,N6,N5,N4)
    ENDDO
!
!     MICROPORE TO MACROPORE TRANSFER
!
  ELSEIF(FINHM(M,N6,N5,N4).LT.0.0)THEN
    IF(VOLWM(M,N6,N5,N4).GT.ZEROS2(NY,NX))THEN
      VFLW=AZMIN1(AMAX1(-VFLWX,FINHM(M,N6,N5,N4)/VOLWM(M,N6,N5,N4)))
    ELSE
      VFLW=-VFLWX
    ENDIF
    DO NTSA=idsa_beg,idsa_psoil_beg-1
      trcsa_RFL(NTSA)=VFLW*AZMAX1(trcsa_solml2(NTSA,N6,N5,N4))
    ENDDO

    DO NTSA=idsa_psoil_beg,idsa_psoil_end
      trcsa_RFL(NTSA)=VFLW*AZMAX1(trcsa_solml2(NTSA,N6,N5,N4))*trcs_VLN(ids_H1PO4,N6,N5,N4)
    ENDDO

    DO NTSA=idsa_pband_beg,idsa_pband_end
      trcsa_RFL(NTSA)=VFLW*AZMAX1(trcsa_solml2(NTSA,N6,N5,N4))*trcs_VLN(ids_H1PO4,N6,N5,N4)
    ENDDO
!
!     NO MACROPORE TO MICROPORE TRANSFER
!
  ELSE
    trcsa_RFL(idsa_beg:idsab_end)=0.0_r8
  ENDIF
  DO NTSA=idsa_beg,idsab_end
    trcsa_RFXS(NTSA,N6,N5,N4)=trcsa_RFL(NTSA)
  ENDDO
  end subroutine SoluteAdvExchMicMacpores
!----------------------------------------------------------------------
  subroutine SoluteDifsExchMicMacpores(M,N,NY,NX,N1,N2,N3,N4,N5,N6)
  implicit none
  integer, intent(in) :: M,N,NY,NX,N1,N2,N3,N4,N5,N6
  real(r8) :: VOLWHS,VOLWT
  integer :: NTSA
  real(r8) :: trcsa_DFV(idsa_beg:idsab_end)
!     DIFFUSIVE FLUXES OF SOLUTES BETWEEN MICROPORES AND
!     MACROPORES FROM AQUEOUS DIFFUSIVITIES AND CONCENTRATION DIFFERENCES
!
!     VOLWM,VOLWHM=micropore,macropore water-filled porosity from watsub.f
!     DFV*S,DFV*B=diffusive solute flux between macro- and micropore in non-band,band
!     XNPH=1/no. of cycles h-1 for water, heat and solute flux calculations
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
!     *2,*H2=solute content of micropores,macropores
!
  IF(VOLWHM(M,N6,N5,N4).GT.ZEROS2(NY,NX))THEN
    VOLWHS=AMIN1(XFRS*VOLT(N6,N5,N4),VOLWHM(M,N6,N5,N4))
    VOLWT=VOLWM(M,N6,N5,N4)+VOLWHS
    DO NTSA=idsa_beg,idsa_psoil_beg-1
      trcsa_DFV(NTSA)=XNPH*(AZMAX1(trcsa_soHml2(NTSA,N6,N5,N4))*VOLWM(M,N6,N5,N4) &
        -AZMAX1(trcsa_solml2(NTSA,N6,N5,N4))*VOLWHS)/VOLWT
    ENDDO

    DO NTSA=idsa_psoil_beg,idsa_psoil_end
      trcsa_DFV(NTSA)=XNPH*(AZMAX1(trcsa_soHml2(NTSA,N6,N5,N4))*VOLWM(M,N6,N5,N4) &
        -AZMAX1(trcsa_solml2(NTSA,N6,N5,N4))*VOLWHS)/VOLWT*trcs_VLN(ids_H1PO4,N6,N5,N4)
    ENDDO

    DO NTSA=idsa_pband_beg,idsa_pband_end
      trcsa_DFV(NTSA)=XNPH*(AZMAX1(trcsa_soHml2(NTSA,N6,N5,N4))*VOLWM(M,N6,N5,N4) &
        -AZMAX1(trcsa_solml2(NTSA,N6,N5,N4))*VOLWHS)/VOLWT*trcs_VLN(ids_H1PO4B,N6,N5,N4)
    ENDDO
  ELSE
    trcsa_DFV(idsa_beg:idsab_end)=0.0_r8
  ENDIF
  DO NTSA=idsa_beg,idsab_end
    trcsa_RFXS(NTSA,N6,N5,N4)=trcsa_RFXS(NTSA,N6,N5,N4)+trcsa_DFV(NTSA)
  ENDDO
  end subroutine SoluteDifsExchMicMacpores
!----------------------------------------------------------------------
  subroutine SoluteAdvDifsMicMacpore(M,N,N1,N2,N3,N4,N5,N6,THETW1)
  implicit none
  integer, intent(in) :: M,N,N1,N2,N3,N4,N5,N6
  REAL(R8), INTENT(IN):: THETW1(JZ,JY,JX)
  integer :: NTSA
!     SOLUTE TRANSPORT IN MICROPORES

  call SoluteAdvMicropore(M,N,N1,N2,N3,N4,N5,N6)
!
  call SoluteDifsMicropore(M,N,N1,N2,N3,N4,N5,N6,THETW1)

!
!     SOLUTE TRANSPORT IN MACROPORES
  call SoluteAdvMacropore(M,N,N1,N2,N3,N4,N5,N6)
!
  call SoluteDifsMacropore(M,N,N1,N2,N3,N4,N5,N6)
!     TOTAL MICROPORE AND MACROPORE SOLUTE TRANSPORT
!     FLUXES = CONVECTIVE + DIFFUSIVE FLUXES
!
!     R*FLS=convective + diffusive solute flux through micropores
!     R*FLW,R*FLB=convective + diffusive solute flux through micropores in non-band,band
!     R*FHS=convective + diffusive solute flux through macropores
!     R*FHW,R*FHB=convective + diffusive solute flux through macropores in non-band,band
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
!
!     ACCUMULATE HOURLY FLUXES FOR USE IN REDIST.F
!
!     X*FLS=hourly convective + diffusive solute flux through micropores
!     X*FLW,X*FLB= hourly convective + diffusive solute flux through micropores in non-band,band
!     X*FHS=hourly convective + diffusive solute flux through macropores
!     X*FHW,X*FHB= hourly convective + diffusive solute flux through macropores in non-band,band
!     R*FLS=convective + diffusive solute flux through micropores
!     R*FLW,X*FLB=convective + diffusive solute flux through micropores in non-band,band
!     R*FHS=convective + diffusive solute flux through macropores
!     R*FHW,X*FHB=convective + diffusive solute flux through macropores in non-band,band
!

  DO NTSA=idsa_beg,idsab_end
    trcsa_XFLS(NTSA,N,N6,N5,N4)=trcsa_XFLS(NTSA,N,N6,N5,N4) &
      +trcsa_RFLS(NTSA,N,N6,N5,N4)
    trcsa_XFHS(NTSA,N,N6,N5,N4)=trcsa_XFHS(NTSA,N,N6,N5,N4) &
      +trcsa_RFHS(NTSA,N,N6,N5,N4)
  ENDDO
  end subroutine SoluteAdvDifsMicMacpore
!----------------------------------------------------------------------
  subroutine SoluteAdvDifsExchMicMacpore(M,N,NY,NX,N1,N2,N3,N4,N5,N6)
  implicit none
  integer, intent(in) :: M,N,NY,NX,N1,N2,N3,N4,N5,N6

  integer :: NTSA

  call SoluteAdvExchMicMacpores(M,N,NY,NX,N1,N2,N3,N4,N5,N6)
!
  call SoluteDifsExchMicMacpores(M,N,NY,NX,N1,N2,N3,N4,N5,N6)


!
!     TOTAL CONVECTIVE +DIFFUSIVE TRANSFER BETWEEN MACROPOES AND MICROPORES
!
!     R*FXS,R*FXB=total convective + diffusive solute flux between macro- and micropore in non-band,band
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
!     RFL*=convective flux between macro- and micropore
!     DFV*=diffusive solute flux between macro- and micropore
!

!
!     ACCUMULATE HOURLY FLUXES FOR USE IN REDIST.F
!
!     X*FXS,X*FXB= hourly convective + diffusive solute flux between macro- and micropore in non-band,band
!     R*FXS,R*FXB=convective + diffusive solute flux between macro- and micropore in non-band,band
!
  DO NTSA=idsa_beg,idsab_end
    trcsa_XFXS(NTSA,N6,N5,N4)=trcsa_XFXS(NTSA,N6,N5,N4)+trcsa_RFXS(NTSA,N6,N5,N4)
  ENDDO

  end subroutine SoluteAdvDifsExchMicMacpore
!----------------------------------------------------------------------
  subroutine UpdateSoluteInSubsurfNeighbors(M,NY,NX,NHW,NHE,NVN,NVS)
!
!     Description:
!
  implicit none
  integer, intent(in) :: M,NY,NX,NHW,NHE,NVN,NVS
  integer :: IFLGB,L,N1,N2,N3,N4,N5,N6
  integer :: LL,N
  real(r8) :: VOLWHS,VOLWT
  real(r8) :: VFLW,VLWPA1,VLWPB1,VLWPA2,VLWPB2
  real(r8) :: THETW1(JZ,JY,JX)

!     begin_execution
!     N3,N2,N1=L,NY,NX of source grid cell
!     N6,N5,N4=L,NY,NX of destination grid cell
!
  IFLGB=0
  DO  L=1,NL(NY,NX)
    N1=NX
    N2=NY
    N3=L
!
!     LOCATE INTERNAL BOUNDARIES BETWEEN ADJACENT GRID CELLS
!
    DO  N=NCN(N2,N1),3
      IF(N.EQ.1)THEN
        IF(NX.EQ.NHE)THEN
          cycle
        ELSE
          N4=NX+1
          N5=NY
          N6=L
        ENDIF
      ELSEIF(N.EQ.2)THEN
        IF(NY.EQ.NVS)THEN
          cycle
        ELSE
          N4=NX
          N5=NY+1
          N6=L
        ENDIF
      ELSEIF(N.EQ.3)THEN
        IF(L.EQ.NL(NY,NX))THEN
          cycle
        ELSE
          N4=NX
          N5=NY
          N6=L+1
        ENDIF
      ENDIF
      DO LL=N6,NL(NY,NX)
        IF(VOLX(LL,N5,N4).GT.ZEROS2(N5,N4))THEN
          N6=LL
          exit
        ENDIF
      enddo
!
!     SOLUTE FLUXES BETWEEN ADJACENT GRID CELLS FROM
!     WATER CONTENTS AND WATER FLUXES 'FLQM' FROM 'WATSUB'
!
!     VOLX,VOLY=soil volume excluding rock, macropore
!     VLNH4,VLNO3,VLPO4=non-band NH4,NO3,PO4 volume fraction
!     VLNHB,VLNOB,VLPOB=band NH4,NO3,PO4 volume fraction
!     VOLWM,VOLWHM=micropore,macropore water-filled porosity from watsub.f
!     THETW=volumetric water content
!     FLPM=change in air volume
!     XNPT=1/number of cycles NPH-1 for gas flux calculations
!
      IF(VOLX(N3,N2,N1).GT.ZEROS2(NY,NX))THEN
        IF(N3.GE.NUM(N2,N1).AND.N6.GE.NUM(N5,N4).AND.N3.LE.NL(N2,N1).AND.N6.LE.NL(N5,N4))THEN
          THETW1(N3,N2,N1)=AZMAX1(VOLWM(M,N3,N2,N1)/VOLY(N3,N2,N1))
          THETW1(N6,N5,N4)=AZMAX1(VOLWM(M,N6,N5,N4)/VOLY(N6,N5,N4))
!
          call SoluteAdvDifsMicMacpore(M,N,N1,N2,N3,N4,N5,N6,THETW1)

!
!     MACROPORE-MICROPORE SOLUTE EXCHANGE WITHIN SOIL
!     LAYER FROM WATER EXCHANGE IN 'WATSUB' AND
!     FROM MACROPORE OR MICROPORE SOLUTE CONCENTRATIONS
!
          IF(N.EQ.3)THEN
             call SoluteAdvDifsExchMicMacpore(M,N,NY,NX,N1,N2,N3,N4,N5,N6)
          ENDIF
        ELSEIF(N.NE.3)THEN
          THETW1(N3,N2,N1)=0.0_r8
          THETW1(N6,N5,N4)=0.0_r8
          trcsa_RFLS(idsa_beg:idsab_end,N,N6,N5,N4)=0.0_r8

          trcsa_RFHS(idsa_beg:idsab_end,N,N6,N5,N4)=0.0_r8
        ENDIF
      ELSE
        THETW1(N3,N2,N1)=0.0_r8
        THETW1(N6,N5,N4)=0.0_r8

        trcsa_RFLS(idsa_beg:idsab_end,N,N6,N5,N4)=0.0_r8
        trcsa_RFHS(idsa_beg:idsab_end,N,N6,N5,N4)=0.0_r8
      ENDIF
    enddo
  enddo
  end subroutine UpdateSoluteInSubsurfNeighbors
!------------------------------------------------------------------------------------------

  subroutine MacMicPoreFluxAdvPlusDifus(NY,NX,trcsa_DFV,trcsa_RFL)
!
!     Description:
!
  implicit none
  integer, intent(in) :: NY,NX
  real(r8), intent(in) :: trcsa_DFV(idsa_beg:idsab_end)
  real(r8), intent(in) :: trcsa_RFL(idsa_beg:idsab_end)
  integer :: NTSA
!     begin_execution
!
!     R*FXS=convective + diffusive solute flux between macropores and micropores
!     RFL*=convective flux between macropores and micropores
!     DFV*=diffusive solute flux between macropores and micropores
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
!
  DO NTSA=idsa_beg,idsab_end
    trcsa_RFXS(NTSA,NU(NY,NX),NY,NX)=trcsa_RFL(NTSA)+trcsa_DFV(NTSA)
  ENDDO
  end subroutine MacMicPoreFluxAdvPlusDifus
end module IngridTranspMod
