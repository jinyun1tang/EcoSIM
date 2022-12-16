module IngridTranspMod
  use data_kind_mod, only : r8 => SHR_KIND_R8
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

  real(r8) :: RFEFLS0,RHYFLS0,RCAFLS0,RMGFLS0,RNAFLS0,RKAFLS0
  real(r8) :: ROHFLS0,RSOFLS0,RCLFLS0,RC3FLS0,RHCFLS0,RAL1FS0
  real(r8) :: RAL2FS0,RAL3FS0,RAL4FS0,RALSFS0,RFE1FS0,RFE2FS0
  real(r8) :: RFE3FS0,RFE4FS0,RFESFS0,RCAOFS0,RCACFS0,RCAHFS0
  real(r8) :: RCASFS0,RMGOFS0,RMGCFS0,RMGHFS0,RMGSFS0,RNACFS0
  real(r8) :: RNASFS0,RKASFS0,RH0PFS0,RH3PFS0,RF1PFS0,RF2PFS0
  real(r8) :: RC0PFS0,RC1PFS0,RC2PFS0,RM1PFS0,RALFLS1,RFEFLS1
  real(r8) :: RHYFLS1,RCAFLS1,RMGFLS1,RNAFLS1,RKAFLS1,ROHFLS1
  real(r8) :: RSOFLS1,RCLFLS1,RC3FLS1,RHCFLS1,RAL1FS1,RAL2FS1
  real(r8) :: RAL3FS1,RAL4FS1,RALSFS1,RFE1FS1,RFE2FS1,RFE3FS1
  real(r8) :: RFE4FS1,RFESFS1,RCAOFS1,RCACFS1,RCAHFS1,RCASFS1
  real(r8) :: RMGOFS1,RMGCFS1,RMGHFS1,RMGSFS1,RNACFS1,RNASFS1
  real(r8) :: RKASFS1,RH0PFS1,RH3PFS1,RF1PFS1,RF2PFS1,RC0PFS1
  real(r8) :: RC1PFS1,RC2PFS1,RM1PFS1,RH0BFB1,RH3BFB1,RF1BFB1
  real(r8) :: RF2BFB1,RC0BFB1,RC1BFB1,RC2BFB1,RM1BFB1
  real(r8) :: RALFLS0

  real(r8) :: C0PO41,C3PO41,CFE1P1,CFE2P1,CCA0P1,CCA1P1,CCA2P1
  real(r8) :: CAL1,CFE1,CHY1,CCA1,CMG1,CNA1,CKA1,COH1,CSO41
  real(r8) :: CKAS2,C0PO42,C3PO42,CFE1P2,CFE2P2,CCA0P2,CCA1P2
  real(r8) :: CALS2,CFE12,CFE22,CFE32,CFE42,CFES2,CCAO2,CCAC2
  real(r8) :: CCL1,CCO31,CHCO31,CAL11,CAL21,CAL31,CAL41,CALS1
  real(r8) :: CFE11,CFE21,CFE31,CFE41,CFES1,CCAO1,CCAC1,CCAH1
  real(r8) :: CCAS1,CMGO1,CMGC1,CMGH1,CMGS1,CNAC1,CNAS1,CKAS1
  real(r8) :: CMG1P1,CAL2,CFE2,CHY2,CCA2,CMG2,CNA2,CKA2,COH2
  real(r8) :: CSO42,CCL2,CCO32,CHCO32,CAL12,CAL22,CAL32,CAL42
  real(r8) :: CCAH2,CCAS2,CMGO2,CMGC2,CMGH2,CMGS2,CNAC2,CNAS2
  real(r8) :: CCA2P2,CMG1P2,C0POB2,C3POB2,CF1PB2,CF2PB2,CC0PB2
  real(r8) :: CC1PB2,CC2PB2,CM1PB2,C0POB1,C3POB1,CF1PB1,CF2PB1
  real(r8) :: CC0PB1,CC1PB1,CC2PB1,CM1PB1
  real(r8) :: DFHC0B,DFHC1B,DFHC2B,DFHM1B,DLYR2,DFHAL,DFHFE
  real(r8) :: DFHHY,DFHCA,DFHMG,DFHNA,DFHKA,DFHOH,DFHSO,DFHCL
  real(r8) :: DFHC3,DFHHC,DFHAL1,DFHAL2,DFHAL3,DFHAL4,DFHALS
  real(r8) :: DFHFE1,DFHFE2,DFHFE3,DFHFE4,DFHFES,DFHCAO,DFHCAC
  real(r8) :: DFHCAH,DFHCAS,DFHMGO,DFHMGC,DFHMGH,DFHMGS,DFHNAC
  real(r8) :: DFHNAS,DFHKAS,DFHH0P,DFHH3P,DFHF1P,DFHF2P,DFHC0P
  real(r8) :: DFHC1P,DFHC2P,DFHM1P,DFHH0B,DFHH3B,DFHF1B,DFHF2B
  real(r8) :: DIFPO,DIFAL,DIFFE,DIFHY,DIFCA,DIFMG,DIFNA,DIFKA
  real(r8) :: DIFOH,DIFSO,DIFCL,DIFC3,DIFHC
  real(r8) :: DLYR0,DLYR1,DISPN

  real(r8) :: RFHAL,RFHFE
  real(r8) :: RFHHY,RFHCA,RFHMG,RFHNA,RFHKA,RFHOH,RFHSO,RFHCL
  real(r8) :: RFHC3,RFHHC,RFHAL1,RFHAL2,RFHAL3,RFHAL4,RFHALS
  real(r8) :: RFHFE1,RFHFE2,RFHFE3,RFHFE4,RFHFES,RFHCAO,RFHCAC
  real(r8) :: RFHCAH,RFHCAS,RFHMGO,RFHMGC,RFHMGH,RFHMGS,RFHNAC
  real(r8) :: RFHNAS,RFHKAS,RFHH0P,RFHH3P,RFHF1P,RFHF2P,RFHC0P
  real(r8) :: RFHC1P,RFHC2P,RFHM1P,RFHH0B,RFHH3B,RFHF1B,RFHF2B
  real(r8) :: RFHC0B,RFHC1B,RFHC2B,RFHM1B
  real(r8) :: TORT0,TORT1,TORTL

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
      call SoluteFluxInSnowpack(M,NY,NX)
!
!     CONVECTIVE SOLUTE EXCHANGE BETWEEN RESIDUE AND SOIL SURFACE
!
      FLWRM1=FLWRM(M,NY,NX)
!
!     FLWRM=litter-soil water flux from watsub.f

      IF(FLWRM1.GT.0.0_r8)THEN
        call Residue2TopsoilSoluteAdvExch(M,NY,NX,FLWRM1)
!
      ELSE
        call Topsoil2ResidueSoluteAdvExch(M,NY,NX,FLWRM1)
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

        call TopsoilResidueSolutedifusExch(M,NY,NX,FLWRM1)
      ELSE
        call ZeroDifusTopsoilResidueExch
      ENDIF
!
!     TOTAL MICROPORE AND MACROPORE SOLUTE TRANSPORT FLUXES BETWEEN
!     ADJACENT GRID CELLS = CONVECTIVE + DIFFUSIVE FLUXES

      call TopsoilResidueFluxAdvPlusDifus(NY,NX)
!
      call AccumHourlyTopsoilReisdueFlux(NY,NX)
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
        call ZeroMicMacPoreExchange
      ENDIF
!
!     DIFFUSIVE FLUXES OF SOLUTES BETWEEN MICROPORES AND
!     MACROPORES FROM AQUEOUS DIFFUSIVITIES AND CONCENTRATION DIFFERENCES
!
      IF(VOLWHM(M,NU(NY,NX),NY,NX).GT.ZEROS2(NY,NX))THEN
        call MacMicPoreSoluteDifusExchange(M,NY,NX)
      ELSE
        call ZeroMacMicPoreSoluteDifusExch
      ENDIF
!
!     TOTAL CONVECTIVE +DIFFUSIVE TRANSFER BETWEEN MACROPOES AND MICROPORES
!      call MacMicPoreFluxAdvPlusDifus(NY,NX)
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
        call ZeroSoluteFluxbySurfaceflow(N1,N2)
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

  subroutine SoluteFluxInSnowpack(M,NY,NX)
!
!     Description:
!
  implicit none
  integer, intent(in) :: M,NY,NX
  integer :: L,ICHKL,L2,NTSA
  real(r8) :: VFLWW,VFLWR,VFLWS,VFLWPO4,VFLWPOB

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
          trcsa_RBLS(NTSA,L2,NY,NX)=trcs_solsml2(NTSA,L,NY,NX)*VFLWW
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

          RALFLS0=trcs_solsml2(idsa_Al,L,NY,NX)*VFLWR
          RFEFLS0=trcs_solsml2(idsa_Fe,L,NY,NX)*VFLWR
          RHYFLS0=trcs_solsml2(idsa_Hp,L,NY,NX)*VFLWR
          RCAFLS0=trcs_solsml2(idsa_Ca,L,NY,NX)*VFLWR
          RMGFLS0=trcs_solsml2(idsa_Mg,L,NY,NX)*VFLWR
          RNAFLS0=trcs_solsml2(idsa_Na,L,NY,NX)*VFLWR
          RKAFLS0=trcs_solsml2(idsa_K,L,NY,NX)*VFLWR
          ROHFLS0=trcs_solsml2(idsa_OH,L,NY,NX)*VFLWR
          RSOFLS0=trcs_solsml2(idsa_SO4,L,NY,NX)*VFLWR
          RCLFLS0=trcs_solsml2(idsa_Cl,L,NY,NX)*VFLWR
          RC3FLS0=trcs_solsml2(idsa_CO3,L,NY,NX)*VFLWR
          RHCFLS0=trcs_solsml2(idsa_HCO3,L,NY,NX)*VFLWR
          RAL1FS0=trcs_solsml2(idsa_AlOH,L,NY,NX)*VFLWR
          RAL2FS0=trcs_solsml2(idsa_AlOH2,L,NY,NX)*VFLWR
          RAL3FS0=trcs_solsml2(idsa_AlOH3,L,NY,NX)*VFLWR
          RAL4FS0=trcs_solsml2(idsa_AlOH4,L,NY,NX)*VFLWR
          RALSFS0=trcs_solsml2(idsa_AlSO4,L,NY,NX)*VFLWR
          RFE1FS0=trcs_solsml2(idsa_FeOH,L,NY,NX)*VFLWR
          RFE2FS0=trcs_solsml2(idsa_FeOH2,L,NY,NX)*VFLWR
          RFE3FS0=trcs_solsml2(idsa_FeOH3,L,NY,NX)*VFLWR
          RFE4FS0=trcs_solsml2(idsa_FeOH4,L,NY,NX)*VFLWR
          RFESFS0=trcs_solsml2(idsa_FeSO4,L,NY,NX)*VFLWR
          RCAOFS0=trcs_solsml2(idsa_CaOH2,L,NY,NX)*VFLWR
          RCACFS0=trcs_solsml2(idsa_CaCO3,L,NY,NX)*VFLWR
          RCAHFS0=trcs_solsml2(idsa_CaHCO3,L,NY,NX)*VFLWR
          RCASFS0=trcs_solsml2(idsa_CaSO4,L,NY,NX)*VFLWR
          RMGOFS0=trcs_solsml2(idsa_MgOH2,L,NY,NX)*VFLWR
          RMGCFS0=trcs_solsml2(idsa_MgCO3,L,NY,NX)*VFLWR
          RMGHFS0=trcs_solsml2(idsa_MgHCO3,L,NY,NX)*VFLWR
          RMGSFS0=trcs_solsml2(idsa_MgSO4,L,NY,NX)*VFLWR
          RNACFS0=trcs_solsml2(idsa_NaCO3,L,NY,NX)*VFLWR
          RNASFS0=trcs_solsml2(idsa_NaSO4,L,NY,NX)*VFLWR
          RKASFS0=trcs_solsml2(idsa_KSO4,L,NY,NX)*VFLWR
          RH0PFS0=trcs_solsml2(idsa_H0PO4,L,NY,NX)*VFLWR
          RH3PFS0=trcs_solsml2(idsa_H3PO4,L,NY,NX)*VFLWR
          RF1PFS0=trcs_solsml2(idsa_FeHPO4,L,NY,NX)*VFLWR
          RF2PFS0=trcs_solsml2(idsa_FeH2PO4,L,NY,NX)*VFLWR
          RC0PFS0=trcs_solsml2(idsa_CaPO4,L,NY,NX)*VFLWR
          RC1PFS0=trcs_solsml2(idsa_CaHPO4,L,NY,NX)*VFLWR
          RC2PFS0=trcs_solsml2(idsa_CaH2PO4,L,NY,NX)*VFLWR
          RM1PFS0=trcs_solsml2(idsa_MgHPO4,L,NY,NX)*VFLWR

          RALFLS1=trcs_solsml2(idsa_Al,L,NY,NX)*VFLWS
          RFEFLS1=trcs_solsml2(idsa_Fe,L,NY,NX)*VFLWS
          RHYFLS1=trcs_solsml2(idsa_Hp,L,NY,NX)*VFLWS
          RCAFLS1=trcs_solsml2(idsa_Ca,L,NY,NX)*VFLWS
          RMGFLS1=trcs_solsml2(idsa_Mg,L,NY,NX)*VFLWS
          RNAFLS1=trcs_solsml2(idsa_Na,L,NY,NX)*VFLWS
          RKAFLS1=trcs_solsml2(idsa_K,L,NY,NX)*VFLWS
          ROHFLS1=trcs_solsml2(idsa_OH,L,NY,NX)*VFLWS
          RSOFLS1=trcs_solsml2(idsa_SO4,L,NY,NX)*VFLWS
          RCLFLS1=trcs_solsml2(idsa_Cl,L,NY,NX)*VFLWS
          RC3FLS1=trcs_solsml2(idsa_CO3,L,NY,NX)*VFLWS
          RHCFLS1=trcs_solsml2(idsa_HCO3,L,NY,NX)*VFLWS
          RAL1FS1=trcs_solsml2(idsa_AlOH,L,NY,NX)*VFLWS
          RAL2FS1=trcs_solsml2(idsa_AlOH2,L,NY,NX)*VFLWS
          RAL3FS1=trcs_solsml2(idsa_AlOH3,L,NY,NX)*VFLWS
          RAL4FS1=trcs_solsml2(idsa_AlOH4,L,NY,NX)*VFLWS
          RALSFS1=trcs_solsml2(idsa_AlSO4,L,NY,NX)*VFLWS
          RFE1FS1=trcs_solsml2(idsa_FeOH,L,NY,NX)*VFLWS
          RFE2FS1=trcs_solsml2(idsa_FeOH2,L,NY,NX)*VFLWS
          RFE3FS1=trcs_solsml2(idsa_FeOH3,L,NY,NX)*VFLWS
          RFE4FS1=trcs_solsml2(idsa_FeOH4,L,NY,NX)*VFLWS
          RFESFS1=trcs_solsml2(idsa_FeSO4,L,NY,NX)*VFLWS
          RCAOFS1=trcs_solsml2(idsa_CaOH2,L,NY,NX)*VFLWS
          RCACFS1=trcs_solsml2(idsa_CaCO3,L,NY,NX)*VFLWS
          RCAHFS1=trcs_solsml2(idsa_CaHCO3,L,NY,NX)*VFLWS
          RCASFS1=trcs_solsml2(idsa_CaSO4,L,NY,NX)*VFLWS
          RMGOFS1=trcs_solsml2(idsa_MgOH2,L,NY,NX)*VFLWS
          RMGCFS1=trcs_solsml2(idsa_MgCO3,L,NY,NX)*VFLWS
          RMGHFS1=trcs_solsml2(idsa_MgHCO3,L,NY,NX)*VFLWS
          RMGSFS1=trcs_solsml2(idsa_MgSO4,L,NY,NX)*VFLWS
          RNACFS1=trcs_solsml2(idsa_NaCO3,L,NY,NX)*VFLWS
          RNASFS1=trcs_solsml2(idsa_NaSO4,L,NY,NX)*VFLWS
          RKASFS1=trcs_solsml2(idsa_KSO4,L,NY,NX)*VFLWS

          RH0PFS1=trcs_solsml2(idsa_H0PO4,L,NY,NX)*VFLWPO4
          RH3PFS1=trcs_solsml2(idsa_H3PO4,L,NY,NX)*VFLWPO4
          RF1PFS1=trcs_solsml2(idsa_FeHPO4,L,NY,NX)*VFLWPO4
          RF2PFS1=trcs_solsml2(idsa_FeH2PO4,L,NY,NX)*VFLWPO4
          RC0PFS1=trcs_solsml2(idsa_CaPO4,L,NY,NX)*VFLWPO4
          RC1PFS1=trcs_solsml2(idsa_CaHPO4,L,NY,NX)*VFLWPO4
          RC2PFS1=trcs_solsml2(idsa_CaH2PO4,L,NY,NX)*VFLWPO4
          RM1PFS1=trcs_solsml2(idsa_MgHPO4,L,NY,NX)*VFLWPO4

          RH0BFB1=trcs_solsml2(idsa_H0PO4,L,NY,NX)*VFLWPOB
          RH3BFB1=trcs_solsml2(idsa_H3PO4,L,NY,NX)*VFLWPOB
          RF1BFB1=trcs_solsml2(idsa_FeHPO4,L,NY,NX)*VFLWPOB
          RF2BFB1=trcs_solsml2(idsa_FeH2PO4,L,NY,NX)*VFLWPOB
          RC0BFB1=trcs_solsml2(idsa_CaPO4,L,NY,NX)*VFLWPOB
          RC1BFB1=trcs_solsml2(idsa_CaHPO4,L,NY,NX)*VFLWPOB
          RC2BFB1=trcs_solsml2(idsa_CaH2PO4,L,NY,NX)*VFLWPOB
          RM1BFB1=trcs_solsml2(idsa_MgHPO4,L,NY,NX)*VFLWPOB
          ICHKL=1
        ENDIF
      ENDIF
    ELSE
      RALFLS0=0.0_r8
      RFEFLS0=0.0_r8
      RHYFLS0=0.0_r8
      RCAFLS0=0.0_r8
      RMGFLS0=0.0_r8
      RNAFLS0=0.0_r8
      RKAFLS0=0.0_r8
      ROHFLS0=0.0_r8
      RSOFLS0=0.0_r8
      RCLFLS0=0.0_r8
      RC3FLS0=0.0_r8
      RHCFLS0=0.0_r8
      RAL1FS0=0.0_r8
      RAL2FS0=0.0_r8
      RAL3FS0=0.0_r8
      RAL4FS0=0.0_r8
      RALSFS0=0.0_r8
      RFE1FS0=0.0_r8
      RFE2FS0=0.0_r8
      RFE3FS0=0.0_r8
      RFE4FS0=0.0_r8
      RFESFS0=0.0_r8
      RCAOFS0=0.0_r8
      RCACFS0=0.0_r8
      RCAHFS0=0.0_r8
      RCASFS0=0.0_r8
      RMGOFS0=0.0_r8
      RMGCFS0=0.0_r8
      RMGHFS0=0.0_r8
      RMGSFS0=0.0_r8
      RNACFS0=0.0_r8
      RNASFS0=0.0_r8
      RKASFS0=0.0_r8
      RH0PFS0=0.0_r8
      RH3PFS0=0.0_r8
      RF1PFS0=0.0_r8
      RF2PFS0=0.0_r8
      RC0PFS0=0.0_r8
      RC1PFS0=0.0_r8
      RC2PFS0=0.0_r8
      RM1PFS0=0.0_r8
      RALFLS1=0.0_r8
      RFEFLS1=0.0_r8
      RHYFLS1=0.0_r8
      RCAFLS1=0.0_r8
      RMGFLS1=0.0_r8
      RNAFLS1=0.0_r8
      RKAFLS1=0.0_r8
      ROHFLS1=0.0_r8
      RSOFLS1=0.0_r8
      RCLFLS1=0.0_r8
      RC3FLS1=0.0_r8
      RHCFLS1=0.0_r8
      RAL1FS1=0.0_r8
      RAL2FS1=0.0_r8
      RAL3FS1=0.0_r8
      RAL4FS1=0.0_r8
      RALSFS1=0.0_r8
      RFE1FS1=0.0_r8
      RFE2FS1=0.0_r8
      RFE3FS1=0.0_r8
      RFE4FS1=0.0_r8
      RFESFS1=0.0_r8
      RCAOFS1=0.0_r8
      RCACFS1=0.0_r8
      RCAHFS1=0.0_r8
      RCASFS1=0.0_r8
      RMGOFS1=0.0_r8
      RMGCFS1=0.0_r8
      RMGHFS1=0.0_r8
      RMGSFS1=0.0_r8
      RNACFS1=0.0_r8
      RNASFS1=0.0_r8
      RKASFS1=0.0_r8
      RH0PFS1=0.0_r8
      RH3PFS1=0.0_r8
      RF1PFS1=0.0_r8
      RF2PFS1=0.0_r8
      RC0PFS1=0.0_r8
      RC1PFS1=0.0_r8
      RC2PFS1=0.0_r8
      RM1PFS1=0.0_r8
      RH0BFB1=0.0_r8
      RH3BFB1=0.0_r8
      RF1BFB1=0.0_r8
      RF2BFB1=0.0_r8
      RC0BFB1=0.0_r8
      RC1BFB1=0.0_r8
      RC2BFB1=0.0_r8
      RM1BFB1=0.0_r8
    ENDIF
  ENDDO
  end subroutine SoluteFluxInSnowpack
!------------------------------------------------------------------------------------------

  subroutine Residue2TopsoilSoluteAdvExch(M,NY,NX,FLWRM1)
!
!     Description:
!
  implicit none
  integer, intent(in) :: M,NY,NX
  real(r8),intent(in) :: FLWRM1
  real(r8) :: VFLW

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
  RFLAL=VFLW*AZMAX1(trcsa_solml2(idsa_Al,0,NY,NX))
  RFLFE=VFLW*AZMAX1(trcsa_solml2(idsa_Fe,0,NY,NX))
  RFLHY=VFLW*AZMAX1(trcsa_solml2(idsa_Hp,0,NY,NX))
  RFLCA=VFLW*AZMAX1(trcsa_solml2(idsa_Ca,0,NY,NX))
  RFLMG=VFLW*AZMAX1(trcsa_solml2(idsa_Mg,0,NY,NX))
  RFLNA=VFLW*AZMAX1(trcsa_solml2(idsa_Na,0,NY,NX))
  RFLKA=VFLW*AZMAX1(trcsa_solml2(idsa_K,0,NY,NX))
  RFLOH=VFLW*AZMAX1(trcsa_solml2(idsa_OH,0,NY,NX))
  RFLSO=VFLW*AZMAX1(trcsa_solml2(idsa_SO4,0,NY,NX))
  RFLCL=VFLW*AZMAX1(trcsa_solml2(idsa_Cl,0,NY,NX))
  RFLC3=VFLW*AZMAX1(trcsa_solml2(idsa_CO3,0,NY,NX))
  RFLHC=VFLW*AZMAX1(trcsa_solml2(idsa_HCO3,0,NY,NX))
  RFLAL1=VFLW*AZMAX1(trcsa_solml2(idsa_AlOH,0,NY,NX))
  RFLAL2=VFLW*AZMAX1(trcsa_solml2(idsa_AlOH2,0,NY,NX))
  RFLAL3=VFLW*AZMAX1(trcsa_solml2(idsa_AlOH3,0,NY,NX))
  RFLAL4=VFLW*AZMAX1(trcsa_solml2(idsa_AlOH4,0,NY,NX))
  RFLALS=VFLW*AZMAX1(trcsa_solml2(idsa_AlSO4,0,NY,NX))
  RFLFE1=VFLW*AZMAX1(trcsa_solml2(idsa_FeOH,0,NY,NX))
  RFLFE2=VFLW*AZMAX1(trcsa_solml2(idsa_FeOH2,0,NY,NX))
  RFLFE3=VFLW*AZMAX1(trcsa_solml2(idsa_FeOH3,0,NY,NX))
  RFLFE4=VFLW*AZMAX1(trcsa_solml2(idsa_FeOH4,0,NY,NX))
  RFLFES=VFLW*AZMAX1(trcsa_solml2(idsa_FeSO4,0,NY,NX))
  RFLCAO=VFLW*AZMAX1(trcsa_solml2(idsa_CaOH2,0,NY,NX))
  RFLCAC=VFLW*AZMAX1(trcsa_solml2(idsa_CaCO3,0,NY,NX))
  RFLCAH=VFLW*AZMAX1(trcsa_solml2(idsa_CaHCO3,0,NY,NX))
  RFLCAS=VFLW*AZMAX1(trcsa_solml2(idsa_CaSO4,0,NY,NX))
  RFLMGO=VFLW*AZMAX1(trcsa_solml2(idsa_MgOH2,0,NY,NX))
  RFLMGC=VFLW*AZMAX1(trcsa_solml2(idsa_MgCO3,0,NY,NX))
  RFLMGH=VFLW*AZMAX1(trcsa_solml2(idsa_MgHCO3,0,NY,NX))
  RFLMGS=VFLW*AZMAX1(trcsa_solml2(idsa_MgSO4,0,NY,NX))
  RFLNAC=VFLW*AZMAX1(trcsa_solml2(idsa_NaCO3,0,NY,NX))
  RFLNAS=VFLW*AZMAX1(trcsa_solml2(idsa_NaSO4,0,NY,NX))
  RFLKAS=VFLW*AZMAX1(trcsa_solml2(idsa_KSO4,0,NY,NX))

  RFLH0P=VFLW*AZMAX1(trcsa_solml2(idsa_H0PO4,0,NY,NX))*trcs_VLN(ids_H1PO4,NU(NY,NX),NY,NX)
  RFLH3P=VFLW*AZMAX1(trcsa_solml2(idsa_H3PO4,0,NY,NX))*trcs_VLN(ids_H1PO4,NU(NY,NX),NY,NX)
  RFLF1P=VFLW*AZMAX1(trcsa_solml2(idsa_FeHPO4,0,NY,NX))*trcs_VLN(ids_H1PO4,NU(NY,NX),NY,NX)
  RFLF2P=VFLW*AZMAX1(trcsa_solml2(idsa_FeH2PO4,0,NY,NX))*trcs_VLN(ids_H1PO4,NU(NY,NX),NY,NX)
  RFLC0P=VFLW*AZMAX1(trcsa_solml2(idsa_CaPO4,0,NY,NX))*trcs_VLN(ids_H1PO4,NU(NY,NX),NY,NX)
  RFLC1P=VFLW*AZMAX1(trcsa_solml2(idsa_CaHPO4,0,NY,NX))*trcs_VLN(ids_H1PO4,NU(NY,NX),NY,NX)
  RFLC2P=VFLW*AZMAX1(trcsa_solml2(idsa_CaH2PO4,0,NY,NX))*trcs_VLN(ids_H1PO4,NU(NY,NX),NY,NX)
  RFLM1P=VFLW*AZMAX1(trcsa_solml2(idsa_MgHPO4,0,NY,NX))*trcs_VLN(ids_H1PO4,NU(NY,NX),NY,NX)

  RFLH0B=VFLW*AZMAX1(trcsa_solml2(idsa_H0PO4B,0,NY,NX))*trcs_VLN(ids_H1PO4B,NU(NY,NX),NY,NX)
  RFLH3B=VFLW*AZMAX1(trcsa_solml2(idsa_H3PO4B,0,NY,NX))*trcs_VLN(ids_H1PO4B,NU(NY,NX),NY,NX)
  RFLF1B=VFLW*AZMAX1(trcsa_solml2(idsa_FeHPO4B,0,NY,NX))*trcs_VLN(ids_H1PO4B,NU(NY,NX),NY,NX)
  RFLF2B=VFLW*AZMAX1(trcsa_solml2(idsa_FeH2PO4B,0,NY,NX))*trcs_VLN(ids_H1PO4B,NU(NY,NX),NY,NX)
  RFLC0B=VFLW*AZMAX1(trcsa_solml2(idsa_CaPO4B,0,NY,NX))*trcs_VLN(ids_H1PO4B,NU(NY,NX),NY,NX)
  RFLC1B=VFLW*AZMAX1(trcsa_solml2(idsa_CaHPO4B,0,NY,NX))*trcs_VLN(ids_H1PO4B,NU(NY,NX),NY,NX)
  RFLC2B=VFLW*AZMAX1(trcsa_solml2(idsa_CaH2PO4B,0,NY,NX))*trcs_VLN(ids_H1PO4B,NU(NY,NX),NY,NX)
  RFLM1B=VFLW*AZMAX1(trcsa_solml2(idsa_MgHPO4B,0,NY,NX))*trcs_VLN(ids_H1PO4B,NU(NY,NX),NY,NX)
  end subroutine Residue2TopsoilSoluteAdvExch
!------------------------------------------------------------------------------------------

  subroutine Topsoil2ResidueSoluteAdvExch(M,NY,NX,FLWRM1)
!
!     Description:
!
  implicit none
  integer, intent(in) :: M,NY,NX
  real(r8),intent(in) :: FLWRM1
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
  RFLAL=VFLW*AZMAX1(trcsa_solml2(idsa_Al,NU(NY,NX),NY,NX))
  RFLFE=VFLW*AZMAX1(trcsa_solml2(idsa_Fe,NU(NY,NX),NY,NX))
  RFLHY=VFLW*AZMAX1(trcsa_solml2(idsa_Hp,NU(NY,NX),NY,NX))
  RFLCA=VFLW*AZMAX1(trcsa_solml2(idsa_Ca,NU(NY,NX),NY,NX))
  RFLMG=VFLW*AZMAX1(trcsa_solml2(idsa_Mg,NU(NY,NX),NY,NX))
  RFLNA=VFLW*AZMAX1(trcsa_solml2(idsa_Na,NU(NY,NX),NY,NX))
  RFLKA=VFLW*AZMAX1(trcsa_solml2(idsa_K,NU(NY,NX),NY,NX))
  RFLOH=VFLW*AZMAX1(trcsa_solml2(idsa_OH,NU(NY,NX),NY,NX))
  RFLSO=VFLW*AZMAX1(trcsa_solml2(idsa_SO4,NU(NY,NX),NY,NX))
  RFLCL=VFLW*AZMAX1(trcsa_solml2(idsa_Cl,NU(NY,NX),NY,NX))
  RFLC3=VFLW*AZMAX1(trcsa_solml2(idsa_CO3,NU(NY,NX),NY,NX))
  RFLHC=VFLW*AZMAX1(trcsa_solml2(idsa_HCO3,NU(NY,NX),NY,NX))
  RFLAL1=VFLW*AZMAX1(trcsa_solml2(idsa_AlOH,NU(NY,NX),NY,NX))
  RFLAL2=VFLW*AZMAX1(trcsa_solml2(idsa_AlOH2,NU(NY,NX),NY,NX))
  RFLAL3=VFLW*AZMAX1(trcsa_solml2(idsa_AlOH3,NU(NY,NX),NY,NX))
  RFLAL4=VFLW*AZMAX1(trcsa_solml2(idsa_AlOH4,NU(NY,NX),NY,NX))
  RFLALS=VFLW*AZMAX1(trcsa_solml2(idsa_AlSO4,NU(NY,NX),NY,NX))
  RFLFE1=VFLW*AZMAX1(trcsa_solml2(idsa_FeOH,NU(NY,NX),NY,NX))
  RFLFE2=VFLW*AZMAX1(trcsa_solml2(idsa_FeOH2,NU(NY,NX),NY,NX))
  RFLFE3=VFLW*AZMAX1(trcsa_solml2(idsa_FeOH3,NU(NY,NX),NY,NX))
  RFLFE4=VFLW*AZMAX1(trcsa_solml2(idsa_FeOH4,NU(NY,NX),NY,NX))
  RFLFES=VFLW*AZMAX1(trcsa_solml2(idsa_FeSO4,NU(NY,NX),NY,NX))
  RFLCAO=VFLW*AZMAX1(trcsa_solml2(idsa_CaOH2,NU(NY,NX),NY,NX))
  RFLCAC=VFLW*AZMAX1(trcsa_solml2(idsa_CaCO3,NU(NY,NX),NY,NX))
  RFLCAH=VFLW*AZMAX1(trcsa_solml2(idsa_CaHCO3,NU(NY,NX),NY,NX))
  RFLCAS=VFLW*AZMAX1(trcsa_solml2(idsa_CaSO4,NU(NY,NX),NY,NX))
  RFLMGO=VFLW*AZMAX1(trcsa_solml2(idsa_MgOH2,NU(NY,NX),NY,NX))
  RFLMGC=VFLW*AZMAX1(trcsa_solml2(idsa_MgCO3,NU(NY,NX),NY,NX))
  RFLMGH=VFLW*AZMAX1(trcsa_solml2(idsa_MgHCO3,NU(NY,NX),NY,NX))
  RFLMGS=VFLW*AZMAX1(trcsa_solml2(idsa_MgSO4,NU(NY,NX),NY,NX))
  RFLNAC=VFLW*AZMAX1(trcsa_solml2(idsa_NaCO3,NU(NY,NX),NY,NX))
  RFLNAS=VFLW*AZMAX1(trcsa_solml2(idsa_NaSO4,NU(NY,NX),NY,NX))
  RFLKAS=VFLW*AZMAX1(trcsa_solml2(idsa_KSO4,NU(NY,NX),NY,NX))
  RFLH0P=VFLW*AZMAX1(trcsa_solml2(idsa_H0PO4,NU(NY,NX),NY,NX))
  RFLH3P=VFLW*AZMAX1(trcsa_solml2(idsa_H3PO4,NU(NY,NX),NY,NX))
  RFLF1P=VFLW*AZMAX1(trcsa_solml2(idsa_FeHPO4,NU(NY,NX),NY,NX))
  RFLF2P=VFLW*AZMAX1(trcsa_solml2(idsa_FeH2PO4,NU(NY,NX),NY,NX))
  RFLC0P=VFLW*AZMAX1(trcsa_solml2(idsa_CaPO4,NU(NY,NX),NY,NX))
  RFLC1P=VFLW*AZMAX1(trcsa_solml2(idsa_CaHPO4,NU(NY,NX),NY,NX))
  RFLC2P=VFLW*AZMAX1(trcsa_solml2(idsa_CaH2PO4,NU(NY,NX),NY,NX))
  RFLM1P=VFLW*AZMAX1(trcsa_solml2(idsa_MgHPO4,NU(NY,NX),NY,NX))
  RFLH0B=VFLW*AZMAX1(trcsa_solml2(idsa_H0PO4B,NU(NY,NX),NY,NX))
  RFLH3B=VFLW*AZMAX1(trcsa_solml2(idsa_H3PO4B,NU(NY,NX),NY,NX))
  RFLF1B=VFLW*AZMAX1(trcsa_solml2(idsa_FeHPO4B,NU(NY,NX),NY,NX))
  RFLF2B=VFLW*AZMAX1(trcsa_solml2(idsa_FeH2PO4B,NU(NY,NX),NY,NX))
  RFLC0B=VFLW*AZMAX1(trcsa_solml2(idsa_CaPO4B,NU(NY,NX),NY,NX))
  RFLC1B=VFLW*AZMAX1(trcsa_solml2(idsa_CaHPO4B,NU(NY,NX),NY,NX))
  RFLC2B=VFLW*AZMAX1(trcsa_solml2(idsa_CaH2PO4B,NU(NY,NX),NY,NX))
  RFLM1B=VFLW*AZMAX1(trcsa_solml2(idsa_MgHPO4B,NU(NY,NX),NY,NX))
  end subroutine Topsoil2ResidueSoluteAdvExch
!------------------------------------------------------------------------------------------

  subroutine TopsoilResidueSolutedifusExch(M,NY,NX,FLWRM1)
!
!     Description:
!
  implicit none
  integer, intent(in) :: M,NY,NX
  real(r8),intent(in) :: FLWRM1
  real(r8) :: VOLWPB,VOLWPA
!     begin_execution
!
!     MICROPORE CONCENTRATIONS FROM WATER IN RESIDUE AND SOIL SURFACE
!
!     C*1,C*2=solute concentration in litter, soil
!     Z*1,Z*2=solute content in litter, soil
!
  VOLWPA=VOLWM(M,NU(NY,NX),NY,NX)*trcs_VLN(ids_H1PO4,NU(NY,NX),NY,NX)
  VOLWPB=VOLWM(M,NU(NY,NX),NY,NX)*trcs_VLN(ids_H1PO4B,NU(NY,NX),NY,NX)
  CAL1=AZMAX1(trcsa_solml2(idsa_Al,0,NY,NX)/VOLWM(M,0,NY,NX))
  CFE1=AZMAX1(trcsa_solml2(idsa_Fe,0,NY,NX)/VOLWM(M,0,NY,NX))
  CHY1=AZMAX1(trcsa_solml2(idsa_Hp,0,NY,NX)/VOLWM(M,0,NY,NX))
  CCA1=AZMAX1(trcsa_solml2(idsa_Ca,0,NY,NX)/VOLWM(M,0,NY,NX))
  CMG1=AZMAX1(trcsa_solml2(idsa_Mg,0,NY,NX)/VOLWM(M,0,NY,NX))
  CNA1=AZMAX1(trcsa_solml2(idsa_Na,0,NY,NX)/VOLWM(M,0,NY,NX))
  CKA1=AZMAX1(trcsa_solml2(idsa_K,0,NY,NX)/VOLWM(M,0,NY,NX))
  COH1=AZMAX1(trcsa_solml2(idsa_OH,0,NY,NX)/VOLWM(M,0,NY,NX))
  CSO41=AZMAX1(trcsa_solml2(idsa_SO4,0,NY,NX)/VOLWM(M,0,NY,NX))
  CCL1=AZMAX1(trcsa_solml2(idsa_Cl,0,NY,NX)/VOLWM(M,0,NY,NX))
  CCO31=AZMAX1(trcsa_solml2(idsa_CO3,0,NY,NX)/VOLWM(M,0,NY,NX))
  CHCO31=AZMAX1(trcsa_solml2(idsa_HCO3,0,NY,NX)/VOLWM(M,0,NY,NX))
  CAL11=AZMAX1(trcsa_solml2(idsa_AlOH,0,NY,NX)/VOLWM(M,0,NY,NX))
  CAL21=AZMAX1(trcsa_solml2(idsa_AlOH2,0,NY,NX)/VOLWM(M,0,NY,NX))
  CAL31=AZMAX1(trcsa_solml2(idsa_AlOH3,0,NY,NX)/VOLWM(M,0,NY,NX))
  CAL41=AZMAX1(trcsa_solml2(idsa_AlOH4,0,NY,NX)/VOLWM(M,0,NY,NX))
  CALS1=AZMAX1(trcsa_solml2(idsa_AlSO4,0,NY,NX)/VOLWM(M,0,NY,NX))
  CFE11=AZMAX1(trcsa_solml2(idsa_FeOH,0,NY,NX)/VOLWM(M,0,NY,NX))
  CFE21=AZMAX1(trcsa_solml2(idsa_FeOH2,0,NY,NX)/VOLWM(M,0,NY,NX))
  CFE31=AZMAX1(trcsa_solml2(idsa_FeOH3,0,NY,NX)/VOLWM(M,0,NY,NX))
  CFE41=AZMAX1(trcsa_solml2(idsa_FeOH4,0,NY,NX)/VOLWM(M,0,NY,NX))
  CFES1=AZMAX1(trcsa_solml2(idsa_FeSO4,0,NY,NX)/VOLWM(M,0,NY,NX))
  CCAO1=AZMAX1(trcsa_solml2(idsa_CaOH2,0,NY,NX)/VOLWM(M,0,NY,NX))
  CCAC1=AZMAX1(trcsa_solml2(idsa_CaCO3,0,NY,NX)/VOLWM(M,0,NY,NX))
  CCAH1=AZMAX1(trcsa_solml2(idsa_CaHCO3,0,NY,NX)/VOLWM(M,0,NY,NX))
  CCAS1=AZMAX1(trcsa_solml2(idsa_CaSO4,0,NY,NX)/VOLWM(M,0,NY,NX))
  CMGO1=AZMAX1(trcsa_solml2(idsa_MgOH2,0,NY,NX)/VOLWM(M,0,NY,NX))
  CMGC1=AZMAX1(trcsa_solml2(idsa_MgCO3,0,NY,NX)/VOLWM(M,0,NY,NX))
  CMGH1=AZMAX1(trcsa_solml2(idsa_MgHCO3,0,NY,NX)/VOLWM(M,0,NY,NX))
  CMGS1=AZMAX1(trcsa_solml2(idsa_MgSO4,0,NY,NX)/VOLWM(M,0,NY,NX))
  CNAC1=AZMAX1(trcsa_solml2(idsa_NaCO3,0,NY,NX)/VOLWM(M,0,NY,NX))
  CNAS1=AZMAX1(trcsa_solml2(idsa_NaSO4,0,NY,NX)/VOLWM(M,0,NY,NX))
  CKAS1=AZMAX1(trcsa_solml2(idsa_KSO4,0,NY,NX)/VOLWM(M,0,NY,NX))
  C0PO41=AZMAX1(trcsa_solml2(idsa_H0PO4,0,NY,NX)/VOLWM(M,0,NY,NX))
  C3PO41=AZMAX1(trcsa_solml2(idsa_H3PO4,0,NY,NX)/VOLWM(M,0,NY,NX))
  CFE1P1=AZMAX1(trcsa_solml2(idsa_FeHPO4,0,NY,NX)/VOLWM(M,0,NY,NX))
  CFE2P1=AZMAX1(trcsa_solml2(idsa_FeH2PO4,0,NY,NX)/VOLWM(M,0,NY,NX))
  CCA0P1=AZMAX1(trcsa_solml2(idsa_CaPO4,0,NY,NX)/VOLWM(M,0,NY,NX))
  CCA1P1=AZMAX1(trcsa_solml2(idsa_CaHPO4,0,NY,NX)/VOLWM(M,0,NY,NX))
  CCA2P1=AZMAX1(trcsa_solml2(idsa_CaH2PO4,0,NY,NX)/VOLWM(M,0,NY,NX))
  CMG1P1=AZMAX1(trcsa_solml2(idsa_MgHPO4,0,NY,NX)/VOLWM(M,0,NY,NX))
  CAL2=AZMAX1(trcsa_solml2(idsa_Al,NU(NY,NX),NY,NX)/VOLWM(M,NU(NY,NX),NY,NX))
  CFE2=AZMAX1(trcsa_solml2(idsa_Fe,NU(NY,NX),NY,NX)/VOLWM(M,NU(NY,NX),NY,NX))
  CHY2=AZMAX1(trcsa_solml2(idsa_Hp,NU(NY,NX),NY,NX)/VOLWM(M,NU(NY,NX),NY,NX))
  CCA2=AZMAX1(trcsa_solml2(idsa_Ca,NU(NY,NX),NY,NX)/VOLWM(M,NU(NY,NX),NY,NX))
  CMG2=AZMAX1(trcsa_solml2(idsa_Mg,NU(NY,NX),NY,NX)/VOLWM(M,NU(NY,NX),NY,NX))
  CNA2=AZMAX1(trcsa_solml2(idsa_Na,NU(NY,NX),NY,NX)/VOLWM(M,NU(NY,NX),NY,NX))
  CKA2=AZMAX1(trcsa_solml2(idsa_K,NU(NY,NX),NY,NX)/VOLWM(M,NU(NY,NX),NY,NX))
  COH2=AZMAX1(trcsa_solml2(idsa_OH,NU(NY,NX),NY,NX)/VOLWM(M,NU(NY,NX),NY,NX))
  CSO42=AZMAX1(trcsa_solml2(idsa_SO4,NU(NY,NX),NY,NX)/VOLWM(M,NU(NY,NX),NY,NX))
  CCL2=AZMAX1(trcsa_solml2(idsa_Cl,NU(NY,NX),NY,NX)/VOLWM(M,NU(NY,NX),NY,NX))
  CCO32=AZMAX1(trcsa_solml2(idsa_CO3,NU(NY,NX),NY,NX)/VOLWM(M,NU(NY,NX),NY,NX))
  CHCO32=AZMAX1(trcsa_solml2(idsa_HCO3,NU(NY,NX),NY,NX)/VOLWM(M,NU(NY,NX),NY,NX))
  CAL12=AZMAX1(trcsa_solml2(idsa_AlOH,NU(NY,NX),NY,NX)/VOLWM(M,NU(NY,NX),NY,NX))
  CAL22=AZMAX1(trcsa_solml2(idsa_AlOH2,NU(NY,NX),NY,NX)/VOLWM(M,NU(NY,NX),NY,NX))
  CAL32=AZMAX1(trcsa_solml2(idsa_AlOH3,NU(NY,NX),NY,NX)/VOLWM(M,NU(NY,NX),NY,NX))
  CAL42=AZMAX1(trcsa_solml2(idsa_AlOH4,NU(NY,NX),NY,NX)/VOLWM(M,NU(NY,NX),NY,NX))
  CALS2=AZMAX1(trcsa_solml2(idsa_AlSO4,NU(NY,NX),NY,NX)/VOLWM(M,NU(NY,NX),NY,NX))
  CFE12=AZMAX1(trcsa_solml2(idsa_FeOH,NU(NY,NX),NY,NX)/VOLWM(M,NU(NY,NX),NY,NX))
  CFE22=AZMAX1(trcsa_solml2(idsa_FeOH2,NU(NY,NX),NY,NX)/VOLWM(M,NU(NY,NX),NY,NX))
  CFE32=AZMAX1(trcsa_solml2(idsa_FeOH3,NU(NY,NX),NY,NX)/VOLWM(M,NU(NY,NX),NY,NX))
  CFE42=AZMAX1(trcsa_solml2(idsa_FeOH4,NU(NY,NX),NY,NX)/VOLWM(M,NU(NY,NX),NY,NX))
  CFES2=AZMAX1(trcsa_solml2(idsa_FeSO4,NU(NY,NX),NY,NX)/VOLWM(M,NU(NY,NX),NY,NX))
  CCAO2=AZMAX1(trcsa_solml2(idsa_CaOH2,NU(NY,NX),NY,NX)/VOLWM(M,NU(NY,NX),NY,NX))
  CCAC2=AZMAX1(trcsa_solml2(idsa_CaCO3,NU(NY,NX),NY,NX)/VOLWM(M,NU(NY,NX),NY,NX))
  CCAH2=AZMAX1(trcsa_solml2(idsa_CaHCO3,NU(NY,NX),NY,NX)/VOLWM(M,NU(NY,NX),NY,NX))
  CCAS2=AZMAX1(trcsa_solml2(idsa_CaSO4,NU(NY,NX),NY,NX)/VOLWM(M,NU(NY,NX),NY,NX))
  CMGO2=AZMAX1(trcsa_solml2(idsa_MgOH2,NU(NY,NX),NY,NX)/VOLWM(M,NU(NY,NX),NY,NX))
  CMGC2=AZMAX1(trcsa_solml2(idsa_MgCO3,NU(NY,NX),NY,NX)/VOLWM(M,NU(NY,NX),NY,NX))
  CMGH2=AZMAX1(trcsa_solml2(idsa_MgHCO3,NU(NY,NX),NY,NX)/VOLWM(M,NU(NY,NX),NY,NX))
  CMGS2=AZMAX1(trcsa_solml2(idsa_MgSO4,NU(NY,NX),NY,NX)/VOLWM(M,NU(NY,NX),NY,NX))
  CNAC2=AZMAX1(trcsa_solml2(idsa_NaCO3,NU(NY,NX),NY,NX)/VOLWM(M,NU(NY,NX),NY,NX))
  CNAS2=AZMAX1(trcsa_solml2(idsa_NaSO4,NU(NY,NX),NY,NX)/VOLWM(M,NU(NY,NX),NY,NX))
  CKAS2=AZMAX1(trcsa_solml2(idsa_KSO4,NU(NY,NX),NY,NX)/VOLWM(M,NU(NY,NX),NY,NX))
  IF(VOLWPA.GT.ZEROS2(NY,NX))THEN
    C0PO42=AZMAX1(trcsa_solml2(idsa_H0PO4,NU(NY,NX),NY,NX)/VOLWPA)
    C3PO42=AZMAX1(trcsa_solml2(idsa_H3PO4,NU(NY,NX),NY,NX)/VOLWPA)
    CFE1P2=AZMAX1(trcsa_solml2(idsa_FeHPO4,NU(NY,NX),NY,NX)/VOLWPA)
    CFE2P2=AZMAX1(trcsa_solml2(idsa_FeH2PO4,NU(NY,NX),NY,NX)/VOLWPA)
    CCA0P2=AZMAX1(trcsa_solml2(idsa_CaPO4,NU(NY,NX),NY,NX)/VOLWPA)
    CCA1P2=AZMAX1(trcsa_solml2(idsa_CaHPO4,NU(NY,NX),NY,NX)/VOLWPA)
    CCA2P2=AZMAX1(trcsa_solml2(idsa_CaH2PO4,NU(NY,NX),NY,NX)/VOLWPA)
    CMG1P2=AZMAX1(trcsa_solml2(idsa_MgHPO4,NU(NY,NX),NY,NX)/VOLWPA)
  ELSE
    C0PO42=0.0_r8
    C3PO42=0.0_r8
    CFE1P2=0.0_r8
    CFE2P2=0.0_r8
    CCA0P2=0.0_r8
    CCA1P2=0.0_r8
    CCA2P2=0.0_r8
    CMG1P2=0.0_r8
  ENDIF
  IF(VOLWPB.GT.ZEROS2(NY,NX))THEN
    C0POB2=AZMAX1(trcsa_solml2(idsa_H0PO4B,NU(NY,NX),NY,NX)/VOLWPB)
    C3POB2=AZMAX1(trcsa_solml2(idsa_H3PO4B,NU(NY,NX),NY,NX)/VOLWPB)
    CF1PB2=AZMAX1(trcsa_solml2(idsa_FeHPO4B,NU(NY,NX),NY,NX)/VOLWPB)
    CF2PB2=AZMAX1(trcsa_solml2(idsa_FeH2PO4B,NU(NY,NX),NY,NX)/VOLWPB)
    CC0PB2=AZMAX1(trcsa_solml2(idsa_CaPO4B,NU(NY,NX),NY,NX)/VOLWPB)
    CC1PB2=AZMAX1(trcsa_solml2(idsa_CaHPO4B,NU(NY,NX),NY,NX)/VOLWPB)
    CC2PB2=AZMAX1(trcsa_solml2(idsa_CaH2PO4B,NU(NY,NX),NY,NX)/VOLWPB)
    CM1PB2=AZMAX1(trcsa_solml2(idsa_MgHPO4B,NU(NY,NX),NY,NX)/VOLWPB)
  ELSE
    C0POB2=C0PO42
    C3POB2=C3PO42
    CF1PB2=CFE1P2
    CF2PB2=CFE2P2
    CC0PB2=CCA0P2
    CC1PB2=CCA1P2
    CC2PB2=CCA2P2
    CM1PB2=CMG1P2
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
  DIFFE=(FESGL2(NU(NY,NX),NY,NX)*TORTL+DISPN)*AREA(3,NU(NY,NX),NY,NX)
  DIFHY=(HYSGL2(NU(NY,NX),NY,NX)*TORTL+DISPN)*AREA(3,NU(NY,NX),NY,NX)
  DIFCA=(CASGL2(NU(NY,NX),NY,NX)*TORTL+DISPN)*AREA(3,NU(NY,NX),NY,NX)
  DIFMG=(GMSGL2(NU(NY,NX),NY,NX)*TORTL+DISPN)*AREA(3,NU(NY,NX),NY,NX)
  DIFNA=(ANSGL2(NU(NY,NX),NY,NX)*TORTL+DISPN)*AREA(3,NU(NY,NX),NY,NX)
  DIFKA=(AKSGL2(NU(NY,NX),NY,NX)*TORTL+DISPN)*AREA(3,NU(NY,NX),NY,NX)
  DIFOH=(OHSGL2(NU(NY,NX),NY,NX)*TORTL+DISPN)*AREA(3,NU(NY,NX),NY,NX)
  DIFSO=(SOSGL2(NU(NY,NX),NY,NX)*TORTL+DISPN)*AREA(3,NU(NY,NX),NY,NX)
  DIFCL=(CLSXL2(NU(NY,NX),NY,NX)*TORTL+DISPN)*AREA(3,NU(NY,NX),NY,NX)
  DIFC3=(C3SGL2(NU(NY,NX),NY,NX)*TORTL+DISPN)*AREA(3,NU(NY,NX),NY,NX)
  DIFHC=(HCSGL2(NU(NY,NX),NY,NX)*TORTL+DISPN)*AREA(3,NU(NY,NX),NY,NX)
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
  DFVAL=DIFAL*(CAL1-CAL2)
  DFVFE=DIFFE*(CFE1-CFE2)
  DFVHY=DIFHY*(CHY1-CHY2)
  DFVCA=DIFCA*(CCA1-CCA2)
  DFVMG=DIFMG*(CMG1-CMG2)
  DFVNA=DIFNA*(CNA1-CNA2)
  DFVKA=DIFKA*(CKA1-CKA2)
  DFVOH=DIFOH*(COH1-COH2)
  DFVSO=DIFSO*(CSO41-CSO42)
  DFVCL=DIFCL*(CCL1-CCL2)
  DFVC3=DIFC3*(CCO31-CCO32)
  DFVHC=DIFHC*(CHCO31-CHCO32)
  DFVAL1=DIFAL*(CAL11-CAL12)
  DFVAL2=DIFAL*(CAL21-CAL22)
  DFVAL3=DIFAL*(CAL31-CAL32)
  DFVAL4=DIFAL*(CAL41-CAL42)
  DFVALS=DIFAL*(CALS1-CALS2)
  DFVFE1=DIFFE*(CFE11-CFE12)
  DFVFE2=DIFFE*(CFE21-CFE22)
  DFVFE3=DIFFE*(CFE31-CFE32)
  DFVFE4=DIFFE*(CFE41-CFE42)
  DFVFES=DIFFE*(CFES1-CFES2)
  DFVCAO=DIFCA*(CCAO1-CCAO2)
  DFVCAC=DIFCA*(CCAC1-CCAC2)
  DFVCAH=DIFCA*(CCAH1-CCAH2)
  DFVCAS=DIFCA*(CCAS1-CCAS2)
  DFVMGO=DIFMG*(CMGO1-CMGO2)
  DFVMGC=DIFMG*(CMGC1-CMGC2)
  DFVMGH=DIFMG*(CMGH1-CMGH2)
  DFVMGS=DIFMG*(CMGS1-CMGS2)
  DFVNAC=DIFNA*(CNAC1-CNAC2)
  DFVNAS=DIFNA*(CNAS1-CNAS2)
  DFVKAS=DIFKA*(CKAS1-CKAS2)
  DFVH0P=DIFPO*(C0PO41-C0PO42)*trcs_VLN(ids_H1PO4,NU(NY,NX),NY,NX)
  DFVH3P=DIFPO*(C3PO41-C3PO42)*trcs_VLN(ids_H1PO4,NU(NY,NX),NY,NX)
  DFVF1P=DIFPO*(CFE1P1-CFE1P2)*trcs_VLN(ids_H1PO4,NU(NY,NX),NY,NX)
  DFVF2P=DIFPO*(CFE2P1-CFE2P2)*trcs_VLN(ids_H1PO4,NU(NY,NX),NY,NX)
  DFVC0P=DIFPO*(CCA0P1-CCA0P2)*trcs_VLN(ids_H1PO4,NU(NY,NX),NY,NX)
  DFVC1P=DIFPO*(CCA1P1-CCA1P2)*trcs_VLN(ids_H1PO4,NU(NY,NX),NY,NX)
  DFVC2P=DIFPO*(CCA2P1-CCA2P2)*trcs_VLN(ids_H1PO4,NU(NY,NX),NY,NX)
  DFVM1P=DIFPO*(CMG1P1-CMG1P2)*trcs_VLN(ids_H1PO4,NU(NY,NX),NY,NX)
  DFVH0B=DIFPO*(C0PO41-C0POB2)*trcs_VLN(ids_H1PO4B,NU(NY,NX),NY,NX)
  DFVH3B=DIFPO*(C3PO41-C3POB2)*trcs_VLN(ids_H1PO4B,NU(NY,NX),NY,NX)
  DFVF1B=DIFPO*(CFE1P1-CF1PB2)*trcs_VLN(ids_H1PO4B,NU(NY,NX),NY,NX)
  DFVF2B=DIFPO*(CFE2P1-CF2PB2)*trcs_VLN(ids_H1PO4B,NU(NY,NX),NY,NX)
  DFVC0B=DIFPO*(CCA0P1-CC0PB2)*trcs_VLN(ids_H1PO4B,NU(NY,NX),NY,NX)
  DFVC1B=DIFPO*(CCA1P1-CC1PB2)*trcs_VLN(ids_H1PO4B,NU(NY,NX),NY,NX)
  DFVC2B=DIFPO*(CCA2P1-CC2PB2)*trcs_VLN(ids_H1PO4B,NU(NY,NX),NY,NX)
  DFVM1B=DIFPO*(CMG1P1-CM1PB2)*trcs_VLN(ids_H1PO4B,NU(NY,NX),NY,NX)
  end subroutine TopsoilResidueSolutedifusExch
!------------------------------------------------------------------------------------------

  subroutine ZeroDifusTopsoilResidueExch
!
!     Description:
!
  implicit none
!     begin_execution
  DFVAL=0.0_r8
  DFVFE=0.0_r8
  DFVHY=0.0_r8
  DFVCA=0.0_r8
  DFVMG=0.0_r8
  DFVNA=0.0_r8
  DFVKA=0.0_r8
  DFVOH=0.0_r8
  DFVSO=0.0_r8
  DFVCL=0.0_r8
  DFVC3=0.0_r8
  DFVHC=0.0_r8
  DFVAL1=0.0_r8
  DFVAL2=0.0_r8
  DFVAL3=0.0_r8
  DFVAL4=0.0_r8
  DFVALS=0.0_r8
  DFVFE1=0.0_r8
  DFVFE2=0.0_r8
  DFVFE3=0.0_r8
  DFVFE4=0.0_r8
  DFVFES=0.0_r8
  DFVCAO=0.0_r8
  DFVCAC=0.0_r8
  DFVCAH=0.0_r8
  DFVCAS=0.0_r8
  DFVMGO=0.0_r8
  DFVMGC=0.0_r8
  DFVMGH=0.0_r8
  DFVMGS=0.0_r8
  DFVNAC=0.0_r8
  DFVNAS=0.0_r8
  DFVKAS=0.0_r8
  DFVH0P=0.0_r8
  DFVH3P=0.0_r8
  DFVF1P=0.0_r8
  DFVF2P=0.0_r8
  DFVC0P=0.0_r8
  DFVC1P=0.0_r8
  DFVC2P=0.0_r8
  DFVM1P=0.0_r8
  DFVH0B=0.0_r8
  DFVH3B=0.0_r8
  DFVF1B=0.0_r8
  DFVF2B=0.0_r8
  DFVC0B=0.0_r8
  DFVC1B=0.0_r8
  DFVC2B=0.0_r8
  DFVM1B=0.0_r8
  end subroutine ZeroDifusTopsoilResidueExch
!------------------------------------------------------------------------------------------

  subroutine TopsoilResidueFluxAdvPlusDifus(NY,NX)
!
!     Description:
!
  implicit none
  integer, intent(in) :: NY,NX
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
  trcsa_RFLS(idsa_Al,3,0,NY,NX)=trcsa_RFL0(idsa_Al,NY,NX)+RALFLS0-RFLAL-DFVAL
  trcsa_RFLS(idsa_Fe,3,0,NY,NX)=trcsa_RFL0(idsa_Fe,NY,NX)+RFEFLS0-RFLFE-DFVFE
  trcsa_RFLS(idsa_Hp,3,0,NY,NX)=trcsa_RFL0(idsa_Hp,NY,NX)+RHYFLS0-RFLHY-DFVHY
  trcsa_RFLS(idsa_Ca,3,0,NY,NX)=trcsa_RFL0(idsa_Ca,NY,NX)+RCAFLS0-RFLCA-DFVCA
  trcsa_RFLS(idsa_Mg,3,0,NY,NX)=trcsa_RFL0(idsa_Mg,NY,NX)+RMGFLS0-RFLMG-DFVMG
  trcsa_RFLS(idsa_Na,3,0,NY,NX)=trcsa_RFL0(idsa_Na,NY,NX)+RNAFLS0-RFLNA-DFVNA
  trcsa_RFLS(idsa_K,3,0,NY,NX)=trcsa_RFL0(idsa_K,NY,NX)+RKAFLS0-RFLKA-DFVKA
  trcsa_RFLS(idsa_OH,3,0,NY,NX)=trcsa_RFL0(idsa_OH,NY,NX)+ROHFLS0-RFLOH-DFVOH
  trcsa_RFLS(idsa_SO4,3,0,NY,NX)=trcsa_RFL0(idsa_SO4,NY,NX)+RSOFLS0-RFLSO-DFVSO
  trcsa_RFLS(idsa_Cl,3,0,NY,NX)=trcsa_RFL0(idsa_Cl,NY,NX)+RCLFLS0-RFLCL-DFVCL
  trcsa_RFLS(idsa_CO3,3,0,NY,NX)=trcsa_RFL0(idsa_CO3,NY,NX)+RC3FLS0-RFLC3-DFVC3
  trcsa_RFLS(idsa_HCO3,3,0,NY,NX)=trcsa_RFL0(idsa_HCO3,NY,NX)+RHCFLS0-RFLHC-DFVHC
  trcsa_RFLS(idsa_AlOH,3,0,NY,NX)=trcsa_RFL0(idsa_AlOH,NY,NX)+RAL1FS0-RFLAL1-DFVAL1
  trcsa_RFLS(idsa_AlOH2,3,0,NY,NX)=trcsa_RFL0(idsa_AlOH2,NY,NX)+RAL2FS0-RFLAL2-DFVAL2
  trcsa_RFLS(idsa_AlOH3,3,0,NY,NX)=trcsa_RFL0(idsa_AlOH3,NY,NX)+RAL3FS0-RFLAL3-DFVAL3
  trcsa_RFLS(idsa_AlOH4,3,0,NY,NX)=trcsa_RFL0(idsa_AlOH4,NY,NX)+RAL4FS0-RFLAL4-DFVAL4
  trcsa_RFLS(idsa_AlSO4,3,0,NY,NX)=trcsa_RFL0(idsa_AlSO4,NY,NX)+RALSFS0-RFLALS-DFVALS
  trcsa_RFLS(idsa_FeOH,3,0,NY,NX)=trcsa_RFL0(idsa_FeOH,NY,NX)+RFE1FS0-RFLFE1-DFVFE1
  trcsa_RFLS(idsa_FeOH2,3,0,NY,NX)=trcsa_RFL0(idsa_FeOH2,NY,NX)+RFE2FS0-RFLFE2-DFVFE2
  trcsa_RFLS(idsa_FeOH3,3,0,NY,NX)=trcsa_RFL0(idsa_FeOH3,NY,NX)+RFE3FS0-RFLFE3-DFVFE3
  trcsa_RFLS(idsa_FeOH4,3,0,NY,NX)=trcsa_RFL0(idsa_FeOH4,NY,NX)+RFE4FS0-RFLFE4-DFVFE4
  trcsa_RFLS(idsa_FeSO4,3,0,NY,NX)=trcsa_RFL0(idsa_FeSO4,NY,NX)+RFESFS0-RFLFES-DFVFES
  trcsa_RFLS(idsa_CaOH2,3,0,NY,NX)=trcsa_RFL0(idsa_CaOH2,NY,NX)+RCAOFS0-RFLCAO-DFVCAO
  trcsa_RFLS(idsa_CaCO3,3,0,NY,NX)=trcsa_RFL0(idsa_CaCO3,NY,NX)+RCACFS0-RFLCAC-DFVCAC
  trcsa_RFLS(idsa_CaHCO3,3,0,NY,NX)=trcsa_RFL0(idsa_CaHCO3,NY,NX)+RCAHFS0-RFLCAH-DFVCAH
  trcsa_RFLS(idsa_CaSO4,3,0,NY,NX)=trcsa_RFL0(idsa_CaSO4,NY,NX)+RCASFS0-RFLCAS-DFVCAS
  trcsa_RFLS(idsa_MgOH2,3,0,NY,NX)=trcsa_RFL0(idsa_MgOH2,NY,NX)+RMGOFS0-RFLMGO-DFVMGO
  trcsa_RFLS(idsa_MgCO3,3,0,NY,NX)=trcsa_RFL0(idsa_MgCO3,NY,NX)+RMGCFS0-RFLMGC-DFVMGC
  trcsa_RFLS(idsa_MgHCO3,3,0,NY,NX)=trcsa_RFL0(idsa_MgHCO3,NY,NX)+RMGHFS0-RFLMGH-DFVMGH
  trcsa_RFLS(idsa_MgSO4,3,0,NY,NX)=trcsa_RFL0(idsa_MgSO4,NY,NX)+RMGSFS0-RFLMGS-DFVMGS
  trcsa_RFLS(idsa_NaCO3,3,0,NY,NX)=trcsa_RFL0(idsa_NaCO3,NY,NX)+RNACFS0-RFLNAC-DFVNAC
  trcsa_RFLS(idsa_NaSO4,3,0,NY,NX)=trcsa_RFL0(idsa_NaSO4,NY,NX)+RNASFS0-RFLNAS-DFVNAS
  trcsa_RFLS(idsa_KSO4,3,0,NY,NX)=trcsa_RFL0(idsa_KSO4,NY,NX)+RKASFS0-RFLKAS-DFVKAS
  trcsa_RFLS(idsa_H0PO4,3,0,NY,NX)=trcsa_RFL0(idsa_H0PO4,NY,NX)+RH0PFS0-RFLH0P-DFVH0P-RFLH0B-DFVH0B
  trcsa_RFLS(idsa_H3PO4,3,0,NY,NX)=trcsa_RFL0(idsa_H3PO4,NY,NX)+RH3PFS0-RFLH3P-DFVH3P-RFLH3B-DFVH3B
  trcsa_RFLS(idsa_FeHPO4,3,0,NY,NX)=trcsa_RFL0(idsa_FeHPO4,NY,NX)+RF1PFS0-RFLF1P-DFVF1P-RFLF1B-DFVF1B
  trcsa_RFLS(idsa_FeH2PO4,3,0,NY,NX)=trcsa_RFL0(idsa_FeH2PO4,NY,NX)+RF2PFS0-RFLF2P-DFVF2P-RFLF2B-DFVF2B
  trcsa_RFLS(idsa_CaPO4,3,0,NY,NX)=trcsa_RFL0(idsa_CaPO4,NY,NX)+RC0PFS0-RFLC0P-DFVC0P-RFLC0B-DFVC0B
  trcsa_RFLS(idsa_CaHPO4,3,0,NY,NX)=trcsa_RFL0(idsa_CaHPO4,NY,NX)+RC1PFS0-RFLC1P-DFVC1P-RFLC1B-DFVC1B
  trcsa_RFLS(idsa_CaH2PO4,3,0,NY,NX)=trcsa_RFL0(idsa_CaH2PO4,NY,NX)+RC2PFS0-RFLC2P-DFVC2P-RFLC2B-DFVC2B
  trcsa_RFLS(idsa_MgHPO4,3,0,NY,NX)=trcsa_RFL0(idsa_MgHPO4,NY,NX)+RM1PFS0-RFLM1P-DFVM1P-RFLM1B-DFVM1B

  trcsa_RFLS(idsa_Al,3,NU(NY,NX),NY,NX)=RALFL1(NY,NX)+RALFLS1+RFLAL+DFVAL
  trcsa_RFLS(idsa_Fe,3,NU(NY,NX),NY,NX)=RFEFL1(NY,NX)+RFEFLS1+RFLFE+DFVFE
  trcsa_RFLS(idsa_Hp,3,NU(NY,NX),NY,NX)=RHYFL1(NY,NX)+RHYFLS1+RFLHY+DFVHY
  trcsa_RFLS(idsa_Ca,3,NU(NY,NX),NY,NX)=RCAFL1(NY,NX)+RCAFLS1+RFLCA+DFVCA
  trcsa_RFLS(idsa_Mg,3,NU(NY,NX),NY,NX)=RMGFL1(NY,NX)+RMGFLS1+RFLMG+DFVMG
  trcsa_RFLS(idsa_Na,3,NU(NY,NX),NY,NX)=RNAFL1(NY,NX)+RNAFLS1+RFLNA+DFVNA
  trcsa_RFLS(idsa_K,3,NU(NY,NX),NY,NX)=RKAFL1(NY,NX)+RKAFLS1+RFLKA+DFVKA
  trcsa_RFLS(idsa_OH,3,NU(NY,NX),NY,NX)=ROHFL1(NY,NX)+ROHFLS1+RFLOH+DFVOH
  trcsa_RFLS(idsa_SO4,3,NU(NY,NX),NY,NX)=RSOFL1(NY,NX)+RSOFLS1+RFLSO+DFVSO
  trcsa_RFLS(idsa_Cl,3,NU(NY,NX),NY,NX)=RCLFL1(NY,NX)+RCLFLS1+RFLCL+DFVCL
  trcsa_RFLS(idsa_CO3,3,NU(NY,NX),NY,NX)=RC3FL1(NY,NX)+RC3FLS1+RFLC3+DFVC3
  trcsa_RFLS(idsa_HCO3,3,NU(NY,NX),NY,NX)=RHCFL1(NY,NX)+RHCFLS1+RFLHC+DFVHC
  trcsa_RFLS(idsa_AlOH,3,NU(NY,NX),NY,NX)=RAL1F1(NY,NX)+RAL1FS1+RFLAL1+DFVAL1
  trcsa_RFLS(idsa_AlOH2,3,NU(NY,NX),NY,NX)=RAL2F1(NY,NX)+RAL2FS1+RFLAL2+DFVAL2
  trcsa_RFLS(idsa_AlOH3,3,NU(NY,NX),NY,NX)=RAL3F1(NY,NX)+RAL3FS1+RFLAL3+DFVAL3
  trcsa_RFLS(idsa_AlOH4,3,NU(NY,NX),NY,NX)=RAL4F1(NY,NX)+RAL4FS1+RFLAL4+DFVAL4
  trcsa_RFLS(idsa_AlSO4,3,NU(NY,NX),NY,NX)=RALSF1(NY,NX)+RALSFS1+RFLALS+DFVALS
  trcsa_RFLS(idsa_FeOH,3,NU(NY,NX),NY,NX)=RFE1F1(NY,NX)+RFE1FS1+RFLFE1+DFVFE1
  trcsa_RFLS(idsa_FeOH2,3,NU(NY,NX),NY,NX)=RFE2F1(NY,NX)+RFE2FS1+RFLFE2+DFVFE2
  trcsa_RFLS(idsa_FeOH3,3,NU(NY,NX),NY,NX)=RFE3F1(NY,NX)+RFE3FS1+RFLFE3+DFVFE3
  trcsa_RFLS(idsa_FeOH4,3,NU(NY,NX),NY,NX)=RFE4F1(NY,NX)+RFE4FS1+RFLFE4+DFVFE4
  trcsa_RFLS(idsa_FeSO4,3,NU(NY,NX),NY,NX)=RFESF1(NY,NX)+RFESFS1+RFLFES+DFVFES
  trcsa_RFLS(idsa_CaOH2,3,NU(NY,NX),NY,NX)=RCAOF1(NY,NX)+RCAOFS1+RFLCAO+DFVCAO
  trcsa_RFLS(idsa_CaCO3,3,NU(NY,NX),NY,NX)=RCACF1(NY,NX)+RCACFS1+RFLCAC+DFVCAC
  trcsa_RFLS(idsa_CaHCO3,3,NU(NY,NX),NY,NX)=RCAHF1(NY,NX)+RCAHFS1+RFLCAH+DFVCAH
  trcsa_RFLS(idsa_CaSO4,3,NU(NY,NX),NY,NX)=RCASF1(NY,NX)+RCASFS1+RFLCAS+DFVCAS
  trcsa_RFLS(idsa_MgOH2,3,NU(NY,NX),NY,NX)=RMGOF1(NY,NX)+RMGOFS1+RFLMGO+DFVMGO
  trcsa_RFLS(idsa_MgCO3,3,NU(NY,NX),NY,NX)=RMGCF1(NY,NX)+RMGCFS1+RFLMGC+DFVMGC
  trcsa_RFLS(idsa_MgHCO3,3,NU(NY,NX),NY,NX)=RMGHF1(NY,NX)+RMGHFS1+RFLMGH+DFVMGH
  trcsa_RFLS(idsa_MgSO4,3,NU(NY,NX),NY,NX)=RMGSF1(NY,NX)+RMGSFS1+RFLMGS+DFVMGS
  trcsa_RFLS(idsa_NaCO3,3,NU(NY,NX),NY,NX)=RNACF1(NY,NX)+RNACFS1+RFLNAC+DFVNAC
  trcsa_RFLS(idsa_NaSO4,3,NU(NY,NX),NY,NX)=RNASF1(NY,NX)+RNASFS1+RFLNAS+DFVNAS
  trcsa_RFLS(idsa_KSO4,3,NU(NY,NX),NY,NX)=RKASF1(NY,NX)+RKASFS1+RFLKAS+DFVKAS
  trcsa_RFLS(idsa_H0PO4,3,NU(NY,NX),NY,NX)=RH0PF1(NY,NX)+RH0PFS1+RFLH0P+DFVH0P
  trcsa_RFLS(idsa_H3PO4,3,NU(NY,NX),NY,NX)=RH3PF1(NY,NX)+RH3PFS1+RFLH3P+DFVH3P
  trcsa_RFLS(idsa_FeHPO4,3,NU(NY,NX),NY,NX)=RF1PF1(NY,NX)+RF1PFS1+RFLF1P+DFVF1P
  trcsa_RFLS(idsa_FeH2PO4,3,NU(NY,NX),NY,NX)=RF2PF1(NY,NX)+RF2PFS1+RFLF2P+DFVF2P
  trcsa_RFLS(idsa_CaPO4,3,NU(NY,NX),NY,NX)=RC0PF1(NY,NX)+RC0PFS1+RFLC0P+DFVC0P
  trcsa_RFLS(idsa_CaHPO4,3,NU(NY,NX),NY,NX)=RC1PF1(NY,NX)+RC1PFS1+RFLC1P+DFVC1P
  trcsa_RFLS(idsa_CaH2PO4,3,NU(NY,NX),NY,NX)=RC2PF1(NY,NX)+RC2PFS1+RFLC2P+DFVC2P
  trcsa_RFLS(idsa_MgHPO4,3,NU(NY,NX),NY,NX)=RM1PF1(NY,NX)+RM1PFS1+RFLM1P+DFVM1P
  trcsa_RFLS(idsa_H0PO4B,3,NU(NY,NX),NY,NX)=RH0BF2(NY,NX)+RH0BFB1+RFLH0B+DFVH0B
  trcsa_RFLS(idsa_H3PO4B,3,NU(NY,NX),NY,NX)=RH3BF2(NY,NX)+RH3BFB1+RFLH3B+DFVH3B
  trcsa_RFLS(idsa_FeHPO4B,3,NU(NY,NX),NY,NX)=RF1BF2(NY,NX)+RF1BFB1+RFLF1B+DFVF1B
  trcsa_RFLS(idsa_FeH2PO4B,3,NU(NY,NX),NY,NX)=RF2BF2(NY,NX)+RF2BFB1+RFLF2B+DFVF2B
  trcsa_RFLS(idsa_CaPO4B,3,NU(NY,NX),NY,NX)=RC0BF2(NY,NX)+RC0BFB1+RFLC0B+DFVC0B
  trcsa_RFLS(idsa_CaHPO4B,3,NU(NY,NX),NY,NX)=RC1BF2(NY,NX)+RC1BFB1+RFLC1B+DFVC1B
  trcsa_RFLS(idsa_CaH2PO4B,3,NU(NY,NX),NY,NX)=RC2BF2(NY,NX)+RC2BFB1+RFLC2B+DFVC2B
  trcsa_RFLS(idsa_MgHPO4B,3,NU(NY,NX),NY,NX)=RM1BF2(NY,NX)+RM1BFB1+RFLM1B+DFVM1B
  end subroutine TopsoilResidueFluxAdvPlusDifus
!------------------------------------------------------------------------------------------

  subroutine AccumHourlyTopsoilReisdueFlux(NY,NX)
!
!     Description:
!
      implicit none
      integer, intent(in) :: NY,NX
!     begin_execution
!     ACCUMULATE HOURLY FLUXES FOR USE IN REDIST.F
!
!     X*FLS=hourly convective + diffusive solute flux
!     X*FLW,X*FLB= hourly convective + diffusive solute flux in non-band,band
!
  trcsa_XFLS(idsa_Al,3,0,NY,NX)=trcsa_XFLS(idsa_Al,3,0,NY,NX)+RALFLS0-RFLAL-DFVAL
  trcsa_XFLS(idsa_Fe,3,0,NY,NX)=trcsa_XFLS(idsa_Fe,3,0,NY,NX)+RFEFLS0-RFLFE-DFVFE
  trcsa_XFLS(idsa_Hp,3,0,NY,NX)=trcsa_XFLS(idsa_Hp,3,0,NY,NX)+RHYFLS0-RFLHY-DFVHY
  trcsa_XFLS(idsa_Ca,3,0,NY,NX)=trcsa_XFLS(idsa_Ca,3,0,NY,NX)+RCAFLS0-RFLCA-DFVCA
  trcsa_XFLS(idsa_Mg,3,0,NY,NX)=trcsa_XFLS(idsa_Mg,3,0,NY,NX)+RMGFLS0-RFLMG-DFVMG
  trcsa_XFLS(idsa_Na,3,0,NY,NX)=trcsa_XFLS(idsa_Na,3,0,NY,NX)+RNAFLS0-RFLNA-DFVNA
  trcsa_XFLS(idsa_K,3,0,NY,NX)=trcsa_XFLS(idsa_K,3,0,NY,NX)+RKAFLS0-RFLKA-DFVKA
  trcsa_XFLS(idsa_OH,3,0,NY,NX)=trcsa_XFLS(idsa_OH,3,0,NY,NX)+ROHFLS0-RFLOH-DFVOH
  trcsa_XFLS(idsa_SO4,3,0,NY,NX)=trcsa_XFLS(idsa_SO4,3,0,NY,NX)+RSOFLS0-RFLSO-DFVSO
  trcsa_XFLS(idsa_Cl,3,0,NY,NX)=trcsa_XFLS(idsa_Cl,3,0,NY,NX)+RCLFLS0-RFLCL-DFVCL
  trcsa_XFLS(idsa_CO3,3,0,NY,NX)=trcsa_XFLS(idsa_CO3,3,0,NY,NX)+RC3FLS0-RFLC3-DFVC3
  trcsa_XFLS(idsa_HCO3,3,0,NY,NX)=trcsa_XFLS(idsa_HCO3,3,0,NY,NX)+RHCFLS0-RFLHC-DFVHC
  trcsa_XFLS(idsa_AlOH,3,0,NY,NX)=trcsa_XFLS(idsa_AlOH,3,0,NY,NX)+RAL1FS0-RFLAL1-DFVAL1
  trcsa_XFLS(idsa_AlOH2,3,0,NY,NX)=trcsa_XFLS(idsa_AlOH2,3,0,NY,NX)+RAL2FS0-RFLAL2-DFVAL2
  trcsa_XFLS(idsa_AlOH3,3,0,NY,NX)=trcsa_XFLS(idsa_AlOH3,3,0,NY,NX)+RAL3FS0-RFLAL3-DFVAL3
  trcsa_XFLS(idsa_AlOH4,3,0,NY,NX)=trcsa_XFLS(idsa_AlOH4,3,0,NY,NX)+RAL4FS0-RFLAL4-DFVAL4
  trcsa_XFLS(idsa_AlSO4,3,0,NY,NX)=trcsa_XFLS(idsa_AlSO4,3,0,NY,NX)+RALSFS0-RFLALS-DFVALS
  trcsa_XFLS(idsa_FeOH,3,0,NY,NX)=trcsa_XFLS(idsa_FeOH,3,0,NY,NX)+RFE1FS0-RFLFE1-DFVFE1
  trcsa_XFLS(idsa_FeOH2,3,0,NY,NX)=trcsa_XFLS(idsa_FeOH2,3,0,NY,NX)+RFE2FS0-RFLFE2-DFVFE2
  trcsa_XFLS(idsa_FeOH3,3,0,NY,NX)=trcsa_XFLS(idsa_FeOH3,3,0,NY,NX)+RFE3FS0-RFLFE3-DFVFE3
  trcsa_XFLS(idsa_FeOH4,3,0,NY,NX)=trcsa_XFLS(idsa_FeOH4,3,0,NY,NX)+RFE4FS0-RFLFE4-DFVFE4
  trcsa_XFLS(idsa_FeSO4,3,0,NY,NX)=trcsa_XFLS(idsa_FeSO4,3,0,NY,NX)+RFESFS0-RFLFES-DFVFES
  trcsa_XFLS(idsa_CaOH2,3,0,NY,NX)=trcsa_XFLS(idsa_CaOH2,3,0,NY,NX)+RCAOFS0-RFLCAO-DFVCAO
  trcsa_XFLS(idsa_CaCO3,3,0,NY,NX)=trcsa_XFLS(idsa_CaCO3,3,0,NY,NX)+RCACFS0-RFLCAC-DFVCAC
  trcsa_XFLS(idsa_CaHCO3,3,0,NY,NX)=trcsa_XFLS(idsa_CaHCO3,3,0,NY,NX)+RCAHFS0-RFLCAH-DFVCAH
  trcsa_XFLS(idsa_CaSO4,3,0,NY,NX)=trcsa_XFLS(idsa_CaSO4,3,0,NY,NX)+RCASFS0-RFLCAS-DFVCAS
  trcsa_XFLS(idsa_MgOH2,3,0,NY,NX)=trcsa_XFLS(idsa_MgOH2,3,0,NY,NX)+RMGOFS0-RFLMGO-DFVMGO
  trcsa_XFLS(idsa_MgCO3,3,0,NY,NX)=trcsa_XFLS(idsa_MgCO3,3,0,NY,NX)+RMGCFS0-RFLMGC-DFVMGC
  trcsa_XFLS(idsa_MgHCO3,3,0,NY,NX)=trcsa_XFLS(idsa_MgHCO3,3,0,NY,NX)+RMGHFS0-RFLMGH-DFVMGH
  trcsa_XFLS(idsa_MgSO4,3,0,NY,NX)=trcsa_XFLS(idsa_MgSO4,3,0,NY,NX)+RMGSFS0-RFLMGS-DFVMGS
  trcsa_XFLS(idsa_NaCO3,3,0,NY,NX)=trcsa_XFLS(idsa_NaCO3,3,0,NY,NX)+RNACFS0-RFLNAC-DFVNAC
  trcsa_XFLS(idsa_NaSO4,3,0,NY,NX)=trcsa_XFLS(idsa_NaSO4,3,0,NY,NX)+RNASFS0-RFLNAS-DFVNAS
  trcsa_XFLS(idsa_KSO4,3,0,NY,NX)=trcsa_XFLS(idsa_KSO4,3,0,NY,NX)+RKASFS0-RFLKAS-DFVKAS
  trcsa_XFLS(idsa_H0PO4,3,0,NY,NX)=trcsa_XFLS(idsa_H0PO4,3,0,NY,NX)+RH0PFS0-RFLH0P-DFVH0P-RFLH0B-DFVH0B
  trcsa_XFLS(idsa_H3PO4,3,0,NY,NX)=trcsa_XFLS(idsa_H3PO4,3,0,NY,NX)+RH3PFS0-RFLH3P-DFVH3P-RFLH3B-DFVH3B
  trcsa_XFLS(idsa_FeHPO4,3,0,NY,NX)=trcsa_XFLS(idsa_FeHPO4,3,0,NY,NX)+RF1PFS0-RFLF1P-DFVF1P-RFLF1B-DFVF1B
  trcsa_XFLS(idsa_FeH2PO4,3,0,NY,NX)=trcsa_XFLS(idsa_FeH2PO4,3,0,NY,NX)+RF2PFS0-RFLF2P-DFVF2P-RFLF2B-DFVF2B
  trcsa_XFLS(idsa_CaPO4,3,0,NY,NX)=trcsa_XFLS(idsa_CaPO4,3,0,NY,NX)+RC0PFS0-RFLC0P-DFVC0P-RFLC0B-DFVC0B
  trcsa_XFLS(idsa_CaHPO4,3,0,NY,NX)=trcsa_XFLS(idsa_CaHPO4,3,0,NY,NX)+RC1PFS0-RFLC1P-DFVC1P-RFLC1B-DFVC1B
  trcsa_XFLS(idsa_CaH2PO4,3,0,NY,NX)=trcsa_XFLS(idsa_CaH2PO4,3,0,NY,NX)+RC2PFS0-RFLC2P-DFVC2P-RFLC2B-DFVC2B
  trcsa_XFLS(idsa_MgHPO4,3,0,NY,NX)=trcsa_XFLS(idsa_MgHPO4,3,0,NY,NX)+RM1PFS0-RFLM1P-DFVM1P-RFLM1B-DFVM1B

  trcsa_XFLS(idsa_Al,3,NU(NY,NX),NY,NX)=trcsa_XFLS(idsa_Al,3,NU(NY,NX),NY,NX)+RALFLS1+RFLAL+DFVAL
  trcsa_XFLS(idsa_Fe,3,NU(NY,NX),NY,NX)=trcsa_XFLS(idsa_Fe,3,NU(NY,NX),NY,NX)+RFEFLS1+RFLFE+DFVFE
  trcsa_XFLS(idsa_Hp,3,NU(NY,NX),NY,NX)=trcsa_XFLS(idsa_Hp,3,NU(NY,NX),NY,NX)+RHYFLS1+RFLHY+DFVHY
  trcsa_XFLS(idsa_Ca,3,NU(NY,NX),NY,NX)=trcsa_XFLS(idsa_Ca,3,NU(NY,NX),NY,NX)+RCAFLS1+RFLCA+DFVCA
  trcsa_XFLS(idsa_Mg,3,NU(NY,NX),NY,NX)=trcsa_XFLS(idsa_Mg,3,NU(NY,NX),NY,NX)+RMGFLS1+RFLMG+DFVMG
  trcsa_XFLS(idsa_Na,3,NU(NY,NX),NY,NX)=trcsa_XFLS(idsa_Na,3,NU(NY,NX),NY,NX)+RNAFLS1+RFLNA+DFVNA
  trcsa_XFLS(idsa_K,3,NU(NY,NX),NY,NX)=trcsa_XFLS(idsa_K,3,NU(NY,NX),NY,NX)+RKAFLS1+RFLKA+DFVKA
  trcsa_XFLS(idsa_OH,3,NU(NY,NX),NY,NX)=trcsa_XFLS(idsa_OH,3,NU(NY,NX),NY,NX)+ROHFLS1+RFLOH+DFVOH
  trcsa_XFLS(idsa_SO4,3,NU(NY,NX),NY,NX)=trcsa_XFLS(idsa_SO4,3,NU(NY,NX),NY,NX)+RSOFLS1+RFLSO+DFVSO
  trcsa_XFLS(idsa_Cl,3,NU(NY,NX),NY,NX)=trcsa_XFLS(idsa_Cl,3,NU(NY,NX),NY,NX)+RCLFLS1+RFLCL+DFVCL
  trcsa_XFLS(idsa_CO3,3,NU(NY,NX),NY,NX)=trcsa_XFLS(idsa_CO3,3,NU(NY,NX),NY,NX)+RC3FLS1+RFLC3+DFVC3
  trcsa_XFLS(idsa_HCO3,3,NU(NY,NX),NY,NX)=trcsa_XFLS(idsa_HCO3,3,NU(NY,NX),NY,NX)+RHCFLS1+RFLHC+DFVHC
  trcsa_XFLS(idsa_AlOH,3,NU(NY,NX),NY,NX)=trcsa_XFLS(idsa_AlOH,3,NU(NY,NX),NY,NX)+RAL1FS1+RFLAL1+DFVAL1
  trcsa_XFLS(idsa_AlOH2,3,NU(NY,NX),NY,NX)=trcsa_XFLS(idsa_AlOH2,3,NU(NY,NX),NY,NX)+RAL2FS1+RFLAL2+DFVAL2
  trcsa_XFLS(idsa_AlOH3,3,NU(NY,NX),NY,NX)=trcsa_XFLS(idsa_AlOH3,3,NU(NY,NX),NY,NX)+RAL3FS1+RFLAL3+DFVAL3
  trcsa_XFLS(idsa_AlOH4,3,NU(NY,NX),NY,NX)=trcsa_XFLS(idsa_AlOH4,3,NU(NY,NX),NY,NX)+RAL4FS1+RFLAL4+DFVAL4
  trcsa_XFLS(idsa_AlSO4,3,NU(NY,NX),NY,NX)=trcsa_XFLS(idsa_AlSO4,3,NU(NY,NX),NY,NX)+RALSFS1+RFLALS+DFVALS
  trcsa_XFLS(idsa_FeOH,3,NU(NY,NX),NY,NX)=trcsa_XFLS(idsa_FeOH,3,NU(NY,NX),NY,NX)+RFE1FS1+RFLFE1+DFVFE1
  trcsa_XFLS(idsa_FeOH2,3,NU(NY,NX),NY,NX)=trcsa_XFLS(idsa_FeOH2,3,NU(NY,NX),NY,NX)+RFE2FS1+RFLFE2+DFVFE2
  trcsa_XFLS(idsa_FeOH3,3,NU(NY,NX),NY,NX)=trcsa_XFLS(idsa_FeOH3,3,NU(NY,NX),NY,NX)+RFE3FS1+RFLFE3+DFVFE3

  trcsa_XFLS(idsa_FeOH4,3,NU(NY,NX),NY,NX)=trcsa_XFLS(idsa_FeOH4,3,NU(NY,NX),NY,NX)+RFE4FS1+RFLFE4+DFVFE4
  trcsa_XFLS(idsa_FeSO4,3,NU(NY,NX),NY,NX)=trcsa_XFLS(idsa_FeSO4,3,NU(NY,NX),NY,NX)+RFESFS1+RFLFES+DFVFES
  trcsa_XFLS(idsa_CaOH2,3,NU(NY,NX),NY,NX)=trcsa_XFLS(idsa_CaOH2,3,NU(NY,NX),NY,NX)+RCAOFS1+RFLCAO+DFVCAO
  trcsa_XFLS(idsa_CaCO3,3,NU(NY,NX),NY,NX)=trcsa_XFLS(idsa_CaCO3,3,NU(NY,NX),NY,NX)+RCACFS1+RFLCAC+DFVCAC
  trcsa_XFLS(idsa_CaHCO3,3,NU(NY,NX),NY,NX)=trcsa_XFLS(idsa_CaHCO3,3,NU(NY,NX),NY,NX)+RCAHFS1+RFLCAH+DFVCAH
  trcsa_XFLS(idsa_CaSO4,3,NU(NY,NX),NY,NX)=trcsa_XFLS(idsa_CaSO4,3,NU(NY,NX),NY,NX)+RCASFS1+RFLCAS+DFVCAS
  trcsa_XFLS(idsa_MgOH2,3,NU(NY,NX),NY,NX)=trcsa_XFLS(idsa_MgOH2,3,NU(NY,NX),NY,NX)+RMGOFS1+RFLMGO+DFVMGO
  trcsa_XFLS(idsa_MgCO3,3,NU(NY,NX),NY,NX)=trcsa_XFLS(idsa_MgCO3,3,NU(NY,NX),NY,NX)+RMGCFS1+RFLMGC+DFVMGC
  trcsa_XFLS(idsa_MgHCO3,3,NU(NY,NX),NY,NX)=trcsa_XFLS(idsa_MgHCO3,3,NU(NY,NX),NY,NX)+RMGHFS1+RFLMGH+DFVMGH
  trcsa_XFLS(idsa_MgSO4,3,NU(NY,NX),NY,NX)=trcsa_XFLS(idsa_MgSO4,3,NU(NY,NX),NY,NX)+RMGSFS1+RFLMGS+DFVMGS
  trcsa_XFLS(idsa_NaCO3,3,NU(NY,NX),NY,NX)=trcsa_XFLS(idsa_NaCO3,3,NU(NY,NX),NY,NX)+RNACFS1+RFLNAC+DFVNAC
  trcsa_XFLS(idsa_NaSO4,3,NU(NY,NX),NY,NX)=trcsa_XFLS(idsa_NaSO4,3,NU(NY,NX),NY,NX)+RNASFS1+RFLNAS+DFVNAS
  trcsa_XFLS(idsa_KSO4,3,NU(NY,NX),NY,NX)=trcsa_XFLS(idsa_KSO4,3,NU(NY,NX),NY,NX)+RKASFS1+RFLKAS+DFVKAS
  trcsa_XFLS(idsa_H0PO4,3,NU(NY,NX),NY,NX)=trcsa_XFLS(idsa_H0PO4,3,NU(NY,NX),NY,NX)+RH0PFS1+RFLH0P+DFVH0P
  trcsa_XFLS(idsa_H3PO4,3,NU(NY,NX),NY,NX)=trcsa_XFLS(idsa_H3PO4,3,NU(NY,NX),NY,NX)+RH3PFS1+RFLH3P+DFVH3P
  trcsa_XFLS(idsa_FeHPO4,3,NU(NY,NX),NY,NX)=trcsa_XFLS(idsa_FeHPO4,3,NU(NY,NX),NY,NX)+RF1PFS1+RFLF1P+DFVF1P
  trcsa_XFLS(idsa_FeH2PO4,3,NU(NY,NX),NY,NX)=trcsa_XFLS(idsa_FeH2PO4,3,NU(NY,NX),NY,NX)+RF2PFS1+RFLF2P+DFVF2P
  trcsa_XFLS(idsa_CaPO4,3,NU(NY,NX),NY,NX)=trcsa_XFLS(idsa_CaPO4,3,NU(NY,NX),NY,NX)+RC0PFS1+RFLC0P+DFVC0P
  trcsa_XFLS(idsa_CaHPO4,3,NU(NY,NX),NY,NX)=trcsa_XFLS(idsa_CaHPO4,3,NU(NY,NX),NY,NX)+RC1PFS1+RFLC1P+DFVC1P
  trcsa_XFLS(idsa_CaH2PO4,3,NU(NY,NX),NY,NX)=trcsa_XFLS(idsa_CaH2PO4,3,NU(NY,NX),NY,NX)+RC2PFS1+RFLC2P+DFVC2P
  trcsa_XFLS(idsa_MgHPO4,3,NU(NY,NX),NY,NX)=trcsa_XFLS(idsa_MgHPO4,3,NU(NY,NX),NY,NX)+RM1PFS1+RFLM1P+DFVM1P
  trcsa_XFLS(idsa_H0PO4B,3,NU(NY,NX),NY,NX)=trcsa_XFLS(idsa_H0PO4B,3,NU(NY,NX),NY,NX)+RH0BFB1+RFLH0B+DFVH0B
  trcsa_XFLS(idsa_H3PO4B,3,NU(NY,NX),NY,NX)=trcsa_XFLS(idsa_H3PO4B,3,NU(NY,NX),NY,NX)+RH3BFB1+RFLH3B+DFVH3B
  trcsa_XFLS(idsa_FeHPO4B,3,NU(NY,NX),NY,NX)=trcsa_XFLS(idsa_FeHPO4B,3,NU(NY,NX),NY,NX)+RF1BFB1+RFLF1B+DFVF1B
  trcsa_XFLS(idsa_FeH2PO4B,3,NU(NY,NX),NY,NX)=trcsa_XFLS(idsa_FeH2PO4B,3,NU(NY,NX),NY,NX)+RF2BFB1+RFLF2B+DFVF2B
  trcsa_XFLS(idsa_CaPO4B,3,NU(NY,NX),NY,NX)=trcsa_XFLS(idsa_CaPO4B,3,NU(NY,NX),NY,NX)+RC0BFB1+RFLC0B+DFVC0B
  trcsa_XFLS(idsa_CaHPO4B,3,NU(NY,NX),NY,NX)=trcsa_XFLS(idsa_CaHPO4B,3,NU(NY,NX),NY,NX)+RC1BFB1+RFLC1B+DFVC1B
  trcsa_XFLS(idsa_CaH2PO4B,3,NU(NY,NX),NY,NX)=trcsa_XFLS(idsa_CaH2PO4B,3,NU(NY,NX),NY,NX)+RC2BFB1+RFLC2B+DFVC2B
  trcsa_XFLS(idsa_MgHPO4B,3,NU(NY,NX),NY,NX)=trcsa_XFLS(idsa_MgHPO4B,3,NU(NY,NX),NY,NX)+RM1BFB1+RFLM1B+DFVM1B
  end subroutine AccumHourlyTopsoilReisdueFlux
!------------------------------------------------------------------------------------------

  subroutine MacToMicPoreSoluteAdvExchange(M,NY,NX)
!
!     Description:
!
  implicit none
  integer, intent(in) :: M,NY,NX
  real(r8) :: VFLW
!     begin_execution

  IF(VOLWHM(M,NU(NY,NX),NY,NX).GT.ZEROS2(NY,NX))THEN
    VFLW=AZMAX1(AMIN1(VFLWX,FINHM(M,NU(NY,NX),NY,NX)/VOLWHM(M,NU(NY,NX),NY,NX)))
  ELSE
    VFLW=VFLWX
  ENDIF
  RFLAL=VFLW*AZMAX1(trcsa_soHml2(idsa_Al,NU(NY,NX),NY,NX))
  RFLFE=VFLW*AZMAX1(trcsa_soHml2(idsa_Fe,NU(NY,NX),NY,NX))
  RFLHY=VFLW*AZMAX1(trcsa_soHml2(idsa_Hp,NU(NY,NX),NY,NX))
  RFLCA=VFLW*AZMAX1(trcsa_soHml2(idsa_Ca,NU(NY,NX),NY,NX))
  RFLMG=VFLW*AZMAX1(trcsa_soHml2(idsa_Mg,NU(NY,NX),NY,NX))
  RFLNA=VFLW*AZMAX1(trcsa_soHml2(idsa_Na,NU(NY,NX),NY,NX))
  RFLKA=VFLW*AZMAX1(trcsa_soHml2(idsa_K,NU(NY,NX),NY,NX))
  RFLOH=VFLW*AZMAX1(trcsa_soHml2(idsa_OH,NU(NY,NX),NY,NX))
  RFLSO=VFLW*AZMAX1(trcsa_soHml2(idsa_SO4,NU(NY,NX),NY,NX))
  RFLCL=VFLW*AZMAX1(trcsa_soHml2(idsa_Cl,NU(NY,NX),NY,NX))
  RFLC3=VFLW*AZMAX1(trcsa_soHml2(idsa_CO3,NU(NY,NX),NY,NX))
  RFLHC=VFLW*AZMAX1(trcsa_soHml2(idsa_HCO3,NU(NY,NX),NY,NX))
  RFLAL1=VFLW*AZMAX1(trcsa_soHml2(idsa_AlOH,NU(NY,NX),NY,NX))
  RFLAL2=VFLW*AZMAX1(trcsa_soHml2(idsa_AlOH2,NU(NY,NX),NY,NX))
  RFLAL3=VFLW*AZMAX1(trcsa_soHml2(idsa_AlOH3,NU(NY,NX),NY,NX))
  RFLAL4=VFLW*AZMAX1(trcsa_soHml2(idsa_AlOH4,NU(NY,NX),NY,NX))
  RFLALS=VFLW*AZMAX1(trcsa_soHml2(idsa_AlSO4,NU(NY,NX),NY,NX))
  RFLFE1=VFLW*AZMAX1(trcsa_soHml2(idsa_FeOH,NU(NY,NX),NY,NX))
  RFLFE2=VFLW*AZMAX1(trcsa_soHml2(idsa_FeOH2,NU(NY,NX),NY,NX))
  RFLFE3=VFLW*AZMAX1(trcsa_soHml2(idsa_FeOH3,NU(NY,NX),NY,NX))
  RFLFE4=VFLW*AZMAX1(trcsa_soHml2(idsa_FeOH4,NU(NY,NX),NY,NX))
  RFLFES=VFLW*AZMAX1(trcsa_soHml2(idsa_FeSO4,NU(NY,NX),NY,NX))
  RFLCAO=VFLW*AZMAX1(trcsa_soHml2(idsa_CaOH2,NU(NY,NX),NY,NX))
  RFLCAC=VFLW*AZMAX1(trcsa_soHml2(idsa_CaCO3,NU(NY,NX),NY,NX))
  RFLCAH=VFLW*AZMAX1(trcsa_soHml2(idsa_CaHCO3,NU(NY,NX),NY,NX))
  RFLCAS=VFLW*AZMAX1(trcsa_soHml2(idsa_CaSO4,NU(NY,NX),NY,NX))
  RFLMGO=VFLW*AZMAX1(trcsa_soHml2(idsa_MgOH2,NU(NY,NX),NY,NX))
  RFLMGC=VFLW*AZMAX1(trcsa_soHml2(idsa_MgCO3,NU(NY,NX),NY,NX))
  RFLMGH=VFLW*AZMAX1(trcsa_soHml2(idsa_MgHCO3,NU(NY,NX),NY,NX))
  RFLMGS=VFLW*AZMAX1(trcsa_soHml2(idsa_MgSO4,NU(NY,NX),NY,NX))
  RFLNAC=VFLW*AZMAX1(trcsa_soHml2(idsa_NaCO3,NU(NY,NX),NY,NX))
  RFLNAS=VFLW*AZMAX1(trcsa_soHml2(idsa_NaSO4,NU(NY,NX),NY,NX))
  RFLKAS=VFLW*AZMAX1(trcsa_soHml2(idsa_KSO4,NU(NY,NX),NY,NX))
  RFLH0P=VFLW*AZMAX1(trcsa_soHml2(idsa_H0PO4,NU(NY,NX),NY,NX))*trcs_VLN(ids_H1PO4,NU(NY,NX),NY,NX)
  RFLH3P=VFLW*AZMAX1(trcsa_soHml2(idsa_H3PO4,NU(NY,NX),NY,NX))*trcs_VLN(ids_H1PO4,NU(NY,NX),NY,NX)
  RFLF1P=VFLW*AZMAX1(trcsa_soHml2(idsa_FeHPO4,NU(NY,NX),NY,NX))*trcs_VLN(ids_H1PO4,NU(NY,NX),NY,NX)
  RFLF2P=VFLW*AZMAX1(trcsa_soHml2(idsa_FeH2PO4,NU(NY,NX),NY,NX))*trcs_VLN(ids_H1PO4,NU(NY,NX),NY,NX)
  RFLC0P=VFLW*AZMAX1(trcsa_soHml2(idsa_CaPO4,NU(NY,NX),NY,NX))*trcs_VLN(ids_H1PO4,NU(NY,NX),NY,NX)
  RFLC1P=VFLW*AZMAX1(trcsa_soHml2(idsa_CaHPO4,NU(NY,NX),NY,NX))*trcs_VLN(ids_H1PO4,NU(NY,NX),NY,NX)
  RFLC2P=VFLW*AZMAX1(trcsa_soHml2(idsa_CaH2PO4,NU(NY,NX),NY,NX))*trcs_VLN(ids_H1PO4,NU(NY,NX),NY,NX)
  RFLM1P=VFLW*AZMAX1(trcsa_soHml2(idsa_MgHPO4,NU(NY,NX),NY,NX))*trcs_VLN(ids_H1PO4,NU(NY,NX),NY,NX)
  RFLH0B=VFLW*AZMAX1(trcsa_soHml2(idsa_H0PO4B,NU(NY,NX),NY,NX))*trcs_VLN(ids_H1PO4B,NU(NY,NX),NY,NX)
  RFLH3B=VFLW*AZMAX1(trcsa_soHml2(idsa_H3PO4B,NU(NY,NX),NY,NX))*trcs_VLN(ids_H1PO4B,NU(NY,NX),NY,NX)
  RFLF1B=VFLW*AZMAX1(trcsa_soHml2(idsa_FeHPO4B,NU(NY,NX),NY,NX))*trcs_VLN(ids_H1PO4B,NU(NY,NX),NY,NX)
  RFLF2B=VFLW*AZMAX1(trcsa_soHml2(idsa_FeH2PO4B,NU(NY,NX),NY,NX))*trcs_VLN(ids_H1PO4B,NU(NY,NX),NY,NX)
  RFLC0B=VFLW*AZMAX1(trcsa_soHml2(idsa_CaPO4B,NU(NY,NX),NY,NX))*trcs_VLN(ids_H1PO4B,NU(NY,NX),NY,NX)
  RFLC1B=VFLW*AZMAX1(trcsa_soHml2(idsa_CaHPO4B,NU(NY,NX),NY,NX))*trcs_VLN(ids_H1PO4B,NU(NY,NX),NY,NX)
  RFLC2B=VFLW*AZMAX1(trcsa_soHml2(idsa_CaH2PO4B,NU(NY,NX),NY,NX))*trcs_VLN(ids_H1PO4B,NU(NY,NX),NY,NX)
  RFLM1B=VFLW*AZMAX1(trcsa_soHml2(idsa_MgHPO4B,NU(NY,NX),NY,NX))*trcs_VLN(ids_H1PO4B,NU(NY,NX),NY,NX)

  trcsa_RFXS(idsa_Al,NU(NY,NX),NY,NX)=RFLAL
  trcsa_RFXS(idsa_Fe,NU(NY,NX),NY,NX)=RFLFE
  trcsa_RFXS(idsa_Hp,NU(NY,NX),NY,NX)=RFLHY
  trcsa_RFXS(idsa_Ca,NU(NY,NX),NY,NX)=RFLCA
  trcsa_RFXS(idsa_Mg,NU(NY,NX),NY,NX)=RFLMG
  trcsa_RFXS(idsa_Na,NU(NY,NX),NY,NX)=RFLNA
  trcsa_RFXS(idsa_K,NU(NY,NX),NY,NX)=RFLKA
  trcsa_RFXS(idsa_OH,NU(NY,NX),NY,NX)=RFLOH
  trcsa_RFXS(idsa_SO4,NU(NY,NX),NY,NX)=RFLSO
  trcsa_RFXS(idsa_Cl,NU(NY,NX),NY,NX)=RFLCL
  trcsa_RFXS(idsa_CO3,NU(NY,NX),NY,NX)=RFLC3
  trcsa_RFXS(idsa_HCO3,NU(NY,NX),NY,NX)=RFLHC
  trcsa_RFXS(idsa_AlOH,NU(NY,NX),NY,NX)=RFLAL1
  trcsa_RFXS(idsa_AlOH2,NU(NY,NX),NY,NX)=RFLAL2
  trcsa_RFXS(idsa_AlOH3,NU(NY,NX),NY,NX)=RFLAL3
  trcsa_RFXS(idsa_AlOH4,NU(NY,NX),NY,NX)=RFLAL4
  trcsa_RFXS(idsa_AlSO4,NU(NY,NX),NY,NX)=RFLALS
  trcsa_RFXS(idsa_FeOH,NU(NY,NX),NY,NX)=RFLFE1
  trcsa_RFXS(idsa_FeOH2,NU(NY,NX),NY,NX)=RFLFE2
  trcsa_RFXS(idsa_FeOH3,NU(NY,NX),NY,NX)=RFLFE3
  trcsa_RFXS(idsa_FeOH4,NU(NY,NX),NY,NX)=RFLFE4
  trcsa_RFXS(idsa_FeSO4,NU(NY,NX),NY,NX)=RFLFES
  trcsa_RFXS(idsa_CaOH2,NU(NY,NX),NY,NX)=RFLCAO
  trcsa_RFXS(idsa_CaCO3,NU(NY,NX),NY,NX)=RFLCAC
  trcsa_RFXS(idsa_CaHCO3,NU(NY,NX),NY,NX)=RFLCAH
  trcsa_RFXS(idsa_CaSO4,NU(NY,NX),NY,NX)=RFLCAS
  trcsa_RFXS(idsa_MgOH2,NU(NY,NX),NY,NX)=RFLMGO
  trcsa_RFXS(idsa_MgCO3,NU(NY,NX),NY,NX)=RFLMGC
  trcsa_RFXS(idsa_MgHCO3,NU(NY,NX),NY,NX)=RFLMGH
  trcsa_RFXS(idsa_MgSO4,NU(NY,NX),NY,NX)=RFLMGS
  trcsa_RFXS(idsa_NaCO3,NU(NY,NX),NY,NX)=RFLNAC
  trcsa_RFXS(idsa_NaSO4,NU(NY,NX),NY,NX)=RFLNAS
  trcsa_RFXS(idsa_KSO4,NU(NY,NX),NY,NX)=RFLKAS
  trcsa_RFXS(idsa_H0PO4,NU(NY,NX),NY,NX)=RFLH0P
  trcsa_RFXS(idsa_H3PO4,NU(NY,NX),NY,NX)=RFLH3P
  trcsa_RFXS(idsa_FeHPO4,NU(NY,NX),NY,NX)=RFLF1P
  trcsa_RFXS(idsa_FeH2PO4,NU(NY,NX),NY,NX)=RFLF2P
  trcsa_RFXS(idsa_CaPO4,NU(NY,NX),NY,NX)=RFLC0P
  trcsa_RFXS(idsa_CaHPO4,NU(NY,NX),NY,NX)=RFLC1P
  trcsa_RFXS(idsa_CaH2PO4,NU(NY,NX),NY,NX)=RFLC2P
  trcsa_RFXS(idsa_MgHPO4,NU(NY,NX),NY,NX)=RFLM1P
  trcsa_RFXS(idsa_H0PO4B,NU(NY,NX),NY,NX)=RFLH0B
  trcsa_RFXS(idsa_H3PO4B,NU(NY,NX),NY,NX)=RFLH3B
  trcsa_RFXS(idsa_FeHPO4B,NU(NY,NX),NY,NX)=RFLF1B
  trcsa_RFXS(idsa_FeH2PO4B,NU(NY,NX),NY,NX)=RFLF2B
  trcsa_RFXS(idsa_CaPO4B,NU(NY,NX),NY,NX)=RFLC0B
  trcsa_RFXS(idsa_CaHPO4B,NU(NY,NX),NY,NX)=RFLC1B
  trcsa_RFXS(idsa_CaH2PO4B,NU(NY,NX),NY,NX)=RFLC2B
  trcsa_RFXS(idsa_MgHPO4B,NU(NY,NX),NY,NX)=RFLM1B
  end subroutine MacToMicPoreSoluteAdvExchange
!------------------------------------------------------------------------------------------

  subroutine MicToMacPoreSoluteAdvExchange(M,NY,NX)
!
!     Description:
!
  implicit none
  integer, intent(in) :: M,NY,NX
  real(r8) :: VFLW
!     begin_execution
  IF(VOLWM(M,NU(NY,NX),NY,NX).GT.ZEROS2(NY,NX))THEN
    VFLW=AZMIN1(AMAX1(-VFLWX,FINHM(M,NU(NY,NX),NY,NX)/VOLWM(M,NU(NY,NX),NY,NX)))
  ELSE
    VFLW=-VFLWX
  ENDIF
  RFLAL=VFLW*AZMAX1(trcsa_solml2(idsa_Al,NU(NY,NX),NY,NX))
  RFLFE=VFLW*AZMAX1(trcsa_solml2(idsa_Fe,NU(NY,NX),NY,NX))
  RFLHY=VFLW*AZMAX1(trcsa_solml2(idsa_Hp,NU(NY,NX),NY,NX))
  RFLCA=VFLW*AZMAX1(trcsa_solml2(idsa_Ca,NU(NY,NX),NY,NX))
  RFLMG=VFLW*AZMAX1(trcsa_solml2(idsa_Mg,NU(NY,NX),NY,NX))
  RFLNA=VFLW*AZMAX1(trcsa_solml2(idsa_Na,NU(NY,NX),NY,NX))
  RFLKA=VFLW*AZMAX1(trcsa_solml2(idsa_K,NU(NY,NX),NY,NX))
  RFLOH=VFLW*AZMAX1(trcsa_solml2(idsa_OH,NU(NY,NX),NY,NX))
  RFLSO=VFLW*AZMAX1(trcsa_solml2(idsa_SO4,NU(NY,NX),NY,NX))
  RFLCL=VFLW*AZMAX1(trcsa_solml2(idsa_Cl,NU(NY,NX),NY,NX))
  RFLC3=VFLW*AZMAX1(trcsa_solml2(idsa_CO3,NU(NY,NX),NY,NX))
  RFLHC=VFLW*AZMAX1(trcsa_solml2(idsa_HCO3,NU(NY,NX),NY,NX))
  RFLAL1=VFLW*AZMAX1(trcsa_solml2(idsa_AlOH,NU(NY,NX),NY,NX))
  RFLAL2=VFLW*AZMAX1(trcsa_solml2(idsa_AlOH2,NU(NY,NX),NY,NX))
  RFLAL3=VFLW*AZMAX1(trcsa_solml2(idsa_AlOH3,NU(NY,NX),NY,NX))
  RFLAL4=VFLW*AZMAX1(trcsa_solml2(idsa_AlOH4,NU(NY,NX),NY,NX))
  RFLALS=VFLW*AZMAX1(trcsa_solml2(idsa_AlSO4,NU(NY,NX),NY,NX))
  RFLFE1=VFLW*AZMAX1(trcsa_solml2(idsa_FeOH,NU(NY,NX),NY,NX))
  RFLFE2=VFLW*AZMAX1(trcsa_solml2(idsa_FeOH2,NU(NY,NX),NY,NX))
  RFLFE3=VFLW*AZMAX1(trcsa_solml2(idsa_FeOH3,NU(NY,NX),NY,NX))
  RFLFE4=VFLW*AZMAX1(trcsa_solml2(idsa_FeOH4,NU(NY,NX),NY,NX))
  RFLFES=VFLW*AZMAX1(trcsa_solml2(idsa_FeSO4,NU(NY,NX),NY,NX))
  RFLCAO=VFLW*AZMAX1(trcsa_solml2(idsa_CaOH2,NU(NY,NX),NY,NX))
  RFLCAC=VFLW*AZMAX1(trcsa_solml2(idsa_CaCO3,NU(NY,NX),NY,NX))
  RFLCAH=VFLW*AZMAX1(trcsa_solml2(idsa_CaHCO3,NU(NY,NX),NY,NX))
  RFLCAS=VFLW*AZMAX1(trcsa_solml2(idsa_CaSO4,NU(NY,NX),NY,NX))
  RFLMGO=VFLW*AZMAX1(trcsa_solml2(idsa_MgOH2,NU(NY,NX),NY,NX))
  RFLMGC=VFLW*AZMAX1(trcsa_solml2(idsa_MgCO3,NU(NY,NX),NY,NX))
  RFLMGH=VFLW*AZMAX1(trcsa_solml2(idsa_MgHCO3,NU(NY,NX),NY,NX))
  RFLMGS=VFLW*AZMAX1(trcsa_solml2(idsa_MgSO4,NU(NY,NX),NY,NX))
  RFLNAC=VFLW*AZMAX1(trcsa_solml2(idsa_NaCO3,NU(NY,NX),NY,NX))
  RFLNAS=VFLW*AZMAX1(trcsa_solml2(idsa_NaSO4,NU(NY,NX),NY,NX))
  RFLKAS=VFLW*AZMAX1(trcsa_solml2(idsa_KSO4,NU(NY,NX),NY,NX))
  RFLH0P=VFLW*AZMAX1(trcsa_solml2(idsa_H0PO4,NU(NY,NX),NY,NX))*trcs_VLN(ids_H1PO4,NU(NY,NX),NY,NX)
  RFLH3P=VFLW*AZMAX1(trcsa_solml2(idsa_H3PO4,NU(NY,NX),NY,NX))*trcs_VLN(ids_H1PO4,NU(NY,NX),NY,NX)
  RFLF1P=VFLW*AZMAX1(trcsa_solml2(idsa_FeHPO4,NU(NY,NX),NY,NX))*trcs_VLN(ids_H1PO4,NU(NY,NX),NY,NX)
  RFLF2P=VFLW*AZMAX1(trcsa_solml2(idsa_FeH2PO4,NU(NY,NX),NY,NX))*trcs_VLN(ids_H1PO4,NU(NY,NX),NY,NX)
  RFLC0P=VFLW*AZMAX1(trcsa_solml2(idsa_CaPO4,NU(NY,NX),NY,NX))*trcs_VLN(ids_H1PO4,NU(NY,NX),NY,NX)
  RFLC1P=VFLW*AZMAX1(trcsa_solml2(idsa_CaHPO4,NU(NY,NX),NY,NX))*trcs_VLN(ids_H1PO4,NU(NY,NX),NY,NX)
  RFLC2P=VFLW*AZMAX1(trcsa_solml2(idsa_CaH2PO4,NU(NY,NX),NY,NX))*trcs_VLN(ids_H1PO4,NU(NY,NX),NY,NX)
  RFLM1P=VFLW*AZMAX1(trcsa_solml2(idsa_MgHPO4,NU(NY,NX),NY,NX))*trcs_VLN(ids_H1PO4,NU(NY,NX),NY,NX)
  RFLH0B=VFLW*AZMAX1(trcsa_solml2(idsa_H0PO4B,NU(NY,NX),NY,NX))*trcs_VLN(ids_H1PO4B,NU(NY,NX),NY,NX)
  RFLH3B=VFLW*AZMAX1(trcsa_solml2(idsa_H3PO4B,NU(NY,NX),NY,NX))*trcs_VLN(ids_H1PO4B,NU(NY,NX),NY,NX)
  RFLF1B=VFLW*AZMAX1(trcsa_solml2(idsa_FeHPO4B,NU(NY,NX),NY,NX))*trcs_VLN(ids_H1PO4B,NU(NY,NX),NY,NX)
  RFLF2B=VFLW*AZMAX1(trcsa_solml2(idsa_FeH2PO4B,NU(NY,NX),NY,NX))*trcs_VLN(ids_H1PO4B,NU(NY,NX),NY,NX)
  RFLC0B=VFLW*AZMAX1(trcsa_solml2(idsa_CaPO4B,NU(NY,NX),NY,NX))*trcs_VLN(ids_H1PO4B,NU(NY,NX),NY,NX)
  RFLC1B=VFLW*AZMAX1(trcsa_solml2(idsa_CaHPO4B,NU(NY,NX),NY,NX))*trcs_VLN(ids_H1PO4B,NU(NY,NX),NY,NX)
  RFLC2B=VFLW*AZMAX1(trcsa_solml2(idsa_CaH2PO4B,NU(NY,NX),NY,NX))*trcs_VLN(ids_H1PO4B,NU(NY,NX),NY,NX)
  RFLM1B=VFLW*AZMAX1(trcsa_solml2(idsa_MgHPO4B,NU(NY,NX),NY,NX))*trcs_VLN(ids_H1PO4B,NU(NY,NX),NY,NX)


  trcsa_RFXS(idsa_Al,NU(NY,NX),NY,NX)=RFLAL
  trcsa_RFXS(idsa_Fe,NU(NY,NX),NY,NX)=RFLFE
  trcsa_RFXS(idsa_Hp,NU(NY,NX),NY,NX)=RFLHY
  trcsa_RFXS(idsa_Ca,NU(NY,NX),NY,NX)=RFLCA
  trcsa_RFXS(idsa_Mg,NU(NY,NX),NY,NX)=RFLMG
  trcsa_RFXS(idsa_Na,NU(NY,NX),NY,NX)=RFLNA
  trcsa_RFXS(idsa_K,NU(NY,NX),NY,NX)=RFLKA
  trcsa_RFXS(idsa_OH,NU(NY,NX),NY,NX)=RFLOH
  trcsa_RFXS(idsa_SO4,NU(NY,NX),NY,NX)=RFLSO
  trcsa_RFXS(idsa_Cl,NU(NY,NX),NY,NX)=RFLCL
  trcsa_RFXS(idsa_CO3,NU(NY,NX),NY,NX)=RFLC3
  trcsa_RFXS(idsa_HCO3,NU(NY,NX),NY,NX)=RFLHC
  trcsa_RFXS(idsa_AlOH,NU(NY,NX),NY,NX)=RFLAL1
  trcsa_RFXS(idsa_AlOH2,NU(NY,NX),NY,NX)=RFLAL2
  trcsa_RFXS(idsa_AlOH3,NU(NY,NX),NY,NX)=RFLAL3
  trcsa_RFXS(idsa_AlOH4,NU(NY,NX),NY,NX)=RFLAL4
  trcsa_RFXS(idsa_AlSO4,NU(NY,NX),NY,NX)=RFLALS
  trcsa_RFXS(idsa_FeOH,NU(NY,NX),NY,NX)=RFLFE1
  trcsa_RFXS(idsa_FeOH2,NU(NY,NX),NY,NX)=RFLFE2
  trcsa_RFXS(idsa_FeOH3,NU(NY,NX),NY,NX)=RFLFE3
  trcsa_RFXS(idsa_FeOH4,NU(NY,NX),NY,NX)=RFLFE4
  trcsa_RFXS(idsa_FeSO4,NU(NY,NX),NY,NX)=RFLFES
  trcsa_RFXS(idsa_CaOH2,NU(NY,NX),NY,NX)=RFLCAO
  trcsa_RFXS(idsa_CaCO3,NU(NY,NX),NY,NX)=RFLCAC
  trcsa_RFXS(idsa_CaHCO3,NU(NY,NX),NY,NX)=RFLCAH
  trcsa_RFXS(idsa_CaSO4,NU(NY,NX),NY,NX)=RFLCAS
  trcsa_RFXS(idsa_MgOH2,NU(NY,NX),NY,NX)=RFLMGO
  trcsa_RFXS(idsa_MgCO3,NU(NY,NX),NY,NX)=RFLMGC
  trcsa_RFXS(idsa_MgHCO3,NU(NY,NX),NY,NX)=RFLMGH
  trcsa_RFXS(idsa_MgSO4,NU(NY,NX),NY,NX)=RFLMGS
  trcsa_RFXS(idsa_NaCO3,NU(NY,NX),NY,NX)=RFLNAC
  trcsa_RFXS(idsa_NaSO4,NU(NY,NX),NY,NX)=RFLNAS
  trcsa_RFXS(idsa_KSO4,NU(NY,NX),NY,NX)=RFLKAS
  trcsa_RFXS(idsa_H0PO4,NU(NY,NX),NY,NX)=RFLH0P
  trcsa_RFXS(idsa_H3PO4,NU(NY,NX),NY,NX)=RFLH3P
  trcsa_RFXS(idsa_FeHPO4,NU(NY,NX),NY,NX)=RFLF1P
  trcsa_RFXS(idsa_FeH2PO4,NU(NY,NX),NY,NX)=RFLF2P
  trcsa_RFXS(idsa_CaPO4,NU(NY,NX),NY,NX)=RFLC0P
  trcsa_RFXS(idsa_CaHPO4,NU(NY,NX),NY,NX)=RFLC1P
  trcsa_RFXS(idsa_CaH2PO4,NU(NY,NX),NY,NX)=RFLC2P
  trcsa_RFXS(idsa_MgHPO4,NU(NY,NX),NY,NX)=RFLM1P
  trcsa_RFXS(idsa_H0PO4B,NU(NY,NX),NY,NX)=RFLH0B
  trcsa_RFXS(idsa_H3PO4B,NU(NY,NX),NY,NX)=RFLH3B
  trcsa_RFXS(idsa_FeHPO4B,NU(NY,NX),NY,NX)=RFLF1B
  trcsa_RFXS(idsa_FeH2PO4B,NU(NY,NX),NY,NX)=RFLF2B
  trcsa_RFXS(idsa_CaPO4B,NU(NY,NX),NY,NX)=RFLC0B
  trcsa_RFXS(idsa_CaHPO4B,NU(NY,NX),NY,NX)=RFLC1B
  trcsa_RFXS(idsa_CaH2PO4B,NU(NY,NX),NY,NX)=RFLC2B
  trcsa_RFXS(idsa_MgHPO4B,NU(NY,NX),NY,NX)=RFLM1B
  end subroutine MicToMacPoreSoluteAdvExchange
!------------------------------------------------------------------------------------------

  subroutine ZeroMicMacPoreExchange
!
!     Description:
!
  implicit none

  RFLAL=0.0_r8
  RFLFE=0.0_r8
  RFLHY=0.0_r8
  RFLCA=0.0_r8
  RFLMG=0.0_r8
  RFLNA=0.0_r8
  RFLKA=0.0_r8
  RFLOH=0.0_r8
  RFLSO=0.0_r8
  RFLCL=0.0_r8
  RFLC3=0.0_r8
  RFLHC=0.0_r8
  RFLAL1=0.0_r8
  RFLAL2=0.0_r8
  RFLAL3=0.0_r8
  RFLAL4=0.0_r8
  RFLALS=0.0_r8
  RFLFE1=0.0_r8
  RFLFE2=0.0_r8
  RFLFE3=0.0_r8
  RFLFE4=0.0_r8
  RFLFES=0.0_r8
  RFLCAO=0.0_r8
  RFLCAC=0.0_r8
  RFLCAH=0.0_r8
  RFLCAS=0.0_r8
  RFLMGO=0.0_r8
  RFLMGC=0.0_r8
  RFLMGH=0.0_r8
  RFLMGS=0.0_r8
  RFLNAC=0.0_r8
  RFLNAS=0.0_r8
  RFLKAS=0.0_r8
  RFLH0P=0.0_r8
  RFLH3P=0.0_r8
  RFLF1P=0.0_r8
  RFLF2P=0.0_r8
  RFLC0P=0.0_r8
  RFLC1P=0.0_r8
  RFLC2P=0.0_r8
  RFLM1P=0.0_r8
  RFLH0B=0.0_r8
  RFLH3B=0.0_r8
  RFLF1B=0.0_r8
  RFLF2B=0.0_r8
  RFLC0B=0.0_r8
  RFLC1B=0.0_r8
  RFLC2B=0.0_r8
  RFLM1B=0.0_r8
  end subroutine ZeroMicMacPoreExchange
!------------------------------------------------------------------------------------------

  subroutine MacMicPoreSoluteDifusExchange(M,NY,NX)
!
!     Description:
!
  implicit none
  integer, intent(in) :: M,NY,NX
  real(r8) :: VOLWHS,VOLWT,VOLWMNU
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
  DFVAL=XNPH*(AZMAX1(trcsa_soHml2(idsa_Al,NU(NY,NX),NY,NX))*VOLWMNU-AZMAX1(trcsa_solml2(idsa_Al,NU(NY,NX),NY,NX))*VOLWHS)/VOLWT
  DFVFE=XNPH*(AZMAX1(trcsa_soHml2(idsa_Fe,NU(NY,NX),NY,NX))*VOLWMNU-AZMAX1(trcsa_solml2(idsa_Fe,NU(NY,NX),NY,NX))*VOLWHS)/VOLWT
  DFVHY=XNPH*(AZMAX1(trcsa_soHml2(idsa_Hp,NU(NY,NX),NY,NX))*VOLWMNU-AZMAX1(trcsa_solml2(idsa_Hp,NU(NY,NX),NY,NX))*VOLWHS)/VOLWT
  DFVCA=XNPH*(AZMAX1(trcsa_soHml2(idsa_Ca,NU(NY,NX),NY,NX))*VOLWMNU-AZMAX1(trcsa_solml2(idsa_Ca,NU(NY,NX),NY,NX))*VOLWHS)/VOLWT
  DFVMG=XNPH*(AZMAX1(trcsa_soHml2(idsa_Mg,NU(NY,NX),NY,NX))*VOLWMNU-AZMAX1(trcsa_solml2(idsa_Mg,NU(NY,NX),NY,NX))*VOLWHS)/VOLWT
  DFVNA=XNPH*(AZMAX1(trcsa_soHml2(idsa_Na,NU(NY,NX),NY,NX))*VOLWMNU-AZMAX1(trcsa_solml2(idsa_Na,NU(NY,NX),NY,NX))*VOLWHS)/VOLWT
  DFVKA=XNPH*(AZMAX1(trcsa_soHml2(idsa_K,NU(NY,NX),NY,NX))*VOLWMNU-AZMAX1(trcsa_solml2(idsa_K,NU(NY,NX),NY,NX))*VOLWHS)/VOLWT
  DFVOH=XNPH*(AZMAX1(trcsa_soHml2(idsa_OH,NU(NY,NX),NY,NX))*VOLWMNU-AZMAX1(trcsa_solml2(idsa_OH,NU(NY,NX),NY,NX))*VOLWHS)/VOLWT
  DFVSO=XNPH*(AZMAX1(trcsa_soHml2(idsa_SO4,NU(NY,NX),NY,NX))*VOLWMNU-AZMAX1(trcsa_solml2(idsa_SO4,NU(NY,NX),NY,NX))*VOLWHS)/VOLWT
  DFVCL=XNPH*(AZMAX1(trcsa_soHml2(idsa_Cl,NU(NY,NX),NY,NX))*VOLWMNU-AZMAX1(trcsa_solml2(idsa_Cl,NU(NY,NX),NY,NX))*VOLWHS)/VOLWT
  DFVC3=XNPH*(AZMAX1(trcsa_soHml2(idsa_CO3,NU(NY,NX),NY,NX))*VOLWMNU-AZMAX1(trcsa_solml2(idsa_CO3,NU(NY,NX),NY,NX))*VOLWHS)/VOLWT
  DFVHC=XNPH*(AZMAX1(trcsa_soHml2(idsa_HCO3,NU(NY,NX),NY,NX))*VOLWMNU-AZMAX1(trcsa_solml2(idsa_HCO3,NU(NY,NX),NY,NX))*VOLWHS)/VOLWT
  DFVAL1=XNPH*(AZMAX1(trcsa_soHml2(idsa_AlOH,NU(NY,NX),NY,NX))*VOLWMNU-AZMAX1(trcsa_solml2(idsa_AlOH,NU(NY,NX),NY,NX))*VOLWHS)/VOLWT
  DFVAL2=XNPH*(AZMAX1(trcsa_soHml2(idsa_AlOH2,NU(NY,NX),NY,NX))*VOLWMNU-AZMAX1(trcsa_solml2(idsa_AlOH2,NU(NY,NX),NY,NX))*VOLWHS)/VOLWT
  DFVAL3=XNPH*(AZMAX1(trcsa_soHml2(idsa_AlOH3,NU(NY,NX),NY,NX))*VOLWMNU-AZMAX1(trcsa_solml2(idsa_AlOH3,NU(NY,NX),NY,NX))*VOLWHS)/VOLWT
  DFVAL4=XNPH*(AZMAX1(trcsa_soHml2(idsa_AlOH4,NU(NY,NX),NY,NX))*VOLWMNU-AZMAX1(trcsa_solml2(idsa_AlOH4,NU(NY,NX),NY,NX))*VOLWHS)/VOLWT
  DFVALS=XNPH*(AZMAX1(trcsa_soHml2(idsa_AlSO4,NU(NY,NX),NY,NX))*VOLWMNU-AZMAX1(trcsa_solml2(idsa_AlSO4,NU(NY,NX),NY,NX))*VOLWHS)/VOLWT
  DFVFE1=XNPH*(AZMAX1(trcsa_soHml2(idsa_FeOH,NU(NY,NX),NY,NX))*VOLWMNU-AZMAX1(trcsa_solml2(idsa_FeOH,NU(NY,NX),NY,NX))*VOLWHS)/VOLWT
  DFVF22=XNPH*(AZMAX1(trcsa_soHml2(idsa_FeOH2,NU(NY,NX),NY,NX))*VOLWMNU-AZMAX1(trcsa_solml2(idsa_FeOH2,NU(NY,NX),NY,NX))*VOLWHS)/VOLWT
  DFVFE3=XNPH*(AZMAX1(trcsa_soHml2(idsa_FeOH3,NU(NY,NX),NY,NX))*VOLWMNU-AZMAX1(trcsa_solml2(idsa_FeOH3,NU(NY,NX),NY,NX))*VOLWHS)/VOLWT
  DFVFE4=XNPH*(AZMAX1(trcsa_soHml2(idsa_FeOH4,NU(NY,NX),NY,NX))*VOLWMNU-AZMAX1(trcsa_solml2(idsa_FeOH4,NU(NY,NX),NY,NX))*VOLWHS)/VOLWT
  DFVFES=XNPH*(AZMAX1(trcsa_soHml2(idsa_FeSO4,NU(NY,NX),NY,NX))*VOLWMNU-AZMAX1(trcsa_solml2(idsa_FeSO4,NU(NY,NX),NY,NX))*VOLWHS)/VOLWT
  DFVCAO=XNPH*(AZMAX1(trcsa_soHml2(idsa_CaOH2,NU(NY,NX),NY,NX))*VOLWMNU-AZMAX1(trcsa_solml2(idsa_CaOH2,NU(NY,NX),NY,NX))*VOLWHS)/VOLWT
  DFVCAC=XNPH*(AZMAX1(trcsa_soHml2(idsa_CaCO3,NU(NY,NX),NY,NX))*VOLWMNU-AZMAX1(trcsa_solml2(idsa_CaCO3,NU(NY,NX),NY,NX))*VOLWHS)/VOLWT
  DFVCAH=XNPH*(AZMAX1(trcsa_soHml2(idsa_CaHCO3,NU(NY,NX),NY,NX))*VOLWMNU-AZMAX1(trcsa_solml2(idsa_CaHCO3,NU(NY,NX),NY,NX))*VOLWHS)/VOLWT
  DFVCAS=XNPH*(AZMAX1(trcsa_soHml2(idsa_CaSO4,NU(NY,NX),NY,NX))*VOLWMNU-AZMAX1(trcsa_solml2(idsa_CaSO4,NU(NY,NX),NY,NX))*VOLWHS)/VOLWT
  DFVMGO=XNPH*(AZMAX1(trcsa_soHml2(idsa_MgOH2,NU(NY,NX),NY,NX))*VOLWMNU-AZMAX1(trcsa_solml2(idsa_MgOH2,NU(NY,NX),NY,NX))*VOLWHS)/VOLWT
  DFVMGC=XNPH*(AZMAX1(trcsa_soHml2(idsa_MgCO3,NU(NY,NX),NY,NX))*VOLWMNU-AZMAX1(trcsa_solml2(idsa_MgCO3,NU(NY,NX),NY,NX))*VOLWHS)/VOLWT
  DFVMGH=XNPH*(AZMAX1(trcsa_soHml2(idsa_MgHCO3,NU(NY,NX),NY,NX))*VOLWMNU-AZMAX1(trcsa_solml2(idsa_MgHCO3,NU(NY,NX),NY,NX))*VOLWHS)/VOLWT
  DFVMGS=XNPH*(AZMAX1(trcsa_soHml2(idsa_MgSO4,NU(NY,NX),NY,NX))*VOLWMNU-AZMAX1(trcsa_solml2(idsa_MgSO4,NU(NY,NX),NY,NX))*VOLWHS)/VOLWT
  DFVNAC=XNPH*(AZMAX1(trcsa_soHml2(idsa_NaCO3,NU(NY,NX),NY,NX))*VOLWMNU-AZMAX1(trcsa_solml2(idsa_NaCO3,NU(NY,NX),NY,NX))*VOLWHS)/VOLWT
  DFVNAS=XNPH*(AZMAX1(trcsa_soHml2(idsa_NaSO4,NU(NY,NX),NY,NX))*VOLWMNU-AZMAX1(trcsa_solml2(idsa_NaSO4,NU(NY,NX),NY,NX))*VOLWHS)/VOLWT
  DFVKAS=XNPH*(AZMAX1(trcsa_soHml2(idsa_KSO4,NU(NY,NX),NY,NX))*VOLWMNU-AZMAX1(trcsa_solml2(idsa_KSO4,NU(NY,NX),NY,NX))*VOLWHS)/VOLWT
  DFVNAC=XNPH*(AZMAX1(trcsa_soHml2(idsa_NaCO3,NU(NY,NX),NY,NX))*VOLWMNU-AZMAX1(trcsa_solml2(idsa_NaCO3,NU(NY,NX),NY,NX))*VOLWHS)/VOLWT
  DFVH0P=XNPH*(AZMAX1(trcsa_soHml2(idsa_H0PO4,NU(NY,NX),NY,NX))*VOLWMNU-AZMAX1(trcsa_solml2(idsa_H0PO4,NU(NY,NX),NY,NX))*VOLWHS)/VOLWT*trcs_VLN(ids_H1PO4,NU(NY,NX),NY,NX)
  DFVH3P=XNPH*(AZMAX1(trcsa_soHml2(idsa_H3PO4,NU(NY,NX),NY,NX))*VOLWMNU-AZMAX1(trcsa_solml2(idsa_H3PO4,NU(NY,NX),NY,NX))*VOLWHS)/VOLWT*trcs_VLN(ids_H1PO4,NU(NY,NX),NY,NX)
  DFVF1P=XNPH*(AZMAX1(trcsa_soHml2(idsa_FeHPO4,NU(NY,NX),NY,NX))*VOLWMNU-AZMAX1(trcsa_solml2(idsa_FeHPO4,NU(NY,NX),NY,NX))*VOLWHS)/VOLWT*trcs_VLN(ids_H1PO4,NU(NY,NX),NY,NX)
  DFVF2P=XNPH*(AZMAX1(trcsa_soHml2(idsa_FeH2PO4,NU(NY,NX),NY,NX))*VOLWMNU-AZMAX1(trcsa_solml2(idsa_FeH2PO4,NU(NY,NX),NY,NX))*VOLWHS)/VOLWT*trcs_VLN(ids_H1PO4,NU(NY,NX),NY,NX)
  DFVC0P=XNPH*(AZMAX1(trcsa_soHml2(idsa_CaPO4,NU(NY,NX),NY,NX))*VOLWMNU-AZMAX1(trcsa_solml2(idsa_CaPO4,NU(NY,NX),NY,NX))*VOLWHS)/VOLWT*trcs_VLN(ids_H1PO4,NU(NY,NX),NY,NX)
  DFVC1P=XNPH*(AZMAX1(trcsa_soHml2(idsa_CaHPO4,NU(NY,NX),NY,NX))*VOLWMNU-AZMAX1(trcsa_solml2(idsa_CaHPO4,NU(NY,NX),NY,NX))*VOLWHS)/VOLWT*trcs_VLN(ids_H1PO4,NU(NY,NX),NY,NX)
  DFVC2P=XNPH*(AZMAX1(trcsa_soHml2(idsa_CaH2PO4,NU(NY,NX),NY,NX))*VOLWMNU-AZMAX1(trcsa_solml2(idsa_CaH2PO4,NU(NY,NX),NY,NX))*VOLWHS)/VOLWT*trcs_VLN(ids_H1PO4,NU(NY,NX),NY,NX)
  DFVM1P=XNPH*(AZMAX1(trcsa_soHml2(idsa_MgHPO4,NU(NY,NX),NY,NX))*VOLWMNU-AZMAX1(trcsa_solml2(idsa_MgHPO4,NU(NY,NX),NY,NX))*VOLWHS)/VOLWT*trcs_VLN(ids_H1PO4,NU(NY,NX),NY,NX)
  DFVH0B=XNPH*(AZMAX1(trcsa_soHml2(idsa_H0PO4B,NU(NY,NX),NY,NX))*VOLWMNU-AZMAX1(trcsa_solml2(idsa_H0PO4B,NU(NY,NX),NY,NX))*VOLWHS)/VOLWT*trcs_VLN(ids_H1PO4B,NU(NY,NX),NY,NX)
  DFVH3B=XNPH*(AZMAX1(trcsa_soHml2(idsa_H3PO4B,NU(NY,NX),NY,NX))*VOLWMNU-AZMAX1(trcsa_solml2(idsa_H3PO4B,NU(NY,NX),NY,NX))*VOLWHS)/VOLWT*trcs_VLN(ids_H1PO4B,NU(NY,NX),NY,NX)
  DFVF1B=XNPH*(AZMAX1(trcsa_soHml2(idsa_FeHPO4B,NU(NY,NX),NY,NX))*VOLWMNU-AZMAX1(trcsa_solml2(idsa_FeHPO4B,NU(NY,NX),NY,NX))*VOLWHS)/VOLWT*trcs_VLN(ids_H1PO4B,NU(NY,NX),NY,NX)
  DFVF2B=XNPH*(AZMAX1(trcsa_soHml2(idsa_FeH2PO4B,NU(NY,NX),NY,NX))*VOLWMNU-AZMAX1(trcsa_solml2(idsa_FeH2PO4B,NU(NY,NX),NY,NX))*VOLWHS)/VOLWT*trcs_VLN(ids_H1PO4B,NU(NY,NX),NY,NX)
  DFVC0B=XNPH*(AZMAX1(trcsa_soHml2(idsa_CaPO4B,NU(NY,NX),NY,NX))*VOLWMNU-AZMAX1(trcsa_solml2(idsa_CaPO4B,NU(NY,NX),NY,NX))*VOLWHS)/VOLWT*trcs_VLN(ids_H1PO4B,NU(NY,NX),NY,NX)
  DFVC1B=XNPH*(AZMAX1(trcsa_soHml2(idsa_CaHPO4B,NU(NY,NX),NY,NX))*VOLWMNU-AZMAX1(trcsa_solml2(idsa_CaHPO4B,NU(NY,NX),NY,NX))*VOLWHS)/VOLWT*trcs_VLN(ids_H1PO4B,NU(NY,NX),NY,NX)
  DFVC2B=XNPH*(AZMAX1(trcsa_soHml2(idsa_CaH2PO4B,NU(NY,NX),NY,NX))*VOLWMNU-AZMAX1(trcsa_solml2(idsa_CaH2PO4B,NU(NY,NX),NY,NX))*VOLWHS)/VOLWT*trcs_VLN(ids_H1PO4B,NU(NY,NX),NY,NX)
  DFVM1B=XNPH*(AZMAX1(trcsa_soHml2(idsa_MgHPO4B,NU(NY,NX),NY,NX))*VOLWMNU-AZMAX1(trcsa_solml2(idsa_MgHPO4B,NU(NY,NX),NY,NX))*VOLWHS)/VOLWT*trcs_VLN(ids_H1PO4B,NU(NY,NX),NY,NX)


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
  trcsa_RFXS(idsa_Al,NU(NY,NX),NY,NX)=trcsa_RFXS(idsa_Al,NU(NY,NX),NY,NX)+DFVAL
  trcsa_RFXS(idsa_Fe,NU(NY,NX),NY,NX)=trcsa_RFXS(idsa_Fe,NU(NY,NX),NY,NX)+DFVFE
  trcsa_RFXS(idsa_Hp,NU(NY,NX),NY,NX)=trcsa_RFXS(idsa_Hp,NU(NY,NX),NY,NX)+DFVHY
  trcsa_RFXS(idsa_Ca,NU(NY,NX),NY,NX)=trcsa_RFXS(idsa_Ca,NU(NY,NX),NY,NX)+DFVCA
  trcsa_RFXS(idsa_Mg,NU(NY,NX),NY,NX)=trcsa_RFXS(idsa_Mg,NU(NY,NX),NY,NX)+DFVMG
  trcsa_RFXS(idsa_Na,NU(NY,NX),NY,NX)=trcsa_RFXS(idsa_Na,NU(NY,NX),NY,NX)+DFVNA
  trcsa_RFXS(idsa_K,NU(NY,NX),NY,NX)=trcsa_RFXS(idsa_K,NU(NY,NX),NY,NX)+DFVKA
  trcsa_RFXS(idsa_OH,NU(NY,NX),NY,NX)=trcsa_RFXS(idsa_OH,NU(NY,NX),NY,NX)+DFVOH
  trcsa_RFXS(idsa_SO4,NU(NY,NX),NY,NX)=trcsa_RFXS(idsa_SO4,NU(NY,NX),NY,NX)+DFVSO
  trcsa_RFXS(idsa_Cl,NU(NY,NX),NY,NX)=trcsa_RFXS(idsa_Cl,NU(NY,NX),NY,NX)+DFVCL
  trcsa_RFXS(idsa_CO3,NU(NY,NX),NY,NX)=trcsa_RFXS(idsa_CO3,NU(NY,NX),NY,NX)+DFVC3
  trcsa_RFXS(idsa_HCO3,NU(NY,NX),NY,NX)=trcsa_RFXS(idsa_HCO3,NU(NY,NX),NY,NX)+DFVHC
  trcsa_RFXS(idsa_AlOH,NU(NY,NX),NY,NX)=trcsa_RFXS(idsa_AlOH,NU(NY,NX),NY,NX)+DFVAL1
  trcsa_RFXS(idsa_AlOH2,NU(NY,NX),NY,NX)=trcsa_RFXS(idsa_AlOH2,NU(NY,NX),NY,NX)+DFVAL2
  trcsa_RFXS(idsa_AlOH3,NU(NY,NX),NY,NX)=trcsa_RFXS(idsa_AlOH3,NU(NY,NX),NY,NX)+DFVAL3
  trcsa_RFXS(idsa_AlOH4,NU(NY,NX),NY,NX)=trcsa_RFXS(idsa_AlOH4,NU(NY,NX),NY,NX)+DFVAL4
  trcsa_RFXS(idsa_AlSO4,NU(NY,NX),NY,NX)=trcsa_RFXS(idsa_AlSO4,NU(NY,NX),NY,NX)+DFVALS
  trcsa_RFXS(idsa_FeOH,NU(NY,NX),NY,NX)=trcsa_RFXS(idsa_FeOH,NU(NY,NX),NY,NX)+DFVFE1
  trcsa_RFXS(idsa_FeOH2,NU(NY,NX),NY,NX)=trcsa_RFXS(idsa_FeOH2,NU(NY,NX),NY,NX)+DFVFE2
  trcsa_RFXS(idsa_FeOH3,NU(NY,NX),NY,NX)=trcsa_RFXS(idsa_FeOH3,NU(NY,NX),NY,NX)+DFVFE3
  trcsa_RFXS(idsa_FeOH4,NU(NY,NX),NY,NX)=trcsa_RFXS(idsa_FeOH4,NU(NY,NX),NY,NX)+DFVFE4
  trcsa_RFXS(idsa_FeSO4,NU(NY,NX),NY,NX)=trcsa_RFXS(idsa_FeSO4,NU(NY,NX),NY,NX)+DFVFES
  trcsa_RFXS(idsa_CaOH2,NU(NY,NX),NY,NX)=trcsa_RFXS(idsa_CaOH2,NU(NY,NX),NY,NX)+DFVCAO
  trcsa_RFXS(idsa_CaCO3,NU(NY,NX),NY,NX)=trcsa_RFXS(idsa_CaCO3,NU(NY,NX),NY,NX)+DFVCAC
  trcsa_RFXS(idsa_CaHCO3,NU(NY,NX),NY,NX)=trcsa_RFXS(idsa_CaHCO3,NU(NY,NX),NY,NX)+DFVCAH
  trcsa_RFXS(idsa_CaSO4,NU(NY,NX),NY,NX)=trcsa_RFXS(idsa_CaSO4,NU(NY,NX),NY,NX)+DFVCAS
  trcsa_RFXS(idsa_MgOH2,NU(NY,NX),NY,NX)=trcsa_RFXS(idsa_MgOH2,NU(NY,NX),NY,NX)+DFVMGO
  trcsa_RFXS(idsa_MgCO3,NU(NY,NX),NY,NX)=trcsa_RFXS(idsa_MgCO3,NU(NY,NX),NY,NX)+DFVMGC
  trcsa_RFXS(idsa_MgHCO3,NU(NY,NX),NY,NX)=trcsa_RFXS(idsa_MgHCO3,NU(NY,NX),NY,NX)+DFVMGH
  trcsa_RFXS(idsa_MgSO4,NU(NY,NX),NY,NX)=trcsa_RFXS(idsa_MgSO4,NU(NY,NX),NY,NX)+DFVMGS
  trcsa_RFXS(idsa_NaCO3,NU(NY,NX),NY,NX)=trcsa_RFXS(idsa_NaCO3,NU(NY,NX),NY,NX)+DFVNAC
  trcsa_RFXS(idsa_NaSO4,NU(NY,NX),NY,NX)=trcsa_RFXS(idsa_NaSO4,NU(NY,NX),NY,NX)+DFVNAS
  trcsa_RFXS(idsa_KSO4,NU(NY,NX),NY,NX)=trcsa_RFXS(idsa_KSO4,NU(NY,NX),NY,NX)+DFVKAS
  trcsa_RFXS(idsa_H0PO4,NU(NY,NX),NY,NX)=trcsa_RFXS(idsa_H0PO4,NU(NY,NX),NY,NX)+DFVH0P
  trcsa_RFXS(idsa_H3PO4,NU(NY,NX),NY,NX)=trcsa_RFXS(idsa_H3PO4,NU(NY,NX),NY,NX)+DFVH3P
  trcsa_RFXS(idsa_FeHPO4,NU(NY,NX),NY,NX)=trcsa_RFXS(idsa_FeHPO4,NU(NY,NX),NY,NX)+DFVF1P
  trcsa_RFXS(idsa_FeH2PO4,NU(NY,NX),NY,NX)=trcsa_RFXS(idsa_FeH2PO4,NU(NY,NX),NY,NX)+DFVF2P
  trcsa_RFXS(idsa_CaPO4,NU(NY,NX),NY,NX)=trcsa_RFXS(idsa_CaPO4,NU(NY,NX),NY,NX)+DFVC0P
  trcsa_RFXS(idsa_CaHPO4,NU(NY,NX),NY,NX)=trcsa_RFXS(idsa_CaHPO4,NU(NY,NX),NY,NX)+DFVC1P
  trcsa_RFXS(idsa_CaH2PO4,NU(NY,NX),NY,NX)=trcsa_RFXS(idsa_CaH2PO4,NU(NY,NX),NY,NX)+DFVC2P
  trcsa_RFXS(idsa_MgHPO4,NU(NY,NX),NY,NX)=trcsa_RFXS(idsa_MgHPO4,NU(NY,NX),NY,NX)+DFVM1P
  trcsa_RFXS(idsa_H0PO4B,NU(NY,NX),NY,NX)=trcsa_RFXS(idsa_H0PO4B,NU(NY,NX),NY,NX)+DFVH0B
  trcsa_RFXS(idsa_H3PO4B,NU(NY,NX),NY,NX)=trcsa_RFXS(idsa_H3PO4B,NU(NY,NX),NY,NX)+DFVH3B
  trcsa_RFXS(idsa_FeHPO4B,NU(NY,NX),NY,NX)=trcsa_RFXS(idsa_FeHPO4B,NU(NY,NX),NY,NX)+DFVF1B
  trcsa_RFXS(idsa_FeH2PO4B,NU(NY,NX),NY,NX)=trcsa_RFXS(idsa_FeH2PO4B,NU(NY,NX),NY,NX)+DFVF2B
  trcsa_RFXS(idsa_CaPO4B,NU(NY,NX),NY,NX)=trcsa_RFXS(idsa_CaPO4B,NU(NY,NX),NY,NX)+DFVC0B
  trcsa_RFXS(idsa_CaHPO4B,NU(NY,NX),NY,NX)=trcsa_RFXS(idsa_CaHPO4B,NU(NY,NX),NY,NX)+DFVC1B
  trcsa_RFXS(idsa_CaH2PO4B,NU(NY,NX),NY,NX)=trcsa_RFXS(idsa_CaH2PO4B,NU(NY,NX),NY,NX)+DFVC2B
  trcsa_RFXS(idsa_MgHPO4B,NU(NY,NX),NY,NX)=trcsa_RFXS(idsa_MgHPO4B,NU(NY,NX),NY,NX)+DFVM1B
  end subroutine MacMicPoreSoluteDifusExchange
!------------------------------------------------------------------------------------------

  subroutine ZeroMacMicPoreSoluteDifusExch
!
!     Description:
!
  implicit none
!     begin_execution
  DFVAL=0.0_r8
  DFVFE=0.0_r8
  DFVHY=0.0_r8
  DFVCA=0.0_r8
  DFVMG=0.0_r8
  DFVNA=0.0_r8
  DFVKA=0.0_r8
  DFVOH=0.0_r8
  DFVSO=0.0_r8
  DFVCL=0.0_r8
  DFVC3=0.0_r8
  DFVHC=0.0_r8
  DFVAL1=0.0_r8
  DFVAL2=0.0_r8
  DFVAL3=0.0_r8
  DFVAL4=0.0_r8
  DFVALS=0.0_r8
  DFVFE1=0.0_r8
  DFVFE2=0.0_r8
  DFVFE3=0.0_r8
  DFVFE4=0.0_r8
  DFVFES=0.0_r8
  DFVCAO=0.0_r8
  DFVCAC=0.0_r8
  DFVCAH=0.0_r8
  DFVCAS=0.0_r8
  DFVMGO=0.0_r8
  DFVMGC=0.0_r8
  DFVMGH=0.0_r8
  DFVMGS=0.0_r8
  DFVNAC=0.0_r8
  DFVNAS=0.0_r8
  DFVKAS=0.0_r8
  DFVH0P=0.0_r8
  DFVH3P=0.0_r8
  DFVF1P=0.0_r8
  DFVF2P=0.0_r8
  DFVC0P=0.0_r8
  DFVC1P=0.0_r8
  DFVC2P=0.0_r8
  DFVM1P=0.0_r8
  DFVH0B=0.0_r8
  DFVH3B=0.0_r8
  DFVF1B=0.0_r8
  DFVF2B=0.0_r8
  DFVC0B=0.0_r8
  DFVC1B=0.0_r8
  DFVC2B=0.0_r8
  DFVM1B=0.0_r8
  end subroutine ZeroMacMicPoreSoluteDifusExch
!------------------------------------------------------------------------------------------

  subroutine AccumHourlyMicMacPoreFlux(NY,NX)
!
!     Description:
!
      implicit none
      integer, intent(in) :: NY,NX
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
  trcsa_XFXS(idsa_Al,NU(NY,NX),NY,NX)=trcsa_XFXS(idsa_Al,NU(NY,NX),NY,NX)+trcsa_RFXS(idsa_Al,NU(NY,NX),NY,NX)
  trcsa_XFXS(idsa_Fe,NU(NY,NX),NY,NX)=trcsa_XFXS(idsa_Fe,NU(NY,NX),NY,NX)+trcsa_RFXS(idsa_Fe,NU(NY,NX),NY,NX)
  trcsa_XFXS(idsa_Hp,NU(NY,NX),NY,NX)=trcsa_XFXS(idsa_Hp,NU(NY,NX),NY,NX)+trcsa_RFXS(idsa_Hp,NU(NY,NX),NY,NX)
  trcsa_XFXS(idsa_Ca,NU(NY,NX),NY,NX)=trcsa_XFXS(idsa_Ca,NU(NY,NX),NY,NX)+trcsa_RFXS(idsa_Ca,NU(NY,NX),NY,NX)
  trcsa_XFXS(idsa_Mg,NU(NY,NX),NY,NX)=trcsa_XFXS(idsa_Mg,NU(NY,NX),NY,NX)+trcsa_RFXS(idsa_Mg,NU(NY,NX),NY,NX)
  trcsa_XFXS(idsa_Na,NU(NY,NX),NY,NX)=trcsa_XFXS(idsa_Na,NU(NY,NX),NY,NX)+trcsa_RFXS(idsa_Na,NU(NY,NX),NY,NX)
  trcsa_XFXS(idsa_K,NU(NY,NX),NY,NX)=trcsa_XFXS(idsa_K,NU(NY,NX),NY,NX)+trcsa_RFXS(idsa_K,NU(NY,NX),NY,NX)
  trcsa_XFXS(idsa_OH,NU(NY,NX),NY,NX)=trcsa_XFXS(idsa_OH,NU(NY,NX),NY,NX)+trcsa_RFXS(idsa_OH,NU(NY,NX),NY,NX)
  trcsa_XFXS(idsa_SO4,NU(NY,NX),NY,NX)=trcsa_XFXS(idsa_SO4,NU(NY,NX),NY,NX)+trcsa_RFXS(idsa_SO4,NU(NY,NX),NY,NX)
  trcsa_XFXS(idsa_Cl,NU(NY,NX),NY,NX)=trcsa_XFXS(idsa_Cl,NU(NY,NX),NY,NX)+trcsa_RFXS(idsa_Cl,NU(NY,NX),NY,NX)
  trcsa_XFXS(idsa_CO3,NU(NY,NX),NY,NX)=trcsa_XFXS(idsa_CO3,NU(NY,NX),NY,NX)+trcsa_RFXS(idsa_CO3,NU(NY,NX),NY,NX)
  trcsa_XFXS(idsa_HCO3,NU(NY,NX),NY,NX)=trcsa_XFXS(idsa_HCO3,NU(NY,NX),NY,NX)+trcsa_RFXS(idsa_HCO3,NU(NY,NX),NY,NX)
  trcsa_XFXS(idsa_AlOH,NU(NY,NX),NY,NX)=trcsa_XFXS(idsa_AlOH,NU(NY,NX),NY,NX)+trcsa_RFXS(idsa_AlOH,NU(NY,NX),NY,NX)
  trcsa_XFXS(idsa_AlOH2,NU(NY,NX),NY,NX)=trcsa_XFXS(idsa_AlOH2,NU(NY,NX),NY,NX)+trcsa_RFXS(idsa_AlOH2,NU(NY,NX),NY,NX)
  trcsa_XFXS(idsa_AlOH3,NU(NY,NX),NY,NX)=trcsa_XFXS(idsa_AlOH3,NU(NY,NX),NY,NX)+trcsa_RFXS(idsa_AlOH3,NU(NY,NX),NY,NX)
  trcsa_XFXS(idsa_AlOH4,NU(NY,NX),NY,NX)=trcsa_XFXS(idsa_AlOH4,NU(NY,NX),NY,NX)+trcsa_RFXS(idsa_AlOH4,NU(NY,NX),NY,NX)
  trcsa_XFXS(idsa_AlSO4,NU(NY,NX),NY,NX)=trcsa_XFXS(idsa_AlSO4,NU(NY,NX),NY,NX)+trcsa_RFXS(idsa_AlSO4,NU(NY,NX),NY,NX)
  trcsa_XFXS(idsa_FeOH,NU(NY,NX),NY,NX)=trcsa_XFXS(idsa_FeOH,NU(NY,NX),NY,NX)+trcsa_RFXS(idsa_FeOH,NU(NY,NX),NY,NX)
  trcsa_XFXS(idsa_FeOH2,NU(NY,NX),NY,NX)=trcsa_XFXS(idsa_FeOH2,NU(NY,NX),NY,NX)+trcsa_RFXS(idsa_FeOH2,NU(NY,NX),NY,NX)
  trcsa_XFXS(idsa_FeOH3,NU(NY,NX),NY,NX)=trcsa_XFXS(idsa_FeOH3,NU(NY,NX),NY,NX)+trcsa_RFXS(idsa_FeOH3,NU(NY,NX),NY,NX)
  trcsa_XFXS(idsa_FeOH4,NU(NY,NX),NY,NX)=trcsa_XFXS(idsa_FeOH4,NU(NY,NX),NY,NX)+trcsa_RFXS(idsa_FeOH4,NU(NY,NX),NY,NX)
  trcsa_XFXS(idsa_FeSO4,NU(NY,NX),NY,NX)=trcsa_XFXS(idsa_FeSO4,NU(NY,NX),NY,NX)+trcsa_RFXS(idsa_FeSO4,NU(NY,NX),NY,NX)
  trcsa_XFXS(idsa_CaOH2,NU(NY,NX),NY,NX)=trcsa_XFXS(idsa_CaOH2,NU(NY,NX),NY,NX)+trcsa_RFXS(idsa_CaOH2,NU(NY,NX),NY,NX)
  trcsa_XFXS(idsa_CaCO3,NU(NY,NX),NY,NX)=trcsa_XFXS(idsa_CaCO3,NU(NY,NX),NY,NX)+trcsa_RFXS(idsa_CaCO3,NU(NY,NX),NY,NX)
  trcsa_XFXS(idsa_CaHCO3,NU(NY,NX),NY,NX)=trcsa_XFXS(idsa_CaHCO3,NU(NY,NX),NY,NX)+trcsa_RFXS(idsa_CaHCO3,NU(NY,NX),NY,NX)
  trcsa_XFXS(idsa_CaSO4,NU(NY,NX),NY,NX)=trcsa_XFXS(idsa_CaSO4,NU(NY,NX),NY,NX)+trcsa_RFXS(idsa_CaSO4,NU(NY,NX),NY,NX)
  trcsa_XFXS(idsa_MgOH2,NU(NY,NX),NY,NX)=trcsa_XFXS(idsa_MgOH2,NU(NY,NX),NY,NX)+trcsa_RFXS(idsa_MgOH2,NU(NY,NX),NY,NX)
  trcsa_XFXS(idsa_MgCO3,NU(NY,NX),NY,NX)=trcsa_XFXS(idsa_MgCO3,NU(NY,NX),NY,NX)+trcsa_RFXS(idsa_MgCO3,NU(NY,NX),NY,NX)
  trcsa_XFXS(idsa_MgHCO3,NU(NY,NX),NY,NX)=trcsa_XFXS(idsa_MgHCO3,NU(NY,NX),NY,NX)+trcsa_RFXS(idsa_MgHCO3,NU(NY,NX),NY,NX)
  trcsa_XFXS(idsa_MgSO4,NU(NY,NX),NY,NX)=trcsa_XFXS(idsa_MgSO4,NU(NY,NX),NY,NX)+trcsa_RFXS(idsa_MgSO4,NU(NY,NX),NY,NX)
  trcsa_XFXS(idsa_NaCO3,NU(NY,NX),NY,NX)=trcsa_XFXS(idsa_NaCO3,NU(NY,NX),NY,NX)+trcsa_RFXS(idsa_NaCO3,NU(NY,NX),NY,NX)
  trcsa_XFXS(idsa_NaSO4,NU(NY,NX),NY,NX)=trcsa_XFXS(idsa_NaSO4,NU(NY,NX),NY,NX)+trcsa_RFXS(idsa_NaSO4,NU(NY,NX),NY,NX)
  trcsa_XFXS(idsa_KSO4,NU(NY,NX),NY,NX)=trcsa_XFXS(idsa_KSO4,NU(NY,NX),NY,NX)+trcsa_RFXS(idsa_KSO4,NU(NY,NX),NY,NX)
  trcsa_XFXS(idsa_H0PO4,NU(NY,NX),NY,NX)=trcsa_XFXS(idsa_H0PO4,NU(NY,NX),NY,NX)+trcsa_RFXS(idsa_H0PO4,NU(NY,NX),NY,NX)
  trcsa_XFXS(idsa_H3PO4,NU(NY,NX),NY,NX)=trcsa_XFXS(idsa_H3PO4,NU(NY,NX),NY,NX)+trcsa_RFXS(idsa_H3PO4,NU(NY,NX),NY,NX)
  trcsa_XFXS(idsa_FeHPO4,NU(NY,NX),NY,NX)=trcsa_XFXS(idsa_FeHPO4,NU(NY,NX),NY,NX)+trcsa_RFXS(idsa_FeHPO4,NU(NY,NX),NY,NX)
  trcsa_XFXS(idsa_FeH2PO4,NU(NY,NX),NY,NX)=trcsa_XFXS(idsa_FeH2PO4,NU(NY,NX),NY,NX)+trcsa_RFXS(idsa_FeH2PO4,NU(NY,NX),NY,NX)
  trcsa_XFXS(idsa_CaPO4,NU(NY,NX),NY,NX)=trcsa_XFXS(idsa_CaPO4,NU(NY,NX),NY,NX)+trcsa_RFXS(idsa_CaPO4,NU(NY,NX),NY,NX)
  trcsa_XFXS(idsa_CaHPO4,NU(NY,NX),NY,NX)=trcsa_XFXS(idsa_CaHPO4,NU(NY,NX),NY,NX)+trcsa_RFXS(idsa_CaHPO4,NU(NY,NX),NY,NX)
  trcsa_XFXS(idsa_CaH2PO4,NU(NY,NX),NY,NX)=trcsa_XFXS(idsa_CaH2PO4,NU(NY,NX),NY,NX)+trcsa_RFXS(idsa_CaH2PO4,NU(NY,NX),NY,NX)
  trcsa_XFXS(idsa_MgHPO4,NU(NY,NX),NY,NX)=trcsa_XFXS(idsa_MgHPO4,NU(NY,NX),NY,NX)+trcsa_RFXS(idsa_MgHPO4,NU(NY,NX),NY,NX)
  trcsa_XFXS(idsa_H0PO4B,NU(NY,NX),NY,NX)=trcsa_XFXS(idsa_H0PO4B,NU(NY,NX),NY,NX)+trcsa_RFXS(idsa_H0PO4B,NU(NY,NX),NY,NX)
  trcsa_XFXS(idsa_H3PO4B,NU(NY,NX),NY,NX)=trcsa_XFXS(idsa_H3PO4B,NU(NY,NX),NY,NX)+trcsa_RFXS(idsa_H3PO4B,NU(NY,NX),NY,NX)
  trcsa_XFXS(idsa_FeHPO4B,NU(NY,NX),NY,NX)=trcsa_XFXS(idsa_FeHPO4B,NU(NY,NX),NY,NX)+trcsa_RFXS(idsa_FeHPO4B,NU(NY,NX),NY,NX)
  trcsa_XFXS(idsa_FeH2PO4B,NU(NY,NX),NY,NX)=trcsa_XFXS(idsa_FeH2PO4B,NU(NY,NX),NY,NX)+trcsa_RFXS(idsa_FeH2PO4B,NU(NY,NX),NY,NX)
  trcsa_XFXS(idsa_CaPO4B,NU(NY,NX),NY,NX)=trcsa_XFXS(idsa_CaPO4B,NU(NY,NX),NY,NX)+trcsa_RFXS(idsa_CaPO4B,NU(NY,NX),NY,NX)
  trcsa_XFXS(idsa_CaHPO4B,NU(NY,NX),NY,NX)=trcsa_XFXS(idsa_CaHPO4B,NU(NY,NX),NY,NX)+trcsa_RFXS(idsa_CaHPO4B,NU(NY,NX),NY,NX)
  trcsa_XFXS(idsa_CaH2PO4B,NU(NY,NX),NY,NX)=trcsa_XFXS(idsa_CaH2PO4B,NU(NY,NX),NY,NX)+trcsa_RFXS(idsa_CaH2PO4B,NU(NY,NX),NY,NX)
  trcsa_XFXS(idsa_MgHPO4B,NU(NY,NX),NY,NX)=trcsa_XFXS(idsa_MgHPO4B,NU(NY,NX),NY,NX)+trcsa_RFXS(idsa_MgHPO4B,NU(NY,NX),NY,NX)
  end subroutine AccumHourlyMicMacPoreFlux
!------------------------------------------------------------------------------------------

  subroutine SoluteFluxBySurfaceOutflow(M,N1,N2)
!
!     Description:
!
  implicit none
  integer, intent(in) :: M,N1,N2
  real(r8) :: VFLW

!     begin_execution

  IF(VOLWM(M,0,N2,N1).GT.ZEROS2(N2,N1))THEN
    VFLW=AMIN1(VFLWX,QRM(M,N2,N1)/VOLWM(M,0,N2,N1))
  ELSE
    VFLW=VFLWX
  ENDIF
  trcsa_RQR0(idsa_Al,N2,N1)=VFLW*AZMAX1(trcsa_solml2(idsa_Al,0,N2,N1))
  trcsa_RQR0(idsa_Fe,N2,N1)=VFLW*AZMAX1(trcsa_solml2(idsa_Fe,0,N2,N1))
  trcsa_RQR0(idsa_Hp,N2,N1)=VFLW*AZMAX1(trcsa_solml2(idsa_Hp,0,N2,N1))
  trcsa_RQR0(idsa_Ca,N2,N1)=VFLW*AZMAX1(trcsa_solml2(idsa_Ca,0,N2,N1))
  trcsa_RQR0(idsa_Mg,N2,N1)=VFLW*AZMAX1(trcsa_solml2(idsa_Mg,0,N2,N1))
  trcsa_RQR0(idsa_Na,N2,N1)=VFLW*AZMAX1(trcsa_solml2(idsa_Na,0,N2,N1))
  trcsa_RQR0(idsa_K,N2,N1)=VFLW*AZMAX1(trcsa_solml2(idsa_K,0,N2,N1))
  trcsa_RQR0(idsa_OH,N2,N1)=VFLW*AZMAX1(trcsa_solml2(idsa_OH,0,N2,N1))
  trcsa_RQR0(idsa_SO4,N2,N1)=VFLW*AZMAX1(trcsa_solml2(idsa_SO4,0,N2,N1))
  trcsa_RQR0(idsa_Cl,N2,N1)=VFLW*AZMAX1(trcsa_solml2(idsa_Cl,0,N2,N1))
  trcsa_RQR0(idsa_CO3,N2,N1)=VFLW*AZMAX1(trcsa_solml2(idsa_CO3,0,N2,N1))
  trcsa_RQR0(idsa_HCO3,N2,N1)=VFLW*AZMAX1(trcsa_solml2(idsa_HCO3,0,N2,N1))
  trcsa_RQR0(idsa_AlOH,N2,N1)=VFLW*AZMAX1(trcsa_solml2(idsa_AlOH,0,N2,N1))
  trcsa_RQR0(idsa_AlOH2,N2,N1)=VFLW*AZMAX1(trcsa_solml2(idsa_AlOH2,0,N2,N1))
  trcsa_RQR0(idsa_AlOH3,N2,N1)=VFLW*AZMAX1(trcsa_solml2(idsa_AlOH3,0,N2,N1))
  trcsa_RQR0(idsa_AlOH4,N2,N1)=VFLW*AZMAX1(trcsa_solml2(idsa_AlOH4,0,N2,N1))
  trcsa_RQR0(idsa_AlSO4,N2,N1)=VFLW*AZMAX1(trcsa_solml2(idsa_AlSO4,0,N2,N1))
  trcsa_RQR0(idsa_FeOH,N2,N1)=VFLW*AZMAX1(trcsa_solml2(idsa_FeOH,0,N2,N1))
  trcsa_RQR0(idsa_FeOH2,N2,N1)=VFLW*AZMAX1(trcsa_solml2(idsa_FeOH2,0,N2,N1))
  trcsa_RQR0(idsa_FeOH3,N2,N1)=VFLW*AZMAX1(trcsa_solml2(idsa_FeOH3,0,N2,N1))
  trcsa_RQR0(idsa_FeOH4,N2,N1)=VFLW*AZMAX1(trcsa_solml2(idsa_FeOH4,0,N2,N1))
  trcsa_RQR0(idsa_FeSO4,N2,N1)=VFLW*AZMAX1(trcsa_solml2(idsa_FeSO4,0,N2,N1))
  trcsa_RQR0(idsa_CaOH2,N2,N1)=VFLW*AZMAX1(trcsa_solml2(idsa_CaOH2,0,N2,N1))
  trcsa_RQR0(idsa_CaCO3,N2,N1)=VFLW*AZMAX1(trcsa_solml2(idsa_CaCO3,0,N2,N1))
  trcsa_RQR0(idsa_CaHCO3,N2,N1)=VFLW*AZMAX1(trcsa_solml2(idsa_CaHCO3,0,N2,N1))
  trcsa_RQR0(idsa_CaSO4,N2,N1)=VFLW*AZMAX1(trcsa_solml2(idsa_CaSO4,0,N2,N1))
  trcsa_RQR0(idsa_MgOH2,N2,N1)=VFLW*AZMAX1(trcsa_solml2(idsa_MgOH2,0,N2,N1))
  trcsa_RQR0(idsa_MgCO3,N2,N1)=VFLW*AZMAX1(trcsa_solml2(idsa_MgCO3,0,N2,N1))
  trcsa_RQR0(idsa_MgHCO3,N2,N1)=VFLW*AZMAX1(trcsa_solml2(idsa_MgHCO3,0,N2,N1))
  trcsa_RQR0(idsa_MgSO4,N2,N1)=VFLW*AZMAX1(trcsa_solml2(idsa_MgSO4,0,N2,N1))
  trcsa_RQR0(idsa_NaCO3,N2,N1)=VFLW*AZMAX1(trcsa_solml2(idsa_NaCO3,0,N2,N1))
  trcsa_RQR0(idsa_NaSO4,N2,N1)=VFLW*AZMAX1(trcsa_solml2(idsa_NaSO4,0,N2,N1))
  trcsa_RQR0(idsa_KSO4,N2,N1)=VFLW*AZMAX1(trcsa_solml2(idsa_KSO4,0,N2,N1))
  trcsa_RQR0(idsa_H0PO4,N2,N1)=VFLW*AZMAX1(trcsa_solml2(idsa_H0PO4,0,N2,N1))
  trcsa_RQR0(idsa_H3PO4,N2,N1)=VFLW*AZMAX1(trcsa_solml2(idsa_H3PO4,0,N2,N1))
  trcsa_RQR0(idsa_FeHPO4,N2,N1)=VFLW*AZMAX1(trcsa_solml2(idsa_FeHPO4,0,N2,N1))
  trcsa_RQR0(idsa_FeH2PO4,N2,N1)=VFLW*AZMAX1(trcsa_solml2(idsa_FeH2PO4,0,N2,N1))
  trcsa_RQR0(idsa_CaPO4,N2,N1)=VFLW*AZMAX1(trcsa_solml2(idsa_CaPO4,0,N2,N1))
  trcsa_RQR0(idsa_CaHPO4,N2,N1)=VFLW*AZMAX1(trcsa_solml2(idsa_CaHPO4,0,N2,N1))
  trcsa_RQR0(idsa_CaH2PO4,N2,N1)=VFLW*AZMAX1(trcsa_solml2(idsa_CaH2PO4,0,N2,N1))
  trcsa_RQR0(idsa_MgHPO4,N2,N1)=VFLW*AZMAX1(trcsa_solml2(idsa_MgHPO4,0,N2,N1))
  end subroutine SoluteFluxBySurfaceOutflow
!------------------------------------------------------------------------------------------

  subroutine ZeroSoluteFluxbySurfaceflow(N1,N2)
!
!     Description:
!
  implicit none
  integer, intent(in) :: N1,N2
!     begin_execution

  trcsa_RQR0(idsa_Al,N2,N1)=0.0_r8
  trcsa_RQR0(idsa_Fe,N2,N1)=0.0_r8
  trcsa_RQR0(idsa_Hp,N2,N1)=0.0_r8
  trcsa_RQR0(idsa_Ca,N2,N1)=0.0_r8
  trcsa_RQR0(idsa_Mg,N2,N1)=0.0_r8
  trcsa_RQR0(idsa_Na,N2,N1)=0.0_r8
  trcsa_RQR0(idsa_K,N2,N1)=0.0_r8
  trcsa_RQR0(idsa_OH,N2,N1)=0.0_r8
  trcsa_RQR0(idsa_SO4,N2,N1)=0.0_r8
  trcsa_RQR0(idsa_Cl,N2,N1)=0.0_r8
  trcsa_RQR0(idsa_CO3,N2,N1)=0.0_r8
  trcsa_RQR0(idsa_HCO3,N2,N1)=0.0_r8
  trcsa_RQR0(idsa_AlOH,N2,N1)=0.0_r8
  trcsa_RQR0(idsa_AlOH2,N2,N1)=0.0_r8
  trcsa_RQR0(idsa_AlOH3,N2,N1)=0.0_r8
  trcsa_RQR0(idsa_AlOH4,N2,N1)=0.0_r8
  trcsa_RQR0(idsa_AlSO4,N2,N1)=0.0_r8
  trcsa_RQR0(idsa_FeOH,N2,N1)=0.0_r8
  trcsa_RQR0(idsa_FeOH2,N2,N1)=0.0_r8
  trcsa_RQR0(idsa_FeOH3,N2,N1)=0.0_r8
  trcsa_RQR0(idsa_FeOH4,N2,N1)=0.0_r8
  trcsa_RQR0(idsa_FeSO4,N2,N1)=0.0_r8
  trcsa_RQR0(idsa_CaOH2,N2,N1)=0.0_r8
  trcsa_RQR0(idsa_CaCO3,N2,N1)=0.0_r8
  trcsa_RQR0(idsa_CaHCO3,N2,N1)=0.0_r8
  trcsa_RQR0(idsa_CaSO4,N2,N1)=0.0_r8
  trcsa_RQR0(idsa_MgOH2,N2,N1)=0.0_r8
  trcsa_RQR0(idsa_MgCO3,N2,N1)=0.0_r8
  trcsa_RQR0(idsa_MgHCO3,N2,N1)=0.0_r8
  trcsa_RQR0(idsa_MgSO4,N2,N1)=0.0_r8
  trcsa_RQR0(idsa_NaCO3,N2,N1)=0.0_r8
  trcsa_RQR0(idsa_NaSO4,N2,N1)=0.0_r8
  trcsa_RQR0(idsa_KSO4,N2,N1)=0.0_r8
  trcsa_RQR0(idsa_H0PO4,N2,N1)=0.0_r8
  trcsa_RQR0(idsa_H3PO4,N2,N1)=0.0_r8
  trcsa_RQR0(idsa_FeHPO4,N2,N1)=0.0_r8
  trcsa_RQR0(idsa_FeH2PO4,N2,N1)=0.0_r8
  trcsa_RQR0(idsa_CaPO4,N2,N1)=0.0_r8
  trcsa_RQR0(idsa_CaHPO4,N2,N1)=0.0_r8
  trcsa_RQR0(idsa_CaH2PO4,N2,N1)=0.0_r8
  trcsa_RQR0(idsa_MgHPO4,N2,N1)=0.0_r8
  end subroutine ZeroSoluteFluxbySurfaceflow
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

          trcsa_RQR(idsa_Al,N,2,N5,N4)=trcsa_RQR0(idsa_Al,N2,N1)*FQRM
          trcsa_RQR(idsa_Fe,N,2,N5,N4)=trcsa_RQR0(idsa_Fe,N2,N1)*FQRM
          trcsa_RQR(idsa_Hp,N,2,N5,N4)=trcsa_RQR0(idsa_Hp,N2,N1)*FQRM
          trcsa_RQR(idsa_Ca,N,2,N5,N4)=trcsa_RQR0(idsa_Ca,N2,N1)*FQRM
          trcsa_RQR(idsa_Mg,N,2,N5,N4)=trcsa_RQR0(idsa_Mg,N2,N1)*FQRM
          trcsa_RQR(idsa_Na,N,2,N5,N4)=trcsa_RQR0(idsa_Na,N2,N1)*FQRM
          trcsa_RQR(idsa_K,N,2,N5,N4)=trcsa_RQR0(idsa_K,N2,N1)*FQRM
          trcsa_RQR(idsa_OH,N,2,N5,N4)=trcsa_RQR0(idsa_OH,N2,N1)*FQRM
          trcsa_RQR(idsa_SO4,N,2,N5,N4)=trcsa_RQR0(idsa_SO4,N2,N1)*FQRM
          trcsa_RQR(idsa_Cl,N,2,N5,N4)=trcsa_RQR0(idsa_Cl,N2,N1)*FQRM
          trcsa_RQR(idsa_CO3,N,2,N5,N4)=trcsa_RQR0(idsa_CO3,N2,N1)*FQRM
          trcsa_RQR(idsa_HCO3,N,2,N5,N4)=trcsa_RQR0(idsa_HCO3,N2,N1)*FQRM
          trcsa_RQR(idsa_AlOH,N,2,N5,N4)=trcsa_RQR0(idsa_AlOH,N2,N1)*FQRM
          trcsa_RQR(idsa_AlOH2,N,2,N5,N4)=trcsa_RQR0(idsa_AlOH2,N2,N1)*FQRM
          trcsa_RQR(idsa_AlOH3,N,2,N5,N4)=trcsa_RQR0(idsa_AlOH3,N2,N1)*FQRM
          trcsa_RQR(idsa_AlOH4,N,2,N5,N4)=trcsa_RQR0(idsa_AlOH4,N2,N1)*FQRM
          trcsa_RQR(idsa_AlSO4,N,2,N5,N4)=trcsa_RQR0(idsa_AlSO4,N2,N1)*FQRM
          trcsa_RQR(idsa_FeOH,N,2,N5,N4)=trcsa_RQR0(idsa_FeOH,N2,N1)*FQRM
          trcsa_RQR(idsa_FeOH2,N,2,N5,N4)=trcsa_RQR0(idsa_FeOH2,N2,N1)*FQRM
          trcsa_RQR(idsa_FeOH3,N,2,N5,N4)=trcsa_RQR0(idsa_FeOH3,N2,N1)*FQRM
          trcsa_RQR(idsa_FeOH4,N,2,N5,N4)=trcsa_RQR0(idsa_FeOH4,N2,N1)*FQRM
          trcsa_RQR(idsa_FeSO4,N,2,N5,N4)=trcsa_RQR0(idsa_FeSO4,N2,N1)*FQRM
          trcsa_RQR(idsa_CaOH2,N,2,N5,N4)=trcsa_RQR0(idsa_CaOH2,N2,N1)*FQRM
          trcsa_RQR(idsa_CaCO3,N,2,N5,N4)=trcsa_RQR0(idsa_CaCO3,N2,N1)*FQRM
          trcsa_RQR(idsa_CaHCO3,N,2,N5,N4)=trcsa_RQR0(idsa_CaHCO3,N2,N1)*FQRM
          trcsa_RQR(idsa_CaSO4,N,2,N5,N4)=trcsa_RQR0(idsa_CaSO4,N2,N1)*FQRM
          trcsa_RQR(idsa_MgOH2,N,2,N5,N4)=trcsa_RQR0(idsa_MgOH2,N2,N1)*FQRM
          trcsa_RQR(idsa_MgCO3,N,2,N5,N4)=trcsa_RQR0(idsa_MgCO3,N2,N1)*FQRM
          trcsa_RQR(idsa_MgHCO3,N,2,N5,N4)=trcsa_RQR0(idsa_MgHCO3,N2,N1)*FQRM
          trcsa_RQR(idsa_MgSO4,N,2,N5,N4)=trcsa_RQR0(idsa_MgSO4,N2,N1)*FQRM
          trcsa_RQR(idsa_NaCO3,N,2,N5,N4)=trcsa_RQR0(idsa_NaCO3,N2,N1)*FQRM
          trcsa_RQR(idsa_NaSO4,N,2,N5,N4)=trcsa_RQR0(idsa_NaSO4,N2,N1)*FQRM
          trcsa_RQR(idsa_KSO4,N,2,N5,N4)=trcsa_RQR0(idsa_KSO4,N2,N1)*FQRM
          trcsa_RQR(idsa_H0PO4,N,2,N5,N4)=trcsa_RQR0(idsa_H0PO4,N2,N1)*FQRM
          trcsa_RQR(idsa_H3PO4,N,2,N5,N4)=trcsa_RQR0(idsa_H3PO4,N2,N1)*FQRM
          trcsa_RQR(idsa_FeHPO4,N,2,N5,N4)=trcsa_RQR0(idsa_FeHPO4,N2,N1)*FQRM
          trcsa_RQR(idsa_FeH2PO4,N,2,N5,N4)=trcsa_RQR0(idsa_FeH2PO4,N2,N1)*FQRM
          trcsa_RQR(idsa_CaPO4,N,2,N5,N4)=trcsa_RQR0(idsa_CaPO4,N2,N1)*FQRM
          trcsa_RQR(idsa_CaHPO4,N,2,N5,N4)=trcsa_RQR0(idsa_CaHPO4,N2,N1)*FQRM
          trcsa_RQR(idsa_CaH2PO4,N,2,N5,N4)=trcsa_RQR0(idsa_CaH2PO4,N2,N1)*FQRM
          trcsa_RQR(idsa_MgHPO4,N,2,N5,N4)=trcsa_RQR0(idsa_MgHPO4,N2,N1)*FQRM
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
            trcsa_RQR(idsa_Al,N,1,N5B,N4B)=trcsa_RQR0(idsa_Al,N2,N1)*FQRM
            trcsa_RQR(idsa_Fe,N,1,N5B,N4B)=trcsa_RQR0(idsa_Fe,N2,N1)*FQRM
            trcsa_RQR(idsa_Hp,N,1,N5B,N4B)=trcsa_RQR0(idsa_Hp,N2,N1)*FQRM
            trcsa_RQR(idsa_Ca,N,1,N5B,N4B)=trcsa_RQR0(idsa_Ca,N2,N1)*FQRM
            trcsa_RQR(idsa_Mg,N,1,N5B,N4B)=trcsa_RQR0(idsa_Mg,N2,N1)*FQRM
            trcsa_RQR(idsa_Na,N,1,N5B,N4B)=trcsa_RQR0(idsa_Na,N2,N1)*FQRM
            trcsa_RQR(idsa_K,N,1,N5B,N4B)=trcsa_RQR0(idsa_K,N2,N1)*FQRM
            trcsa_RQR(idsa_OH,N,1,N5B,N4B)=trcsa_RQR0(idsa_OH,N2,N1)*FQRM
            trcsa_RQR(idsa_SO4,N,1,N5B,N4B)=trcsa_RQR0(idsa_SO4,N2,N1)*FQRM
            trcsa_RQR(idsa_Cl,N,1,N5B,N4B)=trcsa_RQR0(idsa_Cl,N2,N1)*FQRM
            trcsa_RQR(idsa_CO3,N,1,N5B,N4B)=trcsa_RQR0(idsa_CO3,N2,N1)*FQRM
            trcsa_RQR(idsa_HCO3,N,1,N5B,N4B)=trcsa_RQR0(idsa_HCO3,N2,N1)*FQRM
            trcsa_RQR(idsa_AlOH,N,1,N5B,N4B)=trcsa_RQR0(idsa_AlOH,N2,N1)*FQRM
            trcsa_RQR(idsa_AlOH2,N,1,N5B,N4B)=trcsa_RQR0(idsa_AlOH2,N2,N1)*FQRM
            trcsa_RQR(idsa_AlOH3,N,1,N5B,N4B)=trcsa_RQR0(idsa_AlOH3,N2,N1)*FQRM
            trcsa_RQR(idsa_AlOH4,N,1,N5B,N4B)=trcsa_RQR0(idsa_AlOH4,N2,N1)*FQRM
            trcsa_RQR(idsa_AlSO4,N,1,N5B,N4B)=trcsa_RQR0(idsa_AlSO4,N2,N1)*FQRM
            trcsa_RQR(idsa_FeOH,N,1,N5B,N4B)=trcsa_RQR0(idsa_FeOH,N2,N1)*FQRM
            trcsa_RQR(idsa_FeOH2,N,1,N5B,N4B)=trcsa_RQR0(idsa_FeOH2,N2,N1)*FQRM
            trcsa_RQR(idsa_FeOH3,N,1,N5B,N4B)=trcsa_RQR0(idsa_FeOH3,N2,N1)*FQRM
            trcsa_RQR(idsa_FeOH4,N,1,N5B,N4B)=trcsa_RQR0(idsa_FeOH4,N2,N1)*FQRM
            trcsa_RQR(idsa_FeSO4,N,1,N5B,N4B)=trcsa_RQR0(idsa_FeSO4,N2,N1)*FQRM
            trcsa_RQR(idsa_CaOH2,N,1,N5B,N4B)=trcsa_RQR0(idsa_CaOH2,N2,N1)*FQRM
            trcsa_RQR(idsa_CaCO3,N,1,N5B,N4B)=trcsa_RQR0(idsa_CaCO3,N2,N1)*FQRM
            trcsa_RQR(idsa_CaHCO3,N,1,N5B,N4B)=trcsa_RQR0(idsa_CaHCO3,N2,N1)*FQRM
            trcsa_RQR(idsa_CaSO4,N,1,N5B,N4B)=trcsa_RQR0(idsa_CaSO4,N2,N1)*FQRM
            trcsa_RQR(idsa_MgOH2,N,1,N5B,N4B)=trcsa_RQR0(idsa_MgOH2,N2,N1)*FQRM
            trcsa_RQR(idsa_MgCO3,N,1,N5B,N4B)=trcsa_RQR0(idsa_MgCO3,N2,N1)*FQRM
            trcsa_RQR(idsa_MgHCO3,N,1,N5B,N4B)=trcsa_RQR0(idsa_MgHCO3,N2,N1)*FQRM
            trcsa_RQR(idsa_MgSO4,N,1,N5B,N4B)=trcsa_RQR0(idsa_MgSO4,N2,N1)*FQRM
            trcsa_RQR(idsa_NaCO3,N,1,N5B,N4B)=trcsa_RQR0(idsa_NaCO3,N2,N1)*FQRM
            trcsa_RQR(idsa_NaSO4,N,1,N5B,N4B)=trcsa_RQR0(idsa_NaSO4,N2,N1)*FQRM
            trcsa_RQR(idsa_KSO4,N,1,N5B,N4B)=trcsa_RQR0(idsa_KSO4,N2,N1)*FQRM
            trcsa_RQR(idsa_H0PO4,N,1,N5B,N4B)=trcsa_RQR0(idsa_H0PO4,N2,N1)*FQRM
            trcsa_RQR(idsa_H3PO4,N,1,N5B,N4B)=trcsa_RQR0(idsa_H3PO4,N2,N1)*FQRM
            trcsa_RQR(idsa_FeHPO4,N,1,N5B,N4B)=trcsa_RQR0(idsa_FeHPO4,N2,N1)*FQRM
            trcsa_RQR(idsa_FeH2PO4,N,1,N5B,N4B)=trcsa_RQR0(idsa_FeH2PO4,N2,N1)*FQRM
            trcsa_RQR(idsa_CaPO4,N,1,N5B,N4B)=trcsa_RQR0(idsa_CaPO4,N2,N1)*FQRM
            trcsa_RQR(idsa_CaHPO4,N,1,N5B,N4B)=trcsa_RQR0(idsa_CaHPO4,N2,N1)*FQRM
            trcsa_RQR(idsa_CaH2PO4,N,1,N5B,N4B)=trcsa_RQR0(idsa_CaH2PO4,N2,N1)*FQRM
            trcsa_RQR(idsa_MgHPO4,N,1,N5B,N4B)=trcsa_RQR0(idsa_MgHPO4,N2,N1)*FQRM
      !
!     ACCUMULATE HOURLY FLUXES FOR USE IN REDIST.F
!
!     XQR*=hourly solute in runoff
!     RQR*=solute in runoff
!
            trcsa_XQR(idsa_Al,N,1,N5B,N4B)=trcsa_XQR(idsa_Al,N,1,N5B,N4B)+trcsa_RQR(idsa_Al,N,1,N5B,N4B)
            trcsa_XQR(idsa_Fe,N,1,N5B,N4B)=trcsa_XQR(idsa_Fe,N,1,N5B,N4B)+trcsa_RQR(idsa_Fe,N,1,N5B,N4B)
            trcsa_XQR(idsa_Hp,N,1,N5B,N4B)=trcsa_XQR(idsa_Hp,N,1,N5B,N4B)+trcsa_RQR(idsa_Hp,N,1,N5B,N4B)
            trcsa_XQR(idsa_Ca,N,1,N5B,N4B)=trcsa_XQR(idsa_Ca,N,1,N5B,N4B)+trcsa_RQR(idsa_Ca,N,1,N5B,N4B)
            trcsa_XQR(idsa_Mg,N,1,N5B,N4B)=trcsa_XQR(idsa_Mg,N,1,N5B,N4B)+trcsa_RQR(idsa_Mg,N,1,N5B,N4B)
            trcsa_XQR(idsa_Na,N,1,N5B,N4B)=trcsa_XQR(idsa_Na,N,1,N5B,N4B)+trcsa_RQR(idsa_Na,N,1,N5B,N4B)
            trcsa_XQR(idsa_K,N,1,N5B,N4B)=trcsa_XQR(idsa_K,N,1,N5B,N4B)+trcsa_RQR(idsa_K,N,1,N5B,N4B)
            trcsa_XQR(idsa_OH,N,1,N5B,N4B)=trcsa_XQR(idsa_OH,N,1,N5B,N4B)+trcsa_RQR(idsa_OH,N,1,N5B,N4B)
            trcsa_XQR(idsa_SO4,N,1,N5B,N4B)=trcsa_XQR(idsa_SO4,N,1,N5B,N4B)+trcsa_RQR(idsa_SO4,N,1,N5B,N4B)
            trcsa_XQR(idsa_Cl,N,1,N5B,N4B)=trcsa_XQR(idsa_Cl,N,1,N5B,N4B)+trcsa_RQR(idsa_Cl,N,1,N5B,N4B)
            trcsa_XQR(idsa_CO3,N,1,N5B,N4B)=trcsa_XQR(idsa_CO3,N,1,N5B,N4B)+trcsa_RQR(idsa_CO3,N,1,N5B,N4B)
            trcsa_XQR(idsa_HCO3,N,1,N5B,N4B)=trcsa_XQR(idsa_HCO3,N,1,N5B,N4B)+trcsa_RQR(idsa_HCO3,N,1,N5B,N4B)
            trcsa_XQR(idsa_AlOH,N,1,N5B,N4B)=trcsa_XQR(idsa_AlOH,N,1,N5B,N4B)+trcsa_RQR(idsa_AlOH,N,1,N5B,N4B)
            trcsa_XQR(idsa_AlOH2,N,1,N5B,N4B)=trcsa_XQR(idsa_AlOH2,N,1,N5B,N4B)+trcsa_RQR(idsa_AlOH2,N,1,N5B,N4B)
            trcsa_XQR(idsa_AlOH3,N,1,N5B,N4B)=trcsa_XQR(idsa_AlOH3,N,1,N5B,N4B)+trcsa_RQR(idsa_AlOH3,N,1,N5B,N4B)
            trcsa_XQR(idsa_AlOH4,N,1,N5B,N4B)=trcsa_XQR(idsa_AlOH4,N,1,N5B,N4B)+trcsa_RQR(idsa_AlOH4,N,1,N5B,N4B)
            trcsa_XQR(idsa_AlSO4,N,1,N5B,N4B)=trcsa_XQR(idsa_AlSO4,N,1,N5B,N4B)+trcsa_RQR(idsa_AlSO4,N,1,N5B,N4B)
            trcsa_XQR(idsa_FeOH,N,1,N5B,N4B)=trcsa_XQR(idsa_FeOH,N,1,N5B,N4B)+trcsa_RQR(idsa_FeOH,N,1,N5B,N4B)
            trcsa_XQR(idsa_FeOH2,N,1,N5B,N4B)=trcsa_XQR(idsa_FeOH2,N,1,N5B,N4B)+trcsa_RQR(idsa_FeOH2,N,1,N5B,N4B)
            trcsa_XQR(idsa_FeOH3,N,1,N5B,N4B)=trcsa_XQR(idsa_FeOH3,N,1,N5B,N4B)+trcsa_RQR(idsa_FeOH3,N,1,N5B,N4B)
            trcsa_XQR(idsa_FeOH4,N,1,N5B,N4B)=trcsa_XQR(idsa_FeOH4,N,1,N5B,N4B)+trcsa_RQR(idsa_FeOH4,N,1,N5B,N4B)
            trcsa_XQR(idsa_FeSO4,N,1,N5B,N4B)=trcsa_XQR(idsa_FeSO4,N,1,N5B,N4B)+trcsa_RQR(idsa_FeSO4,N,1,N5B,N4B)
            trcsa_XQR(idsa_CaOH2,N,1,N5B,N4B)=trcsa_XQR(idsa_CaOH2,N,1,N5B,N4B)+trcsa_RQR(idsa_CaOH2,N,1,N5B,N4B)
            trcsa_XQR(idsa_CaCO3,N,1,N5B,N4B)=trcsa_XQR(idsa_CaCO3,N,1,N5B,N4B)+trcsa_RQR(idsa_CaCO3,N,1,N5B,N4B)
            trcsa_XQR(idsa_CaHCO3,N,1,N5B,N4B)=trcsa_XQR(idsa_CaHCO3,N,1,N5B,N4B)+trcsa_RQR(idsa_CaHCO3,N,1,N5B,N4B)
            trcsa_XQR(idsa_CaSO4,N,1,N5B,N4B)=trcsa_XQR(idsa_CaSO4,N,1,N5B,N4B)+trcsa_RQR(idsa_CaSO4,N,1,N5B,N4B)
            trcsa_XQR(idsa_MgOH2,N,1,N5B,N4B)=trcsa_XQR(idsa_MgOH2,N,1,N5B,N4B)+trcsa_RQR(idsa_MgOH2,N,1,N5B,N4B)
            trcsa_XQR(idsa_MgCO3,N,1,N5B,N4B)=trcsa_XQR(idsa_MgCO3,N,1,N5B,N4B)+trcsa_RQR(idsa_MgCO3,N,1,N5B,N4B)
            trcsa_XQR(idsa_MgHCO3,N,1,N5B,N4B)=trcsa_XQR(idsa_MgHCO3,N,1,N5B,N4B)+trcsa_RQR(idsa_MgHCO3,N,1,N5B,N4B)
            trcsa_XQR(idsa_MgSO4,N,1,N5B,N4B)=trcsa_XQR(idsa_MgSO4,N,1,N5B,N4B)+trcsa_RQR(idsa_MgSO4,N,1,N5B,N4B)
            trcsa_XQR(idsa_NaCO3,N,1,N5B,N4B)=trcsa_XQR(idsa_NaCO3,N,1,N5B,N4B)+trcsa_RQR(idsa_NaCO3,N,1,N5B,N4B)
            trcsa_XQR(idsa_NaSO4,N,1,N5B,N4B)=trcsa_XQR(idsa_NaSO4,N,1,N5B,N4B)+trcsa_RQR(idsa_NaSO4,N,1,N5B,N4B)
            trcsa_XQR(idsa_KSO4,N,1,N5B,N4B)=trcsa_XQR(idsa_KSO4,N,1,N5B,N4B)+trcsa_RQR(idsa_KSO4,N,1,N5B,N4B)
            trcsa_XQR(idsa_H0PO4,N,1,N5B,N4B)=trcsa_XQR(idsa_H0PO4,N,1,N5B,N4B)+trcsa_RQR(idsa_H0PO4,N,1,N5B,N4B)
            trcsa_XQR(idsa_H3PO4,N,1,N5B,N4B)=trcsa_XQR(idsa_H3PO4,N,1,N5B,N4B)+trcsa_RQR(idsa_H3PO4,N,1,N5B,N4B)
            trcsa_XQR(idsa_FeHPO4,N,1,N5B,N4B)=trcsa_XQR(idsa_FeHPO4,N,1,N5B,N4B)+trcsa_RQR(idsa_FeHPO4,N,1,N5B,N4B)
            trcsa_XQR(idsa_FeH2PO4,N,1,N5B,N4B)=trcsa_XQR(idsa_FeH2PO4,N,1,N5B,N4B)+trcsa_RQR(idsa_FeH2PO4,N,1,N5B,N4B)
            trcsa_XQR(idsa_CaPO4,N,1,N5B,N4B)=trcsa_XQR(idsa_CaPO4,N,1,N5B,N4B)+trcsa_RQR(idsa_CaPO4,N,1,N5B,N4B)
            trcsa_XQR(idsa_CaHPO4,N,1,N5B,N4B)=trcsa_XQR(idsa_CaHPO4,N,1,N5B,N4B)+trcsa_RQR(idsa_CaHPO4,N,1,N5B,N4B)
            trcsa_XQR(idsa_CaH2PO4,N,1,N5B,N4B)=trcsa_XQR(idsa_CaH2PO4,N,1,N5B,N4B)+trcsa_RQR(idsa_CaH2PO4,N,1,N5B,N4B)
            trcsa_XQR(idsa_MgHPO4,N,1,N5B,N4B)=trcsa_XQR(idsa_MgHPO4,N,1,N5B,N4B)+trcsa_RQR(idsa_MgHPO4,N,1,N5B,N4B)
          ELSE
            trcsa_RQR(idsa_Al,N,1,N5B,N4B)=0.0_r8
            trcsa_RQR(idsa_Fe,N,1,N5B,N4B)=0.0_r8
            trcsa_RQR(idsa_Hp,N,1,N5B,N4B)=0.0_r8
            trcsa_RQR(idsa_Ca,N,1,N5B,N4B)=0.0_r8
            trcsa_RQR(idsa_Mg,N,1,N5B,N4B)=0.0_r8
            trcsa_RQR(idsa_Na,N,1,N5B,N4B)=0.0_r8
            trcsa_RQR(idsa_K,N,1,N5B,N4B)=0.0_r8
            trcsa_RQR(idsa_OH,N,1,N5B,N4B)=0.0_r8
            trcsa_RQR(idsa_SO4,N,1,N5B,N4B)=0.0_r8
            trcsa_RQR(idsa_Cl,N,1,N5B,N4B)=0.0_r8
            trcsa_RQR(idsa_CO3,N,1,N5B,N4B)=0.0_r8
            trcsa_RQR(idsa_HCO3,N,1,N5B,N4B)=0.0_r8
            trcsa_RQR(idsa_AlOH,N,1,N5B,N4B)=0.0_r8
            trcsa_RQR(idsa_AlOH2,N,1,N5B,N4B)=0.0_r8
            trcsa_RQR(idsa_AlOH3,N,1,N5B,N4B)=0.0_r8
            trcsa_RQR(idsa_AlOH4,N,1,N5B,N4B)=0.0_r8
            trcsa_RQR(idsa_AlSO4,N,1,N5B,N4B)=0.0_r8
            trcsa_RQR(idsa_FeOH,N,1,N5B,N4B)=0.0_r8
            trcsa_RQR(idsa_FeOH2,N,1,N5B,N4B)=0.0_r8
            trcsa_RQR(idsa_FeOH3,N,1,N5B,N4B)=0.0_r8
            trcsa_RQR(idsa_FeOH4,N,1,N5B,N4B)=0.0_r8
            trcsa_RQR(idsa_FeSO4,N,1,N5B,N4B)=0.0_r8
            trcsa_RQR(idsa_CaOH2,N,1,N5B,N4B)=0.0_r8
            trcsa_RQR(idsa_CaCO3,N,1,N5B,N4B)=0.0_r8
            trcsa_RQR(idsa_CaHCO3,N,1,N5B,N4B)=0.0_r8
            trcsa_RQR(idsa_CaSO4,N,1,N5B,N4B)=0.0_r8
            trcsa_RQR(idsa_MgOH2,N,1,N5B,N4B)=0.0_r8
            trcsa_RQR(idsa_MgCO3,N,1,N5B,N4B)=0.0_r8
            trcsa_RQR(idsa_MgHCO3,N,1,N5B,N4B)=0.0_r8
            trcsa_RQR(idsa_MgSO4,N,1,N5B,N4B)=0.0_r8
            trcsa_RQR(idsa_NaCO3,N,1,N5B,N4B)=0.0_r8
            trcsa_RQR(idsa_NaSO4,N,1,N5B,N4B)=0.0_r8
            trcsa_RQR(idsa_KSO4,N,1,N5B,N4B)=0.0_r8
            trcsa_RQR(idsa_H0PO4,N,1,N5B,N4B)=0.0_r8
            trcsa_RQR(idsa_H3PO4,N,1,N5B,N4B)=0.0_r8
            trcsa_RQR(idsa_FeHPO4,N,1,N5B,N4B)=0.0_r8
            trcsa_RQR(idsa_FeH2PO4,N,1,N5B,N4B)=0.0_r8
            trcsa_RQR(idsa_CaPO4,N,1,N5B,N4B)=0.0_r8
            trcsa_RQR(idsa_CaHPO4,N,1,N5B,N4B)=0.0_r8
            trcsa_RQR(idsa_CaH2PO4,N,1,N5B,N4B)=0.0_r8
            trcsa_RQR(idsa_MgHPO4,N,1,N5B,N4B)=0.0_r8
          ENDIF
        ENDIF
      ELSE
        trcsa_RQR(idsa_Al,N,2,N5,N4)=0.0_r8
        trcsa_RQR(idsa_Fe,N,2,N5,N4)=0.0_r8
        trcsa_RQR(idsa_Hp,N,2,N5,N4)=0.0_r8
        trcsa_RQR(idsa_Ca,N,2,N5,N4)=0.0_r8
        trcsa_RQR(idsa_Mg,N,2,N5,N4)=0.0_r8
        trcsa_RQR(idsa_Na,N,2,N5,N4)=0.0_r8
        trcsa_RQR(idsa_K,N,2,N5,N4)=0.0_r8
        trcsa_RQR(idsa_OH,N,2,N5,N4)=0.0_r8
        trcsa_RQR(idsa_SO4,N,2,N5,N4)=0.0_r8
        trcsa_RQR(idsa_Cl,N,2,N5,N4)=0.0_r8
        trcsa_RQR(idsa_CO3,N,2,N5,N4)=0.0_r8
        trcsa_RQR(idsa_HCO3,N,2,N5,N4)=0.0_r8
        trcsa_RQR(idsa_AlOH,N,2,N5,N4)=0.0_r8
        trcsa_RQR(idsa_AlOH2,N,2,N5,N4)=0.0_r8
        trcsa_RQR(idsa_AlOH3,N,2,N5,N4)=0.0_r8
        trcsa_RQR(idsa_AlOH4,N,2,N5,N4)=0.0_r8
        trcsa_RQR(idsa_AlSO4,N,2,N5,N4)=0.0_r8
        trcsa_RQR(idsa_FeOH,N,2,N5,N4)=0.0_r8
        trcsa_RQR(idsa_FeOH2,N,2,N5,N4)=0.0_r8
        trcsa_RQR(idsa_FeOH3,N,2,N5,N4)=0.0_r8
        trcsa_RQR(idsa_FeOH4,N,2,N5,N4)=0.0_r8
        trcsa_RQR(idsa_FeSO4,N,2,N5,N4)=0.0_r8
        trcsa_RQR(idsa_CaOH2,N,2,N5,N4)=0.0_r8
        trcsa_RQR(idsa_CaCO3,N,2,N5,N4)=0.0_r8
        trcsa_RQR(idsa_CaHCO3,N,2,N5,N4)=0.0_r8
        trcsa_RQR(idsa_CaSO4,N,2,N5,N4)=0.0_r8
        trcsa_RQR(idsa_MgOH2,N,2,N5,N4)=0.0_r8
        trcsa_RQR(idsa_MgCO3,N,2,N5,N4)=0.0_r8
        trcsa_RQR(idsa_MgHCO3,N,2,N5,N4)=0.0_r8
        trcsa_RQR(idsa_MgSO4,N,2,N5,N4)=0.0_r8
        trcsa_RQR(idsa_NaCO3,N,2,N5,N4)=0.0_r8
        trcsa_RQR(idsa_NaSO4,N,2,N5,N4)=0.0_r8
        trcsa_RQR(idsa_KSO4,N,2,N5,N4)=0.0_r8
        trcsa_RQR(idsa_H0PO4,N,2,N5,N4)=0.0_r8
        trcsa_RQR(idsa_H3PO4,N,2,N5,N4)=0.0_r8
        trcsa_RQR(idsa_FeHPO4,N,2,N5,N4)=0.0_r8
        trcsa_RQR(idsa_FeH2PO4,N,2,N5,N4)=0.0_r8
        trcsa_RQR(idsa_CaPO4,N,2,N5,N4)=0.0_r8
        trcsa_RQR(idsa_CaHPO4,N,2,N5,N4)=0.0_r8
        trcsa_RQR(idsa_CaH2PO4,N,2,N5,N4)=0.0_r8
        trcsa_RQR(idsa_MgHPO4,N,2,N5,N4)=0.0_r8
        IF(N4B.GT.0.AND.N5B.GT.0)THEN
          trcsa_RQR(idsa_Al,N,1,N5B,N4B)=0.0_r8
          trcsa_RQR(idsa_Fe,N,1,N5B,N4B)=0.0_r8
          trcsa_RQR(idsa_Hp,N,1,N5B,N4B)=0.0_r8
          trcsa_RQR(idsa_Ca,N,1,N5B,N4B)=0.0_r8
          trcsa_RQR(idsa_Mg,N,1,N5B,N4B)=0.0_r8
          trcsa_RQR(idsa_Na,N,1,N5B,N4B)=0.0_r8
          trcsa_RQR(idsa_K,N,1,N5B,N4B)=0.0_r8
          trcsa_RQR(idsa_OH,N,1,N5B,N4B)=0.0_r8
          trcsa_RQR(idsa_SO4,N,1,N5B,N4B)=0.0_r8
          trcsa_RQR(idsa_Cl,N,1,N5B,N4B)=0.0_r8
          trcsa_RQR(idsa_CO3,N,1,N5B,N4B)=0.0_r8
          trcsa_RQR(idsa_HCO3,N,1,N5B,N4B)=0.0_r8
          trcsa_RQR(idsa_AlOH,N,1,N5B,N4B)=0.0_r8
          trcsa_RQR(idsa_AlOH2,N,1,N5B,N4B)=0.0_r8
          trcsa_RQR(idsa_AlOH3,N,1,N5B,N4B)=0.0_r8
          trcsa_RQR(idsa_AlOH4,N,1,N5B,N4B)=0.0_r8
          trcsa_RQR(idsa_AlSO4,N,1,N5B,N4B)=0.0_r8
          trcsa_RQR(idsa_FeOH,N,1,N5B,N4B)=0.0_r8
          trcsa_RQR(idsa_FeOH2,N,1,N5B,N4B)=0.0_r8
          trcsa_RQR(idsa_FeOH3,N,1,N5B,N4B)=0.0_r8
          trcsa_RQR(idsa_FeOH4,N,1,N5B,N4B)=0.0_r8
          trcsa_RQR(idsa_FeSO4,N,1,N5B,N4B)=0.0_r8
          trcsa_RQR(idsa_CaOH2,N,1,N5B,N4B)=0.0_r8
          trcsa_RQR(idsa_CaCO3,N,1,N5B,N4B)=0.0_r8
          trcsa_RQR(idsa_CaHCO3,N,1,N5B,N4B)=0.0_r8
          trcsa_RQR(idsa_CaSO4,N,1,N5B,N4B)=0.0_r8
          trcsa_RQR(idsa_MgOH2,N,1,N5B,N4B)=0.0_r8
          trcsa_RQR(idsa_MgCO3,N,1,N5B,N4B)=0.0_r8
          trcsa_RQR(idsa_MgHCO3,N,1,N5B,N4B)=0.0_r8
          trcsa_RQR(idsa_MgSO4,N,1,N5B,N4B)=0.0_r8
          trcsa_RQR(idsa_NaCO3,N,1,N5B,N4B)=0.0_r8
          trcsa_RQR(idsa_NaSO4,N,1,N5B,N4B)=0.0_r8
          trcsa_RQR(idsa_KSO4,N,1,N5B,N4B)=0.0_r8
          trcsa_RQR(idsa_H0PO4,N,1,N5B,N4B)=0.0_r8
          trcsa_RQR(idsa_H3PO4,N,1,N5B,N4B)=0.0_r8
          trcsa_RQR(idsa_FeHPO4,N,1,N5B,N4B)=0.0_r8
          trcsa_RQR(idsa_FeH2PO4,N,1,N5B,N4B)=0.0_r8
          trcsa_RQR(idsa_CaPO4,N,1,N5B,N4B)=0.0_r8
          trcsa_RQR(idsa_CaHPO4,N,1,N5B,N4B)=0.0_r8
          trcsa_RQR(idsa_CaH2PO4,N,1,N5B,N4B)=0.0_r8
          trcsa_RQR(idsa_MgHPO4,N,1,N5B,N4B)=0.0_r8
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
            trcsa_RQ(NTSA,N,N5,N4)=VFLW*AZMAX1(trcs_solsml2(NTSA,1,N2,N1))
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
            trcsa_RQ(NTSA,N,N5,N4)=VFLW*AZMAX1(trcs_solsml2(NTSA,1,N5,N4))
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
    RFLAL=VFLW*AZMAX1(trcsa_solml2(idsa_Al,N3,N2,N1))
    RFLFE=VFLW*AZMAX1(trcsa_solml2(idsa_Fe,N3,N2,N1))
    RFLHY=VFLW*AZMAX1(trcsa_solml2(idsa_Hp,N3,N2,N1))
    RFLCA=VFLW*AZMAX1(trcsa_solml2(idsa_Ca,N3,N2,N1))
    RFLMG=VFLW*AZMAX1(trcsa_solml2(idsa_Mg,N3,N2,N1))
    RFLNA=VFLW*AZMAX1(trcsa_solml2(idsa_Na,N3,N2,N1))
    RFLKA=VFLW*AZMAX1(trcsa_solml2(idsa_K,N3,N2,N1))
    RFLOH=VFLW*AZMAX1(trcsa_solml2(idsa_OH,N3,N2,N1))
    RFLSO=VFLW*AZMAX1(trcsa_solml2(idsa_SO4,N3,N2,N1))
    RFLCL=VFLW*AZMAX1(trcsa_solml2(idsa_Cl,N3,N2,N1))
    RFLC3=VFLW*AZMAX1(trcsa_solml2(idsa_CO3,N3,N2,N1))
    RFLHC=VFLW*AZMAX1(trcsa_solml2(idsa_HCO3,N3,N2,N1))
    RFLAL1=VFLW*AZMAX1(trcsa_solml2(idsa_AlOH,N3,N2,N1))
    RFLAL2=VFLW*AZMAX1(trcsa_solml2(idsa_AlOH2,N3,N2,N1))
    RFLAL3=VFLW*AZMAX1(trcsa_solml2(idsa_AlOH3,N3,N2,N1))
    RFLAL4=VFLW*AZMAX1(trcsa_solml2(idsa_AlOH4,N3,N2,N1))
    RFLALS=VFLW*AZMAX1(trcsa_solml2(idsa_AlSO4,N3,N2,N1))
    RFLFE1=VFLW*AZMAX1(trcsa_solml2(idsa_FeOH,N3,N2,N1))
    RFLFE2=VFLW*AZMAX1(trcsa_solml2(idsa_FeOH2,N3,N2,N1))
    RFLFE3=VFLW*AZMAX1(trcsa_solml2(idsa_FeOH3,N3,N2,N1))
    RFLFE4=VFLW*AZMAX1(trcsa_solml2(idsa_FeOH4,N3,N2,N1))
    RFLFES=VFLW*AZMAX1(trcsa_solml2(idsa_FeSO4,N3,N2,N1))
    RFLCAO=VFLW*AZMAX1(trcsa_solml2(idsa_CaOH2,N3,N2,N1))
    RFLCAC=VFLW*AZMAX1(trcsa_solml2(idsa_CaCO3,N3,N2,N1))
    RFLCAH=VFLW*AZMAX1(trcsa_solml2(idsa_CaHCO3,N3,N2,N1))
    RFLCAS=VFLW*AZMAX1(trcsa_solml2(idsa_CaSO4,N3,N2,N1))
    RFLMGO=VFLW*AZMAX1(trcsa_solml2(idsa_MgOH2,N3,N2,N1))
    RFLMGC=VFLW*AZMAX1(trcsa_solml2(idsa_MgCO3,N3,N2,N1))
    RFLMGH=VFLW*AZMAX1(trcsa_solml2(idsa_MgHCO3,N3,N2,N1))
    RFLMGS=VFLW*AZMAX1(trcsa_solml2(idsa_MgSO4,N3,N2,N1))
    RFLNAC=VFLW*AZMAX1(trcsa_solml2(idsa_NaCO3,N3,N2,N1))
    RFLNAS=VFLW*AZMAX1(trcsa_solml2(idsa_NaSO4,N3,N2,N1))
    RFLKAS=VFLW*AZMAX1(trcsa_solml2(idsa_KSO4,N3,N2,N1))
    RFLH0P=VFLW*AZMAX1(trcsa_solml2(idsa_H0PO4,N3,N2,N1))*trcs_VLN(ids_H1PO4,N3,N2,N1)
    RFLH3P=VFLW*AZMAX1(trcsa_solml2(idsa_H3PO4,N3,N2,N1))*trcs_VLN(ids_H1PO4,N3,N2,N1)
    RFLF1P=VFLW*AZMAX1(trcsa_solml2(idsa_FeHPO4,N3,N2,N1))*trcs_VLN(ids_H1PO4,N3,N2,N1)
    RFLF2P=VFLW*AZMAX1(trcsa_solml2(idsa_FeH2PO4,N3,N2,N1))*trcs_VLN(ids_H1PO4,N3,N2,N1)
    RFLC0P=VFLW*AZMAX1(trcsa_solml2(idsa_CaPO4,N3,N2,N1))*trcs_VLN(ids_H1PO4,N3,N2,N1)
    RFLC1P=VFLW*AZMAX1(trcsa_solml2(idsa_CaHPO4,N3,N2,N1))*trcs_VLN(ids_H1PO4,N3,N2,N1)
    RFLC2P=VFLW*AZMAX1(trcsa_solml2(idsa_CaH2PO4,N3,N2,N1))*trcs_VLN(ids_H1PO4,N3,N2,N1)
    RFLM1P=VFLW*AZMAX1(trcsa_solml2(idsa_MgHPO4,N3,N2,N1))*trcs_VLN(ids_H1PO4,N3,N2,N1)
    RFLH0B=VFLW*AZMAX1(trcsa_solml2(idsa_H0PO4B,N3,N2,N1))*trcs_VLN(ids_H1PO4B,N3,N2,N1)
    RFLH3B=VFLW*AZMAX1(trcsa_solml2(idsa_H3PO4B,N3,N2,N1))*trcs_VLN(ids_H1PO4B,N3,N2,N1)
    RFLF1B=VFLW*AZMAX1(trcsa_solml2(idsa_FeHPO4B,N3,N2,N1))*trcs_VLN(ids_H1PO4B,N3,N2,N1)
    RFLF2B=VFLW*AZMAX1(trcsa_solml2(idsa_FeH2PO4B,N3,N2,N1))*trcs_VLN(ids_H1PO4B,N3,N2,N1)
    RFLC0B=VFLW*AZMAX1(trcsa_solml2(idsa_CaPO4B,N3,N2,N1))*trcs_VLN(ids_H1PO4B,N3,N2,N1)
    RFLC1B=VFLW*AZMAX1(trcsa_solml2(idsa_CaHPO4B,N3,N2,N1))*trcs_VLN(ids_H1PO4B,N3,N2,N1)
    RFLC2B=VFLW*AZMAX1(trcsa_solml2(idsa_CaH2PO4B,N3,N2,N1))*trcs_VLN(ids_H1PO4B,N3,N2,N1)
    RFLM1B=VFLW*AZMAX1(trcsa_solml2(idsa_MgHPO4B,N3,N2,N1))*trcs_VLN(ids_H1PO4B,N3,N2,N1)
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
    RFLAL=VFLW*AZMAX1(trcsa_solml2(idsa_Al,N6,N5,N4))
    RFLFE=VFLW*AZMAX1(trcsa_solml2(idsa_Fe,N6,N5,N4))
    RFLHY=VFLW*AZMAX1(trcsa_solml2(idsa_Hp,N6,N5,N4))
    RFLCA=VFLW*AZMAX1(trcsa_solml2(idsa_Ca,N6,N5,N4))
    RFLMG=VFLW*AZMAX1(trcsa_solml2(idsa_Mg,N6,N5,N4))
    RFLNA=VFLW*AZMAX1(trcsa_solml2(idsa_Na,N6,N5,N4))
    RFLKA=VFLW*AZMAX1(trcsa_solml2(idsa_K,N6,N5,N4))
    RFLOH=VFLW*AZMAX1(trcsa_solml2(idsa_OH,N6,N5,N4))
    RFLSO=VFLW*AZMAX1(trcsa_solml2(idsa_SO4,N6,N5,N4))
    RFLCL=VFLW*AZMAX1(trcsa_solml2(idsa_Cl,N6,N5,N4))
    RFLC3=VFLW*AZMAX1(trcsa_solml2(idsa_CO3,N6,N5,N4))
    RFLHC=VFLW*AZMAX1(trcsa_solml2(idsa_HCO3,N6,N5,N4))
    RFLAL1=VFLW*AZMAX1(trcsa_solml2(idsa_AlOH,N6,N5,N4))
    RFLAL2=VFLW*AZMAX1(trcsa_solml2(idsa_AlOH2,N6,N5,N4))
    RFLAL3=VFLW*AZMAX1(trcsa_solml2(idsa_AlOH3,N6,N5,N4))
    RFLAL4=VFLW*AZMAX1(trcsa_solml2(idsa_AlOH4,N6,N5,N4))
    RFLALS=VFLW*AZMAX1(trcsa_solml2(idsa_AlSO4,N6,N5,N4))
    RFLFE1=VFLW*AZMAX1(trcsa_solml2(idsa_FeOH,N6,N5,N4))
    RFLFE2=VFLW*AZMAX1(trcsa_solml2(idsa_FeOH2,N6,N5,N4))
    RFLFE3=VFLW*AZMAX1(trcsa_solml2(idsa_FeOH3,N6,N5,N4))
    RFLFE4=VFLW*AZMAX1(trcsa_solml2(idsa_FeOH4,N6,N5,N4))
    RFLFES=VFLW*AZMAX1(trcsa_solml2(idsa_FeSO4,N6,N5,N4))
    RFLCAO=VFLW*AZMAX1(trcsa_solml2(idsa_CaOH2,N6,N5,N4))
    RFLCAC=VFLW*AZMAX1(trcsa_solml2(idsa_CaCO3,N6,N5,N4))
    RFLCAH=VFLW*AZMAX1(trcsa_solml2(idsa_CaHCO3,N6,N5,N4))
    RFLCAS=VFLW*AZMAX1(trcsa_solml2(idsa_CaSO4,N6,N5,N4))
    RFLMGO=VFLW*AZMAX1(trcsa_solml2(idsa_MgOH2,N6,N5,N4))
    RFLMGC=VFLW*AZMAX1(trcsa_solml2(idsa_MgCO3,N6,N5,N4))
    RFLMGH=VFLW*AZMAX1(trcsa_solml2(idsa_MgHCO3,N6,N5,N4))
    RFLMGS=VFLW*AZMAX1(trcsa_solml2(idsa_MgSO4,N6,N5,N4))
    RFLNAC=VFLW*AZMAX1(trcsa_solml2(idsa_NaCO3,N6,N5,N4))
    RFLNAS=VFLW*AZMAX1(trcsa_solml2(idsa_NaSO4,N6,N5,N4))
    RFLKAS=VFLW*AZMAX1(trcsa_solml2(idsa_KSO4,N6,N5,N4))
    RFLH0P=VFLW*AZMAX1(trcsa_solml2(idsa_H0PO4,N6,N5,N4))
    RFLH3P=VFLW*AZMAX1(trcsa_solml2(idsa_H3PO4,N6,N5,N4))
    RFLF1P=VFLW*AZMAX1(trcsa_solml2(idsa_FeHPO4,N6,N5,N4))
    RFLF2P=VFLW*AZMAX1(trcsa_solml2(idsa_FeH2PO4,N6,N5,N4))
    RFLC0P=VFLW*AZMAX1(trcsa_solml2(idsa_CaPO4,N6,N5,N4))
    RFLC1P=VFLW*AZMAX1(trcsa_solml2(idsa_CaHPO4,N6,N5,N4))
    RFLC2P=VFLW*AZMAX1(trcsa_solml2(idsa_CaH2PO4,N6,N5,N4))
    RFLM1P=VFLW*AZMAX1(trcsa_solml2(idsa_MgHPO4,N6,N5,N4))
    RFLH0B=VFLW*AZMAX1(trcsa_solml2(idsa_H0PO4B,N6,N5,N4))
    RFLH3B=VFLW*AZMAX1(trcsa_solml2(idsa_H3PO4B,N6,N5,N4))
    RFLF1B=VFLW*AZMAX1(trcsa_solml2(idsa_FeHPO4B,N6,N5,N4))
    RFLF2B=VFLW*AZMAX1(trcsa_solml2(idsa_FeH2PO4B,N6,N5,N4))
    RFLC0B=VFLW*AZMAX1(trcsa_solml2(idsa_CaPO4B,N6,N5,N4))
    RFLC1B=VFLW*AZMAX1(trcsa_solml2(idsa_CaHPO4B,N6,N5,N4))
    RFLC2B=VFLW*AZMAX1(trcsa_solml2(idsa_CaH2PO4B,N6,N5,N4))
    RFLM1B=VFLW*AZMAX1(trcsa_solml2(idsa_MgHPO4B,N6,N5,N4))
  ENDIF

!     RFL*=convective flux through micropores
!     DFV*=diffusive solute flux through micropores

  trcsa_RFLS(idsa_Al,N,N6,N5,N4)=RFLAL
  trcsa_RFLS(idsa_Fe,N,N6,N5,N4)=RFLFE
  trcsa_RFLS(idsa_Hp,N,N6,N5,N4)=RFLHY
  trcsa_RFLS(idsa_Ca,N,N6,N5,N4)=RFLCA
  trcsa_RFLS(idsa_Mg,N,N6,N5,N4)=RFLMG
  trcsa_RFLS(idsa_Na,N,N6,N5,N4)=RFLNA
  trcsa_RFLS(idsa_K,N,N6,N5,N4)=RFLKA
  trcsa_RFLS(idsa_OH,N,N6,N5,N4)=RFLOH
  trcsa_RFLS(idsa_SO4,N,N6,N5,N4)=RFLSO
  trcsa_RFLS(idsa_Cl,N,N6,N5,N4)=RFLCL
  trcsa_RFLS(idsa_CO3,N,N6,N5,N4)=RFLC3
  trcsa_RFLS(idsa_HCO3,N,N6,N5,N4)=RFLHC
  trcsa_RFLS(idsa_AlOH,N,N6,N5,N4)=RFLAL1
  trcsa_RFLS(idsa_AlOH2,N,N6,N5,N4)=RFLAL2
  trcsa_RFLS(idsa_AlOH3,N,N6,N5,N4)=RFLAL3
  trcsa_RFLS(idsa_AlOH4,N,N6,N5,N4)=RFLAL4
  trcsa_RFLS(idsa_AlSO4,N,N6,N5,N4)=RFLALS
  trcsa_RFLS(idsa_FeOH,N,N6,N5,N4)=RFLFE1
  trcsa_RFLS(idsa_FeOH2,N,N6,N5,N4)=RFLFE2
  trcsa_RFLS(idsa_FeOH3,N,N6,N5,N4)=RFLFE3
  trcsa_RFLS(idsa_FeOH4,N,N6,N5,N4)=RFLFE4
  trcsa_RFLS(idsa_FeSO4,N,N6,N5,N4)=RFLFES
  trcsa_RFLS(idsa_CaOH2,N,N6,N5,N4)=RFLCAO
  trcsa_RFLS(idsa_CaCO3,N,N6,N5,N4)=RFLCAC
  trcsa_RFLS(idsa_CaHCO3,N,N6,N5,N4)=RFLCAH
  trcsa_RFLS(idsa_CaSO4,N,N6,N5,N4)=RFLCAS
  trcsa_RFLS(idsa_MgOH2,N,N6,N5,N4)=RFLMGO
  trcsa_RFLS(idsa_MgCO3,N,N6,N5,N4)=RFLMGC
  trcsa_RFLS(idsa_MgHCO3,N,N6,N5,N4)=RFLMGH
  trcsa_RFLS(idsa_MgSO4,N,N6,N5,N4)=RFLMGS
  trcsa_RFLS(idsa_NaCO3,N,N6,N5,N4)=RFLNAC
  trcsa_RFLS(idsa_NaSO4,N,N6,N5,N4)=RFLNAS
  trcsa_RFLS(idsa_KSO4,N,N6,N5,N4)=RFLKAS
  trcsa_RFLS(idsa_H0PO4,N,N6,N5,N4)=RFLH0P
  trcsa_RFLS(idsa_H3PO4,N,N6,N5,N4)=RFLH3P
  trcsa_RFLS(idsa_FeHPO4,N,N6,N5,N4)=RFLF1P
  trcsa_RFLS(idsa_FeH2PO4,N,N6,N5,N4)=RFLF2P
  trcsa_RFLS(idsa_CaPO4,N,N6,N5,N4)=RFLC0P
  trcsa_RFLS(idsa_CaHPO4,N,N6,N5,N4)=RFLC1P
  trcsa_RFLS(idsa_CaH2PO4,N,N6,N5,N4)=RFLC2P
  trcsa_RFLS(idsa_MgHPO4,N,N6,N5,N4)=RFLM1P
  trcsa_RFLS(idsa_H0PO4B,N,N6,N5,N4)=RFLH0B
  trcsa_RFLS(idsa_H3PO4B,N,N6,N5,N4)=RFLH3B
  trcsa_RFLS(idsa_FeHPO4B,N,N6,N5,N4)=RFLF1B
  trcsa_RFLS(idsa_FeH2PO4B,N,N6,N5,N4)=RFLF2B
  trcsa_RFLS(idsa_CaPO4B,N,N6,N5,N4)=RFLC0B
  trcsa_RFLS(idsa_CaHPO4B,N,N6,N5,N4)=RFLC1B
  trcsa_RFLS(idsa_CaH2PO4B,N,N6,N5,N4)=RFLC2B
  trcsa_RFLS(idsa_MgHPO4B,N,N6,N5,N4)=RFLM1B
  end subroutine SoluteAdvMicropore
!------------------------------------------------------------------------------------------
  subroutine SoluteDifsMicropore(M,N,N1,N2,N3,N4,N5,N6,THETW1)
  implicit none
  integer, intent(in) :: M,N,N1,N2,N3,N4,N5,N6
  REAL(R8), INTENT(IN):: THETW1(JZ,JY,JX)
  REAL(R8) :: VLWPA1,VLWPB1,VLWPA2,VLWPB2
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
    CAL1=AZMAX1(trcsa_solml2(idsa_Al,N3,N2,N1)/VOLWM(M,N3,N2,N1))
    CFE1=AZMAX1(trcsa_solml2(idsa_Fe,N3,N2,N1)/VOLWM(M,N3,N2,N1))
    CHY1=AZMAX1(trcsa_solml2(idsa_Hp,N3,N2,N1)/VOLWM(M,N3,N2,N1))
    CCA1=AZMAX1(trcsa_solml2(idsa_Ca,N3,N2,N1)/VOLWM(M,N3,N2,N1))
    CMG1=AZMAX1(trcsa_solml2(idsa_Mg,N3,N2,N1)/VOLWM(M,N3,N2,N1))
    CNA1=AZMAX1(trcsa_solml2(idsa_Na,N3,N2,N1)/VOLWM(M,N3,N2,N1))
    CKA1=AZMAX1(trcsa_solml2(idsa_K,N3,N2,N1)/VOLWM(M,N3,N2,N1))
    COH1=AZMAX1(trcsa_solml2(idsa_OH,N3,N2,N1)/VOLWM(M,N3,N2,N1))
    CSO41=AZMAX1(trcsa_solml2(idsa_SO4,N3,N2,N1)/VOLWM(M,N3,N2,N1))
    CCL1=AZMAX1(trcsa_solml2(idsa_Cl,N3,N2,N1)/VOLWM(M,N3,N2,N1))
    CCO31=AZMAX1(trcsa_solml2(idsa_CO3,N3,N2,N1)/VOLWM(M,N3,N2,N1))
    CHCO31=AZMAX1(trcsa_solml2(idsa_HCO3,N3,N2,N1)/VOLWM(M,N3,N2,N1))
    CAL11=AZMAX1(trcsa_solml2(idsa_AlOH,N3,N2,N1)/VOLWM(M,N3,N2,N1))
    CAL21=AZMAX1(trcsa_solml2(idsa_AlOH2,N3,N2,N1)/VOLWM(M,N3,N2,N1))
    CAL31=AZMAX1(trcsa_solml2(idsa_AlOH3,N3,N2,N1)/VOLWM(M,N3,N2,N1))
    CAL41=AZMAX1(trcsa_solml2(idsa_AlOH4,N3,N2,N1)/VOLWM(M,N3,N2,N1))
    CALS1=AZMAX1(trcsa_solml2(idsa_AlSO4,N3,N2,N1)/VOLWM(M,N3,N2,N1))
    CFE11=AZMAX1(trcsa_solml2(idsa_FeOH,N3,N2,N1)/VOLWM(M,N3,N2,N1))
    CFE21=AZMAX1(trcsa_solml2(idsa_FeOH2,N3,N2,N1)/VOLWM(M,N3,N2,N1))
    CFE31=AZMAX1(trcsa_solml2(idsa_FeOH3,N3,N2,N1)/VOLWM(M,N3,N2,N1))
    CFE41=AZMAX1(trcsa_solml2(idsa_FeOH4,N3,N2,N1)/VOLWM(M,N3,N2,N1))
    CFES1=AZMAX1(trcsa_solml2(idsa_FeSO4,N3,N2,N1)/VOLWM(M,N3,N2,N1))
    CCAO1=AZMAX1(trcsa_solml2(idsa_CaOH2,N3,N2,N1)/VOLWM(M,N3,N2,N1))
    CCAC1=AZMAX1(trcsa_solml2(idsa_CaCO3,N3,N2,N1)/VOLWM(M,N3,N2,N1))
    CCAH1=AZMAX1(trcsa_solml2(idsa_CaHCO3,N3,N2,N1)/VOLWM(M,N3,N2,N1))
    CCAS1=AZMAX1(trcsa_solml2(idsa_CaSO4,N3,N2,N1)/VOLWM(M,N3,N2,N1))
    CMGO1=AZMAX1(trcsa_solml2(idsa_MgOH2,N3,N2,N1)/VOLWM(M,N3,N2,N1))
    CMGC1=AZMAX1(trcsa_solml2(idsa_MgCO3,N3,N2,N1)/VOLWM(M,N3,N2,N1))
    CMGH1=AZMAX1(trcsa_solml2(idsa_MgHCO3,N3,N2,N1)/VOLWM(M,N3,N2,N1))
    CMGS1=AZMAX1(trcsa_solml2(idsa_MgSO4,N3,N2,N1)/VOLWM(M,N3,N2,N1))
    CNAC1=AZMAX1(trcsa_solml2(idsa_NaCO3,N3,N2,N1)/VOLWM(M,N3,N2,N1))
    CNAS1=AZMAX1(trcsa_solml2(idsa_NaSO4,N3,N2,N1)/VOLWM(M,N3,N2,N1))
    CKAS1=AZMAX1(trcsa_solml2(idsa_KSO4,N3,N2,N1)/VOLWM(M,N3,N2,N1))
    IF(VLWPA1.GT.ZEROS(N2,N1))THEN
      C0PO41=AZMAX1(trcsa_solml2(idsa_H0PO4,N3,N2,N1)/VLWPA1)
      C3PO41=AZMAX1(trcsa_solml2(idsa_H3PO4,N3,N2,N1)/VLWPA1)
      CFE1P1=AZMAX1(trcsa_solml2(idsa_FeHPO4,N3,N2,N1)/VLWPA1)
      CFE2P1=AZMAX1(trcsa_solml2(idsa_FeH2PO4,N3,N2,N1)/VLWPA1)
      CCA0P1=AZMAX1(trcsa_solml2(idsa_CaPO4,N3,N2,N1)/VLWPA1)
      CCA1P1=AZMAX1(trcsa_solml2(idsa_CaHPO4,N3,N2,N1)/VLWPA1)
      CCA2P1=AZMAX1(trcsa_solml2(idsa_CaH2PO4,N3,N2,N1)/VLWPA1)
      CMG1P1=AZMAX1(trcsa_solml2(idsa_MgHPO4,N3,N2,N1)/VLWPA1)
    ELSE
      C0PO41=0.0_r8
      C3PO41=0.0_r8
      CFE1P1=0.0_r8
      CFE2P1=0.0_r8
      CCA0P1=0.0_r8
      CCA1P1=0.0_r8
      CCA2P1=0.0_r8
      CMG1P1=0.0_r8
    ENDIF
    IF(VLWPB1.GT.ZEROS(N2,N1))THEN
      C0POB1=AZMAX1(trcsa_solml2(idsa_H0PO4B,N3,N2,N1)/VLWPB1)
      C3POB1=AZMAX1(trcsa_solml2(idsa_H3PO4B,N3,N2,N1)/VLWPB1)
      CF1PB1=AZMAX1(trcsa_solml2(idsa_FeHPO4B,N3,N2,N1)/VLWPB1)
      CF2PB1=AZMAX1(trcsa_solml2(idsa_FeH2PO4B,N3,N2,N1)/VLWPB1)
      CC0PB1=AZMAX1(trcsa_solml2(idsa_CaPO4B,N3,N2,N1)/VLWPB1)
      CC1PB1=AZMAX1(trcsa_solml2(idsa_CaHPO4B,N3,N2,N1)/VLWPB1)
      CC2PB1=AZMAX1(trcsa_solml2(idsa_CaH2PO4B,N3,N2,N1)/VLWPB1)
      CM1PB1=AZMAX1(trcsa_solml2(idsa_MgHPO4B,N3,N2,N1)/VLWPB1)
    ELSE
      C0POB1=C0PO41
      C3POB1=C3PO41
      CF1PB1=CFE1P1
      CF2PB1=CFE2P1
      CC0PB1=CCA0P1
      CC1PB1=CCA1P1
      CC2PB1=CCA2P1
      CM1PB1=CMG1P1
    ENDIF
    CAL2=AZMAX1(trcsa_solml2(idsa_Al,N6,N5,N4)/VOLWM(M,N6,N5,N4))
    CFE2=AZMAX1(trcsa_solml2(idsa_Fe,N6,N5,N4)/VOLWM(M,N6,N5,N4))
    CHY2=AZMAX1(trcsa_solml2(idsa_Hp,N6,N5,N4)/VOLWM(M,N6,N5,N4))
    CCA2=AZMAX1(trcsa_solml2(idsa_Ca,N6,N5,N4)/VOLWM(M,N6,N5,N4))
    CMG2=AZMAX1(trcsa_solml2(idsa_Mg,N6,N5,N4)/VOLWM(M,N6,N5,N4))
    CNA2=AZMAX1(trcsa_solml2(idsa_Na,N6,N5,N4)/VOLWM(M,N6,N5,N4))
    CKA2=AZMAX1(trcsa_solml2(idsa_K,N6,N5,N4)/VOLWM(M,N6,N5,N4))
    COH2=AZMAX1(trcsa_solml2(idsa_OH,N6,N5,N4)/VOLWM(M,N6,N5,N4))
    CSO42=AZMAX1(trcsa_solml2(idsa_SO4,N6,N5,N4)/VOLWM(M,N6,N5,N4))
    CCL2=AZMAX1(trcsa_solml2(idsa_Cl,N6,N5,N4)/VOLWM(M,N6,N5,N4))
    CCO32=AZMAX1(trcsa_solml2(idsa_CO3,N6,N5,N4)/VOLWM(M,N6,N5,N4))
    CHCO32=AZMAX1(trcsa_solml2(idsa_HCO3,N6,N5,N4)/VOLWM(M,N6,N5,N4))
    CAL12=AZMAX1(trcsa_solml2(idsa_AlOH,N6,N5,N4)/VOLWM(M,N6,N5,N4))
    CAL22=AZMAX1(trcsa_solml2(idsa_AlOH2,N6,N5,N4)/VOLWM(M,N6,N5,N4))
    CAL32=AZMAX1(trcsa_solml2(idsa_AlOH3,N6,N5,N4)/VOLWM(M,N6,N5,N4))
    CAL42=AZMAX1(trcsa_solml2(idsa_AlOH4,N6,N5,N4)/VOLWM(M,N6,N5,N4))
    CALS2=AZMAX1(trcsa_solml2(idsa_AlSO4,N6,N5,N4)/VOLWM(M,N6,N5,N4))
    CFE12=AZMAX1(trcsa_solml2(idsa_FeOH,N6,N5,N4)/VOLWM(M,N6,N5,N4))
    CFE22=AZMAX1(trcsa_solml2(idsa_FeOH2,N6,N5,N4)/VOLWM(M,N6,N5,N4))
    CFE32=AZMAX1(trcsa_solml2(idsa_FeOH3,N6,N5,N4)/VOLWM(M,N6,N5,N4))
    CFE42=AZMAX1(trcsa_solml2(idsa_FeOH4,N6,N5,N4)/VOLWM(M,N6,N5,N4))
    CFES2=AZMAX1(trcsa_solml2(idsa_FeSO4,N6,N5,N4)/VOLWM(M,N6,N5,N4))
    CCAO2=AZMAX1(trcsa_solml2(idsa_CaOH2,N6,N5,N4)/VOLWM(M,N6,N5,N4))
    CCAC2=AZMAX1(trcsa_solml2(idsa_CaCO3,N6,N5,N4)/VOLWM(M,N6,N5,N4))
    CCAH2=AZMAX1(trcsa_solml2(idsa_CaHCO3,N6,N5,N4)/VOLWM(M,N6,N5,N4))
    CCAS2=AZMAX1(trcsa_solml2(idsa_CaSO4,N6,N5,N4)/VOLWM(M,N6,N5,N4))
    CMGO2=AZMAX1(trcsa_solml2(idsa_MgOH2,N6,N5,N4)/VOLWM(M,N6,N5,N4))
    CMGC2=AZMAX1(trcsa_solml2(idsa_MgCO3,N6,N5,N4)/VOLWM(M,N6,N5,N4))
    CMGH2=AZMAX1(trcsa_solml2(idsa_MgHCO3,N6,N5,N4)/VOLWM(M,N6,N5,N4))
    CMGS2=AZMAX1(trcsa_solml2(idsa_MgSO4,N6,N5,N4)/VOLWM(M,N6,N5,N4))
    CNAC2=AZMAX1(trcsa_solml2(idsa_NaCO3,N6,N5,N4)/VOLWM(M,N6,N5,N4))
    CNAS2=AZMAX1(trcsa_solml2(idsa_NaSO4,N6,N5,N4)/VOLWM(M,N6,N5,N4))
    CKAS2=AZMAX1(trcsa_solml2(idsa_KSO4,N6,N5,N4)/VOLWM(M,N6,N5,N4))
    IF(VLWPA2.GT.ZEROS(N5,N4))THEN
      C0PO42=AZMAX1(trcsa_solml2(idsa_H0PO4,N6,N5,N4)/VOLWM(M,N6,N5,N4))
      C3PO42=AZMAX1(trcsa_solml2(idsa_H3PO4,N6,N5,N4)/VOLWM(M,N6,N5,N4))
      CFE1P2=AZMAX1(trcsa_solml2(idsa_FeHPO4,N6,N5,N4)/VOLWM(M,N6,N5,N4))
      CFE2P2=AZMAX1(trcsa_solml2(idsa_FeH2PO4,N6,N5,N4)/VOLWM(M,N6,N5,N4))
      CCA0P2=AZMAX1(trcsa_solml2(idsa_CaPO4,N6,N5,N4)/VOLWM(M,N6,N5,N4))
      CCA1P2=AZMAX1(trcsa_solml2(idsa_CaHPO4,N6,N5,N4)/VOLWM(M,N6,N5,N4))
      CCA2P2=AZMAX1(trcsa_solml2(idsa_CaH2PO4,N6,N5,N4)/VOLWM(M,N6,N5,N4))
      CMG1P2=AZMAX1(trcsa_solml2(idsa_MgHPO4,N6,N5,N4)/VOLWM(M,N6,N5,N4))
    ELSE
      C0PO42=0.0_r8
      C3PO42=0.0_r8
      CFE1P2=0.0_r8
      CFE2P2=0.0_r8
      CCA0P2=0.0_r8
      CCA1P2=0.0_r8
      CCA2P2=0.0_r8
      CMG1P2=0.0_r8
    ENDIF
    IF(VLWPB2.GT.ZEROS(N5,N4))THEN
      C0POB2=AZMAX1(trcsa_solml2(idsa_H0PO4B,N6,N5,N4)/VOLWM(M,N6,N5,N4))
      C3POB2=AZMAX1(trcsa_solml2(idsa_H3PO4B,N6,N5,N4)/VOLWM(M,N6,N5,N4))
      CF1PB2=AZMAX1(trcsa_solml2(idsa_FeHPO4B,N6,N5,N4)/VOLWM(M,N6,N5,N4))
      CF2PB2=AZMAX1(trcsa_solml2(idsa_FeH2PO4B,N6,N5,N4)/VOLWM(M,N6,N5,N4))
      CC0PB2=AZMAX1(trcsa_solml2(idsa_CaPO4B,N6,N5,N4)/VOLWM(M,N6,N5,N4))
      CC1PB2=AZMAX1(trcsa_solml2(idsa_CaHPO4B,N6,N5,N4)/VOLWM(M,N6,N5,N4))
      CC2PB2=AZMAX1(trcsa_solml2(idsa_CaH2PO4B,N6,N5,N4)/VOLWM(M,N6,N5,N4))
      CM1PB2=AZMAX1(trcsa_solml2(idsa_MgHPO4B,N6,N5,N4)/VOLWM(M,N6,N5,N4))
    ELSE
      C0POB2=C0PO42
      C3POB2=C3PO42
      CF1PB2=CFE1P2
      CF2PB2=CFE2P2
      CC0PB2=CCA0P2
      CC1PB2=CCA1P2
      CC2PB2=CCA2P2
      CM1PB2=CMG1P2
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
    DIFFE=(FESGL2(N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)
    DIFHY=(HYSGL2(N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)
    DIFCA=(CASGL2(N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)
    DIFMG=(GMSGL2(N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)
    DIFNA=(ANSGL2(N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)
    DIFKA=(AKSGL2(N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)
    DIFOH=(OHSGL2(N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)
    DIFSO=(SOSGL2(N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)
    DIFCL=(CLSXL2(N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)
    DIFC3=(C3SGL2(N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)
    DIFHC=(HCSGL2(N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)
!
!     DIFFUSIVE FLUXES BETWEEN CURRENT AND ADJACENT GRID CELL
!     MICROPORES
!
    DFVAL=DIFAL*(CAL1-CAL2)
    DFVFE=DIFFE*(CFE1-CFE2)
    DFVHY=DIFHY*(CHY1-CHY2)
    DFVCA=DIFCA*(CCA1-CCA2)
    DFVMG=DIFMG*(CMG1-CMG2)
    DFVNA=DIFNA*(CNA1-CNA2)
    DFVKA=DIFKA*(CKA1-CKA2)
    DFVOH=DIFOH*(COH1-COH2)
    DFVSO=DIFSO*(CSO41-CSO42)
    DFVCL=DIFCL*(CCL1-CCL2)
    DFVC3=DIFC3*(CCO31-CCO32)
    DFVHC=DIFHC*(CHCO31-CHCO32)
    DFVAL1=DIFAL*(CAL11-CAL12)
    DFVAL2=DIFAL*(CAL21-CAL22)
    DFVAL3=DIFAL*(CAL31-CAL32)
    DFVAL4=DIFAL*(CAL41-CAL42)
    DFVALS=DIFAL*(CALS1-CALS2)
    DFVFE1=DIFFE*(CFE11-CFE12)
    DFVFE2=DIFFE*(CFE21-CFE22)
    DFVFE3=DIFFE*(CFE31-CFE32)
    DFVFE4=DIFFE*(CFE41-CFE42)
    DFVFES=DIFFE*(CFES1-CFES2)
    DFVCAO=DIFCA*(CCAO1-CCAO2)
    DFVCAC=DIFCA*(CCAC1-CCAC2)
    DFVCAH=DIFCA*(CCAH1-CCAH2)
    DFVCAS=DIFCA*(CCAS1-CCAS2)
    DFVMGO=DIFMG*(CMGO1-CMGO2)
    DFVMGC=DIFMG*(CMGC1-CMGC2)
    DFVMGH=DIFMG*(CMGH1-CMGH2)
    DFVMGS=DIFMG*(CMGS1-CMGS2)
    DFVNAC=DIFNA*(CNAC1-CNAC2)
    DFVNAS=DIFNA*(CNAS1-CNAS2)
    DFVKAS=DIFKA*(CKAS1-CKAS2)
    DFVH0P=DIFPO*(C0PO41-C0PO42)*trcs_VLN(ids_H1PO4,N6,N5,N4)
    DFVH3P=DIFPO*(C3PO41-C3PO42)*trcs_VLN(ids_H1PO4,N6,N5,N4)
    DFVF1P=DIFPO*(CFE1P1-CFE1P2)*trcs_VLN(ids_H1PO4,N6,N5,N4)
    DFVF2P=DIFPO*(CFE2P1-CFE2P2)*trcs_VLN(ids_H1PO4,N6,N5,N4)
    DFVC0P=DIFPO*(CCA0P1-CCA0P2)*trcs_VLN(ids_H1PO4,N6,N5,N4)
    DFVC1P=DIFPO*(CCA1P1-CCA1P2)*trcs_VLN(ids_H1PO4,N6,N5,N4)
    DFVC2P=DIFPO*(CCA2P1-CCA2P2)*trcs_VLN(ids_H1PO4,N6,N5,N4)
    DFVM1P=DIFPO*(CMG1P1-CMG1P2)*trcs_VLN(ids_H1PO4,N6,N5,N4)
    DFVH0B=DIFPO*(C0POB1-C0POB2)*trcs_VLN(ids_H1PO4B,N6,N5,N4)
    DFVH3B=DIFPO*(C3POB1-C3POB2)*trcs_VLN(ids_H1PO4B,N6,N5,N4)
    DFVF1B=DIFPO*(CF1PB1-CF1PB2)*trcs_VLN(ids_H1PO4B,N6,N5,N4)
    DFVF2B=DIFPO*(CF2PB1-CF2PB2)*trcs_VLN(ids_H1PO4B,N6,N5,N4)
    DFVC0B=DIFPO*(CC0PB1-CC0PB2)*trcs_VLN(ids_H1PO4B,N6,N5,N4)
    DFVC1B=DIFPO*(CC1PB1-CC1PB2)*trcs_VLN(ids_H1PO4B,N6,N5,N4)
    DFVC2B=DIFPO*(CC2PB1-CC2PB2)*trcs_VLN(ids_H1PO4B,N6,N5,N4)
    DFVM1B=DIFPO*(CM1PB1-CM1PB2)*trcs_VLN(ids_H1PO4B,N6,N5,N4)
  ELSE
    DFVAL=0.0_r8
    DFVFE=0.0_r8
    DFVHY=0.0_r8
    DFVCA=0.0_r8
    DFVMG=0.0_r8
    DFVNA=0.0_r8
    DFVKA=0.0_r8
    DFVOH=0.0_r8
    DFVSO=0.0_r8
    DFVCL=0.0_r8
    DFVC3=0.0_r8
    DFVHC=0.0_r8
    DFVAL1=0.0_r8
    DFVAL2=0.0_r8
    DFVAL3=0.0_r8
    DFVAL4=0.0_r8
    DFVALS=0.0_r8
    DFVFE1=0.0_r8
    DFVFE2=0.0_r8
    DFVFE3=0.0_r8
    DFVFE4=0.0_r8
    DFVFES=0.0_r8
    DFVCAO=0.0_r8
    DFVCAC=0.0_r8
    DFVCAH=0.0_r8
    DFVCAS=0.0_r8
    DFVMGO=0.0_r8
    DFVMGC=0.0_r8
    DFVMGH=0.0_r8
    DFVMGS=0.0_r8
    DFVNAC=0.0_r8
    DFVNAS=0.0_r8
    DFVKAS=0.0_r8
    DFVH0P=0.0_r8
    DFVH3P=0.0_r8
    DFVF1P=0.0_r8
    DFVF2P=0.0_r8
    DFVC0P=0.0_r8
    DFVC1P=0.0_r8
    DFVC2P=0.0_r8
    DFVM1P=0.0_r8
    DFVH0B=0.0_r8
    DFVH3B=0.0_r8
    DFVF1B=0.0_r8
    DFVF2B=0.0_r8
    DFVC0B=0.0_r8
    DFVC1B=0.0_r8
    DFVC2B=0.0_r8
    DFVM1B=0.0_r8
  ENDIF
!     RFL*=convective flux through micropores
!     DFV*=diffusive solute flux through micropores

  trcsa_RFLS(idsa_Al,N,N6,N5,N4)=trcsa_RFLS(idsa_Al,N,N6,N5,N4)+DFVAL
  trcsa_RFLS(idsa_Fe,N,N6,N5,N4)=trcsa_RFLS(idsa_Fe,N,N6,N5,N4)+DFVFE
  trcsa_RFLS(idsa_Hp,N,N6,N5,N4)=trcsa_RFLS(idsa_Hp,N,N6,N5,N4)+DFVHY
  trcsa_RFLS(idsa_Ca,N,N6,N5,N4)=trcsa_RFLS(idsa_Ca,N,N6,N5,N4)+DFVCA
  trcsa_RFLS(idsa_Mg,N,N6,N5,N4)=trcsa_RFLS(idsa_Mg,N,N6,N5,N4)+DFVMG
  trcsa_RFLS(idsa_Na,N,N6,N5,N4)=trcsa_RFLS(idsa_Na,N,N6,N5,N4)+DFVNA
  trcsa_RFLS(idsa_K,N,N6,N5,N4)=trcsa_RFLS(idsa_K,N,N6,N5,N4)+DFVKA
  trcsa_RFLS(idsa_OH,N,N6,N5,N4)=trcsa_RFLS(idsa_OH,N,N6,N5,N4)+DFVOH
  trcsa_RFLS(idsa_SO4,N,N6,N5,N4)=trcsa_RFLS(idsa_SO4,N,N6,N5,N4)+DFVSO
  trcsa_RFLS(idsa_Cl,N,N6,N5,N4)=trcsa_RFLS(idsa_Cl,N,N6,N5,N4)+DFVCL
  trcsa_RFLS(idsa_CO3,N,N6,N5,N4)=trcsa_RFLS(idsa_CO3,N,N6,N5,N4)+DFVC3
  trcsa_RFLS(idsa_HCO3,N,N6,N5,N4)=trcsa_RFLS(idsa_HCO3,N,N6,N5,N4)+DFVHC
  trcsa_RFLS(idsa_AlOH,N,N6,N5,N4)=trcsa_RFLS(idsa_AlOH,N,N6,N5,N4)+DFVAL1
  trcsa_RFLS(idsa_AlOH2,N,N6,N5,N4)=trcsa_RFLS(idsa_AlOH2,N,N6,N5,N4)+DFVAL2
  trcsa_RFLS(idsa_AlOH3,N,N6,N5,N4)=trcsa_RFLS(idsa_AlOH3,N,N6,N5,N4)+DFVAL3
  trcsa_RFLS(idsa_AlOH4,N,N6,N5,N4)=trcsa_RFLS(idsa_AlOH4,N,N6,N5,N4)+DFVAL4
  trcsa_RFLS(idsa_AlSO4,N,N6,N5,N4)=trcsa_RFLS(idsa_AlSO4,N,N6,N5,N4)+DFVALS
  trcsa_RFLS(idsa_FeOH,N,N6,N5,N4)=trcsa_RFLS(idsa_FeOH,N,N6,N5,N4)+DFVFE1
  trcsa_RFLS(idsa_FeOH2,N,N6,N5,N4)=trcsa_RFLS(idsa_FeOH2,N,N6,N5,N4)+DFVFE2
  trcsa_RFLS(idsa_FeOH3,N,N6,N5,N4)=trcsa_RFLS(idsa_FeOH3,N,N6,N5,N4)+DFVFE3
  trcsa_RFLS(idsa_FeOH4,N,N6,N5,N4)=trcsa_RFLS(idsa_FeOH4,N,N6,N5,N4)+DFVFE4
  trcsa_RFLS(idsa_FeSO4,N,N6,N5,N4)=trcsa_RFLS(idsa_FeSO4,N,N6,N5,N4)+DFVFES
  trcsa_RFLS(idsa_CaOH2,N,N6,N5,N4)=trcsa_RFLS(idsa_CaOH2,N,N6,N5,N4)+DFVCAO
  trcsa_RFLS(idsa_CaCO3,N,N6,N5,N4)=trcsa_RFLS(idsa_CaCO3,N,N6,N5,N4)+DFVCAC
  trcsa_RFLS(idsa_CaHCO3,N,N6,N5,N4)=trcsa_RFLS(idsa_CaHCO3,N,N6,N5,N4)+DFVCAH
  trcsa_RFLS(idsa_CaSO4,N,N6,N5,N4)=trcsa_RFLS(idsa_CaSO4,N,N6,N5,N4)+DFVCAS
  trcsa_RFLS(idsa_MgOH2,N,N6,N5,N4)=trcsa_RFLS(idsa_MgOH2,N,N6,N5,N4)+DFVMGO
  trcsa_RFLS(idsa_MgCO3,N,N6,N5,N4)=trcsa_RFLS(idsa_MgCO3,N,N6,N5,N4)+DFVMGC
  trcsa_RFLS(idsa_MgHCO3,N,N6,N5,N4)=trcsa_RFLS(idsa_MgHCO3,N,N6,N5,N4)+DFVMGH
  trcsa_RFLS(idsa_MgSO4,N,N6,N5,N4)=trcsa_RFLS(idsa_MgSO4,N,N6,N5,N4)+DFVMGS
  trcsa_RFLS(idsa_NaCO3,N,N6,N5,N4)=trcsa_RFLS(idsa_NaCO3,N,N6,N5,N4)+DFVNAC
  trcsa_RFLS(idsa_NaSO4,N,N6,N5,N4)=trcsa_RFLS(idsa_NaSO4,N,N6,N5,N4)+DFVNAS
  trcsa_RFLS(idsa_KSO4,N,N6,N5,N4)=trcsa_RFLS(idsa_KSO4,N,N6,N5,N4)+DFVKAS
  trcsa_RFLS(idsa_H0PO4,N,N6,N5,N4)=trcsa_RFLS(idsa_H0PO4,N,N6,N5,N4)+DFVH0P
  trcsa_RFLS(idsa_H3PO4,N,N6,N5,N4)=trcsa_RFLS(idsa_H3PO4,N,N6,N5,N4)+DFVH3P
  trcsa_RFLS(idsa_FeHPO4,N,N6,N5,N4)=trcsa_RFLS(idsa_FeHPO4,N,N6,N5,N4)+DFVF1P
  trcsa_RFLS(idsa_FeH2PO4,N,N6,N5,N4)=trcsa_RFLS(idsa_FeH2PO4,N,N6,N5,N4)+DFVF2P
  trcsa_RFLS(idsa_CaPO4,N,N6,N5,N4)=trcsa_RFLS(idsa_CaPO4,N,N6,N5,N4)+DFVC0P
  trcsa_RFLS(idsa_CaHPO4,N,N6,N5,N4)=trcsa_RFLS(idsa_CaHPO4,N,N6,N5,N4)+DFVC1P
  trcsa_RFLS(idsa_CaH2PO4,N,N6,N5,N4)=trcsa_RFLS(idsa_CaH2PO4,N,N6,N5,N4)+DFVC2P
  trcsa_RFLS(idsa_MgHPO4,N,N6,N5,N4)=trcsa_RFLS(idsa_MgHPO4,N,N6,N5,N4)+DFVM1P
  trcsa_RFLS(idsa_H0PO4B,N,N6,N5,N4)=trcsa_RFLS(idsa_H0PO4B,N,N6,N5,N4)+DFVH0B
  trcsa_RFLS(idsa_H3PO4B,N,N6,N5,N4)=trcsa_RFLS(idsa_H3PO4B,N,N6,N5,N4)+DFVH3B
  trcsa_RFLS(idsa_FeHPO4B,N,N6,N5,N4)=trcsa_RFLS(idsa_FeHPO4B,N,N6,N5,N4)+DFVF1B
  trcsa_RFLS(idsa_FeH2PO4B,N,N6,N5,N4)=trcsa_RFLS(idsa_FeH2PO4B,N,N6,N5,N4)+DFVF2B
  trcsa_RFLS(idsa_CaPO4B,N,N6,N5,N4)=trcsa_RFLS(idsa_CaPO4B,N,N6,N5,N4)+DFVC0B
  trcsa_RFLS(idsa_CaHPO4B,N,N6,N5,N4)=trcsa_RFLS(idsa_CaHPO4B,N,N6,N5,N4)+DFVC1B
  trcsa_RFLS(idsa_CaH2PO4B,N,N6,N5,N4)=trcsa_RFLS(idsa_CaH2PO4B,N,N6,N5,N4)+DFVC2B
  trcsa_RFLS(idsa_MgHPO4B,N,N6,N5,N4)=trcsa_RFLS(idsa_MgHPO4B,N,N6,N5,N4)+DFVM1B

  end subroutine SoluteDifsMicropore
!----------------------------------------------------------------------
  subroutine SoluteAdvMacropore(M,N,N1,N2,N3,N4,N5,N6)
  implicit none
  integer, intent(in) :: M,N,N1,N2,N3,N4,N5,N6

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
      RFHAL=VFLW*AZMAX1((trcsa_soHml2(idsa_Al,N3,N2,N1)-AZMIN1(trcsa_RFXS(idsa_Al,NU(N2,N1),N2,N1))))
      RFHFE=VFLW*AZMAX1((trcsa_soHml2(idsa_Fe,N3,N2,N1)-AZMIN1(trcsa_RFXS(idsa_Fe,NU(N2,N1),N2,N1))))
      RFHHY=VFLW*AZMAX1((trcsa_soHml2(idsa_Hp,N3,N2,N1)-AZMIN1(trcsa_RFXS(idsa_Hp,NU(N2,N1),N2,N1))))
      RFHCA=VFLW*AZMAX1((trcsa_soHml2(idsa_Ca,N3,N2,N1)-AZMIN1(trcsa_RFXS(idsa_Ca,NU(N2,N1),N2,N1))))
      RFHMG=VFLW*AZMAX1((trcsa_soHml2(idsa_Mg,N3,N2,N1)-AZMIN1(trcsa_RFXS(idsa_Mg,NU(N2,N1),N2,N1))))
      RFHNA=VFLW*AZMAX1((trcsa_soHml2(idsa_Na,N3,N2,N1)-AZMIN1(trcsa_RFXS(idsa_Na,NU(N2,N1),N2,N1))))
      RFHKA=VFLW*AZMAX1((trcsa_soHml2(idsa_K,N3,N2,N1)-AZMIN1(trcsa_RFXS(idsa_K,NU(N2,N1),N2,N1))))
      RFHOH=VFLW*AZMAX1((trcsa_soHml2(idsa_OH,N3,N2,N1)-AZMIN1(trcsa_RFXS(idsa_OH,NU(N2,N1),N2,N1))))
      RFHSO=VFLW*AZMAX1((trcsa_soHml2(idsa_SO4,N3,N2,N1)-AZMIN1(trcsa_RFXS(idsa_SO4,NU(N2,N1),N2,N1))))
      RFHCL=VFLW*AZMAX1((trcsa_soHml2(idsa_Cl,N3,N2,N1)-AZMIN1(trcsa_RFXS(idsa_Cl,NU(N2,N1),N2,N1))))
      RFHC3=VFLW*AZMAX1((trcsa_soHml2(idsa_CO3,N3,N2,N1)-AZMIN1(trcsa_RFXS(idsa_CO3,NU(N2,N1),N2,N1))))
      RFHHC=VFLW*AZMAX1((trcsa_soHml2(idsa_HCO3,N3,N2,N1)-AZMIN1(trcsa_RFXS(idsa_HCO3,NU(N2,N1),N2,N1))))
      RFHAL1=VFLW*AZMAX1((trcsa_soHml2(idsa_AlOH,N3,N2,N1)-AZMIN1(trcsa_RFXS(idsa_AlOH,NU(N2,N1),N2,N1))))
      RFHAL2=VFLW*AZMAX1((trcsa_soHml2(idsa_AlOH2,N3,N2,N1)-AZMIN1(trcsa_RFXS(idsa_AlOH2,NU(N2,N1),N2,N1))))
      RFHAL3=VFLW*AZMAX1((trcsa_soHml2(idsa_AlOH3,N3,N2,N1)-AZMIN1(trcsa_RFXS(idsa_AlOH3,NU(N2,N1),N2,N1))))
      RFHAL4=VFLW*AZMAX1((trcsa_soHml2(idsa_AlOH4,N3,N2,N1)-AZMIN1(trcsa_RFXS(idsa_AlOH4,NU(N2,N1),N2,N1))))
      RFHALS=VFLW*AZMAX1((trcsa_soHml2(idsa_AlSO4,N3,N2,N1)-AZMIN1(trcsa_RFXS(idsa_AlSO4,NU(N2,N1),N2,N1))))
      RFHFE1=VFLW*AZMAX1((trcsa_soHml2(idsa_FeOH,N3,N2,N1)-AZMIN1(trcsa_RFXS(idsa_FeOH,NU(N2,N1),N2,N1))))
      RFHFE2=VFLW*AZMAX1((trcsa_soHml2(idsa_FeOH2,N3,N2,N1)-AZMIN1(trcsa_RFXS(idsa_FeOH2,NU(N2,N1),N2,N1))))
      RFHFE3=VFLW*AZMAX1((trcsa_soHml2(idsa_FeOH3,N3,N2,N1)-AZMIN1(trcsa_RFXS(idsa_FeOH3,NU(N2,N1),N2,N1))))
      RFHFE4=VFLW*AZMAX1((trcsa_soHml2(idsa_FeOH4,N3,N2,N1)-AZMIN1(trcsa_RFXS(idsa_FeOH4,NU(N2,N1),N2,N1))))
      RFHFES=VFLW*AZMAX1((trcsa_soHml2(idsa_FeSO4,N3,N2,N1)-AZMIN1(trcsa_RFXS(idsa_FeSO4,NU(N2,N1),N2,N1))))
      RFHCAO=VFLW*AZMAX1((trcsa_soHml2(idsa_CaOH2,N3,N2,N1)-AZMIN1(trcsa_RFXS(idsa_CaOH2,NU(N2,N1),N2,N1))))
      RFHCAC=VFLW*AZMAX1((trcsa_soHml2(idsa_CaCO3,N3,N2,N1)-AZMIN1(trcsa_RFXS(idsa_CaCO3,NU(N2,N1),N2,N1))))
      RFHCAH=VFLW*AZMAX1((trcsa_soHml2(idsa_CaHCO3,N3,N2,N1)-AZMIN1(trcsa_RFXS(idsa_CaHCO3,NU(N2,N1),N2,N1))))
      RFHCAS=VFLW*AZMAX1((trcsa_soHml2(idsa_CaSO4,N3,N2,N1)-AZMIN1(trcsa_RFXS(idsa_CaSO4,NU(N2,N1),N2,N1))))
      RFHMGO=VFLW*AZMAX1((trcsa_soHml2(idsa_MgOH2,N3,N2,N1)-AZMIN1(trcsa_RFXS(idsa_MgOH2,NU(N2,N1),N2,N1))))
      RFHMGC=VFLW*AZMAX1((trcsa_soHml2(idsa_MgCO3,N3,N2,N1)-AZMIN1(trcsa_RFXS(idsa_MgCO3,NU(N2,N1),N2,N1))))
      RFHMGH=VFLW*AZMAX1((trcsa_soHml2(idsa_MgHCO3,N3,N2,N1)-AZMIN1(trcsa_RFXS(idsa_MgHCO3,NU(N2,N1),N2,N1))))
      RFHMGS=VFLW*AZMAX1((trcsa_soHml2(idsa_MgSO4,N3,N2,N1)-AZMIN1(trcsa_RFXS(idsa_MgSO4,NU(N2,N1),N2,N1))))
      RFHNAC=VFLW*AZMAX1((trcsa_soHml2(idsa_NaCO3,N3,N2,N1)-AZMIN1(trcsa_RFXS(idsa_NaCO3,NU(N2,N1),N2,N1))))
      RFHNAS=VFLW*AZMAX1((trcsa_soHml2(idsa_NaSO4,N3,N2,N1)-AZMIN1(trcsa_RFXS(idsa_NaSO4,NU(N2,N1),N2,N1))))
      RFHKAS=VFLW*AZMAX1((trcsa_soHml2(idsa_KSO4,N3,N2,N1)-AZMIN1(trcsa_RFXS(idsa_KSO4,NU(N2,N1),N2,N1))))
      RFHH0P=VFLW*AZMAX1((trcsa_soHml2(idsa_H0PO4,N3,N2,N1)-AZMIN1(trcsa_RFXS(idsa_H0PO4,NU(N2,N1),N2,N1))))*trcs_VLN(ids_H1PO4,N3,N2,N1)
      RFHH3P=VFLW*AZMAX1((trcsa_soHml2(idsa_H3PO4,N3,N2,N1)-AZMIN1(trcsa_RFXS(idsa_H3PO4,NU(N2,N1),N2,N1))))*trcs_VLN(ids_H1PO4,N3,N2,N1)
      RFHF1P=VFLW*AZMAX1((trcsa_soHml2(idsa_FeHPO4,N3,N2,N1)-AZMIN1(trcsa_RFXS(idsa_FeHPO4,NU(N2,N1),N2,N1))))*trcs_VLN(ids_H1PO4,N3,N2,N1)
      RFHF2P=VFLW*AZMAX1((trcsa_soHml2(idsa_FeH2PO4,N3,N2,N1)-AZMIN1(trcsa_RFXS(idsa_FeH2PO4,NU(N2,N1),N2,N1))))*trcs_VLN(ids_H1PO4,N3,N2,N1)
      RFHC0P=VFLW*AZMAX1((trcsa_soHml2(idsa_CaPO4,N3,N2,N1)-AZMIN1(trcsa_RFXS(idsa_CaPO4,NU(N2,N1),N2,N1))))*trcs_VLN(ids_H1PO4,N3,N2,N1)
      RFHC1P=VFLW*AZMAX1((trcsa_soHml2(idsa_CaHPO4,N3,N2,N1)-AZMIN1(trcsa_RFXS(idsa_CaHPO4,NU(N2,N1),N2,N1))))*trcs_VLN(ids_H1PO4,N3,N2,N1)
      RFHC2P=VFLW*AZMAX1((trcsa_soHml2(idsa_CaH2PO4,N3,N2,N1)-AZMIN1(trcsa_RFXS(idsa_CaH2PO4,NU(N2,N1),N2,N1))))*trcs_VLN(ids_H1PO4,N3,N2,N1)
      RFHM1P=VFLW*AZMAX1((trcsa_soHml2(idsa_MgHPO4,N3,N2,N1)-AZMIN1(trcsa_RFXS(idsa_MgHPO4,NU(N2,N1),N2,N1))))*trcs_VLN(ids_H1PO4,N3,N2,N1)
      RFHH0B=VFLW*AZMAX1((trcsa_soHml2(idsa_H0PO4B,N3,N2,N1)-AZMIN1(trcsa_RFXS(idsa_H0PO4B,NU(N2,N1),N2,N1))))*trcs_VLN(ids_H1PO4B,N3,N2,N1)
      RFHH3B=VFLW*AZMAX1((trcsa_soHml2(idsa_H3PO4B,N3,N2,N1)-AZMIN1(trcsa_RFXS(idsa_H3PO4B,NU(N2,N1),N2,N1))))*trcs_VLN(ids_H1PO4B,N3,N2,N1)
      RFHF1B=VFLW*AZMAX1((trcsa_soHml2(idsa_FeHPO4B,N3,N2,N1)-AZMIN1(trcsa_RFXS(idsa_FeHPO4B,NU(N2,N1),N2,N1))))*trcs_VLN(ids_H1PO4B,N3,N2,N1)
      RFHF2B=VFLW*AZMAX1((trcsa_soHml2(idsa_FeH2PO4B,N3,N2,N1)-AZMIN1(trcsa_RFXS(idsa_FeH2PO4B,NU(N2,N1),N2,N1))))*trcs_VLN(ids_H1PO4B,N3,N2,N1)
      RFHC0B=VFLW*AZMAX1((trcsa_soHml2(idsa_CaPO4B,N3,N2,N1)-AZMIN1(trcsa_RFXS(idsa_CaPO4B,NU(N2,N1),N2,N1))))*trcs_VLN(ids_H1PO4B,N3,N2,N1)
      RFHC1B=VFLW*AZMAX1((trcsa_soHml2(idsa_CaHPO4B,N3,N2,N1)-AZMIN1(trcsa_RFXS(idsa_CaHPO4B,NU(N2,N1),N2,N1))))*trcs_VLN(ids_H1PO4B,N3,N2,N1)
      RFHC2B=VFLW*AZMAX1((trcsa_soHml2(idsa_CaH2PO4B,N3,N2,N1)-AZMIN1(trcsa_RFXS(idsa_CaH2PO4B,NU(N2,N1),N2,N1))))*trcs_VLN(ids_H1PO4B,N3,N2,N1)
      RFHM1B=VFLW*AZMAX1((trcsa_soHml2(idsa_MgHPO4B,N3,N2,N1)-AZMIN1(trcsa_RFXS(idsa_MgHPO4B,NU(N2,N1),N2,N1))))*trcs_VLN(ids_H1PO4B,N3,N2,N1)
!
!     OTHERWISE
!
    ELSE
      RFHAL=VFLW*AZMAX1(trcsa_soHml2(idsa_Al,N3,N2,N1))
      RFHFE=VFLW*AZMAX1(trcsa_soHml2(idsa_Fe,N3,N2,N1))
      RFHHY=VFLW*AZMAX1(trcsa_soHml2(idsa_Hp,N3,N2,N1))
      RFHCA=VFLW*AZMAX1(trcsa_soHml2(idsa_Ca,N3,N2,N1))
      RFHMG=VFLW*AZMAX1(trcsa_soHml2(idsa_Mg,N3,N2,N1))
      RFHNA=VFLW*AZMAX1(trcsa_soHml2(idsa_Na,N3,N2,N1))
      RFHKA=VFLW*AZMAX1(trcsa_soHml2(idsa_K,N3,N2,N1))
      RFHOH=VFLW*AZMAX1(trcsa_soHml2(idsa_OH,N3,N2,N1))
      RFHSO=VFLW*AZMAX1(trcsa_soHml2(idsa_SO4,N3,N2,N1))
      RFHCL=VFLW*AZMAX1(trcsa_soHml2(idsa_Cl,N3,N2,N1))
      RFHC3=VFLW*AZMAX1(trcsa_soHml2(idsa_CO3,N3,N2,N1))
      RFHHC=VFLW*AZMAX1(trcsa_soHml2(idsa_HCO3,N3,N2,N1))
      RFHAL1=VFLW*AZMAX1(trcsa_soHml2(idsa_AlOH,N3,N2,N1))
      RFHAL2=VFLW*AZMAX1(trcsa_soHml2(idsa_AlOH2,N3,N2,N1))
      RFHAL3=VFLW*AZMAX1(trcsa_soHml2(idsa_AlOH3,N3,N2,N1))
      RFHAL4=VFLW*AZMAX1(trcsa_soHml2(idsa_AlOH4,N3,N2,N1))
      RFHALS=VFLW*AZMAX1(trcsa_soHml2(idsa_AlSO4,N3,N2,N1))
      RFHFE1=VFLW*AZMAX1(trcsa_soHml2(idsa_FeOH,N3,N2,N1))
      RFHFE2=VFLW*AZMAX1(trcsa_soHml2(idsa_FeOH2,N3,N2,N1))
      RFHFE3=VFLW*AZMAX1(trcsa_soHml2(idsa_FeOH3,N3,N2,N1))
      RFHFE4=VFLW*AZMAX1(trcsa_soHml2(idsa_FeOH4,N3,N2,N1))
      RFHFES=VFLW*AZMAX1(trcsa_soHml2(idsa_FeSO4,N3,N2,N1))
      RFHCAO=VFLW*AZMAX1(trcsa_soHml2(idsa_CaOH2,N3,N2,N1))
      RFHCAC=VFLW*AZMAX1(trcsa_soHml2(idsa_CaCO3,N3,N2,N1))
      RFHCAH=VFLW*AZMAX1(trcsa_soHml2(idsa_CaHCO3,N3,N2,N1))
      RFHCAS=VFLW*AZMAX1(trcsa_soHml2(idsa_CaSO4,N3,N2,N1))
      RFHMGO=VFLW*AZMAX1(trcsa_soHml2(idsa_MgOH2,N3,N2,N1))
      RFHMGC=VFLW*AZMAX1(trcsa_soHml2(idsa_MgCO3,N3,N2,N1))
      RFHMGH=VFLW*AZMAX1(trcsa_soHml2(idsa_MgHCO3,N3,N2,N1))
      RFHMGS=VFLW*AZMAX1(trcsa_soHml2(idsa_MgSO4,N3,N2,N1))
      RFHNAC=VFLW*AZMAX1(trcsa_soHml2(idsa_NaCO3,N3,N2,N1))
      RFHNAS=VFLW*AZMAX1(trcsa_soHml2(idsa_NaSO4,N3,N2,N1))
      RFHKAS=VFLW*AZMAX1(trcsa_soHml2(idsa_KSO4,N3,N2,N1))
      RFHH0P=VFLW*AZMAX1(trcsa_soHml2(idsa_H0PO4,N3,N2,N1))*trcs_VLN(ids_H1PO4,N6,N5,N4)
      RFHH3P=VFLW*AZMAX1(trcsa_soHml2(idsa_H3PO4,N3,N2,N1))*trcs_VLN(ids_H1PO4,N6,N5,N4)
      RFHF1P=VFLW*AZMAX1(trcsa_soHml2(idsa_FeHPO4,N3,N2,N1))*trcs_VLN(ids_H1PO4,N6,N5,N4)
      RFHF2P=VFLW*AZMAX1(trcsa_soHml2(idsa_FeH2PO4,N3,N2,N1))*trcs_VLN(ids_H1PO4,N6,N5,N4)
      RFHC0P=VFLW*AZMAX1(trcsa_soHml2(idsa_CaPO4,N3,N2,N1))*trcs_VLN(ids_H1PO4,N6,N5,N4)
      RFHC1P=VFLW*AZMAX1(trcsa_soHml2(idsa_CaHPO4,N3,N2,N1))*trcs_VLN(ids_H1PO4,N6,N5,N4)
      RFHC2P=VFLW*AZMAX1(trcsa_soHml2(idsa_CaH2PO4,N3,N2,N1))*trcs_VLN(ids_H1PO4,N6,N5,N4)
      RFHM1P=VFLW*AZMAX1(trcsa_soHml2(idsa_MgHPO4,N3,N2,N1))*trcs_VLN(ids_H1PO4,N6,N5,N4)
      RFHH0B=VFLW*AZMAX1(trcsa_soHml2(idsa_H0PO4B,N3,N2,N1))*trcs_VLN(ids_H1PO4B,N6,N5,N4)
      RFHH3B=VFLW*AZMAX1(trcsa_soHml2(idsa_H3PO4B,N3,N2,N1))*trcs_VLN(ids_H1PO4B,N6,N5,N4)
      RFHF1B=VFLW*AZMAX1(trcsa_soHml2(idsa_FeHPO4B,N3,N2,N1))*trcs_VLN(ids_H1PO4B,N6,N5,N4)
      RFHF2B=VFLW*AZMAX1(trcsa_soHml2(idsa_FeH2PO4B,N3,N2,N1))*trcs_VLN(ids_H1PO4B,N6,N5,N4)
      RFHC0B=VFLW*AZMAX1(trcsa_soHml2(idsa_CaPO4B,N3,N2,N1))*trcs_VLN(ids_H1PO4B,N6,N5,N4)
      RFHC1B=VFLW*AZMAX1(trcsa_soHml2(idsa_CaHPO4B,N3,N2,N1))*trcs_VLN(ids_H1PO4B,N6,N5,N4)
      RFHC2B=VFLW*AZMAX1(trcsa_soHml2(idsa_CaH2PO4B,N3,N2,N1))*trcs_VLN(ids_H1PO4B,N6,N5,N4)
      RFHM1B=VFLW*AZMAX1(trcsa_soHml2(idsa_MgHPO4B,N3,N2,N1))*trcs_VLN(ids_H1PO4B,N6,N5,N4)
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
    RFHAL=VFLW*AZMAX1(trcsa_soHml2(idsa_Al,N6,N5,N4))
    RFHFE=VFLW*AZMAX1(trcsa_soHml2(idsa_Fe,N6,N5,N4))
    RFHHY=VFLW*AZMAX1(trcsa_soHml2(idsa_Hp,N6,N5,N4))
    RFHCA=VFLW*AZMAX1(trcsa_soHml2(idsa_Ca,N6,N5,N4))
    RFHMG=VFLW*AZMAX1(trcsa_soHml2(idsa_Mg,N6,N5,N4))
    RFHNA=VFLW*AZMAX1(trcsa_soHml2(idsa_Na,N6,N5,N4))
    RFHKA=VFLW*AZMAX1(trcsa_soHml2(idsa_K,N6,N5,N4))
    RFHOH=VFLW*AZMAX1(trcsa_soHml2(idsa_OH,N6,N5,N4))
    RFHSO=VFLW*AZMAX1(trcsa_soHml2(idsa_SO4,N6,N5,N4))
    RFHCL=VFLW*AZMAX1(trcsa_soHml2(idsa_Cl,N6,N5,N4))
    RFHC3=VFLW*AZMAX1(trcsa_soHml2(idsa_CO3,N6,N5,N4))
    RFHHC=VFLW*AZMAX1(trcsa_soHml2(idsa_HCO3,N6,N5,N4))
    RFHAL1=VFLW*AZMAX1(trcsa_soHml2(idsa_AlOH,N6,N5,N4))
    RFHAL2=VFLW*AZMAX1(trcsa_soHml2(idsa_AlOH2,N6,N5,N4))
    RFHAL3=VFLW*AZMAX1(trcsa_soHml2(idsa_AlOH3,N6,N5,N4))
    RFHAL4=VFLW*AZMAX1(trcsa_soHml2(idsa_AlOH4,N6,N5,N4))
    RFHALS=VFLW*AZMAX1(trcsa_soHml2(idsa_AlSO4,N6,N5,N4))
    RFHFE1=VFLW*AZMAX1(trcsa_soHml2(idsa_FeOH,N6,N5,N4))
    RFHFE2=VFLW*AZMAX1(trcsa_soHml2(idsa_FeOH2,N6,N5,N4))
    RFHFE3=VFLW*AZMAX1(trcsa_soHml2(idsa_FeOH3,N6,N5,N4))
    RFHFE4=VFLW*AZMAX1(trcsa_soHml2(idsa_FeOH4,N6,N5,N4))
    RFHFES=VFLW*AZMAX1(trcsa_soHml2(idsa_FeSO4,N6,N5,N4))
    RFHCAO=VFLW*AZMAX1(trcsa_soHml2(idsa_CaOH2,N6,N5,N4))
    RFHCAC=VFLW*AZMAX1(trcsa_soHml2(idsa_CaCO3,N6,N5,N4))
    RFHCAH=VFLW*AZMAX1(trcsa_soHml2(idsa_CaHCO3,N6,N5,N4))
    RFHCAS=VFLW*AZMAX1(trcsa_soHml2(idsa_CaSO4,N6,N5,N4))
    RFHMGO=VFLW*AZMAX1(trcsa_soHml2(idsa_MgOH2,N6,N5,N4))
    RFHMGC=VFLW*AZMAX1(trcsa_soHml2(idsa_MgCO3,N6,N5,N4))
    RFHMGH=VFLW*AZMAX1(trcsa_soHml2(idsa_MgHCO3,N6,N5,N4))
    RFHMGS=VFLW*AZMAX1(trcsa_soHml2(idsa_MgSO4,N6,N5,N4))
    RFHNAC=VFLW*AZMAX1(trcsa_soHml2(idsa_NaCO3,N6,N5,N4))
    RFHNAS=VFLW*AZMAX1(trcsa_soHml2(idsa_NaSO4,N6,N5,N4))
    RFHKAS=VFLW*AZMAX1(trcsa_soHml2(idsa_KSO4,N6,N5,N4))
    RFHH0P=VFLW*AZMAX1(trcsa_soHml2(idsa_H0PO4,N6,N5,N4))*trcs_VLN(ids_H1PO4,N6,N5,N4)
    RFHH3P=VFLW*AZMAX1(trcsa_soHml2(idsa_H3PO4,N6,N5,N4))*trcs_VLN(ids_H1PO4,N6,N5,N4)
    RFHF1P=VFLW*AZMAX1(trcsa_soHml2(idsa_FeHPO4,N6,N5,N4))*trcs_VLN(ids_H1PO4,N6,N5,N4)
    RFHF2P=VFLW*AZMAX1(trcsa_soHml2(idsa_FeH2PO4,N6,N5,N4))*trcs_VLN(ids_H1PO4,N6,N5,N4)
    RFHC0P=VFLW*AZMAX1(trcsa_soHml2(idsa_CaPO4,N6,N5,N4))*trcs_VLN(ids_H1PO4,N6,N5,N4)
    RFHC1P=VFLW*AZMAX1(trcsa_soHml2(idsa_CaHPO4,N6,N5,N4))*trcs_VLN(ids_H1PO4,N6,N5,N4)
    RFHC2P=VFLW*AZMAX1(trcsa_soHml2(idsa_CaH2PO4,N6,N5,N4))*trcs_VLN(ids_H1PO4,N6,N5,N4)
    RFHM1P=VFLW*AZMAX1(trcsa_soHml2(idsa_MgHPO4,N6,N5,N4))*trcs_VLN(ids_H1PO4,N6,N5,N4)
    RFHH0B=VFLW*AZMAX1(trcsa_soHml2(idsa_H0PO4B,N6,N5,N4))*trcs_VLN(ids_H1PO4B,N6,N5,N4)
    RFHH3B=VFLW*AZMAX1(trcsa_soHml2(idsa_H3PO4B,N6,N5,N4))*trcs_VLN(ids_H1PO4B,N6,N5,N4)
    RFHF1B=VFLW*AZMAX1(trcsa_soHml2(idsa_FeHPO4B,N6,N5,N4))*trcs_VLN(ids_H1PO4B,N6,N5,N4)
    RFHF2B=VFLW*AZMAX1(trcsa_soHml2(idsa_FeH2PO4B,N6,N5,N4))*trcs_VLN(ids_H1PO4B,N6,N5,N4)
    RFHC0B=VFLW*AZMAX1(trcsa_soHml2(idsa_CaPO4B,N6,N5,N4))*trcs_VLN(ids_H1PO4B,N6,N5,N4)
    RFHC1B=VFLW*AZMAX1(trcsa_soHml2(idsa_CaHPO4B,N6,N5,N4))*trcs_VLN(ids_H1PO4B,N6,N5,N4)
    RFHC2B=VFLW*AZMAX1(trcsa_soHml2(idsa_CaH2PO4B,N6,N5,N4))*trcs_VLN(ids_H1PO4B,N6,N5,N4)
    RFHM1B=VFLW*AZMAX1(trcsa_soHml2(idsa_MgHPO4B,N6,N5,N4))*trcs_VLN(ids_H1PO4B,N6,N5,N4)
!
!     NO MACROPORE FLUX
!
  ELSE
    RFHAL=0.0_r8
    RFHFE=0.0_r8
    RFHHY=0.0_r8
    RFHCA=0.0_r8
    RFHMG=0.0_r8
    RFHNA=0.0_r8
    RFHKA=0.0_r8
    RFHOH=0.0_r8
    RFHSO=0.0_r8
    RFHCL=0.0_r8
    RFHC3=0.0_r8
    RFHHC=0.0_r8
    RFHAL1=0.0_r8
    RFHAL2=0.0_r8
    RFHAL3=0.0_r8
    RFHAL4=0.0_r8
    RFHALS=0.0_r8
    RFHFE1=0.0_r8
    RFHFE2=0.0_r8
    RFHFE3=0.0_r8
    RFHFE4=0.0_r8
    RFHFES=0.0_r8
    RFHCAO=0.0_r8
    RFHCAC=0.0_r8
    RFHCAH=0.0_r8
    RFHCAS=0.0_r8
    RFHMGO=0.0_r8
    RFHMGC=0.0_r8
    RFHMGH=0.0_r8
    RFHMGS=0.0_r8
    RFHNAC=0.0_r8
    RFHNAS=0.0_r8
    RFHKAS=0.0_r8
    RFHH0P=0.0_r8
    RFHH3P=0.0_r8
    RFHF1P=0.0_r8
    RFHF2P=0.0_r8
    RFHC0P=0.0_r8
    RFHC1P=0.0_r8
    RFHC2P=0.0_r8
    RFHM1P=0.0_r8
    RFHH0B=0.0_r8
    RFHH3B=0.0_r8
    RFHF1B=0.0_r8
    RFHF2B=0.0_r8
    RFHC0B=0.0_r8
    RFHC1B=0.0_r8
    RFHC2B=0.0_r8
    RFHM1B=0.0_r8
  ENDIF
!     RFH*=convective flux through macropores
!     DFH*=diffusive solute flux through macropores
!

  RALFHS(N,N6,N5,N4)=RFHAL
  RFEFHS(N,N6,N5,N4)=RFHFE
  RHYFHS(N,N6,N5,N4)=RFHHY
  RCAFHS(N,N6,N5,N4)=RFHCA
  RMGFHS(N,N6,N5,N4)=RFHMG
  RNAFHS(N,N6,N5,N4)=RFHNA
  RKAFHS(N,N6,N5,N4)=RFHKA
  ROHFHS(N,N6,N5,N4)=RFHOH
  RSOFHS(N,N6,N5,N4)=RFHSO
  RCLFHS(N,N6,N5,N4)=RFHCL
  RC3FHS(N,N6,N5,N4)=RFHC3
  RHCFHS(N,N6,N5,N4)=RFHHC
  RAL1HS(N,N6,N5,N4)=RFHAL1
  RAL2HS(N,N6,N5,N4)=RFHAL2
  RAL3HS(N,N6,N5,N4)=RFHAL3
  RAL4HS(N,N6,N5,N4)=RFHAL4
  RALSHS(N,N6,N5,N4)=RFHALS
  RFE1HS(N,N6,N5,N4)=RFHFE1
  RFE2HS(N,N6,N5,N4)=RFHFE2
  RFE3HS(N,N6,N5,N4)=RFHFE3
  RFE4HS(N,N6,N5,N4)=RFHFE4
  RFESHS(N,N6,N5,N4)=RFHFES
  RCAOHS(N,N6,N5,N4)=RFHCAO
  RCACHS(N,N6,N5,N4)=RFHCAC
  RCAHHS(N,N6,N5,N4)=RFHCAH
  RCASHS(N,N6,N5,N4)=RFHCAS
  RMGOHS(N,N6,N5,N4)=RFHMGO
  RMGCHS(N,N6,N5,N4)=RFHMGC
  RMGHHS(N,N6,N5,N4)=RFHMGH
  RMGSHS(N,N6,N5,N4)=RFHMGS
  RNACHS(N,N6,N5,N4)=RFHNAC
  RNASHS(N,N6,N5,N4)=RFHNAS
  RKASHS(N,N6,N5,N4)=RFHKAS
  RH0PHS(N,N6,N5,N4)=RFHH0P
  RH3PHS(N,N6,N5,N4)=RFHH3P
  RF1PHS(N,N6,N5,N4)=RFHF1P
  RF2PHS(N,N6,N5,N4)=RFHF2P
  RC0PHS(N,N6,N5,N4)=RFHC0P
  RC1PHS(N,N6,N5,N4)=RFHC1P
  RC2PHS(N,N6,N5,N4)=RFHC2P
  RM1PHS(N,N6,N5,N4)=RFHM1P
  RH0BHB(N,N6,N5,N4)=RFHH0B
  RH3BHB(N,N6,N5,N4)=RFHH3B
  RF1BHB(N,N6,N5,N4)=RFHF1B
  RF2BHB(N,N6,N5,N4)=RFHF2B
  RC0BHB(N,N6,N5,N4)=RFHC0B
  RC1BHB(N,N6,N5,N4)=RFHC1B
  RC2BHB(N,N6,N5,N4)=RFHC2B
  RM1BHB(N,N6,N5,N4)=RFHM1B

  end subroutine SoluteAdvMacropore

!----------------------------------------------------------------------
  subroutine SoluteDifsMacropore(M,N,N1,N2,N3,N4,N5,N6)
  implicit none
  integer, intent(in) :: M,N,N1,N2,N3,N4,N5,N6

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
    CAL1=AZMAX1(trcsa_soHml2(idsa_Al,N3,N2,N1)/VOLWHM(M,N3,N2,N1))
    CFE1=AZMAX1(trcsa_soHml2(idsa_Fe,N3,N2,N1)/VOLWHM(M,N3,N2,N1))
    CHY1=AZMAX1(trcsa_soHml2(idsa_Hp,N3,N2,N1)/VOLWHM(M,N3,N2,N1))
    CCA1=AZMAX1(trcsa_soHml2(idsa_Ca,N3,N2,N1)/VOLWHM(M,N3,N2,N1))
    CMG1=AZMAX1(trcsa_soHml2(idsa_Mg,N3,N2,N1)/VOLWHM(M,N3,N2,N1))
    CNA1=AZMAX1(trcsa_soHml2(idsa_Na,N3,N2,N1)/VOLWHM(M,N3,N2,N1))
    CKA1=AZMAX1(trcsa_soHml2(idsa_K,N3,N2,N1)/VOLWHM(M,N3,N2,N1))
    COH1=AZMAX1(trcsa_soHml2(idsa_OH,N3,N2,N1)/VOLWHM(M,N3,N2,N1))
    CSO41=AZMAX1(trcsa_soHml2(idsa_SO4,N3,N2,N1)/VOLWHM(M,N3,N2,N1))
    CCL1=AZMAX1(trcsa_soHml2(idsa_Cl,N3,N2,N1)/VOLWHM(M,N3,N2,N1))
    CCO31=AZMAX1(trcsa_soHml2(idsa_CO3,N3,N2,N1)/VOLWHM(M,N3,N2,N1))
    CHCO31=AZMAX1(trcsa_soHml2(idsa_HCO3,N3,N2,N1)/VOLWHM(M,N3,N2,N1))
    CAL11=AZMAX1(trcsa_soHml2(idsa_AlOH,N3,N2,N1)/VOLWHM(M,N3,N2,N1))
    CAL21=AZMAX1(trcsa_soHml2(idsa_AlOH2,N3,N2,N1)/VOLWHM(M,N3,N2,N1))
    CAL31=AZMAX1(trcsa_soHml2(idsa_AlOH3,N3,N2,N1)/VOLWHM(M,N3,N2,N1))
    CAL41=AZMAX1(trcsa_soHml2(idsa_AlOH4,N3,N2,N1)/VOLWHM(M,N3,N2,N1))
    CALS1=AZMAX1(trcsa_soHml2(idsa_AlSO4,N3,N2,N1)/VOLWHM(M,N3,N2,N1))
    CFE11=AZMAX1(trcsa_soHml2(idsa_FeOH,N3,N2,N1)/VOLWHM(M,N3,N2,N1))
    CFE21=AZMAX1(trcsa_soHml2(idsa_FeOH2,N3,N2,N1)/VOLWHM(M,N3,N2,N1))
    CFE31=AZMAX1(trcsa_soHml2(idsa_FeOH3,N3,N2,N1)/VOLWHM(M,N3,N2,N1))
    CFE41=AZMAX1(trcsa_soHml2(idsa_FeOH4,N3,N2,N1)/VOLWHM(M,N3,N2,N1))
    CFES1=AZMAX1(trcsa_soHml2(idsa_FeSO4,N3,N2,N1)/VOLWHM(M,N3,N2,N1))
    CCAO1=AZMAX1(trcsa_soHml2(idsa_CaOH2,N3,N2,N1)/VOLWHM(M,N3,N2,N1))
    CCAC1=AZMAX1(trcsa_soHml2(idsa_CaCO3,N3,N2,N1)/VOLWHM(M,N3,N2,N1))
    CCAH1=AZMAX1(trcsa_soHml2(idsa_CaHCO3,N3,N2,N1)/VOLWHM(M,N3,N2,N1))
    CCAS1=AZMAX1(trcsa_soHml2(idsa_CaSO4,N3,N2,N1)/VOLWHM(M,N3,N2,N1))
    CMGO1=AZMAX1(trcsa_soHml2(idsa_MgOH2,N3,N2,N1)/VOLWHM(M,N3,N2,N1))
    CMGC1=AZMAX1(trcsa_soHml2(idsa_MgCO3,N3,N2,N1)/VOLWHM(M,N3,N2,N1))
    CMGH1=AZMAX1(trcsa_soHml2(idsa_MgHCO3,N3,N2,N1)/VOLWHM(M,N3,N2,N1))
    CMGS1=AZMAX1(trcsa_soHml2(idsa_MgSO4,N3,N2,N1)/VOLWHM(M,N3,N2,N1))
    CNAC1=AZMAX1(trcsa_soHml2(idsa_NaCO3,N3,N2,N1)/VOLWHM(M,N3,N2,N1))
    CNAS1=AZMAX1(trcsa_soHml2(idsa_NaSO4,N3,N2,N1)/VOLWHM(M,N3,N2,N1))
    CKAS1=AZMAX1(trcsa_soHml2(idsa_KSO4,N3,N2,N1)/VOLWHM(M,N3,N2,N1))
    C0PO41=AZMAX1(trcsa_soHml2(idsa_H0PO4,N3,N2,N1)/VOLWHM(M,N3,N2,N1))
    C3PO41=AZMAX1(trcsa_soHml2(idsa_H3PO4,N3,N2,N1)/VOLWHM(M,N3,N2,N1))
    CFE1P1=AZMAX1(trcsa_soHml2(idsa_FeHPO4,N3,N2,N1)/VOLWHM(M,N3,N2,N1))
    CFE2P1=AZMAX1(trcsa_soHml2(idsa_FeH2PO4,N3,N2,N1)/VOLWHM(M,N3,N2,N1))
    CCA0P1=AZMAX1(trcsa_soHml2(idsa_CaPO4,N3,N2,N1)/VOLWHM(M,N3,N2,N1))
    CCA1P1=AZMAX1(trcsa_soHml2(idsa_CaHPO4,N3,N2,N1)/VOLWHM(M,N3,N2,N1))
    CCA2P1=AZMAX1(trcsa_soHml2(idsa_CaH2PO4,N3,N2,N1)/VOLWHM(M,N3,N2,N1))
    CMG1P1=AZMAX1(trcsa_soHml2(idsa_MgHPO4,N3,N2,N1)/VOLWHM(M,N3,N2,N1))
    C0POB1=AZMAX1(trcsa_soHml2(idsa_H0PO4B,N3,N2,N1)/VOLWHM(M,N3,N2,N1))
    C3POB1=AZMAX1(trcsa_soHml2(idsa_H3PO4B,N3,N2,N1)/VOLWHM(M,N3,N2,N1))
    CF1PB1=AZMAX1(trcsa_soHml2(idsa_FeHPO4B,N3,N2,N1)/VOLWHM(M,N3,N2,N1))
    CF2PB1=AZMAX1(trcsa_soHml2(idsa_FeH2PO4B,N3,N2,N1)/VOLWHM(M,N3,N2,N1))
    CC0PB1=AZMAX1(trcsa_soHml2(idsa_CaPO4B,N3,N2,N1)/VOLWHM(M,N3,N2,N1))
    CC1PB1=AZMAX1(trcsa_soHml2(idsa_CaHPO4B,N3,N2,N1)/VOLWHM(M,N3,N2,N1))
    CC2PB1=AZMAX1(trcsa_soHml2(idsa_CaH2PO4B,N3,N2,N1)/VOLWHM(M,N3,N2,N1))
    CM1PB1=AZMAX1(trcsa_soHml2(idsa_MgHPO4B,N3,N2,N1)/VOLWHM(M,N3,N2,N1))

    CAL2=AZMAX1(trcsa_soHml2(idsa_Al,N6,N5,N4)/VOLWHM(M,N6,N5,N4))
    CFE2=AZMAX1(trcsa_soHml2(idsa_Fe,N6,N5,N4)/VOLWHM(M,N6,N5,N4))
    CHY2=AZMAX1(trcsa_soHml2(idsa_Hp,N6,N5,N4)/VOLWHM(M,N6,N5,N4))
    CCA2=AZMAX1(trcsa_soHml2(idsa_Ca,N6,N5,N4)/VOLWHM(M,N6,N5,N4))
    CMG2=AZMAX1(trcsa_soHml2(idsa_Mg,N6,N5,N4)/VOLWHM(M,N6,N5,N4))
    CNA2=AZMAX1(trcsa_soHml2(idsa_Na,N6,N5,N4)/VOLWHM(M,N6,N5,N4))
    CKA2=AZMAX1(trcsa_soHml2(idsa_K,N6,N5,N4)/VOLWHM(M,N6,N5,N4))
    COH2=AZMAX1(trcsa_soHml2(idsa_OH,N6,N5,N4)/VOLWHM(M,N6,N5,N4))
    CSO42=AZMAX1(trcsa_soHml2(idsa_SO4,N6,N5,N4)/VOLWHM(M,N6,N5,N4))
    CCL2=AZMAX1(trcsa_soHml2(idsa_Cl,N6,N5,N4)/VOLWHM(M,N6,N5,N4))
    CCO32=AZMAX1(trcsa_soHml2(idsa_CO3,N6,N5,N4)/VOLWHM(M,N6,N5,N4))
    CHCO32=AZMAX1(trcsa_soHml2(idsa_HCO3,N6,N5,N4)/VOLWHM(M,N6,N5,N4))
    CAL12=AZMAX1(trcsa_soHml2(idsa_AlOH,N6,N5,N4)/VOLWHM(M,N6,N5,N4))
    CAL22=AZMAX1(trcsa_soHml2(idsa_AlOH2,N6,N5,N4)/VOLWHM(M,N6,N5,N4))
    CAL32=AZMAX1(trcsa_soHml2(idsa_AlOH3,N6,N5,N4)/VOLWHM(M,N6,N5,N4))
    CAL42=AZMAX1(trcsa_soHml2(idsa_AlOH4,N6,N5,N4)/VOLWHM(M,N6,N5,N4))
    CALS2=AZMAX1(trcsa_soHml2(idsa_AlSO4,N6,N5,N4)/VOLWHM(M,N6,N5,N4))
    CFE12=AZMAX1(trcsa_soHml2(idsa_FeOH,N6,N5,N4)/VOLWHM(M,N6,N5,N4))
    CFE22=AZMAX1(trcsa_soHml2(idsa_FeOH2,N6,N5,N4)/VOLWHM(M,N6,N5,N4))
    CFE32=AZMAX1(trcsa_soHml2(idsa_FeOH3,N6,N5,N4)/VOLWHM(M,N6,N5,N4))
    CFE42=AZMAX1(trcsa_soHml2(idsa_FeOH4,N6,N5,N4)/VOLWHM(M,N6,N5,N4))
    CFES2=AZMAX1(trcsa_soHml2(idsa_FeSO4,N6,N5,N4)/VOLWHM(M,N6,N5,N4))
    CCAO2=AZMAX1(trcsa_soHml2(idsa_CaOH2,N6,N5,N4)/VOLWHM(M,N6,N5,N4))
    CCAC2=AZMAX1(trcsa_soHml2(idsa_CaCO3,N6,N5,N4)/VOLWHM(M,N6,N5,N4))
    CCAH2=AZMAX1(trcsa_soHml2(idsa_CaHCO3,N6,N5,N4)/VOLWHM(M,N6,N5,N4))
    CCAS2=AZMAX1(trcsa_soHml2(idsa_CaSO4,N6,N5,N4)/VOLWHM(M,N6,N5,N4))
    CMGO2=AZMAX1(trcsa_soHml2(idsa_MgOH2,N6,N5,N4)/VOLWHM(M,N6,N5,N4))
    CMGC2=AZMAX1(trcsa_soHml2(idsa_MgCO3,N6,N5,N4)/VOLWHM(M,N6,N5,N4))
    CMGH2=AZMAX1(trcsa_soHml2(idsa_MgHCO3,N6,N5,N4)/VOLWHM(M,N6,N5,N4))
    CMGS2=AZMAX1(trcsa_soHml2(idsa_MgSO4,N6,N5,N4)/VOLWHM(M,N6,N5,N4))
    CNAC2=AZMAX1(trcsa_soHml2(idsa_NaCO3,N6,N5,N4)/VOLWHM(M,N6,N5,N4))
    CNAS2=AZMAX1(trcsa_soHml2(idsa_NaSO4,N6,N5,N4)/VOLWHM(M,N6,N5,N4))
    CKAS2=AZMAX1(trcsa_soHml2(idsa_KSO4,N6,N5,N4)/VOLWHM(M,N6,N5,N4))
    C0PO42=AZMAX1(trcsa_soHml2(idsa_H0PO4,N6,N5,N4)/VOLWHM(M,N6,N5,N4))
    C3PO42=AZMAX1(trcsa_soHml2(idsa_H3PO4,N6,N5,N4)/VOLWHM(M,N6,N5,N4))
    CFE1P2=AZMAX1(trcsa_soHml2(idsa_FeHPO4,N6,N5,N4)/VOLWHM(M,N6,N5,N4))
    CFE2P2=AZMAX1(trcsa_soHml2(idsa_FeH2PO4,N6,N5,N4)/VOLWHM(M,N6,N5,N4))
    CCA0P2=AZMAX1(trcsa_soHml2(idsa_CaPO4,N6,N5,N4)/VOLWHM(M,N6,N5,N4))
    CCA1P2=AZMAX1(trcsa_soHml2(idsa_CaHPO4,N6,N5,N4)/VOLWHM(M,N6,N5,N4))
    CCA2P2=AZMAX1(trcsa_soHml2(idsa_CaH2PO4,N6,N5,N4)/VOLWHM(M,N6,N5,N4))
    CMG1P2=AZMAX1(trcsa_soHml2(idsa_MgHPO4,N6,N5,N4)/VOLWHM(M,N6,N5,N4))
    C0POB2=AZMAX1(trcsa_soHml2(idsa_H0PO4B,N6,N5,N4)/VOLWHM(M,N6,N5,N4))
    C3POB2=AZMAX1(trcsa_soHml2(idsa_H3PO4B,N6,N5,N4)/VOLWHM(M,N6,N5,N4))
    CF1PB2=AZMAX1(trcsa_soHml2(idsa_FeHPO4B,N6,N5,N4)/VOLWHM(M,N6,N5,N4))
    CF2PB2=AZMAX1(trcsa_soHml2(idsa_FeH2PO4B,N6,N5,N4)/VOLWHM(M,N6,N5,N4))
    CC0PB2=AZMAX1(trcsa_soHml2(idsa_CaPO4B,N6,N5,N4)/VOLWHM(M,N6,N5,N4))
    CC1PB2=AZMAX1(trcsa_soHml2(idsa_CaHPO4B,N6,N5,N4)/VOLWHM(M,N6,N5,N4))
    CC2PB2=AZMAX1(trcsa_soHml2(idsa_CaH2PO4B,N6,N5,N4)/VOLWHM(M,N6,N5,N4))
    CM1PB2=AZMAX1(trcsa_soHml2(idsa_MgHPO4B,N6,N5,N4)/VOLWHM(M,N6,N5,N4))
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
    DIFFE=(FESGL2(N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)
    DIFHY=(HYSGL2(N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)
    DIFCA=(CASGL2(N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)
    DIFMG=(GMSGL2(N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)
    DIFNA=(ANSGL2(N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)
    DIFKA=(AKSGL2(N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)
    DIFOH=(OHSGL2(N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)
    DIFSO=(SOSGL2(N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)
    DIFCL=(CLSXL2(N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)
    DIFC3=(C3SGL2(N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)
    DIFHC=(HCSGL2(N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)
!
!     DIFFUSIVE FLUXES BETWEEN CURRENT AND ADJACENT GRID CELL
!     MACROPORES
!
    DFHAL=DIFAL*(CAL1-CAL2)
    DFHFE=DIFFE*(CFE1-CFE2)
    DFHHY=DIFHY*(CHY1-CHY2)
    DFHCA=DIFCA*(CCA1-CCA2)
    DFHMG=DIFMG*(CMG1-CMG2)
    DFHNA=DIFNA*(CNA1-CNA2)
    DFHKA=DIFKA*(CKA1-CKA2)
    DFHOH=DIFOH*(COH1-COH2)
    DFHSO=DIFSO*(CSO41-CSO42)
    DFHCL=DIFCL*(CCL1-CCL2)
    DFHC3=DIFC3*(CCO31-CCO32)
    DFHHC=DIFHC*(CHCO31-CHCO32)
    DFHAL1=DIFAL*(CAL11-CAL12)
    DFHAL2=DIFAL*(CAL21-CAL22)
    DFHAL3=DIFAL*(CAL31-CAL32)
    DFHAL4=DIFAL*(CAL41-CAL42)
    DFHALS=DIFAL*(CALS1-CALS2)
    DFHFE1=DIFFE*(CFE11-CFE12)
    DFHFE2=DIFFE*(CFE21-CFE22)
    DFHFE3=DIFFE*(CFE31-CFE32)
    DFHFE4=DIFFE*(CFE41-CFE42)
    DFHFES=DIFFE*(CFES1-CFES2)
    DFHCAO=DIFCA*(CCAO1-CCAO2)
    DFHCAC=DIFCA*(CCAC1-CCAC2)
    DFHCAH=DIFCA*(CCAH1-CCAH2)
    DFHCAS=DIFCA*(CCAS1-CCAS2)
    DFHMGO=DIFMG*(CMGO1-CMGO2)
    DFHMGC=DIFMG*(CMGC1-CMGC2)
    DFHMGH=DIFMG*(CMGH1-CMGH2)
    DFHMGS=DIFMG*(CMGS1-CMGS2)
    DFHNAC=DIFNA*(CNAC1-CNAC2)
    DFHNAS=DIFNA*(CNAS1-CNAS2)
    DFHKAS=DIFKA*(CKAS1-CKAS2)
    DFHH0P=DIFPO*(C0PO41-C0PO42)*trcs_VLN(ids_H1PO4,N6,N5,N4)
    DFHH3P=DIFPO*(C3PO41-C3PO42)*trcs_VLN(ids_H1PO4,N6,N5,N4)
    DFHF1P=DIFPO*(CFE1P1-CFE1P2)*trcs_VLN(ids_H1PO4,N6,N5,N4)
    DFHF2P=DIFPO*(CFE2P1-CFE2P2)*trcs_VLN(ids_H1PO4,N6,N5,N4)
    DFHC0P=DIFPO*(CCA0P1-CCA0P2)*trcs_VLN(ids_H1PO4,N6,N5,N4)
    DFHC1P=DIFPO*(CCA1P1-CCA1P2)*trcs_VLN(ids_H1PO4,N6,N5,N4)
    DFHC2P=DIFPO*(CCA2P1-CCA2P2)*trcs_VLN(ids_H1PO4,N6,N5,N4)
    DFHM1P=DIFPO*(CMG1P1-CMG1P2)*trcs_VLN(ids_H1PO4,N6,N5,N4)
    DFHH0B=DIFPO*(C0POB1-C0POB2)*trcs_VLN(ids_H1PO4B,N6,N5,N4)
    DFHH3B=DIFPO*(C3POB1-C3POB2)*trcs_VLN(ids_H1PO4B,N6,N5,N4)
    DFHF1B=DIFPO*(CF1PB1-CF1PB2)*trcs_VLN(ids_H1PO4B,N6,N5,N4)
    DFHF2B=DIFPO*(CF2PB1-CF2PB2)*trcs_VLN(ids_H1PO4B,N6,N5,N4)
    DFHC0B=DIFPO*(CC0PB1-CC0PB2)*trcs_VLN(ids_H1PO4B,N6,N5,N4)
    DFHC1B=DIFPO*(CC1PB1-CC1PB2)*trcs_VLN(ids_H1PO4B,N6,N5,N4)
    DFHC2B=DIFPO*(CC2PB1-CC2PB2)*trcs_VLN(ids_H1PO4B,N6,N5,N4)
    DFHM1B=DIFPO*(CM1PB1-CM1PB2)*trcs_VLN(ids_H1PO4B,N6,N5,N4)
  ELSE
    DFHAL=0.0_r8
    DFHFE=0.0_r8
    DFHHY=0.0_r8
    DFHCA=0.0_r8
    DFHMG=0.0_r8
    DFHNA=0.0_r8
    DFHKA=0.0_r8
    DFHOH=0.0_r8
    DFHSO=0.0_r8
    DFHCL=0.0_r8
    DFHC3=0.0_r8
    DFHHC=0.0_r8
    DFHAL1=0.0_r8
    DFHAL2=0.0_r8
    DFHAL3=0.0_r8
    DFHAL4=0.0_r8
    DFHALS=0.0_r8
    DFHFE1=0.0_r8
    DFHFE2=0.0_r8
    DFHFE3=0.0_r8
    DFHFE4=0.0_r8
    DFHFES=0.0_r8
    DFHCAO=0.0_r8
    DFHCAC=0.0_r8
    DFHCAH=0.0_r8
    DFHCAS=0.0_r8
    DFHMGO=0.0_r8
    DFHMGC=0.0_r8
    DFHMGH=0.0_r8
    DFHMGS=0.0_r8
    DFHNAC=0.0_r8
    DFHNAS=0.0_r8
    DFHKAS=0.0_r8
    DFHH0P=0.0_r8
    DFHH3P=0.0_r8
    DFHF1P=0.0_r8
    DFHF2P=0.0_r8
    DFHC0P=0.0_r8
    DFHC1P=0.0_r8
    DFHC2P=0.0_r8
    DFHM1P=0.0_r8
    DFHH0B=0.0_r8
    DFHH3B=0.0_r8
    DFHF1B=0.0_r8
    DFHF2B=0.0_r8
    DFHC0B=0.0_r8
    DFHC1B=0.0_r8
    DFHC2B=0.0_r8
    DFHM1B=0.0_r8
  ENDIF
!     RFH*=convective flux through macropores
!     DFH*=diffusive solute flux through macropores
!

  RALFHS(N,N6,N5,N4)=RALFHS(N,N6,N5,N4)+DFHAL
  RFEFHS(N,N6,N5,N4)=RFEFHS(N,N6,N5,N4)+DFHFE
  RHYFHS(N,N6,N5,N4)=RHYFHS(N,N6,N5,N4)+DFHHY
  RCAFHS(N,N6,N5,N4)=RCAFHS(N,N6,N5,N4)+DFHCA
  RMGFHS(N,N6,N5,N4)=RMGFHS(N,N6,N5,N4)+DFHMG
  RNAFHS(N,N6,N5,N4)=RNAFHS(N,N6,N5,N4)+DFHNA
  RKAFHS(N,N6,N5,N4)=RKAFHS(N,N6,N5,N4)+DFHKA
  ROHFHS(N,N6,N5,N4)=ROHFHS(N,N6,N5,N4)+DFHOH
  RSOFHS(N,N6,N5,N4)=RSOFHS(N,N6,N5,N4)+DFHSO
  RCLFHS(N,N6,N5,N4)=RCLFHS(N,N6,N5,N4)+DFHCL
  RC3FHS(N,N6,N5,N4)=RC3FHS(N,N6,N5,N4)+DFHC3
  RHCFHS(N,N6,N5,N4)=RHCFHS(N,N6,N5,N4)+DFHHC
  RAL1HS(N,N6,N5,N4)=RAL1HS(N,N6,N5,N4)+DFHAL1
  RAL2HS(N,N6,N5,N4)=RAL2HS(N,N6,N5,N4)+DFHAL2
  RAL3HS(N,N6,N5,N4)=RAL3HS(N,N6,N5,N4)+DFHAL3
  RAL4HS(N,N6,N5,N4)=RAL4HS(N,N6,N5,N4)+DFHAL4
  RALSHS(N,N6,N5,N4)=RALSHS(N,N6,N5,N4)+DFHALS
  RFE1HS(N,N6,N5,N4)=RFE1HS(N,N6,N5,N4)+DFHFE1
  RFE2HS(N,N6,N5,N4)=RFE2HS(N,N6,N5,N4)+DFHFE2
  RFE3HS(N,N6,N5,N4)=RFE3HS(N,N6,N5,N4)+DFHFE3
  RFE4HS(N,N6,N5,N4)=RFE4HS(N,N6,N5,N4)+DFHFE4
  RFESHS(N,N6,N5,N4)=RFESHS(N,N6,N5,N4)+DFHFES
  RCAOHS(N,N6,N5,N4)=RCAOHS(N,N6,N5,N4)+DFHCAO
  RCACHS(N,N6,N5,N4)=RCACHS(N,N6,N5,N4)+DFHCAC
  RCAHHS(N,N6,N5,N4)=RCAHHS(N,N6,N5,N4)+DFHCAH
  RCASHS(N,N6,N5,N4)=RCASHS(N,N6,N5,N4)+DFHCAS
  RMGOHS(N,N6,N5,N4)=RMGOHS(N,N6,N5,N4)+DFHMGO
  RMGCHS(N,N6,N5,N4)=RMGCHS(N,N6,N5,N4)+DFHMGC
  RMGHHS(N,N6,N5,N4)=RMGHHS(N,N6,N5,N4)+DFHMGH
  RMGSHS(N,N6,N5,N4)=RMGSHS(N,N6,N5,N4)+DFHMGS
  RNACHS(N,N6,N5,N4)=RNACHS(N,N6,N5,N4)+DFHNAC
  RNASHS(N,N6,N5,N4)=RNASHS(N,N6,N5,N4)+DFHNAS
  RKASHS(N,N6,N5,N4)=RKASHS(N,N6,N5,N4)+DFHKAS
  RH0PHS(N,N6,N5,N4)=RH0PHS(N,N6,N5,N4)+DFHH0P
  RH3PHS(N,N6,N5,N4)=RH3PHS(N,N6,N5,N4)+DFHH3P
  RF1PHS(N,N6,N5,N4)=RF1PHS(N,N6,N5,N4)+DFHF1P
  RF2PHS(N,N6,N5,N4)=RF2PHS(N,N6,N5,N4)+DFHF2P
  RC0PHS(N,N6,N5,N4)=RC0PHS(N,N6,N5,N4)+DFHC0P
  RC1PHS(N,N6,N5,N4)=RC1PHS(N,N6,N5,N4)+DFHC1P
  RC2PHS(N,N6,N5,N4)=RC2PHS(N,N6,N5,N4)+DFHC2P
  RM1PHS(N,N6,N5,N4)=RM1PHS(N,N6,N5,N4)+DFHM1P
  RH0BHB(N,N6,N5,N4)=RH0BHB(N,N6,N5,N4)+DFHH0B
  RH3BHB(N,N6,N5,N4)=RH3BHB(N,N6,N5,N4)+DFHH3B
  RF1BHB(N,N6,N5,N4)=RF1BHB(N,N6,N5,N4)+DFHF1B
  RF2BHB(N,N6,N5,N4)=RF2BHB(N,N6,N5,N4)+DFHF2B
  RC0BHB(N,N6,N5,N4)=RC0BHB(N,N6,N5,N4)+DFHC0B
  RC1BHB(N,N6,N5,N4)=RC1BHB(N,N6,N5,N4)+DFHC1B
  RC2BHB(N,N6,N5,N4)=RC2BHB(N,N6,N5,N4)+DFHC2B
  RM1BHB(N,N6,N5,N4)=RM1BHB(N,N6,N5,N4)+DFHM1B
  end subroutine SoluteDifsMacropore
!----------------------------------------------------------------------
  subroutine SoluteAdvExchMicMacpores(M,N,NY,NX,N1,N2,N3,N4,N5,N6)
  implicit none
  integer, intent(in) :: M,N,NY,NX,N1,N2,N3,N4,N5,N6

  real(r8) :: VFLW
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
    RFLAL=VFLW*AZMAX1(trcsa_soHml2(idsa_Al,N6,N5,N4))
    RFLFE=VFLW*AZMAX1(trcsa_soHml2(idsa_Fe,N6,N5,N4))
    RFLHY=VFLW*AZMAX1(trcsa_soHml2(idsa_Hp,N6,N5,N4))
    RFLCA=VFLW*AZMAX1(trcsa_soHml2(idsa_Ca,N6,N5,N4))
    RFLMG=VFLW*AZMAX1(trcsa_soHml2(idsa_Mg,N6,N5,N4))
    RFLNA=VFLW*AZMAX1(trcsa_soHml2(idsa_Na,N6,N5,N4))
    RFLKA=VFLW*AZMAX1(trcsa_soHml2(idsa_K,N6,N5,N4))
    RFLOH=VFLW*AZMAX1(trcsa_soHml2(idsa_OH,N6,N5,N4))
    RFLSO=VFLW*AZMAX1(trcsa_soHml2(idsa_SO4,N6,N5,N4))
    RFLCL=VFLW*AZMAX1(trcsa_soHml2(idsa_Cl,N6,N5,N4))
    RFLC3=VFLW*AZMAX1(trcsa_soHml2(idsa_CO3,N6,N5,N4))
    RFLHC=VFLW*AZMAX1(trcsa_soHml2(idsa_HCO3,N6,N5,N4))
    RFLAL1=VFLW*AZMAX1(trcsa_soHml2(idsa_AlOH,N6,N5,N4))
    RFLAL2=VFLW*AZMAX1(trcsa_soHml2(idsa_AlOH2,N6,N5,N4))
    RFLAL3=VFLW*AZMAX1(trcsa_soHml2(idsa_AlOH3,N6,N5,N4))
    RFLAL4=VFLW*AZMAX1(trcsa_soHml2(idsa_AlOH4,N6,N5,N4))
    RFLALS=VFLW*AZMAX1(trcsa_soHml2(idsa_AlSO4,N6,N5,N4))
    RFLFE1=VFLW*AZMAX1(trcsa_soHml2(idsa_FeOH,N6,N5,N4))
    RFLFE2=VFLW*AZMAX1(trcsa_soHml2(idsa_FeOH2,N6,N5,N4))
    RFLFE3=VFLW*AZMAX1(trcsa_soHml2(idsa_FeOH3,N6,N5,N4))
    RFLFE4=VFLW*AZMAX1(trcsa_soHml2(idsa_FeOH4,N6,N5,N4))
    RFLFES=VFLW*AZMAX1(trcsa_soHml2(idsa_FeSO4,N6,N5,N4))
    RFLCAO=VFLW*AZMAX1(trcsa_soHml2(idsa_CaOH2,N6,N5,N4))
    RFLCAC=VFLW*AZMAX1(trcsa_soHml2(idsa_CaCO3,N6,N5,N4))
    RFLCAH=VFLW*AZMAX1(trcsa_soHml2(idsa_CaHCO3,N6,N5,N4))
    RFLCAS=VFLW*AZMAX1(trcsa_soHml2(idsa_CaSO4,N6,N5,N4))
    RFLMGO=VFLW*AZMAX1(trcsa_soHml2(idsa_MgOH2,N6,N5,N4))
    RFLMGC=VFLW*AZMAX1(trcsa_soHml2(idsa_MgCO3,N6,N5,N4))
    RFLMGH=VFLW*AZMAX1(trcsa_soHml2(idsa_MgHCO3,N6,N5,N4))
    RFLMGS=VFLW*AZMAX1(trcsa_soHml2(idsa_MgSO4,N6,N5,N4))
    RFLNAC=VFLW*AZMAX1(trcsa_soHml2(idsa_NaCO3,N6,N5,N4))
    RFLNAS=VFLW*AZMAX1(trcsa_soHml2(idsa_NaSO4,N6,N5,N4))
    RFLKAS=VFLW*AZMAX1(trcsa_soHml2(idsa_KSO4,N6,N5,N4))
    RFLH0P=VFLW*AZMAX1(trcsa_soHml2(idsa_H0PO4,N6,N5,N4))*trcs_VLN(ids_H1PO4,N6,N5,N4)
    RFLH3P=VFLW*AZMAX1(trcsa_soHml2(idsa_H3PO4,N6,N5,N4))*trcs_VLN(ids_H1PO4,N6,N5,N4)
    RFLF1P=VFLW*AZMAX1(trcsa_soHml2(idsa_FeHPO4,N6,N5,N4))*trcs_VLN(ids_H1PO4,N6,N5,N4)
    RFLF2P=VFLW*AZMAX1(trcsa_soHml2(idsa_FeH2PO4,N6,N5,N4))*trcs_VLN(ids_H1PO4,N6,N5,N4)
    RFLC0P=VFLW*AZMAX1(trcsa_soHml2(idsa_CaPO4,N6,N5,N4))*trcs_VLN(ids_H1PO4,N6,N5,N4)
    RFLC1P=VFLW*AZMAX1(trcsa_soHml2(idsa_CaHPO4,N6,N5,N4))*trcs_VLN(ids_H1PO4,N6,N5,N4)
    RFLC2P=VFLW*AZMAX1(trcsa_soHml2(idsa_CaH2PO4,N6,N5,N4))*trcs_VLN(ids_H1PO4,N6,N5,N4)
    RFLM1P=VFLW*AZMAX1(trcsa_soHml2(idsa_MgHPO4,N6,N5,N4))*trcs_VLN(ids_H1PO4,N6,N5,N4)
    RFLH0B=VFLW*AZMAX1(trcsa_soHml2(idsa_H0PO4B,N6,N5,N4))*trcs_VLN(ids_H1PO4B,N6,N5,N4)
    RFLH3B=VFLW*AZMAX1(trcsa_soHml2(idsa_H3PO4B,N6,N5,N4))*trcs_VLN(ids_H1PO4B,N6,N5,N4)
    RFLF1B=VFLW*AZMAX1(trcsa_soHml2(idsa_FeHPO4B,N6,N5,N4))*trcs_VLN(ids_H1PO4B,N6,N5,N4)
    RFLF2B=VFLW*AZMAX1(trcsa_soHml2(idsa_FeH2PO4B,N6,N5,N4))*trcs_VLN(ids_H1PO4B,N6,N5,N4)
    RFLC0B=VFLW*AZMAX1(trcsa_soHml2(idsa_CaPO4B,N6,N5,N4))*trcs_VLN(ids_H1PO4B,N6,N5,N4)
    RFLC1B=VFLW*AZMAX1(trcsa_soHml2(idsa_CaHPO4B,N6,N5,N4))*trcs_VLN(ids_H1PO4B,N6,N5,N4)
    RFLC2B=VFLW*AZMAX1(trcsa_soHml2(idsa_CaH2PO4B,N6,N5,N4))*trcs_VLN(ids_H1PO4B,N6,N5,N4)
    RFLM1B=VFLW*AZMAX1(trcsa_soHml2(idsa_MgHPO4B,N6,N5,N4))*trcs_VLN(ids_H1PO4B,N6,N5,N4)
!
!     MICROPORE TO MACROPORE TRANSFER
!
  ELSEIF(FINHM(M,N6,N5,N4).LT.0.0)THEN
    IF(VOLWM(M,N6,N5,N4).GT.ZEROS2(NY,NX))THEN
      VFLW=AZMIN1(AMAX1(-VFLWX,FINHM(M,N6,N5,N4)/VOLWM(M,N6,N5,N4)))
    ELSE
      VFLW=-VFLWX
    ENDIF
    RFLAL=VFLW*AZMAX1(trcsa_solml2(idsa_Al,N6,N5,N4))
    RFLFE=VFLW*AZMAX1(trcsa_solml2(idsa_Fe,N6,N5,N4))
    RFLHY=VFLW*AZMAX1(trcsa_solml2(idsa_Hp,N6,N5,N4))
    RFLCA=VFLW*AZMAX1(trcsa_solml2(idsa_Ca,N6,N5,N4))
    RFLMG=VFLW*AZMAX1(trcsa_solml2(idsa_Mg,N6,N5,N4))
    RFLNA=VFLW*AZMAX1(trcsa_solml2(idsa_Na,N6,N5,N4))
    RFLKA=VFLW*AZMAX1(trcsa_solml2(idsa_K,N6,N5,N4))
    RFLOH=VFLW*AZMAX1(trcsa_solml2(idsa_OH,N6,N5,N4))
    RFLSO=VFLW*AZMAX1(trcsa_solml2(idsa_SO4,N6,N5,N4))
    RFLCL=VFLW*AZMAX1(trcsa_solml2(idsa_Cl,N6,N5,N4))
    RFLC3=VFLW*AZMAX1(trcsa_solml2(idsa_CO3,N6,N5,N4))
    RFLHC=VFLW*AZMAX1(trcsa_solml2(idsa_HCO3,N6,N5,N4))
    RFLAL1=VFLW*AZMAX1(trcsa_solml2(idsa_AlOH,N6,N5,N4))
    RFLAL2=VFLW*AZMAX1(trcsa_solml2(idsa_AlOH2,N6,N5,N4))
    RFLAL3=VFLW*AZMAX1(trcsa_solml2(idsa_AlOH3,N6,N5,N4))
    RFLAL4=VFLW*AZMAX1(trcsa_solml2(idsa_AlOH4,N6,N5,N4))
    RFLALS=VFLW*AZMAX1(trcsa_solml2(idsa_AlSO4,N6,N5,N4))
    RFLFE1=VFLW*AZMAX1(trcsa_solml2(idsa_FeOH,N6,N5,N4))
    RFLFE2=VFLW*AZMAX1(trcsa_solml2(idsa_FeOH2,N6,N5,N4))
    RFLFE3=VFLW*AZMAX1(trcsa_solml2(idsa_FeOH3,N6,N5,N4))
    RFLFE4=VFLW*AZMAX1(trcsa_solml2(idsa_FeOH4,N6,N5,N4))
    RFLFES=VFLW*AZMAX1(trcsa_solml2(idsa_FeSO4,N6,N5,N4))
    RFLCAO=VFLW*AZMAX1(trcsa_solml2(idsa_CaOH2,N6,N5,N4))
    RFLCAC=VFLW*AZMAX1(trcsa_solml2(idsa_CaCO3,N6,N5,N4))
    RFLCAH=VFLW*AZMAX1(trcsa_solml2(idsa_CaHCO3,N6,N5,N4))
    RFLCAS=VFLW*AZMAX1(trcsa_solml2(idsa_CaSO4,N6,N5,N4))
    RFLMGO=VFLW*AZMAX1(trcsa_solml2(idsa_MgOH2,N6,N5,N4))
    RFLMGC=VFLW*AZMAX1(trcsa_solml2(idsa_MgCO3,N6,N5,N4))
    RFLMGH=VFLW*AZMAX1(trcsa_solml2(idsa_MgHCO3,N6,N5,N4))
    RFLMGS=VFLW*AZMAX1(trcsa_solml2(idsa_MgSO4,N6,N5,N4))
    RFLNAC=VFLW*AZMAX1(trcsa_solml2(idsa_NaCO3,N6,N5,N4))
    RFLNAS=VFLW*AZMAX1(trcsa_solml2(idsa_NaSO4,N6,N5,N4))
    RFLKAS=VFLW*AZMAX1(trcsa_solml2(idsa_KSO4,N6,N5,N4))

    RFLH0P=VFLW*AZMAX1(trcsa_solml2(idsa_H0PO4,N6,N5,N4))*trcs_VLN(ids_H1PO4,N6,N5,N4)
    RFLH3P=VFLW*AZMAX1(trcsa_solml2(idsa_H3PO4,N6,N5,N4))*trcs_VLN(ids_H1PO4,N6,N5,N4)
    RFLF1P=VFLW*AZMAX1(trcsa_solml2(idsa_FeHPO4,N6,N5,N4))*trcs_VLN(ids_H1PO4,N6,N5,N4)
    RFLF2P=VFLW*AZMAX1(trcsa_solml2(idsa_FeH2PO4,N6,N5,N4))*trcs_VLN(ids_H1PO4,N6,N5,N4)
    RFLC0P=VFLW*AZMAX1(trcsa_solml2(idsa_CaPO4,N6,N5,N4))*trcs_VLN(ids_H1PO4,N6,N5,N4)
    RFLC1P=VFLW*AZMAX1(trcsa_solml2(idsa_CaHPO4,N6,N5,N4))*trcs_VLN(ids_H1PO4,N6,N5,N4)
    RFLC2P=VFLW*AZMAX1(trcsa_solml2(idsa_CaH2PO4,N6,N5,N4))*trcs_VLN(ids_H1PO4,N6,N5,N4)
    RFLM1P=VFLW*AZMAX1(trcsa_solml2(idsa_MgHPO4,N6,N5,N4))*trcs_VLN(ids_H1PO4,N6,N5,N4)
    RFLH0B=VFLW*AZMAX1(trcsa_solml2(idsa_H0PO4B,N6,N5,N4))*trcs_VLN(ids_H1PO4B,N6,N5,N4)
    RFLH3B=VFLW*AZMAX1(trcsa_solml2(idsa_H3PO4B,N6,N5,N4))*trcs_VLN(ids_H1PO4B,N6,N5,N4)
    RFLF1B=VFLW*AZMAX1(trcsa_solml2(idsa_FeHPO4B,N6,N5,N4))*trcs_VLN(ids_H1PO4B,N6,N5,N4)
    RFLF2B=VFLW*AZMAX1(trcsa_solml2(idsa_FeH2PO4B,N6,N5,N4))*trcs_VLN(ids_H1PO4B,N6,N5,N4)
    RFLC0B=VFLW*AZMAX1(trcsa_solml2(idsa_CaPO4B,N6,N5,N4))*trcs_VLN(ids_H1PO4B,N6,N5,N4)
    RFLC1B=VFLW*AZMAX1(trcsa_solml2(idsa_CaHPO4B,N6,N5,N4))*trcs_VLN(ids_H1PO4B,N6,N5,N4)
    RFLC2B=VFLW*AZMAX1(trcsa_solml2(idsa_CaH2PO4B,N6,N5,N4))*trcs_VLN(ids_H1PO4B,N6,N5,N4)
    RFLM1B=VFLW*AZMAX1(trcsa_solml2(idsa_MgHPO4B,N6,N5,N4))*trcs_VLN(ids_H1PO4B,N6,N5,N4)
!
!     NO MACROPORE TO MICROPORE TRANSFER
!
  ELSE
    RFLAL=0.0_r8
    RFLFE=0.0_r8
    RFLHY=0.0_r8
    RFLCA=0.0_r8
    RFLMG=0.0_r8
    RFLNA=0.0_r8
    RFLKA=0.0_r8
    RFLOH=0.0_r8
    RFLSO=0.0_r8
    RFLCL=0.0_r8
    RFLC3=0.0_r8
    RFLHC=0.0_r8
    RFLAL1=0.0_r8
    RFLAL2=0.0_r8
    RFLAL3=0.0_r8
    RFLAL4=0.0_r8
    RFLALS=0.0_r8
    RFLFE1=0.0_r8
    RFLFE2=0.0_r8
    RFLFE3=0.0_r8
    RFLFE4=0.0_r8
    RFLFES=0.0_r8
    RFLCAO=0.0_r8
    RFLCAC=0.0_r8
    RFLCAH=0.0_r8
    RFLCAS=0.0_r8
    RFLMGO=0.0_r8
    RFLMGC=0.0_r8
    RFLMGH=0.0_r8
    RFLMGS=0.0_r8
    RFLNAC=0.0_r8
    RFLNAS=0.0_r8
    RFLKAS=0.0_r8
    RFLH0P=0.0_r8
    RFLH3P=0.0_r8
    RFLF1P=0.0_r8
    RFLF2P=0.0_r8
    RFLC0P=0.0_r8
    RFLC1P=0.0_r8
    RFLC2P=0.0_r8
    RFLM1P=0.0_r8
    RFLH0B=0.0_r8
    RFLH3B=0.0_r8
    RFLF1B=0.0_r8
    RFLF2B=0.0_r8
    RFLC0B=0.0_r8
    RFLC1B=0.0_r8
    RFLC2B=0.0_r8
    RFLM1B=0.0_r8
  ENDIF
  trcsa_RFXS(idsa_Al,N6,N5,N4)=RFLAL
  trcsa_RFXS(idsa_Fe,N6,N5,N4)=RFLFE
  trcsa_RFXS(idsa_Hp,N6,N5,N4)=RFLHY
  trcsa_RFXS(idsa_Ca,N6,N5,N4)=RFLCA
  trcsa_RFXS(idsa_Mg,N6,N5,N4)=RFLMG
  trcsa_RFXS(idsa_Na,N6,N5,N4)=RFLNA
  trcsa_RFXS(idsa_K,N6,N5,N4)=RFLKA
  trcsa_RFXS(idsa_OH,N6,N5,N4)=RFLOH
  trcsa_RFXS(idsa_SO4,N6,N5,N4)=RFLSO
  trcsa_RFXS(idsa_Cl,N6,N5,N4)=RFLCL
  trcsa_RFXS(idsa_CO3,N6,N5,N4)=RFLC3
  trcsa_RFXS(idsa_HCO3,N6,N5,N4)=RFLHC
  trcsa_RFXS(idsa_AlOH,N6,N5,N4)=RFLAL1
  trcsa_RFXS(idsa_AlOH2,N6,N5,N4)=RFLAL2
  trcsa_RFXS(idsa_AlOH3,N6,N5,N4)=RFLAL3
  trcsa_RFXS(idsa_AlOH4,N6,N5,N4)=RFLAL4
  trcsa_RFXS(idsa_AlSO4,N6,N5,N4)=RFLALS
  trcsa_RFXS(idsa_FeOH,N6,N5,N4)=RFLFE1
  trcsa_RFXS(idsa_FeOH2,N6,N5,N4)=RFLFE2
  trcsa_RFXS(idsa_FeOH3,N6,N5,N4)=RFLFE3
  trcsa_RFXS(idsa_FeOH4,N6,N5,N4)=RFLFE4
  trcsa_RFXS(idsa_FeSO4,N6,N5,N4)=RFLFES
  trcsa_RFXS(idsa_CaOH2,N6,N5,N4)=RFLCAO
  trcsa_RFXS(idsa_CaCO3,N6,N5,N4)=RFLCAC
  trcsa_RFXS(idsa_CaHCO3,N6,N5,N4)=RFLCAH
  trcsa_RFXS(idsa_CaSO4,N6,N5,N4)=RFLCAS
  trcsa_RFXS(idsa_MgOH2,N6,N5,N4)=RFLMGO
  trcsa_RFXS(idsa_MgCO3,N6,N5,N4)=RFLMGC
  trcsa_RFXS(idsa_MgHCO3,N6,N5,N4)=RFLMGH
  trcsa_RFXS(idsa_MgSO4,N6,N5,N4)=RFLMGS
  trcsa_RFXS(idsa_NaCO3,N6,N5,N4)=RFLNAC
  trcsa_RFXS(idsa_NaSO4,N6,N5,N4)=RFLNAS
  trcsa_RFXS(idsa_KSO4,N6,N5,N4)=RFLKAS
  trcsa_RFXS(idsa_H0PO4,N6,N5,N4)=RFLH0P
  trcsa_RFXS(idsa_H3PO4,N6,N5,N4)=RFLH3P
  trcsa_RFXS(idsa_FeHPO4,N6,N5,N4)=RFLF1P
  trcsa_RFXS(idsa_FeH2PO4,N6,N5,N4)=RFLF2P
  trcsa_RFXS(idsa_CaPO4,N6,N5,N4)=RFLC0P
  trcsa_RFXS(idsa_CaHPO4,N6,N5,N4)=RFLC1P
  trcsa_RFXS(idsa_CaH2PO4,N6,N5,N4)=RFLC2P
  trcsa_RFXS(idsa_MgHPO4,N6,N5,N4)=RFLM1P
  trcsa_RFXS(idsa_H0PO4B,N6,N5,N4)=RFLH0B
  trcsa_RFXS(idsa_H3PO4B,N6,N5,N4)=RFLH3B
  trcsa_RFXS(idsa_FeHPO4B,N6,N5,N4)=RFLF1B
  trcsa_RFXS(idsa_FeH2PO4B,N6,N5,N4)=RFLF2B
  trcsa_RFXS(idsa_CaPO4B,N6,N5,N4)=RFLC0B
  trcsa_RFXS(idsa_CaHPO4B,N6,N5,N4)=RFLC1B
  trcsa_RFXS(idsa_CaH2PO4B,N6,N5,N4)=RFLC2B
  trcsa_RFXS(idsa_MgHPO4B,N6,N5,N4)=RFLM1B
  end subroutine SoluteAdvExchMicMacpores
!----------------------------------------------------------------------
  subroutine SoluteDifsExchMicMacpores(M,N,NY,NX,N1,N2,N3,N4,N5,N6)
  implicit none
  integer, intent(in) :: M,N,NY,NX,N1,N2,N3,N4,N5,N6
  real(r8) :: VOLWHS,VOLWT

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
    DFVAL=XNPH*(AZMAX1(trcsa_soHml2(idsa_Al,N6,N5,N4))*VOLWM(M,N6,N5,N4)-AZMAX1(trcsa_solml2(idsa_Al,N6,N5,N4))*VOLWHS)/VOLWT
    DFVFE=XNPH*(AZMAX1(trcsa_soHml2(idsa_Fe,N6,N5,N4))*VOLWM(M,N6,N5,N4)-AZMAX1(trcsa_solml2(idsa_Fe,N6,N5,N4))*VOLWHS)/VOLWT
    DFVHY=XNPH*(AZMAX1(trcsa_soHml2(idsa_Hp,N6,N5,N4))*VOLWM(M,N6,N5,N4)-AZMAX1(trcsa_solml2(idsa_Hp,N6,N5,N4))*VOLWHS)/VOLWT
    DFVCA=XNPH*(AZMAX1(trcsa_soHml2(idsa_Ca,N6,N5,N4))*VOLWM(M,N6,N5,N4)-AZMAX1(trcsa_solml2(idsa_Ca,N6,N5,N4))*VOLWHS)/VOLWT
    DFVMG=XNPH*(AZMAX1(trcsa_soHml2(idsa_Mg,N6,N5,N4))*VOLWM(M,N6,N5,N4)-AZMAX1(trcsa_solml2(idsa_Mg,N6,N5,N4))*VOLWHS)/VOLWT
    DFVNA=XNPH*(AZMAX1(trcsa_soHml2(idsa_Na,N6,N5,N4))*VOLWM(M,N6,N5,N4)-AZMAX1(trcsa_solml2(idsa_Na,N6,N5,N4))*VOLWHS)/VOLWT
    DFVKA=XNPH*(AZMAX1(trcsa_soHml2(idsa_K,N6,N5,N4))*VOLWM(M,N6,N5,N4)-AZMAX1(trcsa_solml2(idsa_K,N6,N5,N4))*VOLWHS)/VOLWT
    DFVOH=XNPH*(AZMAX1(trcsa_soHml2(idsa_OH,N6,N5,N4))*VOLWM(M,N6,N5,N4)-AZMAX1(trcsa_solml2(idsa_OH,N6,N5,N4))*VOLWHS)/VOLWT
    DFVSO=XNPH*(AZMAX1(trcsa_soHml2(idsa_SO4,N6,N5,N4))*VOLWM(M,N6,N5,N4)-AZMAX1(trcsa_solml2(idsa_SO4,N6,N5,N4))*VOLWHS)/VOLWT
    DFVCL=XNPH*(AZMAX1(trcsa_soHml2(idsa_Cl,N6,N5,N4))*VOLWM(M,N6,N5,N4)-AZMAX1(trcsa_solml2(idsa_Cl,N6,N5,N4))*VOLWHS)/VOLWT
    DFVC3=XNPH*(AZMAX1(trcsa_soHml2(idsa_CO3,N6,N5,N4))*VOLWM(M,N6,N5,N4)-AZMAX1(trcsa_solml2(idsa_CO3,N6,N5,N4))*VOLWHS)/VOLWT
    DFVHC=XNPH*(AZMAX1(trcsa_soHml2(idsa_HCO3,N6,N5,N4))*VOLWM(M,N6,N5,N4)-AZMAX1(trcsa_solml2(idsa_HCO3,N6,N5,N4))*VOLWHS)/VOLWT
    DFVAL1=XNPH*(AZMAX1(trcsa_soHml2(idsa_AlOH,N6,N5,N4))*VOLWM(M,N6,N5,N4)-AZMAX1(trcsa_solml2(idsa_AlOH,N6,N5,N4))*VOLWHS)/VOLWT
    DFVAL2=XNPH*(AZMAX1(trcsa_soHml2(idsa_AlOH2,N6,N5,N4))*VOLWM(M,N6,N5,N4)-AZMAX1(trcsa_solml2(idsa_AlOH2,N6,N5,N4))*VOLWHS)/VOLWT
    DFVAL3=XNPH*(AZMAX1(trcsa_soHml2(idsa_AlOH3,N6,N5,N4))*VOLWM(M,N6,N5,N4)-AZMAX1(trcsa_solml2(idsa_AlOH3,N6,N5,N4))*VOLWHS)/VOLWT
    DFVAL4=XNPH*(AZMAX1(trcsa_soHml2(idsa_AlOH4,N6,N5,N4))*VOLWM(M,N6,N5,N4)-AZMAX1(trcsa_solml2(idsa_AlOH4,N6,N5,N4))*VOLWHS)/VOLWT
    DFVALS=XNPH*(AZMAX1(trcsa_soHml2(idsa_AlSO4,N6,N5,N4))*VOLWM(M,N6,N5,N4)-AZMAX1(trcsa_solml2(idsa_AlSO4,N6,N5,N4))*VOLWHS)/VOLWT
    DFVFE1=XNPH*(AZMAX1(trcsa_soHml2(idsa_FeOH,N6,N5,N4))*VOLWM(M,N6,N5,N4)-AZMAX1(trcsa_solml2(idsa_FeOH,N6,N5,N4))*VOLWHS)/VOLWT
    DFVF22=XNPH*(AZMAX1(trcsa_soHml2(idsa_FeOH2,N6,N5,N4))*VOLWM(M,N6,N5,N4)-AZMAX1(trcsa_solml2(idsa_FeOH2,N6,N5,N4))*VOLWHS)/VOLWT
    DFVFE3=XNPH*(AZMAX1(trcsa_soHml2(idsa_FeOH3,N6,N5,N4))*VOLWM(M,N6,N5,N4)-AZMAX1(trcsa_solml2(idsa_FeOH3,N6,N5,N4))*VOLWHS)/VOLWT
    DFVFE4=XNPH*(AZMAX1(trcsa_soHml2(idsa_FeOH4,N6,N5,N4))*VOLWM(M,N6,N5,N4)-AZMAX1(trcsa_solml2(idsa_FeOH4,N6,N5,N4))*VOLWHS)/VOLWT
    DFVFES=XNPH*(AZMAX1(trcsa_soHml2(idsa_FeSO4,N6,N5,N4))*VOLWM(M,N6,N5,N4)-AZMAX1(trcsa_solml2(idsa_FeSO4,N6,N5,N4))*VOLWHS)/VOLWT
    DFVCAO=XNPH*(AZMAX1(trcsa_soHml2(idsa_CaOH2,N6,N5,N4))*VOLWM(M,N6,N5,N4)-AZMAX1(trcsa_solml2(idsa_CaOH2,N6,N5,N4))*VOLWHS)/VOLWT
    DFVCAC=XNPH*(AZMAX1(trcsa_soHml2(idsa_CaCO3,N6,N5,N4))*VOLWM(M,N6,N5,N4)-AZMAX1(trcsa_solml2(idsa_CaCO3,N6,N5,N4))*VOLWHS)/VOLWT
    DFVCAH=XNPH*(AZMAX1(trcsa_soHml2(idsa_CaHCO3,N6,N5,N4))*VOLWM(M,N6,N5,N4)-AZMAX1(trcsa_solml2(idsa_CaHCO3,N6,N5,N4))*VOLWHS)/VOLWT
    DFVCAS=XNPH*(AZMAX1(trcsa_soHml2(idsa_CaSO4,N6,N5,N4))*VOLWM(M,N6,N5,N4)-AZMAX1(trcsa_solml2(idsa_CaSO4,N6,N5,N4))*VOLWHS)/VOLWT
    DFVMGO=XNPH*(AZMAX1(trcsa_soHml2(idsa_MgOH2,N6,N5,N4))*VOLWM(M,N6,N5,N4)-AZMAX1(trcsa_solml2(idsa_MgOH2,N6,N5,N4))*VOLWHS)/VOLWT
    DFVMGC=XNPH*(AZMAX1(trcsa_soHml2(idsa_MgCO3,N6,N5,N4))*VOLWM(M,N6,N5,N4)-AZMAX1(trcsa_solml2(idsa_MgCO3,N6,N5,N4))*VOLWHS)/VOLWT
    DFVMGH=XNPH*(AZMAX1(trcsa_soHml2(idsa_MgHCO3,N6,N5,N4))*VOLWM(M,N6,N5,N4)-AZMAX1(trcsa_solml2(idsa_MgHCO3,N6,N5,N4))*VOLWHS)/VOLWT
    DFVMGS=XNPH*(AZMAX1(trcsa_soHml2(idsa_MgSO4,N6,N5,N4))*VOLWM(M,N6,N5,N4)-AZMAX1(trcsa_solml2(idsa_MgSO4,N6,N5,N4))*VOLWHS)/VOLWT
    DFVNAC=XNPH*(AZMAX1(trcsa_soHml2(idsa_NaCO3,N6,N5,N4))*VOLWM(M,N6,N5,N4)-AZMAX1(trcsa_solml2(idsa_NaCO3,N6,N5,N4))*VOLWHS)/VOLWT
    DFVNAS=XNPH*(AZMAX1(trcsa_soHml2(idsa_NaSO4,N6,N5,N4))*VOLWM(M,N6,N5,N4)-AZMAX1(trcsa_solml2(idsa_NaSO4,N6,N5,N4))*VOLWHS)/VOLWT
    DFVKAS=XNPH*(AZMAX1(trcsa_soHml2(idsa_KSO4,N6,N5,N4))*VOLWM(M,N6,N5,N4)-AZMAX1(trcsa_solml2(idsa_KSO4,N6,N5,N4))*VOLWHS)/VOLWT
    DFVNAC=XNPH*(AZMAX1(trcsa_soHml2(idsa_NaCO3,N6,N5,N4))*VOLWM(M,N6,N5,N4)-AZMAX1(trcsa_solml2(idsa_NaCO3,N6,N5,N4))*VOLWHS)/VOLWT
    DFVH0P=XNPH*(AZMAX1(trcsa_soHml2(idsa_H0PO4,N6,N5,N4))*VOLWM(M,N6,N5,N4)-AZMAX1(trcsa_solml2(idsa_H0PO4,N6,N5,N4))*VOLWHS)/VOLWT*trcs_VLN(ids_H1PO4,N6,N5,N4)
    DFVH3P=XNPH*(AZMAX1(trcsa_soHml2(idsa_H3PO4,N6,N5,N4))*VOLWM(M,N6,N5,N4)-AZMAX1(trcsa_solml2(idsa_H3PO4,N6,N5,N4))*VOLWHS)/VOLWT*trcs_VLN(ids_H1PO4,N6,N5,N4)
    DFVF1P=XNPH*(AZMAX1(trcsa_soHml2(idsa_FeHPO4,N6,N5,N4))*VOLWM(M,N6,N5,N4)-AZMAX1(trcsa_solml2(idsa_FeHPO4,N6,N5,N4))*VOLWHS)/VOLWT*trcs_VLN(ids_H1PO4,N6,N5,N4)
    DFVF2P=XNPH*(AZMAX1(trcsa_soHml2(idsa_FeH2PO4,N6,N5,N4))*VOLWM(M,N6,N5,N4)-AZMAX1(trcsa_solml2(idsa_FeH2PO4,N6,N5,N4))*VOLWHS)/VOLWT*trcs_VLN(ids_H1PO4,N6,N5,N4)
    DFVC0P=XNPH*(AZMAX1(trcsa_soHml2(idsa_CaPO4,N6,N5,N4))*VOLWM(M,N6,N5,N4)-AZMAX1(trcsa_solml2(idsa_CaPO4,N6,N5,N4))*VOLWHS)/VOLWT*trcs_VLN(ids_H1PO4,N6,N5,N4)
    DFVC1P=XNPH*(AZMAX1(trcsa_soHml2(idsa_CaHPO4,N6,N5,N4))*VOLWM(M,N6,N5,N4)-AZMAX1(trcsa_solml2(idsa_CaHPO4,N6,N5,N4))*VOLWHS)/VOLWT*trcs_VLN(ids_H1PO4,N6,N5,N4)
    DFVC2P=XNPH*(AZMAX1(trcsa_soHml2(idsa_CaH2PO4,N6,N5,N4))*VOLWM(M,N6,N5,N4)-AZMAX1(trcsa_solml2(idsa_CaH2PO4,N6,N5,N4))*VOLWHS)/VOLWT*trcs_VLN(ids_H1PO4,N6,N5,N4)
    DFVM1P=XNPH*(AZMAX1(trcsa_soHml2(idsa_MgHPO4,N6,N5,N4))*VOLWM(M,N6,N5,N4)-AZMAX1(trcsa_solml2(idsa_MgHPO4,N6,N5,N4))*VOLWHS)/VOLWT*trcs_VLN(ids_H1PO4,N6,N5,N4)
    DFVH0B=XNPH*(AZMAX1(trcsa_soHml2(idsa_H0PO4B,N6,N5,N4))*VOLWM(M,N6,N5,N4)-AZMAX1(trcsa_solml2(idsa_H0PO4B,N6,N5,N4))*VOLWHS)/VOLWT*trcs_VLN(ids_H1PO4B,N6,N5,N4)
    DFVH3B=XNPH*(AZMAX1(trcsa_soHml2(idsa_H3PO4B,N6,N5,N4))*VOLWM(M,N6,N5,N4)-AZMAX1(trcsa_solml2(idsa_H3PO4B,N6,N5,N4))*VOLWHS)/VOLWT*trcs_VLN(ids_H1PO4B,N6,N5,N4)
    DFVF1B=XNPH*(AZMAX1(trcsa_soHml2(idsa_FeHPO4B,N6,N5,N4))*VOLWM(M,N6,N5,N4)-AZMAX1(trcsa_solml2(idsa_FeHPO4B,N6,N5,N4))*VOLWHS)/VOLWT*trcs_VLN(ids_H1PO4B,N6,N5,N4)
    DFVF2B=XNPH*(AZMAX1(trcsa_soHml2(idsa_FeH2PO4B,N6,N5,N4))*VOLWM(M,N6,N5,N4)-AZMAX1(trcsa_solml2(idsa_FeH2PO4B,N6,N5,N4))*VOLWHS)/VOLWT*trcs_VLN(ids_H1PO4B,N6,N5,N4)
    DFVC0B=XNPH*(AZMAX1(trcsa_soHml2(idsa_CaPO4B,N6,N5,N4))*VOLWM(M,N6,N5,N4)-AZMAX1(trcsa_solml2(idsa_CaPO4B,N6,N5,N4))*VOLWHS)/VOLWT*trcs_VLN(ids_H1PO4B,N6,N5,N4)
    DFVC1B=XNPH*(AZMAX1(trcsa_soHml2(idsa_CaHPO4B,N6,N5,N4))*VOLWM(M,N6,N5,N4)-AZMAX1(trcsa_solml2(idsa_CaHPO4B,N6,N5,N4))*VOLWHS)/VOLWT*trcs_VLN(ids_H1PO4B,N6,N5,N4)
    DFVC2B=XNPH*(AZMAX1(trcsa_soHml2(idsa_CaH2PO4B,N6,N5,N4))*VOLWM(M,N6,N5,N4)-AZMAX1(trcsa_solml2(idsa_CaH2PO4B,N6,N5,N4))*VOLWHS)/VOLWT*trcs_VLN(ids_H1PO4B,N6,N5,N4)
    DFVM1B=XNPH*(AZMAX1(trcsa_soHml2(idsa_MgHPO4B,N6,N5,N4))*VOLWM(M,N6,N5,N4)-AZMAX1(trcsa_solml2(idsa_MgHPO4B,N6,N5,N4))*VOLWHS)/VOLWT*trcs_VLN(ids_H1PO4B,N6,N5,N4)
  ELSE
    DFVAL=0.0_r8
    DFVFE=0.0_r8
    DFVHY=0.0_r8
    DFVCA=0.0_r8
    DFVMG=0.0_r8
    DFVNA=0.0_r8
    DFVKA=0.0_r8
    DFVOH=0.0_r8
    DFVSO=0.0_r8
    DFVCL=0.0_r8
    DFVC3=0.0_r8
    DFVHC=0.0_r8
    DFVAL1=0.0_r8
    DFVAL2=0.0_r8
    DFVAL3=0.0_r8
    DFVAL4=0.0_r8
    DFVALS=0.0_r8
    DFVFE1=0.0_r8
    DFVFE2=0.0_r8
    DFVFE3=0.0_r8
    DFVFE4=0.0_r8
    DFVFES=0.0_r8
    DFVCAO=0.0_r8
    DFVCAC=0.0_r8
    DFVCAH=0.0_r8
    DFVCAS=0.0_r8
    DFVMGO=0.0_r8
    DFVMGC=0.0_r8
    DFVMGH=0.0_r8
    DFVMGS=0.0_r8
    DFVNAC=0.0_r8
    DFVNAS=0.0_r8
    DFVKAS=0.0_r8
    DFVH0P=0.0_r8
    DFVH3P=0.0_r8
    DFVF1P=0.0_r8
    DFVF2P=0.0_r8
    DFVC0P=0.0_r8
    DFVC1P=0.0_r8
    DFVC2P=0.0_r8
    DFVM1P=0.0_r8
    DFVH0B=0.0_r8
    DFVH3B=0.0_r8
    DFVF1B=0.0_r8
    DFVF2B=0.0_r8
    DFVC0B=0.0_r8
    DFVC1B=0.0_r8
    DFVC2B=0.0_r8
    DFVM1B=0.0_r8
  ENDIF
  trcsa_RFXS(idsa_Al,N6,N5,N4)=trcsa_RFXS(idsa_Al,N6,N5,N4)+DFVAL
  trcsa_RFXS(idsa_Fe,N6,N5,N4)=trcsa_RFXS(idsa_Fe,N6,N5,N4)+DFVFE
  trcsa_RFXS(idsa_Hp,N6,N5,N4)=trcsa_RFXS(idsa_Hp,N6,N5,N4)+DFVHY
  trcsa_RFXS(idsa_Ca,N6,N5,N4)=trcsa_RFXS(idsa_Ca,N6,N5,N4)+DFVCA
  trcsa_RFXS(idsa_Mg,N6,N5,N4)=trcsa_RFXS(idsa_Mg,N6,N5,N4)+DFVMG
  trcsa_RFXS(idsa_Na,N6,N5,N4)=trcsa_RFXS(idsa_Na,N6,N5,N4)+DFVNA
  trcsa_RFXS(idsa_K,N6,N5,N4)=trcsa_RFXS(idsa_K,N6,N5,N4)+DFVKA
  trcsa_RFXS(idsa_OH,N6,N5,N4)=trcsa_RFXS(idsa_OH,N6,N5,N4)+DFVOH
  trcsa_RFXS(idsa_SO4,N6,N5,N4)=trcsa_RFXS(idsa_SO4,N6,N5,N4)+DFVSO
  trcsa_RFXS(idsa_Cl,N6,N5,N4)=trcsa_RFXS(idsa_Cl,N6,N5,N4)+DFVCL
  trcsa_RFXS(idsa_CO3,N6,N5,N4)=trcsa_RFXS(idsa_CO3,N6,N5,N4)+DFVC3
  trcsa_RFXS(idsa_HCO3,N6,N5,N4)=trcsa_RFXS(idsa_HCO3,N6,N5,N4)+DFVHC
  trcsa_RFXS(idsa_AlOH,N6,N5,N4)=trcsa_RFXS(idsa_AlOH,N6,N5,N4)+DFVAL1
  trcsa_RFXS(idsa_AlOH2,N6,N5,N4)=trcsa_RFXS(idsa_AlOH2,N6,N5,N4)+DFVAL2
  trcsa_RFXS(idsa_AlOH3,N6,N5,N4)=trcsa_RFXS(idsa_AlOH3,N6,N5,N4)+DFVAL3
  trcsa_RFXS(idsa_AlOH4,N6,N5,N4)=trcsa_RFXS(idsa_AlOH4,N6,N5,N4)+DFVAL4
  trcsa_RFXS(idsa_AlSO4,N6,N5,N4)=trcsa_RFXS(idsa_AlSO4,N6,N5,N4)+DFVALS
  trcsa_RFXS(idsa_FeOH,N6,N5,N4)=trcsa_RFXS(idsa_FeOH,N6,N5,N4)+DFVFE1
  trcsa_RFXS(idsa_FeOH2,N6,N5,N4)=trcsa_RFXS(idsa_FeOH2,N6,N5,N4)+DFVFE2
  trcsa_RFXS(idsa_FeOH3,N6,N5,N4)=trcsa_RFXS(idsa_FeOH3,N6,N5,N4)+DFVFE3
  trcsa_RFXS(idsa_FeOH4,N6,N5,N4)=trcsa_RFXS(idsa_FeOH4,N6,N5,N4)+DFVFE4
  trcsa_RFXS(idsa_FeSO4,N6,N5,N4)=trcsa_RFXS(idsa_FeSO4,N6,N5,N4)+DFVFES
  trcsa_RFXS(idsa_CaOH2,N6,N5,N4)=trcsa_RFXS(idsa_CaOH2,N6,N5,N4)+DFVCAO
  trcsa_RFXS(idsa_CaCO3,N6,N5,N4)=trcsa_RFXS(idsa_CaCO3,N6,N5,N4)+DFVCAC
  trcsa_RFXS(idsa_CaHCO3,N6,N5,N4)=trcsa_RFXS(idsa_CaHCO3,N6,N5,N4)+DFVCAH
  trcsa_RFXS(idsa_CaSO4,N6,N5,N4)=trcsa_RFXS(idsa_CaSO4,N6,N5,N4)+DFVCAS
  trcsa_RFXS(idsa_MgOH2,N6,N5,N4)=trcsa_RFXS(idsa_MgOH2,N6,N5,N4)+DFVMGO
  trcsa_RFXS(idsa_MgCO3,N6,N5,N4)=trcsa_RFXS(idsa_MgCO3,N6,N5,N4)+DFVMGC
  trcsa_RFXS(idsa_MgHCO3,N6,N5,N4)=trcsa_RFXS(idsa_MgHCO3,N6,N5,N4)+DFVMGH
  trcsa_RFXS(idsa_MgSO4,N6,N5,N4)=trcsa_RFXS(idsa_MgSO4,N6,N5,N4)+DFVMGS
  trcsa_RFXS(idsa_NaCO3,N6,N5,N4)=trcsa_RFXS(idsa_NaCO3,N6,N5,N4)+DFVNAC
  trcsa_RFXS(idsa_NaSO4,N6,N5,N4)=trcsa_RFXS(idsa_NaSO4,N6,N5,N4)+DFVNAS
  trcsa_RFXS(idsa_KSO4,N6,N5,N4)=trcsa_RFXS(idsa_KSO4,N6,N5,N4)+DFVKAS
  trcsa_RFXS(idsa_H0PO4,N6,N5,N4)=trcsa_RFXS(idsa_H0PO4,N6,N5,N4)+DFVH0P
  trcsa_RFXS(idsa_H3PO4,N6,N5,N4)=trcsa_RFXS(idsa_H3PO4,N6,N5,N4)+DFVH3P
  trcsa_RFXS(idsa_FeHPO4,N6,N5,N4)=trcsa_RFXS(idsa_FeHPO4,N6,N5,N4)+DFVF1P
  trcsa_RFXS(idsa_FeH2PO4,N6,N5,N4)=trcsa_RFXS(idsa_FeH2PO4,N6,N5,N4)+DFVF2P
  trcsa_RFXS(idsa_CaPO4,N6,N5,N4)=trcsa_RFXS(idsa_CaPO4,N6,N5,N4)+DFVC0P
  trcsa_RFXS(idsa_CaHPO4,N6,N5,N4)=trcsa_RFXS(idsa_CaHPO4,N6,N5,N4)+DFVC1P
  trcsa_RFXS(idsa_CaH2PO4,N6,N5,N4)=trcsa_RFXS(idsa_CaH2PO4,N6,N5,N4)+DFVC2P
  trcsa_RFXS(idsa_MgHPO4,N6,N5,N4)=trcsa_RFXS(idsa_MgHPO4,N6,N5,N4)+DFVM1P
  trcsa_RFXS(idsa_H0PO4B,N6,N5,N4)=trcsa_RFXS(idsa_H0PO4B,N6,N5,N4)+DFVH0B
  trcsa_RFXS(idsa_H3PO4B,N6,N5,N4)=trcsa_RFXS(idsa_H3PO4B,N6,N5,N4)+DFVH3B
  trcsa_RFXS(idsa_FeHPO4B,N6,N5,N4)=trcsa_RFXS(idsa_FeHPO4B,N6,N5,N4)+DFVF1B
  trcsa_RFXS(idsa_FeH2PO4B,N6,N5,N4)=trcsa_RFXS(idsa_FeH2PO4B,N6,N5,N4)+DFVF2B
  trcsa_RFXS(idsa_CaPO4B,N6,N5,N4)=trcsa_RFXS(idsa_CaPO4B,N6,N5,N4)+DFVC0B
  trcsa_RFXS(idsa_CaHPO4B,N6,N5,N4)=trcsa_RFXS(idsa_CaHPO4B,N6,N5,N4)+DFVC1B
  trcsa_RFXS(idsa_CaH2PO4B,N6,N5,N4)=trcsa_RFXS(idsa_CaH2PO4B,N6,N5,N4)+DFVC2B
  trcsa_RFXS(idsa_MgHPO4B,N6,N5,N4)=trcsa_RFXS(idsa_MgHPO4B,N6,N5,N4)+DFVM1B
  end subroutine SoluteDifsExchMicMacpores
!----------------------------------------------------------------------
  subroutine SoluteAdvDifsMicMacpore(M,N,N1,N2,N3,N4,N5,N6,THETW1)
  implicit none
  integer, intent(in) :: M,N,N1,N2,N3,N4,N5,N6
  REAL(R8), INTENT(IN):: THETW1(JZ,JY,JX)

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
  trcsa_XFLS(idsa_Al,N,N6,N5,N4)=trcsa_XFLS(idsa_Al,N,N6,N5,N4)+trcsa_RFLS(idsa_Al,N,N6,N5,N4)
  trcsa_XFLS(idsa_Fe,N,N6,N5,N4)=trcsa_XFLS(idsa_Fe,N,N6,N5,N4)+trcsa_RFLS(idsa_Fe,N,N6,N5,N4)
  trcsa_XFLS(idsa_Hp,N,N6,N5,N4)=trcsa_XFLS(idsa_Hp,N,N6,N5,N4)+trcsa_RFLS(idsa_Hp,N,N6,N5,N4)
  trcsa_XFLS(idsa_Ca,N,N6,N5,N4)=trcsa_XFLS(idsa_Ca,N,N6,N5,N4)+trcsa_RFLS(idsa_Ca,N,N6,N5,N4)
  trcsa_XFLS(idsa_Mg,N,N6,N5,N4)=trcsa_XFLS(idsa_Mg,N,N6,N5,N4)+trcsa_RFLS(idsa_Mg,N,N6,N5,N4)
  trcsa_XFLS(idsa_Na,N,N6,N5,N4)=trcsa_XFLS(idsa_Na,N,N6,N5,N4)+trcsa_RFLS(idsa_Na,N,N6,N5,N4)
  trcsa_XFLS(idsa_K,N,N6,N5,N4)=trcsa_XFLS(idsa_K,N,N6,N5,N4)+trcsa_RFLS(idsa_K,N,N6,N5,N4)
  trcsa_XFLS(idsa_OH,N,N6,N5,N4)=trcsa_XFLS(idsa_OH,N,N6,N5,N4)+trcsa_RFLS(idsa_OH,N,N6,N5,N4)
  trcsa_XFLS(idsa_SO4,N,N6,N5,N4)=trcsa_XFLS(idsa_SO4,N,N6,N5,N4)+trcsa_RFLS(idsa_SO4,N,N6,N5,N4)
  trcsa_XFLS(idsa_Cl,N,N6,N5,N4)=trcsa_XFLS(idsa_Cl,N,N6,N5,N4)+trcsa_RFLS(idsa_Cl,N,N6,N5,N4)
  trcsa_XFLS(idsa_CO3,N,N6,N5,N4)=trcsa_XFLS(idsa_CO3,N,N6,N5,N4)+trcsa_RFLS(idsa_CO3,N,N6,N5,N4)
  trcsa_XFLS(idsa_HCO3,N,N6,N5,N4)=trcsa_XFLS(idsa_HCO3,N,N6,N5,N4)+trcsa_RFLS(idsa_HCO3,N,N6,N5,N4)
  trcsa_XFLS(idsa_AlOH,N,N6,N5,N4)=trcsa_XFLS(idsa_AlOH,N,N6,N5,N4)+trcsa_RFLS(idsa_AlOH,N,N6,N5,N4)
  trcsa_XFLS(idsa_AlOH2,N,N6,N5,N4)=trcsa_XFLS(idsa_AlOH2,N,N6,N5,N4)+trcsa_RFLS(idsa_AlOH2,N,N6,N5,N4)
  trcsa_XFLS(idsa_AlOH3,N,N6,N5,N4)=trcsa_XFLS(idsa_AlOH3,N,N6,N5,N4)+trcsa_RFLS(idsa_AlOH3,N,N6,N5,N4)
  trcsa_XFLS(idsa_AlOH4,N,N6,N5,N4)=trcsa_XFLS(idsa_AlOH4,N,N6,N5,N4)+trcsa_RFLS(idsa_AlOH4,N,N6,N5,N4)
  trcsa_XFLS(idsa_AlSO4,N,N6,N5,N4)=trcsa_XFLS(idsa_AlSO4,N,N6,N5,N4)+trcsa_RFLS(idsa_AlSO4,N,N6,N5,N4)
  trcsa_XFLS(idsa_FeOH,N,N6,N5,N4)=trcsa_XFLS(idsa_FeOH,N,N6,N5,N4)+trcsa_RFLS(idsa_FeOH,N,N6,N5,N4)
  trcsa_XFLS(idsa_FeOH2,N,N6,N5,N4)=trcsa_XFLS(idsa_FeOH2,N,N6,N5,N4)+trcsa_RFLS(idsa_FeOH2,N,N6,N5,N4)
  trcsa_XFLS(idsa_FeOH3,N,N6,N5,N4)=trcsa_XFLS(idsa_FeOH3,N,N6,N5,N4)+trcsa_RFLS(idsa_FeOH3,N,N6,N5,N4)
  trcsa_XFLS(idsa_FeOH4,N,N6,N5,N4)=trcsa_XFLS(idsa_FeOH4,N,N6,N5,N4)+trcsa_RFLS(idsa_FeOH4,N,N6,N5,N4)
  trcsa_XFLS(idsa_FeSO4,N,N6,N5,N4)=trcsa_XFLS(idsa_FeSO4,N,N6,N5,N4)+trcsa_RFLS(idsa_FeSO4,N,N6,N5,N4)
  trcsa_XFLS(idsa_CaOH2,N,N6,N5,N4)=trcsa_XFLS(idsa_CaOH2,N,N6,N5,N4)+trcsa_RFLS(idsa_CaOH2,N,N6,N5,N4)
  trcsa_XFLS(idsa_CaCO3,N,N6,N5,N4)=trcsa_XFLS(idsa_CaCO3,N,N6,N5,N4)+trcsa_RFLS(idsa_CaCO3,N,N6,N5,N4)
  trcsa_XFLS(idsa_CaHCO3,N,N6,N5,N4)=trcsa_XFLS(idsa_CaHCO3,N,N6,N5,N4)+trcsa_RFLS(idsa_CaHCO3,N,N6,N5,N4)
  trcsa_XFLS(idsa_CaSO4,N,N6,N5,N4)=trcsa_XFLS(idsa_CaSO4,N,N6,N5,N4)+trcsa_RFLS(idsa_CaSO4,N,N6,N5,N4)
  trcsa_XFLS(idsa_MgOH2,N,N6,N5,N4)=trcsa_XFLS(idsa_MgOH2,N,N6,N5,N4)+trcsa_RFLS(idsa_MgOH2,N,N6,N5,N4)
  trcsa_XFLS(idsa_MgCO3,N,N6,N5,N4)=trcsa_XFLS(idsa_MgCO3,N,N6,N5,N4)+trcsa_RFLS(idsa_MgCO3,N,N6,N5,N4)
  trcsa_XFLS(idsa_MgHCO3,N,N6,N5,N4)=trcsa_XFLS(idsa_MgHCO3,N,N6,N5,N4)+trcsa_RFLS(idsa_MgHCO3,N,N6,N5,N4)
  trcsa_XFLS(idsa_MgSO4,N,N6,N5,N4)=trcsa_XFLS(idsa_MgSO4,N,N6,N5,N4)+trcsa_RFLS(idsa_MgSO4,N,N6,N5,N4)
  trcsa_XFLS(idsa_NaCO3,N,N6,N5,N4)=trcsa_XFLS(idsa_NaCO3,N,N6,N5,N4)+trcsa_RFLS(idsa_NaCO3,N,N6,N5,N4)
  trcsa_XFLS(idsa_NaSO4,N,N6,N5,N4)=trcsa_XFLS(idsa_NaSO4,N,N6,N5,N4)+trcsa_RFLS(idsa_NaSO4,N,N6,N5,N4)
  trcsa_XFLS(idsa_KSO4,N,N6,N5,N4)=trcsa_XFLS(idsa_KSO4,N,N6,N5,N4)+trcsa_RFLS(idsa_KSO4,N,N6,N5,N4)
  trcsa_XFLS(idsa_H0PO4,N,N6,N5,N4)=trcsa_XFLS(idsa_H0PO4,N,N6,N5,N4)+trcsa_RFLS(idsa_H0PO4,N,N6,N5,N4)
  trcsa_XFLS(idsa_H3PO4,N,N6,N5,N4)=trcsa_XFLS(idsa_H3PO4,N,N6,N5,N4)+trcsa_RFLS(idsa_H3PO4,N,N6,N5,N4)
  trcsa_XFLS(idsa_FeHPO4,N,N6,N5,N4)=trcsa_XFLS(idsa_FeHPO4,N,N6,N5,N4)+trcsa_RFLS(idsa_FeHPO4,N,N6,N5,N4)
  trcsa_XFLS(idsa_FeH2PO4,N,N6,N5,N4)=trcsa_XFLS(idsa_FeH2PO4,N,N6,N5,N4)+trcsa_RFLS(idsa_FeH2PO4,N,N6,N5,N4)
  trcsa_XFLS(idsa_CaPO4,N,N6,N5,N4)=trcsa_XFLS(idsa_CaPO4,N,N6,N5,N4)+trcsa_RFLS(idsa_CaPO4,N,N6,N5,N4)
  trcsa_XFLS(idsa_CaHPO4,N,N6,N5,N4)=trcsa_XFLS(idsa_CaHPO4,N,N6,N5,N4)+trcsa_RFLS(idsa_CaHPO4,N,N6,N5,N4)
  trcsa_XFLS(idsa_CaH2PO4,N,N6,N5,N4)=trcsa_XFLS(idsa_CaH2PO4,N,N6,N5,N4)+trcsa_RFLS(idsa_CaH2PO4,N,N6,N5,N4)
  trcsa_XFLS(idsa_MgHPO4,N,N6,N5,N4)=trcsa_XFLS(idsa_MgHPO4,N,N6,N5,N4)+trcsa_RFLS(idsa_MgHPO4,N,N6,N5,N4)

  trcsa_XFLS(idsa_H0PO4B,N,N6,N5,N4)=trcsa_XFLS(idsa_H0PO4B,N,N6,N5,N4)+trcsa_RFLS(idsa_H0PO4B,N,N6,N5,N4)
  trcsa_XFLS(idsa_H3PO4B,N,N6,N5,N4)=trcsa_XFLS(idsa_H3PO4B,N,N6,N5,N4)+trcsa_RFLS(idsa_H3PO4B,N,N6,N5,N4)
  trcsa_XFLS(idsa_FeHPO4B,N,N6,N5,N4)=trcsa_XFLS(idsa_FeHPO4B,N,N6,N5,N4)+trcsa_RFLS(idsa_FeHPO4B,N,N6,N5,N4)
  trcsa_XFLS(idsa_FeH2PO4B,N,N6,N5,N4)=trcsa_XFLS(idsa_FeH2PO4B,N,N6,N5,N4)+trcsa_RFLS(idsa_FeH2PO4B,N,N6,N5,N4)
  trcsa_XFLS(idsa_CaPO4B,N,N6,N5,N4)=trcsa_XFLS(idsa_CaPO4B,N,N6,N5,N4)+trcsa_RFLS(idsa_CaPO4B,N,N6,N5,N4)
  trcsa_XFLS(idsa_CaHPO4B,N,N6,N5,N4)=trcsa_XFLS(idsa_CaHPO4B,N,N6,N5,N4)+trcsa_RFLS(idsa_CaHPO4B,N,N6,N5,N4)
  trcsa_XFLS(idsa_CaH2PO4B,N,N6,N5,N4)=trcsa_XFLS(idsa_CaH2PO4B,N,N6,N5,N4)+trcsa_RFLS(idsa_CaH2PO4B,N,N6,N5,N4)
  trcsa_XFLS(idsa_MgHPO4B,N,N6,N5,N4)=trcsa_XFLS(idsa_MgHPO4B,N,N6,N5,N4)+trcsa_RFLS(idsa_MgHPO4B,N,N6,N5,N4)

  trcsa_XFHS(idsa_Al,N,N6,N5,N4)=trcsa_XFHS(idsa_Al,N,N6,N5,N4)+RALFHS(N,N6,N5,N4)
  trcsa_XFHS(idsa_Fe,N,N6,N5,N4)=trcsa_XFHS(idsa_Fe,N,N6,N5,N4)+RFEFHS(N,N6,N5,N4)
  trcsa_XFHS(idsa_Hp,N,N6,N5,N4)=trcsa_XFHS(idsa_Hp,N,N6,N5,N4)+RHYFHS(N,N6,N5,N4)
  trcsa_XFHS(idsa_Ca,N,N6,N5,N4)=trcsa_XFHS(idsa_Ca,N,N6,N5,N4)+RCAFHS(N,N6,N5,N4)
  trcsa_XFHS(idsa_Mg,N,N6,N5,N4)=trcsa_XFHS(idsa_Mg,N,N6,N5,N4)+RMGFHS(N,N6,N5,N4)
  trcsa_XFHS(idsa_Na,N,N6,N5,N4)=trcsa_XFHS(idsa_Na,N,N6,N5,N4)+RNAFHS(N,N6,N5,N4)
  trcsa_XFHS(idsa_K,N,N6,N5,N4)=trcsa_XFHS(idsa_K,N,N6,N5,N4)+RKAFHS(N,N6,N5,N4)
  trcsa_XFHS(idsa_OH,N,N6,N5,N4)=trcsa_XFHS(idsa_OH,N,N6,N5,N4)+ROHFHS(N,N6,N5,N4)
  trcsa_XFHS(idsa_SO4,N,N6,N5,N4)=trcsa_XFHS(idsa_SO4,N,N6,N5,N4)+RSOFHS(N,N6,N5,N4)
  trcsa_XFHS(idsa_Cl,N,N6,N5,N4)=trcsa_XFHS(idsa_Cl,N,N6,N5,N4)+RCLFHS(N,N6,N5,N4)
  trcsa_XFHS(idsa_CO3,N,N6,N5,N4)=trcsa_XFHS(idsa_CO3,N,N6,N5,N4)+RC3FHS(N,N6,N5,N4)
  trcsa_XFHS(idsa_HCO3,N,N6,N5,N4)=trcsa_XFHS(idsa_HCO3,N,N6,N5,N4)+RHCFHS(N,N6,N5,N4)
  trcsa_XFHS(idsa_AlOH,N,N6,N5,N4)=trcsa_XFHS(idsa_AlOH,N,N6,N5,N4)+RAL1HS(N,N6,N5,N4)
  trcsa_XFHS(idsa_AlOH2,N,N6,N5,N4)=trcsa_XFHS(idsa_AlOH2,N,N6,N5,N4)+RAL2HS(N,N6,N5,N4)
  trcsa_XFHS(idsa_AlOH3,N,N6,N5,N4)=trcsa_XFHS(idsa_AlOH3,N,N6,N5,N4)+RAL3HS(N,N6,N5,N4)
  trcsa_XFHS(idsa_AlOH4,N,N6,N5,N4)=trcsa_XFHS(idsa_AlOH4,N,N6,N5,N4)+RAL4HS(N,N6,N5,N4)
  trcsa_XFHS(idsa_AlSO4,N,N6,N5,N4)=trcsa_XFHS(idsa_AlSO4,N,N6,N5,N4)+RALSHS(N,N6,N5,N4)
  trcsa_XFHS(idsa_FeOH,N,N6,N5,N4)=trcsa_XFHS(idsa_FeOH,N,N6,N5,N4)+RFE1HS(N,N6,N5,N4)
  trcsa_XFHS(idsa_FeOH2,N,N6,N5,N4)=trcsa_XFHS(idsa_FeOH2,N,N6,N5,N4)+RFE2HS(N,N6,N5,N4)
  trcsa_XFHS(idsa_FeOH3,N,N6,N5,N4)=trcsa_XFHS(idsa_FeOH3,N,N6,N5,N4)+RFE3HS(N,N6,N5,N4)
  trcsa_XFHS(idsa_FeOH4,N,N6,N5,N4)=trcsa_XFHS(idsa_FeOH4,N,N6,N5,N4)+RFE4HS(N,N6,N5,N4)
  trcsa_XFHS(idsa_FeSO4,N,N6,N5,N4)=trcsa_XFHS(idsa_FeSO4,N,N6,N5,N4)+RFESHS(N,N6,N5,N4)
  trcsa_XFHS(idsa_CaOH2,N,N6,N5,N4)=trcsa_XFHS(idsa_CaOH2,N,N6,N5,N4)+RCAOHS(N,N6,N5,N4)
  trcsa_XFHS(idsa_CaCO3,N,N6,N5,N4)=trcsa_XFHS(idsa_CaCO3,N,N6,N5,N4)+RCACHS(N,N6,N5,N4)
  trcsa_XFHS(idsa_CaHCO3,N,N6,N5,N4)=trcsa_XFHS(idsa_CaHCO3,N,N6,N5,N4)+RCAHHS(N,N6,N5,N4)
  trcsa_XFHS(idsa_CaSO4,N,N6,N5,N4)=trcsa_XFHS(idsa_CaSO4,N,N6,N5,N4)+RCASHS(N,N6,N5,N4)
  trcsa_XFHS(idsa_MgOH2,N,N6,N5,N4)=trcsa_XFHS(idsa_MgOH2,N,N6,N5,N4)+RMGOHS(N,N6,N5,N4)
  trcsa_XFHS(idsa_MgCO3,N,N6,N5,N4)=trcsa_XFHS(idsa_MgCO3,N,N6,N5,N4)+RMGCHS(N,N6,N5,N4)
  trcsa_XFHS(idsa_MgHCO3,N,N6,N5,N4)=trcsa_XFHS(idsa_MgHCO3,N,N6,N5,N4)+RMGHHS(N,N6,N5,N4)
  trcsa_XFHS(idsa_MgSO4,N,N6,N5,N4)=trcsa_XFHS(idsa_MgSO4,N,N6,N5,N4)+RMGSHS(N,N6,N5,N4)
  trcsa_XFHS(idsa_NaCO3,N,N6,N5,N4)=trcsa_XFHS(idsa_NaCO3,N,N6,N5,N4)+RNACHS(N,N6,N5,N4)
  trcsa_XFHS(idsa_NaSO4,N,N6,N5,N4)=trcsa_XFHS(idsa_NaSO4,N,N6,N5,N4)+RNASHS(N,N6,N5,N4)
  trcsa_XFHS(idsa_KSO4,N,N6,N5,N4)=trcsa_XFHS(idsa_KSO4,N,N6,N5,N4)+RKASHS(N,N6,N5,N4)
  trcsa_XFHS(idsa_H0PO4,N,N6,N5,N4)=trcsa_XFHS(idsa_H0PO4,N,N6,N5,N4)+RH0PHS(N,N6,N5,N4)
  trcsa_XFHS(idsa_H3PO4,N,N6,N5,N4)=trcsa_XFHS(idsa_H3PO4,N,N6,N5,N4)+RH3PHS(N,N6,N5,N4)
  trcsa_XFHS(idsa_FeHPO4,N,N6,N5,N4)=trcsa_XFHS(idsa_FeHPO4,N,N6,N5,N4)+RF1PHS(N,N6,N5,N4)
  trcsa_XFHS(idsa_FeH2PO4,N,N6,N5,N4)=trcsa_XFHS(idsa_FeH2PO4,N,N6,N5,N4)+RF2PHS(N,N6,N5,N4)
  trcsa_XFHS(idsa_CaPO4,N,N6,N5,N4)=trcsa_XFHS(idsa_CaPO4,N,N6,N5,N4)+RC0PHS(N,N6,N5,N4)
  trcsa_XFHS(idsa_CaHPO4,N,N6,N5,N4)=trcsa_XFHS(idsa_CaHPO4,N,N6,N5,N4)+RC1PHS(N,N6,N5,N4)
  trcsa_XFHS(idsa_CaH2PO4,N,N6,N5,N4)=trcsa_XFHS(idsa_CaH2PO4,N,N6,N5,N4)+RC2PHS(N,N6,N5,N4)
  trcsa_XFHS(idsa_MgHPO4,N,N6,N5,N4)=trcsa_XFHS(idsa_MgHPO4,N,N6,N5,N4)+RM1PHS(N,N6,N5,N4)

  trcsa_XFHS(idsa_H0PO4B,N,N6,N5,N4)=trcsa_XFHS(idsa_H0PO4B,N,N6,N5,N4)+RH0BHB(N,N6,N5,N4)
  trcsa_XFHS(idsa_H3PO4B,N,N6,N5,N4)=trcsa_XFHS(idsa_H3PO4B,N,N6,N5,N4)+RH3BHB(N,N6,N5,N4)
  trcsa_XFHS(idsa_FeHPO4B,N,N6,N5,N4)=trcsa_XFHS(idsa_FeHPO4B,N,N6,N5,N4)+RF1BHB(N,N6,N5,N4)
  trcsa_XFHS(idsa_FeH2PO4B,N,N6,N5,N4)=trcsa_XFHS(idsa_FeH2PO4B,N,N6,N5,N4)+RF2BHB(N,N6,N5,N4)
  trcsa_XFHS(idsa_CaPO4B,N,N6,N5,N4)=trcsa_XFHS(idsa_CaPO4B,N,N6,N5,N4)+RC0BHB(N,N6,N5,N4)
  trcsa_XFHS(idsa_CaHPO4B,N,N6,N5,N4)=trcsa_XFHS(idsa_CaHPO4B,N,N6,N5,N4)+RC1BHB(N,N6,N5,N4)
  trcsa_XFHS(idsa_CaH2PO4B,N,N6,N5,N4)=trcsa_XFHS(idsa_CaH2PO4B,N,N6,N5,N4)+RC2BHB(N,N6,N5,N4)
  trcsa_XFHS(idsa_MgHPO4B,N,N6,N5,N4)=trcsa_XFHS(idsa_MgHPO4B,N,N6,N5,N4)+RM1BHB(N,N6,N5,N4)
  end subroutine SoluteAdvDifsMicMacpore
!----------------------------------------------------------------------
  subroutine SoluteAdvDifsExchMicMacpore(M,N,NY,NX,N1,N2,N3,N4,N5,N6)
  implicit none
  integer, intent(in) :: M,N,NY,NX,N1,N2,N3,N4,N5,N6


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
  trcsa_XFXS(idsa_Al,N6,N5,N4)=trcsa_XFXS(idsa_Al,N6,N5,N4)+trcsa_RFXS(idsa_Al,N6,N5,N4)
  trcsa_XFXS(idsa_Fe,N6,N5,N4)=trcsa_XFXS(idsa_Fe,N6,N5,N4)+trcsa_RFXS(idsa_Fe,N6,N5,N4)
  trcsa_XFXS(idsa_Hp,N6,N5,N4)=trcsa_XFXS(idsa_Hp,N6,N5,N4)+trcsa_RFXS(idsa_Hp,N6,N5,N4)
  trcsa_XFXS(idsa_Ca,N6,N5,N4)=trcsa_XFXS(idsa_Ca,N6,N5,N4)+trcsa_RFXS(idsa_Ca,N6,N5,N4)
  trcsa_XFXS(idsa_Mg,N6,N5,N4)=trcsa_XFXS(idsa_Mg,N6,N5,N4)+trcsa_RFXS(idsa_Mg,N6,N5,N4)
  trcsa_XFXS(idsa_Na,N6,N5,N4)=trcsa_XFXS(idsa_Na,N6,N5,N4)+trcsa_RFXS(idsa_Na,N6,N5,N4)
  trcsa_XFXS(idsa_K,N6,N5,N4)=trcsa_XFXS(idsa_K,N6,N5,N4)+trcsa_RFXS(idsa_K,N6,N5,N4)
  trcsa_XFXS(idsa_OH,N6,N5,N4)=trcsa_XFXS(idsa_OH,N6,N5,N4)+trcsa_RFXS(idsa_OH,N6,N5,N4)
  trcsa_XFXS(idsa_SO4,N6,N5,N4)=trcsa_XFXS(idsa_SO4,N6,N5,N4)+trcsa_RFXS(idsa_SO4,N6,N5,N4)
  trcsa_XFXS(idsa_Cl,N6,N5,N4)=trcsa_XFXS(idsa_Cl,N6,N5,N4)+trcsa_RFXS(idsa_Cl,N6,N5,N4)
  trcsa_XFXS(idsa_CO3,N6,N5,N4)=trcsa_XFXS(idsa_CO3,N6,N5,N4)+trcsa_RFXS(idsa_CO3,N6,N5,N4)
  trcsa_XFXS(idsa_HCO3,N6,N5,N4)=trcsa_XFXS(idsa_HCO3,N6,N5,N4)+trcsa_RFXS(idsa_HCO3,N6,N5,N4)
  trcsa_XFXS(idsa_AlOH,N6,N5,N4)=trcsa_XFXS(idsa_AlOH,N6,N5,N4)+trcsa_RFXS(idsa_AlOH,N6,N5,N4)
  trcsa_XFXS(idsa_AlOH2,N6,N5,N4)=trcsa_XFXS(idsa_AlOH2,N6,N5,N4)+trcsa_RFXS(idsa_AlOH2,N6,N5,N4)
  trcsa_XFXS(idsa_AlOH3,N6,N5,N4)=trcsa_XFXS(idsa_AlOH3,N6,N5,N4)+trcsa_RFXS(idsa_AlOH3,N6,N5,N4)
  trcsa_XFXS(idsa_AlOH4,N6,N5,N4)=trcsa_XFXS(idsa_AlOH4,N6,N5,N4)+trcsa_RFXS(idsa_AlOH4,N6,N5,N4)
  trcsa_XFXS(idsa_AlSO4,N6,N5,N4)=trcsa_XFXS(idsa_AlSO4,N6,N5,N4)+trcsa_RFXS(idsa_AlSO4,N6,N5,N4)
  trcsa_XFXS(idsa_FeOH,N6,N5,N4)=trcsa_XFXS(idsa_FeOH,N6,N5,N4)+trcsa_RFXS(idsa_FeOH,N6,N5,N4)
  trcsa_XFXS(idsa_FeOH2,N6,N5,N4)=trcsa_XFXS(idsa_FeOH2,N6,N5,N4)+trcsa_RFXS(idsa_FeOH2,N6,N5,N4)
  trcsa_XFXS(idsa_FeOH3,N6,N5,N4)=trcsa_XFXS(idsa_FeOH3,N6,N5,N4)+trcsa_RFXS(idsa_FeOH3,N6,N5,N4)
  trcsa_XFXS(idsa_FeOH4,N6,N5,N4)=trcsa_XFXS(idsa_FeOH4,N6,N5,N4)+trcsa_RFXS(idsa_FeOH4,N6,N5,N4)
  trcsa_XFXS(idsa_FeSO4,N6,N5,N4)=trcsa_XFXS(idsa_FeSO4,N6,N5,N4)+trcsa_RFXS(idsa_FeSO4,N6,N5,N4)
  trcsa_XFXS(idsa_CaOH2,N6,N5,N4)=trcsa_XFXS(idsa_CaOH2,N6,N5,N4)+trcsa_RFXS(idsa_CaOH2,N6,N5,N4)
  trcsa_XFXS(idsa_CaCO3,N6,N5,N4)=trcsa_XFXS(idsa_CaCO3,N6,N5,N4)+trcsa_RFXS(idsa_CaCO3,N6,N5,N4)
  trcsa_XFXS(idsa_CaHCO3,N6,N5,N4)=trcsa_XFXS(idsa_CaHCO3,N6,N5,N4)+trcsa_RFXS(idsa_CaHCO3,N6,N5,N4)
  trcsa_XFXS(idsa_CaSO4,N6,N5,N4)=trcsa_XFXS(idsa_CaSO4,N6,N5,N4)+trcsa_RFXS(idsa_CaSO4,N6,N5,N4)
  trcsa_XFXS(idsa_MgOH2,N6,N5,N4)=trcsa_XFXS(idsa_MgOH2,N6,N5,N4)+trcsa_RFXS(idsa_MgOH2,N6,N5,N4)
  trcsa_XFXS(idsa_MgCO3,N6,N5,N4)=trcsa_XFXS(idsa_MgCO3,N6,N5,N4)+trcsa_RFXS(idsa_MgCO3,N6,N5,N4)
  trcsa_XFXS(idsa_MgHCO3,N6,N5,N4)=trcsa_XFXS(idsa_MgHCO3,N6,N5,N4)+trcsa_RFXS(idsa_MgHCO3,N6,N5,N4)
  trcsa_XFXS(idsa_MgSO4,N6,N5,N4)=trcsa_XFXS(idsa_MgSO4,N6,N5,N4)+trcsa_RFXS(idsa_MgSO4,N6,N5,N4)
  trcsa_XFXS(idsa_NaCO3,N6,N5,N4)=trcsa_XFXS(idsa_NaCO3,N6,N5,N4)+trcsa_RFXS(idsa_NaCO3,N6,N5,N4)
  trcsa_XFXS(idsa_NaSO4,N6,N5,N4)=trcsa_XFXS(idsa_NaSO4,N6,N5,N4)+trcsa_RFXS(idsa_NaSO4,N6,N5,N4)
  trcsa_XFXS(idsa_KSO4,N6,N5,N4)=trcsa_XFXS(idsa_KSO4,N6,N5,N4)+trcsa_RFXS(idsa_KSO4,N6,N5,N4)
  trcsa_XFXS(idsa_H0PO4,N6,N5,N4)=trcsa_XFXS(idsa_H0PO4,N6,N5,N4)+trcsa_RFXS(idsa_H0PO4,N6,N5,N4)
  trcsa_XFXS(idsa_H3PO4,N6,N5,N4)=trcsa_XFXS(idsa_H3PO4,N6,N5,N4)+trcsa_RFXS(idsa_H3PO4,N6,N5,N4)
  trcsa_XFXS(idsa_FeHPO4,N6,N5,N4)=trcsa_XFXS(idsa_FeHPO4,N6,N5,N4)+trcsa_RFXS(idsa_FeHPO4,N6,N5,N4)
  trcsa_XFXS(idsa_FeH2PO4,N6,N5,N4)=trcsa_XFXS(idsa_FeH2PO4,N6,N5,N4)+trcsa_RFXS(idsa_FeH2PO4,N6,N5,N4)
  trcsa_XFXS(idsa_CaPO4,N6,N5,N4)=trcsa_XFXS(idsa_CaPO4,N6,N5,N4)+trcsa_RFXS(idsa_CaPO4,N6,N5,N4)
  trcsa_XFXS(idsa_CaHPO4,N6,N5,N4)=trcsa_XFXS(idsa_CaHPO4,N6,N5,N4)+trcsa_RFXS(idsa_CaHPO4,N6,N5,N4)
  trcsa_XFXS(idsa_CaH2PO4,N6,N5,N4)=trcsa_XFXS(idsa_CaH2PO4,N6,N5,N4)+trcsa_RFXS(idsa_CaH2PO4,N6,N5,N4)
  trcsa_XFXS(idsa_MgHPO4,N6,N5,N4)=trcsa_XFXS(idsa_MgHPO4,N6,N5,N4)+trcsa_RFXS(idsa_MgHPO4,N6,N5,N4)
  trcsa_XFXS(idsa_H0PO4B,N6,N5,N4)=trcsa_XFXS(idsa_H0PO4B,N6,N5,N4)+trcsa_RFXS(idsa_H0PO4B,N6,N5,N4)
  trcsa_XFXS(idsa_H3PO4B,N6,N5,N4)=trcsa_XFXS(idsa_H3PO4B,N6,N5,N4)+trcsa_RFXS(idsa_H3PO4B,N6,N5,N4)
  trcsa_XFXS(idsa_FeHPO4B,N6,N5,N4)=trcsa_XFXS(idsa_FeHPO4B,N6,N5,N4)+trcsa_RFXS(idsa_FeHPO4B,N6,N5,N4)
  trcsa_XFXS(idsa_FeH2PO4B,N6,N5,N4)=trcsa_XFXS(idsa_FeH2PO4B,N6,N5,N4)+trcsa_RFXS(idsa_FeH2PO4B,N6,N5,N4)
  trcsa_XFXS(idsa_CaPO4B,N6,N5,N4)=trcsa_XFXS(idsa_CaPO4B,N6,N5,N4)+trcsa_RFXS(idsa_CaPO4B,N6,N5,N4)
  trcsa_XFXS(idsa_CaHPO4B,N6,N5,N4)=trcsa_XFXS(idsa_CaHPO4B,N6,N5,N4)+trcsa_RFXS(idsa_CaHPO4B,N6,N5,N4)
  trcsa_XFXS(idsa_CaH2PO4B,N6,N5,N4)=trcsa_XFXS(idsa_CaH2PO4B,N6,N5,N4)+trcsa_RFXS(idsa_CaH2PO4B,N6,N5,N4)
  trcsa_XFXS(idsa_MgHPO4B,N6,N5,N4)=trcsa_XFXS(idsa_MgHPO4B,N6,N5,N4)+trcsa_RFXS(idsa_MgHPO4B,N6,N5,N4)
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
          trcsa_RFLS(idsa_Al,N,N6,N5,N4)=0.0_r8
          trcsa_RFLS(idsa_Fe,N,N6,N5,N4)=0.0_r8
          trcsa_RFLS(idsa_Hp,N,N6,N5,N4)=0.0_r8
          trcsa_RFLS(idsa_Ca,N,N6,N5,N4)=0.0_r8
          trcsa_RFLS(idsa_Mg,N,N6,N5,N4)=0.0_r8
          trcsa_RFLS(idsa_Na,N,N6,N5,N4)=0.0_r8
          trcsa_RFLS(idsa_K,N,N6,N5,N4)=0.0_r8
          trcsa_RFLS(idsa_OH,N,N6,N5,N4)=0.0_r8
          trcsa_RFLS(idsa_SO4,N,N6,N5,N4)=0.0_r8
          trcsa_RFLS(idsa_Cl,N,N6,N5,N4)=0.0_r8
          trcsa_RFLS(idsa_CO3,N,N6,N5,N4)=0.0_r8
          trcsa_RFLS(idsa_HCO3,N,N6,N5,N4)=0.0_r8
          trcsa_RFLS(idsa_AlOH,N,N6,N5,N4)=0.0_r8
          trcsa_RFLS(idsa_AlOH2,N,N6,N5,N4)=0.0_r8
          trcsa_RFLS(idsa_AlOH3,N,N6,N5,N4)=0.0_r8
          trcsa_RFLS(idsa_AlOH4,N,N6,N5,N4)=0.0_r8
          trcsa_RFLS(idsa_AlSO4,N,N6,N5,N4)=0.0_r8
          trcsa_RFLS(idsa_FeOH,N,N6,N5,N4)=0.0_r8
          trcsa_RFLS(idsa_FeOH2,N,N6,N5,N4)=0.0_r8
          trcsa_RFLS(idsa_FeOH3,N,N6,N5,N4)=0.0_r8
          trcsa_RFLS(idsa_FeOH4,N,N6,N5,N4)=0.0_r8
          trcsa_RFLS(idsa_FeSO4,N,N6,N5,N4)=0.0_r8
          trcsa_RFLS(idsa_CaOH2,N,N6,N5,N4)=0.0_r8
          trcsa_RFLS(idsa_CaCO3,N,N6,N5,N4)=0.0_r8
          trcsa_RFLS(idsa_CaHCO3,N,N6,N5,N4)=0.0_r8
          trcsa_RFLS(idsa_CaSO4,N,N6,N5,N4)=0.0_r8
          trcsa_RFLS(idsa_MgOH2,N,N6,N5,N4)=0.0_r8
          trcsa_RFLS(idsa_MgCO3,N,N6,N5,N4)=0.0_r8
          trcsa_RFLS(idsa_MgHCO3,N,N6,N5,N4)=0.0_r8
          trcsa_RFLS(idsa_MgSO4,N,N6,N5,N4)=0.0_r8
          trcsa_RFLS(idsa_NaCO3,N,N6,N5,N4)=0.0_r8
          trcsa_RFLS(idsa_NaSO4,N,N6,N5,N4)=0.0_r8
          trcsa_RFLS(idsa_KSO4,N,N6,N5,N4)=0.0_r8
          trcsa_RFLS(idsa_H0PO4,N,N6,N5,N4)=0.0_r8
          trcsa_RFLS(idsa_H3PO4,N,N6,N5,N4)=0.0_r8
          trcsa_RFLS(idsa_FeHPO4,N,N6,N5,N4)=0.0_r8
          trcsa_RFLS(idsa_FeH2PO4,N,N6,N5,N4)=0.0_r8
          trcsa_RFLS(idsa_CaPO4,N,N6,N5,N4)=0.0_r8
          trcsa_RFLS(idsa_CaHPO4,N,N6,N5,N4)=0.0_r8
          trcsa_RFLS(idsa_CaH2PO4,N,N6,N5,N4)=0.0_r8
          trcsa_RFLS(idsa_MgHPO4,N,N6,N5,N4)=0.0_r8
          trcsa_RFLS(idsa_H0PO4B,N,N6,N5,N4)=0.0_r8
          trcsa_RFLS(idsa_H3PO4B,N,N6,N5,N4)=0.0_r8
          trcsa_RFLS(idsa_FeHPO4B,N,N6,N5,N4)=0.0_r8
          trcsa_RFLS(idsa_FeH2PO4B,N,N6,N5,N4)=0.0_r8
          trcsa_RFLS(idsa_CaPO4B,N,N6,N5,N4)=0.0_r8
          trcsa_RFLS(idsa_CaHPO4B,N,N6,N5,N4)=0.0_r8
          trcsa_RFLS(idsa_CaH2PO4B,N,N6,N5,N4)=0.0_r8
          trcsa_RFLS(idsa_MgHPO4B,N,N6,N5,N4)=0.0_r8
          RALFHS(N,N6,N5,N4)=0.0_r8
          RFEFHS(N,N6,N5,N4)=0.0_r8
          RHYFHS(N,N6,N5,N4)=0.0_r8
          RCAFHS(N,N6,N5,N4)=0.0_r8
          RMGFHS(N,N6,N5,N4)=0.0_r8
          RNAFHS(N,N6,N5,N4)=0.0_r8
          RKAFHS(N,N6,N5,N4)=0.0_r8
          ROHFHS(N,N6,N5,N4)=0.0_r8
          RSOFHS(N,N6,N5,N4)=0.0_r8
          RCLFHS(N,N6,N5,N4)=0.0_r8
          RC3FHS(N,N6,N5,N4)=0.0_r8
          RHCFHS(N,N6,N5,N4)=0.0_r8
          RAL1HS(N,N6,N5,N4)=0.0_r8
          RAL2HS(N,N6,N5,N4)=0.0_r8
          RAL3HS(N,N6,N5,N4)=0.0_r8
          RAL4HS(N,N6,N5,N4)=0.0_r8
          RALSHS(N,N6,N5,N4)=0.0_r8
          RFE1HS(N,N6,N5,N4)=0.0_r8
          RFE2HS(N,N6,N5,N4)=0.0_r8
          RFE3HS(N,N6,N5,N4)=0.0_r8
          RFE4HS(N,N6,N5,N4)=0.0_r8
          RFESHS(N,N6,N5,N4)=0.0_r8
          RCAOHS(N,N6,N5,N4)=0.0_r8
          RCACHS(N,N6,N5,N4)=0.0_r8
          RCAHHS(N,N6,N5,N4)=0.0_r8
          RCASHS(N,N6,N5,N4)=0.0_r8
          RMGOHS(N,N6,N5,N4)=0.0_r8
          RMGCHS(N,N6,N5,N4)=0.0_r8
          RMGHHS(N,N6,N5,N4)=0.0_r8
          RMGSHS(N,N6,N5,N4)=0.0_r8
          RNACHS(N,N6,N5,N4)=0.0_r8
          RNASHS(N,N6,N5,N4)=0.0_r8
          RKASHS(N,N6,N5,N4)=0.0_r8
          RH0PHS(N,N6,N5,N4)=0.0_r8
          RH3PHS(N,N6,N5,N4)=0.0_r8
          RF1PHS(N,N6,N5,N4)=0.0_r8
          RF2PHS(N,N6,N5,N4)=0.0_r8
          RC0PHS(N,N6,N5,N4)=0.0_r8
          RC1PHS(N,N6,N5,N4)=0.0_r8
          RC2PHS(N,N6,N5,N4)=0.0_r8
          RM1PHS(N,N6,N5,N4)=0.0_r8
          RH0BHB(N,N6,N5,N4)=0.0_r8
          RH3BHB(N,N6,N5,N4)=0.0_r8
          RF1BHB(N,N6,N5,N4)=0.0_r8
          RF2BHB(N,N6,N5,N4)=0.0_r8
          RC0BHB(N,N6,N5,N4)=0.0_r8
          RC1BHB(N,N6,N5,N4)=0.0_r8
          RC2BHB(N,N6,N5,N4)=0.0_r8
          RM1BHB(N,N6,N5,N4)=0.0_r8
        ENDIF
      ELSE
        THETW1(N3,N2,N1)=0.0_r8
        THETW1(N6,N5,N4)=0.0_r8
        trcsa_RFLS(idsa_Al,N,N6,N5,N4)=0.0_r8
        trcsa_RFLS(idsa_Fe,N,N6,N5,N4)=0.0_r8
        trcsa_RFLS(idsa_Hp,N,N6,N5,N4)=0.0_r8
        trcsa_RFLS(idsa_Ca,N,N6,N5,N4)=0.0_r8
        trcsa_RFLS(idsa_Mg,N,N6,N5,N4)=0.0_r8
        trcsa_RFLS(idsa_Na,N,N6,N5,N4)=0.0_r8
        trcsa_RFLS(idsa_K,N,N6,N5,N4)=0.0_r8
        trcsa_RFLS(idsa_OH,N,N6,N5,N4)=0.0_r8
        trcsa_RFLS(idsa_SO4,N,N6,N5,N4)=0.0_r8
        trcsa_RFLS(idsa_Cl,N,N6,N5,N4)=0.0_r8
        trcsa_RFLS(idsa_CO3,N,N6,N5,N4)=0.0_r8
        trcsa_RFLS(idsa_HCO3,N,N6,N5,N4)=0.0_r8
        trcsa_RFLS(idsa_AlOH,N,N6,N5,N4)=0.0_r8
        trcsa_RFLS(idsa_AlOH2,N,N6,N5,N4)=0.0_r8
        trcsa_RFLS(idsa_AlOH3,N,N6,N5,N4)=0.0_r8
        trcsa_RFLS(idsa_AlOH4,N,N6,N5,N4)=0.0_r8
        trcsa_RFLS(idsa_AlSO4,N,N6,N5,N4)=0.0_r8
        trcsa_RFLS(idsa_FeOH,N,N6,N5,N4)=0.0_r8
        trcsa_RFLS(idsa_FeOH2,N,N6,N5,N4)=0.0_r8
        trcsa_RFLS(idsa_FeOH3,N,N6,N5,N4)=0.0_r8
        trcsa_RFLS(idsa_FeOH4,N,N6,N5,N4)=0.0_r8
        trcsa_RFLS(idsa_FeSO4,N,N6,N5,N4)=0.0_r8
        trcsa_RFLS(idsa_CaOH2,N,N6,N5,N4)=0.0_r8
        trcsa_RFLS(idsa_CaCO3,N,N6,N5,N4)=0.0_r8
        trcsa_RFLS(idsa_CaHCO3,N,N6,N5,N4)=0.0_r8
        trcsa_RFLS(idsa_CaSO4,N,N6,N5,N4)=0.0_r8
        trcsa_RFLS(idsa_MgOH2,N,N6,N5,N4)=0.0_r8
        trcsa_RFLS(idsa_MgCO3,N,N6,N5,N4)=0.0_r8
        trcsa_RFLS(idsa_MgHCO3,N,N6,N5,N4)=0.0_r8
        trcsa_RFLS(idsa_MgSO4,N,N6,N5,N4)=0.0_r8
        trcsa_RFLS(idsa_NaCO3,N,N6,N5,N4)=0.0_r8
        trcsa_RFLS(idsa_NaSO4,N,N6,N5,N4)=0.0_r8
        trcsa_RFLS(idsa_KSO4,N,N6,N5,N4)=0.0_r8
        trcsa_RFLS(idsa_H0PO4,N,N6,N5,N4)=0.0_r8
        trcsa_RFLS(idsa_H3PO4,N,N6,N5,N4)=0.0_r8
        trcsa_RFLS(idsa_FeHPO4,N,N6,N5,N4)=0.0_r8
        trcsa_RFLS(idsa_FeH2PO4,N,N6,N5,N4)=0.0_r8
        trcsa_RFLS(idsa_CaPO4,N,N6,N5,N4)=0.0_r8
        trcsa_RFLS(idsa_CaHPO4,N,N6,N5,N4)=0.0_r8
        trcsa_RFLS(idsa_CaH2PO4,N,N6,N5,N4)=0.0_r8
        trcsa_RFLS(idsa_MgHPO4,N,N6,N5,N4)=0.0_r8
        trcsa_RFLS(idsa_H0PO4B,N,N6,N5,N4)=0.0_r8
        trcsa_RFLS(idsa_H3PO4B,N,N6,N5,N4)=0.0_r8
        trcsa_RFLS(idsa_FeHPO4B,N,N6,N5,N4)=0.0_r8
        trcsa_RFLS(idsa_FeH2PO4B,N,N6,N5,N4)=0.0_r8
        trcsa_RFLS(idsa_CaPO4B,N,N6,N5,N4)=0.0_r8
        trcsa_RFLS(idsa_CaHPO4B,N,N6,N5,N4)=0.0_r8
        trcsa_RFLS(idsa_CaH2PO4B,N,N6,N5,N4)=0.0_r8
        trcsa_RFLS(idsa_MgHPO4B,N,N6,N5,N4)=0.0_r8
        RALFHS(N,N6,N5,N4)=0.0_r8
        RFEFHS(N,N6,N5,N4)=0.0_r8
        RHYFHS(N,N6,N5,N4)=0.0_r8
        RCAFHS(N,N6,N5,N4)=0.0_r8
        RMGFHS(N,N6,N5,N4)=0.0_r8
        RNAFHS(N,N6,N5,N4)=0.0_r8
        RKAFHS(N,N6,N5,N4)=0.0_r8
        ROHFHS(N,N6,N5,N4)=0.0_r8
        RSOFHS(N,N6,N5,N4)=0.0_r8
        RCLFHS(N,N6,N5,N4)=0.0_r8
        RC3FHS(N,N6,N5,N4)=0.0_r8
        RHCFHS(N,N6,N5,N4)=0.0_r8
        RAL1HS(N,N6,N5,N4)=0.0_r8
        RAL2HS(N,N6,N5,N4)=0.0_r8
        RAL3HS(N,N6,N5,N4)=0.0_r8
        RAL4HS(N,N6,N5,N4)=0.0_r8
        RALSHS(N,N6,N5,N4)=0.0_r8
        RFE1HS(N,N6,N5,N4)=0.0_r8
        RFE2HS(N,N6,N5,N4)=0.0_r8
        RFE3HS(N,N6,N5,N4)=0.0_r8
        RFE4HS(N,N6,N5,N4)=0.0_r8
        RFESHS(N,N6,N5,N4)=0.0_r8
        RCAOHS(N,N6,N5,N4)=0.0_r8
        RCACHS(N,N6,N5,N4)=0.0_r8
        RCAHHS(N,N6,N5,N4)=0.0_r8
        RCASHS(N,N6,N5,N4)=0.0_r8
        RMGOHS(N,N6,N5,N4)=0.0_r8
        RMGCHS(N,N6,N5,N4)=0.0_r8
        RMGHHS(N,N6,N5,N4)=0.0_r8
        RMGSHS(N,N6,N5,N4)=0.0_r8
        RNACHS(N,N6,N5,N4)=0.0_r8
        RNASHS(N,N6,N5,N4)=0.0_r8
        RKASHS(N,N6,N5,N4)=0.0_r8
        RH0PHS(N,N6,N5,N4)=0.0_r8
        RH3PHS(N,N6,N5,N4)=0.0_r8
        RF1PHS(N,N6,N5,N4)=0.0_r8
        RF2PHS(N,N6,N5,N4)=0.0_r8
        RC0PHS(N,N6,N5,N4)=0.0_r8
        RC1PHS(N,N6,N5,N4)=0.0_r8
        RC2PHS(N,N6,N5,N4)=0.0_r8
        RM1PHS(N,N6,N5,N4)=0.0_r8
        RH0BHB(N,N6,N5,N4)=0.0_r8
        RH3BHB(N,N6,N5,N4)=0.0_r8
        RF1BHB(N,N6,N5,N4)=0.0_r8
        RF2BHB(N,N6,N5,N4)=0.0_r8
        RC0BHB(N,N6,N5,N4)=0.0_r8
        RC1BHB(N,N6,N5,N4)=0.0_r8
        RC2BHB(N,N6,N5,N4)=0.0_r8
        RM1BHB(N,N6,N5,N4)=0.0_r8
      ENDIF
    enddo
  enddo
  end subroutine UpdateSoluteInSubsurfNeighbors
!------------------------------------------------------------------------------------------

  subroutine MacMicPoreFluxAdvPlusDifus(NY,NX)
!
!     Description:
!
  implicit none
  integer, intent(in) :: NY,NX
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
  trcsa_RFXS(idsa_Al,NU(NY,NX),NY,NX)=RFLAL+DFVAL
  trcsa_RFXS(idsa_Fe,NU(NY,NX),NY,NX)=RFLFE+DFVFE
  trcsa_RFXS(idsa_Hp,NU(NY,NX),NY,NX)=RFLHY+DFVHY
  trcsa_RFXS(idsa_Ca,NU(NY,NX),NY,NX)=RFLCA+DFVCA
  trcsa_RFXS(idsa_Mg,NU(NY,NX),NY,NX)=RFLMG+DFVMG
  trcsa_RFXS(idsa_Na,NU(NY,NX),NY,NX)=RFLNA+DFVNA
  trcsa_RFXS(idsa_K,NU(NY,NX),NY,NX)=RFLKA+DFVKA
  trcsa_RFXS(idsa_OH,NU(NY,NX),NY,NX)=RFLOH+DFVOH
  trcsa_RFXS(idsa_SO4,NU(NY,NX),NY,NX)=RFLSO+DFVSO
  trcsa_RFXS(idsa_Cl,NU(NY,NX),NY,NX)=RFLCL+DFVCL
  trcsa_RFXS(idsa_CO3,NU(NY,NX),NY,NX)=RFLC3+DFVC3
  trcsa_RFXS(idsa_HCO3,NU(NY,NX),NY,NX)=RFLHC+DFVHC
  trcsa_RFXS(idsa_AlOH,NU(NY,NX),NY,NX)=RFLAL1+DFVAL1
  trcsa_RFXS(idsa_AlOH2,NU(NY,NX),NY,NX)=RFLAL2+DFVAL2
  trcsa_RFXS(idsa_AlOH3,NU(NY,NX),NY,NX)=RFLAL3+DFVAL3
  trcsa_RFXS(idsa_AlOH4,NU(NY,NX),NY,NX)=RFLAL4+DFVAL4
  trcsa_RFXS(idsa_AlSO4,NU(NY,NX),NY,NX)=RFLALS+DFVALS
  trcsa_RFXS(idsa_FeOH,NU(NY,NX),NY,NX)=RFLFE1+DFVFE1
  trcsa_RFXS(idsa_FeOH2,NU(NY,NX),NY,NX)=RFLFE2+DFVFE2
  trcsa_RFXS(idsa_FeOH3,NU(NY,NX),NY,NX)=RFLFE3+DFVFE3
  trcsa_RFXS(idsa_FeOH4,NU(NY,NX),NY,NX)=RFLFE4+DFVFE4
  trcsa_RFXS(idsa_FeSO4,NU(NY,NX),NY,NX)=RFLFES+DFVFES
  trcsa_RFXS(idsa_CaOH2,NU(NY,NX),NY,NX)=RFLCAO+DFVCAO
  trcsa_RFXS(idsa_CaCO3,NU(NY,NX),NY,NX)=RFLCAC+DFVCAC
  trcsa_RFXS(idsa_CaHCO3,NU(NY,NX),NY,NX)=RFLCAH+DFVCAH
  trcsa_RFXS(idsa_CaSO4,NU(NY,NX),NY,NX)=RFLCAS+DFVCAS
  trcsa_RFXS(idsa_MgOH2,NU(NY,NX),NY,NX)=RFLMGO+DFVMGO
  trcsa_RFXS(idsa_MgCO3,NU(NY,NX),NY,NX)=RFLMGC+DFVMGC
  trcsa_RFXS(idsa_MgHCO3,NU(NY,NX),NY,NX)=RFLMGH+DFVMGH
  trcsa_RFXS(idsa_MgSO4,NU(NY,NX),NY,NX)=RFLMGS+DFVMGS
  trcsa_RFXS(idsa_NaCO3,NU(NY,NX),NY,NX)=RFLNAC+DFVNAC
  trcsa_RFXS(idsa_NaSO4,NU(NY,NX),NY,NX)=RFLNAS+DFVNAS
  trcsa_RFXS(idsa_KSO4,NU(NY,NX),NY,NX)=RFLKAS+DFVKAS
  trcsa_RFXS(idsa_H0PO4,NU(NY,NX),NY,NX)=RFLH0P+DFVH0P
  trcsa_RFXS(idsa_H3PO4,NU(NY,NX),NY,NX)=RFLH3P+DFVH3P
      trcsa_RFXS(idsa_FeHPO4,NU(NY,NX),NY,NX)=RFLF1P+DFVF1P
      trcsa_RFXS(idsa_FeH2PO4,NU(NY,NX),NY,NX)=RFLF2P+DFVF2P
      trcsa_RFXS(idsa_CaPO4,NU(NY,NX),NY,NX)=RFLC0P+DFVC0P
      trcsa_RFXS(idsa_CaHPO4,NU(NY,NX),NY,NX)=RFLC1P+DFVC1P
      trcsa_RFXS(idsa_CaH2PO4,NU(NY,NX),NY,NX)=RFLC2P+DFVC2P
      trcsa_RFXS(idsa_MgHPO4,NU(NY,NX),NY,NX)=RFLM1P+DFVM1P
      trcsa_RFXS(idsa_H0PO4B,NU(NY,NX),NY,NX)=RFLH0B+DFVH0B
      trcsa_RFXS(idsa_H3PO4B,NU(NY,NX),NY,NX)=RFLH3B+DFVH3B
      trcsa_RFXS(idsa_FeHPO4B,NU(NY,NX),NY,NX)=RFLF1B+DFVF1B
      trcsa_RFXS(idsa_FeH2PO4B,NU(NY,NX),NY,NX)=RFLF2B+DFVF2B
      trcsa_RFXS(idsa_CaPO4B,NU(NY,NX),NY,NX)=RFLC0B+DFVC0B
      trcsa_RFXS(idsa_CaHPO4B,NU(NY,NX),NY,NX)=RFLC1B+DFVC1B
      trcsa_RFXS(idsa_CaH2PO4B,NU(NY,NX),NY,NX)=RFLC2B+DFVC2B
      trcsa_RFXS(idsa_MgHPO4B,NU(NY,NX),NY,NX)=RFLM1B+DFVM1B
      end subroutine MacMicPoreFluxAdvPlusDifus
end module IngridTranspMod
