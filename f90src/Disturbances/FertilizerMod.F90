module FertilizerMod
  use data_kind_mod, only: r8 => DAT_KIND_R8
  use EcoSIMCtrlMod, only: lverb
  use DebugToolMod,  only: PrintInfo
  use SurfLitterDataType
  use EcoSimSumDataType
  use FertilizerDataType
  use ClimForcDataType
  use SoilBGCDataType
  use FlagDataType
  use GridDataType
  use EcosimBGCFluxType
  use SOMDataType
  use EcoSiMParDataMod
  use EcoSIMCtrlDataType
  use SoilPropertyDataType
  use MicrobialDataType
  use AqueChemDatatype
  use MiniMathMod
  use EcosimConst
  use fileUtil
implicit none
  private

  CHARACTER(LEN=*), PARAMETER :: MOD_FILENAME=&
  __FILE__

  public :: ApplyFertilizerAtNoon

  contains
!------------------------------------------------------------------------------------------

  subroutine ApplyFertilizerAtNoon(I,J,NHW,NHE,NVN,NVS)
!
  implicit none
  integer, intent(in) :: I,J,NHW,NHE,NVN,NVS

  integer :: NX,NY
  real(r8) :: OFC(2),OFN(2),OFP(2)
  integer :: LFDPTH = 0
!     begin_execution

  D8990: DO NX=NHW,NHE
    D8995: DO NY=NVN,NVS
      IF(J.EQ.INT(SolarNoonHour_col(NY,NX)))THEN

        call ApplyMineralFertilizer(I,J,NY,NX,LFDPTH,OFC,OFN,OFP)
!
!     SOIL LAYER NUMBER IN WHICH PLANT OR ANIMAL RESIDUES ARE APPLIED
!
       if(lverb)write(iulog,*) 'ApplyManure'
        call ApplyManure(I,J,NY,NX,LFDPTH,OFC,OFN,OFP)
!
!     FERTILIZER UREA, NITRIFICATION INHIBITORS
        if(lverb)write(iulog,*)'ApplyUreaNitrifierInhibitor'
        call ApplyUreaNitrifierInhibitor(I,J,NY,NX,LFDPTH)

      ENDIF

    ENDDO D8995
  ENDDO D8990
  end subroutine ApplyFertilizerAtNoon
!------------------------------------------------------------------------------------------

  subroutine ApplyUreaNitrifierInhibitor(I,J,NY,NX,LFDPTH)

  implicit none
  integer, intent(in) :: I,J,NY,NX
  integer, intent(in) :: LFDPTH
  integer :: L
!     begin_execution
!
!     IYTYP=fertilizer release type from fertilizer input file
!     FERT=fertilizer type from fertilizer input file
!     IUTYP=urea hydrolysis inhibitor type (1=no,2=yes)
!     ZNHU0,ZNHUI=initial,current urea hydrolysis inhibition activity
!     ZNFN0,ZNFNI=initial,current nitrification inhibition activity
! urea application
  IF(FERT(ifert_urea,I,NY,NX).GT.0._r8 .OR. FERT(ifert_urea_band,I,NY,NX).GT.0._r8)THEN  
    IF(IYTYP(0,I,NY,NX).EQ.0)THEN
      IUTYP(NY,NX)=0
    ELSEIF(IYTYP(0,I,NY,NX).EQ.1 .OR. IYTYP(0,I,NY,NX).EQ.3)THEN
      IUTYP(NY,NX)=1
    ELSE
      !urea hydrolysis is on
      IUTYP(NY,NX)=2
    ENDIF

    D9964: DO L=0,NL(NY,NX)
      IF(L.EQ.LFDPTH)THEN
        ZNHU0_vr(L,NY,NX)=1.0_r8
        ZNHUI_vr(L,NY,NX)=1.0_r8
      ELSE
        ZNHU0_vr(L,NY,NX)=0._r8
        ZNHUI_vr(L,NY,NX)=0._r8
      ENDIF
    ENDDO D9964
  ENDIF

  IF(IYTYP(0,I,NY,NX).EQ.3 .OR. IYTYP(0,I,NY,NX).EQ.4)THEN
    D9965: DO L=0,NL(NY,NX)
      IF(L.EQ.LFDPTH)THEN
        ZNFN0_vr(L,NY,NX)=1.0_r8
        ZNFNI_vr(L,NY,NX)=1.0_r8
      ELSE
        ZNFN0_vr(L,NY,NX)=0._r8
        ZNFNI_vr(L,NY,NX)=0._r8
      ENDIF
    ENDDO D9965
  ENDIF
  end subroutine ApplyUreaNitrifierInhibitor
!------------------------------------------------------------------------------------------

  subroutine ApplyManure(I,J,NY,NX,LFDPTH,OFC,OFN,OFP)
  implicit none
  integer, intent(in) :: I,J,NY,NX
  integer, intent(inout) :: LFDPTH
  real(r8), intent(in) :: OFC(2),OFN(2),OFP(2)
  real(r8) :: CNOF(4),CPOF(4)
  real(r8) :: CORGCX
  real(r8) :: CNOFT
  real(r8) :: CPOFT
  real(r8) :: FDPTHM
  real(r8) :: FRNT,FRPT
  REAL(R8) :: OSCI,OSNI,OSPI
  REAL(R8) :: OSCX,OSNX,OSPX
  REAL(R8) :: OMC1,OMN1,OMP1
  real(r8) :: OQC1,OQN1,OQP1
  real(r8) :: OSC1,OSN1,OSP1
  REAL(R8) :: RNT,RPT
  integer  :: L,K,M,N,NN,NGL,MID
  real(r8) :: tglds
  real(r8) :: OMC1g,OMN1g,OMP1g
!     begin_execution
  associate(                           &
    k_fine_litr => micpar%k_fine_litr, &
    k_manure    => micpar%k_manure   , &
    ilignin     => micpar%ilignin    , &
    icellulos   => micpar%icellulos  , &
    icarbhyro   => micpar%icarbhyro  , &
    iprotein    => micpar%iprotein     &
  )
!     LFDPTH=layer number
!
  IF(OFC(1)+OFC(2).GT.0._r8)THEN
    DO  L=0,JZ
      FDPTHM=FDPTH(I,NY,NX)+CumDepz2LayBottom_vr(NU(NY,NX)-1,NY,NX)
      IF(FDPTHM.LE.0._r8)THEN
        LFDPTH=0
        exit
      ELSEIF(CumDepz2LayBottom_vr(L,NY,NX).GE.FDPTHM)THEN
        LFDPTH=L
        exit
      ENDIF
    ENDDO
!
!     ALLOCATION OF PLANT RESIDUE APPLICATION TO
!     RESIDUE PROTEIN, CH2O, CELLULOSE, LIGNIN
!
!     CFOSC=fraction of litter allocated to protein(1)
!     soluble CH2O(2), cellulose(3) and lignin(4)
!     ITYPE=litter type entered in fertilizer input file
!
!     MAIZE
!
    IF(IYTYP(1,I,NY,NX).EQ.1)THEN
      CFOSC(iprotein,k_fine_litr,LFDPTH,NY,NX)=0.080_r8
      CFOSC(icarbhyro,k_fine_litr,LFDPTH,NY,NX)=0.245_r8
      CFOSC(icellulos,k_fine_litr,LFDPTH,NY,NX)=0.613_r8
      CFOSC(ilignin,k_fine_litr,LFDPTH,NY,NX)=0.062_r8
!
!     WHEAT
!
    ELSEIF(IYTYP(1,I,NY,NX).EQ.2)THEN
      CFOSC(iprotein,k_fine_litr,LFDPTH,NY,NX)=0.125_r8
      CFOSC(icarbhyro,k_fine_litr,LFDPTH,NY,NX)=0.171_r8
      CFOSC(icellulos,k_fine_litr,LFDPTH,NY,NX)=0.560_r8
      CFOSC(ilignin,k_fine_litr,LFDPTH,NY,NX)=0.144_r8
!
!     SOYBEAN
!
    ELSEIF(IYTYP(1,I,NY,NX).EQ.3)THEN
      CFOSC(iprotein,k_fine_litr,LFDPTH,NY,NX)=0.138_r8
      CFOSC(icarbhyro,k_fine_litr,LFDPTH,NY,NX)=0.426_r8
      CFOSC(icellulos,k_fine_litr,LFDPTH,NY,NX)=0.316_r8
      CFOSC(ilignin,k_fine_litr,LFDPTH,NY,NX)=0.120_r8
!
!     OLD STRAW
!
    ELSEIF(IYTYP(1,I,NY,NX).EQ.4)THEN
      CFOSC(iprotein,k_fine_litr,LFDPTH,NY,NX)=0.075_r8
      CFOSC(icarbhyro,k_fine_litr,LFDPTH,NY,NX)=0.125_r8
      CFOSC(icellulos,k_fine_litr,LFDPTH,NY,NX)=0.550_r8
      CFOSC(ilignin,k_fine_litr,LFDPTH,NY,NX)=0.250_r8
!
!     STRAW
!
    ELSEIF(IYTYP(1,I,NY,NX).EQ.5)THEN
      CFOSC(iprotein,k_fine_litr,LFDPTH,NY,NX)=0.036_r8
      CFOSC(icarbhyro,k_fine_litr,LFDPTH,NY,NX)=0.044_r8
      CFOSC(icellulos,k_fine_litr,LFDPTH,NY,NX)=0.767_r8
      CFOSC(ilignin,k_fine_litr,LFDPTH,NY,NX)=0.153_r8
!
!     COMPOST
!
    ELSEIF(IYTYP(1,I,NY,NX).EQ.6)THEN
      CFOSC(iprotein,k_fine_litr,LFDPTH,NY,NX)=0.143_r8
      CFOSC(icarbhyro,k_fine_litr,LFDPTH,NY,NX)=0.015_r8
      CFOSC(icellulos,k_fine_litr,LFDPTH,NY,NX)=0.640_r8
      CFOSC(ilignin,k_fine_litr,LFDPTH,NY,NX)=0.202_r8
!
!     GREEN MANURE
!
    ELSEIF(IYTYP(1,I,NY,NX).EQ.7)THEN
      CFOSC(iprotein,k_fine_litr,LFDPTH,NY,NX)=0.202_r8
      CFOSC(icarbhyro,k_fine_litr,LFDPTH,NY,NX)=0.013_r8
      CFOSC(icellulos,k_fine_litr,LFDPTH,NY,NX)=0.560_r8
      CFOSC(ilignin,k_fine_litr,LFDPTH,NY,NX)=0.225_r8
!
!     SIMPLE SUBSTRATE
!
    ELSEIF(IYTYP(1,I,NY,NX).EQ.10)THEN
      CFOSC(iprotein,k_fine_litr,LFDPTH,NY,NX)=0.000_r8
      CFOSC(icarbhyro,k_fine_litr,LFDPTH,NY,NX)=1.000_r8
      CFOSC(icellulos,k_fine_litr,LFDPTH,NY,NX)=0.000_r8
      CFOSC(ilignin,k_fine_litr,LFDPTH,NY,NX)=0.000_r8
    ELSE
      CFOSC(iprotein,k_fine_litr,LFDPTH,NY,NX)=0.075_r8
      CFOSC(icarbhyro,k_fine_litr,LFDPTH,NY,NX)=0.125_r8
      CFOSC(icellulos,k_fine_litr,LFDPTH,NY,NX)=0.550_r8
      CFOSC(ilignin,k_fine_litr,LFDPTH,NY,NX)=0.250_r8
    ENDIF
!
!     ALLOCATION OF ANIMAL MANURE APPLICATION TO
!     RESIDUE PROTEIN, CH2O, CELLULOSE, LIGNIN
!
!     RUMINANT
!
    IF(IYTYP(2,I,NY,NX).EQ.1)THEN
      CFOSC(iprotein,k_manure,LFDPTH,NY,NX)=0.036_r8
      CFOSC(icarbhyro,k_manure,LFDPTH,NY,NX)=0.044_r8
      CFOSC(icellulos,k_manure,LFDPTH,NY,NX)=0.630_r8
      CFOSC(ilignin,k_manure,LFDPTH,NY,NX)=0.290_r8
!
!     NON-RUMINANT
!
    ELSEIF(IYTYP(2,I,NY,NX).EQ.2)THEN
      CFOSC(iprotein,k_manure,LFDPTH,NY,NX)=0.138_r8
      CFOSC(icarbhyro,k_manure,LFDPTH,NY,NX)=0.401_r8
      CFOSC(icellulos,k_manure,LFDPTH,NY,NX)=0.316_r8
      CFOSC(ilignin,k_manure,LFDPTH,NY,NX)=0.145_r8
!
!     GRAZING
!
    ELSEIF(IYTYP(2,I,NY,NX).EQ.3)THEN
      CFOSC(iprotein,k_manure,LFDPTH,NY,NX)=0.036_r8
      CFOSC(icarbhyro,k_manure,LFDPTH,NY,NX)=0.044_r8
      CFOSC(icellulos,k_manure,LFDPTH,NY,NX)=0.630_r8
      CFOSC(ilignin,k_manure,LFDPTH,NY,NX)=0.290_r8
!
!     OTHER
!
    ELSE
      CFOSC(iprotein,k_manure,LFDPTH,NY,NX)=0.138_r8
      CFOSC(icarbhyro,k_manure,LFDPTH,NY,NX)=0.401_r8
      CFOSC(icellulos,k_manure,LFDPTH,NY,NX)=0.316_r8
      CFOSC(ilignin,k_manure,LFDPTH,NY,NX)=0.145_r8
    ENDIF
!
!     DISTRIBUTE RESIDUE APPLICATION AMONG COMPONENTS OF RESIDUE complex
!
!     OFC,OFN,OFP=litter C,N,P application from fertilizer file
!

    D2965: DO K=1,2
      OSCI=OFC(K)*AREA(3,LFDPTH,NY,NX)
      OSNI=OFN(K)*AREA(3,LFDPTH,NY,NX)
      OSPI=OFP(K)*AREA(3,LFDPTH,NY,NX)
      IF(VLSoilMicPMass_vr(LFDPTH,NY,NX).GT.ZEROS(NY,NX))THEN
        CORGCX=OSCI/VLSoilMicPMass_vr(LFDPTH,NY,NX)
      ELSE
        CORGCX=orgcden
      ENDIF
      OSCX=0._r8
      OSNX=0._r8
      OSPX=0._r8
!
!     BIOMASSES OF MICROBIAL POPULATIONS IN RESIDUE
!
!     OMC,OMN,OMP=microbial biomass in litter application
!     OMCI=microbial biomass content in litter
!     OMCF,OMCA=hetero,autotrophic biomass composition in litter
!
      D2960: DO N=1,NumMicbFunGrupsPerCmplx
        tglds=JGnfo(N)-JGnfo(N)+1
        D2961: DO M=1,nlbiomcp
          OMC1=AZMAX1(AMIN1(OSCI*micpar%OMCI(M,K)*micpar%OMCF(N),OSCI-OSCX))
          OMN1=AZMAX1(AMIN1(OMC1*micpar%rNCOMCa(M,N,K),OSNI-OSNX))
          OMP1=AZMAX1(AMIN1(OMC1*micpar%rPCOMCa(M,N,K),OSPI-OSPX))
          DO NGL=JGnio(N),JGnfo(N)
            MID=micpar%get_micb_id(M,NGL)
            OMC1g=OMC1/tglds
            OMN1g=OMN1/tglds
            OMP1g=OMP1/tglds
            mBiomeHeter_vr(ielmc,MID,K,LFDPTH,NY,NX)=mBiomeHeter_vr(ielmc,MID,K,LFDPTH,NY,NX)+OMC1g
            mBiomeHeter_vr(ielmn,MID,K,LFDPTH,NY,NX)=mBiomeHeter_vr(ielmn,MID,K,LFDPTH,NY,NX)+OMN1g
            mBiomeHeter_vr(ielmp,MID,K,LFDPTH,NY,NX)=mBiomeHeter_vr(ielmp,MID,K,LFDPTH,NY,NX)+OMP1g
          ENDDO
          OSCX=OSCX+OMC1
          OSNX=OSNX+OMN1
          OSPX=OSPX+OMP1
          D2962: DO NN=1,NumMicbFunGrupsPerCmplx
            tglds=JGnfA(N)-JGniA(N)+1
            DO NGL=JGniA(NN),JGnfA(NN)
              MID=micpar%get_micb_id(M,NGL)
              OMC1g=OMC1/tglds
              OMN1g=OMN1/tglds
              OMP1g=OMP1/tglds
              mBiomeAutor_vr(ielmc,MID,LFDPTH,NY,NX)=mBiomeAutor_vr(ielmc,MID,LFDPTH,NY,NX)+OMC1g*micpar%OMCA(NN)
              mBiomeAutor_vr(ielmn,MID,LFDPTH,NY,NX)=mBiomeAutor_vr(ielmn,MID,LFDPTH,NY,NX)+OMN1g*micpar%OMCA(NN)
              mBiomeAutor_vr(ielmp,MID,LFDPTH,NY,NX)=mBiomeAutor_vr(ielmp,MID,LFDPTH,NY,NX)+OMP1g*micpar%OMCA(NN)
            ENDDO
            OSCX=OSCX+OMC1*micpar%OMCA(NN)
            OSNX=OSNX+OMN1*micpar%OMCA(NN)
            OSPX=OSPX+OMP1*micpar%OMCA(NN)
          ENDDO D2962
        ENDDO D2961
      ENDDO D2960
!
!     DOC, DON AND DOP IN RESIDUE
!
!     OQC,OQN,OQP=DOC,DON,DOP in litter
!
      OQC1=AMIN1(0.1_r8*OSCX,OSCI-OSCX)
      OQN1=AMIN1(0.1_r8*OSNX,OSNI-OSNX)
      OQP1=AMIN1(0.1_r8*OSPX,OSPI-OSPX)

      DOM_vr(idom_doc,K,LFDPTH,NY,NX)=DOM_vr(idom_doc,K,LFDPTH,NY,NX)+OQC1
      DOM_vr(idom_don,K,LFDPTH,NY,NX)=DOM_vr(idom_don,K,LFDPTH,NY,NX)+OQN1
      DOM_vr(idom_dop,K,LFDPTH,NY,NX)=DOM_vr(idom_dop,K,LFDPTH,NY,NX)+OQP1
!
!     REMAINDER DISTRIBUTED TO RESIDUE FRACTIONS
!
!     OSC,OSN,OSP,OSA=SOC,SON,SOP,colonized SOC in litter
!     VOLT=litter volume
!     AmendCFlx_CumYr_col,FertNFlx_CumYr_col,FerPFlx_CumYr_col=accumulated litter C,N,P application
!     Eco_NBP_CumYr_col=accumulated net biome productivity
!
      OSCX=OSCX+OQC1
      OSNX=OSNX+OQN1
      OSPX=OSPX+OQP1
      CNOFT=0._r8
      CPOFT=0._r8
      IF(OSCI-OSCX.GT.ZEROS(NY,NX))THEN
        RNT=0._r8
        RPT=0._r8
        D965: DO M=1,jsken
          RNT=RNT+(OSCI-OSCX)*CFOSC(M,K,LFDPTH,NY,NX)*micpar%CNOFC(M,K)
          RPT=RPT+(OSCI-OSCX)*CFOSC(M,K,LFDPTH,NY,NX)*micpar%CPOFC(M,K)
        ENDDO D965
        FRNT=(OSNI-OSNX)/RNT
        FRPT=(OSPI-OSPX)/RPT
        D970: DO M=1,jsken
          CNOF(M)=micpar%CNOFC(M,K)*FRNT
          CPOF(M)=micpar%CPOFC(M,K)*FRPT
          CNOFT=CNOFT+CFOSC(M,K,LFDPTH,NY,NX)*CNOF(M)
          CPOFT=CPOFT+CFOSC(M,K,LFDPTH,NY,NX)*CPOF(M)
        ENDDO D970
      ELSE
        D975: DO M=1,jsken
          CNOF(M)=0._r8
          CPOF(M)=0._r8
        ENDDO D975
      ENDIF
      D2970: DO M=1,jsken
        OSC1=CFOSC(M,K,LFDPTH,NY,NX)*(OSCI-OSCX)
        IF(CNOFT.GT.ZERO)THEN
          OSN1=CFOSC(M,K,LFDPTH,NY,NX)*CNOF(M)/CNOFT*(OSNI-OSNX)
        ELSE
          OSN1=0._r8
        ENDIF
        IF(CPOFT.GT.ZERO)THEN
          OSP1=CFOSC(M,K,LFDPTH,NY,NX)*CPOF(M)/CPOFT*(OSPI-OSPX)
        ELSE
          OSP1=0._r8
        ENDIF
        SolidOM_vr(ielmc,M,K,LFDPTH,NY,NX)=SolidOM_vr(ielmc,M,K,LFDPTH,NY,NX)+OSC1
        SolidOM_vr(ielmn,M,K,LFDPTH,NY,NX)=SolidOM_vr(ielmn,M,K,LFDPTH,NY,NX)+OSN1
        SolidOM_vr(ielmp,M,K,LFDPTH,NY,NX)=SolidOM_vr(ielmp,M,K,LFDPTH,NY,NX)+OSP1
        SolidOMAct_vr(M,K,LFDPTH,NY,NX)=SolidOMAct_vr(M,K,LFDPTH,NY,NX)+OSC1*micpar%OMCI(1,K)
        IF(LFDPTH.EQ.0)THEN
          VGeomLayer_vr(LFDPTH,NY,NX)=VGeomLayer_vr(LFDPTH,NY,NX)+OSC1*ppmc/BulkDensLitR(micpar%k_fine_litr)
        ENDIF
      ENDDO D2970
      tAmendOrgC_lnd=tAmendOrgC_lnd+OSCI
      TORGN=TORGN+OSNI
      TORGP=TORGP+OSPI
      AmendCFlx_CumYr_col(NY,NX)=AmendCFlx_CumYr_col(NY,NX)+OSCI
      FertNFlx_CumYr_col(NY,NX)=FertNFlx_CumYr_col(NY,NX)+OSNI
      FerPFlx_CumYr_col(NY,NX)=FerPFlx_CumYr_col(NY,NX)+OSPI
      IF(IYTYP(2,I,NY,NX).LT.3)THEN
        Eco_NBP_CumYr_col(NY,NX)=Eco_NBP_CumYr_col(NY,NX)+OSCI
      ENDIF
    ENDDO D2965
  ENDIF
  end associate
  end subroutine ApplyManure
!------------------------------------------------------------------------------------------

  subroutine ApplyMineralFertilizer(I,J,NY,NX,LFDPTH,OFC,OFN,OFP)
  implicit none
  integer, intent(in) :: I,J,NY,NX
  real(r8), intent(out) :: OFC(2),OFN(2),OFP(2)
  integer, intent(out) :: LFDPTH

  character(len=*), parameter :: subname='ApplyMineralFertilizer'
  real(r8) :: BAREF
  real(r8) :: CVRDF
  real(r8) :: CAC
  real(r8) :: CAS
  real(r8) :: CACX
  real(r8) :: CASX
  real(r8) :: FDPTHF
  real(r8) :: H0PO4T,H1PO4T,H2PO4T,H3PO4T
  real(r8) :: PMA,PMB,PHA
  real(r8) :: PALPOT,PFEPOT
  real(r8) :: PCAPDT,PCAPHT,PCAPMT
  real(r8) :: PMAX,PMBX,PHAX
  REAL(R8) :: XN4T
  real(r8) :: XOH0T,XOH1T,XOH2T,XH1PT,XH2PT
  real(r8) :: Z4A,Z3A,ZUA,ZOA,Z4B,Z3B
  REAL(R8) :: ZUB,ZOB
  real(r8) :: ZNH4T,ZNH3T,ZNO3T,ZNO2T
  real(r8) :: ZFE1PT,ZFE2PT
  real(r8) :: ZCA0PT,ZCA1PT,ZCA2PT
  real(r8) :: ZMG1PT,Z4AX,Z3AX,ZUAX,ZOAX
  real(r8) :: Z4BX,Z3BX,ZUBX,ZOBX
  integer :: L

!     begin_execution
!
!     NH4,NH3,UREA,NO3 FERTILIZER APPLICATION
!
!     *A,*B=broadcast,banded
!     Z4,Z3,ZU,ZO=NH4,NH3,urea,NO3
!
  call PrintInfo('beg '//subname)
  Z4A=FERT(1,I,NY,NX)
  Z3A=FERT(2,I,NY,NX)
  ZUA=FERT(3,I,NY,NX)
  ZOA=FERT(4,I,NY,NX)
  Z4B=FERT(5,I,NY,NX)
  Z3B=FERT(6,I,NY,NX)
  ZUB=FERT(7,I,NY,NX)
  ZOB=FERT(8,I,NY,NX)
!
!     MONOCALCIUM PHOSPHATE OR HYDROXYAPATITE
!
!     PM*,PH*=Ca(H2PO4)2,apatite
!
  PMA=FERT(9,I,NY,NX)
  PMB=FERT(10,I,NY,NX)
  PHA=FERT(11,I,NY,NX)
!
!     LIME AND GYPSUM
!
!     CAC,CAS=CaCO3,CaSO4
!
  CAC=FERT(12,I,NY,NX)
  CAS=FERT(13,I,NY,NX)
!
!     PLANT(1) AND ANIMAL(2) RESIDUE C, N AND P
!
  OFC(1)=FERT(14,I,NY,NX)
  OFN(1)=FERT(15,I,NY,NX)
  OFP(1)=FERT(16,I,NY,NX)
  OFC(2)=FERT(17,I,NY,NX)
  OFN(2)=FERT(18,I,NY,NX)
  OFP(2)=FERT(19,I,NY,NX)
!
!     SOIL LAYER NUMBER AT DEPTH OF FERTILIZER APPLICATION
!
!     LFDPTH=layer number
!     CVRDF=fraction of fertilizer applied to surface litter
!
  IF(Z4A+Z3A+ZUA+ZOA+Z4B+Z3B+ZUB+ZOB+PMA+PMB+PHA+CAC+CAS.GT.0._r8)THEN
    FDPTHF=FDPTH(I,NY,NX)+CumDepz2LayBottom_vr(NU(NY,NX)-1,NY,NX)
    IF(FDPTHF.LE.0._r8.AND.isclose(Z4B+Z3B+ZUB+ZOB+PMB,0._r8))THEN
      LFDPTH=0
      CVRDF=1.0_r8-EXP(-0.8E-02_r8*(SoilOrgM_vr(ielmc,0,NY,NX)/AREA(3,0,NY,NX)))
    ELSE
      D65: DO L=NUI(NY,NX),JZ
        IF(CumDepz2LayBottom_vr(L,NY,NX).GE.FDPTHF)THEN
          LFDPTH=L
          CVRDF=1.0_r8
          exit
        ENDIF
      ENDDO D65
    ENDIF

    BAREF=1.0_r8-CVRDF
!
!     RESET WIDTH AND DEPTH OF NH4 FERTILIZER BAND IF NEW BAND
!     AND ADD REMAINS OF ANY EXISTING FERTILIZER BAND TO NEW BAND
!
!     ROWN=width of NH4 band row
!     DPNHB,WDNHB=depth,width of NH4 band
!     VLNHB,VLNH4=soil volume in NH4 band,non-band
!
    IF((Z4B+Z3B+ZUB.GT.0._r8).OR.((trcs_solml_vr(ids_NH4B,LFDPTH,NY,NX).GT.0._r8 &
      .OR.trcs_solml_vr(idg_NH3B,LFDPTH,NY,NX).GT.0._r8).AND.IFNHB(NY,NX).EQ.0))THEN
      IFNHB(NY,NX)=1
      ROWSpaceNH4_col(NY,NX)=ROWI(I,NY,NX)
      D50: DO L=NUI(NY,NX),JZ
        IF(L.LT.LFDPTH)THEN
          BandThicknessNH4_vr(L,NY,NX)=DLYR_3D(3,L,NY,NX)
          BandWidthNH4_vr(L,NY,NX)=0._r8
        ELSEIF(L.EQ.LFDPTH)THEN
          BandThicknessNH4_vr(L,NY,NX)=AMAX1(0.025_r8,FDPTHF-CumDepz2LayBottom_vr(L-1,NY,NX))
          BandWidthNH4_vr(L,NY,NX)=AMIN1(0.025_r8,ROWSpaceNH4_col(NY,NX))
        ELSE
          BandThicknessNH4_vr(L,NY,NX)=0._r8
          BandWidthNH4_vr(L,NY,NX)=0._r8
        ENDIF
        IF(DLYR_3D(3,L,NY,NX).GT.ZERO2)THEN
          trcs_VLN_vr(ids_NH4B,L,NY,NX)=AMIN1(0.999_r8,BandWidthNH4_vr(L,NY,NX)/ROWSpaceNH4_col(NY,NX) &
            *BandThicknessNH4_vr(L,NY,NX)/DLYR_3D(3,L,NY,NX))
        ELSE
          trcs_VLN_vr(ids_NH4B,L,NY,NX)=0._r8
        ENDIF
        trcs_VLN_vr(ids_NH4,L,NY,NX)    = 1.0_r8-trcs_VLN_vr(ids_NH4B,L,NY,NX)
        trcs_VLN_vr(idg_NH3B,L,NY,NX)   = trcs_VLN_vr(ids_NH4B,L,NY,NX)
        trcs_VLN_vr(idg_NH3,L,NY,NX)    = trcs_VLN_vr(ids_NH4,L,NY,NX)
        ZNH4T                           = trcs_solml_vr(ids_NH4,L,NY,NX)+trcs_solml_vr(ids_NH4B,L,NY,NX)
        ZNH3T                           = trcs_solml_vr(idg_NH3,L,NY,NX)+trcs_solml_vr(idg_NH3B,L,NY,NX)
        XN4T                            = trcx_solml_vr(idx_NH4,L,NY,NX)+trcx_solml_vr(idx_NH4B,L,NY,NX)
        trcs_solml_vr(ids_NH4,L,NY,NX)  = ZNH4T*trcs_VLN_vr(ids_NH4,L,NY,NX)
        trcs_solml_vr(idg_NH3,L,NY,NX)  = ZNH3T*trcs_VLN_vr(idg_NH3,L,NY,NX)
        trcs_solml_vr(ids_NH4B,L,NY,NX) = ZNH4T*trcs_VLN_vr(ids_NH4B,L,NY,NX)
        trcs_solml_vr(idg_NH3B,L,NY,NX) = ZNH3T*trcs_VLN_vr(idg_NH3B,L,NY,NX)
        trcx_solml_vr(idx_NH4,L,NY,NX)  = XN4T*trcs_VLN_vr(ids_NH4,L,NY,NX)
        trcx_solml_vr(idx_NH4B,L,NY,NX) = XN4T*trcs_VLN_vr(ids_NH4B,L,NY,NX)
      ENDDO D50
      BandDepthNH4_col(NY,NX)=BandThicknessNH4_vr(LFDPTH,NY,NX)+CumDepz2LayBottom_vr(LFDPTH-1,NY,NX)
    ENDIF
!
!     RESET WIDTH AND DEPTH OF NO3 FERTILIZER BAND IF NEW BAND
!     AND ADD REMAINS OF ANY EXISTING FERTILIZER BAND TO NEW BAND
!
!     ROWO=width of NO3 band row
!     DPNOB,WDNOB=depth,width of NO3 band
!     VLNOB,VLNO3=soil volume in NO3 band,non-band
!
    IF((Z4B+Z3B+ZUB+ZOB.GT.0._r8).OR.((trcs_solml_vr(ids_NO3B,LFDPTH,NY,NX).GT.0._r8 &
      .OR.trcs_solml_vr(ids_NO2B,LFDPTH,NY,NX).GT.0._r8).AND.IFNOB(NY,NX).EQ.0))THEN
      IFNOB(NY,NX)=1
      ROWSpaceNO3_col(NY,NX)=ROWI(I,NY,NX)
      D45: DO L=NUI(NY,NX),JZ
        IF(L.LT.LFDPTH)THEN
          BandThicknessNO3_vr(L,NY,NX)=DLYR_3D(3,L,NY,NX)
          BandWidthNO3_vr(L,NY,NX)=0._r8
        ELSEIF(L.EQ.LFDPTH)THEN
          BandThicknessNO3_vr(L,NY,NX)=AMAX1(0.01_r8,FDPTHF-CumDepz2LayBottom_vr(L-1,NY,NX))
          BandWidthNO3_vr(L,NY,NX)=AMIN1(0.01_r8,ROWSpaceNO3_col(NY,NX))
        ELSE
          BandThicknessNO3_vr(L,NY,NX)=0._r8
          BandWidthNO3_vr(L,NY,NX)=0._r8
        ENDIF
        IF(DLYR_3D(3,L,NY,NX).GT.ZERO2)THEN
          trcs_VLN_vr(ids_NO3B,L,NY,NX)=AMIN1(0.999_r8,BandWidthNO3_vr(L,NY,NX)/ROWSpaceNO3_col(NY,NX) &
            *BandThicknessNO3_vr(L,NY,NX)/DLYR_3D(3,L,NY,NX))
        ELSE
          trcs_VLN_vr(ids_NO3B,L,NY,NX)=0._r8
        ENDIF

        trcs_VLN_vr(ids_NO3,L,NY,NX)=1.0_r8-trcs_VLN_vr(ids_NO3B,L,NY,NX)
        trcs_VLN_vr(ids_NO2B,L,NY,NX)=trcs_VLN_vr(ids_NO3B,L,NY,NX)
        trcs_VLN_vr(ids_NO2,L,NY,NX)=trcs_VLN_vr(ids_NO3,L,NY,NX)
        ZNO3T=trcs_solml_vr(ids_NO3,L,NY,NX)+trcs_solml_vr(ids_NO3B,L,NY,NX)
        ZNO2T=trcs_solml_vr(ids_NO2,L,NY,NX)+trcs_solml_vr(ids_NO2B,L,NY,NX)

        trcs_solml_vr(ids_NO3,L,NY,NX)=ZNO3T*trcs_VLN_vr(ids_NO3,L,NY,NX)
        trcs_solml_vr(ids_NO2,L,NY,NX)=ZNO2T*trcs_VLN_vr(ids_NO2,L,NY,NX)
        trcs_solml_vr(ids_NO3B,L,NY,NX)=ZNO3T*trcs_VLN_vr(ids_NO3B,L,NY,NX)
        trcs_solml_vr(ids_NO2B,L,NY,NX)=ZNO2T*trcs_VLN_vr(ids_NO2B,L,NY,NX)
      ENDDO D45
      BandDepthNO3_col(NY,NX)=BandThicknessNO3_vr(LFDPTH,NY,NX)+CumDepz2LayBottom_vr(LFDPTH-1,NY,NX)
    ENDIF
!
!     RESET WIDTH AND DEPTH OF PO4 FERTILIZER BAND IF NEW BAND
!     AND ADD REMAINS OF ANY EXISTING FERTILIZER BAND TO NEW BAND
!
!     ROWP=width of H2PO4 band row
!     DPPOB,WDPOB=depth,width of H2PO4 band
!     VLPOB,VLPO4=soil volume in H2PO4 band,non-band
!
    IF((PMB.GT.0.0).OR.(trcs_solml_vr(ids_H2PO4B,LFDPTH,NY,NX).GT.0._r8.AND.IFPOB(NY,NX).EQ.0))THEN
      IFPOB(NY,NX)=1
      ROWSpacePO4_col(NY,NX)=ROWI(I,NY,NX)
      DO  L=NUI(NY,NX),JZ
        IF(L.LT.LFDPTH)THEN
          BandThicknessPO4_vr(L,NY,NX) = DLYR_3D(3,L,NY,NX)
          BandWidthPO4_vr(L,NY,NX)     = AMIN1(0.01,ROWSpacePO4_col(NY,NX))
        ELSEIF(L.EQ.LFDPTH)THEN
          BandThicknessPO4_vr(L,NY,NX) = AMAX1(0.01,FDPTHF-CumDepz2LayBottom_vr(L-1,NY,NX))
          BandWidthPO4_vr(L,NY,NX)     = AMIN1(0.01,ROWSpacePO4_col(NY,NX))
        ELSE
          BandThicknessPO4_vr(L,NY,NX) = 0._r8
          BandWidthPO4_vr(L,NY,NX)     = 0._r8
        ENDIF
        IF(DLYR_3D(3,L,NY,NX).GT.ZERO2)THEN
          trcs_VLN_vr(ids_H1PO4B,L,NY,NX)=AMIN1(0.999_r8,BandWidthPO4_vr(L,NY,NX)/ROWSpacePO4_col(NY,NX) &
          *BandThicknessPO4_vr(L,NY,NX)/DLYR_3D(3,L,NY,NX))
        ELSE
          trcs_VLN_vr(ids_H1PO4B,L,NY,NX)=0._r8
        ENDIF
        trcs_VLN_vr(ids_H1PO4,L,NY,NX)  = 1.0-trcs_VLN_vr(ids_H1PO4B,L,NY,NX)
        trcs_VLN_vr(ids_H2PO4B,L,NY,NX) = trcs_VLN_vr(ids_H1PO4B,L,NY,NX)
        trcs_VLN_vr(ids_H2PO4,L,NY,NX)  = trcs_VLN_vr(ids_H1PO4,L,NY,NX)

        H1PO4T = trcs_solml_vr(ids_H1PO4,L,NY,NX)+trcs_solml_vr(ids_H1PO4B,L,NY,NX)
        H2PO4T = trcs_solml_vr(ids_H2PO4,L,NY,NX)+trcs_solml_vr(ids_H2PO4B,L,NY,NX)
        XOH0T  = trcx_solml_vr(idx_OHe,L,NY,NX)+trcx_solml_vr(idx_OHeB,L,NY,NX)
        XOH1T  = trcx_solml_vr(idx_OH,L,NY,NX)+trcx_solml_vr(idx_OHB,L,NY,NX)
        XOH2T  = trcx_solml_vr(idx_OHp,L,NY,NX)+trcx_solml_vr(idx_OHpB,L,NY,NX)
        XH1PT  = trcx_solml_vr(idx_HPO4,L,NY,NX)+trcx_solml_vr(idx_HPO4B,L,NY,NX)
        XH2PT  = trcx_solml_vr(idx_H2PO4,L,NY,NX)+trcx_solml_vr(idx_H2PO4B,L,NY,NX)
        PALPOT = trcp_saltpml_vr(idsp_AlPO4,L,NY,NX)+trcp_saltpml_vr(idsp_AlPO4B,L,NY,NX)
        PFEPOT = trcp_saltpml_vr(idsp_FePO4,L,NY,NX)+trcp_saltpml_vr(idsp_FePO4B,L,NY,NX)
        PCAPDT = trcp_saltpml_vr(idsp_CaHPO4,L,NY,NX)+trcp_saltpml_vr(idsp_CaHPO4B,L,NY,NX)
        PCAPHT = trcp_saltpml_vr(idsp_HA,L,NY,NX)+trcp_saltpml_vr(idsp_HAB,L,NY,NX)
        PCAPMT = trcp_saltpml_vr(idsp_CaH4P2O8,L,NY,NX)+trcp_saltpml_vr(idsp_CaH4P2O8B,L,NY,NX)

        trcs_solml_vr(ids_H1PO4,L,NY,NX)  = H1PO4T*trcs_VLN_vr(ids_H1PO4,L,NY,NX)
        trcs_solml_vr(ids_H2PO4,L,NY,NX)  = H2PO4T*trcs_VLN_vr(ids_H1PO4,L,NY,NX)
        trcs_solml_vr(ids_H1PO4B,L,NY,NX) = H1PO4T*trcs_VLN_vr(ids_H1PO4B,L,NY,NX)
        trcs_solml_vr(ids_H2PO4B,L,NY,NX) = H2PO4T*trcs_VLN_vr(ids_H1PO4B,L,NY,NX)

        trcx_solml_vr(idx_OHe,L,NY,NX)          = XOH0T*trcs_VLN_vr(ids_H1PO4,L,NY,NX)
        trcx_solml_vr(idx_OH,L,NY,NX)           = XOH1T*trcs_VLN_vr(ids_H1PO4,L,NY,NX)
        trcx_solml_vr(idx_OHp,L,NY,NX)          = XOH2T*trcs_VLN_vr(ids_H1PO4,L,NY,NX)
        trcx_solml_vr(idx_HPO4,L,NY,NX)         = XH1PT*trcs_VLN_vr(ids_H1PO4,L,NY,NX)
        trcx_solml_vr(idx_H2PO4,L,NY,NX)        = XH2PT*trcs_VLN_vr(ids_H1PO4,L,NY,NX)
        trcx_solml_vr(idx_OHeB,L,NY,NX)         = XOH0T*trcs_VLN_vr(ids_H1PO4B,L,NY,NX)
        trcx_solml_vr(idx_OHB,L,NY,NX)          = XOH1T*trcs_VLN_vr(ids_H1PO4B,L,NY,NX)
        trcx_solml_vr(idx_OHpB,L,NY,NX)         = XOH2T*trcs_VLN_vr(ids_H1PO4B,L,NY,NX)
        trcx_solml_vr(idx_HPO4B,L,NY,NX)        = XH1PT*trcs_VLN_vr(ids_H1PO4B,L,NY,NX)
        trcx_solml_vr(idx_H2PO4B,L,NY,NX)       = XH2PT*trcs_VLN_vr(ids_H1PO4B,L,NY,NX)
        trcp_saltpml_vr(idsp_AlPO4,L,NY,NX)     = PALPOT*trcs_VLN_vr(ids_H1PO4,L,NY,NX)
        trcp_saltpml_vr(idsp_FePO4,L,NY,NX)     = PFEPOT*trcs_VLN_vr(ids_H1PO4,L,NY,NX)
        trcp_saltpml_vr(idsp_CaHPO4,L,NY,NX)    = PCAPDT*trcs_VLN_vr(ids_H1PO4,L,NY,NX)
        trcp_saltpml_vr(idsp_HA,L,NY,NX)        = PCAPHT*trcs_VLN_vr(ids_H1PO4,L,NY,NX)
        trcp_saltpml_vr(idsp_CaH4P2O8,L,NY,NX)  = PCAPMT*trcs_VLN_vr(ids_H1PO4,L,NY,NX)
        trcp_saltpml_vr(idsp_AlPO4B,L,NY,NX)    = PALPOT*trcs_VLN_vr(ids_H1PO4B,L,NY,NX)
        trcp_saltpml_vr(idsp_FePO4B,L,NY,NX)    = PFEPOT*trcs_VLN_vr(ids_H1PO4B,L,NY,NX)
        trcp_saltpml_vr(idsp_CaHPO4B,L,NY,NX)   = PCAPDT*trcs_VLN_vr(ids_H1PO4B,L,NY,NX)
        trcp_saltpml_vr(idsp_HAB,L,NY,NX)       = PCAPHT*trcs_VLN_vr(ids_H1PO4B,L,NY,NX)
        trcp_saltpml_vr(idsp_CaH4P2O8B,L,NY,NX) = PCAPMT*trcs_VLN_vr(ids_H1PO4B,L,NY,NX)

        if(salt_model)then
          H0PO4T = trcSalt_solml_vr(idsalt_H0PO4,L,NY,NX)+trcSalt_solml_vr(idsalt_H0PO4B,L,NY,NX)
          H3PO4T = trcSalt_solml_vr(idsalt_H3PO4,L,NY,NX)+trcSalt_solml_vr(idsalt_H3PO4B,L,NY,NX)
          ZFE1PT = trcSalt_solml_vr(idsalt_FeHPO4,L,NY,NX)+trcSalt_solml_vr(idsalt_FeHPO4B,L,NY,NX)
          ZFE2PT = trcSalt_solml_vr(idsalt_FeH2PO4,L,NY,NX)+trcSalt_solml_vr(idsalt_FeH2PO4B,L,NY,NX)
          ZCA0PT = trcSalt_solml_vr(idsalt_CaPO4,L,NY,NX)+trcSalt_solml_vr(idsalt_CaPO4B,L,NY,NX)
          ZCA1PT = trcSalt_solml_vr(idsalt_CaHPO4,L,NY,NX)+trcSalt_solml_vr(idsalt_CaHPO4B,L,NY,NX)
          ZCA2PT = trcSalt_solml_vr(idsalt_CaH4P2O8,L,NY,NX)+trcSalt_solml_vr(idsalt_CaH4P2O8B,L,NY,NX)
          ZMG1PT = trcSalt_solml_vr(idsalt_MgHPO4,L,NY,NX)+trcSalt_solml_vr(idsalt_MgHPO4B,L,NY,NX)
          
          trcSalt_solml_vr(idsalt_H0PO4,L,NY,NX)=H0PO4T*trcs_VLN_vr(ids_H1PO4,L,NY,NX)
          trcSalt_solml_vr(idsalt_H3PO4,L,NY,NX)=H3PO4T*trcs_VLN_vr(ids_H1PO4,L,NY,NX)
          trcSalt_solml_vr(idsalt_FeHPO4,L,NY,NX)=ZFE1PT*trcs_VLN_vr(ids_H1PO4,L,NY,NX)
          trcSalt_solml_vr(idsalt_FeH2PO4,L,NY,NX)=ZFE2PT*trcs_VLN_vr(ids_H1PO4,L,NY,NX)
          trcSalt_solml_vr(idsalt_CaPO4,L,NY,NX)=ZCA0PT*trcs_VLN_vr(ids_H1PO4,L,NY,NX)
          trcSalt_solml_vr(idsalt_CaHPO4,L,NY,NX)=ZCA1PT*trcs_VLN_vr(ids_H1PO4,L,NY,NX)
          trcSalt_solml_vr(idsalt_CaH4P2O8,L,NY,NX)=ZCA2PT*trcs_VLN_vr(ids_H1PO4,L,NY,NX)
          trcSalt_solml_vr(idsalt_MgHPO4,L,NY,NX)=ZMG1PT*trcs_VLN_vr(ids_H1PO4,L,NY,NX)
          trcSalt_solml_vr(idsalt_H0PO4B,L,NY,NX)=H0PO4T*trcs_VLN_vr(ids_H1PO4B,L,NY,NX)

          trcSalt_solml_vr(idsalt_H3PO4B,L,NY,NX)=H3PO4T*trcs_VLN_vr(ids_H1PO4B,L,NY,NX)
          trcSalt_solml_vr(idsalt_FeHPO4B,L,NY,NX)=ZFE1PT*trcs_VLN_vr(ids_H1PO4B,L,NY,NX)
          trcSalt_solml_vr(idsalt_FeH2PO4B,L,NY,NX)=ZFE2PT*trcs_VLN_vr(ids_H1PO4B,L,NY,NX)
          trcSalt_solml_vr(idsalt_CaPO4B,L,NY,NX)=ZCA0PT*trcs_VLN_vr(ids_H1PO4B,L,NY,NX)
          trcSalt_solml_vr(idsalt_CaHPO4B,L,NY,NX)=ZCA1PT*trcs_VLN_vr(ids_H1PO4B,L,NY,NX)
          trcSalt_solml_vr(idsalt_CaH4P2O8B,L,NY,NX)=ZCA2PT*trcs_VLN_vr(ids_H1PO4B,L,NY,NX)
          trcSalt_solml_vr(idsalt_MgHPO4B,L,NY,NX)=ZMG1PT*trcs_VLN_vr(ids_H1PO4B,L,NY,NX)
        endif

      ENDDO
      BandDepthPO4_col(NY,NX)=BandThicknessPO4_vr(LFDPTH,NY,NX)+CumDepz2LayBottom_vr(LFDPTH-1,NY,NX)
    ENDIF
!
!     UPDATE STATE VARIABLES FOR BROADCAST AND BANDED FERTILIZER
!     NH4, NH3, UREA, NO3, PO4, LIME AND GYPSUM IN SOIL
!     AND CONVERT FROM G TO MOLE
!
!     ZNH4FA,ZNH3FA,ZNHUFA,ZNO3FA=bdcast NH4,NH3,urea,NO3 fertilizer
!     ZNH4FB,ZNH3FB,ZNHUFB,ZNO3FB=banded NH4,NH3,urea,NO3 fertilizer
!     PCAPM1,PCAPD1,PCAPH1=concn of precip CaH4P2O8,CaHPO4,apatite in non-band
!     PCAPMB,PCAPDB,PCAPHB=concn of precip CaH4P2O8,CaHPO4,apatite in band
!     PCACO,PCASO=precipitated CaCO3,CaSO4
!
    Z4AX = Z4A*AREA(3,LFDPTH,NY,NX)/natomw
    Z3AX = Z3A*AREA(3,LFDPTH,NY,NX)/natomw
    ZUAX = ZUA*AREA(3,LFDPTH,NY,NX)/natomw
    ZOAX = ZOA*AREA(3,LFDPTH,NY,NX)/natomw
    Z4BX = Z4B*AREA(3,LFDPTH,NY,NX)/natomw
    Z3BX = Z3B*AREA(3,LFDPTH,NY,NX)/natomw
    ZUBX = ZUB*AREA(3,LFDPTH,NY,NX)/natomw
    ZOBX = ZOB*AREA(3,LFDPTH,NY,NX)/natomw
    PMAX = PMA*AREA(3,LFDPTH,NY,NX)/(2.0_r8*patomw)
    PMBX = PMB*AREA(3,LFDPTH,NY,NX)/(2.0_r8*patomw)
    PHAX = PHA*AREA(3,LFDPTH,NY,NX)/(3.0_r8*patomw)
    CACX = CAC*AREA(3,LFDPTH,NY,NX)/40._r8
    CASX = CAS*AREA(3,LFDPTH,NY,NX)/40._r8

    FertN_mole_soil_vr(ifert_nh4,LFDPTH,NY,NX)  = AZMAX1(FertN_mole_soil_vr(ifert_nh4,LFDPTH,NY,NX)+Z4AX*CVRDF)
    FertN_mole_soil_vr(ifert_urea,LFDPTH,NY,NX) = AZMAX1(FertN_mole_soil_vr(ifert_urea,LFDPTH,NY,NX)+ZUAX*CVRDF)
    FertN_mole_soil_vr(ifert_no3,LFDPTH,NY,NX)  = AZMAX1(FertN_mole_soil_vr(ifert_no3,LFDPTH,NY,NX)+ZOAX*CVRDF)

    FertN_mole_Band_vr(ifert_nh4_band,LFDPTH,NY,NX)  = AZMAX1(FertN_mole_Band_vr(ifert_nh4_band,LFDPTH,NY,NX)+Z4BX*CVRDF)
    FertN_mole_Band_vr(ifert_urea_band,LFDPTH,NY,NX) = AZMAX1(FertN_mole_Band_vr(ifert_urea_band,LFDPTH,NY,NX)+ZUBX*CVRDF)
    FertN_mole_Band_vr(ifert_no3_band,LFDPTH,NY,NX)  = AZMAX1(FertN_mole_Band_vr(ifert_no3_band,LFDPTH,NY,NX)+ZOBX*CVRDF)

    trcp_saltpml_vr(idsp_CaH4P2O8,LFDPTH,NY,NX)=trcp_saltpml_vr(idsp_CaH4P2O8,LFDPTH,NY,NX)+PMAX*trcs_VLN_vr(ids_H1PO4,LFDPTH,NY,NX)*CVRDF
    trcp_saltpml_vr(idsp_CaH4P2O8B,LFDPTH,NY,NX)=trcp_saltpml_vr(idsp_CaH4P2O8B,LFDPTH,NY,NX)+PMAX*trcs_VLN_vr(ids_H1PO4B,LFDPTH,NY,NX)*CVRDF+PMBX*CVRDF
    trcp_saltpml_vr(idsp_HA,LFDPTH,NY,NX)=trcp_saltpml_vr(idsp_HA,LFDPTH,NY,NX)+PHAX*trcs_VLN_vr(ids_H1PO4,LFDPTH,NY,NX)*CVRDF
    trcp_saltpml_vr(idsp_HAB,LFDPTH,NY,NX)=trcp_saltpml_vr(idsp_HAB,LFDPTH,NY,NX)+PHAX*trcs_VLN_vr(ids_H1PO4B,LFDPTH,NY,NX)*CVRDF
    IF(LFDPTH.EQ.0)THEN
      FertN_mole_soil_vr(ifert_nh4,NU(NY,NX),NY,NX)=AZMAX1(FertN_mole_soil_vr(ifert_nh4,NU(NY,NX),NY,NX)+Z4AX*BAREF)
      FertN_mole_soil_vr(ifert_nh3,NU(NY,NX),NY,NX)=AZMAX1(FertN_mole_soil_vr(ifert_nh3,NU(NY,NX),NY,NX)+Z3AX)
      FertN_mole_soil_vr(ifert_urea,NU(NY,NX),NY,NX)=AZMAX1(FertN_mole_soil_vr(ifert_urea,NU(NY,NX),NY,NX)+ZUAX*BAREF)
      FertN_mole_soil_vr(ifert_no3,NU(NY,NX),NY,NX)=AZMAX1(FertN_mole_soil_vr(ifert_no3,NU(NY,NX),NY,NX)+ZOAX*BAREF)

      FertN_mole_Band_vr(ifert_nh4_band,NU(NY,NX),NY,NX)=AZMAX1(FertN_mole_Band_vr(ifert_nh4_band,NU(NY,NX),NY,NX)+Z4BX*BAREF)
      FertN_mole_Band_vr(ifert_nh3_band,NU(NY,NX),NY,NX)=AZMAX1(FertN_mole_Band_vr(ifert_nh3_band,NU(NY,NX),NY,NX)+Z3BX)
      FertN_mole_Band_vr(ifert_urea_band,NU(NY,NX),NY,NX)=AZMAX1(FertN_mole_Band_vr(ifert_urea_band,NU(NY,NX),NY,NX)+ZUBX*BAREF)
      FertN_mole_Band_vr(ifert_no3_band,NU(NY,NX),NY,NX)=AZMAX1(FertN_mole_Band_vr(ifert_no3_band,NU(NY,NX),NY,NX)+ZOBX*BAREF)

      trcp_saltpml_vr(idsp_CaH4P2O8,NU(NY,NX),NY,NX)=trcp_saltpml_vr(idsp_CaH4P2O8,NU(NY,NX),NY,NX)+PMAX*trcs_VLN_vr(ids_H1PO4,NU(NY,NX),NY,NX)*BAREF
      trcp_saltpml_vr(idsp_CaH4P2O8B,NU(NY,NX),NY,NX)=trcp_saltpml_vr(idsp_CaH4P2O8B,NU(NY,NX),NY,NX)+PMAX*trcs_VLN_vr(ids_H1PO4B,NU(NY,NX),NY,NX)*BAREF+PMBX*BAREF
      trcp_saltpml_vr(idsp_HA,NU(NY,NX),NY,NX)=trcp_saltpml_vr(idsp_HA,NU(NY,NX),NY,NX)+PHAX*trcs_VLN_vr(ids_H1PO4,NU(NY,NX),NY,NX)*BAREF
      trcp_saltpml_vr(idsp_HAB,NU(NY,NX),NY,NX)=trcp_saltpml_vr(idsp_HAB,NU(NY,NX),NY,NX)+PHAX*trcs_VLN_vr(ids_H1PO4B,NU(NY,NX),NY,NX)*BAREF
    ELSE
      FertN_mole_soil_vr(ifert_nh3,LFDPTH,NY,NX)      = AZMAX1(FertN_mole_soil_vr(ifert_nh3,LFDPTH,NY,NX)+Z3AX*CVRDF)
      FertN_mole_Band_vr(ifert_nh3_band,LFDPTH,NY,NX) = AZMAX1(FertN_mole_Band_vr(ifert_nh3_band,LFDPTH,NY,NX)+Z3BX*CVRDF)
    ENDIF
    trcp_saltpml_vr(idsp_CaCO3,NU(NY,NX),NY,NX)=trcp_saltpml_vr(idsp_CaCO3,NU(NY,NX),NY,NX)+CACX
    trcp_saltpml_vr(idsp_CaSO4,NU(NY,NX),NY,NX)=trcp_saltpml_vr(idsp_CaSO4,NU(NY,NX),NY,NX)+CASX
    TZIN   = TZIN+natomw*(Z4AX+Z3AX+ZUAX+ZOAX+Z4BX+Z3BX+ZUBX+ZOBX)
    TPIN   = TPIN+62.0_r8*(PMAX+PMBX)+93.0_r8*PHAX
    TIONIN = TIONIN+2.0_r8*(CACX+CASX)
    FertNFlx_CumYr_col(NY,NX)=FertNFlx_CumYr_col(NY,NX)+natomw*(Z4AX+Z4BX+Z3AX+Z3BX+ZUAX+ZUBX+ZOAX+ZOBX)
    FerPFlx_CumYr_col(NY,NX)=FerPFlx_CumYr_col(NY,NX)+62.0_r8*(PMAX+PMBX)+93.0_r8*PHAX
  ENDIF

  call PrintInfo('end '//subname)  
  end subroutine ApplyMineralFertilizer

end module FertilizerMod