module PlantDisturbMod
!
!! Description:
! code to apply distance to plants
  use data_kind_mod, only : r8 => SHR_KIND_R8
  use minimathmod, only : test_aeqb
  use SOMDataType
  use GrosubPars
  use PlantTraitDataType
  use GridConsts
  use FlagDataType
  use SoilWaterDataType
  use EcoSIMCtrlDataType
  use PlantDataRateType
  use FertilizerDataType
  use ClimForcDataType
  use EcosimConst
  use PlantMngmtDataType
  use RootDataType
  use CanopyDataType
  use EcoSimSumDataType
  use SoilBGCDataType
  use EcosimBGCFluxType
  use GridDataType
  implicit none
  private
  character(len=*),private, parameter :: mod_filename = __FILE__
! end_include_section

  public :: PrepLandscapeGrazing
  contains

!------------------------------------------------------------------------------------------

  subroutine PrepLandscapeGrazing(I,J,NHW,NHE,NVN,NVS)
  implicit none
  integer, intent(in) :: I, J
  integer, intent(in) :: NHW,NHE,NVN,NVS
  integer :: NX,NY,NZ,nn,nx1,NY1
  real(r8) :: WTSHTZ
!     begin_execution

  DO 2995 NX=NHW,NHE
    DO 2990 NY=NVN,NVS
      DO 2985 NZ=1,NP(NY,NX)
!
!     IHVST=harvest type:0=none,1=grain,2=all above-ground
!                       ,3=pruning,4=grazing,5=fire,6=herbivory
!     LSG=landscape grazing section number
!     WTSHTZ,WTSHTA=total,average biomass in landscape grazing section
!
        IF(IHVST(NZ,I,NY,NX).EQ.4.OR.IHVST(NZ,I,NY,NX).EQ.6)THEN
          WTSHTZ=0
          NN=0
          DO 1995 NX1=NHW,NHE
            DO 1990 NY1=NVN,NVS
              IF(LSG(NZ,NY1,NX1).EQ.LSG(NZ,NY,NX))THEN
                IF(IFLGC(NZ,NY1,NX1).EQ.1)THEN
                  WTSHTZ=WTSHTZ+WTSHT(NZ,NY1,NX1)
                  NN=NN+1
                ENDIF
              ENDIF
1990        CONTINUE
1995      CONTINUE
          IF(NN.GT.0)THEN
            WTSHTA(NZ,NY,NX)=WTSHTZ/NN
          ELSE
            WTSHTA(NZ,NY,NX)=WTSHT(NZ,NY,NX)
          ENDIF
        ENDIF
2985  CONTINUE
2990  CONTINUE
2995  CONTINUE
  end subroutine PrepLandscapeGrazing
end module PlantDisturbMod
