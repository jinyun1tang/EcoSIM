module ForcTypeMod
  use data_kind_mod     , only : r8 => SHR_KIND_R8
implicit none

  character(len=*),private, parameter :: mod_filename = __FILE__
  type, public :: forc_type
!primary variables
    real(r8) :: CCH4E       !atmospheric CH4 concentration, [g m-3]
    real(r8) :: COXYE       !atmospheric O2 concentration, [g m-3]
    real(r8) :: COXQ        !surface irrigation  O2 concentration, [g m-3]
    real(r8) :: COXR        !precipitation  O2 concentration, [g m-3]
    real(r8) :: FLQRI       !irrigation flux into surface litter, [m3 d-2 h-1]
    real(r8) :: FLQRQ       !precipitation flux into surface litter, [m3 d-2 h-1]
    real(r8) :: OFFSET      !offset for calculating temperature in Arrhenius curves, [oC]

    real(r8) :: THETY       !air-dry water content, [m3 m-3]
    real(r8) :: TKS         !temperature in kelvin, [K]
    real(r8) :: THETW       !volumetric water content [m3 m-3]
    real(r8) :: pH          !pH value
    real(r8) :: BKVL        !mass of soil layer	Mg d-2
    real(r8) :: VOLX        !volume of soil layer	m3 d-2
    real(r8) :: VOLW        !soil micropore water content [m3 d-2]

    real(r8) :: VLNOB       !NO3 band volume fracrion, [-]
    real(r8) :: VLNO3       !NO3 non-band volume fraction,[-] := 1-VLNOB
    real(r8) :: VLNHB       !NH4 band volume fraction, [-]
    real(r8) :: VLNH4       !NH4 non-band volume fraction, [-]:=1-VLNHB
    real(r8) :: VLPOB       !PO4 band volume fracrion, [-]
    real(r8) :: VLPO4       !PO4 non-band volume fraction,[-]:=1-VLPOB

    real(r8) :: PSISM       !soil micropore matric water potential, [MPa]
    real(r8) :: OLSGL       !aqueous O2 diffusivity, [m2 h-1], set in hour1
    real(r8) :: ORGC        !total soil organic C [g d-2]
    real(r8) :: CFOMC(1:2)  !allocation coefficient to humus fractions

    real(r8) :: CNH4B       !NH4 concentration band micropore	[g m-3], derived from ZNH4B
    real(r8) :: ZNH4B       !NH4 band micropore, [g d-2]
    real(r8) :: CNH4S       !NH4 concentration non-band micropore	[g m-3], derived from ZNH4S
    real(r8) :: ZNH4S       !NH4 non-band micropore, [g d-2]
    real(r8) :: CNO3B       !NO3 concentration band micropore	[g m-3], derived from ZNO3B
    real(r8) :: ZNO3B       !NO3 band micropore, [g d-2]
    real(r8) :: CNO3S       !NO3 concentration non-band micropore	[g m-3], derived from ZNO3S
    real(r8) :: ZNO3S       !NO3 band micropore, [g d-2]
    real(r8) :: CH2P4       !aqueous PO4 concentration non-band	[g m-3], derived from H2PO4
    real(r8) :: H2PO4       !PO4 non-band micropore, [g d-2]
    real(r8) :: CH2P4B      !aqueous PO4 concentration band	[g m-3], derived from H2POB
    real(r8) :: H2POB       !PO4 band micropore, [g d-2]
    real(r8) :: CH1P4       !aqueous H1PO4 concentration non-band [g m-3], derived from H1PO4
    real(r8) :: H1PO4       !soil aqueous HPO4 content micropore non-band, [mol d-2]
    real(r8) :: CH1P4B      !aqueous H1PO4 concentration band [g m-3], derived from H1POB
    real(r8) :: H1POB       !soil aqueous HPO4 content micropore band, [mol d-2]
    real(r8) :: CNO2B       !aqueous HNO2 concentration band [g m-3], derived from ZNO2B
    real(r8) :: ZNO2B       !NO2  band micropore, [g d-2]
    real(r8) :: CNO2S       !NO2 concentration non-band micropore	[g m-3], derived from ZNO2S
    real(r8) :: ZNO2S       !NO2  non-band micropore, [g d-2]
    real(r8) :: DFGS        !coefficient for dissolution - volatilization, []
    real(r8) :: THETPM      !soil air-filled porosity, [m3 m-3]

!derived variables
    real(r8) :: FILM        !soil water film thickness , [m]
    real(r8) :: TORT        !soil tortuosity, []
    real(r8) :: VOLPM       !soil air content, [m3 d-2]

    real(r8) :: EPOC      !partitioning coefficient between POC and litter, [], hour1.f
    real(r8) :: EHUM      !partitioning coefficient between humus and microbial residue, [], hour1.f
    real(r8) :: CCO2S     !aqueous CO2 concentration micropore	[g m-3]
    real(r8) :: CZ2OS     !aqueous N2O concentration micropore	[g m-3]
    real(r8) :: Z2OS      !aqueous N2O micropore, [g d-2]
    real(r8) :: COXYS     !aqueous O2 concentration micropore	[g m-3]
    real(r8) :: OXYS      !aqueous O2  micropore	[g d-2]
    real(r8) :: COXYG     !gaseous O2 concentration	[g m-3]
    real(r8) :: CZ2GS     !aqueous N2 concentration micropore	[g m-3]
    real(r8) :: CH2GS     !aqueous H2 concentration	[g m-3]
    real(r8) :: H2GS      !aqueous H2 	[g d-2]
    real(r8) :: CCH4G     !gaseous CH4 concentration	[g m-3]
    real(r8) :: CH4S      !aqueous CO2  micropore	[g d-2]
    real(r8) :: SCH4L     !solubility of CH4, [m3 m-3]
    real(r8) :: SOXYL     !solubility of O2, [m3 m-3]
    real(r8) :: ZNFN0     !initial nitrification inhibition activity
    real(r8) :: ZNFNI     !current nitrification inhibition activity

    real(r8) :: VOLY      !micropore volume, [m3 d-2]
    real(r8) :: POROS     !soil porosity, [m3 m-3]
    real(r8) :: FC        !water contents at field capacity, [m3 m-3]
    real(r8) :: TFND      !temperature effect on aqueous diffusivity

 !litter layer
    real(r8) :: VOLR      !surface litter volume, [m3 d-2]
    real(r8) :: VOLWRX    !surface litter water holding capacity, [m3 d-2]
 !non litter layer
!    real(r8) :: ROXYY       !total root + microbial O2 uptake from previous hour, [g d-2 h-1], updated in hour1
!    real(r8) :: RN2OY       !total root + microbial N2O uptake from previous hour, [g d-2 h-1]
!    real(r8) :: RNO2Y       !total root + microbial NO2 uptake non-band from previous hour, [g d-2 h-1]
!    real(r8) :: RN2BY       !total root + microbial NO2 uptake band from previous hour, [g d-2 h-1]
!    real(r8) :: RNHBY       !total root + microbial NH4 uptake band from previous hour, [g d-2 h-1]
!    real(r8) :: RN3BY       !total root + microbial NO3 uptake band from previous hour, [g d-2 h-1]
!    real(r8) :: RPOBY       !total root + microbial PO4 uptake band from previous hour, [g d-2 h-1]
!    real(r8) :: RP1BY       !HPO4 demand in band by all microbial, root, myco populations from previous hour
!    real(r8) :: RNH4Y       !total root + microbial NH4 uptake non-band from previous hour, [g d-2 h-1]
!    real(r8) :: RNO3Y       !total root + microbial NO3 uptake non-band from previous hour, [g d-2 h-1]
!    real(r8) :: RPO4Y       !total root + microbial PO4 uptake non-band from previous hour, [g d-2 h-1]
!    real(r8) :: RP14Y       !HPO4 demand in non-band by all microbial, root, myco populations from previous hour
!    real(r8) :: ROXYF       !net gaseous O2 flux, [g d-2 h-1], updated in redist.f
!    real(r8) :: RCH4L       !net aqueous CH4 flux, [g d-2 h-1], updated in redist.f
!    real(r8) :: ROXYL       !net aqueous O2 flux from previous hour, [g d-2 h-1], updated in redist.f
  end type forc_type

  contains


!------------------------------------------------------------------------------------------

  subroutine ReadForc(forc)

  implicit none
  type(forc_type), intent(inout) :: forc



  end subroutine ReadForc

end module ForcTypeMod
