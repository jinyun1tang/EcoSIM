module GridConsts
  use data_kind_mod, only : r8 => DAT_KIND_R8
implicit none
  public
  save
  character(len=*), private, parameter :: mod_filename = &
  __FILE__

  integer :: JX=6
  integer :: JY=4
  integer :: JX0,JY0
  integer :: JZ=20
  integer :: JH
  integer :: JV
  integer :: JD
  integer :: MaxNumBranches=10      !number of plant branches
  integer :: NumGrowthStages     !number of plant growth states
  integer :: NumOfPlantMorphUnits !
  integer :: MaxNumRootAxes     !maximum number of root layers
  integer :: NumLitterGroups !number of liter groups
  integer, PARAMETER :: JP=5    !maximum pft in a given topgraphic column
  integer, PARAMETER :: NumOfCanopyLayers=10   !# of canopy layers
  integer, PARAMETER :: JS=5
  integer, PARAMETER :: NumOfLeafZenithSectors=4   !# of sectors for the leaf zenith [0,pi/2]
  integer, PARAMETER :: NumOfLeafAzimuthSectors=4   !# of sectors for the leaf azimuth, [0,pi]
  integer, PARAMETER :: NumOfSkyAzimuSects=4   !# of sectors for the sky azimuth  [0,2*pi]
  integer, parameter :: MaxNodesPerBranch=25!# of nodes for plant canopy
  integer, pointer :: JGnio(:)   !guid indices for organic-microbial complex
  integer, pointer :: JGnfo(:)   !guid indices for organic-microbial complex
  integer, pointer :: JGniA(:)   !guid indices for autotrophic-microbial complex
  integer, pointer :: JGnfA(:)   !guid indices for autotrophic-microbial complex
  integer  :: NumMicrobAutrophCmplx            !total number of microbial guilds in the autotrophic complex
  integer  :: NumHetetrMicCmplx             !total number of microbial guilds in one organic-microbial complex
  integer  :: NumLiveHeterBioms         !total number of live biomass component in one heterotroph organo-microbe complex
  integer  :: NumLiveAutoBioms          !total number of live biomass component in one autotrophic complex
  type, public :: bounds_type
   integer :: NHW
   integer :: NVN
   integer :: NHE
   integer :: NVS
   integer :: begg,endg
   integer :: begt,endt
   integer :: begc,endc
   integer :: begp,endp
   integer :: ngrid
   integer :: ntopou
   integer :: ncols
   integer :: npfts
   integer, pointer :: icol(:,:)    => null()   !column id
   integer, pointer :: ipft(:,:,:)  => null()   !pft id
  end type bounds_type

  type(bounds_type) :: bounds

end module GridConsts
