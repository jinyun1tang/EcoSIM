
! This code is adapted from Alquimia with the credits below:
! Alquimia Copyright (c) 2013-2016, The Regents of the University of California,
! through Lawrence Berkeley National Laboratory (subject to receipt of any
! required approvals from the U.S. Dept. of Energy).  All rights reserved.
!
! Alquimia is available under a BSD license. See LICENSE.txt for more
! information.
!
! If you have questions about your rights to use or distribute this software,
! please contact Berkeley Lab's Technology Transfer and Intellectual Property
! Management at TTD@lbl.gov referring to Alquimia (LBNL Ref. 2013-119).
!
! NOTICE.  This software was developed under funding from the U.S. Department
! of Energy.  As such, the U.S. Government has been granted for itself and
! others acting on its behalf a paid-up, nonexclusive, irrevocable, worldwide
! license in the Software to reproduce, prepare derivative works, and perform
! publicly and display publicly.  Beginning five (5) years after the date
! permission to assert copyright is obtained from the U.S. Department of Energy,
! and subject to any subsequent five (5) year renewals, the U.S. Government is
! granted for itself and others acting on its behalf a paid-up, nonexclusive,
! irrevocable, worldwide license in the Software to reproduce, prepare derivative
! works, distribute copies to the public, perform publicly and display publicly,
! and to permit others to do so.
!
! Authors: Benjamin Andre <bandre@lbl.gov>
!        : Andrew Graus <agraus@lbl.gov>
!

! **************************************************************************** !
!
! ATS-EcoSIM Containers module
!
! Author: Andrew Graus
!
! WARNINGS:
!
!   * The data structures defined here are dictated by
!     the API! Do NOT change them unless you make
!     corresponding changes to the c containers (and doc).
!
!   * The order of the data members matters! If num_primary is the
!     first member and num_minerals is the third, then they must be in
!     those positions on both the c and fortran side of the interface!
!
!   * The names of the containers and their members should be the same
!     for both c and fortran. The language interface doesn't require
!     this, but it makes it easier for the reader to understand what is
!     going on.
!
! **************************************************************************** !

module BGCContainers_module

  use, intrinsic :: iso_c_binding

  implicit none

  integer (c_int), parameter :: kBGCMaxStringLength = 512
  integer (c_int), parameter :: kBGCMaxWordLength = 32

  integer (c_int), parameter :: kBGCNoError = 0
  integer (c_int), parameter :: kBGCErrorInvalidEngine = 1
  integer (c_int), parameter :: kBGCErrorUnknownConstraintName = 2
  integer (c_int), parameter :: kBGCErrorUnsupportedFunctionality = 3
  integer (c_int), parameter :: kBGCErrorEngineIntegrity = 4577

  character (13), parameter :: kBGCStringTotalAqueous = 'total_aqueous'
  character (12), parameter :: kBGCStringTotalSorbed = 'total_sorbed'
  character (25), parameter :: kBGCStringTotalAqueousPlusSorbed = 'total_aqueous_plus_sorbed'
  character (4), parameter :: kBGCStringFree = 'free'
  character (2), parameter :: kBGCStringPH = 'pH'
  character (7), parameter :: kBGCStringMineral = 'mineral'
  character (3), parameter :: kBGCStringGas = 'gas'
  character (6), parameter :: kBGCStringCharge = 'charge'

  type, public, bind(c) :: BGCVectorDouble
     integer (c_int) :: size
     integer (c_int) :: capacity
     type (c_ptr) :: data
  end type BGCVectorDouble

  type, public, bind(c) :: BGCVectorInt
     integer (c_int) :: size
     integer (c_int) :: capacity
     type (c_ptr) :: data
  end type BGCVectorInt

  type, public, bind(c) :: BGCVectorString
     integer (c_int) :: size
     integer (c_int) :: capacity
     type (c_ptr) :: data
  end type BGCVectorString

  type, public, bind(c) :: BGCMatrixDouble
    integer (c_int) :: rows
    integer (c_int) :: cols
    integer (c_int) :: cap_rows
    integer (c_int) :: cap_cols
    type (c_ptr) :: data
  end type BGCMatrixDouble

  type, public, bind(c) :: BGCMatrixInt
    integer (c_int) :: rows
    integer (c_int) :: cols
    integer (c_int) :: cap_rows
    integer (c_int) :: cap_cols
    type (c_ptr) :: data
  end type BGCMatrixInt

  type, public, bind(c) :: BGCMatrixString
    integer (c_int) :: rows
    integer (c_int) :: cols
    integer (c_int) :: cap_rows
    integer (c_int) :: cap_cols
    type (c_ptr) :: data
  end type BGCMatrixString

  type, public, bind(c) :: BGCTensorDouble
    integer (c_int) :: rows
    integer (c_int) :: cols
    integer (c_int) :: procs
    integer (c_int) :: cap_rows
    integer (c_int) :: cap_cols
    integer (c_int) :: cap_procs
    type (c_ptr) :: data
  end type BGCTensorDouble

  type, public, bind(c) :: BGCTensorInt
    integer (c_int) :: rows
    integer (c_int) :: cols
    integer (c_int) :: procs
    integer (c_int) :: cap_rows
    integer (c_int) :: cap_cols
    integer (c_int) :: cap_procs
    type (c_ptr) :: data
  end type BGCTensorInt

  type, public, bind(c) :: BGCSizes
     integer (c_int) :: ncells_per_col_
     integer (c_int) :: num_components
     integer (c_int) :: num_columns
  end type BGCSizes

  type, public, bind(c) :: BGCState
     ! I think I have to write the data as vector doubles
     type (BGCMatrixDouble) :: liquid_density
     type (BGCMatrixDouble) :: gas_density
     type (BGCMatrixDouble) :: ice_density
     type (BGCMatrixDouble) :: porosity
     type (BGCMatrixDouble) :: water_content
     type (BGCMatrixDouble) :: matric_pressure
     type (BGCMatrixDouble) :: temperature
     type (BGCMatrixDouble) :: hydraulic_conductivity
     type (BGCMatrixDouble) :: bulk_density
     type (BGCMatrixDouble) :: subsurface_water_source
     type (BGCMatrixDouble) :: subsurface_energy_source
     type (BGCVectorDouble) :: surface_energy_source
     type (BGCVectorDouble) :: surface_water_source
     type (BGCTensorDouble) :: total_component_concentration
  end type BGCState

  type, public, bind(c) :: BGCProperties
     type (BGCMatrixDouble) :: liquid_saturation
     type (BGCMatrixDouble) :: gas_saturation
     type (BGCMatrixDouble) :: ice_saturation
     type (BGCMatrixDouble) :: relative_permeability
     type (BGCMatrixDouble) :: thermal_conductivity
     type (BGCMatrixDouble) :: volume
     type (BGCMatrixDouble) :: depth
     type (BGCMatrixDouble) :: depth_c
     type (BGCMatrixDouble) :: dz
     type (BGCMatrixDouble) :: plant_wilting_factor
     type (BGCMatrixDouble) :: rooting_depth_fraction
     type (BGCVectorDouble)  :: shortwave_radiation
     type (BGCVectorDouble) :: longwave_radiation
     type (BGCVectorDouble) :: air_temperature
     type (BGCVectorDouble) :: vapor_pressure_air
     type (BGCVectorDouble) :: wind_speed
     type (BGCVectorDouble) :: precipitation
     type (BGCVectorDouble) :: precipitation_snow
     type (BGCVectorDouble) :: elevation
     type (BGCVectorDouble) :: aspect
     type (BGCVectorDouble) :: slope
     real (c_double) :: atm_n2
     real (c_double) :: atm_o2
     real (c_double) :: atm_co2
     real (c_double) :: atm_ch4
     real (c_double) :: atm_n2o
     real (c_double) :: atm_h2
     real (c_double) :: atm_nh3
     real (c_double) :: heat_capacity
     real (c_double) :: field_capacity
     real (c_double) :: wilting_point
  end type BGCProperties

  type, public, bind(c) :: BGCAuxiliaryData
     type (BGCVectorInt) :: aux_ints
     type (BGCVectorDouble) :: aux_doubles
  end type BGCAuxiliaryData

end module BGCContainers_module
