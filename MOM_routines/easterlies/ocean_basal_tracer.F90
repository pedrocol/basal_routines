module ocean_basal_tracer_mod
!</CONTACT>
!
!<CONTACT EMAIL="paul.spence@gmail.com"> Paul Spence
!</CONTACT>
!
!<CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> S. M. Griffies 
!</CONTACT>
!
!<OVERVIEW>
! Add basal fresh water flux to ocean.
!</OVERVIEW>
!
!<DESCRIPTION>
! This module applies a 3D fresh water flux, presumably from basal 
! melting (but the source of the FW is arbitrary). The code is based
! on ocean_rivermix.F90. The basal meltwater input file is from 
! Merino et al., 2018: https://doi.org/10.1016/j.ocemod.2017.12.006
! It should apply to both salt and temp tracers. The water
! can occur at any location and with any distribution in the domain.  
!
! The user is responsible for providing (and registering) the data on
! the model grid. Check the land/sea mask in your input file.
!</DESCRIPTION>
!
!<NAMELIST NAME="ocean_basal_tracer_nml">
!  <DATA NAME="use_this_module" TYPE="logical">
!  For using this module.  Default use_this_module=.false.
!  </DATA> 
!  <DATA NAME="test_nml" TYPE="logical">
!  Testing input.nml vars for this moduile.  Default test_nml=.false.
!  Still need to test this.
!  </DATA> 
!</NAMELIST>
!
use diag_manager_mod,         only: register_diag_field
use fms_mod,                  only: write_version_number, open_namelist_file, close_file
use fms_mod,                  only: file_exist
use fms_mod,                  only: open_namelist_file, check_nml_error, close_file
use fms_mod,                  only: read_data, lowercase, FATAL, WARNING, stdout, stdlog
use mpp_mod,                  only: input_nml_file, mpp_sum, mpp_error, mpp_max
use time_interp_external_mod, only: init_external_field, time_interp_external
use time_manager_mod,         only: time_type, set_date, get_time
use time_manager_mod,         only: operator( + ), operator( - ), operator( // )
use time_manager_mod,         only: operator( > ), operator( == ), operator( <= )

use ocean_domains_mod,        only: get_local_indices
use ocean_parameters_mod,     only: missing_value 
use ocean_types_mod,          only: ocean_domain_type, ocean_grid_type, ocean_thickness_type
use ocean_types_mod,          only: ocean_prog_tracer_type, ocean_options_type, ocean_time_type
use ocean_types_mod,          only: ocean_density_type
use ocean_workspace_mod,      only: wrk1, wrk2, wrk3, wrk4, wrk5, wrk1_2d
use ocean_util_mod,           only: diagnose_3d
use constants_mod,            only: epsln
use axis_utils_mod,           only: frac_index, nearest_index
use ocean_types_mod,          only: ocean_time_type, ocean_time_steps_type, ocean_options_type
use ocean_util_mod,           only: diagnose_2d, diagnose_3d, diagnose_sum
use ocean_tracer_util_mod,    only: diagnose_3d_rho
use ocean_tempsalt_mod,       only: contemp_from_pottemp

implicit none

private

private watermass_diag_init_ba
private watermass_diag_basal

#include <ocean_memory.h>

type ocean_basal_type
   integer :: id                                             ! time_interp_external index
   character(len=32) :: name                                 ! tracer name corresponding to basal
   real, dimension(:,:,:), pointer :: damp_coeff   => NULL() ! 3d inverse damping rate (tracer units/ sec)
   real, dimension(:,:)  , pointer :: damp_coeff2d => NULL() ! 2d inverse damping rate (tracer units/ sec)
   real, dimension(:,:)  , pointer :: field        => NULL() ! runoff
end type ocean_basal_type


type(ocean_basal_type), allocatable, dimension(:) :: Basal
type(ocean_domain_type), pointer :: Dom => NULL()
type(ocean_grid_type),   pointer :: Grd => NULL()

public ocean_basal_tracer_init
public basal_tracer_source

character(len=126)  :: version = '$Id: ocean_basal_tracer.F90,v 20.0 2013/12/14 00:16:24 fms Exp $'
character (len=128) :: tagname = '$Name: tikal $'

! for diagnostics 
logical :: used
integer, dimension(:), allocatable :: id_basal_tend
integer :: id_basal_fwflx        =-1
integer :: id_basal_fwflx2d        =-1

integer, dimension(:), allocatable :: id_basalfw
integer, dimension(:), allocatable :: id_basalmix_on_nrho

integer :: id_neut_rho_basalmix    =-1
integer :: id_wdian_rho_basalmix   =-1
integer :: id_tform_rho_basalmix   =-1

integer :: id_neut_rho_basalmix_on_nrho   =-1
integer :: id_wdian_rho_basalmix_on_nrho   =-1
integer :: id_tform_rho_basalmix_on_nrho   =-1

integer :: id_neut_rho_pbl_ba_kn          =-1
integer :: id_neut_rho_pbl_ba_kn_on_nrho  =-1
integer :: id_wdian_rho_pbl_ba_kn         =-1
integer :: id_wdian_rho_pbl_ba_kn_on_nrho =-1
integer :: id_tform_rho_pbl_ba_kn         =-1
integer :: id_tform_rho_pbl_ba_kn_on_nrho =-1

integer :: id_neut_temp_pbl_ba_kn          =-1
integer :: id_neut_temp_pbl_ba_kn_on_nrho  =-1
integer :: id_wdian_temp_pbl_ba_kn         =-1
integer :: id_wdian_temp_pbl_ba_kn_on_nrho =-1
integer :: id_tform_temp_pbl_ba_kn         =-1
integer :: id_tform_temp_pbl_ba_kn_on_nrho =-1

integer :: id_neut_salt_pbl_ba_kn          =-1
integer :: id_neut_salt_pbl_ba_kn_on_nrho  =-1
integer :: id_wdian_salt_pbl_ba_kn         =-1
integer :: id_wdian_salt_pbl_ba_kn_on_nrho =-1
integer :: id_tform_salt_pbl_ba_kn         =-1
integer :: id_tform_salt_pbl_ba_kn_on_nrho =-1

integer :: id_neut_rho_pbl_ba_pr          =-1
integer :: id_neut_rho_pbl_ba_pr_on_nrho  =-1
integer :: id_wdian_rho_pbl_ba_pr         =-1
integer :: id_wdian_rho_pbl_ba_pr_on_nrho =-1
integer :: id_tform_rho_pbl_ba_pr         =-1
integer :: id_tform_rho_pbl_ba_pr_on_nrho =-1

integer :: id_neut_temp_pbl_ba_pr          =-1
integer :: id_neut_temp_pbl_ba_pr_on_nrho  =-1
integer :: id_wdian_temp_pbl_ba_pr         =-1
integer :: id_wdian_temp_pbl_ba_pr_on_nrho =-1
integer :: id_tform_temp_pbl_ba_pr         =-1
integer :: id_tform_temp_pbl_ba_pr_on_nrho =-1

integer :: id_neut_salt_pbl_ba_pr          =-1
integer :: id_neut_salt_pbl_ba_pr_on_nrho  =-1
integer :: id_wdian_salt_pbl_ba_pr         =-1
integer :: id_wdian_salt_pbl_ba_pr_on_nrho =-1
integer :: id_tform_salt_pbl_ba_pr         =-1
integer :: id_tform_salt_pbl_ba_pr_on_nrho =-1

integer :: id_eta_tend_basalmix        =-1
integer :: id_eta_tend_basalmix_glob   =-1

integer :: id_tbasal=-1

integer :: num_prog_tracers      = 0
logical :: module_is_initialized = .FALSE.
logical :: damp_coeff_3d         = .false. 
logical :: use_basal_module       = .false. 
logical :: test_nml              = .false. 

! internally set for computing watermass diagnostics
logical :: compute_watermass_diag_ba = .false.

! for global area normalization
real    :: cellarea_r

! Adaptive restoring
real    :: athresh               = 0.5
real    :: npower                = 1.0
real    :: lambda                = 0.0083
real    :: taumin                = 720
logical :: use_adaptive_restore  = .false.
logical :: use_basal_after_init = .false.
logical :: use_normalising       = .false.
logical :: use_hard_thump        = .false.
logical :: deflate               = .false.
real    :: deflate_fraction      = 0.6
integer :: secs_to_restore       = 0
integer :: days_to_restore       = 1

logical :: limit_temp            = .true.
real    :: limit_temp_min        = -1.8
real    :: limit_temp_restore    = 10800.

logical :: limit_salt            = .false.
real    :: limit_salt_min        = 0.01
real    :: limit_salt_restore    = 3600.

integer :: secs_restore
integer :: initial_day
integer :: initial_secs
real, allocatable :: sdiffo(:)
logical :: debug_all_in_top_cell = .false.


namelist /ocean_basal_tracer_nml/ use_basal_module, test_nml, damp_coeff_3d

contains


!#######################################################################
! <SUBROUTINE NAME="ocean_basal_tracer_init">
!
! <DESCRIPTION>
! This subroutine is intended to be used to initialize addition of basal melt
! water fluxes from an input file. 
! Everything in this subroutine is a user prototype, and should be replacable.
! </DESCRIPTION>
!
subroutine ocean_basal_tracer_init(Grid, Domain, Time, T_prog, dtime, Ocean_options, Dens)

  type(ocean_grid_type),        intent(in), target :: Grid
  type(ocean_domain_type),      intent(in), target :: Domain
  type(ocean_time_type),        intent(in)         :: Time
  type(ocean_prog_tracer_type), intent(in)         :: T_prog(:)
  real,                         intent(in)         :: dtime
  type(ocean_options_type),     intent(inout)      :: Ocean_options
  type(ocean_density_type),     intent(in)         :: Dens


  integer :: i, j, k, n
  integer :: ioun, io_status, ierr
  integer :: index_temp
  integer :: secs, days
  real    :: dtimer

  character(len=128) :: name

  integer :: stdoutunit,stdlogunit 
  stdoutunit=stdout();stdlogunit=stdlog() 

  if ( module_is_initialized ) then 
    call mpp_error(FATAL, &
    '==>Error in ocean_basal_tracer_mod (ocean_basal_tracer_init): module already initialized')
  endif 

  module_is_initialized = .TRUE.

  num_prog_tracers = 3

  allocate ( Basal(num_prog_tracers) )

  call write_version_number(version, tagname)

  ! provide for namelist over-ride of default values
#ifdef INTERNAL_FILE_NML
  read(input_nml_file, nml=ocean_basal_tracer_nml, iostat=io_status)
  ierr = check_nml_error(io_status, 'ocean_basal_tracer_nml')
#else
  ioun =  open_namelist_file()
  read(ioun, ocean_basal_tracer_nml, IOSTAT=io_status)
  ierr = check_nml_error(io_status, 'ocean_basal_tracer_nml')
  call close_file(ioun)
#endif
  write(stdlogunit, ocean_basal_tracer_nml)
  write(stdoutunit, '(/)')
  write(stdoutunit, ocean_basal_tracer_nml)

  Dom => Domain
  Grd => Grid

#ifndef MOM_STATIC_ARRAYS    
  call get_local_indices(Domain, isd, ied, jsd, jed, isc, iec, jsc, jec)
  nk = Grid%nk
  !bottom level Grid%kmt(i,j)
#endif

  do n=1,num_prog_tracers
     Basal(n)%id = -1
  enddo

  dtimer = 1.0/dtime

  if(use_basal_module) then
      write(stdoutunit,*)'==>Note from ocean_basal_tracer_mod: Using this module.'
      Ocean_options%ocean_basal_tracer= 'Used ocean tracer basal.'
  else
      write(stdoutunit,*)'==>Note from ocean_basal_tracer_mod: NOT using ocean tracer basal.'
      Ocean_options%ocean_basal_tracer= 'Did NOT use ocean tracer basal.'
      return
  endif

  do n = 1, num_prog_tracers
    ! Read basal fw flux  data
    name = 'INPUT/basal_fw.nc'
    if ( n == 1 ) then
       Basal(n)%id = init_external_field(name,'meltingFlux',domain=Domain%domain2d)
       if (Basal(n)%id < 1) then
         call mpp_error(FATAL,&
         '==>Error: in ocean_basal_tracer_mod:  basal fw values are not specified')
       endif
       Basal(n)%name = 'basal_fw' 
    elseif ( n == 2 ) then
       Basal(n)%id = init_external_field(name,'sodepmin_isf',domain=Domain%domain2d)
       if (Basal(n)%id < 1) then
         call mpp_error(FATAL,&
         '==>Error: in ocean_basal_tracer_mod:  basal front values are not specified')
       endif
       Basal(n)%name = 'depmin' 
    elseif ( n == 3 ) then
       Basal(n)%id = init_external_field(name,'sodepmax_isf',domain=Domain%domain2d)
       if (Basal(n)%id < 1) then
         call mpp_error(FATAL,&
         '==>Error: in ocean_basal_tracer_mod:  basal GL values are not specified')
       endif
       Basal(n)%name = 'depmax' 
    endif
    write(stdoutunit,*) '==> Using basal freshwater flux data specified from file '//trim(name)
  enddo

  ! register diagnostic outputs
  id_tbasal = register_diag_field ('ocean_model','temp_basal',Grd%tracer_axes(1:2),          &
              Time%model_time, 'temperature of basal melt runoff entering the ocean',        &
              'degC' , missing_value=missing_value,range=(/-1.e10,1.e10/))
  id_basal_fwflx = register_diag_field ('ocean_model','basal_fwflx', Grd%tracer_axes(1:3),   &
              Time%model_time, 'mass flux of liquid basal meltwater entering ocean ',        &
              '(kg/m^3)*(m/sec)', missing_value=missing_value,range=(/-1e6,1e6/),            &
              standard_name='water_flux_into_sea_water_from_basal_melting')
  id_basal_fwflx2d = register_diag_field ('ocean_model','basal_fwflx2d', Grd%tracer_axes(1:2),   &
              Time%model_time, '2d mass flux of liquid basal meltwater entering ocean ',        &
              '(kg/m^3)*(m/sec)', missing_value=missing_value,range=(/-1e6,1e6/),            &
              standard_name='water_flux_into_sea_water_from_basal_melting')

  
  allocate (id_basal_tend(num_prog_tracers))
  id_basal_tend = -1

  do n=1,num_prog_tracers
     id_basal_tend(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_basal_fw_tend',&
               Grd%tracer_axes(1:3), Time%model_time, 'rho*dzt*tendency due to basal freshwater flux',          &
               trim(T_prog(n)%flux_units), missing_value=missing_value, range=(/-1.e10,1.e10/))
  enddo

  ! register for diag_manager
  num_prog_tracers = size(T_prog(:))
  allocate (id_basalfw(num_prog_tracers))
  allocate (id_basalmix_on_nrho(num_prog_tracers))
  id_basalfw   = -1
  id_basalmix_on_nrho   = -1
  do n=1,num_prog_tracers
     if(T_prog(n)%name == 'temp') then
        id_basalfw(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_basalmix', &
                         Grd%tracer_axes(1:3), Time%model_time,                                 &
                         'cp*basalmix*rho_dzt*temp', 'Watt/m^2',                                &
                         missing_value=missing_value, range=(/-1.e10,1.e10/))
        id_basalmix_on_nrho(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_basalmix_on_nrho', &
                         Dens%neutralrho_axes(1:3), Time%model_time,                                 &
                         'cp*basalmix*rho_dzt*temp binned to neutral density', 'Watt/m^2',       &
                         missing_value=missing_value, range=(/-1.e20,1.e20/))
     else
        id_basalfw(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_basalmix', &
                         Grd%tracer_axes(1:3), Time%model_time,                                 &
                         'basalmix*rho_dzt*tracer for '//trim(T_prog(n)%name),                  &
                         trim(T_prog(n)%flux_units),                                            &
                         missing_value=missing_value, range=(/-1.e10,1.e10/))
        id_basalmix_on_nrho(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_basalmix_on_nrho', &
                         Dens%neutralrho_axes(1:3), Time%model_time,                                 &
                         'basalmix*rho_dzt*tracer for '//trim(T_prog(n)%name)//' binned to neutral density',&
                         trim(T_prog(n)%flux_units),                                            &
                         missing_value=missing_value, range=(/-1.e10,1.e10/))
    endif
  enddo


  
  !try to write out basal fwflx data to output for test
  ! get sponge value for current time
  wrk1=0.0
  call time_interp_external(Basal(1)%id, Time%model_time, wrk1)
  call time_interp_external(Basal(2)%id, Time%model_time, wrk1)
  call time_interp_external(Basal(3)%id, Time%model_time, wrk1)
  !call diagnose_3d(Time, Grd, id_basal_fwflx,wrk1(:,:,:))

  call watermass_diag_init_ba(Time, Dens)


end subroutine ocean_basal_tracer_init
! </SUBROUTINE> NAME="ocean_basal_tracer_init"


!#######################################################################
! <SUBROUTINE NAME="basal_tracer_source">
!
! <DESCRIPTION>
! This subroutine calculates thickness weighted and density weighted
! time tendencies of tracers due to damping by basal.
! </DESCRIPTION>
!
subroutine basal_tracer_source(Time, Time_steps, Thickness, Dens, T_prog, basal, diff_cbt, index_temp, &
                               index_salt, basal3d)

  type(ocean_time_type),          intent(in)      :: Time
  type(ocean_time_steps_type),    intent(in)      :: Time_steps
  type(ocean_thickness_type),     intent(in)      :: Thickness
  type(ocean_density_type),       intent(in)      :: Dens
  type(ocean_prog_tracer_type),   intent(inout)   :: T_prog(:)
  real, dimension(isd:,jsd:),     intent(inout)   :: basal
  real, dimension(isd:,jsd:,:,:), intent(inout)   :: diff_cbt
  real, dimension(isd:,jsd:,:)  , intent(inout)   :: basal3d
  integer,                        intent(in)      :: index_temp
  integer,                        intent(in)      :: index_salt
  integer :: param_choice,n
  integer, allocatable, dimension(:,:) :: misfkt,misfkb ! Top and bottom input depths

  if(.not. use_basal_module) return

  allocate ( misfkt(isd:ied,jsd:jed), misfkb(isd:ied,jsd:jed) )
  misfkt(:,:) = 0
  misfkb(:,:) = 0
  !basal(:,:) = 0.0

  param_choice = 1
  cellarea_r = 1.0/(epsln + Grd%tcellsurf)
  num_prog_tracers = size(T_prog(:))

  IF ( param_choice == 1 ) THEN
    CALL basal_tracer_source_1(Time, Time_steps, Thickness, T_prog(1:num_prog_tracers), basal, diff_cbt, index_temp, &
                               index_salt, misfkt, misfkb, Dens, basal3d)
  ELSEIF ( param_choice == 3 ) THEN
    !Do nothing for the moment
  ENDIF

  call watermass_diag_basal(Time, Dens, T_prog, basal, &
  T_prog(index_temp)%wrk1(:,:,:),T_prog(index_salt)%wrk1(:,:,:), misfkt, misfkb)

  deallocate ( misfkt, misfkb )

end subroutine basal_tracer_source
! </SUBROUTINE> NAME="basal_tracer_source"

!#######################################################################
! <SUBROUTINE NAME="basal_tracer_source_1">
!
! <DESCRIPTION>
! This subroutine calculates thickness weighted and density weighted
! time tendencies of tracers due to damping by basal.
! </DESCRIPTION>
!
subroutine basal_tracer_source_1(Time, Time_steps, Thickness, T_prog, basal_i,diff_cbt,index_temp, index_salt, &
                                 misfkt, misfkb, Dens, basal3d)
  ! Case: specified fwf and heat flux forcing beneath the ice shelf
  type(ocean_time_type),          intent(in)    :: Time
  type(ocean_time_steps_type),    intent(in)    :: Time_steps
  type(ocean_thickness_type),     intent(in)    :: Thickness
  type(ocean_prog_tracer_type),   intent(inout) :: T_prog(:)
  real, dimension(isd:,jsd:),     intent(inout)    :: basal_i
  integer,                        intent(in)     :: index_temp
  integer,                        intent(in)     :: index_salt
  real, dimension(isd:,jsd:,:,:), intent(inout)  :: diff_cbt
  integer, dimension(isd:,jsd:),  intent(inout)  :: misfkt,misfkb ! Top and bottom input depths
  type(ocean_density_type),       intent(in)     :: Dens
  real, dimension(isd:,jsd:,:),   intent(inout)  :: basal3d
  real    :: dtime
  integer :: taum1, tau, taup1
  integer :: i, j, k, n, nz, nn
  real, allocatable, dimension(:,:)    :: misfzt,misfzb ! Top and bottom input depths
  real, allocatable, dimension(:,:) :: fwfisf        ! fresh water flux from the isf (fwfisf <0 mean melting) [Kg/m2/s]
  real, allocatable, dimension(:,:) :: qisf          ! heat flux
  real                              :: rau0          ! volumic mass of reference     [kg/m3]
  real                              :: rcp           ! heat capacity     [J/K]
  real                              :: rau0_rcp, r1_rau0_rcp, r1_rau0, soce, rLfusisf
  real, allocatable, dimension(:,:) :: stbl          ! Salinity top boundary layer
  real, allocatable, dimension(:,:) :: zt_frz        ! Freezing point temperature
  real, allocatable, dimension(:,:) :: rhisf_tbl, r1_hisf_tbl   !: thickness of tbl  [m]
  real, allocatable, dimension(:,:,:) :: risf_tsc_b , risf_tsc ! before and now T & S isf contents [K.m/s & PSU.m/s]
  real, allocatable, dimension(:,:,:) :: zt_frz3d        ! Freezing point temperature
  real, allocatable, dimension(:,:,:,:) :: risf_tsc_3d_b , risf_tsc_3d ! before and now T & S isf contents [K.m/s & PSU.m/s]
  integer ::   ikt, ikb, import_file, method, firstlev, gade_line     ! local integers
  integer :: nvars
  !rivermix
  real    :: depth, thkocean
  real    :: delta(nk), delta_rho_tocean(nk), delta_rho0_triver(nk)
  real    :: zextra, zinsert, tracerextra, tracernew(nk)
  real    :: tracer_input, tfreezing, tbasal, tbasal_sum
  real    :: maxinsertiondepth,mininsertiondepth
  logical :: river_diffuse_temp=.true.    ! to enhance diffusivity of temp at river mouths over river_thickness column
  logical :: river_diffuse_salt=.true.    ! to enhance diffusivity of salt at river mouths over river_thickness column
  integer :: stdoutunit,stdlogunit
  character(len=128) :: name
  real    :: sal, press
  real    :: c_p = 3974.0
  real    :: L_f = 3.34*10**5
  stdoutunit=stdout();stdlogunit=stdlog()
!#######################################################################

  allocate ( misfzt(isd:ied,jsd:jed), misfzb(isd:ied,jsd:jed) )
  allocate ( fwfisf(isd:ied,jsd:jed),   qisf(isd:ied,jsd:jed) )
  allocate (   stbl(isd:ied,jsd:jed), zt_frz(isd:ied,jsd:jed) )
  allocate (   zt_frz3d(isd:ied,jsd:jed,nk) )
  allocate ( rhisf_tbl(isd:ied,jsd:jed), r1_hisf_tbl(isd:ied,jsd:jed) )
  allocate ( risf_tsc_b(isd:ied,jsd:jed,2), risf_tsc(isd:ied,jsd:jed,2) ) !1=temp 2=sal
  allocate ( risf_tsc_3d_b(isd:ied,jsd:jed,nk,2), risf_tsc_3d(isd:ied,jsd:jed,nk,2) ) !1=temp 2=sal

  taum1 = Time%taum1
  tau   = Time%tau
  taup1 = Time%taup1

  wrk1  = 0.0 !Asumes T=0 and S=0
  dtime = Time_steps%dtts
  num_prog_tracers = size(T_prog(:))
  nvars = 3
  misfzt = 0.0
  misfzb = 0.0
  sal = 0.0
  press = 0.0
  tfreezing = 0.0
  tbasal = 0.0
  tbasal_sum = 0.0

  do n=1,num_prog_tracers
    T_prog(n)%wrk1(:,:,:) = 0.0
  enddo
  delta             = 0.0
  delta_rho_tocean  = 0.0
  delta_rho0_triver = 0.0

  import_file = 2

  if ( import_file == 2) then
    if ( Basal(1)%id > 0 ) then
       ! get basal value for current time
       call time_interp_external(Basal(1)%id, Time%model_time, wrk1)
       do j=jsd,jed
         do i=isd,ied
            fwfisf(i,j)=wrk1(i,j,1)
         enddo
       enddo
    endif
    if (Basal(2)%id > 0) then
       wrk1  = 0.0
       ! get basal value for current time
       call time_interp_external(Basal(2)%id, Time%model_time, wrk1)
       do j=jsd,jed
         do i=isd,ied
            misfzt(i,j)=wrk1(i,j,1)
         enddo
       enddo
    endif
    if (Basal(3)%id > 0) then
       wrk1  = 0.0 
       ! get basal value for current time
       call time_interp_external(Basal(3)%id, Time%model_time, wrk1)
       do j=jsd,jed
         do i=isd,ied
            misfzb(i,j)=wrk1(i,j,1)
         enddo
       enddo
    endif
  else
    misfzt(:,:)=100
    misfzb(:,:)=400
    fwfisf(:,:) = 0.0
    fwfisf(:,:) = basal_i(:,:)
  endif

  !For this first test we use the original river values
  fwfisf(:,:) = 0.0
  fwfisf(:,:) = basal_i(:,:)

  wrk2  = 0.0

  do j=jsc,jec
     do i=isc,iec

        tracernew(:) = 0.0

        if (fwfisf(i,j) > 0.0 .and. Grd%kmt(i,j) > 0) then

           if ( misfzt(i,j) == 0 .or. misfzb(i,j) == 0) then
              misfzt(i,j)=200
              misfzb(i,j)=400
           endif
           if ( misfzb(i,j) < misfzt(i,j) ) then
              misfzt(i,j)=200
              misfzb(i,j)=400
           endif
           
           !Test like BG03
           !misfzb(:,:) = misfzt(:,:) + 50

        ! the array "river" contains the volume rate (m/s) or mass
        ! rate (kg/m2/sec) of fluid with tracer 
        ! that is to be distributed in the vertical.
           mininsertiondepth = misfzt(i,j)
           depth       = min(Grd%ht(i,j),mininsertiondepth) !min between bathy and value
           misfkt(i,j) = min(Grd%kmt(i,j),floor(frac_index(depth,Grd%zw))) !min(mbathy,
           misfkt(i,j) = max(1,misfkt(i,j))
           !maxinsertiondepth = Grd%zt(misfkb(i,j))
           maxinsertiondepth = misfzb(i,j)
           depth       = min(Grd%ht(i,j),maxinsertiondepth)     ! be sure not to discharge river content into rock, ht = ocean topography
           misfkb(i,j) = min(Grd%kmt(i,j),floor(frac_index(depth,Grd%zw))) ! max number of k-levels into which discharge rivers
           misfkb(i,j) = max(1,misfkb(i,j))                                         ! make sure have at least one cell to discharge into

           ! determine fractional thicknesses of grid cells 
           thkocean = 0.0
           do k=misfkt(i,j),misfkb(i,j)
              thkocean = thkocean + Thickness%rho_dzt(i,j,k,tau)
           enddo
           do k=misfkt(i,j),misfkb(i,j)
              delta(k) = Thickness%rho_dzt(i,j,k,tau)/(epsln+thkocean)
           enddo

           do n=1,num_prog_tracers
              gade_line = 0

              if (gade_line == 0 ) then
                 if ( n == 1 ) then
                    nn = index_salt       
                 elseif ( n == 2 ) then
                    nn = index_temp
                 else
                    nn = n
                 end if

                 tbasal = 0.0
                 tbasal_sum = 0.0

                 if ( trim(T_prog(nn)%name) /= 'temp' ) then !salt and age
                    zextra=0.0
                    do k=misfkb(i,j),misfkt(i,j),-1
                       tracernew(k) = 0.0

                       if (k.eq.misfkb(i,j)) then
                          tracerextra=0.0
                       else
                          tracerextra = tracernew(k+1)
                       endif

                       basal3d(i,j,k) = fwfisf(i,j)*delta(k)
                       zinsert = basal3d(i,j,k)*dtime
                       tracernew(k) = (tracerextra*zextra + T_prog(nn)%field(i,j,k,tau)*Thickness%rho_dzt(i,j,k,tau) + &
                                       tbasal*zinsert) / (zextra+Thickness%rho_dzt(i,j,k,tau)+zinsert)

                       zextra=zextra+zinsert

                       tbasal_sum = tbasal_sum + tbasal*delta(k)
                       wrk2(i,j,k) = tracernew(k) - T_prog(nn)%field(i,j,k,tau) !Delta S
                    enddo

                    T_prog(nn)%tbasal(i,j) = tbasal_sum !average for ocean_diagnostics

                    firstlev = misfkt(i,j)

                 else    !temp
                    do k=misfkt(i,j),misfkb(i,j)
                       tracernew(k) = 0.0
                       tracernew(k) = T_prog(index_temp)%field(i,j,k,tau) + &
                                    wrk2(i,j,k) / T_prog(index_salt)%field(i,j,k,tau) * L_f / c_p

                       zinsert = fwfisf(i,j)*dtime*delta(k)
                       tbasal = ( (tracernew(k) * (Thickness%rho_dzt(i,j,k,tau)+zinsert)) - &
                               T_prog(index_temp)%field(i,j,k,tau)*Thickness%rho_dzt(i,j,k,tau) ) / zinsert 
                       tbasal_sum = tbasal_sum + tbasal*delta(k)
                    enddo
        
                    T_prog(nn)%tbasal(i,j) = tbasal_sum !average for ocean_diagnostics

                 endif !gade


                 k=1
                 if ( misfkt(i,j) > 1 ) tracernew(k) = Thickness%rho_dzt(i,j,k,tau)
                 T_prog(nn)%wrk1(i,j,k) = (tracernew(k)*(Thickness%rho_dzt(i,j,k,tau)+fwfisf(i,j)*dtime) -&
                                           T_prog(nn)%field(i,j,k,tau)*Thickness%rho_dzt(i,j,k,tau))/dtime
                 do k=misfkt(i,j),misfkb(i,j)
                    T_prog(nn)%wrk1(i,j,k) = Thickness%rho_dzt(i,j,k,tau)*(tracernew(k) - T_prog(nn)%field(i,j,k,tau))/dtime !Tendency
                 enddo

                 !Update tendency
                 do k=misfkt(i,j),misfkb(i,j)
                    T_prog(nn)%th_tendency(i,j,k) = T_prog(nn)%th_tendency(i,j,k) + T_prog(nn)%wrk1(i,j,k)
                 enddo
              elseif (gade_line == 1 ) then

                 zextra=0.0

                 tbasal_sum = 0.0

                 if ( trim(T_prog(n)%name) == 'temp' ) then
                    sal = T_prog(index_salt)%field(i,j,k,tau)
                    !sal = 0.0
                    press = Dens%pressure_at_depth(i,j,k)
                    call frz_preteos10(tfreezing, press, sal)
                    !tbasal = tfreezing
                    tfreezing = contemp_from_pottemp(sal,tfreezing)
                    tbasal = min(tfreezing,T_prog(index_temp)%field(i,j,k,tau))
                 endif
                 if ( trim(T_prog(n)%name) == 'salt' ) tbasal = 0.0
                 !tbasal = T_prog(n)%triver(i,j)

                 method = 1

                 if ( method == 0 ) then
 
                    do k=misfkb(i,j),misfkt(i,j),-1
  
                       tracernew(k) = 0.0

                       if (k.eq.misfkb(i,j)) then
                          tracerextra=0.0
                       else
                          tracerextra = tracernew(k+1)
                       endif

                       zinsert = fwfisf(i,j)*dtime*delta(k)
                       tracernew(k) = (tracerextra*zextra + T_prog(n)%field(i,j,k,tau)*Thickness%rho_dzt(i,j,k,tau) + &
                                      tbasal*zinsert) / (zextra+Thickness%rho_dzt(i,j,k,tau)+zinsert)
                       basal3d(i,j,k) = ( zextra + zinsert ) / ( dtime*delta(k) ) !Equivalente basal flux
                       zextra=zextra+zinsert

                       tbasal_sum = tbasal_sum + tbasal*delta(k)
                    enddo

                    k=misfkt(i,j) !Treatment at the first level
                    T_prog(n)%wrk1(i,j,k) = (tracernew(k)*(Thickness%rho_dzt(i,j,k,tau)+fwfisf(i,j)*dtime) -&
                                            T_prog(n)%field(i,j,k,tau)*Thickness%rho_dzt(i,j,k,tau))/dtime

                    firstlev = misfkt(i,j)+1
                 elseif ( method == 1 ) then
                    do k=misfkt(i,j),misfkb(i,j)
                       basal3d(i,j,k) = fwfisf(i,j)*delta(k)
                       zinsert = fwfisf(i,j)*dtime*delta(k)
                       tracernew(k) = (T_prog(n)%field(i,j,k,tau)*Thickness%rho_dzt(i,j,k,tau) + tbasal*zinsert) / &
                                      (Thickness%rho_dzt(i,j,k,tau)+zinsert)
                       tbasal_sum = tbasal_sum + tbasal*delta(k)
                    enddo
                    firstlev = misfkt(i,j)                 
                 endif !method

                 T_prog(n)%tbasal(i,j) = tbasal_sum !average for ocean_diagnostics
  
                 do k=firstlev,misfkb(i,j)
                    T_prog(n)%wrk1(i,j,k) = Thickness%rho_dzt(i,j,k,tau)*(tracernew(k) - T_prog(n)%field(i,j,k,tau))/dtime !Tendency 
                 enddo

                 if(debug_all_in_top_cell) then
                     k=1
                     T_prog(n)%wrk1(i,j,:) = 0.0
                     T_prog(n)%wrk1(i,j,k) = Grd%tmask(i,j,k)*fwfisf(i,j)*T_prog(n)%triver(i,j)
                 endif

                 do k=misfkt(i,j),misfkb(i,j)
                    T_prog(n)%th_tendency(i,j,k) = T_prog(n)%th_tendency(i,j,k) + T_prog(n)%wrk1(i,j,k)
                 enddo
              endif !gade line    
           enddo !n
        endif ! fwfisf > 0
     enddo !i
  enddo !j

  !!Diags
  if (id_tbasal > 0) then !basal temp
     wrk1_2d=0.0
     do j=jsc,jec
        do i=isc,iec
           if ( fwfisf(i,j) /= 0.0 ) then
              wrk1_2d(i,j) = T_prog(index_temp)%tbasal(i,j)
           endif
        enddo
     enddo
     call diagnose_2d(Time, Grd, id_tbasal, wrk1_2d(:,:))
  endif

  if (id_basal_fwflx > 0) then !basal flux
     ! basal entering the ocean (kg/m^3)*(m/s)
     call diagnose_3d(Time, Grd, id_basal_fwflx, basal3d(:,:,:))
  endif

  if (id_basal_fwflx2d> 0) then !basal flux
     ! basal entering the ocean (kg/m^3)*(m/s)
     call diagnose_2d(Time, Grd, id_basal_fwflx2d, fwfisf(:,:))
  endif


  do n=1,num_prog_tracers
     if(id_basalfw(n) > 0) then !temp/salt flux
        call diagnose_3d(Time, Grd, id_basalfw(n), T_prog(n)%wrk1(:,:,:)*T_prog(n)%conversion)
     endif
     if(id_basalmix_on_nrho(n) > 0) then
        call diagnose_3d_rho(Time, Dens, id_basalmix_on_nrho(n), T_prog(n)%wrk1*T_prog(n)%conversion)
     endif
  enddo

  deallocate ( misfzt, misfzb )
  deallocate ( fwfisf,   qisf )
  deallocate (   stbl, zt_frz )
  deallocate ( zt_frz3d )
  deallocate ( rhisf_tbl, r1_hisf_tbl )
  deallocate ( risf_tsc_b, risf_tsc )
  deallocate ( risf_tsc_3d_b, risf_tsc_3d )

!  if(river_diffuse_temp) then
!     call river_kappa(Time, Thickness, T_prog(index_temp), diff_cbt(isd:ied,jsd:jed,:,1))
!  endif
!  if(river_diffuse_salt) then
!     call river_kappa(Time, Thickness, T_prog(index_salt), diff_cbt(isd:ied,jsd:jed,:,2))
!  endif


end subroutine basal_tracer_source_1
! </SUBROUTINE> NAME="basal_tracer_source_1"

!#######################################################################
! <SUBROUTINE NAME="frz_preteos10">
!
! <DESCRIPTION>
! To use the accurate freezing point temperature of seawater,
! which is a nonlinear function of salinity and pressure.
! air_saturated_water = true
! </DESCRIPTION>
!
subroutine frz_preteos10(tfreezing, press, sal)
   real, intent(inout)    :: tfreezing
   real, intent(in)       :: press
   real, intent(in)       :: sal
   real     :: tf_num, tf_den, tfreeze, sqrts
   real     :: tfreeze_check
   ! coefficients in freezing temperature with preTEOS10
   real :: a0, a1, a2, a3, a4, a5, a6
   real :: b0, b1, b2, b3
   real :: c1, c2

   ! coefficients in the preTEOS10 freezing point of seawater
   ! assume prognostic temperature variable is potential temp
   ! note that TEOS10 coefficients set as parameters above. 

   a0 =  2.5180516744541290e-03
   a1 = -5.8545863698926184e-02
   a2 =  2.2979985780124325e-03
   a3 = -3.0086338218235500e-04
   a4 = -7.0023530029351803e-04
   a5 =  8.4149607219833806e-09
   a6 =  1.1845857563107403e-11

   b0 =  1.0000000000000000e+00
   b1 = -3.8493266309172074e-05
   b2 =  9.1686537446749641e-10
   b3 =  1.3632481944285909e-06

   c1 = -2.5180516744541290e-03
   c2 =  1.428571428571429e-05
   tfreeze_check = -2.076426227617581
   sqrts   = sqrt(sal)
   tf_num  = a0 + sal*(a1 + sqrts*(a2 + sqrts*a3)) + press*(a4 + press*(a5 + sal*a6))
   tf_den  = b0 + press*(b1 + press*b2) + sal*sal*sqrts*b3
   tfreeze = tf_num/tf_den + (c1 + sal*c2)
   tfreezing = tfreeze

end subroutine frz_preteos10
! </SUBROUTINE> NAME="frz_preteos10"


!#######################################################################
! <SUBROUTINE NAME="watermass_diag_init_ba">
!
! <DESCRIPTION>
! Initialization of watermass diagnostic output files.
! </DESCRIPTION>
!
subroutine watermass_diag_init_ba(Time, Dens)

  type(ocean_time_type),    intent(in) :: Time
  type(ocean_density_type), intent(in) :: Dens
  integer :: stdoutunit
  stdoutunit=stdout()
  compute_watermass_diag_ba = .false.
  ! runoff plus calving mixing
  id_neut_rho_basalmix = register_diag_field ('ocean_model', 'neut_rho_basalmix',&
     Grd%tracer_axes(1:3), Time%model_time,                                      &
     'update of locally ref potrho from basalmix scheme',                        &
     '(kg/m^3)/sec', missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_neut_rho_basalmix > 0) compute_watermass_diag_ba = .true.

  id_neut_rho_basalmix_on_nrho = register_diag_field ('ocean_model',                    &
    'neut_rho_basalmix_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,            &
    'update of locally ref potrho from basalmix scheme as binned to neutral rho layers',&
    '(kg/m^3)/sec', missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_neut_rho_basalmix_on_nrho > 0) compute_watermass_diag_ba = .true.

  id_wdian_rho_basalmix = register_diag_field ('ocean_model', 'wdian_rho_basalmix',&
     Grd%tracer_axes(1:3), Time%model_time,                                        &
     'dianeutral mass transport due to basalmix scheme',                           &
     'kg/sec', missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_wdian_rho_basalmix > 0) compute_watermass_diag_ba = .true.

  id_wdian_rho_basalmix_on_nrho = register_diag_field ('ocean_model',                  &
    'wdian_rho_basalmix_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,          &
    'dianeutral mass transport due to basalmix scheme as binned to neutral rho layers',&
    'kg/sec', missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_wdian_rho_basalmix_on_nrho > 0) compute_watermass_diag_ba = .true.

  id_tform_rho_basalmix = register_diag_field ('ocean_model', 'tform_rho_basalmix',&
     Grd%tracer_axes(1:3), Time%model_time,                                        &
     'transform due to basalmix scheme pre-binning to neutral rho',                &
     'kg/sec', missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_tform_rho_basalmix > 0) compute_watermass_diag_ba = .true.

  id_tform_rho_basalmix_on_nrho = register_diag_field ('ocean_model',          &
    'tform_rho_basalmix_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,  &
    'watermass transform from basalmix scheme as binned to neutral rho layers',&
    'kg/sec', missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_tform_rho_basalmix_on_nrho > 0) compute_watermass_diag_ba = .true.

  ! process advective form of runoff + calving  
  id_neut_rho_pbl_ba_pr = register_diag_field ('ocean_model',   &
   'neut_rho_pbl_ba_pr', Grd%tracer_axes(1:3), Time%model_time, &
   'process advective-form material time derivative from basal',&
   'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_neut_rho_pbl_ba_pr > 0) compute_watermass_diag_ba=.true.

  id_neut_rho_pbl_ba_pr_on_nrho = register_diag_field ('ocean_model',                        &
   'neut_rho_pbl_ba_pr_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,                 &
   'process advective-form material time derivative from basal as binned to neutral density',&
   'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_neut_rho_pbl_ba_pr_on_nrho > 0) compute_watermass_diag_ba=.true.

  id_wdian_rho_pbl_ba_pr = register_diag_field ('ocean_model',  &
   'wdian_rho_pbl_ba_pr', Grd%tracer_axes(1:3), Time%model_time,&
   'process advective-form dianeutral transport from basal',    &
   'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_wdian_rho_pbl_ba_pr > 0) compute_watermass_diag_ba=.true.

  id_wdian_rho_pbl_ba_pr_on_nrho = register_diag_field ('ocean_model',                    &
    'wdian_rho_pbl_ba_pr_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,            &
    'process advective-form dianeutral transport from basal as binned to neutral density',&
    'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_wdian_rho_pbl_ba_pr_on_nrho > 0) compute_watermass_diag_ba=.true.

  id_tform_rho_pbl_ba_pr = register_diag_field ('ocean_model',                              &
    'tform_rho_pbl_ba_pr', Grd%tracer_axes(1:3), Time%model_time,                           &
    'process advective-form water mass transform from basal pre-binning to neutral density',&
    'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_tform_rho_pbl_ba_pr > 0) compute_watermass_diag_ba=.true.

   id_tform_rho_pbl_ba_pr_on_nrho = register_diag_field ('ocean_model',                &
   'tform_rho_pbl_ba_pr_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,         &
   'process advective-form water mass transform from basal binned to neutral density',&
   'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_tform_rho_pbl_ba_pr_on_nrho > 0) compute_watermass_diag_ba=.true.

  ! process advective form of runoff + calving from temperature pieces
  id_neut_temp_pbl_ba_pr = register_diag_field ('ocean_model',               &
   'neut_temp_pbl_ba_pr', Grd%tracer_axes(1:3), Time%model_time,             &
   'temp related process advective-form material time derivative from basal',&
   'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_neut_temp_pbl_ba_pr > 0) compute_watermass_diag_ba=.true.

  id_neut_temp_pbl_ba_pr_on_nrho = register_diag_field ('ocean_model',                                 &
   'neut_temp_pbl_ba_pr_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,                          &
   'temp related process advective-form material time derivative from basal binned to neutral density',&
   'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_neut_temp_pbl_ba_pr_on_nrho > 0) compute_watermass_diag_ba=.true.

  id_wdian_temp_pbl_ba_pr = register_diag_field ('ocean_model',          &
   'wdian_temp_pbl_ba_pr', Grd%tracer_axes(1:3), Time%model_time,        &
   'temp related process advective-form dianeutral transport from basal',&
   'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_wdian_temp_pbl_ba_pr > 0) compute_watermass_diag_ba=.true.

  id_wdian_temp_pbl_ba_pr_on_nrho = register_diag_field ('ocean_model',                             &
    'wdian_temp_pbl_ba_pr_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,                     &
    'temp related process advective-form dianeutral transport from basal binned to neutral density',&
    'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_wdian_temp_pbl_ba_pr_on_nrho > 0) compute_watermass_diag_ba=.true.

  id_tform_temp_pbl_ba_pr = register_diag_field ('ocean_model',           &
    'tform_temp_pbl_ba_pr', Grd%tracer_axes(1:3), Time%model_time,        &
    'temp related process advective-form water mass transform from basal',&
    'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_tform_temp_pbl_ba_pr > 0) compute_watermass_diag_ba=.true.

  id_tform_temp_pbl_ba_pr_on_nrho = register_diag_field ('ocean_model',                            &
   'tform_temp_pbl_ba_pr_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,                     &
   'temp related process advective-form water mass transform from basal binned to neutral density',&
   'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_tform_temp_pbl_ba_pr_on_nrho > 0) compute_watermass_diag_ba=.true.

  ! process advective form of runoff + calving from salinity pieces
  id_neut_salt_pbl_ba_pr = register_diag_field ('ocean_model',               &
   'neut_salt_pbl_ba_pr', Grd%tracer_axes(1:3), Time%model_time,             &
   'salt related process advective-form material time derivative from river',&
   'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_neut_salt_pbl_ba_pr > 0) compute_watermass_diag_ba=.true.

  id_neut_salt_pbl_ba_pr_on_nrho = register_diag_field ('ocean_model',                                 &
   'neut_salt_pbl_ba_pr_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,                          &
   'salt related process advective-form material time derivative from basal binned to neutral density',&
   'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_neut_salt_pbl_ba_pr_on_nrho > 0) compute_watermass_diag_ba=.true.

  id_wdian_salt_pbl_ba_pr = register_diag_field ('ocean_model',          &
   'wdian_salt_pbl_ba_pr', Grd%tracer_axes(1:3), Time%model_time,        &
   'salt related process advective-form dianeutral transport from basal',&
   'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_wdian_salt_pbl_ba_pr > 0) compute_watermass_diag_ba=.true.

  id_wdian_salt_pbl_ba_pr_on_nrho = register_diag_field ('ocean_model',                             &
    'wdian_salt_pbl_ba_pr_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,                     &
    'salt related process advective-form dianeutral transport from basal binned to neutral density',&
    'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_wdian_salt_pbl_ba_pr_on_nrho > 0) compute_watermass_diag_ba=.true.

  id_tform_salt_pbl_ba_pr = register_diag_field ('ocean_model',           &
    'tform_salt_pbl_ba_pr', Grd%tracer_axes(1:3), Time%model_time,        &
    'salt related process advective-form water mass transform from basal',&
    'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_tform_salt_pbl_ba_pr > 0) compute_watermass_diag_ba=.true.

  id_tform_salt_pbl_ba_pr_on_nrho = register_diag_field ('ocean_model',                            &
   'tform_salt_pbl_ba_pr_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,                     &
   'salt related process advective-form water mass transform from basal binned to neutral density',&
   'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_tform_salt_pbl_ba_pr_on_nrho > 0) compute_watermass_diag_ba=.true.

  ! kinematic advective form of runoff + calving  
  id_neut_rho_pbl_ba_kn = register_diag_field ('ocean_model',     &
   'neut_rho_pbl_ba_kn', Grd%tracer_axes(1:3), Time%model_time,   &
   'kinematic advective-form material time derivative from basal',&
   'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_neut_rho_pbl_ba_kn > 0) compute_watermass_diag_ba=.true.

  id_neut_rho_pbl_ba_kn_on_nrho = register_diag_field ('ocean_model',                           &
    'neut_rho_pbl_ba_kn_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,                   &
    'kinematic advective-form material time derivative from basal as binned to neutral density',&
    'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_neut_rho_pbl_ba_kn_on_nrho > 0) compute_watermass_diag_ba=.true.

  id_wdian_rho_pbl_ba_kn = register_diag_field ('ocean_model',   &
    'wdian_rho_pbl_ba_kn', Grd%tracer_axes(1:3), Time%model_time,&
    'kinematic advective-form dianeutral transport from basal',  &
    'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_wdian_rho_pbl_ba_kn > 0) compute_watermass_diag_ba=.true.

  id_wdian_rho_pbl_ba_kn_on_nrho = register_diag_field ('ocean_model',                      &
    'wdian_rho_pbl_ba_kn_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,              &
    'kinematic advective-form dianeutral transport from basal as binned to neutral density',&
    'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_wdian_rho_pbl_ba_kn_on_nrho > 0) compute_watermass_diag_ba=.true.

  id_tform_rho_pbl_ba_kn = register_diag_field ('ocean_model',                                &
    'tform_rho_pbl_ba_kn', Grd%tracer_axes(1:3), Time%model_time,                             &
    'kinematic advective-form water mass transform from basal pre-binning to neutral density',&
    'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_tform_rho_pbl_ba_kn > 0) compute_watermass_diag_ba=.true.

  id_tform_rho_pbl_ba_kn_on_nrho = register_diag_field ('ocean_model',                   &
    'tform_rho_pbl_ba_kn_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,           &
    'kinematic advective-form water mass transform from basal binned to neutral density',&
    'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_tform_rho_pbl_ba_kn_on_nrho > 0) compute_watermass_diag_ba=.true.

  ! kinematic advective form of runoff + calving from temperature contributions
  id_neut_temp_pbl_ba_kn = register_diag_field ('ocean_model',                 &
   'neut_temp_pbl_ba_kn', Grd%tracer_axes(1:3), Time%model_time,               &
   'temp related kinematic advective-form material time derivative from basal',&
   'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_neut_temp_pbl_ba_kn > 0) compute_watermass_diag_ba=.true.

  id_neut_temp_pbl_ba_kn_on_nrho = register_diag_field ('ocean_model',                                    &
    'neut_temp_pbl_ba_kn_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,                            &
    'temp related kinematic advective-form material time derivative from basal binned to neutral density',&
    'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_neut_temp_pbl_ba_kn_on_nrho > 0) compute_watermass_diag_ba=.true.

  id_wdian_temp_pbl_ba_kn = register_diag_field ('ocean_model',             &
    'wdian_temp_pbl_ba_kn', Grd%tracer_axes(1:3), Time%model_time,          &
    'temp related kinematic advective-form dianeutral transport from basal',&
    'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_wdian_temp_pbl_ba_kn > 0) compute_watermass_diag_ba=.true.

  id_wdian_temp_pbl_ba_kn_on_nrho = register_diag_field ('ocean_model',                               &
    'wdian_temp_pbl_ba_kn_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,                       &
    'temp related kinematic advective-form dianeutral transport from basal binned to neutral density',&
    'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_wdian_temp_pbl_ba_kn_on_nrho > 0) compute_watermass_diag_ba=.true.

  id_tform_temp_pbl_ba_kn = register_diag_field ('ocean_model',             &
    'tform_temp_pbl_ba_kn', Grd%tracer_axes(1:3), Time%model_time,          &
    'temp related kinematic advective-form water mass transform from basal',&
    'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_tform_temp_pbl_ba_kn > 0) compute_watermass_diag_ba=.true.

  id_tform_temp_pbl_ba_kn_on_nrho = register_diag_field ('ocean_model',                               &
    'tform_temp_pbl_ba_kn_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,                       &
    'temp related kinematic advective-form water mass transform from basal binned to neutral density',&
    'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_tform_temp_pbl_ba_kn_on_nrho > 0) compute_watermass_diag_ba=.true.

  ! kinematic advective form of runoff + calving from salinity contributions
  id_neut_salt_pbl_ba_kn = register_diag_field ('ocean_model',                 &
   'neut_salt_pbl_ba_kn', Grd%tracer_axes(1:3), Time%model_time,               &
   'salt related kinematic advective-form material time derivative from basal',&
   'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_neut_salt_pbl_ba_kn > 0) compute_watermass_diag_ba=.true.

  id_neut_salt_pbl_ba_kn_on_nrho = register_diag_field ('ocean_model',                                    &
    'neut_salt_pbl_ba_kn_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,                            &
    'salt related kinematic advective-form material time derivative from basal binned to neutral density',&
    'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_neut_salt_pbl_ba_kn_on_nrho > 0) compute_watermass_diag_ba=.true.

  id_wdian_salt_pbl_ba_kn = register_diag_field ('ocean_model',             &
    'wdian_salt_pbl_ba_kn', Grd%tracer_axes(1:3), Time%model_time,          &
    'salt related kinematic advective-form dianeutral transport from basal',&
    'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_wdian_salt_pbl_ba_kn > 0) compute_watermass_diag_ba=.true.

  id_wdian_salt_pbl_ba_kn_on_nrho = register_diag_field ('ocean_model',                               &
    'wdian_salt_pbl_ba_kn_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,                       &
    'salt related kinematic advective-form dianeutral transport from basal binned to neutral density',&
    'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_wdian_salt_pbl_ba_kn_on_nrho > 0) compute_watermass_diag_ba=.true.

  id_tform_salt_pbl_ba_kn = register_diag_field ('ocean_model',             &
    'tform_salt_pbl_ba_kn', Grd%tracer_axes(1:3), Time%model_time,          &
    'salt related kinematic advective-form water mass transform from basal',&
    'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_tform_salt_pbl_ba_kn > 0) compute_watermass_diag_ba=.true.

  id_tform_salt_pbl_ba_kn_on_nrho = register_diag_field ('ocean_model',                               &
    'tform_salt_pbl_ba_kn_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,                       &
    'salt related kinematic advective-form water mass transform from basal binned to neutral density',&
    'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_tform_salt_pbl_ba_kn_on_nrho > 0) compute_watermass_diag_ba=.true.

  ! eta tendency terms 
  id_eta_tend_basalmix= -1
  id_eta_tend_basalmix= register_diag_field ('ocean_model','eta_tend_basalmix', &
       Grd%tracer_axes(1:2), Time%model_time,                                   &
       'non-Bouss steric sea level tendency from basalmix', 'm/s',              &
       missing_value=missing_value, range=(/-1e10,1.e10/))
  if(id_eta_tend_basalmix > 0) compute_watermass_diag_ba = .true.

  id_eta_tend_basalmix_glob= -1
  id_eta_tend_basalmix_glob= register_diag_field ('ocean_model', 'eta_tend_basalmix_glob',&
       Time%model_time,                                                                   &
       'global mean non-bouss steric sea level tendency from basalmix',                   &
       'm/s', missing_value=missing_value, range=(/-1e10,1.e10/))
  if(id_eta_tend_basalmix_glob > 0) compute_watermass_diag_ba = .true.

    if(compute_watermass_diag_ba) then
    write(stdoutunit,'(/a/)') &
    '==>Note: running ocean_basalmix_mod w/ compute_watermass_diag_ba=.true.'
  endif


end subroutine watermass_diag_init_ba
! </SUBROUTINE> NAME="watermass_diag_init_ba"

!#######################################################################
! <SUBROUTINE NAME="watermass_diag_basal">
!
! <DESCRIPTION>
! watermass diagnostics for river = runoff + calving.
! </DESCRIPTION>
!
subroutine watermass_diag_basal(Time, Dens, T_prog, basal, &
                          temp_wrk, salt_wrk, misfkt, misfkb)

  type(ocean_time_type),        intent(in)  :: Time
  type(ocean_density_type),     intent(in)  :: Dens
  type(ocean_prog_tracer_type), intent(in)  :: T_prog(:)
  real, dimension(isd:,jsd:),   intent(in)  :: basal
  real, dimension(isd:,jsd:,:), intent(in)  :: temp_wrk
  real, dimension(isd:,jsd:,:), intent(in)  :: salt_wrk
  integer, dimension(isd:,jsd:),  intent(in) :: misfkt,misfkb ! Top and bottom input depths
  integer :: index_temp=-1
  integer :: index_salt=-1
  real, dimension(isd:ied,jsd:jed) :: eta_tend
  integer :: i,j,k,tau,n

  if (.not. compute_watermass_diag_ba) return

  tau = Time%tau

  num_prog_tracers = size(T_prog(:))
  do n=1,num_prog_tracers
     if (T_prog(n)%name == 'temp') index_temp = n
     if (T_prog(n)%name == 'salt') index_salt = n
  enddo

  wrk1(:,:,:) = 0.0
  wrk2(:,:,:) = 0.0
  wrk3(:,:,:) = 0.0
  wrk4(:,:,:) = 0.0
  wrk5(:,:,:) = 0.0

  ! flux-form contributions to material time derivative and dianeutral transport
  do k=1,nk
     do j=jsc,jec
        do i=isc,iec
           wrk1(i,j,k) = Grd%tmask(i,j,k)*(Dens%drhodT(i,j,k)*temp_wrk(i,j,k)+Dens%drhodS(i,j,k)*salt_wrk(i,j,k))
           wrk2(i,j,k) = wrk1(i,j,k)*Dens%rho_dztr_tau(i,j,k)
           wrk3(i,j,k) = wrk2(i,j,k)*Dens%stratification_factor(i,j,k)
           wrk4(i,j,k) = wrk1(i,j,k)*Grd%dat(i,j)*Dens%watermass_factor(i,j,k)
           wrk5(i,j,k) =-wrk1(i,j,k)/(epsln+Dens%rho(i,j,k,tau)**2) ! for eta_tend
        enddo
     enddo
  enddo


  call diagnose_3d(Time, Grd, id_neut_rho_basalmix, wrk2(:,:,:))
  call diagnose_3d(Time, Grd, id_wdian_rho_basalmix, wrk3(:,:,:))
  call diagnose_3d(Time, Grd, id_tform_rho_basalmix, wrk4(:,:,:))
  call diagnose_3d_rho(Time, Dens, id_neut_rho_basalmix_on_nrho, wrk2)
  call diagnose_3d_rho(Time, Dens, id_wdian_rho_basalmix_on_nrho, wrk3)
  call diagnose_3d_rho(Time, Dens, id_tform_rho_basalmix_on_nrho, wrk4)

  if(id_eta_tend_basalmix > 0 .or. id_eta_tend_basalmix_glob > 0) then
      eta_tend(:,:) = 0.0
      do k=1,nk
         do j=jsc,jec
            do i=isc,iec
               eta_tend(i,j) = eta_tend(i,j) + wrk5(i,j,k)
            enddo
         enddo
      enddo
      call diagnose_2d(Time, Grd, id_eta_tend_basalmix, eta_tend(:,:))
      call diagnose_sum(Time, Grd, Dom, id_eta_tend_basalmix_glob, eta_tend, cellarea_r)
  endif


  !-----------------------------------------------------------------------------
  ! advective-form material time derivative and associated dianeutral transport
  !
  ! assumes all river enters to the top grid cell; this is not correct when
  ! have insertion of river water into depths.

  ! kinematic form of the contributions
  wrk1(:,:,:) = 0.0
  wrk2(:,:,:) = 0.0
  wrk3(:,:,:) = 0.0
  wrk4(:,:,:) = 0.0

  !k=1
  do j=jsc,jec
     do i=isc,iec
        do k=misfkt(i,j),misfkb(i,j)
           wrk1(i,j,k) = Grd%tmask(i,j,k)*basal(i,j)                            &
                *(Dens%drhodT(i,j,k)*(0.0-T_prog(index_temp)%field(i,j,k,tau))  &
                +Dens%drhodS(i,j,k)*(0.0-T_prog(index_salt)%field(i,j,k,tau)))
           wrk2(i,j,k) = wrk1(i,j,k)*Dens%rho_dztr_tau(i,j,k)
           wrk3(i,j,k) = wrk2(i,j,k)*Dens%stratification_factor(i,j,k)
           wrk4(i,j,k) = wrk1(i,j,k)*Grd%dat(i,j)*Dens%watermass_factor(i,j,k)
        enddo
     enddo
  enddo

  call diagnose_3d(Time, Grd, id_neut_rho_pbl_ba_kn, wrk2(:,:,:))
  call diagnose_3d(Time, Grd, id_wdian_rho_pbl_ba_kn, wrk3(:,:,:))
  call diagnose_3d(Time, Grd, id_tform_rho_pbl_ba_kn, wrk4(:,:,:))
  call diagnose_3d_rho(Time, Dens, id_neut_rho_pbl_ba_kn_on_nrho, wrk2)
  call diagnose_3d_rho(Time, Dens, id_wdian_rho_pbl_ba_kn_on_nrho, wrk3)
  call diagnose_3d_rho(Time, Dens, id_tform_rho_pbl_ba_kn_on_nrho, wrk4)

  ! kinematic form of the contributions from temperature
  wrk1(:,:,:) = 0.0
  wrk2(:,:,:) = 0.0
  wrk3(:,:,:) = 0.0
  wrk4(:,:,:) = 0.0

  !k=1
  do j=jsc,jec
     do i=isc,iec
        do k=misfkt(i,j),misfkb(i,j)
           wrk1(i,j,k) = Grd%tmask(i,j,k)*basal(i,j)                            &
                *(Dens%drhodT(i,j,k)*(0.0-T_prog(index_temp)%field(i,j,k,tau))  &
                 +0.0)
           wrk2(i,j,k) = wrk1(i,j,k)*Dens%rho_dztr_tau(i,j,k)
           wrk3(i,j,k) = wrk2(i,j,k)*Dens%stratification_factor(i,j,k)
           wrk4(i,j,k) = wrk1(i,j,k)*Grd%dat(i,j)*Dens%watermass_factor(i,j,k)
        enddo
     enddo
  enddo

  call diagnose_3d(Time, Grd, id_neut_temp_pbl_ba_kn, wrk2(:,:,:))
  call diagnose_3d(Time, Grd, id_wdian_temp_pbl_ba_kn, wrk3(:,:,:))
  call diagnose_3d(Time, Grd, id_tform_temp_pbl_ba_kn, wrk4(:,:,:))
  call diagnose_3d_rho(Time, Dens, id_neut_temp_pbl_ba_kn_on_nrho, wrk2)
  call diagnose_3d_rho(Time, Dens, id_wdian_temp_pbl_ba_kn_on_nrho, wrk3)
  call diagnose_3d_rho(Time, Dens, id_tform_temp_pbl_ba_kn_on_nrho, wrk4)

   ! kinematic form of the contributions from salinity 
  wrk1(:,:,:) = 0.0
  wrk2(:,:,:) = 0.0
  wrk3(:,:,:) = 0.0
  wrk4(:,:,:) = 0.0

  !k=1
  do j=jsc,jec
     do i=isc,iec
        do k=misfkt(i,j),misfkb(i,j)
           wrk1(i,j,k) = Grd%tmask(i,j,k)*basal(i,j)                            &
                *(0.0                                                           &
                +Dens%drhodS(i,j,k)*(0.0-T_prog(index_salt)%field(i,j,k,tau)))
           wrk2(i,j,k) = wrk1(i,j,k)*Dens%rho_dztr_tau(i,j,k)
           wrk3(i,j,k) = wrk2(i,j,k)*Dens%stratification_factor(i,j,k)
           wrk4(i,j,k) = wrk1(i,j,k)*Grd%dat(i,j)*Dens%watermass_factor(i,j,k)
        enddo
     enddo
  enddo

  call diagnose_3d(Time, Grd, id_neut_salt_pbl_ba_kn, wrk2(:,:,:))
  call diagnose_3d(Time, Grd, id_wdian_salt_pbl_ba_kn, wrk3(:,:,:))
  call diagnose_3d(Time, Grd, id_tform_salt_pbl_ba_kn, wrk4(:,:,:))
  call diagnose_3d_rho(Time, Dens, id_neut_salt_pbl_ba_kn_on_nrho, wrk2)
  call diagnose_3d_rho(Time, Dens, id_wdian_salt_pbl_ba_kn_on_nrho, wrk3)
  call diagnose_3d_rho(Time, Dens, id_tform_salt_pbl_ba_kn_on_nrho, wrk4)

   ! process contributions
  wrk1(:,:,:) = 0.0
  wrk2(:,:,:) = 0.0
  wrk3(:,:,:) = 0.0
  wrk4(:,:,:) = 0.0

  !k=1
  do j=jsc,jec
     do i=isc,iec
        do k=misfkt(i,j),misfkb(i,j)
           wrk1(i,j,k) = Grd%tmask(i,j,k)*basal(i,j)                                                 &
           *(Dens%drhodT(i,j,k)*(T_prog(index_temp)%triver(i,j)-T_prog(index_temp)%field(i,j,k,tau))&
           +Dens%drhodS(i,j,k)*(T_prog(index_salt)%triver(i,j)-T_prog(index_salt)%field(i,j,k,tau)))
           wrk2(i,j,k) = wrk1(i,j,k)*Dens%rho_dztr_tau(i,j,k)
           wrk3(i,j,k) = wrk2(i,j,k)*Dens%stratification_factor(i,j,k)
           wrk4(i,j,k) = wrk1(i,j,k)*Grd%dat(i,j)*Dens%watermass_factor(i,j,k)
        enddo
     enddo
  enddo

  call diagnose_3d(Time, Grd, id_neut_rho_pbl_ba_pr, wrk2(:,:,:))
  call diagnose_3d(Time, Grd, id_wdian_rho_pbl_ba_pr, wrk3(:,:,:))
  call diagnose_3d(Time, Grd, id_tform_rho_pbl_ba_pr, wrk4(:,:,:))
  call diagnose_3d_rho(Time, Dens, id_neut_rho_pbl_ba_pr_on_nrho, wrk2)
  call diagnose_3d_rho(Time, Dens, id_wdian_rho_pbl_ba_pr_on_nrho, wrk3)
  call diagnose_3d_rho(Time, Dens, id_tform_rho_pbl_ba_pr_on_nrho, wrk4)

  ! process contribution associated with temperature 
  wrk1(:,:,:) = 0.0
  wrk2(:,:,:) = 0.0
  wrk3(:,:,:) = 0.0
  wrk4(:,:,:) = 0.0

  !k=1
  do j=jsc,jec
     do i=isc,iec
        do k=misfkt(i,j),misfkb(i,j)
           wrk1(i,j,k) = Grd%tmask(i,j,k)*basal(i,j)                                                 &
            *(Dens%drhodT(i,j,k)*(T_prog(index_temp)%triver(i,j)-T_prog(index_temp)%field(i,j,k,tau))&
            +0.0)
           wrk2(i,j,k) = wrk1(i,j,k)*Dens%rho_dztr_tau(i,j,k)
           wrk3(i,j,k) = wrk2(i,j,k)*Dens%stratification_factor(i,j,k)
           wrk4(i,j,k) = wrk1(i,j,k)*Grd%dat(i,j)*Dens%watermass_factor(i,j,k)
        enddo
     enddo
  enddo

  call diagnose_3d(Time, Grd, id_neut_temp_pbl_ba_pr, wrk2(:,:,:))
  call diagnose_3d(Time, Grd, id_wdian_temp_pbl_ba_pr, wrk3(:,:,:))
  call diagnose_3d(Time, Grd, id_tform_temp_pbl_ba_pr, wrk4(:,:,:))
  call diagnose_3d_rho(Time, Dens, id_neut_temp_pbl_ba_pr_on_nrho, wrk2)
  call diagnose_3d_rho(Time, Dens, id_wdian_temp_pbl_ba_pr_on_nrho, wrk3)
  call diagnose_3d_rho(Time, Dens, id_tform_temp_pbl_ba_pr_on_nrho, wrk4)

  ! process contribution associated with salinity 
  wrk1(:,:,:) = 0.0
  wrk2(:,:,:) = 0.0
  wrk3(:,:,:) = 0.0
  wrk4(:,:,:) = 0.0
  !k=1
  do j=jsc,jec
     do i=isc,iec
        do k=misfkt(i,j),misfkb(i,j)
           wrk1(i,j,k) = Grd%tmask(i,j,k)*basal(i,j)                                                 &
            *(0.0                                                                                    &
            +Dens%drhodS(i,j,k)*(T_prog(index_salt)%triver(i,j)-T_prog(index_salt)%field(i,j,k,tau)))
           wrk2(i,j,k) = wrk1(i,j,k)*Dens%rho_dztr_tau(i,j,k)
           wrk3(i,j,k) = wrk2(i,j,k)*Dens%stratification_factor(i,j,k)
           wrk4(i,j,k) = wrk1(i,j,k)*Grd%dat(i,j)*Dens%watermass_factor(i,j,k)
        enddo
     enddo
  enddo

  call diagnose_3d(Time, Grd, id_neut_salt_pbl_ba_pr, wrk2(:,:,:))
  call diagnose_3d(Time, Grd, id_wdian_salt_pbl_ba_pr, wrk3(:,:,:))
  call diagnose_3d(Time, Grd, id_tform_salt_pbl_ba_pr, wrk4(:,:,:))
  call diagnose_3d_rho(Time, Dens, id_neut_salt_pbl_ba_pr_on_nrho, wrk2)
  call diagnose_3d_rho(Time, Dens, id_wdian_salt_pbl_ba_pr_on_nrho, wrk3)
  call diagnose_3d_rho(Time, Dens, id_tform_salt_pbl_ba_pr_on_nrho, wrk4)

end subroutine watermass_diag_basal
! </SUBROUTINE> NAME="watermass_diag_basal"


!#######################################################################
! <SUBROUTINE NAME="river_kappa">
!
! <DESCRIPTION>
! This subroutine enhances the vertical diffusivity kappa over
! a vertical column whose thickness is set by river_diffusion_thickness
! and whose horizontal location is given by the rmask array.
! Note that rmask can be > 0 even if river=0 in the case when
! use virtual salt flux.
! The enhanced diffusivity is maximum at the top cell and is linearly
! interpolated to the normal diffusivity at the depth set by
! river_diffusion_thickness
! </DESCRIPTION>
!
!subroutine river_kappa (Time, Thickness, Tracer, kappa)

!  type(ocean_time_type),         intent(in)    :: Time
!  type(ocean_thickness_type),    intent(in)    :: Thickness
!  type(ocean_prog_tracer_type),  intent(in)    :: Tracer
!  real, dimension(isd:,jsd:,:),  intent(inout) :: kappa
!  real    :: river_diffusivity=0.0      ! enhancement to the vertical diffusivity (m^2/s) at river mouths
!  real    :: river_diffusion_thickness=0.0 ! static thickness (m) of ocean column where diffuse tracer at river mouths.
!                                         ! actual thickness is based on model grid spacing. min thickness=dtz(1).  integer :: i, j, k, nz
!  real    :: depth
!  real, dimension(nk) :: zw_ij
!
!  if(.not. module_is_initialized ) then
!    call mpp_error(FATAL, &
!    '==>Error in ocean_rivermix_mod (river_kappa): module must be initialized')
!  endif
!
!  wrk1=0.0
!
!  do j=jsc,jec
!     do i=isc,iec
!
!        do k=1,nk
!           zw_ij(k) = Thickness%depth_zwt(i,j,k)
!        enddo
!
!        if (Tracer%riverdiffuse(i,j) > 0.0 .and. Grd%kmt(i,j) > 0) then
!            depth = min(Grd%ht(i,j),river_diffusion_thickness)
!            nz    = min(Grd%kmt(i,j),floor(frac_index(depth,zw_ij)))
!            nz    = max(1,nz)
!            do k=1,nz
!              wrk1(i,j,k)  = river_diffusivity*(1.0 - Thickness%depth_zt(i,j,k)/Thickness%depth_zt(i,j,nk))
!              kappa(i,j,k) = kappa(i,j,k) + wrk1(i,j,k)
!            enddo
!        endif
!     enddo
!  enddo
!
!  !if (id_diff_cbt_river_t > 0 .and. Tracer%name=='temp') then
!  !   call diagnose_3d(Time, Grd, id_diff_cbt_river_t, wrk1(:,:,:))
!  !endif
!  !if (id_diff_cbt_river_s > 0 .and. Tracer%name=='salt') then
!  !   call diagnose_3d(Time, Grd, id_diff_cbt_river_s, wrk1(:,:,:))
!  !endif
!
!end subroutine river_kappa
!! </SUBROUTINE> NAME="river_kappa"


!#######################################################################
! <SUBROUTINE NAME="eos_fzp">
!
! <DESCRIPTION>
! Freezing point temperature for sea water
! </DESCRIPTION>
!

!! </SUBROUTINE> NAME="eos_fzp"
!   SUBROUTINE  eos_fzp( psal, ptf, pdep )
!      !!----------------------------------------------------------------------
!      !!                 ***  ROUTINE eos_fzp  ***
!      !!
!      !! ** Purpose :   Compute the freezing point temperature [Celsius]
!      !!
!      !! ** Method  :   UNESCO freezing point (ptf) in Celsius is given by
!      !!       ptf(t,z) = (-.0575+1.710523e-3*sqrt(abs(s))-2.154996e-4*s)*s - 7.53e-4*z
!      !!       checkvalue: tf=-2.588567 Celsius for s=40psu, z=500m
!      !!
!      !! Reference  :   UNESCO tech. papers in the marine science no. 28. 1978
!      !!----------------------------------------------------------------------
!      REAL(wp), DIMENSION(jpi,jpj), INTENT(in   )           ::   psal   ! salinity   [psu]
!      REAL(wp), DIMENSION(jpi,jpj), INTENT(in   ), OPTIONAL ::   pdep   ! depth      [m]
!      REAL(wp), DIMENSION(jpi,jpj), INTENT(out  )           ::   ptf    ! freezing temperature [Celsius]
!      !
!      INTEGER  ::   ji, jj, jformulation          ! dummy loop indices
!      REAL(wp) ::   zt, zs, z1_S0   ! local scalars
!      !!----------------------------------------------------------------------
!      !
!      jformulation = 1
!      !
!      IF ( jformulation == 1 ) THEN      !==  CT,SA (TEOS-10 and S-EOS formulations) ==!
!         !
!         z1_S0 = 1._wp / 35.16504_wp
!         DO jj = 1, jpj
!            DO ji = 1, jpi
!               zs= SQRT( ABS( psal(ji,jj) ) * z1_S0 )           ! square root salinity
!               ptf(ji,jj) = ((((1.46873e-03_wp*zs-9.64972e-03_wp)*zs+2.28348e-02_wp)*zs &
!                  &          - 3.12775e-02_wp)*zs+2.07679e-02_wp)*zs-5.87701e-02_wp
!            END DO
!         END DO
!         ptf(:,:) = ptf(:,:) * psal(:,:)
!         !
!         IF( PRESENT( pdep ) )   ptf(:,:) = ptf(:,:) - 7.53e-4 * pdep(:,:)
!         !
!      ELSEIF ( jformulation == 2) THEN                !==  PT,SP (UNESCO formulation)  ==!
!         !
!         ptf(:,:) = ( - 0.0575_wp + 1.710523e-3_wp * SQRT( psal(:,:) )   &
!            &                     - 2.154996e-4_wp *       psal(:,:)   ) * psal(:,:)
!            !
!         IF( PRESENT( pdep ) )   ptf(:,:) = ptf(:,:) - 7.53e-4 * pdep(:,:)
!         !
!      ENDIF
!      !
!  END SUBROUTINE eos_fzp
!! </SUBROUTINE> NAME="eos_fzp"


end module ocean_basal_tracer_mod
