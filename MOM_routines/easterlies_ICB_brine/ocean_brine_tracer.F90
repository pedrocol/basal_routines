module ocean_brine_tracer_mod
!</CONTACT>
!
!<CONTACT EMAIL="paul.spence@gmail.com"> Paul Spence
!</CONTACT>
!
!<CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> S. M. Griffies 
!</CONTACT>
!
!<OVERVIEW>
! Distribute brine rejection at depth.
!</OVERVIEW>
!
!<DESCRIPTION>
! This module distributed at depth brine rejection. Brine rejection
! happens as a change in concentration (and not as a salt flux).
! The array brine(i,j) equals the array wfi_form of the sea-ice model
! (ocean_sbc).
! Different strategies can be followed to distribute brine at depth.
! In this module vertical homogeneous and Barthelemy et al., 2015
! http://dx.doi.org/10.1016/j.ocemod.2014.12.009
!
!<NAMELIST NAME="ocean_brine_tracer_nml">
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

implicit none

private

#include <ocean_memory.h>

type(ocean_domain_type), pointer :: Dom => NULL()
type(ocean_grid_type),   pointer :: Grd => NULL()

public ocean_brine_tracer_init
public brine_tracer_source

character(len=126)  :: version = '$Id: ocean_brine_tracer.F90,v 20.0 2013/12/14 00:16:24 fms Exp $'
character (len=128) :: tagname = '$Name: tikal $'

! for diagnostics
logical :: used
integer, dimension(:), allocatable :: id_basal_tend
integer :: id_brine_fwflx        =-1
integer :: id_brine_fwflx2d      =-1

integer :: num_prog_tracers      = 0
logical :: module_is_initialized = .FALSE.
logical :: use_basal_module       = .true.
logical :: use_brine_module       = .true.
logical :: use_icb_module       = .true.
logical :: test_nml = .false.

namelist /ocean_basal_tracer_nml/ use_basal_module, use_icb_module, use_brine_module
          

contains

!#######################################################################
! <SUBROUTINE NAME="ocean_brine_tracer_init">
!
! <DESCRIPTION>
! This subroutine is intended to be used to initialize addition of basal melt
! water fluxes from an input file. 
! Everything in this subroutine is a user prototype, and should be replacable.
! </DESCRIPTION>
!
subroutine ocean_brine_tracer_init(Grid, Domain, Time, T_prog, dtime, Ocean_options, Dens)

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
    '==>Error in ocean_brine_tracer_mod (ocean_brine_tracer_init): module already initialized')
  endif

  module_is_initialized = .TRUE.

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

  dtimer = 1.0/dtime

  if(use_brine_module) then
      write(stdoutunit,*)'==>Note from ocean_brine_tracer_mod: Using this module.'
      Ocean_options%ocean_brine_tracer= 'Used ocean tracer brine.'
  else
      write(stdoutunit,*)'==>Note from ocean_brine_tracer_mod: NOT using ocean tracer brine.'
      Ocean_options%ocean_brine_tracer= 'Did NOT use ocean tracer brine.'
      return
  endif

  ! register diagnostic outputs
  id_brine_fwflx = register_diag_field ('ocean_model','brine_fwflx', Grd%tracer_axes(1:3),   &
              Time%model_time, '3d mass flux of liquid brine meltwater leaving ocean ',        &
              '(kg/m^3)*(m/sec)', missing_value=missing_value,range=(/-1e6,1e6/),            &
              standard_name='water_flux3d_out_sea_water_from_brine_melting')
  id_brine_fwflx2d = register_diag_field ('ocean_model','brine_fwflx2d', Grd%tracer_axes(1:2),   &
              Time%model_time, '2d mass flux of liquid brine meltwater leaving ocean ',        &
              '(kg/m^3)*(m/sec)', missing_value=missing_value,range=(/-1e6,1e6/),            &
              standard_name='water_flux2d_out_sea_water_from_brine_melting')

end subroutine ocean_brine_tracer_init
! </SUBROUTINE> NAME="ocean_brine_tracer_init"

!#######################################################################
! <SUBROUTINE NAME="brine_tracer_source">
!
! <DESCRIPTION>
! This subroutine calculates thickness weighted and density weighted
! time tendencies of tracers due to damping by brine.
! </DESCRIPTION>
!
subroutine brine_tracer_source(Time, Time_steps, Thickness, Dens, T_prog, brine, index_temp, &
                               index_salt, brine3d)

  type(ocean_time_type),          intent(in)      :: Time
  type(ocean_time_steps_type),    intent(in)      :: Time_steps
  type(ocean_thickness_type),     intent(inout)   :: Thickness
  type(ocean_density_type),       intent(in)      :: Dens
  type(ocean_prog_tracer_type),   intent(inout)   :: T_prog(:)
  real, dimension(isd:,jsd:),     intent(inout)   :: brine
  real, dimension(isd:,jsd:,:)  , intent(inout)   :: brine3d
  integer,                        intent(in)      :: index_temp
  integer,                        intent(in)      :: index_salt
  integer :: i, j, k, n
  integer :: param_choice
  integer :: taum1, tau, taup1
  real    :: maxinsertiondepth, depth, thkocean
  real    :: delta(nk)
  integer :: max_nk


  if(.not. use_brine_module) return

  param_choice = 1

  IF ( param_choice == 1 ) THEN !Uniform dstribution 200m
      do j=jsc,jec
         do i=isc,iec

            if (brine(i,j) < 0.0 .and. Grd%kmt(i,j) > 0) then
               thkocean = 0.0

               maxinsertiondepth = 200.0
               depth  = min(Grd%ht(i,j),maxinsertiondepth)                ! be sure not to discharge river content into rock, ht = ocean topography
               max_nk = min(Grd%kmt(i,j),floor(frac_index(depth,Grd%zw))) ! max number of k-levels into which discharge rivers
               do k=1,max_nk
                  thkocean = thkocean + Thickness%rho_dzt(i,j,k,tau)
               enddo

               do k=1,max_nk
                  !Brine rejected is performed via change in the concentration (it's not a salt flux)
                  delta(k) = Thickness%rho_dzt(i,j,k,tau)/(epsln+thkocean)
                  brine3d(i,j,k) = brine(i,j)*delta(k)
                  Thickness%mass_source(i,j,k) = Thickness%mass_source(i,j,k) + brine3d(i,j,k)
               enddo
            endif

         enddo
      enddo

          
  ELSEIF ( param_choice == 2 ) THEN !Barthelemy et al., 2015
          !Do nothing for the moment
  ENDIF
  

end subroutine brine_tracer_source
! </SUBROUTINE> NAME="brine_tracer_source"









end module ocean_brine_tracer_mod