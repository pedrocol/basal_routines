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
! on ocean_sponges_tracer.F90. The basal meltwater input file is from 
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
use ocean_workspace_mod,      only: wrk1, wrk2
use ocean_util_mod,           only: diagnose_3d
use constants_mod,            only: epsln
use axis_utils_mod,           only: frac_index, nearest_index

implicit none

private

#include <ocean_memory.h>

type ocean_basal_type
   integer :: id                                             ! time_interp_external index
   character(len=32) :: name                                 ! tracer name corresponding to basal
   real, dimension(:,:,:), pointer :: damp_coeff   => NULL() ! 3d inverse damping rate (tracer units/ sec)
   real, dimension(:,:)  , pointer :: damp_coeff2d => NULL() ! 2d inverse damping rate (tracer units/ sec)
   real, dimension(:,:)  , pointer :: basal_fw     => NULL() ! runoff
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

integer :: num_prog_tracers      = 0
logical :: module_is_initialized = .FALSE.
logical :: damp_coeff_3d         = .false. 
logical :: use_this_module       = .false. 
logical :: test_nml              = .false. 

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


namelist /ocean_basal_tracer_nml/ use_this_module, test_nml, damp_coeff_3d

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
subroutine ocean_basal_tracer_init(Grid, Domain, Time, T_prog, dtime, Ocean_options)

  type(ocean_grid_type),        intent(in), target :: Grid
  type(ocean_domain_type),      intent(in), target :: Domain
  type(ocean_time_type),        intent(in)         :: Time
  type(ocean_prog_tracer_type), intent(in)         :: T_prog(:)
  real,                         intent(in)         :: dtime
  type(ocean_options_type),     intent(inout)      :: Ocean_options

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

  !num_prog_tracers = size(T_prog(:))
  num_prog_tracers = 1

  allocate( Basal(1) )

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

  if(use_this_module) then
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
    Basal(n)%id = init_external_field(name,'basal_fw',domain=Domain%domain2d)
    if (Basal(n)%id < 1) then
      call mpp_error(FATAL,&
      '==>Error: in ocean_basal_tracer_mod:  basal fw values are not specified')
    endif
    write(stdoutunit,*) '==> Using basal freshwater flux data specified from file '//trim(name)
  enddo

  ! register diagnostic outputs
  id_basal_fwflx = register_diag_field('ocean_model','basal_fwflx', Grd%tracer_axes(1:3),&
       Time%model_time, 'mass flux of liquid basal meltwater entering ocean ',    &
       '(kg/m^3)*(m/sec)', missing_value=missing_value,range=(/-1e6,1e6/),     &
       standard_name='water_flux_into_sea_water_from_basal_melting')
  
  allocate (id_basal_tend(num_prog_tracers))
  id_basal_tend = -1

  do n=1,num_prog_tracers
     id_basal_tend(n) = register_diag_field ('ocean_model', 'basal_fw_tend',&
               Grd%tracer_axes(1:3), Time%model_time, 'rho*dzt*tendency due to basal freshwater flux',          &
               trim(T_prog(n)%flux_units), missing_value=missing_value, range=(/-1.e10,1.e10/))
  enddo
  
  !try to write out basal fwflx data to output for test
  ! get sponge value for current time
  wrk1=0.0
  call time_interp_external(Basal(1)%id, Time%model_time, wrk1)
  call diagnose_3d(Time, Grd, id_basal_fwflx,wrk1(:,:,:))

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
subroutine basal_tracer_source(Time, Thickness, T_prog, dtime)

  type(ocean_time_type),        intent(in)    :: Time
  type(ocean_thickness_type),   intent(in)    :: Thickness
  type(ocean_prog_tracer_type), intent(inout) :: T_prog(:)
  real,                         intent(in)         :: dtime

  integer :: param_choice

  param_choice = 0

  IF ( param_choice == 1 ) THEN
    CALL basal_tracer_source_0(Time, Thickness, T_prog)
  ELSEIF ( param_choice == 2 ) THEN
    CALL basal_tracer_source_1(Time, Thickness, T_prog, dtime)
  ELSEIF ( param_choice == 3 ) THEN
    CALL basal_tracer_source_paul(Time, Thickness, T_prog)
  ENDIF

end subroutine basal_tracer_source
! </SUBROUTINE> NAME="basal_tracer_source"

!#######################################################################
! <SUBROUTINE NAME="basal_tracer_source_0">
!
! <DESCRIPTION>
! This subroutine calculates thickness weighted and density weighted
! time tendencies of tracers due to damping by basal.
! </DESCRIPTION>
!
subroutine basal_tracer_source_0(Time, Thickness, T_prog)
  ! Case: specified fwf and heat flux forcing beneath the ice shelf
  type(ocean_time_type),        intent(in)    :: Time
  type(ocean_thickness_type),   intent(in)    :: Thickness
  type(ocean_prog_tracer_type), intent(inout) :: T_prog(:)
  integer :: taum1, tau
  integer :: i, j, k, n
  integer, allocatable, dimension(:,:) :: misfkt,misfkb ! Top and bottom input depths
  integer, allocatable, dimension(:,:) :: fwfisf        ! fresh water flux from the isf (fwfisf <0 mean melting) [Kg/m2/s]
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
  integer ::   ikt, ikb     ! local integers
!#######################################################################


  allocate ( misfkt(isd:ied,jsd:jed), misfkb(isd:ied,jsd:jed) )
  allocate ( fwfisf(isd:ied,jsd:jed),   qisf(isd:ied,jsd:jed) )
  allocate (   stbl(isd:ied,jsd:jed), zt_frz(isd:ied,jsd:jed) )
  allocate (   zt_frz3d(isd:ied,jsd:jed,nk) )
  allocate ( rhisf_tbl(isd:ied,jsd:jed), r1_hisf_tbl(isd:ied,jsd:jed) )
  allocate ( risf_tsc_b(isd:ied,jsd:jed,2), risf_tsc(isd:ied,jsd:jed,2) ) !1=temp 2=sal
  allocate ( risf_tsc_3d_b(isd:ied,jsd:jed,nk,2), risf_tsc_3d(isd:ied,jsd:jed,nk,2) ) !1=temp 2=sal

  taum1 = Time%taum1
  tau   = Time%tau
  wrk1  = 0.0 !Asumes T=0 and S=0

  !do n=1,size(T_prog(:))
  do n=1,num_prog_tracers

    ! Need to reinitialise wrk2 here due to limiting of temperature.
    wrk2  = 0.0

    fwfisf(:,:) = 1 ! fresh water flux from the isf (fwfisf <0 mean melting)

    rLfusisf    = 0.334e6    !: latent heat of fusion of ice shelf     [J/kg]
    qisf(:,:)   = fwfisf(:,:) * rLfusisf               ! heat flux

    rau0        = 1026                 !: volumic mass of reference     [kg/m3]
    rcp         = 3991.86795711963      !: heat capacity     [J/K]
    rau0_rcp    = rau0 * rcp
    r1_rau0_rcp = 1. / rau0_rcp
    r1_rau0     = 1. / rau0

    soce        =   34.7        ! Sea salinity, approximation Pedro
    stbl(:,:)   = soce          ! Salinity top boundary layer

    !Freezing point temp
    !DO jj = 1,jpj
    !   DO ji = 1,jpi
    !      zdep(ji,jj)=gdepw_n(ji,jj,misfkt(ji,jj))
    !   END DO
    !END DO

    !CALL eos_fzp( stbl(:,:), zt_frz(:,:), zdep(:,:) ) ! freezing point temperature at depth z
    zt_frz(:,:) = 273.15 
    zt_frz3d(:,:,:) = 273.15 

    ! Before and now values
    risf_tsc(:,:,1) = qisf(:,:) * r1_rau0_rcp - fwfisf(:,:) * zt_frz(:,:) * r1_rau0 !:before and now T & S isf contents [K.m/s & PSU.m/s]
    risf_tsc(:,:,2) = 0.0
    risf_tsc_b(:,:,:)= risf_tsc(:,:,:) !Equal for the moment, constant source
    ! 10 *  ((0.334e6 * (1/(1026*3991.86795711963)) ) - ( 272 * (1/1026)))

    !Some dummy values for ice shelfs geometry
    misfkt(:,:) = 4
    misfkb(:,:) = 6
    !rhisf_tbl(:,:) = SUM(Thickness%dzt(:,:,misfkt:misfkb))
    rhisf_tbl(:,:) = 50 !Thickness%dzt(:,:,4)
    r1_hisf_tbl(:,:) = 1. / rhisf_tbl(:,:)


    do j=jsc,jec
       do i=isc,iec
          ikt = misfkt(i,j)
          ikb = misfkb(i,j)
          ! level fully include in the ice shelf boundary layer
          ! sign - because fwf sign of evapo (rnf sign of precip)
          do k = ikt, ikb - 1
            risf_tsc_3d(i,j,k,1) = (qisf(i,j) * r1_rau0_rcp - fwfisf(i,j) * zt_frz3d(i,j,k) * r1_rau0) * r1_hisf_tbl(i,j) ! K/s
            risf_tsc_3d(i,j,k,2) = 0.0
          enddo
       enddo
    enddo
    risf_tsc_3d_b(:,:,:,:)= risf_tsc_3d(:,:,:,:) !Equal for the moment, constant source

!    ! compute bottom level of isf tbl and thickness of tbl below the ice shelf
!    DO jj = 1,jpj
!        DO ji = 1,jpi
!           ikt = misfkt(ji,jj)
!           ikb = misfkt(ji,jj)
!           ! thickness of boundary layer at least the top level thickness
!           rhisf_tbl(ji,jj) = MAX(rhisf_tbl_0(ji,jj), e3t_n(ji,jj,ikt))
!            ! determine the deepest level influenced by the boundary layer
!           DO jk = ikt+1, mbkt(ji,jj)
!              IF( (SUM(e3t_n(ji,jj,ikt:jk-1)) < rhisf_tbl(ji,jj)) .AND. (tmask(ji,jj,jk) == 1) )   ikb = jk
!           END DO
!           rhisf_tbl(ji,jj) = MIN(rhisf_tbl(ji,jj), SUM(e3t_n(ji,jj,ikt:ikb))) ! limit the tbl to water thickness.
!           misfkb(ji,jj) = ikb                                                   ! last wet level of the tbl
!           r1_hisf_tbl(ji,jj) = 1._wp / rhisf_tbl(ji,jj)
!
!           zhk           = SUM( e3t_n(ji, jj, ikt:ikb - 1)) * r1_hisf_tbl(ji,jj) ! proportion of tbl cover by cell from ikt to ikb - 1
!           ralpha(ji,jj) = rhisf_tbl(ji,jj) * (1._wp - zhk ) / e3t_n(ji,jj,ikb)  ! proportion of bottom cell influenced by boundary layer
!        END DO
!     END DO

    if (trim(T_prog(n)%name) == 'temp' ) then
       do j=jsc,jec
          do i=isc,iec
             ikt = misfkt(i,j)
             ikb = misfkb(i,j)
             ! level fully include in the ice shelf boundary layer
             ! sign - because fwf sign of evapo (rnf sign of precip)
             do k = ikt, ikb - 1
                wrk2(i,j,k) = wrk2(i,j,k) + 0.5 * ( risf_tsc_3d_b(i,j,k,1) + risf_tsc_3d(i,j,k,1) )
             enddo
          enddo
       enddo
    end if

    if ( trim(T_prog(n)%name) == 'salt' ) then
       do j=jsc,jec
          do i=isc,iec
             ikt = misfkt(i,j)
             ikb = misfkb(i,j)
             ! level fully include in the ice shelf boundary layer
             ! sign - because fwf sign of evapo (rnf sign of precip)
             do k = ikt, ikb - 1
                wrk2(i,j,k) = wrk2(i,j,k) + 0.5 * ( risf_tsc_3d_b(i,j,k,2) + risf_tsc_3d(i,j,k,2) )
             enddo
          enddo
       enddo   
    end if

    ! Update tendency
    do k = 1, nk
        do j = jsc, jec
            do i = isc, iec
                T_prog(n)%th_tendency(i,j,k) = T_prog(n)%th_tendency(i,j,k) &
                                              + wrk2(i,j,k)
            end do
        end do
    end do


    !if (id_basal_tend(n) > 0) call diagnose_3d(Time, Grd, id_basal_tend(n), &
    !     T_prog(n)%conversion*wrk2(:,:,:))

  enddo !n

  deallocate ( misfkt, misfkb )
  deallocate ( fwfisf,   qisf )
  deallocate (   stbl, zt_frz )
  deallocate ( zt_frz3d )
  deallocate ( rhisf_tbl, r1_hisf_tbl )
  deallocate ( risf_tsc_b, risf_tsc ) !1=temp 2=sal
  deallocate ( risf_tsc_3d_b, risf_tsc_3d )
end subroutine basal_tracer_source_0
! </SUBROUTINE> NAME="basal_tracer_source_0"

!#######################################################################
! <SUBROUTINE NAME="basal_tracer_source_1">
!
! <DESCRIPTION>
! This subroutine calculates thickness weighted and density weighted
! time tendencies of tracers due to damping by basal.
! </DESCRIPTION>
!
subroutine basal_tracer_source_1(Time, Thickness, T_prog, dtime)
  ! Case: specified fwf and heat flux forcing beneath the ice shelf
  type(ocean_time_type),        intent(in)    :: Time
  type(ocean_thickness_type),   intent(in)    :: Thickness
  type(ocean_prog_tracer_type), intent(inout) :: T_prog(:)
  real,                         intent(in)         :: dtime
  integer :: taum1, tau
  integer :: i, j, k, n, nz
  integer, allocatable, dimension(:,:) :: misfkt,misfkb ! Top and bottom input depths
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
  integer ::   ikt, ikb     ! local integers
  !rivermix
  real    :: depth, thkocean
  real    :: delta(nk), delta_rho_tocean(nk), delta_rho0_triver(nk)
  real    :: zextra, zinsert, tracerextra, tracernew(nk)
  real    :: tracer_input
  real    :: maxinsertiondepth
  real, dimension(:,:), allocatable :: tracer_flux

!#######################################################################

  allocate ( misfkt(isd:ied,jsd:jed), misfkb(isd:ied,jsd:jed) )
  allocate ( fwfisf(isd:ied,jsd:jed),   qisf(isd:ied,jsd:jed) )
  allocate (   stbl(isd:ied,jsd:jed), zt_frz(isd:ied,jsd:jed) )
  allocate (   zt_frz3d(isd:ied,jsd:jed,nk) )
  allocate ( rhisf_tbl(isd:ied,jsd:jed), r1_hisf_tbl(isd:ied,jsd:jed) )
  allocate ( risf_tsc_b(isd:ied,jsd:jed,2), risf_tsc(isd:ied,jsd:jed,2) ) !1=temp 2=sal
  allocate ( risf_tsc_3d_b(isd:ied,jsd:jed,nk,2), risf_tsc_3d(isd:ied,jsd:jed,nk,2) ) !1=temp 2=sal
  allocate ( tracer_flux(isd:ied,jsd:jed) )

  taum1 = Time%taum1
  tau   = Time%tau
  wrk1  = 0.0 !Asumes T=0 and S=0

  do n=1,num_prog_tracers
    T_prog(n)%wrk1(:,:,:) = 0.0
  enddo
  delta             = 0.0
  delta_rho_tocean  = 0.0
  delta_rho0_triver = 0.0
  tracer_flux       = 0.0

  if (Basal(1)%id > 0) then
     ! get basal value for current time
     call time_interp_external(Basal(1)%id, Time%model_time, wrk1)
  endif

  fwfisf(:,:) = 0 ! fresh water flux from the isf (fwfisf <0 mean melting)
  !fwfisf(:,:) = 1 ! fresh water flux from the isf (fwfisf <0 mean melting)
  !do k = 1,nk
    do j=jsd,jed
       do i=isd,ied
          fwfisf(i,j)=wrk1(i,j,1)
       enddo
    enddo
  !enddo

  rLfusisf    = 0.334e6    !: latent heat of fusion of ice shelf     [J/kg]
  qisf(:,:)   = fwfisf(:,:) * rLfusisf               ! heat flux

  rau0        = 1026                 !: volumic mass of reference     [kg/m3]
  rcp         = 3991.86795711963      !: heat capacity     [J/K]
  rau0_rcp    = rau0 * rcp
  r1_rau0_rcp = 1. / rau0_rcp
  r1_rau0     = 1. / rau0

  soce        =   34.7        ! Sea salinity, approximation Pedro
  stbl(:,:)   = soce          ! Salinity top boundary layer


  !CALL eos_fzp( stbl(:,:), zt_frz(:,:), zdep(:,:) ) ! freezing point temperature at depth z
  zt_frz(:,:) = 273.15
  zt_frz3d(:,:,:) = 273.15

  ! Before and now values
  risf_tsc(:,:,1) = qisf(:,:) * r1_rau0_rcp - fwfisf(:,:) * zt_frz(:,:) * r1_rau0 !:before and now T & S isf contents [K.m/s & PSU.m/s]
  risf_tsc(:,:,2) = 0.0
  risf_tsc_b(:,:,:)= risf_tsc(:,:,:) !Equal for the moment, constant source
  ! 10 *  ((0.334e6 * (1/(1026*3991.86795711963)) ) - ( 272 * (1/1026)))

  !Some dummy values for ice shelfs geometry
  misfkt(:,:) = 4
  misfkb(:,:) = 8
  !rhisf_tbl(:,:) = SUM(Thickness%dzt(:,:,misfkt:misfkb))
  rhisf_tbl(:,:) = 50 !Thickness%dzt(:,:,4)
  r1_hisf_tbl(:,:) = 1. / rhisf_tbl(:,:)

  do j=jsc,jec
     do i=isc,iec
        ikt = misfkt(i,j)
        ikb = misfkb(i,j)
        ! level fully include in the ice shelf boundary layer
        ! sign - because fwf sign of evapo (rnf sign of precip)
        do k = ikt, ikb - 1
           risf_tsc_3d(i,j,k,1) = (qisf(i,j) * r1_rau0_rcp - fwfisf(i,j) * zt_frz3d(i,j,k) * r1_rau0) * r1_hisf_tbl(i,j) ! K/s
           risf_tsc_3d(i,j,k,2) = 0.0
        enddo
     enddo
  enddo
  risf_tsc_3d_b(:,:,:,:)= risf_tsc_3d(:,:,:,:) !Equal for the moment, constant source


  ! Need to reinitialise wrk2 here due to limiting of temperature.
  wrk2  = 0.0

  do j=jsc,jec
     do i=isc,iec

        if (fwfisf(i,j) > 0.0 .and. Grd%kmt(i,j) > 0) then


        ! the array "river" contains the volume rate (m/s) or mass
        ! rate (kg/m2/sec) of fluid with tracer 
        ! that is to be distributed in the vertical. 
           maxinsertiondepth = Grd%zt(misfkb(i,j))
           depth = min(Grd%ht(i,j),maxinsertiondepth)     ! be sure not to discharge river content into rock, ht = ocean topography
           misfkb(i,j) = min(Grd%kmt(i,j),floor(frac_index(depth,Grd%zw))) ! max number of k-levels into which discharge rivers
           misfkb(i,j)    = max(1,misfkb(i,j))                                         ! make sure have at least one cell to discharge into

           ! determine fractional thicknesses of grid cells 
           thkocean = 0.0
           do k=misfkt(i,j),misfkb(i,j)
              thkocean = thkocean + Thickness%rho_dzt(i,j,k,tau)
           enddo
           do k=misfkt(i,j),misfkb(i,j)
              delta(k) = Thickness%rho_dzt(i,j,k,tau)/(epsln+thkocean)
           enddo

           do n=1,num_prog_tracers
              tracer_flux(i,j) = fwfisf(i,j)

              zextra=0.0
              do k=misfkb(i,j),misfkt(i,j),-1
                 tracernew(k) = 0.0

                 if (k.eq.misfkb(i,j)) then
                     tracerextra=0.0
                 else
                     tracerextra = tracernew(k+1)
                 endif

                 zinsert = tracer_flux(i,j)*dtime*delta(k)
                 tracernew(k) = (tracerextra*zextra + T_prog(n)%field(i,j,k,tau)*Thickness%rho_dzt(i,j,k,tau) + &
                                 & T_prog(n)%triver(i,j)*zinsert) / (zextra+Thickness%rho_dzt(i,j,k,tau)+zinsert)
                 zextra=zextra+zinsert
              enddo

              k=1 !Treatment at the surface
              if ( misfkt(i,j) == k ) then
                 T_prog(n)%wrk1(i,j,k) = (tracernew(k)*(Thickness%rho_dzt(i,j,k,tau)+tracer_flux(i,j)*dtime) -&
                                       T_prog(n)%field(i,j,k,tau)*Thickness%rho_dzt(i,j,k,tau))/dtime
              endif 

              do k=misfkt(i,j),misfkb(i,j)
                 T_prog(n)%wrk1(i,j,k) = Thickness%rho_dzt(i,j,k,tau)*(tracernew(k) - T_prog(n)%field(i,j,k,tau))/dtime !Tendency 
              enddo

              if(debug_all_in_top_cell) then
                  k=1
                  T_prog(n)%wrk1(i,j,:) = 0.0
                  T_prog(n)%wrk1(i,j,k) = Grd%tmask(i,j,k)*tracer_flux(i,j)*T_prog(n)%triver(i,j)
              endif

              do k=misfkt(i,j),misfkb(i,j)
                 T_prog(n)%th_tendency(i,j,k) = T_prog(n)%th_tendency(i,j,k) + T_prog(n)%wrk1(i,j,k)
              enddo

           enddo !n
        endif ! fwfisf > 0
     enddo !i
  enddo !j

  deallocate ( misfkt, misfkb )
  deallocate ( fwfisf,   qisf )
  deallocate (   stbl, zt_frz )
  deallocate ( zt_frz3d )
  deallocate ( rhisf_tbl, r1_hisf_tbl )
  deallocate ( risf_tsc_b, risf_tsc ) !1=temp 2=sal
  deallocate ( risf_tsc_3d_b, risf_tsc_3d )
  deallocate ( tracer_flux )

end subroutine basal_tracer_source_1
! </SUBROUTINE> NAME="basal_tracer_source_1"


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


!#######################################################################
! <SUBROUTINE NAME="basal_tracer_source_paul">
!
! <DESCRIPTION>
! This subroutine calculates thickness weighted and density weighted
! time tendencies of tracers due to damping by basal.
! </DESCRIPTION>
!
subroutine basal_tracer_source_paul(Time, Thickness, T_prog)

  type(ocean_time_type),        intent(in)    :: Time
  type(ocean_thickness_type),   intent(in)    :: Thickness
  type(ocean_prog_tracer_type), intent(inout) :: T_prog(:)

  integer :: taum1, tau
  integer :: i, j, k, n

  real    :: sdiff, sum_val, adaptive_coeff
  integer :: secs, days, numsecs
  logical :: do_adaptive_restore = .false.  ! Only restore in specified time period.
  logical :: do_normal_restore = .false.  ! Only restore in specified time period.
  logical, save :: first_pass = .true.

  if(.not. use_this_module) return

  taum1 = Time%taum1
  tau   = Time%tau
  wrk1  = 0.0

  if (use_adaptive_restore) then
    secs_restore = days_to_restore * 86400 + secs_to_restore
    sdiff = 0.0
    call get_time(Time%model_time, secs, days)
    numsecs = (days - initial_day) * 86400 + secs - initial_secs

    do_adaptive_restore = (numsecs < secs_restore)
    do_normal_restore = (.not. do_adaptive_restore) .and. use_basal_after_init
  else
    do_normal_restore = .true.
  endif

  do n=1,size(T_prog(:))

    ! Need to reinitialise wrk2 here due to limiting of temperature.
    wrk2  = 0.0

    if (Basal(n)%id > 0) then

        ! get basal value for current time
        call time_interp_external(Basal(n)%id, Time%model_time, wrk1)
        if (do_adaptive_restore) then
           if (first_pass) then
              if (use_hard_thump) then
                 do k = 1, nk
                    do j = jsd, jed
                       do i = isd, ied
                          T_prog(n)%field(i,j,k,taum1) = wrk1(i,j,k)
                          T_prog(n)%field(i,j,k,tau) = wrk1(i,j,k)
                       end do
                    end do
                 end do
              end if
              wrk2 = abs(T_prog(n)%field(:,:,:,taum1) - wrk1(:,:,:)) &
                     * Grd%tmask(:,:,:)
              sum_val = sum(wrk2(isc:iec,jsc:jec,:))
              call mpp_sum(sum_val)
              sdiffo(n) = sum_val / Grd%wet_t_points
              if (.not. use_normalising) then
                sdiffo(n) = 1.0
              else
                sdiffo(n) = maxval(wrk2)
                call mpp_max(sdiffo(n))
                sdiffo(n) = 1.01 * sdiffo(n)
              endif
           endif

           do k = 1, nk
             do j = jsd, jed
               do i = isd, ied
                 sdiff = abs(T_prog(n)%field(i,j,k,taum1) - wrk1(i,j,k))
                 sdiff = max(sdiff, 1e-9) ! Minimum tolerance
                 adaptive_coeff = 1.0 / (taumin - &
                                         (lambda * real(secs_restore)) &
                                         / (((sdiff / sdiffo(n))**npower) & 
                                            * log(1 - athresh)))
                 wrk2(i,j,k) = Thickness%rho_dzt(i,j,k,tau) * adaptive_coeff &
                                  * (wrk1(i,j,k) - T_prog(n)%field(i,j,k,taum1))
               enddo
             enddo
           enddo
        else if (do_normal_restore) then
            do k = 1, nk
                do j = jsc, jec
                    do i = isc, iec
                        wrk2(i,j,k) = Thickness%rho_dzt(i,j,k,tau) &
                                     * Basal(n)%damp_coeff(i,j,k) &
                                     * (wrk1(i,j,k) - T_prog(n)%field(i,j,k,taum1))
                    end do
                end do
             end do
        end if

    end if

    ! These tracer relaxation schemes were in the original OFAM and AusCOM
    ! source code, and may be needed by other OFAM-based experiments.
    !
    ! For now, let's keep them in the same place, but they are now turned off
    ! by default.  They may need to be moved or deleted in the future.

    ! Limit to -1.8 deg with 3hour restoring (got other limits here for initialising. Not sure if needed (namelist use?).
    if (limit_temp .and. trim(T_prog(n)%name) == 'temp') then
       do k=1,nk
          do j=jsc,jec
             do i=isc,iec
                wrk2(i,j,k) = wrk2(i,j,k) + Thickness%rho_dzt(i,j,k,tau) &
                    * max(limit_temp_min - T_prog(n)%field(i,j,k,taum1), 0.0) &
                    / limit_temp_restore
!                wrk2(i,j,k) = wrk2(i,j,k)+Thickness%rho_dzt(i,j,k,tau)*min(43.0-T_prog(n)%field(i,j,k,taum1),0.0)/3600.0
             enddo
          enddo
       enddo
    end if

    if (limit_salt .and. trim(T_prog(n)%name) == 'salt') then
       do k=1,nk
          do j=jsc,jec
             do i=isc,iec
                wrk2(i,j,k) = wrk2(i,j,k) + Thickness%rho_dzt(i,j,k,tau) &
                    * max(limit_salt_min - T_prog(n)%field(i,j,k,taum1), 0.0) &
                    / limit_salt_restore
             enddo
          enddo
       enddo
    end if

    ! Update tendency
    do k = 1, nk
        do j = jsc, jec
            do i = isc, iec
                T_prog(n)%th_tendency(i,j,k) = T_prog(n)%th_tendency(i,j,k) &
                                              + wrk2(i,j,k)
            end do
        end do
    end do

    if (id_basal_tend(n) > 0) call diagnose_3d(Time, Grd, id_basal_tend(n), &
         T_prog(n)%conversion*wrk2(:,:,:))

  enddo

  first_pass = .false.

  return


end subroutine basal_tracer_source_paul
! </SUBROUTINE> NAME="basal_tracer_source_paul"

end module ocean_basal_tracer_mod
