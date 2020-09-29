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

implicit none

private

#include <ocean_memory.h>

type ocean_basal_type
   integer :: id                                             ! time_interp_external index
   character(len=32) :: name                                 ! tracer name corresponding to basal
   real, dimension(:,:,:), pointer :: damp_coeff   => NULL() ! 3d inverse damping rate (tracer units/ sec)
   real, dimension(:,:)  , pointer :: damp_coeff2d => NULL() ! 2d inverse damping rate (tracer units/ sec)
end type ocean_basal_type


type(ocean_basal_type), allocatable, dimension(:) :: Basal
!Pedro
type(ocean_basal_type), allocatable, dimension(:,:) :: misfkt,misfkb ! Top and bottom input depths
!Pedro
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

logical :: limit_temp            = .false.
real    :: limit_temp_min        = -1.8
real    :: limit_temp_restore    = 10800.

logical :: limit_salt            = .false.
real    :: limit_salt_min        = 0.01
real    :: limit_salt_restore    = 3600.

integer :: secs_restore
integer :: initial_day
integer :: initial_secs
real, allocatable :: sdiffo(:)

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
        PRINT *, "ocean_basal_tracer_init"
  stdoutunit=stdout();stdlogunit=stdlog() 

  if ( module_is_initialized ) then 
    call mpp_error(FATAL, &
    '==>Error in ocean_basal_tracer_mod (ocean_basal_tracer_init): module already initialized')
  endif 

  module_is_initialized = .TRUE.

  !num_prog_tracers = size(T_prog(:))
  num_prog_tracers = 1

  allocate( Basal(num_prog_tracers) )

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

  !Top and bottom Fw input
  misfkt(:,:) = 200
  misfkb(:,:) = 600

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
subroutine basal_tracer_source(Time, Thickness, T_prog)

!Pedro
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
PRINT *, "basal_tracer_source"
!
!  if(.not. use_this_module) return
!
!  taum1 = Time%taum1
!  tau   = Time%tau
!  wrk1  = 0.0
!
!  if (use_adaptive_restore) then
!    secs_restore = days_to_restore * 86400 + secs_to_restore
!    sdiff = 0.0
!    call get_time(Time%model_time, secs, days)
!    numsecs = (days - initial_day) * 86400 + secs - initial_secs
!
!    do_adaptive_restore = (numsecs < secs_restore)
!    do_normal_restore = (.not. do_adaptive_restore) .and. use_basal_after_init
!  else
!    do_normal_restore = .true.
!  endif
!
!  do n=1,size(T_prog(:))
!
!    ! Need to reinitialise wrk2 here due to limiting of temperature.
!    wrk2  = 0.0
!
!    if (Basal(n)%id > 0) then
!
!        ! get basal value for current time
!        call time_interp_external(Basal(n)%id, Time%model_time, wrk1)
!        if (do_adaptive_restore) then
!           if (first_pass) then
!              if (use_hard_thump) then
!                 do k = 1, nk
!                    do j = jsd, jed
!                       do i = isd, ied
!                          T_prog(n)%field(i,j,k,taum1) = wrk1(i,j,k)
!                          T_prog(n)%field(i,j,k,tau) = wrk1(i,j,k)
!                       end do
!                    end do
!                 end do
!              end if
!              wrk2 = abs(T_prog(n)%field(:,:,:,taum1) - wrk1(:,:,:)) &
!                     * Grd%tmask(:,:,:)
!              sum_val = sum(wrk2(isc:iec,jsc:jec,:))
!              call mpp_sum(sum_val)
!              sdiffo(n) = sum_val / Grd%wet_t_points
!              if (.not. use_normalising) then
!                sdiffo(n) = 1.0
!              else
!                sdiffo(n) = maxval(wrk2)
!                call mpp_max(sdiffo(n))
!                sdiffo(n) = 1.01 * sdiffo(n)
!              endif
!           endif
!
!           do k = 1, nk
!             do j = jsd, jed
!               do i = isd, ied
!                 sdiff = abs(T_prog(n)%field(i,j,k,taum1) - wrk1(i,j,k))
!                 sdiff = max(sdiff, 1e-9) ! Minimum tolerance
!                 adaptive_coeff = 1.0 / (taumin - &
!                                         (lambda * real(secs_restore)) &
!                                         / (((sdiff / sdiffo(n))**npower) & 
!                                            * log(1 - athresh)))
!                 wrk2(i,j,k) = Thickness%rho_dzt(i,j,k,tau) * adaptive_coeff &
!                                  * (wrk1(i,j,k) - T_prog(n)%field(i,j,k,taum1))
!               enddo
!             enddo
!           enddo
!        else if (do_normal_restore) then
!            do k = 1, nk
!                do j = jsc, jec
!                    do i = isc, iec
!                        wrk2(i,j,k) = Thickness%rho_dzt(i,j,k,tau) &
!                                     * Basal(n)%damp_coeff(i,j,k) &
!                                     * (wrk1(i,j,k) - T_prog(n)%field(i,j,k,taum1))
!                    end do
!                end do
!             end do
!        end if
!
!    end if
!
!    ! These tracer relaxation schemes were in the original OFAM and AusCOM
!    ! source code, and may be needed by other OFAM-based experiments.
!    !
!    ! For now, let's keep them in the same place, but they are now turned off
!    ! by default.  They may need to be moved or deleted in the future.
!
!    ! Limit to -1.8 deg with 3hour restoring (got other limits here for initialising. Not sure if needed (namelist use?).
!    if (limit_temp .and. trim(T_prog(n)%name) == 'temp') then
!       do k=1,nk
!          do j=jsc,jec
!             do i=isc,iec
!                wrk2(i,j,k) = wrk2(i,j,k) + Thickness%rho_dzt(i,j,k,tau) &
!                    * max(limit_temp_min - T_prog(n)%field(i,j,k,taum1), 0.0) &
!                    / limit_temp_restore
!!                wrk2(i,j,k) = wrk2(i,j,k)+Thickness%rho_dzt(i,j,k,tau)*min(43.0-T_prog(n)%field(i,j,k,taum1),0.0)/3600.0
!             enddo
!          enddo
!       enddo
!    end if
!
!    if (limit_salt .and. trim(T_prog(n)%name) == 'salt') then
!       do k=1,nk
!          do j=jsc,jec
!             do i=isc,iec
!                wrk2(i,j,k) = wrk2(i,j,k) + Thickness%rho_dzt(i,j,k,tau) &
!                    * max(limit_salt_min - T_prog(n)%field(i,j,k,taum1), 0.0) &
!                    / limit_salt_restore
!             enddo
!          enddo
!       enddo
!    end if
!
!    ! Update tendency
!    do k = 1, nk
!        do j = jsc, jec
!            do i = isc, iec
!                T_prog(n)%th_tendency(i,j,k) = T_prog(n)%th_tendency(i,j,k) &
!                                              + wrk2(i,j,k)
!            end do
!        end do
!    end do
!
!    if (id_basal_tend(n) > 0) call diagnose_3d(Time, Grd, id_basal_tend(n), &
!         T_prog(n)%conversion*wrk2(:,:,:))
!
!  enddo
!
!  first_pass = .false.
!
!  return


end subroutine basal_tracer_source
! </SUBROUTINE> NAME="basal_tracer_source"

  SUBROUTINE sbc_isf_bg03(kt)
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE sbc_isf_bg03  ***
      !!
      !! ** Purpose : add net heat and fresh water flux from ice shelf melting
      !!          into the adjacent ocean
      !!
      !! ** Method  :   See reference, just ISOMIP
      !!
      !! ** Reference : Beckmann and Goosse (2003), "A parameterization of ice shelf-ocean
      !!         interaction for climate models", Ocean Modelling 5(2003) 157-170.
      !!         (hereafter BG)
      !! History :  06-02  (C. Wang) Original code
      !!----------------------------------------------------------------------
      INTEGER, INTENT ( in ) :: kt
      !
      INTEGER  :: i, j, k ! dummy loop index
      INTEGER  :: ik         ! current level
      real :: zt_sum     ! sum of the temperature between 200m and 600m
      real :: zt_ave     ! averaged temperature between 200m and 600m
      real :: zt_frz     ! freezing point temperature at depth z
      real :: zpress     ! pressure to compute the freezing point in depth
      real :: grav = 9.81 !gravity value
      !!----------------------------------------------------------------------
      !
      DO i = isc, iec
         DO j = jsc, jec
            ik = misfkt(ji,jj)
            !! Initialize arrays to 0 (each step)
            zt_sum = 0.e0_wp
            IF ( ik > 1 ) THEN
               ! 1. -----------the average temperature between 200m and 600m ---------------------
               DO jk = misfkt(ji,jj),misfkb(ji,jj)
                  ! Calculate freezing temperature
                  zpress = grav*rau0*gdept_n(ji,jj,ik)*1.e-04
                  CALL eos_fzp(stbl(ji,jj), zt_frz, zpress)
                  zt_sum = zt_sum + (tsn(ji,jj,jk,jp_tem)-zt_frz) * e3t_n(ji,jj,jk) * tmask(ji,jj,jk)  ! sum temp
               END DO
               zt_ave = zt_sum/rhisf_tbl(ji,jj) ! calcul mean value
               ! 2. ------------Net heat flux and fresh water flux due to the ice shelf
               ! For those corresponding to zonal boundary
               qisf(ji,jj) = - rau0 * rcp * rn_gammat0 * risfLeff(ji,jj) * e1t(ji,jj) * zt_ave  &
                           & * r1_e1e2t(ji,jj) * tmask(ji,jj,jk)

               fwfisf(ji,jj) = qisf(ji,jj) / rLfusisf          !fresh water flux kg/(m2s)
               fwfisf(ji,jj) = fwfisf(ji,jj) * ( soce / stbl(ji,jj) )
               !add to salinity trend
            ELSE
               qisf(ji,jj) = 0._wp   ;   fwfisf(ji,jj) = 0._wp
            END IF
         END DO
      END DO
      !
  END SUBROUTINE sbc_isf_bg03

  SUBROUTINE sbc_isf_cav( kt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE sbc_isf_cav  ***
      !!
      !! ** Purpose :   handle surface boundary condition under ice shelf
      !!
      !! ** Method  : Mathiot et al., 2014 with two options of eq.
      !!
      !! ** Action  :   utau, vtau : remain unchanged
      !!                taum, wndm : remain unchanged
      !!                qns        : update heat flux below ice shelf
      !!                emp, emps  : update freshwater flux below ice shelf
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! ocean time step
      !
      INTEGER  ::   ji, jj     ! dummy loop indices
      INTEGER  ::   nit
      LOGICAL  ::   lit
      REAL(wp) ::   zlamb1, zlamb2, zlamb3
      REAL(wp) ::   zeps1,zeps2,zeps3,zeps4,zeps6,zeps7
      REAL(wp) ::   zaqe,zbqe,zcqe,zaqer,zdis,zsfrz,zcfac
      REAL(wp) ::   zeps = 1.e-20_wp
      REAL(wp) ::   zerr
      REAL(wp), DIMENSION(jpi,jpj) ::   zfrz
      REAL(wp), DIMENSION(jpi,jpj) ::   zgammat, zgammas
      REAL(wp), DIMENSION(jpi,jpj) ::   zfwflx, zhtflx, zhtflx_b
      !!---------------------------------------------------------------------
      !
      ! coeficient for linearisation of potential tfreez
      ! Crude approximation for pressure (but commonly used)
      IF ( l_useCT ) THEN   ! linearisation from Jourdain et al. (2017)
         zlamb1 =-0.0564_wp
         zlamb2 = 0.0773_wp
         zlamb3 =-7.8633e-8 * grav * rau0
      ELSE                  ! linearisation from table 4 (Asay-Davis et al., 2015)
         zlamb1 =-0.0573_wp
         zlamb2 = 0.0832_wp
         zlamb3 =-7.53e-8 * grav * rau0
      ENDIF
      !
      ! initialisation
      zgammat(:,:) = rn_gammat0 ; zgammas (:,:) = rn_gammas0
      zhtflx (:,:) = 0.0_wp     ; zhtflx_b(:,:) = 0.0_wp
      zfwflx (:,:) = 0.0_wp

      ! compute ice shelf melting
      nit = 1 ; lit = .TRUE.
      DO WHILE ( lit )    ! maybe just a constant number of iteration as in blk_core is fine
         SELECT CASE ( nn_isfblk )
         CASE ( 1 )   !  ISOMIP formulation (2 equations) for volume flux (Hunter et al., 2006)
            ! Calculate freezing temperature
            CALL eos_fzp( stbl(:,:), zfrz(:,:), risfdep(:,:) )

            ! compute gammat every where (2d)
            CALL sbc_isf_gammats(zgammat, zgammas, zhtflx, zfwflx)

            ! compute upward heat flux zhtflx and upward water flux zwflx
            DO jj = 1, jpj
               DO ji = 1, jpi
                  zhtflx(ji,jj) =   zgammat(ji,jj)*rcp*rau0*(ttbl(ji,jj)-zfrz(ji,jj))
                  zfwflx(ji,jj) = - zhtflx(ji,jj)/rLfusisf
               END DO
            END DO

            ! Compute heat flux and upward fresh water flux
            qisf  (:,:) = - zhtflx(:,:) * (1._wp - tmask(:,:,1)) * ssmask(:,:)
            fwfisf(:,:) =   zfwflx(:,:) * (1._wp - tmask(:,:,1)) * ssmask(:,:)
         CASE ( 2 )  ! ISOMIP+ formulation (3 equations) for volume flux (Asay-Davis et al., 2015)
            ! compute gammat every where (2d)
            CALL sbc_isf_gammats(zgammat, zgammas, zhtflx, zfwflx)

            ! compute upward heat flux zhtflx and upward water flux zwflx
            ! Resolution of a 2d equation from equation 21, 22 and 23 to find Sb (Asay-Davis et al., 2015)
            DO jj = 1, jpj
               DO ji = 1, jpi
                  ! compute coeficient to solve the 2nd order equation
                  zeps1 = rcp*rau0*zgammat(ji,jj)
                  zeps2 = rLfusisf*rau0*zgammas(ji,jj)
                  zeps3 = rhoisf*rcpisf*rkappa/MAX(risfdep(ji,jj),zeps)
                  zeps4 = zlamb2+zlamb3*risfdep(ji,jj)
                  zeps6 = zeps4-ttbl(ji,jj)
                  zeps7 = zeps4-tsurf
                  zaqe  = zlamb1 * (zeps1 + zeps3)
                  zaqer = 0.5_wp/MIN(zaqe,-zeps)
                  zbqe  = zeps1*zeps6+zeps3*zeps7-zeps2
                  zcqe  = zeps2*stbl(ji,jj)
                  zdis  = zbqe*zbqe-4.0_wp*zaqe*zcqe

                  ! Presumably zdis can never be negative because gammas is very small compared to gammat
                  ! compute s freeze
                  zsfrz=(-zbqe-SQRT(zdis))*zaqer
                  IF ( zsfrz < 0.0_wp ) zsfrz=(-zbqe+SQRT(zdis))*zaqer

                  ! compute t freeze (eq. 22)
                  zfrz(ji,jj)=zeps4+zlamb1*zsfrz

                  ! zfwflx is upward water flux
                  ! zhtflx is upward heat flux (out of ocean)
                  ! compute the upward water and heat flux (eq. 28 and eq. 29)
                  zfwflx(ji,jj) = rau0 * zgammas(ji,jj) * (zsfrz-stbl(ji,jj)) / MAX(zsfrz,zeps)
                  zhtflx(ji,jj) = zgammat(ji,jj) * rau0 * rcp * (ttbl(ji,jj) - zfrz(ji,jj) )
               END DO
            END DO

            ! compute heat and water flux
            qisf  (:,:) = - zhtflx(:,:) * (1._wp - tmask(:,:,1)) * ssmask(:,:)
            fwfisf(:,:) =   zfwflx(:,:) * (1._wp - tmask(:,:,1)) * ssmask(:,:)

         END SELECT

         ! define if we need to iterate (nn_gammablk 0/1 do not need iteration)
         IF ( nn_gammablk <  2 ) THEN ; lit = .FALSE.
         ELSE
            ! check total number of iteration
            IF (nit >= 100) THEN ; CALL ctl_stop( 'STOP', 'sbc_isf_hol99 : too many iteration ...' )
            ELSE                 ; nit = nit + 1
            END IF

            ! compute error between 2 iterations
            ! if needed save gammat and compute zhtflx_b for next iteration
            zerr = MAXVAL(ABS(zhtflx-zhtflx_b))
            IF ( zerr <= 0.01_wp ) THEN ; lit = .FALSE.
            ELSE                        ; zhtflx_b(:,:) = zhtflx(:,:)
            END IF
         END IF
      END DO
      !
      CALL iom_put('isfgammat', zgammat)
      CALL iom_put('isfgammas', zgammas)
      ! 
   END SUBROUTINE sbc_isf_cav

   SUBROUTINE sbc_isf_gammats(pgt, pgs, pqhisf, pqwisf )
      !!----------------------------------------------------------------------
      !! ** Purpose    : compute the coefficient echange for heat flux
      !!
      !! ** Method     : gamma assume constant or depends of u* and stability
      !!
      !! ** References : Holland and Jenkins, 1999, JPO, p1787-1800, eq 14
      !!                Jenkins et al., 2010, JPO, p2298-2312
      !!---------------------------------------------------------------------
      REAL(wp), DIMENSION(:,:), INTENT(  out) ::   pgt   , pgs      ! 
      REAL(wp), DIMENSION(:,:), INTENT(in   ) ::   pqhisf, pqwisf   ! 
      !
      INTEGER  :: ji, jj                     ! loop index
      INTEGER  :: ikt                        ! local integer
      REAL(wp) :: zdku, zdkv                 ! U, V shear 
      REAL(wp) :: zPr, zSc, zRc              ! Prandtl, Scmidth and Richardson number 
      REAL(wp) :: zmob, zmols                ! Monin Obukov length, coriolis factor at T point
      REAL(wp) :: zbuofdep, zhnu             ! Bouyancy length scale, sublayer tickness
      REAL(wp) :: zhmax                      ! limitation of mol
      REAL(wp) :: zetastar                   ! stability parameter
      REAL(wp) :: zgmolet, zgmoles, zgturb   ! contribution of modelecular sublayer and turbulence 
      REAL(wp) :: zcoef                      ! temporary coef
      REAL(wp) :: zdep
      REAL(wp) :: zeps = 1.0e-20_wp
      REAL(wp), PARAMETER :: zxsiN = 0.052_wp   ! dimensionless constant
      REAL(wp), PARAMETER :: znu   = 1.95e-6_wp ! kinamatic viscosity of sea water (m2.s-1)
      REAL(wp), DIMENSION(2) :: zts, zab
      REAL(wp), DIMENSION(jpi,jpj) :: zustar   ! U, V at T point and friction velocity
      !!---------------------------------------------------------------------
      !
      SELECT CASE ( nn_gammablk )
      CASE ( 0 ) ! gamma is constant (specified in namelist)
         !! ISOMIP formulation (Hunter et al, 2006)
         pgt(:,:) = rn_gammat0
         pgs(:,:) = rn_gammas0

      CASE ( 1 ) ! gamma is assume to be proportional to u*
         !! Jenkins et al., 2010, JPO, p2298-2312
         !! Adopted by Asay-Davis et al. (2015)
         !! compute ustar (eq. 24)
!!gm  NB  use pCdU here so that it will incorporate local boost of Cd0 and log layer case :
!!         zustar(:,:) = SQRT( rCdU_top(:,:) * SQRT(utbl(:,:) * utbl(:,:) + vtbl(:,:) * vtbl(:,:) + r_ke0_top) )
!! or better :  compute ustar in zdfdrg  and use it here as well as in TKE, GLS and Co
!!
!!     ===>>>>    GM  to be done this chrismas
!!         
!!gm end
         zustar(:,:) = SQRT( r_Cdmin_top * (utbl(:,:) * utbl(:,:) + vtbl(:,:) * vtbl(:,:) + r_ke0_top) )

         !! Compute gammats
         pgt(:,:) = zustar(:,:) * rn_gammat0
         pgs(:,:) = zustar(:,:) * rn_gammas0

      CASE ( 2 ) ! gamma depends of stability of boundary layer
         !! Holland and Jenkins, 1999, JPO, p1787-1800, eq 14
         !! as MOL depends of flux and flux depends of MOL, best will be iteration (TO DO)
         !! compute ustar
         zustar(:,:) = SQRT( r_Cdmin_top * (utbl(:,:) * utbl(:,:) + vtbl(:,:) * vtbl(:,:) + r_ke0_top) )

         !! compute Pr and Sc number (can be improved)
         zPr =   13.8_wp
         zSc = 2432.0_wp

         !! compute gamma mole
         zgmolet = 12.5_wp * zPr ** (2.0/3.0) - 6.0_wp
         zgmoles = 12.5_wp * zSc ** (2.0/3.0) - 6.0_wp

         !! compute gamma
         DO ji = 2, jpi
            DO jj = 2, jpj
               ikt = mikt(ji,jj)

               IF( zustar(ji,jj) == 0._wp ) THEN           ! only for kt = 1 I think
                  pgt = rn_gammat0
                  pgs = rn_gammas0
               ELSE
                  !! compute Rc number (as done in zdfric.F90)
!!gm better to do it like in the new zdfric.F90   i.e. avm weighted Ri computation
!!gm moreover, use Max(rn2,0) to take care of static instabilities....
                  zcoef = 0.5_wp / e3w_n(ji,jj,ikt+1)
                  !                                            ! shear of horizontal velocity
                  zdku = zcoef * (  un(ji-1,jj  ,ikt  ) + un(ji,jj,ikt  )  &
                     &             -un(ji-1,jj  ,ikt+1) - un(ji,jj,ikt+1)  )
                  zdkv = zcoef * (  vn(ji  ,jj-1,ikt  ) + vn(ji,jj,ikt  )  &
                     &             -vn(ji  ,jj-1,ikt+1) - vn(ji,jj,ikt+1)  )
                  !                                            ! richardson number (minimum value set to zero)
                  zRc = rn2(ji,jj,ikt+1) / MAX( zdku*zdku + zdkv*zdkv, zeps )

                  !! compute bouyancy 
                  zts(jp_tem) = ttbl(ji,jj)
                  zts(jp_sal) = stbl(ji,jj)
                  zdep        = gdepw_n(ji,jj,ikt)
                  !
                  CALL eos_rab( zts, zdep, zab )
                  !
                  !! compute length scale 
                  zbuofdep = grav * ( zab(jp_tem) * pqhisf(ji,jj) - zab(jp_sal) * pqwisf(ji,jj) )  !!!!!!!!!!!!!!!!!!!!!!!!!!!!

                  !! compute Monin Obukov Length
                  ! Maximum boundary layer depth
                  zhmax = gdept_n(ji,jj,mbkt(ji,jj)) - gdepw_n(ji,jj,mikt(ji,jj)) -0.001_wp
                  ! Compute Monin obukhov length scale at the surface and Ekman depth:
                  zmob   = zustar(ji,jj) ** 3 / (vkarmn * (zbuofdep + zeps))
                  zmols  = SIGN(1._wp, zmob) * MIN(ABS(zmob), zhmax) * tmask(ji,jj,ikt)

                  !! compute eta* (stability parameter)
                  zetastar = 1._wp / ( SQRT(1._wp + MAX(zxsiN * zustar(ji,jj) / ( ABS(ff_f(ji,jj)) * zmols * zRc ), 0._wp)))

                  !! compute the sublayer thickness
                  zhnu = 5 * znu / zustar(ji,jj)

                  !! compute gamma turb
                  zgturb = 1._wp / vkarmn * LOG(zustar(ji,jj) * zxsiN * zetastar * zetastar / ( ABS(ff_f(ji,jj)) * zhnu )) &
                  &      + 1._wp / ( 2 * zxsiN * zetastar ) - 1._wp / vkarmn

                  !! compute gammats
                  pgt(ji,jj) = zustar(ji,jj) / (zgturb + zgmolet)
                  pgs(ji,jj) = zustar(ji,jj) / (zgturb + zgmoles)
               END IF
            END DO
         END DO
         CALL lbc_lnk_multi( 'sbcisf', pgt, 'T', 1., pgs, 'T', 1.)
      END SELECT
      !
   END SUBROUTINE sbc_isf_gammats

   SUBROUTINE sbc_isf_tbl( pvarin, pvarout, cd_ptin )
      !!----------------------------------------------------------------------
      !!                  ***  SUBROUTINE sbc_isf_tbl  ***
      !!
      !! ** Purpose : compute mean T/S/U/V in the boundary layer at T- point
      !!
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(:,:,:), INTENT(in   ) :: pvarin
      REAL(wp), DIMENSION(:,:)  , INTENT(  out) :: pvarout
      CHARACTER(len=1),           INTENT(in   ) :: cd_ptin ! point of variable in/out
      !
      INTEGER ::   ji, jj, jk                   ! loop index
      INTEGER ::   ikt, ikb                     ! top and bottom index of the tbl
      REAL(wp) ::   ze3, zhk
      REAL(wp), DIMENSION(jpi,jpj) :: zhisf_tbl ! thickness of the tbl
      REAL(wp), DIMENSION(jpi,jpj) :: zvarout
      !!----------------------------------------------------------------------

      ! initialisation
      pvarout(:,:)=0._wp

      SELECT CASE ( cd_ptin )
      CASE ( 'U' ) ! compute U in the top boundary layer at T- point 
         !
         zvarout(:,:)=0._wp
         !
         DO jj = 1,jpj
            DO ji = 1,jpi
               ikt = miku(ji,jj) ; ikb = miku(ji,jj)
               ! thickness of boundary layer at least the top level thickness
               zhisf_tbl(ji,jj) = MAX( rhisf_tbl_0(ji,jj) , e3u_n(ji,jj,ikt) )

               ! determine the deepest level influenced by the boundary layer
               DO jk = ikt+1, mbku(ji,jj)
                  IF ( (SUM(e3u_n(ji,jj,ikt:jk-1)) < zhisf_tbl(ji,jj)) .AND. (umask(ji,jj,jk) == 1) ) ikb = jk
                  IF ( (SUM(e3u_n(ji,jj,ikt:jk-1)) < zhisf_tbl(ji,jj)) .AND. (umask(ji,jj,jk) == 1) ) ikb = jk
               END DO
               zhisf_tbl(ji,jj) = MIN(zhisf_tbl(ji,jj), SUM(e3u_n(ji,jj,ikt:ikb)))  ! limit the tbl to water thickness.

               ! level fully include in the ice shelf boundary layer
               DO jk = ikt, ikb - 1
                  ze3 = e3u_n(ji,jj,jk)
                  zvarout(ji,jj) = zvarout(ji,jj) + pvarin(ji,jj,jk) / zhisf_tbl(ji,jj) * ze3
               END DO

               ! level partially include in ice shelf boundary layer 
               zhk = SUM( e3u_n(ji, jj, ikt:ikb - 1)) / zhisf_tbl(ji,jj)
               zvarout(ji,jj) = zvarout(ji,jj) + pvarin(ji,jj,ikb) * (1._wp - zhk)
            END DO
         END DO
         DO jj = 2, jpj
            DO ji = 2, jpi
!!gm a wet-point only average should be used here !!!
               pvarout(ji,jj) = 0.5_wp * (zvarout(ji,jj) + zvarout(ji-1,jj))
            END DO
         END DO
         CALL lbc_lnk('sbcisf', pvarout,'T',-1.)

      CASE ( 'V' ) ! compute V in the top boundary layer at T- point 
         !
         zvarout(:,:)=0._wp
         !
         DO jj = 1,jpj
            DO ji = 1,jpi
               ikt = mikv(ji,jj) ; ikb = mikv(ji,jj)
               ! thickness of boundary layer at least the top level thickness
               zhisf_tbl(ji,jj) = MAX(rhisf_tbl_0(ji,jj), e3v_n(ji,jj,ikt))

               ! determine the deepest level influenced by the boundary layer
               DO jk = ikt+1, mbkv(ji,jj)
                  IF ( (SUM(e3v_n(ji,jj,ikt:jk-1)) < zhisf_tbl(ji,jj)) .AND. (vmask(ji,jj,jk) == 1) ) ikb = jk
                  IF ( (SUM(e3u_n(ji,jj,ikt:jk-1)) < zhisf_tbl(ji,jj)) .AND. (umask(ji,jj,jk) == 1) ) ikb = jk
               END DO
               zhisf_tbl(ji,jj) = MIN(zhisf_tbl(ji,jj), SUM(e3u_n(ji,jj,ikt:ikb)))  ! limit the tbl to water thickness.

               ! level fully include in the ice shelf boundary layer
               DO jk = ikt, ikb - 1
                  ze3 = e3u_n(ji,jj,jk)
                  zvarout(ji,jj) = zvarout(ji,jj) + pvarin(ji,jj,jk) / zhisf_tbl(ji,jj) * ze3
               END DO

               ! level partially include in ice shelf boundary layer 
               zhk = SUM( e3u_n(ji, jj, ikt:ikb - 1)) / zhisf_tbl(ji,jj)
               zvarout(ji,jj) = zvarout(ji,jj) + pvarin(ji,jj,ikb) * (1._wp - zhk)
            END DO
         END DO
         DO jj = 2, jpj
            DO ji = 2, jpi
!!gm a wet-point only average should be used here !!!
               pvarout(ji,jj) = 0.5_wp * (zvarout(ji,jj) + zvarout(ji-1,jj))
            END DO
         END DO
         CALL lbc_lnk('sbcisf', pvarout,'T',-1.)

      CASE ( 'V' ) ! compute V in the top boundary layer at T- point 
         !
         zvarout(:,:)=0._wp
         !
         DO jj = 1,jpj
            DO ji = 1,jpi
               ikt = mikv(ji,jj) ; ikb = mikv(ji,jj)
               ! thickness of boundary layer at least the top level thickness
               zhisf_tbl(ji,jj) = MAX(rhisf_tbl_0(ji,jj), e3v_n(ji,jj,ikt))

               ! determine the deepest level influenced by the boundary layer
               DO jk = ikt+1, mbkv(ji,jj)
                  IF ( (SUM(e3v_n(ji,jj,ikt:jk-1)) < zhisf_tbl(ji,jj)) .AND. (vmask(ji,jj,jk) == 1) ) ikb = jk
               END DO
               zhisf_tbl(ji,jj) = MIN(zhisf_tbl(ji,jj), SUM(e3v_n(ji,jj,ikt:ikb)))  ! limit the tbl to water thickness.

               ! level fully include in the ice shelf boundary layer
               DO jk = ikt, ikb - 1
                  ze3 = e3v_n(ji,jj,jk)
                  zvarout(ji,jj) = zvarout(ji,jj) + pvarin(ji,jj,jk) / zhisf_tbl(ji,jj) * ze3
               END DO

               ! level partially include in ice shelf boundary layer 
               zhk = SUM( e3v_n(ji, jj, ikt:ikb - 1)) / zhisf_tbl(ji,jj)
               zvarout(ji,jj) = zvarout(ji,jj) + pvarin(ji,jj,ikb) * (1._wp - zhk)
            END DO
         END DO
         DO jj = 2, jpj
            DO ji = 2, jpi
!!gm a wet-point only average should be used here !!!
               pvarout(ji,jj) = 0.5_wp * (zvarout(ji,jj) + zvarout(ji,jj-1))
            END DO
         END DO
         CALL lbc_lnk('sbcisf', pvarout,'T',-1.)

      CASE ( 'T' ) ! compute T in the top boundary layer at T- point 
         DO jj = 1,jpj
            DO ji = 1,jpi
               ikt = misfkt(ji,jj)
               ikb = misfkb(ji,jj)

               ! level fully include in the ice shelf boundary layer
               DO jk = ikt, ikb - 1
                  ze3 = e3t_n(ji,jj,jk)
                  pvarout(ji,jj) = pvarout(ji,jj) + pvarin(ji,jj,jk) * r1_hisf_tbl(ji,jj) * ze3
               END DO

               ! level partially include in ice shelf boundary layer 
               zhk = SUM( e3t_n(ji, jj, ikt:ikb - 1)) * r1_hisf_tbl(ji,jj)
              pvarout(ji,jj) = pvarout(ji,jj) + pvarin(ji,jj,ikb) * (1._wp - zhk)
            END DO
         END DO
      END SELECT
      !
      ! mask mean tbl value
      pvarout(:,:) = pvarout(:,:) * ssmask(:,:)
      !
   END SUBROUTINE sbc_isf_tbl


end module ocean_basal_tracer_mod
