

module sas_conv_intr
!---------------------------------------------------------------------------------
! Purpose:
!
! CAM interface to the Simplified Arakawa Schubert deep convection scheme
!
! Author: Yi-Chi Wang
! Nov 25, 2013 - output heo and heos for determining kbconv and kb
! Nov 14, 2013 - add output for SAS triggering mask flag
! Oct 12, 2013 - remove zm evap
! Sept, 2013 - no convective transport
! Nov 8, 2012  - channel the flag of activated deep cumulus into shallow.
! Oct 1, 2012  - remove the output of SASQCKO into netcdf file
!                configuration with 60 levels would cause problem.
! July 5, 2012 - remove the utility of cloud liquid/ice partition
!                    and direct detrainment of cloud water from deep cps
!                cloudice = .false.
! Apr 23, 2012 - change the time step for sas_conv as 0.5*ztodt (=delt)
!                however, the tendency is still for ztodt (=2*delt) <= reverse
! Apr 1,2 2012 - add subroutine for convective transport
! Mar 14, 2012 - add tendency to liquid water
! Mar 7,  2012 - add new output arguments to subroutine sascnvn
! Feb 21, 2012 - to the 
!      Q: whether to include a module for the format and physical constants of convection
! January 2010 modified by J. Kay to add COSP simulator fields to physics buffer
!---------------------------------------------------------------------------------
! Yi-Chi : need to check the physconst and select the constant
   use shr_kind_mod, only: r8=>shr_kind_r8
   use physconst,    only: cpair,latvap 
   use ppgrid,       only: pver, pcols, pverp, begchunk, endchunk
   use zm_conv,      only: zm_conv_evap
   use sas_conv,     only: sascnvn, convtran_sas
   use cam_history,  only: outfld, addfld, add_default, phys_decomp
   use perf_mod
   use cam_logfile,  only: iulog
   !use phys_buffer,  only: pbuf_add
   implicit none
   private
   save

   ! Public methods

   public ::&
      sas_conv_register,           &! register fields in physics buffer
      sas_conv_init,               &! initialize donner_deep module
      sas_conv_tend,               &! return tendencies
      sas_conv_tend_2!,               &! return tendencies

   ! Private module data
   ! Yi-Chi: Apr 2012 : add downdraft detrainment
   real(r8), allocatable, dimension(:,:,:) :: dd  !(pcols,pver,begchunk:endchunk)

   real(r8), allocatable, dimension(:,:,:) :: mu  !(pcols,pver,begchunk:endchunk)
   real(r8), allocatable, dimension(:,:,:) :: eu  !(pcols,pver,begchunk:endchunk)
   real(r8), allocatable, dimension(:,:,:) :: du  !(pcols,pver,begchunk:endchunk)
   real(r8), allocatable, dimension(:,:,:) :: md  !(pcols,pver,begchunk:endchunk)
   real(r8), allocatable, dimension(:,:,:) :: ed  !(pcols,pver,begchunk:endchunk)
   real(r8), allocatable, dimension(:,:,:) :: dp  !(pcols,pver,begchunk:endchunk) 
	! wg layer thickness in mbs (between upper/lower interface).
   real(r8), allocatable, dimension(:,:)   :: dsubcld  !(pcols,begchunk:endchunk)
	! wg layer thickness in mbs between lcl and maxi.

   integer, allocatable, dimension(:,:) :: jt   !(pcols,begchunk:endchunk)
        ! wg top  level index of deep cumulus convection.
   integer, allocatable, dimension(:,:) :: maxg !(pcols,begchunk:endchunk)
        ! wg gathered values of maxi.
   integer, allocatable, dimension(:,:) :: ideep !(pcols,begchunk:endchunk)               
	! w holds position of gathered points vs longitude index

   integer, allocatable, dimension(:) :: lengath !(begchunk:endchunk)

   integer ::& ! indices for fields in the physics buffer
      dp_flxprc_idx, &
      dp_flxsnw_idx, &
      dp_cldliq_idx, &
      dp_cldice_idx, &
      prec_dp_idx,   &
      snow_dp_idx


!  indices for fields in the physics buffer
   integer  ::    cld_idx          = 0    
   integer  ::    icwmrdp_idx      = 0     
   integer  ::    rprddp_idx       = 0    
   integer  ::    fracis_idx       = 0   
   integer  ::    nevapr_dpcu_idx  = 0    

!=========================================================================================
contains
!=========================================================================================

subroutine sas_conv_register

!----------------------------------------
! Purpose: register fields with the physics buffer
!----------------------------------------

  use physics_buffer, only : pbuf_add_field, dtype_r8

  implicit none

  integer idx

! Flux of precipitation from deep convection (kg/m2/s)
   call pbuf_add_field('DP_FLXPRC','global',dtype_r8,(/pcols,pverp/),dp_flxprc_idx)

! Flux of snow from deep convection (kg/m2/s)
   call pbuf_add_field('DP_FLXSNW','global',dtype_r8,(/pcols,pverp/),dp_flxsnw_idx)

! deep gbm cloud liquid water (kg/kg)
   call pbuf_add_field('DP_CLDLIQ','global',dtype_r8,(/pcols,pver/), dp_cldliq_idx)

! deep gbm cloud liquid water (kg/kg)
   call pbuf_add_field('DP_CLDICE','global',dtype_r8,(/pcols,pver/), dp_cldice_idx)

end subroutine sas_conv_register

!=========================================================================================

subroutine sas_conv_init(pref_edge)

!----------------------------------------
! Purpose:  declare output fields, initialize variables needed by convection
!----------------------------------------

  use cam_history,    only: outfld, addfld, add_default, phys_decomp
  use ppgrid,         only: pcols, pver
  use zm_conv,        only: zm_convi
  use pmgrid,         only: plev,plevp
  use spmd_utils,     only: masterproc
  use error_messages, only: alloc_err	
  use phys_control,   only: phys_deepconv_pbl, phys_getopts, cam_physpkg_is
  use physics_buffer, only: pbuf_get_index

  implicit none

  real(r8),intent(in) :: pref_edge(plevp)        ! reference pressures at interfaces

  logical :: no_deep_pbl    ! if true, no deep convection in PBL
  integer  limcnv           ! top interface level limit for convection
  integer k, istat
  logical :: history_budget ! output tendencies and state variables for CAM4
                            ! temperature, water vapor, cloud ice and cloud
                            ! liquid budgets.
  integer :: history_budget_histfile_num ! output history file number for budget fields

!
! Allocate space for arrays private to this module
!
     allocate( mu(pcols,pver,begchunk:endchunk), stat=istat )
      call alloc_err( istat, 'sas_conv_tend', 'mu', &
                      pcols*pver*((endchunk-begchunk)+1) )
     allocate( eu(pcols,pver,begchunk:endchunk), stat=istat )
      call alloc_err( istat, 'sas_conv_tend', 'eu', &
                      pcols*pver*((endchunk-begchunk)+1) )
     allocate( du(pcols,pver,begchunk:endchunk), stat=istat )
      call alloc_err( istat, 'sas_conv_tend', 'du', &
                      pcols*pver*((endchunk-begchunk)+1) )
     allocate( md(pcols,pver,begchunk:endchunk), stat=istat )
      call alloc_err( istat, 'sas_conv_tend', 'md', &
                      pcols*pver*((endchunk-begchunk)+1) )
     allocate( ed(pcols,pver,begchunk:endchunk), stat=istat )
      call alloc_err( istat, 'sas_conv_tend', 'ed', &
                      pcols*pver*((endchunk-begchunk)+1) )
!------ Yi-Chi:Apr 1, 2012: for downdraft detrainment in SAS
     allocate( dd(pcols,pver,begchunk:endchunk), stat=istat )
      call alloc_err( istat, 'sas_conv_tend', 'dd', &
                      pcols*pver*((endchunk-begchunk)+1) )
!------
     allocate( dp(pcols,pver,begchunk:endchunk), stat=istat )
      call alloc_err( istat, 'sas_conv_tend', 'dp', &
                      pcols*pver*((endchunk-begchunk)+1) )
     allocate( dsubcld(pcols,begchunk:endchunk), stat=istat )
      call alloc_err( istat, 'sas_conv_tend', 'dsubcld', &
                      pcols*((endchunk-begchunk)+1) )
     allocate( jt(pcols,begchunk:endchunk), stat=istat )
      call alloc_err( istat, 'sas_conv_tend', 'jt', &
                      pcols*((endchunk-begchunk)+1) )
     allocate( maxg(pcols,begchunk:endchunk), stat=istat )
      call alloc_err( istat, 'sas_conv_tend', 'maxg', &
                      pcols*((endchunk-begchunk)+1) )
     allocate( ideep(pcols,begchunk:endchunk), stat=istat )
      call alloc_err( istat, 'sas_conv_tend', 'ideep', &
                      pcols*((endchunk-begchunk)+1) )
     allocate( lengath(begchunk:endchunk), stat=istat )
      call alloc_err( istat, 'sas_conv_tend', 'lengath', &
                      ((endchunk-begchunk)+1) )


! 
! Register fields with the output buffer
!


    call addfld ('PRECZ   ','m/s     ',1,    'A','total precipitation from ZM convection',        phys_decomp)
    call addfld ('SASDT    ','K/s     ',pver, 'A','T tendency - SAS moist convection', phys_decomp)
    call addfld ('SASDQ    ','kg/kg/s ',pver, 'A','Q tendency - SAS moist convection', phys_decomp)
    call addfld ('SASDICE ','kg/kg/s ',pver, 'A','Cloud ice tendency - SAS convection',phys_decomp)
    call addfld ('SASDLIQ ','kg/kg/s ',pver, 'A','Cloud liq tendency - SAS convection',phys_decomp)
    
    call addfld ('EVAPTSAS ','K/s     ',pver, 'A','T tendency - Evaporation/snow prod from SAS convection',phys_decomp)
    call addfld ('FZSNTSAS ','K/s     ',pver, 'A','T tendency - Rain to snow conversion from SAS convection',phys_decomp)
    call addfld ('EVSNTSAS ','K/s     ',pver, 'A','T tendency - Snow to rain prod from SAS convection',phys_decomp)
    call addfld ('EVAPQSAS ','kg/kg/s ',pver, 'A','Q tendency - Evaporation from SAS moist convection',phys_decomp)
    call add_default('EVAPTSAS  ', 1, ' ')
    call add_default('EVAPQSAS  ', 1, ' ')    

    call addfld ('SASFLXPRC','kg/m2/s ',pverp, 'A','Flux of precipitation from SAS convection'       ,phys_decomp)
    call addfld ('SASFLXSNW','kg/m2/s ',pverp, 'A','Flux of snow from SAS convection'                ,phys_decomp)
    call addfld ('SASNTPRPD','kg/kg/s ',pver , 'A','Net precipitation production from SAS convection',phys_decomp)
    call addfld ('SASNTSNPD','kg/kg/s ',pver , 'A','Net snow production from SAS convection'         ,phys_decomp)
    call addfld ('SASEIHEAT','W/kg'    ,pver , 'A','Heating by ice and evaporation in SAS convection',phys_decomp)
    
    call addfld ('CMFMCDSAS','kg/m2/s ',pverp,'A','Convection mass flux from SAS deep ',phys_decomp)
    call addfld ('PRECCDSAS','m/s     ',1,    'A','Convective precipitation rate from SAS deep',phys_decomp)
    call add_default ('CMFMCDSAS', 1, ' ')
    call add_default ('PRECCDSAS', 1, ' ')

    call addfld ('PCONVB','Pa'    ,1 , 'A','convection base pressure',phys_decomp)
    call addfld ('PCONVT','Pa'    ,1 , 'A','convection top  pressure',phys_decomp)
    call add_default ('PCONVB', 1, ' ')
    call add_default ('PCONVT', 1, ' ')
    ! +++ Yi-Chi : add CIN flag output
    
    call addfld ('SAS_flgmask',''    ,1 , 'A','flags for SAS triggering: 1:CIN, 2:dry PBL,3:cld depth',phys_decomp)
    !call add_default ('SAS_cinflg', 1, ' ')
    call addfld ('SAS_cinmask',''    ,1 , 'A','distance between LFC and original level when CIN flag is on',phys_decomp)
    !call add_default ('SAS_cinmask', 1, ' ')
    call addfld ('SAS_heo ', 'J',    pverp, 'A', 'SAS moist static energy for kbconv',phys_decomp)
    call addfld ('SAS_heso', 'J',    pverp, 'A', 'SAS saturated moist static energy for kbconv',  phys_decomp)

    call addfld ('SAS_xmb ', '',    1, 'A', 'SAS cloud-base mass flux',phys_decomp)
    call addfld ('SAS_xk', 'J',     1, 'A', 'SAS saturated moist static energy for kbconv',  phys_decomp)
    call addfld ('SAS_fld', 'J',    1, 'A', 'SAS large-scale forcing',phys_decomp)
    call addfld ('SAS_dtconv', 'minutes', 1, 'A', 'SAS convective time scale',  phys_decomp)
    call addfld ('SAS_aa1crit', 'J', 1, 'A', 'SAS critical CAPE',  phys_decomp)
    call addfld ('SAS_xomega', 'Pa/s', 1, 'A', 'SAS cloud base vertical velocity',  phys_decomp)
    call addfld ('SAS_qcko', 'kg/kg', 1, 'A', 'SAS saturated updraft moisture',  phys_decomp)
    call addfld ('SAS_duq', 'kg/kg/s', 1, 'A', 'SAS detrained updraft',  phys_decomp)
    ! ---

    call addfld ('CAPE',   'J/kg',       1, 'A', 'Convectively available potential energy', phys_decomp)
   call add_default ('CAPE', 1, ' ')
    call addfld ('FREQSAS ','fraction  ',1  ,'A', 'Fractional occurance of SAS convection',phys_decomp) 
    call add_default ('FREQSAS', 1, ' ')

    call addfld ('SASMTT ', 'K/s',     pver, 'A', 'T tendency - SAS convective momentum transport',phys_decomp)
    call addfld ('SASMTU',  'm/s2',    pver, 'A', 'U tendency - SAS convective momentum transport',  phys_decomp)
    call addfld ('SASMTV',  'm/s2',    pver, 'A', 'V tendency - SAS convective momentum transport',  phys_decomp)

    call addfld ('SASMU',   'kg/m2/s', pver, 'A', 'SAS convection updraft mass flux',   phys_decomp)
    call add_default ('SASMU', 1, ' ')
    call addfld ('SASMD',   'kg/m2/s', pver, 'A', 'SAS convection downdraft mass flux', phys_decomp)
    call add_default ('SASMD', 1, ' ')
    call addfld ('SASMED',   '1/s', pver, 'A', 'SAS convection entrained downdraft',   phys_decomp)
    call add_default ('SASMED', 1, ' ')
    call addfld ('SASMDU',   '1/s', pver, 'A', 'SAS convection detrained updraft mass', phys_decomp)
    call add_default ('SASMDU', 1, ' ')
    call addfld ('SASMDD',   '1/s', pver, 'A', 'SAS convection entrained downdraft',   phys_decomp)
    call add_default ('SASMDD', 1, ' ')
    call addfld ('SASMEU',   '1/s', pver, 'A', 'SAS convection detrained updraft mass', phys_decomp)
    call add_default ('SASMEU', 1, ' ')
    !updraft properties!-------------
    !call addfld ('SASQCKO',   '', pver, 'A', 'SAS convection normalized updraft moisture', phys_decomp)
    !call add_default ('SASQCKO', 1, ' ')
    call addfld ('SASDUQ',   '', pver, 'A', 'SAS updraft detrained moisture', phys_decomp)
    call add_default ('SASDUQ', 1, ' ')
    !----------------------------
    call addfld ('SASDQL',   'kg/kg/s', pver, 'A', 'SAS convection cloud water change before convective transport', phys_decomp)
    call add_default ('SASDQL', 1, ' ')
    call addfld ('SASDQICE',   'kg/kg/s', pver, 'A', 'SAS convection cloud ice change before convective transport', phys_decomp)
    call add_default ('SASDQICE', 1, ' ')

    
    call phys_getopts(history_budget_out = history_budget, history_budget_histfile_num_out = history_budget_histfile_num)
!    if ( history_budget ) then
       call add_default('EVAPTSAS  ', history_budget_histfile_num, ' ')
       call add_default('EVAPQSAS  ', history_budget_histfile_num, ' ')
       call add_default('SASDT     ', history_budget_histfile_num, ' ')
       call add_default('SASDQ     ', history_budget_histfile_num, ' ')
       call add_default('SASDLIQ   ', history_budget_histfile_num, ' ')
       call add_default('SASDICE   ', history_budget_histfile_num, ' ')

       if( cam_physpkg_is('cam4') .or. cam_physpkg_is('cam5') ) then
          call add_default('SASMTT    ', history_budget_histfile_num, ' ')
       end if

!    end if
!
! Limit deep convection to regions below 40 mb
! Note this calculation is repeated in the shallow convection interface
!
    limcnv = 0   ! null value to check against below
    if (pref_edge(1) >= 4.e3_r8) then
       limcnv = 1
    else
       do k=1,plev
          if (pref_edge(k) < 4.e3_r8 .and. pref_edge(k+1) >= 4.e3_r8) then
             limcnv = k
             exit
          end if
       end do
       if ( limcnv == 0 ) limcnv = plevp
    end if
    
    if (masterproc) then
       write(iulog,*)'SAS_CONV_INIT: Deep convection will be capped at intfc ',limcnv, &
            ' which is ',pref_edge(limcnv),' pascals'
    end if
        
    no_deep_pbl = phys_deepconv_pbl()
    call zm_convi(limcnv,no_deep_pbl_in = no_deep_pbl)

     cld_idx         = pbuf_get_index('CLD')
     icwmrdp_idx     = pbuf_get_index('ICWMRDP')
     rprddp_idx      = pbuf_get_index('RPRDDP')
     fracis_idx      = pbuf_get_index('FRACIS')
     nevapr_dpcu_idx = pbuf_get_index('NEVAPR_DPCU')

    prec_dp_idx     = pbuf_get_index('PREC_DP')
    snow_dp_idx     = pbuf_get_index('SNOW_DP')

end subroutine sas_conv_init
!=========================================================================================
!subroutine zm_conv_tend(state, ptend, tdt, pbuf)
! Yi-Chi : Modify the output of sas_conv_tend
!subroutine zm_conv_tend(prec    , &
!     pblh    ,mcon    ,cme     ,          &
!     tpert   ,dlf     ,pflx    ,zdu      , &
!     rliq    , &
!     ztodt   ,snow    ,&
!     jctop   ,jcbot , &
!     state   ,ptend_all   ,landfrac   ,pbuf  )

subroutine sas_conv_tend( mcon    ,     &
     dlf     ,zdu      , &
     rliq    , &
     ztodt   , &
     jctop   ,jcbot , &
     state   ,ptend_all   ,landfrac   ,pbuf  ,ideepout)
! (in)pblh, tperb, (out)cme, pflx are removed.  
!** dlf, zdu, and rliq are all related to detraining cloud water**
!   need to be further refined.

   use module_gfs_funcphys, only  : gfuncphys ! for saturated vapor computation
   use module_gfs_machine ,  only : kind_phys ! for input: Yi-Chi
   use cam_history,   only: outfld
   use physics_types, only: physics_state, physics_ptend
   use physics_types, only: physics_ptend_init, physics_update
   use physics_types, only: physics_state_copy, physics_state_dealloc
   use physics_types, only: physics_ptend_sum, physics_ptend_dealloc

   use phys_grid,     only: get_lat_p, get_lon_p
   use time_manager,  only: get_nstep, is_first_step
!   use phys_buffer,   only: pbuf_size_max, pbuf_fld, pbuf_old_tim_idx
   use physics_buffer, only : pbuf_get_field, physics_buffer_desc, pbuf_old_tim_idx

   use constituents,  only: pcnst, cnst_get_ind, cnst_is_convtran1
   use check_energy,  only: check_energy_chng
   use physconst,     only: gravit
   use phys_control,  only: cam_physpkg_is

   ! Arguments

   type(physics_state), intent(in ) :: state          ! Physics state variables
   type(physics_ptend), intent(out) :: ptend_all          ! indivdual parameterization tendencies
   type(physics_buffer_desc), pointer :: pbuf(:)

   real(r8), intent(in) :: ztodt                       ! 2 delta t (model time increment)
   real(r8), intent(in) :: landfrac(pcols)             ! RBN - Landfrac 

   real(r8), intent(out) :: mcon(pcols,pverp)  ! Convective mass flux--m sub c
   real(r8), intent(out) :: dlf(pcols,pver)    ! scattrd version of the detraining cld h2o tend
!   real(r8), intent(out) :: pflx(pcols,pverp)  ! scattered precip flux at each level
!   real(r8), intent(out) :: cme(pcols,pver)    ! cmf condensation - evaporation
   real(r8), intent(out) :: zdu(pcols,pver)    ! detraining mass flux

!   real(r8), intent(out) :: prec(pcols)   ! total precipitation
!   real(r8), intent(out) :: snow(pcols)   ! snow from SAS convection 
   real(r8), intent(out) :: rliq(pcols) ! reserved liquid (not yet in cldliq) for energy integrals
   integer, intent(out) :: ideepout(pcols) ! flag for shallow convection

   ! Local variables

   integer :: i,k,m
   integer :: ilon                      ! global longitude index of a column
   integer :: ilat                      ! global latitude index of a column
   integer :: nstep
   integer :: ixcldice, ixcldliq      ! constituent indices for cloud liquid and ice water.
   integer :: lchnk                   ! chunk identifier
   integer :: ncol                    ! number of atmospheric columns
   integer :: itim                    ! for physics buffer fields

   real(r8) :: ftem(pcols,pver)              ! Temporary workspace for outfld variables
   real(r8) :: ntprprd(pcols,pver)    ! evap outfld: net precip production in layer
   real(r8) :: ntsnprd(pcols,pver)    ! evap outfld: net snow production in layer
   real(r8) :: tend_s_snwprd  (pcols,pver) ! Heating rate of snow production
   real(r8) :: tend_s_snwevmlt(pcols,pver) ! Heating rate of evap/melting of snow
   real(r8) :: fake_dpdry(pcols,pver) ! used in convtran call

   ! physics types
   type(physics_state) :: state1        ! locally modify for evaporation to use, not returned
   type(physics_ptend) :: ptend_loc     ! package tendencies

   ! physics buffer fields
   real(r8), pointer, dimension(:)   :: prec         ! total precipitation
   real(r8), pointer, dimension(:)   :: snow         ! snow from ZM convection
   real(r8), pointer, dimension(:,:) :: cld
   real(r8), pointer, dimension(:,:) :: ql           ! wg grid slice of cloud liquid water.
   real(r8), pointer, dimension(:,:) :: rprd         ! rain production rate
   real(r8), pointer, dimension(:,:,:) :: fracis  ! fraction of transported species that are insoluble
   real(r8), pointer, dimension(:,:) :: evapcdp      ! Evaporation of deep convective precipitation
   real(r8), pointer, dimension(:,:) :: flxprec      ! Convective-scale flux of precip at interfaces (kg/m2/s)
   real(r8), pointer, dimension(:,:) :: flxsnow      ! Convective-scale flux of snow   at interfaces (kg/m2/s)
   real(r8), pointer, dimension(:,:) :: dp_cldliq
   real(r8), pointer, dimension(:,:) :: dp_cldice
! Yi-Chi : whether this should change to integer
   real(r8) :: jctop(pcols)  ! o row of top-of-deep-convection indices passed out.
   real(r8) :: jcbot(pcols)  ! o row of base of cloud indices passed out.
   real(r8) :: pcont(pcols), pconb(pcols), freqzm(pcols)

   ! history output fields
   real(r8) :: cape(pcols)        ! w  convective available potential energy.
   real(r8) :: mu_out(pcols,pver)
   real(r8) :: md_out(pcols,pver)

   ! +++ Yi-Chi : remove all variables for momentum transport
   ! used in momentum transport calculation
   ! ---
   integer  :: ii
   logical  :: lq(pcnst)

   !---------------------------------------
   ! Yi-Chi : variables used for SAS calculation
   logical  :: convtran_flag, zmevap_flag
   integer  :: iconst, igpvs
   integer  :: kbot(pcols), ktop(pcols)
   integer  :: ncloud, jcap
   real(r8)  :: rcs(pcols)         ! = 1 as in WRF
   real(r8) :: slimsk(pcols)      ! land mask : land =1, sea = otherwise
   real(r8) :: dpp(pcols,pver)    ! difference of pressure level
   real(r8) :: prsl(pcols,pver)   ! mid-level pressure
   real(r8) :: ps(pcols)          ! surface pressure
   real(r8) :: dot(pcols,pver)    ! vertical pressure velocity
   real(r8) :: phil(pcols,pver)   ! surface geopotential
   real(r8) :: qlsas(pcols,pver,2)   ! If only cloud water 1:cloud water; If two sepicies: ice; 2: cloud water
   real(r8) :: q1sas(pcols,pver)  ! water vapor for SAS
   real(r8) :: q1(pcols,pver,pcnst)     ! water vapor of t+dt; pcnst=1 default
   real(r8) :: t1(pcols,pver)     ! temperature of t+dt
   real(r8) :: s1(pcols,pver)     ! dry static energy at t+dt
   real(r8) :: u1(pcols,pver)     ! zonal velocity of t+dt
   real(r8) :: v1(pcols,pver)     ! meridional velocity of t+dt
   real(r8) :: ud_mf(pcols,pver)  ! mass flux of updraft
   real(r8) :: dd_mf(pcols,pver)  ! mass flux of downdraft
   ! Yi-Chi: temporary variables for vertical reversal
   real(r8) :: dpptmp(pcols,pver)
   real(r8) :: prsltmp(pcols,pver)
   real(r8) :: dottmp(pcols,pver)
   real(r8) :: philtmp(pcols,pver)
   real(r8) :: q1sastmp(pcols,pver)
   real(r8) :: qlsastmp(pcols,pver,2)! environmental cloud condensate mixing ratio
   ! +++ Yi-Chi (Sept 2013)
   real(r8) :: qlicsastmp(pcols,pver)! in-cloud mixing ratio
   ! ---
   real(r8) :: t1tmp(pcols,pver)     ! temperature of t+dt
   real(r8) :: u1tmp(pcols,pver)     ! zonal velocity of t+dt
   real(r8) :: v1tmp(pcols,pver)     ! meridional velocity of t+dt
   real(r8) :: ud_mftmp(pcols,pver)  ! mass flux of updraft
   real(r8) :: dd_mftmp(pcols,pver)  ! mass flux of downdraft
   real(r8) :: dlftmp(pcols,pver)    ! scattrd version of the detraining cld h2o tend
   real(r8) :: zdutmp(pcols,pver)    ! detraining mass flux
   real(r8) :: rprdtmp(pcols,pver)         ! rain production rate
   real(r8) :: evapcdptmp(pcols,pver)     ! Evaporation of deep convective precipitation
!  Yi-Chi : Apr 2012
   real(r8) :: dptmp(pcols,pver)  ! tmp variable for dp(i,k)
   real(r8) :: eutmp(pcols,pver)   ! updraft entrainment
   real(r8) :: dutmp(pcols,pver)   ! updraft detrainment
   real(r8) :: edtmp(pcols,pver)   ! downdraft entrainment
   real(r8) :: ddtmp(pcols,pver)   ! downdraft detrainment
   real(r8) :: euout(pcols,pver)   ! updraft entrainment
   real(r8) :: duout(pcols,pver)   ! updraft detrainment
   real(r8) :: edout(pcols,pver)   ! downdraft entrainment
   real(r8) :: ddout(pcols,pver)   ! downdraft detrainment
! Yi-Chi : Apr 2012 cloud ice
   logical  :: cloudice           ! treatment of cloud ice
   integer  :: k_inv                    !  Vertical index for incoming fields [ no ]
   !real(r8) :: qice(pcols,pver)   ! cloud ice from SAS scheme
   !------------------
   ! for log files
   !LOGICAL :: ex,op
   !character(10) :: myfile 
   !----------------------------------------------------------------------

   ! initialize
   lchnk = state%lchnk
   ncol  = state%ncol
   nstep = get_nstep()

   ftem = 0._r8   
   mu_out(:,:) = 0._r8
   md_out(:,:) = 0._r8
   ! Yi-Chi : initialization
   jt(:,:)    = pver
   maxg(:,:)  = 0
!   wind_tends(:ncol,:pver,:) = 0.0_r8

   call physics_state_copy(state,state1)   ! copy state to local state1.
   lq(:) = .FALSE.
   lq(1) = .TRUE.

   call physics_ptend_init(ptend_loc, state%psetcols, 'sas_convr', ls=.true., lq=lq, lu=.true., lv=.true.)! initialize local ptend type
!   call physics_ptend_init(ptend_loc, state%psetcols, 'sas_convr', ls=.true., lq=lq)

!   call physics_ptend_init(ptend_loc)  ! initialize local ptend type
!   call physics_ptend_init(ptend_all)  ! initialize output ptend type
!   call physics_tend_init(tend)        ! tend type here is a null place holder

!
! Associate pointers with physics buffer fields
!
   itim = pbuf_old_tim_idx()
   call pbuf_get_field(pbuf, cld_idx,         cld,    start=(/1,1,itim/), kount=(/pcols,pver,1/) )
   call pbuf_get_field(pbuf, icwmrdp_idx,     ql )
   call pbuf_get_field(pbuf, rprddp_idx,      rprd )
   call pbuf_get_field(pbuf, fracis_idx,      fracis, start=(/1,1,1/),    kount=(/pcols, pver, pcnst/) )
    call pbuf_get_field(pbuf, nevapr_dpcu_idx, evapcdp )
    call pbuf_get_field(pbuf, prec_dp_idx,     prec )
    call pbuf_get_field(pbuf, snow_dp_idx,     snow )


!
! Begin with Zhang-McFarlane (1996) convection parameterization
!
   call t_startf ('sas_convr')

!---------------------------------------------
! Yi-Chi: replace this part with the SAS core
!   call sascnvn(im,ix,km,jcap,delt,del,prsl,ps,phil,ql,
!     &     q1,t1,u1,v1,rcs,cldwrk,rn,kbot,ktop,kcnv,slimsk,
!     &     dot,ncloud,ud_mf,dd_mf,dt_mf) 
! 1. initialization
! 2. reverse vertical levels

   ncloud = 1
   jcap   = 126
   ! Yi-Chi : follow WRF, all Pa->bar (*0.00001) ->cb (*100)
   dpp  = state%pdel*.001
   prsl = state%pmid*.001
   ps   = state%ps*.001
   dot  = state%omega*0.001 !Pa/s -> bar/s(*.00001) -> cb(*100)
   !write(iulog,*) 'Yi-Chi(sas_conv_intr) state%omega:',state%omega   
   phil = state%zm*gravit
   ! Yi-Chi (Apr27, 2012) : cloud ice
   !cloudice = .true.
   ! Yi-Chi (July 27, 2012) : cloud ice
   cloudice = .false.  

!   write(iulog,*) 'Yi-Chi(sas_conv_intr before sascnvn),ql,qlsas:',ql,qlsas
   eutmp(:,:) = 0.
   dutmp(:,:) = 0.
   edtmp(:,:) = 0.
   ddtmp(:,:) = 0.
   euout(:,:) = 0.
   duout(:,:) = 0.
   edout(:,:) = 0.
   ddout(:,:) = 0.
   !qckotmp(:,:) = 0.
   !qckoout(:,:) = 0.
   !duqtmp(:,:)  = 0.
   !duqout(:,:)  = 0.
   dlftmp(:,:)  = 0.

   q1sas(:,:)    = 0.
   q1sastmp(:,:) = 0.

   qlsastmp(:,:,:) = 0.
   qlicsastmp(:,:) = 0.

   do i = 1, ncol
      rcs(i)    = 1
     if(landfrac(i).gt.0.5) then
      slimsk(i) = 1. ! for land
     else
      slimsk(i) = 0. ! for water
     endif
   do k = 1, pver
     ! Yi-Chi : Apr27, 2012
     if(cloudice) then
      qlsas(i,k,1) = state%q(i,k,3) ! (Apr 2012) cloud ice
      qlsas(i,k,2) = state%q(i,k,2)   ! cld water
     else
      qlsas(i,k,1) = state%q(i,k,2)  ! (Mar 2012) cloud water
      qlsas(i,k,2) = -1000.   ! set to show only cld water is considered
     endif
      q1sas(i,k)   = state%q(i,k,1) ! save water vapor to SAS format
      do iconst = 1, pcnst
      q1(i,k,iconst)   = 0.         ! initialize CAM format
      enddo
      rprd(i,k)    = 0.
      evapcdp(i,k) = 0.
      zdu(i,k)     = 0.
      dlf(i,k)     = 0.
   enddo
      rliq(i)      = 0.
      prec(i)      = 0.
   enddo

!   write(iulog,*) 'Yi-Chi(sas_conv_intr before sascnvn),pcnst:',pcnst
!   write(iulog,*) 'Yi-Chi(sas_conv_intr before sascnvn),ql,qlsas:',ql,qlsas
!   q1      = state%q !(pcols,pver,pcnst):1 as water vapor, 2 as cloud water
   t1      = state%t
   u1      = state%u
   v1      = state%v

!   write(iulog,*) 'Yi-Chi(sas_conv_intr): pressure level prsl before reverse:',prsl
   ! reverse the vertical level: 
   ! CAM: from top to bottom; while GFS: from bottom to top
   ! Yi-Chi : use a temporary matrix for it
!   do i = 1, ncol
   do k = 1, pver
      k_inv               = pver - k + 1
      dpptmp (:ncol,k)    = dpp (:ncol,k_inv)
      prsltmp(:ncol,k)    = prsl(:ncol,k_inv)
      dottmp (:ncol,k)    = dot (:ncol,k_inv)
      philtmp(:ncol,k)    = phil(:ncol,k_inv)
      !if(cloudice) then ! Yi-Chi : Apr27 2012
        qlsastmp(:ncol,k,1) = qlsas(:ncol,k_inv,1)! cloud ice
        qlsastmp(:ncol,k,2) = qlsas(:ncol,k_inv,2)! cloud water
      !else
      !  qlsastmp(i,k,1) = qlsas(i,pver-k+1,1)! only for cloud water case
      !endif
      q1sastmp(:ncol,k)   = q1sas(:ncol,k_inv)  ! mixing ratio
      t1tmp(:ncol,k)      = t1(:ncol,k_inv)
      u1tmp(:ncol,k)      = u1(:ncol,k_inv)
      v1tmp(:ncol,k)      = v1(:ncol,k_inv)
   enddo
!   enddo
   ! write(iulog,*) 'Yi-Chi(sas_conv_intr): cloud water qlsas before SAS:',qlsastmp(i,k,1)
!    write(iulog,*) 'Yi-Chi(sas_conv_intr): water vapor q1sas before SAS:',q1sas
!   write(iulog,*) 'Yi-Chi(sas_conv_intr): pressure increment  before SAS:',dpp
!   write(iulog,*) 'Yi-Chi(sas_conv_intr): pressure level prsl before SAS:',prsl
!   write(iulog,*) 'Yi-Chi(sas_conv_intr): pressure level pver before SAS:',pver
   ! compute all tables for saturated vapor pressure
   call gfuncphys

   ! replace cldwrf with cape
   !call sascnvn(ncol,ncol,pver,jcap,.5_r8*ztodt,dpptmp,prsltmp, &
   !call sascnvn(ncol,ncol,pver,jcap,ztodt,dpptmp,prsltmp, &
   ! +++ ycw add lchnk
   call sascnvn(lchnk, ncol,ncol,pver,jcap,ztodt,dpptmp,prsltmp, &
   ! --- ycw add END
          ! +++ Yi-Chi : add qlicsastmp as in-cloud mixing ratio
          !ps,philtmp,qlsastmp, &
          ps,philtmp,qlsastmp, qlicsastmp,  &
          ! ---
          q1sastmp,t1tmp,u1tmp,v1tmp,rcs,cape,prec,kbot,ktop,slimsk, &
          dottmp,ncloud,ud_mftmp,dd_mftmp, & !)
          rprdtmp, evapcdptmp, zdutmp, dlftmp, & !) ! Yi-Chi: Mar 2012
!          ideep(:,lchnk), lengath(lchnk))
          ideep(:,lchnk), lengath(lchnk), & !)
          eutmp,dutmp,edtmp,ddtmp)
          ! +++ ycw move arguments to outfld
          !,qckotmp,duqtmp,flgmaskout,cinmaskout, &!) !Yi-Chi: Apr 1, 2012
          !heotmp, hesotmp, heo1tmp, heso1tmp, &!)
          !xmbout, xkout, fldout,dtconvout, aa1crit,xomega)
          ! ---
!     &     dot,ncloud,ud_mf,dd_mf,dt_mf)
   ! remove: kcnv, dt_mf

   ! reverse the vertical level: back to CAM format
   ! CAM: from top to bottom; while GFS: from bottom to top
!   do i = 1, ncol
   do k = 1, pver
      k_inv            = pver-k+1
      t1(:ncol,k)      = t1tmp (:ncol,k_inv)
      u1(:ncol,k)      = u1tmp (:ncol,k_inv)
      v1(:ncol,k)      = v1tmp (:ncol,k_inv)
      ud_mf(:ncol,k)   = ud_mftmp(:ncol,k_inv)
      dd_mf(:ncol,k)   = -1.*dd_mftmp(:ncol,k_inv) ! upward as positive
!      qlsas(i,k,2) = qlsas(i,pver-k+1,2)
!      q1sas(i,k)   = q1sas(i,pver-k+1)
      rprd(:ncol,k)    = rprdtmp(:ncol,k_inv)
      evapcdp(:ncol,k) = evapcdptmp(:ncol,k_inv) ! from kg/kg to kg/(kg s)
      zdu(:ncol,k)     = zdutmp(:ncol,k_inv)
      dlf(:ncol,k)     = dlftmp(:ncol,k_inv) ! detrain/entrain mass
      dptmp(:ncol,k)= state%pdel(:ncol,k_inv)
      euout(:ncol,k)= eutmp(:ncol,k_inv)!/dptmp
      duout(:ncol,k)= dutmp(:ncol,k_inv)!/dptmp
      edout(:ncol,k)= edtmp(:ncol,k_inv)!/dptmp
      ddout(:ncol,k)= ddtmp(:ncol,k_inv)!/dptmp
   ! updraft properties
   !   qckoout(:ncol,k) = qckotmp(:ncol,k_inv)
   !   duqout (:ncol,k) = duqtmp (:ncol,k_inv)
   ! moist static energy heo and saturated MSE heso
   !   on mid-layer : km
      !heo1out (i,k) = heo1tmp (:ncol,k_inv)
      !heso1out(i,k) = heso1tmp(:ncol,k_inv)
   enddo
   prec(:ncol)      = prec(:ncol)/ztodt ! from m to m/s for CAM
   jctop(:ncol)     = pver-min(max(ktop(:ncol),1),pver)+1
   jcbot(:ncol)     = pver-max(min(kbot(:ncol),pver),1)+1

!   enddo
!   write(iulog,*) 'Yi-Chi(sas_conv_intr):pver, jcbot',pver,jcbot,'jctop',jctop
!   write(iulog,*) 'Yi-Chi(sas_conv_intr):rprd, evapcdp',rprd,evapcdp
!   write(iulog,*) 'Yi-Chi(sas_conv_intr): dlf',dlf 
!   write(iulog,*) 'Yi-Chi(sas_conv_intr): prec(m/s)',prec
!   write(iulog,*) 'Yi-Chi(sas_conv_intr):ztodt',ztodt
!    write(iulog,*) 'Yi-Chi(sas_conv_intr): cloud water qlsas after SAS:',qlsastmp(i,k,1)

   !  Yi-Chi : caution: need further study on the constituents in CAM
   !do i = 1, ncol
   do k = 1, pver
         k_inv       = pver-k+1
         q1(:ncol,k,1)   = q1sastmp  (:ncol,k_inv)    ! save water vapor to q1 CAM format; pcnst = 1 for water vapor
         ql(:ncol,k)     = qlicsastmp(:ncol,k_inv)
       if(cloudice) then ! Yi-Chi : apr27
         ! Yi-Chi: here only consider the cloud water case in sascnvn
         q1(:ncol,k,3)   = qlsastmp  (:ncol,k_inv,1)      ! cloud ice
         q1(:ncol,k,2)   = qlsastmp  (:ncol,k_inv,2)      ! grid-scale liquid water
       else
         q1(:ncol,k,2)   = qlsastmp  (:ncol,k_inv,1)      ! grid-scale liquid water
       endif
   enddo
   !enddo


   call outfld('SASMED', edout, pcols, lchnk)
   call outfld('SASMDU', duout, pcols, lchnk)
   call outfld('SASMDD', ddout, pcols, lchnk)
   call outfld('SASMEU', euout, pcols, lchnk)
!   write(iulog,*) 'Yi-Chi(sas_conv_intr):ql, q1',ql,q1
   ! calculate snow as zero
   do i = 1, ncol
    snow(i) = 0.
   end do
   ! calculate mcon = ud_mf + dd_mf
   do k = 1, pver
   do i = 1, ncol
      mcon(i,k) = ud_mf(i,k) + dd_mf(i,k)
      if(mcon(i,k).lt.0.) then 
         ! Yi-Chi: let levels below cloud base as zero
         mcon(i,k) = 0.
      !else ! Yi-Chi: even downdraft is larger than updraft, let it pass
      !   mcon(i,k) = abs(ud_mf(i,k) + dd_mf(i,k))
      endif
   enddo
   enddo
   ! calculate dlf(i,k) : detraining cloud water tendency
   ! calculate reserved cloud water
   do k = 1, pver
   do i = 1, ncol
    ! du_mf(i,k)*qlsas(i,k,2) ! in CAM, it's du(i,k)
    !dpp (kPa)
    rliq(i)  = rliq(i) + dlf(i,k)*(dpp(i,k)*1000.)/gravit
   end do
   end do
   rliq(:ncol) = rliq(:ncol) /1000._r8
   !write(iulog,*) 'Yi-Chi(sas_conv_intr):rliq ',rliq*1000.*latvap
   !write(iulog,*) 'Yi-Chi(sas_conv_intr):prec ',prec*1000.*latvap
   !write(iulog,*) 'Yi-Chi(sas_conv_intr):latvap,cpair: ',latvap,cpair

   ! Yi-Chi (July4 2012): remove direct detrainment of cloud water
   !ptend_loc%lq(2) = .TRUE.
   ! Yi-Chi (Apr 2012) : cloud ice
   if(cloudice) then 
      ptend_loc%lq(3) = .TRUE.
   endif

   ptend_loc%u(:ncol,:pver) = (u1(:ncol,:pver)-state%u(:ncol,:pver))/ztodt
   ptend_loc%v(:ncol,:pver) = (v1(:ncol,:pver)-state%v(:ncol,:pver))/ztodt
   ptend_loc%s(:ncol,:pver) = (t1(:ncol,:pver)-state%t(:ncol,:pver))*cpair/ztodt
   ! in zm_convr: heat = ds/dt * cpres
   ! in ZM_convr: only water vapor is changed; cloud water is passed through dlf.
   ptend_loc%q(:ncol,:pver,1) =  &
           (q1(:ncol,:pver,1)-state%q(:ncol,:pver,1))/ztodt
   ! Yi-Chi (July4 2012)
   !   remove the direct detrainment of cloud water
   !! include the change of ql
   !!ptend_loc%q(:ncol,:pver,2) =  &
   !!        (ql(:ncol,:pver)-state%q(:ncol,:pver,2))/ztodt
   if(cloudice) then ! Yi-Chi : apr27 for cloud ice
      ptend_loc%q(:ncol,:pver,3) =  &
           (q1(:ncol,:pver,3)-state%q(:ncol,:pver,3))/ztodt
   endif

!   do iconst = 1, pcnst
!      ptend_loc%q(:ncol,:pver,iconst) =  &
!           (q1(:ncol,:pver,iconst)-state%q(:ncol,:pver,iconst))/ztodt
!   enddo

   call outfld('CAPE', cape, pcols, lchnk)        ! RBN - CAPE output
!
! Output fractional occurance of SAS convection
!
   ideepout(:) = 0
   freqzm  (:) = 0._r8
   do i = 1,lengath(lchnk)
      freqzm  (ideep(i,lchnk)) = 1.0_r8
      ideepout(ideep(i,lchnk)) = 1
   end do
   call outfld('FREQSAS  ',freqzm          ,pcols   ,lchnk   )

!
! Convert mass flux from reported mb/s to kg/m^2/s
!
   mcon(:ncol,:pver) = mcon(:ncol,:pver) * 100._r8/gravit

   ! Store upward and downward mass fluxes in un-gathered arrays
   ! + convert from mb/s to kg/m^2/s
   do i=1,lengath(lchnk) 
      do k=1,pver
         ii = ideep(i,lchnk)
         ! Yi-Chi (Apr25)
         eu(i,k,lchnk)= euout(ii,k)
         du(i,k,lchnk)= duout(ii,k)
         ed(i,k,lchnk)= edout(ii,k)
         dd(i,k,lchnk)= ddout(ii,k)
         ! Yi-Chi (Apr12) : ud_mf and dd_mf are in unit of mb/s
         !   the transformation is in sas_conv.F90
         mu(i,k,lchnk) = ud_mf(ii,k)
         md(i,k,lchnk) = dd_mf(ii,k)
         !!Yi-Chi
         !mu_out(ii,k) =  ud_mf(ii,k)!* 100._r8/gravit
         !md_out(ii,k) =  dd_mf(ii,k)!* 100._r8/gravit
         mu_out(ii,k) = mu(i,k,lchnk) * 100._r8/gravit
         md_out(ii,k) = md(i,k,lchnk) * 100._r8/gravit
      end do
   end do

!   write(iulog,*) 'Yi-Chi(sas_conv_intr):mcon',mcon
!   write(iulog,*) 'Yi-Chi(sas_conv_intr): ud_mf',ud_mf
!   write(iulog,*) 'Yi-Chi(sas_conv_intr): dd_mf',dd_mf

   call outfld('SASMU', mu_out(1,1), pcols, lchnk)
   call outfld('SASMD', md_out(1,1), pcols, lchnk)

!   ptend_loc%name  = 'sas_convr'
!   ptend_loc%ls    = .TRUE.
!   ptend_loc%lq(1) = .TRUE.

   ftem(:ncol,:pver) = ptend_loc%s(:ncol,:pver)/cpair
   call outfld('SASDT    ',ftem           ,pcols   ,lchnk   )
   call outfld('SASDQ    ',ptend_loc%q(1,1,1) ,pcols   ,lchnk   )
   call outfld('SASDQL   ',ptend_loc%q(1,1,2) ,pcols   ,lchnk   )
   call outfld('SASDQICE ',ptend_loc%q(1,1,3) ,pcols   ,lchnk   )
   call t_stopf ('sas_convr')

!    do i = 1,pcols
!    do i = 1,nco
! Yi-Chi (Apr 2012) : prepare for convective transport
! (1) save jctop and jcbot into gathered arrays
! (2) output dp
   do i = 1,lengath(lchnk)
     jt(i,lchnk)   = int(jctop(ideep(i,lchnk)))
     maxg(i,lchnk) = int(jcbot(ideep(i,lchnk)))
     dp(i,:,lchnk) = state%pdel(ideep(i,lchnk),:)/100._r8
   end do

   pcont(:ncol) = state%ps(:ncol)
   pconb(:ncol) = state%ps(:ncol)
   do i = 1,lengath(lchnk)
       if (maxg(i,lchnk).gt.jt(i,lchnk)) then
          pcont(ideep(i,lchnk)) = state%pmid(ideep(i,lchnk),jt(i,lchnk))  ! gathered array (or jctop ungathered)
          pconb(ideep(i,lchnk)) = state%pmid(ideep(i,lchnk),maxg(i,lchnk))! gathered array
       endif
       !     write(iulog,*) ' pcont, pconb ', pcont(i), pconb(i), cnt(i), cnb(i)
    end do
    call outfld('PCONVT  ',pcont          ,pcols   ,lchnk   )
    call outfld('PCONVB  ',pconb          ,pcols   ,lchnk   )

  ! This name triggers a special case in physics_types.F90:physics_update()
  call physics_ptend_init(ptend_all, state%psetcols, 'convect_deep')

  ! add tendency from this process to tendencies from other processes
  call physics_ptend_sum(ptend_loc,ptend_all, ncol)

  ! update physics state type state1 with ptend_loc
  call physics_update(state1, ptend_loc, ztodt)

!  +++++
!  Yi-Chi :
!   start of rainfall reevaporation with zm_conv_evap subroutine
!
! As SAS has already considered rainfall evaporation,
!   a flag is set to
!  ----
  zmevap_flag = .false.

  if(zmevap_flag) then

  ! initialize ptend for next process
  lq(:) = .FALSE.
  lq(1) = .TRUE.
  call physics_ptend_init(ptend_loc, state1%psetcols, 'zm_conv_evap', ls=.true., lq=lq)

   call t_startf ('zm_conv_evap')
!
! Determine the phase of the precipitation produced and add latent heat of fusion
! Evaporate some of the precip directly into the environment (Sundqvist)
! Allow this to use the updated state1 and the fresh ptend_loc type
! heating and specific humidity tendencies produced
!

    call pbuf_get_field(pbuf, dp_flxprc_idx, flxprec    )
    call pbuf_get_field(pbuf, dp_flxsnw_idx, flxsnow    )
    call pbuf_get_field(pbuf, dp_cldliq_idx, dp_cldliq  )
    call pbuf_get_field(pbuf, dp_cldice_idx, dp_cldice  )
    dp_cldliq(:ncol,:) = 0._r8
    dp_cldice(:ncol,:) = 0._r8

    call zm_conv_evap(state1%ncol,state1%lchnk, &
         state1%t,state1%pmid,state1%pdel,state1%q(:pcols,:pver,1), &
         ptend_loc%s, tend_s_snwprd, tend_s_snwevmlt, ptend_loc%q(:pcols,:pver,1), &
         rprd, cld, ztodt, &
         prec, snow, ntprprd, ntsnprd , flxprec, flxsnow)

    evapcdp(:ncol,:pver) = ptend_loc%q(:ncol,:pver,1)
!
! Write out variables from zm_conv_evap
!
   ftem(:ncol,:pver) = ptend_loc%s(:ncol,:pver)/cpair
   call outfld('EVAPTSAS ',ftem           ,pcols   ,lchnk   )
   ftem(:ncol,:pver) = tend_s_snwprd  (:ncol,:pver)/cpair
   call outfld('FZSNTSAS ',ftem           ,pcols   ,lchnk   )
   ftem(:ncol,:pver) = tend_s_snwevmlt(:ncol,:pver)/cpair
   call outfld('EVSNTSAS ',ftem           ,pcols   ,lchnk   )
   call outfld('EVAPQSAS ',ptend_loc%q(1,1,1) ,pcols   ,lchnk   )
   call outfld('SASFLXPRC', flxprec, pcols, lchnk)
   call outfld('SASFLXSNW', flxsnow, pcols, lchnk)
   call outfld('SASNTPRPD', ntprprd, pcols, lchnk)
   call outfld('SASNTSNPD', ntsnprd, pcols, lchnk)
   call outfld('SASEIHEAT', ptend_loc%s, pcols, lchnk)
   call outfld('CMFMCDSAS   ',mcon ,  pcols   ,lchnk   )
   call outfld('PRECCDSAS   ',prec,  pcols   ,lchnk   )


   call t_stopf ('zm_conv_evap')

   call outfld('PRECZ   ', prec   , pcols, lchnk)

  ! add tendency from this process to tend from other processes here
  call physics_ptend_sum(ptend_loc,ptend_all, ncol)

  ! update physics state type state1 with ptend_loc
  call physics_update(state1, ptend_loc, ztodt)

endif ! Yi-Chi : end of zmevaporation flag

! -----------end of zmevap part-------

! add convtran flag for debugging
! +++
   convtran_flag = .false.
!   convtran_flag = .true.
! ---
   if(convtran_flag) then
  ! Transport cloud water and ice only
   call cnst_get_ind('CLDLIQ', ixcldliq)
   call cnst_get_ind('CLDICE', ixcldice)

   lq(:)  = .FALSE.
   lq(2:) = cnst_is_convtran1(2:)
   call physics_ptend_init(ptend_loc, state1%psetcols, 'convtran1', lq=lq)


   ! dpdry is not used in this call to convtran since the cloud liquid and ice mixing
   ! ratios are moist
   fake_dpdry(:,:) = 0._r8

   call t_startf ('convtran1')
   call convtran_sas (lchnk,                                        &
                  ptend_loc%lq,state1%q, pcnst,  mu(:,:,lchnk), md(:,:,lchnk),   &
                  du(:,:,lchnk), eu(:,:,lchnk), ed(:,:,lchnk), dd(:,:,lchnk),  dp(:,:,lchnk), &
                  jt(:,lchnk),maxg(:,lchnk), ideep(:,lchnk), 1, lengath(lchnk),  &
                  nstep,   fracis,  ptend_loc%q, fake_dpdry)

   call t_stopf ('convtran1')

   call outfld('SASDICE ',ptend_loc%q(1,1,ixcldice) ,pcols   ,lchnk   )
   call outfld('SASDLIQ ',ptend_loc%q(1,1,ixcldliq) ,pcols   ,lchnk   )

   ! add tendency from this process to tend from other processes here
   call physics_ptend_sum(ptend_loc,ptend_all, ncol)


   end if

   call physics_state_dealloc(state1)
   call physics_ptend_dealloc(ptend_loc)


end subroutine sas_conv_tend



!==========================================================================
! Yi-Chi : 
!  * add back convective transport subroutine (April, 2012) 
!  * remove convection transport subroutine March, 2012
!       
! See original subdirectory
!--------
subroutine sas_conv_tend_2( state,  ptend,  ztodt, pbuf  )

   use physics_types, only: physics_state, physics_ptend, physics_ptend_init
   use time_manager,  only: get_nstep
   use physics_buffer, only: pbuf_get_index, pbuf_get_field, physics_buffer_desc
   use constituents,  only: pcnst, cnst_get_ind, cnst_is_convtran1
   use error_messages, only: alloc_err

! Arguments
   type(physics_state), intent(in )   :: state          ! Physics state variables
   type(physics_ptend), intent(out)   :: ptend          ! indivdual parameterization tendencies

   type(physics_buffer_desc), pointer :: pbuf(:)

   real(r8), intent(in) :: ztodt                          ! 2 delta t (model time increment)

! Local variables
   integer :: i, lchnk, istat
   integer :: nstep
   real(r8), dimension(pcols,pver) :: dpdry

! physics buffer fields
   integer itim, ifld
   real(r8), pointer, dimension(:,:,:) :: fracis  ! fraction of transported species that are insoluble
   logical   :: lq(pcnst)

! Initialize
!
  lq(:) = .FALSE.
  lq(:) = .not. cnst_is_convtran1(:)
  call physics_ptend_init(ptend, state%psetcols, 'convtran2', lq=lq )

!
! Associate pointers with physics buffer fields
!
   ifld = pbuf_get_index('FRACIS')
   call pbuf_get_field(pbuf, fracis_idx, fracis, start=(/1,1,1/), kount=(/pcols, pver, pcnst/) )

!
! Transport all constituents except cloud water and ice
!

  lchnk = state%lchnk

   nstep = get_nstep()

   if (any(ptend%lq(:))) then
      ! initialize dpdry for call to convtran
      ! it is used for tracers of dry mixing ratio type
      dpdry = 0._r8
      do i = 1,lengath(lchnk)
         dpdry(i,:) = state%pdeldry(ideep(i,lchnk),:)/100._r8
      end do

      call t_startf ('convtran2')
! +++ Yi-Chi :
   call convtran_sas (lchnk,                                        &
                  ptend%lq,state%q, pcnst,  mu(:,:,lchnk), md(:,:,lchnk),   &
                  du(:,:,lchnk), eu(:,:,lchnk), ed(:,:,lchnk), dd(:,:,lchnk),  dp(:,:,lchnk), &
                  jt(:,lchnk),maxg(:,lchnk), ideep(:,lchnk), 1, lengath(lchnk),  &
                  nstep,   fracis,  ptend%q, dpdry)

   end if

end subroutine sas_conv_tend_2

end module sas_conv_intr
