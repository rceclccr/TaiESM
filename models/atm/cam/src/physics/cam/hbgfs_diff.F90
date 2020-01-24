module hbgfs_diff

! ----------------------------------------
! ycwang Nov 2014 Port the code to CAM5.3
! After counseling with Chao-An, add some comments and remove several arguments.
! ----------------------------------------
!reivsed on 2013/01/28
!joanne chen Nov 2012 modified to accommadte H&P 2011 PBL scheme
!This subroutine is based on the structure of HB scheme,
!    the evaluation of PBL height, diffusivity of local and non-local scheme are also
!    the same as HB scheme, 
!    only plus cloud-top-driven diffusivity into local and
!    non-local diffusivity.
!The calculation of cloud-top-driven diffusivity is in subroutine austausch_atm,
!    and combined with local diffusivity in subroutine austausch_atm, and non-local
!    diffusivity in subroutine austausch_pbl.
!    Furthermore, vertical-diffusion.F90 , eddy_diff.F90 and diffusion_solver.F90
!    also need to be modified due to the additional cloud-top-driven diffusivity.

  !---------------------------------------------------------------------------------
  ! Module to compute mixing coefficients associated with turbulence in the 
  ! planetary boundary layer and elsewhere.  PBL coefficients are based on Holtslag 
  ! and Boville, 1991.
  !
  ! Public interfaces:
  !    init_hb_diff     initializes time independent coefficients
  !    compute_hb_diff  computes eddy diffusivities and counter-gradient fluxes
  !
  ! Private methods:
  !       trbintd         initializes time dependent variables
  !       pblintd         initializes time dependent variables that depend pbl depth
  !       austausch_atm   computes free atmosphere exchange coefficients
  !       austausch_pbl   computes pbl exchange coefficients
  !
  !---------------------------Code history--------------------------------
  ! Standardized:      J. Rosinski, June 1992
  ! Reviewed:          P. Rasch, B. Boville, August 1992
  ! Reviewed:          P. Rasch, April 1996
  ! Reviewed:          B. Boville, April 1996
  ! rewritten:         B. Boville, May 2000
  ! rewritten:         B. Stevens, August 2000
  ! modularized:       J. McCaa, September 2004
  !---------------------------------------------------------------------------------
  use shr_kind_mod, only: r8 => shr_kind_r8
  use spmd_utils,   only: masterproc          ! output from hb_init should be eliminated
  use ppgrid,       only: pver, pverp, pcols  ! these should be passed in
  use cam_logfile,  only: iulog
  !joanne use radiation to calculate solar zenith angle
  use constituents, only: pcnst
  use radiation
  !joanne
  ! Yi-Chi
  !use hb_diff, only:pblintd
  ! end Yi-Chi

  implicit none
  private
  save

  ! +++ Yi-Chi +++
  ! Public interfaces
  public init_hbgfs_diff
  public compute_hbgfs_diff
  !public init_hb_diff
  !public compute_hb_diff
  public pblintd
  ! ---------Yi-Chi---
  !
  ! PBL limits
  !
  real(r8), parameter :: pblmaxp   = 4.e4_r8        ! pbl max depth in pressure units
  real(r8), parameter :: zkmin     = 0.01_r8        ! Minimum kneutral*f(ri)
  !
  ! PBL Parameters
  !
  real(r8), parameter :: onet  = 1._r8/3._r8 ! 1/3 power in wind gradient expression
  real(r8), parameter :: betam = 15.0_r8  ! Constant in wind gradient expression
  real(r8), parameter :: betas =  5.0_r8  ! Constant in surface layer gradient expression
  real(r8), parameter :: betah = 15.0_r8  ! Constant in temperature gradient expression 
  real(r8), parameter :: fakn  =  7.2_r8  ! Constant in turbulent prandtl number
  real(r8), parameter :: fak   =  8.5_r8  ! Constant in surface temperature excess         
  real(r8), parameter :: ricr  =  0.3_r8  ! Critical richardson number
  real(r8), parameter :: sffrac=  0.1_r8  ! Surface layer fraction of boundary layer
  real(r8), parameter :: binm  = betam*sffrac       ! betam * sffrac
  real(r8), parameter :: binh  = betah*sffrac       ! betah * sffrac

  ! Pbl constants set using values from other parts of code

  real(r8) :: cpair      ! Specific heat of dry air
  real(r8) :: g          ! Gravitational acceleration
  real(r8) :: ml2(pverp) ! Mixing lengths squared
  real(r8) :: vk         ! Von Karman's constant
  real(r8) :: ccon       ! fak * sffrac * vk

  integer :: npbl       ! Maximum number of levels in pbl from surface
  integer :: ntop_turb  ! Top level to which turbulent vertical diffusion is applied.
  integer :: nbot_turb  ! Bottom level to which turbulent vertical diff is applied.

!===============================================================================
CONTAINS
!===============================================================================

subroutine init_hbgfs_diff(gravx, cpairx, ntop_eddy, nbot_eddy, pref_mid, &
                        vkx, eddy_scheme)

   !----------------------------------------------------------------------- 
   ! 
   ! Initialize time independent variables of turbulence/pbl package.
   ! 
   !-----------------------------------------------------------------------

   !------------------------------Arguments--------------------------------
   real(r8), intent(in) :: gravx     ! acceleration of gravity
   real(r8), intent(in) :: cpairx    ! specific heat of dry air
   real(r8), intent(in) :: pref_mid(pver)! reference pressures at midpoints
   real(r8), intent(in) :: vkx       ! Von Karman's constant
   integer, intent(in)  :: ntop_eddy ! Top level to which eddy vert diff is applied.
   integer, intent(in)  :: nbot_eddy ! Bottom level to which eddy vert diff is applied.
   character(len=16),  intent(in) :: eddy_scheme

   !---------------------------Local workspace-----------------------------
   integer :: k                     ! vertical loop index
   !-----------------------------------------------------------------------

   ! Basic constants
   cpair = cpairx
   g     = gravx
   vk    = vkx
   ccon  = fak*sffrac*vk
   ntop_turb = ntop_eddy
   nbot_turb = nbot_eddy

   ! Set the square of the mixing lengths.
   ml2(ntop_turb) = 0._r8
   do k = ntop_turb+1, nbot_turb
      ml2(k) = 30.0_r8**2                 ! HB scheme: length scale = 30m  
      if  ( eddy_scheme .eq. 'HBR' ) then      
         ml2(k) = 1.0_r8**2               ! HBR scheme: length scale = 1m  
      end if
   end do
   ml2(nbot_turb+1) = 0._r8

   ! Limit pbl height to regions below 400 mb
   ! npbl = max number of levels (from bottom) in pbl

   npbl = 0
   do k=nbot_turb,ntop_turb,-1
      if (pref_mid(k) >= pblmaxp) then
         npbl = npbl + 1
      end if
   end do
   npbl = max(npbl,1)

   if (masterproc) then
   ! +++ Yi-Chi
      write(iulog,*)'INIT_HBGFS_DIFF: PBL height will be limited to bottom ',npbl, &
   ! write(iulog,*)'INIT_HB_DIFF: PBL height will be limited to bottom ',npbl, &
   ! ---Yi-Chi
         ' model levels. Top is ',pref_mid(pverp-npbl),' pascals'
   end if

end subroutine init_hbgfs_diff

!===============================================================================
!joanne add pint to calculate cloud-top cooling
  subroutine compute_hbgfs_diff(lchnk, ncol, pint,           &
!joanne
! subroutine compute_hb_diff(lchnk, ncol,           &
       th      ,t       ,q       ,z       ,zi      , &
       pmid    ,u       ,v       ,taux    ,tauy    , &
       shflx   ,qflx    ,obklen  ,ustar   ,pblh    , &
       kvm     ,kvh     ,kvq     ,cgh     ,cgs     , &
       tpert   ,qpert   ,cldn    ,ocnfrac ,tke     , &
       ri      , &
!joanne add swh, hlw to calculate cloud-top cooling, and other diagnostic variables
       !eddy_scheme )
       eddy_scheme ,swh, hlw, kbfs, khfs, kqfs, kvf, &
       ksct, kscu, dktx, dkux)
!joanne

    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: 
    !  Interface routines for calcualtion and diatnostics of turbulence related
    !  coefficients
    !
    ! Author: B. Stevens (rewrite August 2000)
    ! 
    !-----------------------------------------------------------------------

    use pbl_utils, only: virtem, calc_ustar, calc_obklen

    !------------------------------Arguments--------------------------------
    !
    ! Input arguments
    !
    integer, intent(in) :: lchnk                      ! chunk index (for debug only)
    integer, intent(in) :: ncol                       ! number of atmospheric columns
    ! +++ Yi-Chi : add for CAM5.3
    real(r8), intent(in) :: pint(pcols,pverp)
    ! +++ Yi-Chi

    real(r8), intent(in)  :: th(pcols,pver)           ! potential temperature [K]
    real(r8), intent(in)  :: t(pcols,pver)            ! temperature (used for density)
    ! +++ Yi-Chi
    real(r8), intent(in)  :: q(pcols,pver,pcnst)            ! specific humidity [kg/kg]
    !real(r8), intent(in)  :: q(pcols,pver)            ! specific humidity [kg/kg]
    ! --- Yi-Chi
    real(r8), intent(in)  :: z(pcols,pver)            ! height above surface [m]
    real(r8), intent(in)  :: zi(pcols,pverp)          ! height above surface [m]
    real(r8), intent(in)  :: u(pcols,pver)            ! zonal velocity
    real(r8), intent(in)  :: v(pcols,pver)            ! meridional velocity
    real(r8), intent(in)  :: taux(pcols)              ! zonal stress [N/m2]
    real(r8), intent(in)  :: tauy(pcols)              ! meridional stress [N/m2]
    real(r8), intent(in)  :: shflx(pcols)             ! sensible heat flux
    real(r8), intent(in)  :: qflx(pcols)              ! water vapor flux
    real(r8), intent(in)  :: pmid(pcols,pver)         ! midpoint pressures
    real(r8), intent(in)  :: cldn(pcols,pver)         ! new cloud fraction
    real(r8), intent(in)  :: ocnfrac(pcols)           ! Land fraction
!joanne add swh, hlw to calculate cloud-top-driven diffusivity, and other variables used in GFS
    real(r8), intent(in)  :: swh(pcols,pverp)                                 ! short wave heating
    real(r8), intent(in)  :: hlw(pcols,pverp)                                 ! long wave heating
    real(r8), intent(out) :: dkux(pcols,pverp)       ! GFS
    real(r8), intent(out) :: dktx(pcols,pverp)       ! GFS
    !real(r8) :: kqfs(pcols,pcnst)                     ! kinematic surf constituent flux (kg/m2/s)
    real(r8) :: cku(pcols,pverp)       ! GFS
    real(r8) :: ckt(pcols,pverp)       ! GFS
!joanne
    character(len=16), intent(in) :: eddy_scheme

    !
    ! Output arguments
    !
    real(r8), intent(out) :: kvm(pcols,pverp)         ! eddy diffusivity for momentum [m2/s]
    real(r8), intent(out) :: kvh(pcols,pverp)         ! eddy diffusivity for heat [m2/s]
    real(r8), intent(out) :: kvq(pcols,pverp)         ! eddy diffusivity for constituents [m2/s]
    real(r8), intent(out) :: cgh(pcols,pverp)         ! counter-gradient term for heat [J/kg/m]
    real(r8), intent(out) :: cgs(pcols,pverp)         ! counter-gradient star (cg/flux)
    real(r8), intent(out) :: tpert(pcols)             ! convective temperature excess
    real(r8), intent(out) :: qpert(pcols)             ! convective humidity excess
    real(r8), intent(out) :: ustar(pcols)             ! surface friction velocity [m/s]
    real(r8), intent(out) :: obklen(pcols)            ! Obukhov length
    real(r8), intent(out) :: pblh(pcols)              ! boundary-layer height [m]
    real(r8), intent(out) :: tke(pcols,pverp)         ! turbulent kinetic energy (estimated)
    real(r8), intent(out) :: ri(pcols,pver)           ! richardson number: n2/s2
    !
    !---------------------------Local workspace-----------------------------
    !
    real(r8) :: thv(pcols,pver)         ! virtual temperature
    real(r8) :: rrho(pcols)             ! 1./bottom level density
    real(r8) :: wstar(pcols)            ! convective velocity scale [m/s]
    real(r8) :: kqfs(pcols)             ! kinematic surf constituent flux (kg/m2/s)
    real(r8) :: khfs(pcols)             ! kinimatic surface heat flux 
    real(r8) :: kbfs(pcols)             ! surface buoyancy flux 
    real(r8) :: kvf(pcols,pverp)        ! free atmospheric eddy diffsvty [m2/s]
    real(r8) :: s2(pcols,pver)          ! shear squared
    real(r8) :: n2(pcols,pver)          ! brunt vaisaila frequency
    real(r8) :: bge(pcols)              ! buoyancy gradient enhancment
    integer  :: ktopbl(pcols)           ! index of first midpoint inside pbl

    !joanne add for GFS
    real(r8) :: ksct(pcols,pverp)        ! check GFS diff.[m2/s]
    real(r8) :: kscu(pcols,pverp)        ! check GFS diff.[m2/s]
    !real(r8) :: ksft(pcols,pverp)        ! check GFS diff. [m2/s]
    !real(r8) :: ksfu(pcols,pverp)        ! check GFS diff. [m2/s]
    !joanne

    !
    ! Initialize time dependent variables that do not depend on pbl height
    !

    ! virtual temperature
    ! +++ Yi-Chi
    thv(:ncol,ntop_turb:) = virtem(th(:ncol,ntop_turb:),q(:ncol,ntop_turb:,1))
    !thv(:ncol,ntop_turb:) = virtem(th(:ncol,ntop_turb:),q(:ncol,ntop_turb:))
    ! --- Yi-Chi

    ! Compute ustar, Obukhov length, and kinematic surface fluxes.
    call calc_ustar(t(:ncol,pver),pmid(:ncol,pver),taux(:ncol),tauy(:ncol), &
         rrho(:ncol),ustar(:ncol))
    call calc_obklen(th(:ncol,pver), thv(:ncol,pver), qflx(:ncol),  &
                     shflx(:ncol),   rrho(:ncol),     ustar(:ncol), &
                     khfs(:ncol),    kqfs(:ncol),     kbfs(:ncol),  &
                     obklen(:ncol))
    ! Calculate s2, n2, and Richardson number.
    call trbintd(ncol    ,                            &
         thv     ,z       ,u       ,v       , &
         s2      ,n2      ,ri      )
    !
    ! Initialize time dependent variables that do depend on pbl height
    !
    ! +++ Yi-Chi : Nov 2014
    !write (iulog,*) 'Yi-Chi:(hbgfs_diff) npbl',npbl
    call  pblintd(ncol    ,                            &
         thv     ,z       ,u       ,v       , &
         ustar   ,obklen  ,kbfs    ,pblh    ,wstar   , &
         zi      ,cldn    ,ocnfrac ,bge     )
    !
    ! Get free atmosphere exchange coefficients
    !
!joanne add variables for GFS
    call austausch_atm(lchnk, ncol ,pint ,ri,s2 ,n2 ,z ,zi ,kvf ,&
                       th,q,t ,swh, hlw , u, v, &
                       ksct, kscu, ckt,cku)
!    call austausch_atm(ncol    ,ri      ,s2      ,kvf     )
!end joanne
    ! 
    ! Get pbl exchange coefficients
    !
!joanne add variables for GFS
    call austausch_pbl(lchnk, ncol,                    &
         z       ,kvf     ,kqfs    ,khfs    ,kbfs    , &
         obklen  ,ustar   ,wstar   ,pblh    ,kvm     , &
         kvh     ,cgh     ,cgs     ,tpert   ,qpert   , &
         ktopbl  ,tke     ,bge     ,eddy_scheme, &
!         ksct, kscu, dktx, dkux, ckt, cku)
         dktx, dkux, ckt, cku) ! Yi-Chi removes ksct, kscu
!    call austausch_pbl(lchnk, ncol,                    &
!         z       ,kvf     ,kqfs    ,khfs    ,kbfs    , &
!         obklen  ,ustar   ,wstar   ,pblh    ,kvm     , &
!         kvh     ,cgh     ,cgs     ,tpert   ,qpert   , &
!         ktopbl  ,tke     ,bge     ,eddy_scheme)
!end joanne
    !

    kvq(:ncol,:) = kvh(:ncol,:)

    return
  end subroutine compute_hbgfs_diff
  !
  !===============================================================================
  subroutine trbintd(ncol    ,                            &
       thv     ,z       ,u       ,v       , &
       s2      ,n2      ,ri      )

    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: 
    !  Time dependent initialization
    ! 
    ! Method: 
    !  Diagnosis of variables that do not depend on mixing assumptions or
    !  PBL depth
    !
    ! Author: B. Stevens (extracted from pbldiff, August, 2000)
    ! 
    !-----------------------------------------------------------------------
    !------------------------------Arguments--------------------------------
    !
    ! Input arguments
    !
    integer, intent(in) :: ncol                      ! number of atmospheric columns

    real(r8), intent(in)  :: thv(pcols,pver)         ! virtual temperature
    real(r8), intent(in)  :: z(pcols,pver)           ! height above surface [m]
    real(r8), intent(in)  :: u(pcols,pver)           ! windspeed x-direction [m/s]
    real(r8), intent(in)  :: v(pcols,pver)           ! windspeed y-direction [m/s]

    !
    ! Output arguments
    !
    real(r8), intent(out) :: s2(pcols,pver)          ! shear squared
    real(r8), intent(out) :: n2(pcols,pver)          ! brunt vaisaila frequency
    real(r8), intent(out) :: ri(pcols,pver)          ! richardson number: n2/s2
    !
    !---------------------------Local workspace-----------------------------
    !
    integer  :: i                        ! longitude index
    integer  :: k                        ! level index

    real(r8) :: vvk                      ! velocity magnitude squared
    real(r8) :: dvdz2                    ! velocity shear squared
    real(r8) :: dz                       ! delta z between midpoints
    !
    ! Compute shear squared (s2), brunt vaisaila frequency (n2) and related richardson
    ! number (ri). Use virtual temperature to compute n2.
    !

    do k=ntop_turb,nbot_turb-1
       do i=1,ncol
          dvdz2   = (u(i,k)-u(i,k+1))**2 + (v(i,k)-v(i,k+1))**2
          dvdz2   = max(dvdz2,1.e-30_r8)
          dz      = z(i,k) - z(i,k+1)
          s2(i,k) = dvdz2/(dz**2)
          n2(i,k) = g*2.0_r8*( thv(i,k) - thv(i,k+1))/((thv(i,k) + thv(i,k+1))*dz)
          ri(i,k) = n2(i,k)/s2(i,k)
       end do
    end do

    return
  end subroutine trbintd
  !
  !===============================================================================
  ! +++ Yi-Chi : remove repetitive subroutines
       
  subroutine pblintd(ncol    ,                            &
       thv     ,z       ,u       ,v       , &
       ustar   ,obklen  ,kbfs    ,pblh    ,wstar   , &
       zi      ,cldn    ,ocnfrac ,bge     )
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: 
    ! Diagnose standard PBL variables
    ! 
    ! Method: 
    ! Diagnosis of PBL depth and related variables.  In this case only wstar.
    ! The PBL depth follows:
    !    Holtslag, A.A.M., and B.A. Boville, 1993:
    !    Local versus Nonlocal Boundary-Layer Diffusion in a Global Climate
    !    Model. J. Clim., vol. 6., p. 1825--1842.
    !
    ! Updated by Holtslag and Hack to exclude the surface layer from the
    ! definition of the boundary layer Richardson number. Ri is now defined
    ! across the outer layer of the pbl (between the top of the surface
    ! layer and the pbl top) instead of the full pbl (between the surface and
    ! the pbl top). For simiplicity, the surface layer is assumed to be the
    ! region below the first model level (otherwise the boundary layer depth
    ! determination would require iteration).
    !
    ! Modified for boundary layer height diagnosis: Bert Holtslag, june 1994
    ! >>>>>>>>>  (Use ricr = 0.3 in this formulation)
    ! 
    ! Author: B. Stevens (extracted from pbldiff, August 2000)
    ! 
    !-----------------------------------------------------------------------
    !------------------------------Arguments--------------------------------
    !
    ! Input arguments
    !
    integer, intent(in) :: ncol                      ! number of atmospheric columns

    real(r8), intent(in)  :: thv(pcols,pver)         ! virtual temperature
    real(r8), intent(in)  :: z(pcols,pver)           ! height above surface [m]
    real(r8), intent(in)  :: u(pcols,pver)           ! windspeed x-direction [m/s]
    real(r8), intent(in)  :: v(pcols,pver)           ! windspeed y-direction [m/s]
    real(r8), intent(in)  :: ustar(pcols)            ! surface friction velocity [m/s]
    real(r8), intent(in)  :: obklen(pcols)           ! Obukhov length
    real(r8), intent(in)  :: kbfs(pcols)             ! sfc kinematic buoyancy flux [m^2/s^3]
    real(r8), intent(in)  :: zi(pcols,pverp)         ! height above surface [m]
    real(r8), intent(in)  :: cldn(pcols,pver)        ! new cloud fraction
    real(r8), intent(in)  :: ocnfrac(pcols)          ! Land fraction

    !
    ! Output arguments
    !
    real(r8), intent(out) :: wstar(pcols)            ! convective sclae velocity [m/s]
    real(r8), intent(out) :: pblh(pcols)             ! boundary-layer height [m]
    real(r8), intent(out) :: bge(pcols)              ! buoyancy gradient enhancment
    !
    !---------------------------Local parameters----------------------------
    !
    real(r8), parameter   :: tiny = 1.e-30_r8           ! lower bound for wind magnitude
    real(r8), parameter   :: fac  = 100._r8             ! ustar parameter in height diagnosis
    !
    !---------------------------Local workspace-----------------------------
    !
    integer  :: i                       ! longitude index
    integer  :: k                       ! level index

    real(r8) :: phiminv(pcols)          ! inverse phi function for momentum
    real(r8) :: phihinv(pcols)          ! inverse phi function for heat
    real(r8) :: rino(pcols,pver)        ! bulk Richardson no. from level to ref lev
    real(r8) :: tlv(pcols)              ! ref. level pot tmp + tmp excess
    real(r8) :: vvk                     ! velocity magnitude squared

    logical  :: unstbl(pcols)           ! pts w/unstbl pbl (positive virtual ht flx)
    logical  :: check(pcols)            ! True=>chk if Richardson no.>critcal
    logical  :: ocncldcheck(pcols)      ! True=>if ocean surface and cloud in lowest layer
    !
    ! Compute Obukhov length virtual temperature flux and various arrays for use later:
    !
    do i=1,ncol
       check(i)     = .true.
       rino(i,pver) = 0.0_r8
       pblh(i)      = z(i,pver)
    end do
    !
    !
    ! PBL height calculation:  Scan upward until the Richardson number between
    ! the first level and the current level exceeds the "critical" value.
    !
    do k=pver-1,pver-npbl+1,-1
       do i=1,ncol
          if (check(i)) then
             vvk = (u(i,k) - u(i,pver))**2 + (v(i,k) - v(i,pver))**2 + fac*ustar(i)**2
             vvk = max(vvk,tiny)
             rino(i,k) = g*(thv(i,k) - thv(i,pver))*(z(i,k)-z(i,pver))/(thv(i,pver)*vvk)
             if (rino(i,k) >= ricr) then
                pblh(i) = z(i,k+1) + (ricr - rino(i,k+1))/(rino(i,k) - rino(i,k+1)) * &
                     (z(i,k) - z(i,k+1))
                check(i) = .false.
             end if
          end if
       end do
    end do
    !
    ! Estimate an effective surface temperature to account for surface fluctuations
    !
    do i=1,ncol
       !write (iulog,*) 'Yi-Chi:(hbgfs_diff) npbl',npbl
       !write (iulog,*) 'Yi-Chi:(hbgfs_diff) pverp',pverp
       if (check(i)) pblh(i) = z(i,pverp-npbl)
       unstbl(i) = (kbfs(i) > 0._r8)
       check(i)  = (kbfs(i) > 0._r8)
       if (check(i)) then
          phiminv(i)   = (1._r8 - binm*pblh(i)/obklen(i))**onet
          rino(i,pver) = 0.0_r8
          tlv(i)       = thv(i,pver) + kbfs(i)*fak/( ustar(i)*phiminv(i) )
       end if
    end do
    !
    ! Improve pblh estimate for unstable conditions using the convective temperature excess:
    !
    do i = 1,ncol
       bge(i) = 1.e-8_r8
    end do
    do k=pver-1,pver-npbl+1,-1
       do i=1,ncol
          if (check(i)) then
             vvk = (u(i,k) - u(i,pver))**2 + (v(i,k) - v(i,pver))**2 + fac*ustar(i)**2
             vvk = max(vvk,tiny)
             rino(i,k) = g*(thv(i,k) - tlv(i))*(z(i,k)-z(i,pver))/(thv(i,pver)*vvk)
             if (rino(i,k) >= ricr) then
                pblh(i) = z(i,k+1) + (ricr - rino(i,k+1))/(rino(i,k) - rino(i,k+1))* &
                     (z(i,k) - z(i,k+1))
                bge(i) = 2._r8*g/(thv(i,k)+thv(i,k+1))*(thv(i,k)-thv(i,k+1))/(z(i,k)-z(i,k+1))*pblh(i)
                if (bge(i).lt.0._r8) then
                   bge(i) = 1.e-8_r8
                endif
                check(i) = .false.
             end if
          end if
       end do
    end do
    !
    ! PBL height must be greater than some minimum mechanical mixing depth
    ! Several investigators have proposed minimum mechanical mixing depth
    ! relationships as a function of the local friction velocity, u*.  We
    ! make use of a linear relationship of the form h = c u* where c=700.
    ! The scaling arguments that give rise to this relationship most often
    ! represent the coefficient c as some constant over the local coriolis
    ! parameter.  Here we make use of the experimental results of Koracin
    ! and Berkowicz (1988) [BLM, Vol 43] for wich they recommend 0.07/f
    ! where f was evaluated at 39.5 N and 52 N.  Thus we use a typical mid
    ! latitude value for f so that c = 0.07/f = 700.  Also, do not allow
    ! PBL to exceed some maximum (npbl) number of allowable points
    !
    do i=1,ncol
       if (check(i)) pblh(i) = z(i,pverp-npbl)
       pblh(i) = max(pblh(i),700.0_r8*ustar(i))
       wstar(i) = (max(0._r8,kbfs(i))*g*pblh(i)/thv(i,pver))**onet
    end do
    !
    ! Final requirement on PBL heightis that it must be greater than the depth
    ! of the lowest model level over ocean if there is any cloud diagnosed in
    ! the lowest model level.  This is to deal with the inadequacies of the
    ! current "dry" formulation of the boundary layer, where this test is
    ! used to identify circumstances where there is marine stratus in the
    ! lowest level, and to provide a weak ventilation of the layer to avoid
    ! a pathology in the cloud scheme (locking in low-level stratiform cloud)
    ! If over an ocean surface, and any cloud is diagnosed in the
    ! lowest level, set pblh to 50 meters higher than top interface of lowest level
    !
    !  jrm This is being applied everywhere (not just ocean)!
    do i=1,ncol
       ocncldcheck(i) = .false.
       if (cldn(i,pver).ge.0.0_r8) ocncldcheck(i) = .true.
       if (ocncldcheck(i)) pblh(i) = max(pblh(i),zi(i,pver) + 50._r8)
    end do
    !
    return
  end subroutine pblintd
  ! end Yi-Chi
  !
  !===============================================================================
!joanne add variables for GFS
  subroutine austausch_atm(lchnk, ncol   ,pint ,ri,s2,n2, z ,zi ,kvf , &
                           th, q ,t ,swh, hlw ,u ,v , &
                           ksct, kscu, ckt, cku)
! subroutine austausch_atm(ncol    ,ri      ,s2      ,kvf     )
!end joanne
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: 
    !  Computes exchange coefficients for free turbulent flows. 
    ! 
    ! Method: 
    !
    ! The free atmosphere diffusivities are based on standard mixing length
    ! forms for the neutral diffusivity multiplied by functns of Richardson
    ! number. K = l^2 * |dV/dz| * f(Ri). The same functions are used for
    ! momentum, potential temperature, and constitutents. 
    !
    ! The stable Richardson num function (Ri>0) is taken from Holtslag and
    ! Beljaars (1989), ECMWF proceedings. f = 1 / (1 + 10*Ri*(1 + 8*Ri))
    ! The unstable Richardson number function (Ri<0) is taken from  CCM1.
    ! f = sqrt(1 - 18*Ri)
    ! 
    ! Author: B. Stevens (rewrite, August 2000)
    ! 
    !-----------------------------------------------------------------------
    !joanne to calculate cloud-top radiative cooling
      use radiation
      use phys_grid,       only: get_rlat_all_p, get_rlon_all_p
      use time_manager, only: get_curr_calday
    !joanne
    !------------------------------Arguments--------------------------------
    !
    ! Input arguments
    !
    integer, intent(in) :: ncol                     ! number of atmospheric columns
!joanne add for GFS
    integer, intent(in) :: lchnk                    ! local chunk index (for debug only)
    real(r8), intent(in)  :: th(pcols,pver)          ! potential temperature [K]
    real(r8), intent(in)  :: q(pcols,pver,pcnst)     ! specific humidity [kg/kg]
    real(r8), intent(in)  :: t(pcols,pver)            ! temperature (used for density)
    real(r8), intent(in)  :: u(pcols,pver)           ! windspeed x-direction [m/s]
    real(r8), intent(in)  :: v(pcols,pver)           ! windspeed y-direction [m/s]
    real(r8),intent(in)   :: n2(pcols,pver)          ! brunt vaisaila frequency
    real(r8), intent(in)  :: zi(pcols,pverp)          ! height above surface [m]
    real(r8), intent(in)  :: z(pcols,pver)            ! height above surface [m]
    real(r8), intent(in) :: pint(pcols,pverp)   ! pressure at layer interfaces (Pa)
    real(r8), intent(in) :: swh(pcols,pverp)                                 ! short wave heating
    real(r8), intent(in) :: hlw(pcols,pverp)                                 ! long wave heating
    real(r8) :: xkzmo(pcols,pverp)       ! GFS
    real(r8) :: xkzo(pcols,pverp)       ! GFS
    real(r8) :: test(pcols)             ! test coszrs
    real(r8) :: coszrs(pcols)           ! coszrs
    logical  :: flg(pcols)           ! GFS
    logical  :: scuflg(pcols)           ! GFS
    real(r8) :: rdzt(pcols,pver)        ! GFS
!joanne

    real(r8), intent(in)  ::  s2(pcols,pver)        ! shear squared
    real(r8), intent(in)  ::  ri(pcols,pver)        ! richardson no
    !
    ! Output arguments
    !
    real(r8), intent(out) :: kvf(pcols,pverp)       ! coefficient for heat and tracers
!joanne add for GFS
    real(r8) :: camkvf(pcols,pverp)       ! coefficient for heat and tracers
!   real(r8) :: dku(pcols,pverp)       ! GFS
!   real(r8) :: dkt(pcols,pverp)       ! GFS
    real(r8), intent(out) :: cku(pcols,pverp)       ! GFS
    real(r8), intent(out) :: ckt(pcols,pverp)       ! GFS
    real(r8),intent(out) :: ksct(pcols,pverp)       ! check GFS diff.[m2/s]
    real(r8),intent(out) :: kscu(pcols,pverp)       ! check GFS diff.[m2/s]
    !real(r8) :: ksft(pcols,pverp)        ! check GFS diff. [m2/s]
    !real(r8) :: ksfu(pcols,pverp)        ! check GFS diff. [m2/s]

    real(r8) :: qlx(pcols,pver)       ! GFS
    real(r8) :: qtx(pcols,pver)       ! GFS
    real(r8) :: thetae(pcols,pverp)       ! GFS
    real(r8) :: thvx(pcols,pverp)       ! GFS
    real(r8) :: thlvx(pcols,pverp)       ! GFS
    real(r8) :: thlvx1(pcols)       ! GFS
    real(r8) :: radx(pcols,pver)       !GFS
    real(r8) :: bf(pcols,pver)        !GFS
!joanne
    !
    !---------------------------Local workspace-----------------------------
    !
    real(r8) :: fofri                  ! f(ri)
    real(r8) :: kvn                    ! neutral Kv
!joanne add for GFS
!
    real(r8) :: ti                     ! GFS dt
    real(r8) :: rdz                     ! GFS rdz
    real(r8) :: dw2                    ! GFS
    real(r8) :: shr2                    ! GFS
    real(r8) :: bvf2                    ! GFS

    real(r8) :: rigfs                     ! GFS
    real(r8) :: zk                     ! GFS
    real(r8) :: rl2                     ! GFS
    real(r8) :: dk                     ! GFS
    real(r8) :: sri                     ! GFS
    real(r8) :: tem2                     ! GFS
    real(r8) :: tem1                     ! GFS
    real(r8) :: tem                     ! GFS
    real(r8) :: prnum                     ! GFS
    real(r8) :: ptem                     ! GFS
    real(r8) :: ptem1                     ! GFS
    real(r8) :: ptem2                     ! GFS

    real(r8) :: cteit(pcols)           ! GFS
    real(r8) :: rent(pcols)           ! GFS
    real(r8) :: hrad(pcols)           ! GFS
    real(r8) :: zd(pcols)           ! GFS
    real(r8) :: zdd(pcols)           ! GFS
    real(r8) :: vrad(pcols)           ! GFS
    real(r8) :: scheck(pcols,pver)           ! GFS
    real(r8) :: radmin(pcols)           ! GFS
    integer :: kcld(pcols)           ! GFS
    integer :: icld(pcols)           ! GFS
    integer :: krad(pcols)           ! GFS
    integer :: lcld(pcols)           ! GFS

    integer :: pver2gfs                 !GFS

    real(r8) :: calday       ! current calendar day
    real(r8) :: clat(pcols)                   ! current latitudes(radians)
    real(r8) :: clon(pcols)                   ! current longitudes(radians)

    integer  :: kk                     ! GFS
!
!joanne

    integer  :: i                      ! longitude index
    integer  :: k                      ! vertical index
    !
    !-----------------------------------------------------------------------
    !
    ! The surface diffusivity is always zero
    !
    kvf(:ncol,pverp) = 0.0_r8
    !
    ! Set the vertical diffusion coefficient above the top diffusion level
    ! Note that nbot_turb != pver is not supported
    !
    kvf(:ncol,1:ntop_turb) = 0.0_r8

    !
    !joanne : Main part for GFS cloud-top diffusivity calculation
    ! 1. to calculate Cosine solar zenith angle
    ! 2. 
    !  
    !joanne to calculate Cosine solar zenith angle
      calday = get_curr_calday()
! Cosine solar zenith angle for current time step
      call get_rlat_all_p(lchnk, ncol, clat)
      call get_rlon_all_p(lchnk, ncol, clon)
      call zenith (calday, clat, clon, coszrs, ncol)
!joanne

     do k=pver+1,1,-1
     do i=1,ncol
        xkzo(i,k)=0.
        xkzmo(i,k)=0.
     enddo
     enddo
     do k=pver,1,-1
     do i=1,ncol
        rdzt(i,k)=0.
     enddo
     enddo

!joanne main part of evaluating cloud-top-driven diffusivity "ckt" and "cku",
!and plus them in local diffusivity kvf.
      do k = pver,2,-1
        do i=1,ncol
          rdzt(i,k) = 1.0 / (z(i,k-1) - z(i,k))
        enddo
      enddo

!   ! vertical background diffusivity
!1.=xkzm
      do k = pver,2,-1
        do i=1,ncol
          tem1      = 1.0 - pint(i,k) / pint(i,pver+1)
          tem1      = tem1 * tem1 * 10.0
          xkzo(i,k) = 1. * min(1.0, exp(-tem1))
        enddo
      enddo
!   !vertical background diffusivity for momentum
!3.=xkzmu
      do k = pver,2,-1
        do i=1,ncol
          ptem = pint(i,k) / pint(i,pver+1)
          if(ptem.ge.0.2) then
            xkzmo(i,k) = 3.
            ptem1 = pint(i,k)
          else
            tem1 = 1.0 - pint(i,k) / ptem1
            tem1 = tem1 * tem1 * 5.0
            xkzmo(i,k) = 3. * min(1.0, exp(-tem1))
          endif
        enddo
      enddo

!   !diffusivity in the inversion layer is set to be xkzminv (m^2/s)
!0.3=xkzminv
      pver2gfs=pver/2
      do k = pver,pver2gfs,-1
        do i=1,ncol
          if(zi(i,k).gt.250.) then
            tem1 = (t(i,k-1)-t(i,k)) * rdzt(i,k)
            if(tem1 .gt. 1.e-5) then
               xkzo(i,k) = min(xkzo(i,k),0.3)
            endif
          endif
        enddo
      enddo

!1.e-12=qlmin,1.e-8=qmin

      do k = pver,1,-1
        do i = 1,ncol
          tem=q(i,k,2)+q(i,k,3)
          qlx(i,k)   = max(tem,1.e-12)
          qtx(i,k)   = max(q(i,k,1),1.e-8)+qlx(i,k)
          ptem       = qlx(i,k)
          ptem1      = 2.44*1e+6*max(q(i,k,1),1.e-12)/(1004.*t(i,k))
          thetae(i,k)= th(i,k)*(1.+ptem1)
          thvx(i,k)  = th(i,k)*(1.+0.61_r8*max(q(i,k,1),1.e-8)-ptem)
          ptem2      = th(i,k)-((2.44*1e+6)/1004.)*ptem
          thlvx(i,k) = ptem2*(1.+0.61_r8*qtx(i,k))
        enddo
      enddo

      do k = pver+1,1,-1
        do i = 1,ncol
!         dku(i,k)  = 0.
!         dkt(i,k)  = 0.
          cku(i,k)  = 0.
          ckt(i,k)  = 0.
          !ksfu(i,k)=0.
          !ksft(i,k)=0.
          kscu(i,k)=0.
          ksct(i,k)=0.
        enddo
      enddo
      do k = pver,2,-1
        do i = 1,ncol
          tem       = zi(i,k)-zi(i,k+1)
          radx(i,k) = tem*(swh(i,k)*coszrs(i)+hlw(i,k))
        enddo
      enddo
!      write (iulog,*) 'swh=',swh(:ncol,:pver)
!      write (iulog,*) 'hlw=',hlw(:ncol,:pver)
!0.2=rentf1
      do i=1,ncol
         scuflg(i)=.true.
         cteit(i)=0.
         rent(i)=0.2
         icld(i)=0
         hrad(i)=zi(i,pver+1)
         krad(i)=pver
         kcld(i)=2
         lcld(i)=2
         zd(i)=0.
         radmin(i)=0.
         flg(i)  = scuflg(i)
      enddo
!2500.=zstblmax
      do k = pver,2,-1
        do i=1,ncol
          if(flg(i).and.z(i,k).ge.2500.) then
             lcld(i)=k
             flg(i)=.false.
          endif
      enddo
      enddo

!
!  compute buoyancy flux
!
      do k = pver, 2,-1
      do i = 1, ncol
         bf(i,k) = (thvx(i,k-1)-thvx(i,k))*rdzt(i,k)
      enddo
      enddo
!   look for stratocumulus
!   3.5e-5=qlcr
      do i = 1, ncol
        flg(i)=scuflg(i)
      enddo
      pver2gfs=pver/2
      do k = pver2gfs,pver
      do i = 1, ncol
        if(flg(i).and.k.ge.lcld(i)) then
          if(qlx(i,k).ge.3.5e-5) then
             kcld(i)=k
             flg(i)=.false.
          endif
        endif
      enddo
      enddo
      do i = 1, ncol
        if(scuflg(i).and.kcld(i).eq.2) scuflg(i)=.false.
      enddo


      do i = 1, ncol
        flg(i)=scuflg(i)
      enddo
      do k = pver2gfs,pver
      do i = 1, ncol
        if(flg(i).and.k.ge.kcld(i)) then
          if(qlx(i,k).ge.3.5e-5) then
            if(radx(i,k).lt.radmin(i)) then
              radmin(i)=radx(i,k)
              krad(i)=k
            endif
          else
            flg(i)=.false.
          endif
        endif
      enddo
      enddo
      do i = 1, ncol
        if(scuflg(i).and.krad(i).ge.pver) scuflg(i)=.false.
        if(scuflg(i).and.radmin(i).ge.0.) scuflg(i)=.false.
      enddo

      do i = 1, ncol
        flg(i)=scuflg(i)
      enddo
      do k = pver2gfs,pver-1
      do i = 1, ncol
        if(flg(i).and.k.ge.krad(i)) then
          if(qlx(i,k).ge.3.5e-5) then
            icld(i)=icld(i)+1
          else
            flg(i)=.false.
          endif
        endif
      enddo
      enddo
      do i = 1, ncol
!       if(scuflg(i).and.icld(i).gt.pver) scuflg(i)=.false.
        if(scuflg(i).and.icld(i).lt.1) scuflg(i)=.false.
      enddo
!
      do i = 1, ncol
        if(scuflg(i)) then
           hrad(i) = zi(i,krad(i))
        endif
      enddo
!
      do i = 1, ncol
        if(scuflg(i).and.hrad(i).lt.zi(i,pver)) scuflg(i)=.false.
      enddo
!500.=cldtime
      do i = 1, ncol
        if(scuflg(i)) then
          k    = krad(i)
          tem  = zi(i,k)-zi(i,k+1)
          tem1 = 500.*radmin(i)/tem
          thlvx1(i) = thlvx(i,k)+tem1
!         if(thlvx1(i).gt.thlvx(i,k-1)) scuflg(i)=.false.
        endif
      enddo
!
      do i = 1, ncol
         flg(i)=scuflg(i)
      enddo
      do k = pver2gfs,pver
      do i = 1, ncol
        if(flg(i).and.k.ge.krad(i))then
          if(thlvx1(i).le.thlvx(i,k))then
             tem=zi(i,k)-zi(i,k+1)
             zd(i)=zd(i)+tem
          else
             flg(i)=.false.
          endif
        endif
      enddo
      enddo

      do i = 1, ncol
        if(scuflg(i))then
!         kk = max(1, krad(i)+1-icld(i))
          kk = min(pver, krad(i)-1+icld(i))
          zdd(i) = hrad(i)-zi(i,kk+1)
        endif
      enddo
!1/3=h1=onet
      do i = 1, ncol
        if(scuflg(i))then
          zd(i) = max(zd(i),zdd(i))
          zd(i) = min(zd(i),hrad(i))
          tem   = (g/th(i,pver))*zd(i)*(-radmin(i))
          vrad(i)= tem**onet
        endif
      enddo

!  compute diffusion coefficients for cloud-top driven diffusion
!  if the condition for cloud-top instability is met,
!  increase entrainment flux at cloud top
!2.44*1e+6=hvap,1004.=cp,0.7=actei,1.=rentf2
      do i = 1, ncol
        if(scuflg(i)) then
           k = krad(i)
           tem = thetae(i,k) - thetae(i,k-1)
           tem1 = qtx(i,k) - qtx(i,k-1)
           if (tem.gt.0..and.tem1.gt.0.) then
             cteit(i)= 1004.*tem/((2.44*1e+6)*tem1)
             if(cteit(i).gt.0.7) rent(i) = 1.
           endif
        endif
      enddo
!1e-3=tdzmin
      do i = 1, ncol
        if(scuflg(i)) then
           k = krad(i)
           tem1  = max(bf(i,k),1e-3)
           ckt(i,k) = -rent(i)*radmin(i)/tem1
           cku(i,k) = ckt(i,k)
!           write (iulog,*) 'ckt=',ckt(i,k),'cku=',cku(i,k)
!joanne ksct,kscu is just for analysis, not for evaluation.
           ksct(i,k)=ckt(i,k)
           kscu(i,k)=cku(i,k)
!joanne
        endif
      enddo
!0.85=radfac,0.4=vk,0.=dkmin,1000.=dxmax
      do k = pver, pver2gfs,-1
         do i=1,ncol
            if(scuflg(i).and.k.gt.krad(i)) then
               tem1=hrad(i)-zd(i)
               tem2=zi(i,k)-tem1
               if(tem2.gt.0.) then
                  ptem= tem2/zd(i)
                  if(ptem.ge.1.) ptem= 1.
                  ptem= tem2*ptem*sqrt(1.-ptem)
                  scheck(i,k)=ptem
                  ckt(i,k) = 0.85*0.4*vrad(i)*ptem
                  cku(i,k) = 0.75*ckt(i,k)
                  ckt(i,k) = max(ckt(i,k),0.)
                  ckt(i,k) = min(ckt(i,k),1000.)
                  cku(i,k) = max(cku(i,k),0.)
                  cku(i,k) = min(cku(i,k),1000.)
                  ksct(i,k)=ckt(i,k)
                  kscu(i,k)=cku(i,k)
               endif
            endif
         enddo
      enddo
!joanne finish the process of evaluating cloud-top-driven diffusivity "ckt" and "cku",
!and then plus them in local diffusivity kvf.
    ! joanne end : Main part for GFS cloud-top diffusivity

    !
    ! Compute the free atmosphere vertical diffusion coefficients: kvh = kvq = kvm. 
    !
   do k = ntop_turb, nbot_turb-1
       do i=1,ncol
          if (ri(i,k) < 0.0_r8) then
             fofri = sqrt(max(1._r8 - 18._r8*ri(i,k),0._r8))
          else 
             fofri = 1.0_r8/(1.0_r8 + 10.0_r8*ri(i,k)*(1.0_r8 + 8.0_r8*ri(i,k)))    
          end if
          kvn = ml2(k)*sqrt(s2(i,k))
          kvf(i,k+1) = max(zkmin,kvn*fofri)
       end do
    end do

    return
  end subroutine austausch_atm
  !
  !===============================================================================
!joanne add variables for GFS
  subroutine austausch_pbl(lchnk ,ncol    ,          &
       z       ,kvf     ,kqfs    ,khfs    ,kbfs    , &
       obklen  ,ustar   ,wstar   ,pblh    ,kvm     , &
       kvh     ,cgh     ,cgs     ,tpert   ,qpert   , &
!       ktopbl  ,tke     ,bge     ,eddy_scheme)
       ktopbl  ,tke     ,bge     ,eddy_scheme, &
!       ksct, kscu, dktx, dkux, ckt, cku)
       dktx, dkux, ckt, cku)
!joanne
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: 
    ! Atmospheric Boundary Layer Computation
    ! 
    ! Method: 
    ! Nonlocal scheme that determines eddy diffusivities based on a
    ! specified boundary layer height and a turbulent velocity scale;
    ! also, countergradient effects for heat and moisture, and constituents
    ! are included, along with temperature and humidity perturbations which
    ! measure the strength of convective thermals in the lower part of the
    ! atmospheric boundary layer.
    !
    ! For more information, see Holtslag, A.A.M., and B.A. Boville, 1993:
    ! Local versus Nonlocal Boundary-Layer Diffusion in a Global Climate
    ! Model. J. Clim., vol. 6., p. 1825--1842.
    !
    ! Updated by Holtslag and Hack to exclude the surface layer from the
    ! definition of the boundary layer Richardson number. Ri is now defined
    ! across the outer layer of the pbl (between the top of the surface
    ! layer and the pbl top) instead of the full pbl (between the surface and
    ! the pbl top). For simiplicity, the surface layer is assumed to be the
    ! region below the first model level (otherwise the boundary layer depth
    ! determination would require iteration).
    !
    ! Author: B. Boville, B. Stevens (rewrite August 2000)
    ! 
    !-----------------------------------------------------------------------
!++ debug code to be removed after validation of PBL codes
     use phys_debug, only: phys_debug_hbdiff1
!++ debug code to be removed after validation of PBL codes
    !------------------------------Arguments--------------------------------
    !
    ! Input arguments
    !
    integer, intent(in) :: lchnk                    ! local chunk index (for debug only)
    integer, intent(in) :: ncol                     ! number of atmospheric columns

    real(r8), intent(in) :: z(pcols,pver)           ! height above surface [m]
    real(r8), intent(in) :: kvf(pcols,pverp)        ! free atmospheric eddy diffsvty [m2/s]
    real(r8), intent(in) :: kqfs(pcols)             ! kinematic surf cnstituent flux (kg/m2/s)
    real(r8), intent(in) :: khfs(pcols)             ! kinimatic surface heat flux 
    real(r8), intent(in) :: kbfs(pcols)             ! surface buoyancy flux 
    real(r8), intent(in) :: pblh(pcols)             ! boundary-layer height [m]
    real(r8), intent(in) :: obklen(pcols)           ! Obukhov length
    real(r8), intent(in) :: ustar(pcols)            ! surface friction velocity [m/s]
    real(r8), intent(in) :: wstar(pcols)            ! convective velocity scale [m/s]
    real(r8), intent(in) :: bge(pcols)              ! buoyancy gradient enhancment
    character(len=16), intent(in) :: eddy_scheme

!joanne add variables for GFS
    !real(r8) :: ksct(pcols,pverp)        ! check GFS diff.[m2/s]
   ! real(r8) :: kscu(pcols,pverp)        ! check GFS diff.[m2/s]
    real(r8), intent(in) :: cku(pcols,pverp)       ! GFS
    real(r8), intent(in) :: ckt(pcols,pverp)       ! GFS
!joanne

    !
    ! Output arguments
    !
    real(r8), intent(out) :: kvm(pcols,pverp)       ! eddy diffusivity for momentum [m2/s]
    real(r8), intent(out) :: kvh(pcols,pverp)       ! eddy diffusivity for heat [m2/s]
    real(r8), intent(out) :: cgh(pcols,pverp)       ! counter-gradient term for heat [J/kg/m]
    real(r8), intent(out) :: cgs(pcols,pverp)       ! counter-gradient star (cg/flux)
    real(r8), intent(out) :: tpert(pcols)           ! convective temperature excess
    real(r8), intent(out) :: qpert(pcols)           ! convective humidity excess

    integer,  intent(out) :: ktopbl(pcols)          ! index of first midpoint inside pbl
    real(r8), intent(out) :: tke(pcols,pverp)       ! turbulent kinetic energy (estimated)
!joanne add for store surface diffusivity
    real(r8), intent(out) :: dkux(pcols,pverp)       ! GFS
    real(r8), intent(out) :: dktx(pcols,pverp)       ! GFS
!joanne
    !
    !---------------------------Local workspace-----------------------------
    !
    integer  :: i                       ! longitude index
    integer  :: k                       ! level index

    real(r8) :: phiminv(pcols)          ! inverse phi function for momentum
    real(r8) :: phihinv(pcols)          ! inverse phi function for heat
    real(r8) :: wm(pcols)               ! turbulent velocity scale for momentum
    real(r8) :: zp(pcols)               ! current level height + one level up 
    real(r8) :: fak1(pcols)             ! k*ustar*pblh     
    real(r8) :: fak2(pcols)             ! k*wm*pblh
    real(r8) :: fak3(pcols)             ! fakn*wstar/wm
    real(r8) :: pblk(pcols)             ! level eddy diffusivity for momentum
    real(r8) :: pr(pcols)               ! Prandtl number for eddy diffusivities
    real(r8) :: zl(pcols)               ! zmzp / Obukhov length
    real(r8) :: zh(pcols)               ! zmzp / pblh
    real(r8) :: zzh(pcols)              ! (1-(zmzp/pblh))**2
    real(r8) :: zmzp                    ! level height halfway between zm and zp
    real(r8) :: term                    ! intermediate calculation
    real(r8) :: kve                     ! diffusivity at entrainment layer in unstable cases 

    logical  :: unstbl(pcols)           ! pts w/unstbl pbl (positive virtual ht flx)
    logical  :: pblpt(pcols)            ! pts within pbl
    !
    ! Initialize height independent arrays
    !

    !drb initialize variables for runtime error checking
    kvm = 0._r8	
    kvh = 0._r8
    kve = 0._r8
    cgh = 0._r8
    cgs = 0._r8
    tpert = 0._r8
    qpert = 0._r8
    ktopbl = 0._r8
    tke = 0._r8
    ! Yi-Chi : initialization
    dkux = 0._r8
    dktx = 0._r8
    ! Yi-Chi end.

    do i=1,ncol
       unstbl(i) = (kbfs(i) > 0._r8)
       pblk(i) = 0.0_r8
       fak1(i) = ustar(i)*pblh(i)*vk
       if (unstbl(i)) then
          phiminv(i) = (1._r8 - binm*pblh(i)/obklen(i))**onet
          phihinv(i) = sqrt(1._r8 - binh*pblh(i)/obklen(i))
          wm(i)      = ustar(i)*phiminv(i)
          fak2(i)    = wm(i)*pblh(i)*vk
          fak3(i)    = fakn*wstar(i)/wm(i)
          tpert(i)   = max(khfs(i)*fak/wm(i),0._r8)
          qpert(i)   = max(kqfs(i)*fak/wm(i),0._r8)
       else
          tpert(i)   = max(khfs(i)*fak/ustar(i),0._r8)
          qpert(i)   = max(kqfs(i)*fak/ustar(i),0._r8)
       end if
    end do
    !
    ! Initialize output arrays with free atmosphere values
    !
    do k=1,pverp
       do i=1,ncol
          kvm(i,k) = kvf(i,k)
          kvh(i,k) = kvf(i,k)
          cgh(i,k) = 0.0_r8
          cgs(i,k) = 0.0_r8
       end do
    end do
    !
    ! Main level loop to compute the diffusivities and counter-gradient terms. These terms are 
    ! only calculated at points determined to be in the interior of the pbl (pblpt(i)==.true.),
    ! and then calculations are directed toward regime: stable vs unstable, surface vs outer 
    ! layer.
    !
    do k=pver,pver-npbl+2,-1
       do i=1,ncol
          pblpt(i) = (z(i,k) < pblh(i))
          if (pblpt(i)) then
             ktopbl(i) = k
             zp(i)  = z(i,k-1)
             if (zkmin == 0.0_r8 .and. zp(i) > pblh(i)) zp(i) = pblh(i)
             zmzp    = 0.5_r8*(z(i,k) + zp(i)) ! we think this is an approximation to the interface height (where KVs are calculated)
             zh(i)   = zmzp/pblh(i)
             zl(i)   = zmzp/obklen(i)
             zzh(i)  = zh(i)*max(0._r8,(1._r8 - zh(i)))**2
             if (unstbl(i)) then
                if (zh(i) < sffrac) then
                   term     = (1._r8 - betam*zl(i))**onet
                   pblk(i)  = fak1(i)*zzh(i)*term
                   pr(i)    = term/sqrt(1._r8 - betah*zl(i))
                else
                   pblk(i)  = fak2(i)*zzh(i)
                   pr(i)    = phiminv(i)/phihinv(i) + ccon*fak3(i)/fak
                   cgs(i,k) = fak3(i)/(pblh(i)*wm(i))
                   cgh(i,k) = khfs(i)*cgs(i,k)*cpair
                end if
             else
                if (zl(i) <= 1._r8) then
                   pblk(i) = fak1(i)*zzh(i)/(1._r8 + betas*zl(i))
                else
                   pblk(i) = fak1(i)*zzh(i)/(betas + zl(i))
                end if
                pr(i)    = 1._r8
             end if
!joanne store surface diffusivity in dkux, dktx (are needed in vertical_diffusion.F90 and diffusion_solver.F90 ).
!and combine cloud-top-driven cku, ckt as kvm,kvh
!              dkux(i,k) = max(pblk(i),kvf(i,k))
!              dktx(i,k) = max(pblk(i)/pr(i),kvf(i,k))
!0115
               dkux(i,k) = max(pblk(i),kvf(i,k)-ckt(i,k))
               dktx(i,k) = max(pblk(i)/pr(i),kvf(i,k)-ckt(i,k))
               kvm(i,k) = dkux(i,k)+cku(i,k)
               kvh(i,k) = dktx(i,k)+ckt(i,k)
!0115
!!0115         kvm(i,k) = dkux(i,k)+kscu(i,k)
!!0115         kvh(i,k) = dktx(i,k)+ksct(i,k)
!              kvm(i,k) = max(pblk(i),kvf(i,k))+kscu(i,k)
!              kvh(i,k) = max(pblk(i)/pr(i),kvf(i,k))+ksct(i,k)
! -----default setup ------------
!             kvm(i,k) = max(pblk(i),kvf(i,k))
!             kvh(i,k) = max(pblk(i)/pr(i),kvf(i,k))
! end joanne
          end if
       end do
    end do

!++ debug code to be removed after validation of PBL codes
    call phys_debug_hbdiff1(lchnk, pblh, zl, zh)
!++ debug code to be removed after validation of PBL codes

    !
    ! Check whether last allowed midpoint is within pbl
    !

    if  ( eddy_scheme .eq. 'HBR' ) then  
       ! apply new diffusivity at entrainment zone 
       do i = 1,ncol
          if (bge(i) > 1.e-7_r8) then
             k = ktopbl(i)
             kve = 0.2_r8*(wstar(i)**3+5._r8*ustar(i)**3)/bge(i)
             kvm(i,k) = kve
             kvh(i,k) = kve
          end if
       end do
    end if

    ! Crude estimate of tke (tke=0 above boundary layer)
    do k = max(pverp-npbl,2),pverp
       do i = 1, ncol
          if (z(i,k-1) < pblh(i)) then
             tke(i,k) = ( kvm(i,k) / pblh(i) ) ** 2
          endif
       end do
    end do
    return
  end subroutine austausch_pbl

end module hbgfs_diff
