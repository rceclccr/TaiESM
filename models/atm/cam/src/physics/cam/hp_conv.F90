
   module hp_conv
! Yi-Chi : (June, 2013)
!  - apply modified budget code from Dr.Pan with shallow
! Yi-Chi : (Feb , 2013)
!  - correct the vertical location of icwmr
!    transform qlk from interface level to mid-point qlkh
! Yi-Chi : (Dec 22, 2012) 
!  - correct collocation of cloud liquid
! Yi-Chi (Dec 10, 2012):
!  - correct the evap and heat based on the formulation in PBL scheme
!  - put ix = pcols, im = ncol 
!  - correct indices of tendencies
! Yi-Chi (Dec 7, 2012):
!  - combine ktop(out), kbot(out) and cnt, cnb to save space
!  - check 
! Yi-Chi (Dec 4, 2012):
!  - modify the in-cloud total water mixing ration icwmr
! Yi-Chi (July 19, 2012):
!  - change the calculation of temperature and moisture tendency 
!     to subroutine shalcnv_inv
!  - add momentum transport u and v from shallow
!  - correct the units of cmfmc, qc2, rprdsh, evapcsh
!  - remove redundant variables ud_mf, dt_mf, 
! cjshiu Mar 2012 Modified to accommadte H&P 2011 shallow convection scheme
! Moist convection. Primarily data used by both Zhang-McFarlane convection
! and Hack shallow convective schemes.
!
! $Id$
!
   use shr_kind_mod, only: r8 => shr_kind_r8
   use cam_logfile,  only: iulog
   use spmd_utils,   only: masterproc
   use abortutils,   only: endrun
!cjshiu   use funcphys, only: gpvs,fpvs,fpvsq,fpvsx
!   use funcphys
   implicit none

   private
   save
!
! Public interfaces
!
   public  &
     hpconv_readnl,   & ! read hpconv_nl namelist
     shalcnv_init,    & !  Initialization of data for HP moist shallow convection
     shalcnv_inv,     &   !  HP Shallow convection interface for vertical transpose
     shalcnv             !  HP Shallow convection

!cjs !
!cjs ! Private data used for Hack shallow convection
!cjs !
!cjs ! Private data used for HP shallow convection
   real(r8), parameter :: unset_r8 = huge(1.0_r8)

  ! Namelist variables
   real(r8) :: hpconv_terr = unset_r8    
   real(r8) :: hpconv_c0 = unset_r8 
   real(r8) :: hpconv_c1 = unset_r8 

!cjshiu add
   real(r8) :: terr     
   real(r8) :: c0  
   real(r8) :: c1 
   
!cjshiu  for H&P 2011
   real(r8) :: grav   ! add descrip
   real(r8) :: cp   
   real(r8) :: hvap 
   real(r8) :: rv   
   real(r8) :: fv   
   real(r8) :: t0c  
   real(r8) :: rd   
   real(r8) :: cvap 
   real(r8) :: cliq 
   real(r8) :: eps  
   real(r8) :: epsm1 

contains

subroutine hpconv_readnl(nlfile)

   use namelist_utils,  only: find_group_name
   use units,           only: getunit, freeunit
   use mpishorthand

   character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

   ! Local variables
   integer :: unitn, ierr
   character(len=*), parameter :: subname = 'hpconv_readnl'

!cjshiu need to check if there few parameters needed to be set in namelist for HP
!cjs set  parameter(terr=0.,c0=.002,c1=5.e-4,delta=fv) of H&P2011 to input from namelist
   namelist /hpconv_nl/ hpconv_terr, hpconv_c0, hpconv_c1
   !-----------------------------------------------------------------------------

   if (masterproc) then
      unitn = getunit()
      open( unitn, file=trim(nlfile), status='old' )
      call find_group_name(unitn, 'hpconv_nl', status=ierr)
      if (ierr == 0) then
         read(unitn, hpconv_nl, iostat=ierr)
         if (ierr /= 0) then
            call endrun(subname // ':: ERROR reading namelist')
         end if
      end if
      close(unitn)
      call freeunit(unitn)

      ! set local variables
      terr   = hpconv_terr
      c0     = hpconv_c0
      c1     = hpconv_c1

   end if

#ifdef SPMD
   ! Broadcast namelist variables
   !call mpibcast(cmftau,            1, mpir8,  0, mpicom)
   call mpibcast(c0,                1, mpir8,  0, mpicom)
#endif

end subroutine hpconv_readnl

!================================================================================================

subroutine shalcnv_init (con_g , con_cp , con_hvap ,   &
                         con_rv, con_fvirt, con_t0c,   &
                         con_rd, con_cvap, con_cliq,   &
                         con_eps, con_epsm1)
                           

!cjshiu modified to initialization of H&P2011
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Initialize moist convective mass flux procedure common block, cmfmca
! 
! Method: 
! <Describe the algorithm(s) used in the routine.> 
! <Also include any applicable external references.> 
! 
! Author: J. Hack
! 
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use dycore, only: dycore_is, get_resolution
   use spmd_utils, only: masterproc
!------------------------------Arguments--------------------------------
!
! Input arguments
!
  
   real(r8), intent(in) :: con_g            ! description
   real(r8), intent(in) :: con_cp            ! description
   real(r8), intent(in) :: con_hvap            ! description
   real(r8), intent(in) :: con_rv            ! description
   real(r8), intent(in) :: con_fvirt            ! description
   real(r8), intent(in) :: con_t0c            ! description
   real(r8), intent(in) :: con_rd            ! description
   real(r8), intent(in) :: con_cvap            ! description
   real(r8), intent(in) :: con_cliq            ! description
   real(r8), intent(in) :: con_eps            ! description
   real(r8), intent(in) :: con_epsm1            ! description

!cjs   real(r8), intent(in) :: rair              ! gas constant for dry air
!cjs   real(r8), intent(in) :: cpair             ! specific heat of dry air
!cjs   real(r8), intent(in) :: gravit            ! acceleration due to gravity
!cjs   real(r8), intent(in) :: latvap            ! latent heat of vaporization
!cjs   real(r8), intent(in) :: rhowtr            ! density of liquid water (STP)
!cjs   integer,  intent(in) :: limcnv_in         ! top interface level limit for convection

!cjs   ! local variables
!cjs   character(len=32)    :: hgrid             ! horizontal grid specifier
!
!-----------------------------------------------------------------------
!
! Initialize physical constants for moist convective mass flux procedure
!

   grav = con_g
   cp   = con_cp 
   hvap = con_hvap
   rv   = con_rv 
   fv   = con_fvirt 
   t0c  = con_t0c
   rd   = con_rd 
   cvap = con_cvap 
   cliq = con_cliq
   eps  = con_eps
   epsm1 = con_epsm1

!cjs   cp     = cpair         ! specific heat of dry air
!cjs   hlat   = latvap        ! latent heat of vaporization
!cjs   grav   = gravit        ! gravitational constant
!cjs   rgas   = rair          ! gas constant for dry air
!cjs   rhoh2o = rhowtr        ! density of liquid water (STP)

end subroutine shalcnv_init

  subroutine shalcnv_inv( lchnk   , im , ix , km , pcnst , jcap_inv ,     &
                         delt_inv , del_inv , prsl_inv , ps_inv ,         &
                         zm_inv   , ql_inv  ,                             &
                         q1_inv   , t1_inv  , u1_inv ,v1_inv ,            & 
! +++ Yi-Chi
                        rcs_inv , rn_inv , kcnv_inv , slimsk_inv , &
!                         rcs_inv , rn_inv , cnb_inv , cnt_inv ,       &
! --- Yi-Chi
                         dot_inv , ncloud_inv , hpbl_inv ,              &
                         heat_inv , evap_inv , &
!                         heat_inv , evap_inv , ud_mf_inv , dt_mf_inv,   & 
                         sten_inv , qten_inv , uten_inv , vten_inv ,& 
                         cmfmc_inv, rprdsh_inv , cmfsl_inv , cmflq_inv ,&
                         qc_inv , icwmr_inv , rliq_inv ,                &
! +++ Yi-Chi
!                         evapcv_inv)
                         cnb_inv , cnt_inv, evapcv_inv)
! --- Yi-Chi

!cjshiu more or less following method of uwshcn to do swap of vertical coordinate

    use shr_kind_mod, only: r8 => shr_kind_r8
    use physconst,       only: gravit
! +++ YiChi
!    use funcphys, only: gfuncphys
    use module_gfs_funcphys, only: gfuncphys
!---- YiChi
    implicit none

    integer , intent(in)  :: lchnk                !later need to add description and unit  
    integer , intent(in)  :: im
    integer , intent(in)  :: ix
    integer , intent(in)  :: km
    integer , intent(in)  :: pcnst
    integer , intent(in)  :: jcap_inv
    integer , intent(in)  :: ncloud_inv
    ! +++ Yi-Chi : passed to cnt and cnb
    !integer , intent(out)  :: kbot_inv(im)
    !integer , intent(out)  :: ktop_inv(im)
    ! --- Yi-Chi
    integer , intent(in)  :: kcnv_inv(im)

    real(r8), intent(in)  :: delt_inv
    real(r8), intent(in)  :: del_inv(ix,km)
    real(r8), intent(in)  :: prsl_inv(ix,km)
    real(r8), intent(in)  :: ps_inv(im)
    real(r8), intent(in)  :: zm_inv(ix,km)
!cjshiu need to confirm the 3rd dimension
    !real(r8), intent(in)  :: ql_inv(ix,km,2)
    ! Yi-Chi
    real(r8), intent(in)  :: ql_inv(ix,km)
    ! ---
    real(r8), intent(in)  :: q1_inv(ix,km)
    real(r8), intent(in)  :: t1_inv(ix,km)
    real(r8), intent(in)  :: u1_inv(ix,km)
    real(r8), intent(in)  :: v1_inv(ix,km)
    real(r8), intent(in)  :: rcs_inv(im)
    real(r8), intent(out)  :: rn_inv(im)
!cjshiu    real(r8), intent(out)  :: rn_inv(ix)
    real(r8), intent(in)  :: slimsk_inv(im)
    real(r8), intent(in)  :: dot_inv(ix,km)
    real(r8), intent(in)  :: hpbl_inv(im)
!cjshiu note the input from CAM is 2-D now assign to 1-D probably need to do a temp parameters to assign value to avoid possible segmentation fault
    real(r8), intent(in)  :: heat_inv(im)
    real(r8), intent(in)  :: evap_inv(im)
!cjshiu need to check the difference between ix and im used for some of the parameters.
! Yi-Chi: remove "in" from tendency
    real(r8), intent(out)  :: uten_inv(ix,km)
    real(r8), intent(out)  :: vten_inv(ix,km)
    real(r8), intent(out)  :: sten_inv(ix,km)
    real(r8), intent(out)  :: qten_inv(ix,km,pcnst)
!cjshiu    real(r8), intent(inout)  :: qten_inv(ix,km,2)
    real(r8), intent(out)  :: cmfmc_inv(ix,km)
    real(r8), intent(out)  :: rprdsh_inv(ix,km)
    real(r8), intent(out)  :: cmfsl_inv(ix,km)
    real(r8), intent(out)  :: cmflq_inv(ix,km)
    real(r8), intent(out)  :: qc_inv(ix,km)
    real(r8), intent(out)  :: icwmr_inv(ix,km)
    real(r8), intent(out)  :: rliq_inv(ix)
! +++ Yi-Chi : combined with ktop and kbot
    real(r8), intent(out)  :: cnt_inv(im)
    real(r8), intent(out)  :: cnb_inv(im)
! --- Yi-Chi
    real(r8), intent(out)  :: evapcv_inv(ix,km)

!cjshiu local variables for vertical transpose
    integer    :: kbot_tmp(im)
    integer    :: ktop_tmp(im)
    integer    :: kcnv_tmp(im)
    real(r8)   :: del_tmp(ix,km)
    real(r8)   :: prsl_tmp(ix,km)
    real(r8)   :: ps_tmp(im)
    real(r8)   :: phis_tmp(ix,km)
    real(r8)   :: ql_tmp(ix,km,2)
    real(r8)   :: q1_tmp(ix,km)
    real(r8)   :: t1_tmp(ix,km)
    real(r8)   :: u1_tmp(ix,km)
    real(r8)   :: v1_tmp(ix,km)
    real(r8)   :: rn_tmp(im)
    real(r8)   :: dot_tmp(ix,km)    

!cjshiu    real(r8)   :: qten_tmp(ix,km,2)
    real(r8)   :: cmfmc_tmp(ix,km)
    real(r8)   :: rprdsh_tmp(ix,km)
    real(r8)   :: cmfsl_tmp(ix,km)
    real(r8)   :: cmflq_tmp(ix,km)
    real(r8)   :: qc_tmp(ix,km)
    real(r8)   :: icwmr_tmp(ix,km)
    real(r8)   :: cnt_tmp(im)
    real(r8)   :: cnb_tmp(im)
    real(r8)   :: evapcv_tmp(ix,km)


    integer                :: k                        !  Vertical index for local fields [ no ] 
    integer                :: k_inv                    !  Vertical index for incoming fields [ no ]
    integer                :: m                    ! cjshiu check if we really need it?? 

!cjshiu try to move x 0.001 in convect_shallow to here
!  Yi-Chi : ix=>im
     ps_tmp(:im)=ps_inv(:im)*0.001                       !cjshiu Pa to hPa

    do k = 1, km
       k_inv         = km + 1 - k
       del_tmp(:ix,k)    = del_inv(:ix,k_inv)*0.001      !cjshiu Pa to hPa(or mb)
       prsl_tmp(:ix,k)   = prsl_inv(:ix,k_inv)*0.001     !cjshiu Pa to hPa (HP use cb)
       phis_tmp(:ix,k)   = zm_inv(:ix,k_inv)*gravit
       ! Yi-Chi
       ql_tmp(:ix,k,1)  = ql_inv(:ix,k_inv)
       ql_tmp(:ix,k,2)  = -1000. ! no ice, only liquid
       !ql_tmp(:ix,k,2)  = ql_inv(:ix,k_inv,2)
       ! ----
       q1_tmp(:ix,k)    = q1_inv(:ix,k_inv)
       t1_tmp(:ix,k)    = t1_inv(:ix,k_inv)
       u1_tmp(:ix,k)    = u1_inv(:ix,k_inv)
       v1_tmp(:ix,k)    = v1_inv(:ix,k_inv)
       dot_tmp(:ix,k)   = dot_inv(:ix,k_inv)*0.001       !cjshiu Pa to hPa
    enddo

    !--------------------
    ! Yi-Chi : add initialization for output variables  
    uten_inv(:ix,:km)         = 0.0
    vten_inv(:ix,:km)         = 0.0
    sten_inv(:ix,:km)         = 0.0
    do m = 1, pcnst
    qten_inv(:ix,:km,m)         = 0.0
    !enddo
    enddo
    ! --------------------


!+++ YiChi : This is the flag for whether deep convection is active.
    kcnv_tmp(:) = kcnv_inv(:)
!     write(iulog,*)'YiChi(hp_conv):kcnv_tmp:',kcnv_tmp
!---YiChi
    do m=1, im
      rn_tmp(m)=rn_inv(m)
    enddo


!cjs       cnt_tmp(:im) = 1._r8
!cjs       cnb_tmp(:im) = real(km, r8)
       cnt_tmp(:im) = real(km, r8)
       cnb_tmp(:im) = 1._r8

!cjs       write(iulog,*) 'in cnt_tmp=',cnt_tmp(:ix)
!cjs       write(iulog,*) 'in cnb_tmp=',cnb_tmp(:ix)
!cjshiu       write(iulog,*) 'before call shalcnv pcnst=',pcnst
    
!cjshiu check
    call gfuncphys

    call shalcnv(lchnk    , im , ix , km , pcnst , jcap_inv,            &
                 delt_inv , del_tmp , prsl_tmp , ps_tmp ,               &
                 phis_tmp , ql_tmp  ,                                   &
                 q1_tmp   , t1_tmp  , u1_tmp ,v1_tmp ,                  &
                 rcs_inv  , rn_tmp  , kbot_tmp , ktop_tmp ,             &
                 kcnv_tmp , slimsk_inv ,                                &
                 dot_tmp  , ncloud_inv , hpbl_inv ,                     &
                 heat_inv , evap_inv ,                                  &
                 cmfmc_tmp, rprdsh_tmp ,                                &
                 cmfsl_tmp, cmflq_tmp ,                                 &
                 qc_tmp   , icwmr_tmp , rliq_inv ,                      &
                 evapcv_tmp)

   ! Reverse cloud top/base interface indices          
!cjs    do m=1, im
!cjs       write(iulog,*) 'm= out cnt_tmp= ',m,cnt_tmp(m)
!cjs       write(iulog,*) 'm= out cnb_tmp= ',m,cnb_tmp(m)
 !cjs    enddo
! Yi-Chi :
      !jctop(i)     = pver-min(max(ktop(i),1),pver)+1
      !jcbot(i)     = pver-max(min(kbot(i),pver),1)+1
       cnt_inv(:im) = km + 1 - min(max(ktop_tmp(:im),1),km)
       cnb_inv(:im) = km + 1 - max(min(kbot_tmp(:im),km),1)

    do k = 1, km
       k_inv                       = km + 1 - k
!       ud_mf_inv(:ix,k_inv)        = ud_mf_tmp(:ix,k)
!       dt_mf_inv(:ix,k_inv)        = dt_mf_tmp(:ix,k)
!------Yi-Chi-------------
       sten_inv(:im,k_inv)         = (t1_tmp(:im,k)-t1_inv(:im,k_inv))*cp/delt_inv
       qten_inv(:im,k_inv,1)       = (q1_tmp(:im,k)-q1_inv(:im,k_inv))/delt_inv

       uten_inv(:im,k_inv)       = (u1_tmp(:im,k)-u1_inv(:im,k_inv))/delt_inv
       vten_inv(:im,k_inv)       = (v1_tmp(:im,k)-v1_inv(:im,k_inv))/delt_inv
       cmfmc_inv(:im,k_inv)      = cmfmc_tmp(:im,k) * 100./gravit ! from Pa/s ->  kg m-2 s
       rprdsh_inv(:im,k_inv)     = rprdsh_tmp(:im,k)
       cmfsl_inv(:im,k_inv)      = cmfsl_tmp(:im,k)
       cmflq_inv(:im,k_inv)      = cmflq_tmp(:im,k)
       qc_inv(:im,k_inv)         = qc_tmp(:im,k)
       icwmr_inv(:im,k_inv)      = icwmr_tmp(:im,k)
       evapcv_inv(:im,k_inv)     = evapcv_tmp(:im,k)
    enddo

!+++ YC debug for CMFDT++++   
! print out whenever shallow conv is activated.
!   do m = 1, im
!      if( maxval(cmfmc_inv(m,:km)) > 0._r8 ) then
!          write(iulog,*)'Yi-Chi:hp_conv466:sten_inv= ',sten_inv(m,:)
!          write(iulog,*)'Yi-Chi:hp_conv466:qten_inv= ',qten_inv(m,:,1)
!      end if
!   end do
!     write(iulog,*)'Yi-Chi:hp_conv452:t1_inv,t1_tmp= ',t1_tmp,t1_inv
!     write(iulog,*)'Yi-Chi:hp_conv452:cp,delt_inv= ',cp,delt_inv
!cjshiu specify smallest value to avoid floating invalid
!cjs    do k=1, km
!cjs     do m=1, im
!cjs      if (qc_inv(m,k) < 1.e-25_r8) then
!cjs         qc_inv(m,k) = 1.e-25_r8
!cjs      end if
!cjs    enddo
!cjs    enddo 
!cjshiu end specify    


!cjshiu note need to confirm the unit of precc in m or m/s if m/s need to do converting 
    do m =1,im
      rn_inv(m)=rn_tmp(m)/delt_inv
    enddo 

  end subroutine shalcnv_inv

  subroutine shalcnv(lchnk, im, ix, km, pcnst, jcap,                    &
                         delt, del, prsl, ps,                           &
                         phil, ql,                                      &
                         q1, t1, u1, v1,                                &
                         rcs, rn, kbot, ktop, kcnv, slimsk,             &
                         dot, ncloud, hpbl,                             &
!                         heat, evap, ud_mf, dt_mf,                      &
                         heat, evap, &
                         cmfmc_out , rprdsh_out , &
!                         sten_out , qten_out , cmfmc_out , rprdsh_out , &
                         cmfsl_out , cmflq_out ,                        &
                         qc_out , icwmr_out , rliq_out ,                &
!                         cnt_out , cnb_out , evapcv_out)
! +++ Yi-Chi
                         evapcv_out)
! --- Yi-Chi

   use shr_kind_mod, only: r8 => shr_kind_r8
!cjs      use cam_history,    only : outfld, addfld, phys_decomp
!cjs      use constituents,   only : qmin, cnst_get_type_byind, cnst_get_ind
!
!cjshiu direct assign r8
!cjs      use machine , only : kind_phys
!cjshiu use function under this subroutine
use module_gfs_funcphys, only: fpvs

    implicit none
!cjshiu added for additional input from CAM model
!cjshiu      integer , parameter :: r8 = selected_real_kind(12)    !  8 byte real
    integer , intent(in)  :: lchnk
    integer , intent(in)  :: pcnst
    integer   :: m
!cjshiu addend
    integer            im, ix,  km, jcap, ncloud,                     &
                         kbot(im), ktop(im), kcnv(im)
!    &,                  me
    real(r8)             delt
    real(r8)             ps(im),     del(ix,km),  prsl(ix,km),        &
                           ql(ix,km,2),q1(ix,km),   t1(ix,km),          &
                           u1(ix,km),  v1(ix,km),   rcs(im),            &
!cjshiu move rn to outputs     &                     rn(im),     slimsk(im), 
                           slimsk(im),                                  &
                           dot(ix,km), phil(ix,km), hpbl(im),           &
                           heat(im),   evap(im)

! hchuang code change mass flux output
      real(r8) ud_mf(im,km), dt_mf(im,km)
    real(r8), intent(out)   :: cmfmc_out(ix,km)       ! moist convection cloud mass flux [ kg/m2/s ]
    real(r8), intent(out)   :: rprdsh_out(ix,km)      ! dq/dt due to convective rainout  [ kg/m2/s ]
    real(r8), intent(out)   :: cmfsl_out(ix,km)       ! convective lw static energy flux [ J/kg * kg/m2/s = W/m2 ]
    real(r8), intent(out)   :: cmflq_out(ix,km)       ! convective total water flux [ kg/kg * kg/m2/s ]
    real(r8), intent(out)   :: rn(im)                 ! convective precipitation rate [m] (note: CAM uses m/s)
    real(r8), intent(out)   :: qc_out(ix,km)          ! dq/dt due to export of cloud water [ kg/kg/s ]
    real(r8), intent(out)   :: icwmr_out(ix,km)       ! cloud water mixing ratio (??)
    real(r8), intent(out)   :: rliq_out(im)           ! Vertical integral of qc_out [ m/s ]
! +++ Yi-Chi : combine with ktop and kbot
!    real(r8), intent(out)   :: cnt_out(im)            ! top level of convective activity
!    real(r8), intent(out)   :: cnb_out(im)            ! bottom level of convective activity
! --- Yi-Chi
    real(r8), intent(out)   :: evapcv_out(ix,km)       ! Tendency of evaporation of precipitation [ kg/kg/s ]
!cjshiu addend
!cjshiu need to check if there is unit issue

!
    integer              i,j,indx, jmn, k, kk, latd, lond, km1
    integer              kpbl(im)
!
    real(r8)             c0,      cpoel,   dellat,  delta,            &
                           desdt,   deta,    detad,   dg,               &
                           dh,      dhh,     dlnsig,  dp,               &
                           dq,      dqsdp,   dqsdt,   dt,               &
                           dt2,     dtmax,   dtmin,   dv1h,             &
                           dv1q,    dv2h,    dv2q,    dv1u,             &
                           dv1v,    dv2u,    dv2v,    dv3q,             &
                           dv3h,    dv3u,    dv3v,    clam,             &
                           dz,      dz1,     e1,                        &
                           el2orc,  elocp,   aafac,                     &
                           es,      etah,    h1,      dthk,             &
                           evef,    evfact,  evfactl, fact1,            &
                           fact2,   factor,  fjcap,                     &
                           g,       gamma,   pprime,  betaw,            &
                           qlk,     qrch,    qs,      c1,               &
                           rain,    rfact,   shear,   tem1,             &
                           tem2,    terr,    val,     val1,             &
                           val2,    w1,      w1l,     w1s,              &
                           w2,      w2l,     w2s,     w3,               &
                           w3l,     w3s,     w4,      w4l,              &
                           w4s,     tem,     ptem,    ptem1,            &
                           pgcon
!
    integer              kb(im), kbcon(im), kbcon1(im),               &
                           ktcon(im), ktcon1(im),                       &
                           kbm(im), kmax(im)
!
    real(r8)             aa1(im),                                     &
                           delhbar(im), delq(im),   delq2(im),          &
                           delqbar(im), delqev(im), deltbar(im),        &
                           deltv(im),   edt(im),                        &
                           wstar(im),   sflx(im),                       &
                           pdot(im),    po(im,km),                      &
                           qcond(im),   qevap(im),  hmax(im),           &
                           rntot(im),   vshear(im),                     &
                           xlamud(im),  xmb(im),    xmbmax(im),         &
                           delubar(im), delvbar(im)
    ! +++ Yi-Chi
      real(r8) delqlbar(im)
    ! ---
!c
    real(r8) cincr, cincrmax, cincrmin
    real(r8),parameter:: grav   =9.80665e+0     ! gravity           (m/s2)
    real(r8),parameter:: cp     =1.0046e+3      ! spec heat air @p    (J/kg/K)
    real(r8),parameter:: hvap   =2.5000e+6      ! lat heat H2O cond   (J/kg)
    real(r8),parameter:: rv     =4.6150e+2      ! gas constant H2O    (J/kg/K)
    real(r8),parameter:: rd     =2.8705e+2      ! gas constant air    (J/kg/K)
    real(r8),parameter:: t0c    =2.7315e+2      ! temp at 0C          (K)
    real(r8),parameter:: cvap   =1.8460e+3      ! spec heat H2O gas   (J/kg/K)
    real(r8),parameter:: cliq   =4.1855e+3      ! spec heat H2O liq   (J/kg/K)
    real(r8),parameter:: fv  =rv/rd-1.
    real(r8),parameter:: eps    =rd/rv
    real(r8),parameter:: epsm1  =rd/rv-1.

!cc
!cjshiu!c  physical parameters
    parameter(g=grav)
    parameter(cpoel=cp/hvap,elocp=hvap/cp,                            &
                el2orc=hvap*hvap/(rv*cp))
!cjshiu git some problem from input namelist file (restore to do debug) 
    parameter(terr=0.,c0=.002,c1=5.e-4) ! default
!    parameter(terr=0.,c0=.002,c1=1.e-2)

    parameter(delta=fv)
    parameter(fact1=(cvap-cliq)/rv,fact2=hvap/rv-fact1*t0c)
    parameter(cincrmax=180.,cincrmin=120.,dthk=25.)
    parameter(h1=0.33333333)

! Yi-Chi +++
    real(r8)    rrho(im) ! air density used for calculating sflx
    real(r8)    icwmr_tmp1(im,km) !
! Yi-Chi ---

!c  local variables and arrays
    real(r8)             pfld(im,km),    to(im,km),     qo(im,km),    &
                           uo(im,km),      vo(im,km),     qeso(im,km)
!c  cloud water
!     real(kind=kind_phys) qlko_ktcon(im), dellal(im,km), tvo(im,km),
    real(r8)             qlko_ktcon(im), dellal(im,km),               &
                           dbyo(im,km),    zo(im,km),     xlamue(im,km),&
                           heo(im,km),     heso(im,km),                 &
                           dellah(im,km),  dellaq(im,km),               &
                           dellau(im,km),  dellav(im,km), hcko(im,km),  &
                           ucko(im,km),    vcko(im,km),   qcko(im,km),  &
                           eta(im,km),     zi(im,km),     pwo(im,km),   &
! +++ Yi-Chi : June, 2013
                           qrcko(im,km), tx1(im)
!                           tx1(im)
! ------------
!
    logical totflg, cnvflg(im), flg(im)
!
    real(r8) tf, tcr, tcrf
    parameter (tf=233.16, tcr=263.16, tcrf=1.0/(tcr-tf))

!cjshiu add initial values for the following parameters
      do k = 1, km
        do i = 1, im
         ud_mf(i,k) = 0.0_r8
         dt_mf(i,k) = 0.0_r8
         cmfmc_out(i,k) = 0._r8      ! moist convection cloud mass flux [ kg/m2/s ]
         rprdsh_out(i,k) = 0._r8     ! dq/dt due to convective rainout  [ kg/m2/s ]
         cmfsl_out(i,k) = 0._r8      ! convective lw static energy flux [ J/kg * kg/m2/s = W/m2 ]
         cmflq_out(i,k) = 0._r8      ! convective total water flux [ kg/kg * kg/m2/s ]
!cj    rn(:im) = 0.0_r8                ! convective precipitation rate [m] (note: CAM uses m/s)
         qc_out(i,k) = 0.0_r8         ! dq/dt due to export of cloud water [ kg/kg/s ]
         icwmr_tmp1(i,k) = 0.0_r8
         icwmr_out(i,k) = 0.0_r8      ! cloud water mixing ratio (??)
         evapcv_out(i,k) = 0.0_r8      ! Tendency of evaporation of precipitation [ kg/kg/s ]
        enddo
      enddo
      do i = 1, im
         rn(i) = 0.0_r8                ! convective precipitation rate [m] (note: CAM uses m/s)
         rliq_out(i) = 0.
      enddo
     
      !check
      !print *, 'run after initialization in shalcnv'

!c-----------------------------------------------------------------------
!
      km1 = km - 1

!cjshiu  check
!cjs      print *, 'run before compute surface buoyancy flux in shalcnv'
!cjs      print *, 'km1= km= ',km1,km
!c
!c  compute surface buoyancy flux
!c  
      do i=1,im
        ! Yi-Chi +++
        ! correct surface buoyancy flux with air density
        !print*, '(Yi-Chi (hp_conv): rd, t1, prsl)', rd, t1(i,1), prsl(i,1)
        !print*, '(Yi-Chi (hp_conv): prsl)', prsl
        rrho(i)    = rd*t1(i,1)/(prsl(i,1)*1000.) ! inverse of density (m3/kg); rair->rd
        !print*, '(Yi-Chi (hp_conv): sflx,cflx, rrho', heat(i), evap(i), rrho(i)
        heat(i)    = heat(i)*rrho(i)/cp  ! heat : shflx->heat flux
        evap(i)    = evap(i)*rrho(i)
        sflx(i)    = heat(i) + fv*(t1(i,1)*(100000./(prsl(i,1)*1000.))**(rd/cp))*evap(i)
        !print*, '(Yi-Chi (hp_conv): surf t1,theta1',t1(i,1), t1(i,1)*(100000./(prsl(i,1)*1000.))**(rd/cp)
        !print*, '(Yi-Chi (hp_conv): cp, heat,evap, sflx,fv',cp, heat(i), evap(i), sflx(i),fv
        ! sflx(i) = heat(i)+fv*t1(i,1)*evap(i)
        ! Yi-Chi ---
!cjshiu  check
!cjs      print *, 'run after compute buoyancy flux i= sflx(i)= fv= ',i,sflx(i),fv
!cjs      print *, 'run after compute buoyancy flux i= heat(i)= t1(i,1)= evap(i)= ',i,heat(i),t1(i,1),evap(i)
      enddo
!c
!c  initialize arrays
!c
      do i=1,im
        cnvflg(i) = .true.
!cjshiu  check
!cjs      print *, 'run in initialize arrays i= kcnv(i)=',i,kcnv
        if(kcnv(i).eq.1) cnvflg(i) = .false.
        ! write(iulog,*) 'YiChi(hp_conv): kcnv',kcnv(i)
        if(sflx(i).le.0.) cnvflg(i) = .false.
        if(cnvflg(i)) then
          kbot(i)=km+1
          ktop(i)=0
        endif
        rn(i)=0.
        kbcon(i)=km
        ktcon(i)=1
        kb(i)=km
        pdot(i) = 0.
        qlko_ktcon(i) = 0.
        edt(i)  = 0.
        aa1(i)  = 0.
        vshear(i) = 0.
      enddo
! hchuang code change
      do k = 1, km
        do i = 1, im
          ud_mf(i,k) = 0.
          dt_mf(i,k) = 0.
        enddo
      enddo
!cjshiu  check
!      print *, 'run before totflg1'
!!
      totflg = .true.
      do i=1,im
        totflg = totflg .and. (.not. cnvflg(i))
      enddo
!cjshiu  check
!      print *, 'run before totflg2'

      if(totflg) return
!cjshiu  check
!      print *, 'run after if(totflg) return'
!!
!c
      dt2   = delt
      val   =         1200.
      dtmin = max(dt2, val )
      val   =         3600.
      dtmax = max(dt2, val )
!c  model tunable parameters are all here
      clam    = .3
      aafac   = .1
      betaw   = .03
!c     evef    = 0.07
      evfact  = 0.3
      evfactl = 0.3
!
!     pgcon   = 0.7     ! Gregory et al. (1997, QJRMS)
      pgcon   = 0.55    ! Zhang & Wu (2003,JAS)
      fjcap   = (float(jcap) / 126.) ** 2
      val     =           1.
      fjcap   = max(fjcap,val)
      w1l     = -8.e-3
      w2l     = -4.e-2
      w3l     = -5.e-3
      w4l     = -5.e-4
      w1s     = -2.e-4
      w2s     = -2.e-3
      w3s     = -1.e-3
      w4s     = -2.e-5
!c
!c  define top layer for search of the downdraft originating layer
!c  and the maximum thetae for updraft
!c
      do i=1,im
        kbm(i)   = km
        kmax(i)  = km
        tx1(i)   = 1.0 / ps(i)
      enddo
!     
      do k = 1, km
        do i=1,im
          if (prsl(i,k)*tx1(i) .gt. 0.70) kbm(i)   = k + 1
          if (prsl(i,k)*tx1(i) .gt. 0.60) kmax(i)  = k + 1
        enddo
      enddo
      do i=1,im
        kbm(i)   = min(kbm(i),kmax(i))
      enddo

!c
!c  hydrostatic height assume zero terr and compute
!c  updraft entrainment rate as an inverse function of height
!c
      do k = 1, km
        do i=1,im
          zo(i,k) = phil(i,k) / g
        enddo
      enddo
      do k = 1, km1
        do i=1,im
          zi(i,k) = 0.5*(zo(i,k)+zo(i,k+1))
          xlamue(i,k) = clam / zi(i,k)
        enddo
      enddo
      do i=1,im
        xlamue(i,km) = xlamue(i,km1)
      enddo
!c
!c  pbl height
!c
      do i=1,im
        flg(i) = cnvflg(i)
        kpbl(i)= 1
      enddo
      do k = 2, km1
        do i=1,im
          if (flg(i).and.zo(i,k).le.hpbl(i)) then
            kpbl(i) = k
          else
            flg(i) = .false.
          endif
        enddo
      enddo
      do i=1,im
        kpbl(i)= min(kpbl(i),kbm(i))
      enddo

!c
!c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!c   convert surface pressure to mb from cb
!c
      do k = 1, km
        do i = 1, im
          if (cnvflg(i) .and. k .le. kmax(i)) then
            pfld(i,k) = prsl(i,k) * 10.0
            eta(i,k)  = 1.
            hcko(i,k) = 0.
            qcko(i,k) = 0.
            ! +++ Yi-Chi : June 22, 2013
            qrcko(i,k)= 0.
            ! ---
            ucko(i,k) = 0.
            vcko(i,k) = 0.
            dbyo(i,k) = 0.
            pwo(i,k)  = 0.
            dellal(i,k) = 0.
            to(i,k)   = t1(i,k)
            qo(i,k)   = q1(i,k)
            uo(i,k)   = u1(i,k) * rcs(i)
            vo(i,k)   = v1(i,k) * rcs(i)
          endif
        enddo
      enddo
!c
!c  column variables
!c  p is pressure of the layer (mb)
!c  t is temperature at t-dt (k)..tn
!c  q is mixing ratio at t-dt (kg/kg)..qn
!c  to is temperature at t+dt (k)... this is after advection and turbulan
!c  qo is mixing ratio at t+dt (kg/kg)..q1
!c
      do k = 1, km
        do i=1,im
          if (cnvflg(i) .and. k .le. kmax(i)) then
            qeso(i,k) = 0.01 * fpvs(to(i,k))      ! fpvs is in pa
            qeso(i,k) = eps * qeso(i,k) / (pfld(i,k) + epsm1*qeso(i,k))
            val1      =             1.e-8
            qeso(i,k) = max(qeso(i,k), val1)
            val2      =           1.e-10
            qo(i,k)   = max(qo(i,k), val2 )
!           qo(i,k)   = min(qo(i,k),qeso(i,k))
!           tvo(i,k)  = to(i,k) + delta * to(i,k) * qo(i,k)
          endif
        enddo
      enddo

!cjshiu  check
!      print *, 'run before compute moist static energy in shalcnv'
!      print *, 'val1= val2= ',val1,val2
!c
!c  compute moist static energy
!c
      do k = 1, km
        do i=1,im
          if (cnvflg(i) .and. k .le. kmax(i)) then
!           tem       = g * zo(i,k) + cp * to(i,k)
            tem       = phil(i,k) + cp * to(i,k)
            heo(i,k)  = tem  + hvap * qo(i,k)
            heso(i,k) = tem  + hvap * qeso(i,k)
!c           heo(i,k)  = min(heo(i,k),heso(i,k))
          endif
        enddo
      enddo
!c
!c  determine level with largest moist static energy within pbl
!c  this is the level where updraft starts
!c
      do i=1,im
         if (cnvflg(i)) then
            hmax(i) = heo(i,1)
            kb(i) = 1
         endif
      enddo
      do k = 2, km
        do i=1,im
          if (cnvflg(i).and.k.le.kpbl(i)) then
            if(heo(i,k).gt.hmax(i)) then
              kb(i)   = k
              hmax(i) = heo(i,k)
            endif
          endif
        enddo
      enddo
!c
      do k = 1, km1
        do i=1,im
          if (cnvflg(i) .and. k .le. kmax(i)-1) then
            dz      = .5 * (zo(i,k+1) - zo(i,k))
            dp      = .5 * (pfld(i,k+1) - pfld(i,k))
!cjshiu check floating invalid
!cjs            write(iulog,*)'i= k= to(i,k+1)= dz= dp= ',i,k,to(i,k+1),dz,dp 
            es      = 0.01 * fpvs(to(i,k+1))      ! fpvs is in pa
            pprime  = pfld(i,k+1) + epsm1 * es
            qs      = eps * es / pprime
            dqsdp   = - qs / pprime
            desdt   = es * (fact1 / to(i,k+1) + fact2 / (to(i,k+1)**2))
!cjshiu check floating invalid
!cjs            write(iulog,*)'i= k= qs= pfld(i,k+1)= desdt= es= pprime= ',i,k,qs,pfld(i,k+1),desdt,es,pprime 
            dqsdt   = qs * pfld(i,k+1) * desdt / (es * pprime)
            gamma   = el2orc * qeso(i,k+1) / (to(i,k+1)**2)
            dt      = (g * dz + hvap * dqsdp * dp) / (cp * (1. + gamma))
            dq      = dqsdt * dt + dqsdp * dp
            to(i,k) = to(i,k+1) + dt
            qo(i,k) = qo(i,k+1) + dq
            po(i,k) = .5 * (pfld(i,k) + pfld(i,k+1))
          endif
        enddo
      enddo
!
      do k = 1, km1
        do i=1,im
          if (cnvflg(i) .and. k .le. kmax(i)-1) then
            qeso(i,k) = 0.01 * fpvs(to(i,k))      ! fpvs is in pa
            qeso(i,k) = eps * qeso(i,k) / (po(i,k) + epsm1*qeso(i,k))
            val1      =             1.e-8
            qeso(i,k) = max(qeso(i,k), val1)
            val2      =           1.e-10
            qo(i,k)   = max(qo(i,k), val2 )
!           qo(i,k)   = min(qo(i,k),qeso(i,k))
            heo(i,k)  = .5 * g * (zo(i,k) + zo(i,k+1)) +                &
                        cp * to(i,k) + hvap * qo(i,k)
            heso(i,k) = .5 * g * (zo(i,k) + zo(i,k+1)) +                &
                        cp * to(i,k) + hvap * qeso(i,k)
            uo(i,k)   = .5 * (uo(i,k) + uo(i,k+1))
            vo(i,k)   = .5 * (vo(i,k) + vo(i,k+1))
          endif
        enddo
      enddo

!c
!c  look for the level of free convection as cloud base
!c
      do i=1,im
        flg(i)   = cnvflg(i)
        if(flg(i)) kbcon(i) = kmax(i)
      enddo
      do k = 2, km1
        do i=1,im
          if (flg(i).and.k.lt.kbm(i)) then
            if(k.gt.kb(i).and.heo(i,kb(i)).gt.heso(i,k)) then
              kbcon(i) = k
              flg(i)   = .false.
            endif
          endif
        enddo
      enddo
!c
      do i=1,im
        if(cnvflg(i)) then
          if(kbcon(i).eq.kmax(i)) cnvflg(i) = .false.
        endif
      enddo
!!
      totflg = .true.
      do i=1,im
        totflg = totflg .and. (.not. cnvflg(i))
      enddo
      if(totflg) return
!!
!c
!c  determine critical convective inhibition
!c  as a function of vertical velocity at cloud base.
!c
      do i=1,im
        if(cnvflg(i)) then
          pdot(i)  = 10.* dot(i,kbcon(i))
        endif
      enddo
      do i=1,im
        if(cnvflg(i)) then
          if(slimsk(i).eq.1.) then
            w1 = w1l
            w2 = w2l
            w3 = w3l
            w4 = w4l
          else
            w1 = w1s
            w2 = w2s
            w3 = w3s
            w4 = w4s
          endif
          if(pdot(i).le.w4) then
            ptem = (pdot(i) - w4) / (w3 - w4)
          elseif(pdot(i).ge.-w4) then
            ptem = - (pdot(i) + w4) / (w4 - w3)
          else
            ptem = 0.
          endif
          val1    =             -1.
          ptem = max(ptem,val1)
          val2    =             1.
          ptem = min(ptem,val2)
          ptem = 1. - ptem
          ptem1= .5*(cincrmax-cincrmin)
          cincr = cincrmax - ptem * ptem1
          tem1 = pfld(i,kb(i)) - pfld(i,kbcon(i))
          if(tem1.gt.cincr) then
             cnvflg(i) = .false.
          endif
        endif
      enddo
!!
      totflg = .true.
      do i=1,im
        totflg = totflg .and. (.not. cnvflg(i))
      enddo
      if(totflg) return

!!
!c
!c  assume the detrainment rate for the updrafts to be same as 
!c  the entrainment rate at cloud base
!c
      do i = 1, im
        if(cnvflg(i)) then
          xlamud(i) = xlamue(i,kbcon(i))
        endif
      enddo

!cjshiu check
!      write(iulog,*)'run to determine updraft mass flux for the subcloud layers'
!c
!c  determine updraft mass flux for the subcloud layers
!c
      do k = km1, 1, -1
        do i = 1, im
          if (cnvflg(i)) then
            if(k.lt.kbcon(i).and.k.ge.kb(i)) then
              dz       = zi(i,k+1) - zi(i,k)
              ptem     = 0.5*(xlamue(i,k)+xlamue(i,k+1))-xlamud(i)
              eta(i,k) = eta(i,k+1) / (1. + ptem * dz)
            endif
          endif
        enddo
      enddo
!c
!c  compute mass flux above cloud base
!c
      do k = 2, km1
        do i = 1, im
         if(cnvflg(i))then
           if(k.gt.kbcon(i).and.k.lt.kmax(i)) then
              dz       = zi(i,k) - zi(i,k-1)
              ptem     = 0.5*(xlamue(i,k)+xlamue(i,k-1))-xlamud(i)
              eta(i,k) = eta(i,k-1) * (1 + ptem * dz)
           endif
         endif
        enddo
      enddo
!c
!c  compute updraft cloud property
!c
      do i = 1, im
        if(cnvflg(i)) then
          indx         = kb(i)
          hcko(i,indx) = heo(i,indx)
!cjshiu output this
          cmfsl_out(i,indx) = hcko(i,indx)
!cjshiu add end
          ucko(i,indx) = uo(i,indx)
          vcko(i,indx) = vo(i,indx)
        endif
      enddo
!c
      do k = 2, km1
        do i = 1, im
          if (cnvflg(i)) then
            if(k.gt.kb(i).and.k.lt.kmax(i)) then
              dz   = zi(i,k) - zi(i,k-1)
              tem  = 0.5 * (xlamue(i,k)+xlamue(i,k-1)) * dz
              tem1 = 0.5 * xlamud(i) * dz
              factor = 1. + tem - tem1
              ptem = 0.5 * tem + pgcon
              ptem1= 0.5 * tem - pgcon
              hcko(i,k) = ((1.-tem1)*hcko(i,k-1)+tem*0.5*               &
                           (heo(i,k)+heo(i,k-1)))/factor
!cjshiu output this
              cmfsl_out(i,k) = hcko(i,k)
!cjshiu add end
              ucko(i,k) = ((1.-tem1)*ucko(i,k-1)+ptem*uo(i,k)           &
                           +ptem1*uo(i,k-1))/factor
              vcko(i,k) = ((1.-tem1)*vcko(i,k-1)+ptem*vo(i,k)           &
                           +ptem1*vo(i,k-1))/factor
              dbyo(i,k) = hcko(i,k) - heso(i,k)
            endif
          endif
        enddo
      enddo
!c
!c   taking account into convection inhibition due to existence of
!c    dry layers below cloud base
!c
      do i=1,im
        flg(i) = cnvflg(i)
        kbcon1(i) = kmax(i)
      enddo
      do k = 2, km1
      do i=1,im
        if (flg(i).and.k.lt.kbm(i)) then
          if(k.ge.kbcon(i).and.dbyo(i,k).gt.0.) then
            kbcon1(i) = k
            flg(i)    = .false.
          endif
        endif
      enddo
      enddo
      do i=1,im
        if(cnvflg(i)) then
          if(kbcon1(i).eq.kmax(i)) cnvflg(i) = .false.
        endif
      enddo
      do i=1,im
        if(cnvflg(i)) then
          tem = pfld(i,kbcon(i)) - pfld(i,kbcon1(i))
          if(tem.gt.dthk) then
             cnvflg(i) = .false.
          endif
        endif
      enddo
!!
      totflg = .true.
      do i = 1, im
        totflg = totflg .and. (.not. cnvflg(i))
      enddo
      if(totflg) return
!!
!c
!c  determine first guess cloud top as the level of zero buoyancy
!c    limited to the level of sigma=0.7
!c
      do i = 1, im
        flg(i) = cnvflg(i)
        if(flg(i)) ktcon(i) = kbm(i)
      enddo
      do k = 2, km1
      do i=1,im
        if (flg(i).and.k .lt. kbm(i)) then
          if(k.gt.kbcon1(i).and.dbyo(i,k).lt.0.) then
             ktcon(i) = k
             flg(i)   = .false.
          endif
        endif
      enddo
      enddo
!c
!c  turn off shallow convection if cloud top is less than pbl top
!c
!     do i=1,im
!       if(cnvflg(i)) then
!         kk = kpbl(i)+1
!         if(ktcon(i).le.kk) cnvflg(i) = .false.
!       endif
!     enddo
!!
!     totflg = .true.
!     do i = 1, im
!       totflg = totflg .and. (.not. cnvflg(i))
!     enddo
!     if(totflg) return
!!
!c
!c  specify upper limit of mass flux at cloud base
!c
      do i = 1, im
        if(cnvflg(i)) then
!         xmbmax(i) = .1
!
          k = kbcon(i)
          dp = 1000. * del(i,k)
          xmbmax(i) = dp / (g * dt2)
!
!         tem = dp / (g * dt2)
!         xmbmax(i) = min(tem, xmbmax(i))
        endif
      enddo
!c
!c  compute cloud moisture property and precipitation
!c
      do i = 1, im
        if (cnvflg(i)) then
          aa1(i) = 0.
          qcko(i,kb(i)) = qo(i,kb(i))
          ! +++ Yi-Chi : June 22, 2013
          qrcko(i,kb(i)) = qo(i,kb(i))
          ! ---
        endif
      enddo
      do k = 2, km1
        do i = 1, im
          if (cnvflg(i)) then
            if(k.gt.kb(i).and.k.lt.ktcon(i)) then
              dz    = zi(i,k) - zi(i,k-1)
              gamma = el2orc * qeso(i,k) / (to(i,k)**2)
              qrch = qeso(i,k)                                          &
                  + gamma * dbyo(i,k) / (hvap * (1. + gamma))
!cj
              tem  = 0.5 * (xlamue(i,k)+xlamue(i,k-1)) * dz
              tem1 = 0.5 * xlamud(i) * dz
              factor = 1. + tem - tem1
              qcko(i,k) = ((1.-tem1)*qcko(i,k-1)+tem*0.5*               &
                          (qo(i,k)+qo(i,k-1)))/factor
          ! +++ Yi-Chi : June 22, 2013
          qrcko(i,k) = qcko(i,k)
          ! --
              dq = eta(i,k) * (qcko(i,k) - qrch)
!c
!             rhbar(i) = rhbar(i) + qo(i,k) / qeso(i,k)
!c
!c  below lfc check if there is excess moisture to release latent heat
!c
              if(k.ge.kbcon(i).and.dq.gt.0.) then
                etah = .5 * (eta(i,k) + eta(i,k-1))
                if(ncloud.gt.0.) then
                  dp = 1000. * del(i,k)
                  qlk = dq / (eta(i,k) + etah * (c0 + c1) * dz)
                  dellal(i,k) = etah * c1 * dz * qlk * g / dp
!!cjshiu output this
!                  qc_out(i,k)=dellal(i,k)
!!cjshiu add end
                else
                  qlk = dq / (eta(i,k) + etah * c0 * dz)
                endif
                aa1(i) = aa1(i) - dz * g * qlk
                qcko(i,k)= qlk + qrch
                pwo(i,k) = etah * c0 * dz * qlk
!cjshiu output this
!                cmfmc_out(i,k) = qcko(i,k)
!cjs check 
!cjs                write(iulog,*)'i= k= qcko(i,k)= pwo(i,k)=',i,k,qcko(i,k),pwo(i,k)
!cjs                print *, 'i= k= qcko(i,k)= pwo(i,k)=',i,k,qcko(i,k),pwo(i,k)
! Yi-Chi----------
                 icwmr_tmp1(i,k) = eta(i,k) * qlk
!                icwmr_out(i,k) = etah * qlk * g / (1000. * del(i,k))
!                icwmr_out(i,k) = cmfmc_out(i,k)*g/(1000. * del(i,k))*dt2
                 !icwmr_out(i,k) = qcko(i,k)*g/(1000. * del(i,k))*dt2
!                rprdsh_out(i,k) = pwo(i,k) * g / (1000. * del(i,k))
!cjshiu add end
              endif
            endif
          endif
        enddo
      enddo

!cjshiu check
!     write(iulog,*)'run to before calculate cloud work function'

!c
!c  calculate cloud work function
!c
      do k = 2, km1
        do i = 1, im
          if (cnvflg(i)) then
            if(k.ge.kbcon(i).and.k.lt.ktcon(i)) then
              dz1 = zo(i,k+1) - zo(i,k)
              gamma = el2orc * qeso(i,k) / (to(i,k)**2)
              rfact =  1. + delta * cp * gamma                          &
                      * to(i,k) / hvap
              aa1(i) = aa1(i) +                                         &
                       dz1 * (g / (cp * to(i,k)))                       &
                       * dbyo(i,k) / (1. + gamma)                       &
                      * rfact
              val = 0.
              aa1(i)=aa1(i)+                                            &
                       dz1 * g * delta *                                &
                      max(val,(qeso(i,k) - qo(i,k)))
            endif
          endif
        enddo
      enddo
      do i = 1, im
        if(cnvflg(i).and.aa1(i).le.0.) cnvflg(i) = .false.
      enddo
!!
      totflg = .true.
      do i=1,im
        totflg = totflg .and. (.not. cnvflg(i))
      enddo
      if(totflg) return

!!
!c
!c  estimate the onvective overshooting as the level
!c    where the [aafac * cloud work function] becomes zero,
!c    which is the final cloud top
!c    limited to the level of sigma=0.7
!c
      do i = 1, im
        if (cnvflg(i)) then
          aa1(i) = aafac * aa1(i)
        endif
      enddo
!c
      do i = 1, im
        flg(i) = cnvflg(i)
        ktcon1(i) = kbm(i)
      enddo
      do k = 2, km1
        do i = 1, im
          if (flg(i)) then
            if(k.ge.ktcon(i).and.k.lt.kbm(i)) then
              dz1 = zo(i,k+1) - zo(i,k)
              gamma = el2orc * qeso(i,k) / (to(i,k)**2)
              rfact =  1. + delta * cp * gamma                          &
                      * to(i,k) / hvap
              aa1(i) = aa1(i) +                                         &
                       dz1 * (g / (cp * to(i,k)))                       &
                       * dbyo(i,k) / (1. + gamma)                       &
                       * rfact
              if(aa1(i).lt.0.) then
                ktcon1(i) = k
                flg(i) = .false.
              endif
            endif
          endif
        enddo
      enddo
!c
!c  compute cloud moisture property, detraining cloud water
!c    and precipitation in overshooting layers
!c
      do k = 2, km1
        do i = 1, im
          if (cnvflg(i)) then
            if(k.ge.ktcon(i).and.k.lt.ktcon1(i)) then
              dz    = zi(i,k) - zi(i,k-1)
              gamma = el2orc * qeso(i,k) / (to(i,k)**2)                 
              qrch = qeso(i,k)                                          &
                  + gamma * dbyo(i,k) / (hvap * (1. + gamma))
!cj
              tem  = 0.5 * (xlamue(i,k)+xlamue(i,k-1)) * dz
              tem1 = 0.5 * xlamud(i) * dz
              factor = 1. + tem - tem1
              qcko(i,k) = ((1.-tem1)*qcko(i,k-1)+tem*0.5*               &
                          (qo(i,k)+qo(i,k-1)))/factor
! Yi-Chi
              ! +++ Yi-Chi: June 2013
              qrcko(i,k) = qcko(i,k)
              ! ----
!!cjshiu output this
!              cmfmc_out(i,k) = qcko(i,k)
!!cjshiu add end
!cj
              dq = eta(i,k) * (qcko(i,k) - qrch)
!c
!c  check if there is excess moisture to release latent heat
!c
              if(dq.gt.0.) then
                etah = .5 * (eta(i,k) + eta(i,k-1))
                if(ncloud.gt.0.) then
                  dp = 1000. * del(i,k)
                  qlk = dq / (eta(i,k) + etah * (c0 + c1) * dz)
                  dellal(i,k) = etah * c1 * dz * qlk * g / dp
!cjshiu output this
!                  qc_out(i,k)=dellal(i,k)
!cjshiu add end
                else
                  qlk = dq / (eta(i,k) + etah * c0 * dz)
                endif
                qcko(i,k) = qlk + qrch
                pwo(i,k) = etah * c0 * dz * qlk
! Yi-Chi
!!cjshiu output this
!                cmfmc_out(i,k) = qcko(i,k)
                 icwmr_tmp1(i,k) = eta(i,k) * qlk
!                 icwmr_out(i,k) = etah * qlk * g / (1000. * del(i,k))
!                icwmr_out(i,k) = cmfmc_out(i,k)*g/(1000. * del(i,k))*dt2
! Yi-Chi
!                rprdsh_out(i,k) = pwo(i,k) * g / (1000. * del(i,k))
!cjshiu add end
              endif
            endif
          endif
        enddo
      enddo
!c
!c exchange ktcon with ktcon1
!c
      do i = 1, im
        if(cnvflg(i)) then
          kk = ktcon(i)
          ktcon(i) = ktcon1(i)
          ktcon1(i) = kk
        endif
      enddo
!c
!c  this section is ready for cloud water
!c
      if(ncloud.gt.0) then
!c
!c  compute liquid and vapor separation at cloud top
!c
      do i = 1, im
        if(cnvflg(i)) then
          k = ktcon(i) - 1
          gamma = el2orc * qeso(i,k) / (to(i,k)**2)
          qrch = qeso(i,k)                                              &
               + gamma * dbyo(i,k) / (hvap * (1. + gamma))
          dq = qcko(i,k) - qrch
!c
!c  check if there is excess moisture to release latent heat
!c
          if(dq.gt.0.) then
            qlko_ktcon(i) = dq
            qcko(i,k) = qrch
          endif
        endif
      enddo
      endif
!c
!c--- compute precipitation efficiency in terms of windshear
!c
      do i = 1, im
        if(cnvflg(i)) then
          vshear(i) = 0.
        endif
      enddo
      do k = 2, km
        do i = 1, im
          if (cnvflg(i)) then
            if(k.gt.kb(i).and.k.le.ktcon(i)) then
              shear= sqrt((uo(i,k)-uo(i,k-1)) ** 2                     &
                       + (vo(i,k)-vo(i,k-1)) ** 2)
              vshear(i) = vshear(i) + shear
            endif
          endif
        enddo
      enddo
      do i = 1, im
        if(cnvflg(i)) then
          vshear(i) = 1.e3 * vshear(i) / (zi(i,ktcon(i))-zi(i,kb(i)))  
          e1=1.591-.639*vshear(i)                                      &
            +.0953*(vshear(i)**2)-.00496*(vshear(i)**3)
          edt(i)=1.-e1
          val =         .9
          edt(i) = min(edt(i),val)
          val =         .0
          edt(i) = max(edt(i),val)
        endif
      enddo
!c
!c--- what would the change be, that a cloud with unit mass
!c--- will do to the environment?
!c
      do k = 1, km
        do i = 1, im
          if(cnvflg(i) .and. k .le. kmax(i)) then
            dellah(i,k) = 0.
            dellaq(i,k) = 0.
            dellau(i,k) = 0.
            dellav(i,k) = 0.
          endif
        enddo
      enddo
!c
!c--- changed due to subsidence and entrainment
!c
      do k = 2, km1
        do i = 1, im
          if (cnvflg(i)) then
            if(k.gt.kb(i).and.k.lt.ktcon(i)) then
              dp = 1000. * del(i,k)
              dz = zi(i,k) - zi(i,k-1)
!c
              dv1h = heo(i,k)
              dv2h = .5 * (heo(i,k) + heo(i,k-1))
              dv3h = heo(i,k-1)
              dv1q = qo(i,k)
              dv2q = .5 * (qo(i,k) + qo(i,k-1))
              dv3q = qo(i,k-1)
              dv1u = uo(i,k)
              dv2u = .5 * (uo(i,k) + uo(i,k-1))
              dv3u = uo(i,k-1)
              dv1v = vo(i,k)
              dv2v = .5 * (vo(i,k) + vo(i,k-1))
              dv3v = vo(i,k-1)
!c
              tem  = 0.5 * (xlamue(i,k)+xlamue(i,k-1))
              tem1 = xlamud(i)
!cj
              dellah(i,k) = dellah(i,k) +                               &
           ( eta(i,k)*dv1h - eta(i,k-1)*dv3h                            &
          -  tem*eta(i,k-1)*dv2h*dz                                     &
          +  tem1*eta(i,k-1)*.5*(hcko(i,k)+hcko(i,k-1))*dz              &
               ) *g/dp
!cj
              dellaq(i,k) = dellaq(i,k) +                               &
           ( eta(i,k)*dv1q - eta(i,k-1)*dv3q                            &
          -  tem*eta(i,k-1)*dv2q*dz                                     &
!+++ Yi-Chi : June 2013
           +  tem1*eta(i,k-1)*.5*(qrcko(i,k)+qcko(i,k-1))*dz            &
!          +  tem1*eta(i,k-1)*.5*(qcko(i,k)+qcko(i,k-1))*dz              &
! -----
               ) *g/dp
! Yi-Chi : This two quantities are not correct!
!cjshiu output this
!              sten_out(i,k) = dellah(i,k)
!              cmflq_out(i,k) = dellaq(i,k) * dp / g
!cjshiu add end
!cj
              dellau(i,k) = dellau(i,k) +                               &
           ( eta(i,k)*dv1u - eta(i,k-1)*dv3u                            &
          -  tem*eta(i,k-1)*dv2u*dz                                     &
          +  tem1*eta(i,k-1)*.5*(ucko(i,k)+ucko(i,k-1))*dz              &
          -  pgcon*eta(i,k-1)*(dv1u-dv3u)                               &
               ) *g/dp
!cj
!cj
              dellav(i,k) = dellav(i,k) +                               &
           ( eta(i,k)*dv1v - eta(i,k-1)*dv3v                            &
          -  tem*eta(i,k-1)*dv2v*dz                                     &
          +  tem1*eta(i,k-1)*.5*(vcko(i,k)+vcko(i,k-1))*dz              &
          -  pgcon*eta(i,k-1)*(dv1v-dv3v)                               &
               ) *g/dp
!cj
            endif
          endif
        enddo
      enddo
!c
!c------- cloud top
!c
      do i = 1, im
        if(cnvflg(i)) then
          indx = ktcon(i)
          dp = 1000. * del(i,indx)
          dv1h = heo(i,indx-1)
          dellah(i,indx) = eta(i,indx-1) *                              &
                          (hcko(i,indx-1) - dv1h) * g / dp
          dv1q = qo(i,indx-1)
          dellaq(i,indx) = eta(i,indx-1) *                              &
                          (qcko(i,indx-1) - dv1q) * g / dp
! Yi-Chi : These two are not stend_out
!cjshiu output this
!          sten_out(i,indx) = dellah(i,indx)
!          cmflq_out(i,indx) = dellaq(i,indx) * dp / g
!cjshiu add end
          dv1u = uo(i,indx-1)
          dellau(i,indx) = eta(i,indx-1) *                              &
                          (ucko(i,indx-1) - dv1u) * g / dp
          dv1v = vo(i,indx-1)
          dellav(i,indx) = eta(i,indx-1) *                              &
                          (vcko(i,indx-1) - dv1v) * g / dp
!c
!c  cloud water
!c
          dellal(i,indx) = eta(i,indx-1) *                              &
                           qlko_ktcon(i) * g / dp
!cjshiu output this
!          qc_out(i,indx)=dellal(i,indx)
!cjshiu add end
        endif
      enddo


!cjshiu add output for rliq_out
!
!      do k = 1, km
!        do i = 1, im
!         dp = 1000. * del(i,k)
!         rliq_out(i) = rliq_out(i) +  qc_out(i,k) * dp / g
!        enddo
!      enddo
!      rliq_out(:im) = rliq_out(:im) /1000._r8
!cjshiu add end

!cjshiu check
!     write(iulog,*)'run to before mass flux at cloud base for shallow convection'

!c
!c  mass flux at cloud base for shallow convection
!c  (Grant, 2001)
!c
      do i= 1, im
        if(cnvflg(i)) then
          k = kbcon(i)
!         ptem = g*sflx(i)*zi(i,k)/t1(i,1)
          ptem = g*sflx(i)*hpbl(i)/t1(i,1)
          wstar(i) = ptem**h1
          tem = po(i,k)*100. / (rd*t1(i,k))
          xmb(i) = betaw*tem*wstar(i)
          xmb(i) = min(xmb(i),xmbmax(i))
      !    print*, '(Yi-Chi (hp_conv): hpbl, wstar,xmb:',hpbl(i), wstar(i), xmb(i)
        endif
      enddo
!c
!c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!c
      do k = 1, km
        do i = 1, im
          if (cnvflg(i) .and. k .le. kmax(i)) then
            qeso(i,k) = 0.01 * fpvs(t1(i,k))      ! fpvs is in pa
            qeso(i,k) = eps * qeso(i,k) / (pfld(i,k) + epsm1*qeso(i,k))
            val     =             1.e-8
            qeso(i,k) = max(qeso(i,k), val )
          endif
        enddo
      enddo
!c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!c
      do i = 1, im
        delhbar(i) = 0.
        delqbar(i) = 0.
        deltbar(i) = 0.
        delubar(i) = 0.
        delvbar(i) = 0.
        qcond(i) = 0.
      enddo
      do k = 1, km
        do i = 1, im
          if (cnvflg(i)) then
            if(k.gt.kb(i).and.k.le.ktcon(i)) then
              dellat = (dellah(i,k) - hvap * dellaq(i,k)) / cp
              t1(i,k) = t1(i,k) + dellat * xmb(i) * dt2
              q1(i,k) = q1(i,k) + dellaq(i,k) * xmb(i) * dt2
              tem = 1./rcs(i)
              u1(i,k) = u1(i,k) + dellau(i,k) * xmb(i) * dt2 * tem
              v1(i,k) = v1(i,k) + dellav(i,k) * xmb(i) * dt2 * tem
              dp = 1000. * del(i,k)
              delhbar(i) = delhbar(i) + dellah(i,k)*xmb(i)*dp/g
              delqbar(i) = delqbar(i) + dellaq(i,k)*xmb(i)*dp/g
              deltbar(i) = deltbar(i) + dellat*xmb(i)*dp/g
              delubar(i) = delubar(i) + dellau(i,k)*xmb(i)*dp/g
              delvbar(i) = delvbar(i) + dellav(i,k)*xmb(i)*dp/g
            endif
          endif
        enddo
      enddo
      do k = 1, km
        do i = 1, im
          if (cnvflg(i)) then
            if(k.gt.kb(i).and.k.le.ktcon(i)) then
              qeso(i,k) = 0.01 * fpvs(t1(i,k))      ! fpvs is in pa
              qeso(i,k) = eps * qeso(i,k)/(pfld(i,k) + epsm1*qeso(i,k))
              val     =             1.e-8
              qeso(i,k) = max(qeso(i,k), val )
            endif
          endif
        enddo
      enddo
!c
      do i = 1, im
        rntot(i) = 0.
        delqev(i) = 0.
        delq2(i) = 0.
        flg(i) = cnvflg(i)
      enddo
      do k = km, 1, -1
        do i = 1, im
          if (cnvflg(i)) then
            if(k.lt.ktcon(i).and.k.gt.kb(i)) then
              rntot(i) = rntot(i) + pwo(i,k) * xmb(i) * .001 * dt2
            endif
          endif
        enddo
      enddo


!cjshiu check
!     write(iulog,*)'run to before evaporation of rain'
!c
!c evaporating rain
!c
      do k = km, 1, -1
        do i = 1, im
          if (k .le. kmax(i)) then
            deltv(i) = 0.
            delq(i) = 0.
            qevap(i) = 0.
            if(cnvflg(i)) then
              if(k.lt.ktcon(i).and.k.gt.kb(i)) then
                rn(i) = rn(i) + pwo(i,k) * xmb(i) * .001 * dt2
              ! Yi-Chi : for rainout
              rprdsh_out(i,k) = pwo(i,k) * xmb(i)*g/(1000. * del(i,k))
              endif
            endif
            if(flg(i).and.k.lt.ktcon(i)) then
              evef = edt(i) * evfact
              if(slimsk(i).eq.1.) evef=edt(i) * evfactl
!             if(slimsk(i).eq.1.) evef=.07
!c             if(slimsk(i).ne.1.) evef = 0.
              qcond(i) = evef * (q1(i,k) - qeso(i,k))                   &
                       / (1. + el2orc * qeso(i,k) / t1(i,k)**2)
              dp = 1000. * del(i,k)
              if(rn(i).gt.0..and.qcond(i).lt.0.) then
                qevap(i) = -qcond(i) * (1.-exp(-.32*sqrt(dt2*rn(i))))
                qevap(i) = min(qevap(i), rn(i)*1000.*g/dp)
                delq2(i) = delqev(i) + .001 * qevap(i) * dp / g
              endif
              if(rn(i).gt.0..and.qcond(i).lt.0..and.                    &
                delq2(i).gt.rntot(i)) then
                qevap(i) = 1000.* g * (rntot(i) - delqev(i)) / dp
                flg(i) = .false.
              endif
              if(rn(i).gt.0..and.qevap(i).gt.0.) then
                tem  = .001 * dp / g
                tem1 = qevap(i) * tem
                if(tem1.gt.rn(i)) then
                  qevap(i) = rn(i) / tem
                  rn(i) = 0.
                else
                  rn(i) = rn(i) - tem1
                endif
                q1(i,k) = q1(i,k) + qevap(i)
                t1(i,k) = t1(i,k) - elocp * qevap(i)
                deltv(i) = - elocp*qevap(i)/dt2
                delq(i) =  + qevap(i)/dt2
                delqev(i) = delqev(i) + .001*dp*qevap(i)/g
              endif
! Yi-Chi:
              rprdsh_out(i,k)  = rprdsh_out(i,k) - qevap(i)/dt2
!cjshiu add output for tendency of evaporation of precipitation
              evapcv_out(i,k) = qevap(i)/dt2
!cjs              write(iulog,*)'i= k= evapcv_out= qevap(i)= ',i,k,evapcv_out(i,k),qevap(i)
!cjs              print *,'i= k= evapcv_out= qevap(i)= ',i,k,evapcv_out(i,k),qevap(i)
!cjshiu addend              
              dellaq(i,k) = dellaq(i,k) + delq(i) / xmb(i)
              delqbar(i) = delqbar(i) + delq(i)*dp/g
              deltbar(i) = deltbar(i) + deltv(i)*dp/g
            endif
          endif
        enddo
      enddo
!cj
!     do i = 1, im
!!     if(me.eq.31.and.cnvflg(i)) then
!     if(cnvflg(i)) then
!       print *, ' shallow delhbar, delqbar, deltbar = ', &
!                 delhbar(i),hvap*delqbar(i),cp*deltbar(i) 
!       print *, ' shallow delubar, delvbar = ',delubar(i),delvbar(i)
!       print *, ' precip =', hvap*rn(i)*1000./dt2
!       print*,'pdif= ',pfld(i,kbcon(i))-pfld(i,ktcon(i))
!     endif
!     enddo
!cj

      do i = 1, im
        if(cnvflg(i)) then
          if(rn(i).lt.0..or..not.flg(i)) rn(i) = 0.
          ktop(i) = ktcon(i)
          kbot(i) = kbcon(i)
          kcnv(i) = 0
!cjshiu given value to cnt and cnb
!          cnt_out(i)=ktop(i)*1.
!          cnb_out(i)=kbot(i)*1.
!cjs          write(iulog,*)' i= ktop(i)= kbot(i)= cnt_out(i)= cnb_out(i)= ',i,ktop(i),kbot(i),cnt_out(i),cnb_out(i)
!cjshiu add end
        endif
      enddo

      delqlbar = 0.0
!c
!c  cloud water
!c
      if (ncloud.gt.0) then
!
      do k = 1, km1
        do i = 1, im
          if (cnvflg(i)) then
            if (k.gt.kb(i).and.k.le.ktcon(i)) then
              ! Yi-Chi +++
              icwmr_out(i,k) = 0.5*(icwmr_tmp1(i,k)+icwmr_tmp1(i,k-1)) &
                          * g / (1000. * del(i,k))
              icwmr_out(i,k) = icwmr_out(i,k) * xmb(i) * dt2
              ! Yi-Chi ---
              tem  = dellal(i,k) * xmb(i) * dt2
              tem1 = max(0.0, min(1.0, (tcr-t1(i,k))*tcrf))
              if (ql(i,k,2) .gt. -999.0) then
                ql(i,k,1) = ql(i,k,1) + tem * tem1            ! ice
                ql(i,k,2) = ql(i,k,2) + tem *(1.0-tem1)       ! water
                ! Yi-Chi : detrained liquid water
                qc_out(i,k) = tem*(1.0-tem1)/dt2
                !-------------------------
              else
                ! +++ Yi-Chi add for debugging
                dp = 1000. * del(i,k)
                delqlbar(i) = delqlbar(i)+tem*dp/g
                ! ---
                ql(i,k,1) = ql(i,k,1) + tem
                qc_out(i,k) = tem/dt2
              endif
!cjshiu add constituent tendency
! note the present formulation is wrong (need to use sort of Ptend%q= ql - xxx
!cjs              do m = 1, pcnst
!cjs              do m = 1, 2
!cjshiu              do m = 1, 2
!cjs               qten_out(i,k,m)  = ql(i,k,m) / dt2
!cjs check this for debugging float invalid
!cjs                write(iulog,*)'k= i= ql(i,k,m)= qten_out(i,k,m)=',k,i,m,ql(i,k,m),qten_out(i,k,m)
!cjshiu                write(iulog,*)'k= i= ql(i,k,m)= tem= tem1= ',k,i,m,ql(i,k,m),tem,tem1
!cjshiu              enddo
!cjshiu add end
            endif
          endif
        enddo
      enddo
!
      endif

!cjshiu add output for rliq_out
!
      do k = 1, km
        do i = 1, im
         dp = 1000. * del(i,k)
         rliq_out(i) = rliq_out(i) +  qc_out(i,k) * dp / g
        enddo
      enddo
      rliq_out(:im) = rliq_out(:im) /1000._r8
!cjshiu add end


! +++ ycw +++
!    do i = 1, im
!    if(cnvflg(i)) then
!      write(iulog,*) 'dt2',dt2
!      write(iulog,*) ' shallow delqlbar ', hvap*delqlbar(i)/dt2
!      write(iulog,*) ' precip =', hvap*rn(i)*1000./dt2
!      write(iulog,*) ' precip + delqlbar =', hvap*rn(i)*1000./dt2 + hvap*delqlbar(i)/dt2
!      write(iulog,*) ' difq =', hvap*rn(i)*1000./dt2 + hvap*delqlbar(i)/dt2 + hvap*delqbar(i)
!      !write(iulog,*) 'dlfsas =', dlfsas
!      !write(iulog,*) 'ktcon:',ktcon(i)
!      !write(iulog,*) 'qlko_ktcon:',qlko_ktcon(i)
!    endif
!    enddo
! --- ycw ---

!
! hchuang code change
!
      do k = 1, km
        do i = 1, im
          if(cnvflg(i)) then
            if(k.ge.kb(i) .and. k.lt.ktop(i)) then
              ud_mf(i,k) = eta(i,k) * xmb(i) * dt2
              ! Yi-Chi : mass flux of updraft
              cmfmc_out(i,k) = ud_mf(i,k)*g/(100.*dt2)
            endif
          endif
        enddo
      enddo
      do i = 1, im
        if(cnvflg(i)) then
           k = ktop(i)-1
           dt_mf(i,k) = ud_mf(i,k)
        endif
      enddo
!!
      return
  end subroutine shalcnv

end module hp_conv
