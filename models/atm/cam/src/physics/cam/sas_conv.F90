
module sas_conv

!---------------------------------------------------------------------------------
! Purpose:
!
! Interface from Zhang-McFarlane convection scheme, includes evaporation of convective
! precip from the ZM scheme
! Oct, 2014:
! - port to CAM5.3
! Nov 26, 2013:
! - output heo1 and heso1 on mid-layer for comparison
! Nov 25, 2013:
! - output heo and heso on interface
! Nov 21, 2013:
! - add CIN for distance
! Nov 14, 2013:
! - add cnvflg mask for SAS triggering
!   flag = 1 for 
! Oct 12, 2013:
! - remove zmevap subroutine
! Oct 9, 2013:
!  - as qlk at interface level, qlic is interprelated onto mid-layer level.
! Sept 30, 2013:
!  - correct ql 
! Sept 23, 2013:
!  - no overshoot + c1half
! August 19, 2013 :
!  - recalculate the cloud work function with new approximation
!    of si - sbar and qi - qbar.
!  - correct an error of 2nd derivative calculation
! August 15, 2013 :
!  - add qrch 2nd derivative into the calculation of updraft and downdraft
! June 22, 2013 :
!  - test for the new modification of GFS
!  - /nuwa_cluster/home/ychwang/PROGRAMS/CODE_gfs_cpscheme_june2013/
! March 22, 2013:
! - revise the method for calculating cloud top overshooting.
!   add one condition for 
! March 20, 2013:
!   to test the sensitivity of function for saturated vapor pressure,
!   rewrite the function to alternatives.
! Oct 2, 2012  : SAS only do calculations for 
! Sept 3, 2012 : add temporary variables to ensure detrain/entrain output 
!                are only shown as SAS is activated.
! Sept 1, 2012 : add output for the detrained moisture
! August 31, 2012: test entrainment rate by adding a coefficient as in EC 
!        
! August 27, 2012: test the sensitivity 
!                   cxlamu = 0.5*10e-4
! August 14,2012: Yi-Chi : add output of updraft moisture properties
! July 9 2012: Yi-Chi : change the declaration of ideep from ncol to pcols
!                       for all I/O variables of subroutine SASCNVN
! July 4 2012: Yi-Chi : use Bolton's formulation for saturated vapor pres 
!                    to replace subroutine fpvs
! July2012: Yi-Chi : remove evaporation(qevap) from rainout at all levels
! Apr 2012: Yi-Chi : include the partition of ice and cloud water
!                      for the detrained cloud water
! Apr 2012: Yi-Chi : include convective transport (subroutine convtran_sas)
! Mar 2012: Yi-Chi : most ZM related variables are kept to
!                    Future work can include namelist variables and develop
!                    input in namelist
! Feb 2012: Yi-Chi : The scheme is modified to accomodate SAS scheme.
!
!
!---------------------------------------------------------------------------------
  use shr_kind_mod,    only: r8 => shr_kind_r8
  use spmd_utils,      only: masterproc
  use ppgrid,          only: pcols, pver, pverp
!  use cldwat,          only: cldwat_fice
  use physconst,       only: cpair, epsilo, gravit, latice, latvap, tmelt, rair, &
                             cpwv, cpliq, rh2o
  use abortutils,      only: endrun
  use cam_logfile,     only: iulog
  use cam_history,  only: outfld
  implicit none

  save
  private                         ! Make default type private to the module
!
! PUBLIC: interfaces
!
  public zmconv_readnl            ! read zmconv_nl namelist
  public sas_convi                 ! initialization of SAS scheme
! Yi-Chi
  public sascnvn                    ! core of SAS scheme
!  public sas_convr                 ! core of SAS scheme
   public convtran_sas             ! convective transport for SAS scheme

!
! Private data
!
   real(r8), parameter :: unset_r8 = huge(1.0_r8)
   real(r8) :: zmconv_c0_lnd = unset_r8
   real(r8) :: zmconv_c0_ocn = unset_r8
   real(r8) :: zmconv_ke     = unset_r8

   real(r8) rl         ! wg latent heat of vaporization.
   real(r8) cpres      ! specific heat at constant pressure in j/kg-degk.
   real(r8), parameter :: capelmt = 70._r8  ! threshold value for cape for deep convection.
   real(r8) :: ke           ! Tunable evaporation efficiency set from namelist input zmconv_ke
   real(r8) :: c0_lnd       ! set from namelist input zmconv_c0_lnd
   real(r8) :: c0_ocn       ! set from namelist input zmconv_c0_ocn
   real(r8) tau   ! convective time scale
   real(r8),parameter :: a = 21.656_r8
   real(r8),parameter :: b = 5418._r8
   real(r8),parameter :: c1e = 6.112_r8
   real(r8),parameter :: c2e = 17.67_r8
   real(r8),parameter :: c3e = 243.5_r8
   real(r8) :: tfreez
   real(r8) :: eps1


   logical :: no_deep_pbl ! default = .false.
                          ! no_deep_pbl = .true. eliminates deep convection entirely within PBL


!moved from moistconvection.F90
   real(r8) :: rgrav       ! reciprocal of grav
   real(r8) :: rgas        ! gas constant for dry air
   real(r8) :: grav        ! = gravit
   real(r8) :: cp          ! = cpres = cpair

   integer  limcnv       ! top interface level limit for convection

   real(r8),parameter ::  tiedke_add = 0.5_r8

contains

subroutine zmconv_readnl(nlfile)

   use namelist_utils,  only: find_group_name
   use units,           only: getunit, freeunit
   use mpishorthand

   character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

   ! Local variables
   integer :: unitn, ierr
   character(len=*), parameter :: subname = 'zmconv_readnl'

   namelist /zmconv_nl/ zmconv_c0_lnd, zmconv_c0_ocn, zmconv_ke
   !-----------------------------------------------------------------------------

   if (masterproc) then
      unitn = getunit()
      open( unitn, file=trim(nlfile), status='old' )
      call find_group_name(unitn, 'zmconv_nl', status=ierr)
      if (ierr == 0) then
         read(unitn, zmconv_nl, iostat=ierr)
         if (ierr /= 0) then
            call endrun(subname // ':: ERROR reading namelist')
         end if
      end if
      close(unitn)
      call freeunit(unitn)

      ! set local variables
      c0_lnd = zmconv_c0_lnd
      c0_ocn = zmconv_c0_ocn
      ke = zmconv_ke

   end if

#ifdef SPMD
   ! Broadcast namelist variables
   call mpibcast(c0_lnd,            1, mpir8,  0, mpicom)
   call mpibcast(c0_ocn,            1, mpir8,  0, mpicom)
   call mpibcast(ke,                1, mpir8,  0, mpicom)
#endif

end subroutine zmconv_readnl

!-------------------------------------------------------
!-------------------------------------------------------
subroutine sas_convi(limcnv_in, no_deep_pbl_in)

   use dycore,       only: dycore_is, get_resolution

   integer, intent(in)           :: limcnv_in       ! top interface level limit for convection
   logical, intent(in), optional :: no_deep_pbl_in  ! no_deep_pbl = .true. eliminates ZM convection entirely within PBL

   ! local variables
   character(len=32)   :: hgrid           ! horizontal grid specifier

   ! Initialization of ZM constants
   limcnv = limcnv_in
   tfreez = tmelt
   eps1   = epsilo
   rl     = latvap
   cpres  = cpair
   rgrav  = 1.0_r8/gravit
   rgas   = rair
   grav   = gravit
   cp     = cpres

   if ( present(no_deep_pbl_in) )  then
      no_deep_pbl = no_deep_pbl_in
   else
      no_deep_pbl = .false.
   endif

   ! tau=4800. were used in canadian climate center. however, in echam3 t42,
   ! convection is too weak, thus adjusted to 2400.

   hgrid = get_resolution()
   tau = 3600._r8

   if ( masterproc ) then
      write(iulog,*) 'tuning parameters zm_convi: tau',tau
      write(iulog,*) 'tuning parameters zm_convi: c0_lnd',c0_lnd, ', c0_ocn', c0_ocn
      write(iulog,*) 'tuning parameters zm_convi: ke',ke
      write(iulog,*) 'tuning parameters zm_convi: no_deep_pbl',no_deep_pbl
   endif

   if ( masterproc ) write(iulog,*)'**** SAS: After initialization ****'

end subroutine sas_convi

!-----------------------------------------------------------------
! Yi-Chi : March 7, 2012
!          add rainfall, evaporation, detraining mass flux
!-----------------------------------------------------------------
subroutine sascnvn(lchnk, im  ,ix  ,km  , &
                   jcap , delt,del ,prsl,ps ,phil  , &
                   ql   , qlic, &!ql,  &
                   q1   , t1  ,u1  ,v1  ,rcs,cldwrk, &
                   rn   , kbot,ktop,slimsk, &
           dot,ncloud,ud_mf,dd_mf, &
           rainout,qevapout,zdusas,dlfsas, &!)
           ideep,lengath, & !)
           euout,duout,edout,ddout)
!           ,qckoout,duqout,flgmask,cinmask, &!)
!           heoout, hesoout,heoout1,hesoout1, &!)
!           xmbout, xkout, fldout,dtconvout, aa1crit,xomega)
!           dot,ncloud,ud_mf,dd_mf)
!    &     dot,ncloud,ud_mf,dd_mf,dt_mf,me)
!
!--------------------------------------
! YI-CHI: FOLLOW WRF, add three files for the data type and physical constant
!-------------------------------------------------------------------
      use ppgrid,          only: pcols
      use module_gfs_machine ,  only : kind_phys
!      use module_gfs_funcphys, only : fpvs
      use physconst,       only: grav=>gravit,cp=>cpair, hvap=> latvap, &
                                 rv=>rh2o, rair, t0c=>tmelt, &
                                 cvap=> cpwv, cliq=>cpliq
                                 
!-------------------------------------------------------
      implicit none
!-------------------------------------------------------
!
      integer            im, ix,  km, jcap, ncloud,   &
                         kbot(pcols), ktop(pcols), kcnv(pcols)
!                        kbot(im), ktop(im), kcnv(im) 
!    &,                  me
      ! ycw add for diagnostic output
      integer            lchnk
      ! ycw add END
      real(kind=kind_phys) delt
      real(kind=kind_phys) ps(pcols),     del(pcols,km),  prsl(pcols,km), &
                          ql(pcols,km,2),q1(pcols,km),   t1(pcols,km),   &
                          u1(pcols,km),  v1(pcols,km),   rcs(pcols),     &
                          cldwrk(pcols), & ! output
                          rn(pcols),      slimsk(pcols), &
                          dot(pcols,km), phil(pcols,km), &
! hchuang code change mass flux output
                          ud_mf(pcols,km),dd_mf(pcols,km),dt_mf(pcols,km), &
! Yi-Chi : output for CAM5 Mar 2012
                          rainout(pcols,km), qevapout(pcols,km), &
                          dlfsas(pcols,km),  zdusas(pcols,km), &
! Yi-Chi : detrain, entrain mass flux
                          euout(pcols,km), duout(pcols,km),qckoout(pcols,km), &  ! updraft entrainment
                          edout(pcols,km), ddout(pcols,km), &  ! downdraft detrainment
                          duqout(pcols,km), &
                          qlic(pcols,km),qlic_tmp(pcols,km) ! in-cloud mixing ratio
! Yi-Chi (August 2012)
      real(kind=kind_phys) eutmp(pcols,km), dutmp(pcols,km), &  ! updraft entrainment
                           edtmp(pcols,km), ddtmp(pcols,km), &  ! downdraft detrainment
                           duqtmp(pcols,km)
!
      integer              i, j, indx, jmn, k, kk, latd, lond, km1
!
      real(kind=kind_phys) clam, cxlamu, xlamde, xlamdd, &
                           c0xlamu
! 
      real(kind=kind_phys) adw,     aup,     aafac,  &
                          beta,    betal,   betas,  &
                          c0,      cpoel,   dellat,  delta, &
                          desdt,   deta,    detad,   dg,    &
                          dh,      dhh,     dlnsig,  dp,    &
                          dq,      dqsdp,   dqsdt,   dt,    &
                          dt2,     dtmax,   dtmin,   dv1h,  &
                          dv1q,    dv2h,    dv2q,    dv1u,  &
                          dv1v,    dv2u,    dv2v,    dv3q,  &
                          dv3h,    dv3u,    dv3v,           &
                          dz,      dz1,     e1,      edtmax,&
                          edtmaxl, edtmaxs, el2orc,  elocp, &
                          es,      etah,    cthk,    dthk,  &
                          evef,    evfact,  evfactl, fact1, &
                          fact2,   factor,  fjcap,   fkm,   &
                          g,       gamma,   pprime,         &
                          qlk,     qrch,    qs,      c1,    &
                          rain,    rfact,   shear,   tem1,  &
                          tem2,    terr,    val,     val1,  &
                          val2,    w1,      w1l,     w1s,   &
                          w2,      w2l,     w2s,     w3,    &
                          w3l,     w3s,     w4,      w4l,   &
                          w4s,     xdby,    xpw,     xpwd,  &
                          xqrch,   mbdt,    tem,            &
                          ptem,    ptem1,   pgcon
!
      integer              kb(im), kbcon(im), kbcon1(im),    &
                          ktcon(im), ktcon1(im),            &
                          jmin(im), lmin(im), kbmax(im),    &
                          kbm(im), kmax(im)
!
      real(kind=kind_phys) aa1(im),     acrt(im),   acrtfct(im), &
                          delhbar(im), delq(im),   delq2(im),   &
                          delqbar(im), delqev(im), deltbar(im), &
                          deltv(im),   dtconv(im), edt(im),     &
                          edto(im),    edtx(im),   fld(im),     &
                          hcdo(im,km), hmax(im),   hmin(im),    &
                          ucdo(im,km), vcdo(im,km),aa2(im),     &
                          pbcdif(im),  pdot(im),   po(im,km),   &
                          pwavo(im),   pwevo(im),  xlamud(im),  &
                          qcdo(im,km), qcond(im),  qevap(im),   &
                          rntot(im),   vshear(im), xaa0(im),    &
                          xk(im),      xlamd(im),               &
                          xmb(im),     xmbmax(im), xpwav(im),   &
                          xpwev(im),   delubar(im),delvbar(im)
! +++ Yi-Chi
      real(kind=kind_phys) delqlbar(im)
! ---

!j
      real(kind=kind_phys) cincr, cincrmax, cincrmin
!j

!  physical parameters
      ! +++ Yi-Chi
     real(kind=kind_phys) eps,epsm1,fv
!      real(kind=kind_phys) grav,cp,hvap,rv,t0c,cvap,cliq,eps,epsm1,fv
!                   rv => con_rv, fv => con_fvirt, t0c => con_t0c, &
!                   cvap => con_cvap, cliq => con_cliq, &
!                   eps => con_eps, epsm1 => con_epsm1

!!      parameter(grav=gravit,cp=cpair, hvap=latvap)
!!      parameter(rv=rh2o, t0c=tmelt, cvap= cpwv, cliq=cpliq)
!      parameter(g=grav)
!      parameter(cpoel=cp/hvap,elocp=hvap/cp, &
!               el2orc=hvap*hvap/(rv*cp))
      ! +++ Yi-Chi
!      parameter(fv=(rh2o/rair)-1, eps=rair/rh2o, epsm1=(rair/rh2o)-1)
      !parameter(terr=0.,c0=.002,c1=.001,delta=fv)
      !parameter(terr=0.,c0=.002,c1=0.0)
      parameter(terr=0.,c0=.002,c1=.002) ! default
      !parameter(terr=0.,c0=.002,c1=.002,delta=fv) ! default
      ! ---
!      parameter(fact1=(cvap-cliq)/rv,fact2=hvap/rv-fact1*t0c)
      parameter(cthk=150.,cincrmax=180.,cincrmin=120.,dthk=25.)
!  local variables and arrays
      real(kind=kind_phys) pfld(im,km),    to(im,km),     qo(im,km), &
                          uo(im,km),      vo(im,km),     qeso(im,km)
!  cloud water
      real(kind=kind_phys) qlko_ktcon(im), dellal(im,km), tvo(im,km), &
                          dbyo(im,km),    zo(im,km),     xlamue(im,km), &
                          fent1(im,km),   fent2(im,km),  frh(im,km),    &
                          heo(im,km),     heso(im,km),                  &
                          qrcd(im,km),    dellah(im,km), dellaq(im,km), &
                          dellau(im,km),  dellav(im,km), hcko(im,km),   &
                          ucko(im,km),    vcko(im,km),   qcko(im,km),   &
                          eta(im,km),     etad(im,km),   zi(im,km),     &
                          qrcdo(im,km),   pwo(im,km),    pwdo(im,km),   &
! +++ Yi-Chi : June 22, 2013
                          tx1(im),        sumx(im),      qrcko(im,km)
!                          tx1(im),        sumx(im)
! ---
!    &,                    rhbar(im)
!
      logical totflg, cnvflg(im), flg(im)
!     Yi-Chi
      real(kind=kind_phys) dqaa1, dsaa1 ! temp variables for cloud work function calculation
      logical flgfvps
      real(kind=kind_phys) tfreez !July 2012 to calculate Bolton's saturated vapor pressure
!      parameter(tfreez=tmelt)
! July 9, 2012: change the index of ideep to pcols
!      integer ideep(im) ! for CAM5: indices of grids with active deep convect 
      integer ideep(pcols) ! for CAM5: indices of grids with active deep convect
      integer lengath   ! for CAM5: number of grids with active deep convect
! Nov 14, 2013 : write out flag mask for triggering
      real(kind=kind_phys) flgmask(im), cinmask(im)
      ! heo and heso interpolated on interface
      real(kind=kind_phys) heoout(im,km), hesoout(im,km)
      ! heo and heso on mid-layer
      real(kind=kind_phys) heoout1(im,km), hesoout1(im,km)
      ! closure-related variables
      real(kind=kind_phys) xmbout(im), xkout(im), fldout(im),dtconvout(im), aa1crit(im),xomega(im)


      real(kind=kind_phys) pcrit(15), acritt(15), acrit(15)
!     save pcrit, acritt
      data pcrit/850.,800.,750.,700.,650.,600.,550.,500.,450.,400., &
                350.,300.,250.,200.,150./
      data acritt/.0633,.0445,.0553,.0664,.075,.1082,.1521,.2216, &
                .3151,.3677,.41,.5255,.7663,1.1686,1.6851/
!  gdas derived acrit
!     data acritt/.203,.515,.521,.566,.625,.665,.659,.688,
!    &            .743,.813,.886,.947,1.138,1.377,1.896/
      real(kind=kind_phys) tf, tcr, tcrf
      parameter (tf=233.16, tcr=263.16, tcrf=1.0/(tcr-tf))

! ++++++++++++++++++++
      g = grav
      cpoel=cp/hvap
      elocp=hvap/cp
      el2orc=hvap*hvap/(rv*cp)
      fv=(rv/rair)-1 
      eps=rair/rv
      epsm1=(rair/rv)-1
      delta = fv
      fact1=(cvap-cliq)/rv
      fact2=hvap/rv-fact1*t0c

!
!-----------------------------------------------------------------------
!
      km1 = km - 1
! +++ Yi-Chi : check all inputs : Nov 2013
!     write(iulog,*) 'Yi-Chi(sas_conv) land mask (slimsk)',slimsk
!     write(iulog,*) 'Yi-Chi(sas_conv) mid-layer pres levels (prsl)',prsl
!     write(iulog,*) 'Yi-Chi(sas_conv) pressure depth (del)',del
!     write(iulog,*) 'Yi-Chi(sas_conv) omega (dot)',dot
!     write(iulog,*) 'Yi-Chi(sas_conv) geopotential (phil)',phil
!     write(iulog,*) 'Yi-Chi(sas_conv) surface prs (ps)',ps
!     write(iulog,*) 'Yi-Chi(sas_conv) temp input (t1)',t1
!     write(iulog,*) 'Yi-Chi(sas_conv) moisture input (q1)',q1
!     write(iulog,*) 'Yi-Chi(sas_conv) zonal (u1)',u1
!     write(iulog,*) 'Yi-Chi(sas_conv) meridional (v1)',v1

! +++Yi-Chi (March 2013)
!  set a flag for vapor pressure calculation
      flgfvps = .true.
!      tfreez = tmelt
!
!  initialize arrays
!
      do i=1,im
        cnvflg(i) = .true.
        rn(i)=0.
        kbot(i)=km+1
        ktop(i)=0
        kbcon(i)=km
        ktcon(i)=1
        dtconv(i) = 3600.
        cldwrk(i) = 0.
        pdot(i) = 0.
        pbcdif(i)= 0.
        lmin(i) = 1
        jmin(i) = 1
        qlko_ktcon(i) = 0.
        edt(i)  = 0.
        edto(i) = 0.
        edtx(i) = 0.
        acrt(i) = 0.
        acrtfct(i) = 1.
        aa1(i)  = 0.
        aa2(i)  = 0.
        xaa0(i) = 0.
        pwavo(i)= 0.
        pwevo(i)= 0.
        xpwav(i)= 0.
        xpwev(i)= 0.
        vshear(i) = 0.
      enddo
! Yi-Chi : initialize arrays
    do k = 1, km
      do i = 1, im
          rainout(i,k) = 0. ! rain production at each level
          qevapout(i,k)= 0. ! rain evaporation in 
          zdusas(i,k)  = 0. ! detraining mass flux(xlam)
          dlfsas(i,k)  = 0. ! detraining liquid water(dellal*xmb*dt2)
          !  output entrainment/detrainment rate for convective transport
          !  unit : fraction of unit mass flux
          euout(i,k) = 0.  ! updraft entrainment
          duout(i,k) = 0.  ! updraft detrainment
          edout(i,k) = 0.  ! downdraft entrainment
          ddout(i,k) = 0.  ! downdraft detrainment
          duqout(i,k)= 0.  ! output detrained moisture
          ! updraft moisture
          qckoout(i,k) = 0. ! updraft specific humidity
          ! temporary arrays
          eutmp(i,k) = 0.  ! updraft entrainment
          dutmp(i,k) = 0.  ! updraft detrainment
          edtmp(i,k) = 0.  ! downdraft entrainment
          ddtmp(i,k) = 0.  ! downdraft detrainment
          duqtmp(i,k)= 0.  ! output detrained moisture
          ! updraft moisture
          !qckotmp(i,k) = 0. ! updraft specific humidity
          ! in-cloud liquid water : Sept 2013
          qlic(i,k)= 0.
          qlic_tmp(i,k)= 0.
          heoout(i,k) = 0.
          hesoout(i,k) = 0.
          heoout1(i,k) = 0.
          hesoout1(i,k) = 0.
      enddo
    enddo
    do i = 1, im
      ideep(i) = 0
      ! +++ Yi-Chi :
      flgmask(i) = 0. ! initialize triggering flag mask
      cinmask(i) = 0. ! initialize cin mask
      ! ---
      xmbout(i)   = 0.
      xkout(i)    = 0.
      fldout(i)   = 0.
      dtconvout(i)= 0.
      aa1crit(i)  = 0.
      xomega(i)   = 0.
    enddo
      lengath  = 0

! hchuang code change
      do k = 1, km
        do i = 1, im
          ud_mf(i,k) = 0.
          dd_mf(i,k) = 0.
          dt_mf(i,k) = 0.
        enddo
      enddo
!
      do k = 1, 15
        acrit(k) = acritt(k) * (975. - pcrit(k))
      enddo
      dt2 = delt
      val   =         1200.
      dtmin = max(dt2, val )
      val   =         3600.
      dtmax = max(dt2, val )
!  model tunable parameters are all here
      mbdt    = 10.
      edtmaxl = .3
      edtmaxs = .3
      clam    = .1
      ! +++ Yi-Chi : Aug 2013
      !aafac   = 0.  ! noovershoot
      aafac   = .1 ! default
      ! ---
!     betal   = .15
!     betas   = .15
      betal   = .05
      betas   = .05
!     evef    = 0.07
      evfact  = 0.3
      evfactl = 0.3
!     Yi-Chi : test the sensitivity of cloud cover to
!              the entrainment rate;
!              cxlamu becomes half implying less environmental
!              influence on entrainment rate
      cxlamu  = 1.0e-4
      !cxlamu  = 0.5*1.0e-4
      !cxlamu  = 1.0e-5
      !cxlamu   = 0.0*1.0e-4
! Yi-Chi :
!     add one coefficient for scaling down
!      c0xlamu   = 1.0e-4
!----------------------------
      xlamde  = 1.0e-4
      xlamdd  = 1.0e-4
!
!     pgcon   = 0.7     ! Gregory et al. (1997, QJRMS)
      pgcon   = 0.55    ! Zhang & Wu (2003,JAS)
      fjcap   = (float(jcap) / 126.) ** 2
      val     =           1.
      fjcap   = max(fjcap,val)
      fkm     = (float(km) / 28.) ** 2
      fkm     = max(fkm,val)
      w1l     = -8.e-3 
      w2l     = -4.e-2
      w3l     = -5.e-3 
      w4l     = -5.e-4
      w1s     = -2.e-4
      w2s     = -2.e-3
      w3s     = -1.e-3
      w4s     = -2.e-5
!
!  define top layer for search of the downdraft originating layer
!  and the maximum thetae for updraft
!
      do i=1,im
        kbmax(i) = km
        kbm(i)   = km
        kmax(i)  = km
        tx1(i)   = 1.0 / ps(i)
      enddo
!     
      do k = 1, km
        do i=1,im
          if (prsl(i,k)*tx1(i) .gt. 0.04) kmax(i)  = k + 1
          if (prsl(i,k)*tx1(i) .gt. 0.45) kbmax(i) = k + 1
          if (prsl(i,k)*tx1(i) .gt. 0.70) kbm(i)   = k + 1
        enddo
      enddo
      do i=1,im
      ! +++ Yi-Chi : June 22, 2013
        kmax(i)  = min(km,kmax(i))
      ! ---
        kbmax(i) = min(kbmax(i),kmax(i))
        kbm(i)   = min(kbm(i),kmax(i))
      enddo
      !write(iulog,*) 'Yi-Chi:kmax,kbmax,kbm=',kmax,kbmax,kbm
!
!  hydrostatic height assume zero terr and initially assume
!    updraft entrainment rate as an inverse function of height 
!
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
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   convert surface pressure to mb from cb
!
      do k = 1, km
        do i = 1, im
          if (k .le. kmax(i)) then
            pfld(i,k) = prsl(i,k) * 10.0
            eta(i,k)  = 1.
            fent1(i,k)= 1.
            fent2(i,k)= 1.
            frh(i,k)  = 0.
            hcko(i,k) = 0.
            qcko(i,k) = 0.
            ! +++ Yi-Chi : June 22, 2013
            qrcko(i,k)= 0.
            ! ---
            ucko(i,k) = 0.
            vcko(i,k) = 0.
            etad(i,k) = 1.
            hcdo(i,k) = 0.
            qcdo(i,k) = 0.
            ucdo(i,k) = 0.
            vcdo(i,k) = 0.
            qrcd(i,k) = 0.
            qrcdo(i,k)= 0.
            dbyo(i,k) = 0.
            pwo(i,k)  = 0.
            pwdo(i,k) = 0.
            dellal(i,k) = 0.
            to(i,k)   = t1(i,k)
            qo(i,k)   = q1(i,k)
            uo(i,k)   = u1(i,k) * rcs(i)
            vo(i,k)   = v1(i,k) * rcs(i)
          endif
        enddo
      enddo

      !print *, 'Yi-Chi: sas_conv586:to',to
      !print *, 'Yi-Chi: sas_conv586:qo',qo
      !print *, 'Yi-Chi: sas_conv586:uo',uo
      !print *, 'Yi-Chi: sas_conv586:vo',vo

!
!  column variables
!  p is pressure of the layer (mb)
!  t is temperature at t-dt (k)..tn
!  q is mixing ratio at t-dt (kg/kg)..qn
!  to is temperature at t+dt (k)... this is after advection and turbulan
!  qo is mixing ratio at t+dt (kg/kg)..q1
!
      do k = 1, km
        do i=1,im
          if (k .le. kmax(i)) then
            ! Yi-Chi (March 2013)
            qeso(i,k) = esat(to(i,k),flgfvps)
            ! Yi-Chi: change the saturated vapor pressure (July 2012)
            !esat = c1e*exp(c2e*(to(i,k)-tfreez)/(c3e+to(i,k)-tfreez))       ! esat(T) in mb
            !qeso(i,k) = 0.01 * fpvs(to(i,k))      ! fpvs is in pa
            !qeso(i,k) = esat
            !-----Yi-Chi--------
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
!
!  compute moist static energy
!
      do k = 1, km
        do i=1,im
          if (k .le. kmax(i)) then
!           tem       = g * zo(i,k) + cp * to(i,k)
            tem       = phil(i,k) + cp * to(i,k)
            heo(i,k)  = tem  + hvap * qo(i,k)
            heso(i,k) = tem  + hvap * qeso(i,k)
            ! +++ Yi-Chi
            heoout1(i,k)   = heo(i,k)
            hesoout1(i,k)  = heso(i,k)
            ! ---
!           heo(i,k)  = min(heo(i,k),heso(i,k))
          endif
        enddo
      enddo
!
!  determine level with largest moist static energy
!  this is the level where updraft starts
!
      do i=1,im
        hmax(i) = heo(i,1)
        kb(i)   = 1
      enddo
      do k = 2, km
        do i=1,im
          if (k .le. kbm(i)) then
            if(heo(i,k).gt.hmax(i)) then
              kb(i)   = k
              hmax(i) = heo(i,k)
            endif
          endif
        enddo
      enddo
!    Yi-Chi: calculate to,qo,po at mid-layer
      do k = 1, km1
        do i=1,im
          if (k .le. kmax(i)-1) then
            dz      = .5 * (zo(i,k+1) - zo(i,k))
            dp      = .5 * (pfld(i,k+1) - pfld(i,k))
            ! Yi-Chi (March 2013)
            es = esat(to(i,k+1),flgfvps)
            ! Yi-Chi (July 2012)
            !esat = c1e*exp(c2e*(to(i,k+1)-tfreez)/(c3e+to(i,k+1)-tfreez))       ! esat(T) in mb
            !es   = esat ! mb->Pa
            !write(iulog,*) 'Yi-Chi(sas_conv)esat,to',esat,to(i,k+1)-tfreez
            !es      = 0.01 * fpvs(to(i,k+1))      ! fpvs is in pa
            !---------------------------
            !write(iulog,*) 'fpvs,',fpvs(to(i,k+1))
            pprime  = pfld(i,k+1) + epsm1 * es
            qs      = eps * es / pprime
            dqsdp   = - qs / pprime
            desdt   = es * (fact1 / to(i,k+1) + fact2 / (to(i,k+1)**2))
            !write(iulog,*) 'pfld,',pfld(i,k+1),'desdt,',desdt,'pprime',pprime
            !write(iulog,*) 'es,',es,'to',to(i,k+1),'qs',qs
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
!     Yi-Chi: calculate q0,heo,heso,uo,vo at the mid-layers in subcloud
      do k = 1, km1
        do i=1,im
          if (k .le. kmax(i)-1) then
            ! Yi-Chi (March 2013)
            qeso(i,k) = esat(to(i,k),flgfvps)
            ! Yi-Chi (July 2012)
            !esat = c1e*exp(c2e*(to(i,k)-tfreez)/(c3e+to(i,k)-tfreez))       ! esat(T) in mb
            !qeso(i,k) = 0.01*esat*100 ! mb->Pa
            !write(iulog,*) 'Yi-Chi(sas_conv606)tfreez',tfreez
            !write(iulog,*) 'Yi-Chi(sas_conv606)esat,to',esat,to(i,k)-tfreez
            !qeso(i,k) = 0.01 * fpvs(to(i,k))      ! fpvs is in pa
            !------------------------
            !write(iulog,*) 'Yi-Chi(sas_conv732)esat',qeso(i,k)
            qeso(i,k) = eps * qeso(i,k) / (po(i,k) + epsm1*qeso(i,k))
            !write(iulog,*) 'Yi-Chi(sas_conv732)qeso,po',qeso(i,k),po(i,k)
            val1      =             1.e-8
            qeso(i,k) = max(qeso(i,k), val1)
            val2      =           1.e-10
            qo(i,k)   = max(qo(i,k), val2 )
!           qo(i,k)   = min(qo(i,k),qeso(i,k))
            frh(i,k)  = 1. - min(qo(i,k)/qeso(i,k), 1.)
            heo(i,k)  = .5 * g * (zo(i,k) + zo(i,k+1)) + &
                       cp * to(i,k) + hvap * qo(i,k)
            heso(i,k) = .5 * g * (zo(i,k) + zo(i,k+1)) + &
                       cp * to(i,k) + hvap * qeso(i,k)
            uo(i,k)   = .5 * (uo(i,k) + uo(i,k+1))
            vo(i,k)   = .5 * (vo(i,k) + vo(i,k+1))
          endif
        enddo
      enddo
!
!  look for the level of free convection as cloud base
!
      do i=1,im
        flg(i)   = .true.
        kbcon(i) = kmax(i)
      enddo
      do k = 1, km1
        do i=1,im
          if (flg(i).and.k.le.kbmax(i)) then
            if(k.gt.kb(i).and.heo(i,kb(i)).gt.heso(i,k)) then
              kbcon(i) = k
              flg(i)   = .false.
            endif
          endif
        enddo
      enddo
      !++++ Yi-Chi : Nov 25
      do k = 1, km1
      do i=1,im
        if (k .le. kmax(i)-1) then
        heoout(i,k)   = heo(i,k)
        hesoout(i,k)  = heso(i,k)
        endif
        !write(iulog,*) 'Yi-Chi(sas_conv) k,kmax,gz,cpT',k,kmax(i),.5*g*(zo(i,k)+zo(i,k+1)), cp*to(i,k)
        !write(iulog,*) 'Yi-Chi(sas_conv) k,Lq,Lq*',k,hvap * qo(i,k),hvap * qeso(i,k)
      enddo
      enddo
      !----
!
      do i=1,im
        if(kbcon(i).eq.kmax(i)) cnvflg(i) = .false.
!        write(iulog,*) 'Yi-Chi 588:kbcon v.s. kmax,cnvflg',kbcon(i),cnvflg
      enddo
!!
      totflg = .true.
      do i=1,im
        totflg = totflg .and. (.not. cnvflg(i))
      enddo
      if(totflg) return
!!
!
!  determine critical convective inhibition
!  as a function of vertical velocity at cloud base.
!
      do i=1,im
        if(cnvflg(i)) then
          pdot(i)  = 10.* dot(i,kbcon(i))
        endif
      enddo
!      write(iulog,*) 'Yi-Chi(sas_conv):dot',dot
!      write(iulog,*) 'Yi-Chi(sas_conv):kbcon, pdot',kbcon,pdot
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
            tem = (pdot(i) - w4) / (w3 - w4)
          elseif(pdot(i).ge.-w4) then
            tem = - (pdot(i) + w4) / (w4 - w3)
          else
            tem = 0.
          endif
          val1    =             -1.
          tem = max(tem,val1)
          val2    =             1.
          tem = min(tem,val2)
          tem = 1. - tem
          tem1= .5*(cincrmax-cincrmin)
          cincr = cincrmax - tem * tem1
          pbcdif(i) = pfld(i,kb(i)) - pfld(i,kbcon(i))
          ! +++ Yi-Chi
          cinmask(i) = pbcdif(i)
          !write(iulog,*) 'Yi-Chi(sas_conv) slimsk,w3,w4',slimsk,w3,w4
          !write(iulog,*) 'Yi-Chi(sas_conv) pbcdif',cinmask(i)
          ! ---
          if(pbcdif(i).gt.cincr) then
             cnvflg(i) = .false.
             ! +++ Yi-Chi
             flgmask(i) = 1.
             !cinmask(i) = pbcdif(i)
             ! --- Yi-Chi
          endif
        endif
!         write(iulog,*) 'Yi-Chi638:pbcdif,kb,kbcon,cincr,cnvflg',pbcdif(i),kb(i),kbcon(i),cincr,cnvflg(i)
      enddo
!       write(iulog,*) 'Yi-Chi638:cincr,cnvflg',cincr,cnvflg
!      write(iulog,*) 'Yi-Chi638:kbcon,kb,ktcon,cnvflg',kbcon,kb,ktcon,cnvflg
!      write(iulog,*) 'Yi-Chi(LFC-Base):cincr,slimsk,pdot,pbcdif',cincr,slimsk,pdot,pbcdif
!!
      totflg = .true.
      do i=1,im
        totflg = totflg .and. (.not. cnvflg(i))
      enddo
      if(totflg) return

!!
!
!  assume that updraft entrainment rate above cloud base is
!    same as that at cloud base
!
      do k = 2, km1
        do i=1,im
          if(cnvflg(i).and.  &
           (k.gt.kbcon(i).and.k.lt.kmax(i))) then
              xlamue(i,k) = xlamue(i,kbcon(i))
          endif
        enddo
      enddo
!
!  assume the detrainment rate for the updrafts to be same as
!  the entrainment rate at cloud base
!
      do i = 1, im
        if(cnvflg(i)) then
          xlamud(i) = xlamue(i,kbcon(i))
        endif
      enddo
!
!  functions rapidly decreasing with height, mimicking a cloud ensemble
!    (Bechtold et al., 2008)
!
      do k = 2, km1
        do i=1,im
          if(cnvflg(i).and.  &
           (k.gt.kbcon(i).and.k.lt.kmax(i))) then
              tem = qeso(i,k)/qeso(i,kbcon(i))
              fent1(i,k) = tem**2
              fent2(i,k) = tem**3
          endif
        enddo
      enddo
!
!  final entrainment rate as the sum of turbulent part and organized entrainment
!    depending on the environmental relative humidity
!    (Bechtold et al., 2008)
!
      do k = 2, km1
        do i=1,im
          if(cnvflg(i).and.  &
           (k.ge.kbcon(i).and.k.lt.kmax(i))) then
              tem = cxlamu * frh(i,k) * fent2(i,k)
              xlamue(i,k) = xlamue(i,k)*fent1(i,k) + tem
          endif
        enddo
      enddo
!
!  determine updraft mass flux for the subcloud layers
!
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
!
!  compute mass flux above cloud base
!
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
!
!  compute updraft cloud properties
!
      do i = 1, im
        if(cnvflg(i)) then
          indx         = kb(i)
          hcko(i,indx) = heo(i,indx)
          ucko(i,indx) = uo(i,indx)
          vcko(i,indx) = vo(i,indx)
          pwavo(i)     = 0.
        endif
      enddo
!
!  cloud property is modified by the entrainment process
!
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
              hcko(i,k) = ((1.-tem1)*hcko(i,k-1)+tem*0.5*  &
                          (heo(i,k)+heo(i,k-1)))/factor
              ucko(i,k) = ((1.-tem1)*ucko(i,k-1)+ptem*uo(i,k) &
                          +ptem1*uo(i,k-1))/factor
              vcko(i,k) = ((1.-tem1)*vcko(i,k-1)+ptem*vo(i,k) &
                          +ptem1*vo(i,k-1))/factor
              dbyo(i,k) = hcko(i,k) - heso(i,k)
            endif
          endif
        enddo
      enddo
!
!   taking account into convection inhibition due to existence of
!    dry layers below cloud base
!
      do i=1,im
        flg(i) = cnvflg(i)
        kbcon1(i) = kmax(i)
      enddo
      do k = 2, km1
      do i=1,im
        if (flg(i).and.k.lt.kmax(i)) then
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
             ! +++ Yi-Chi
             flgmask(i) = 2. ! for dry PBL causing large CIN
             ! --- Yi-Chi
          endif
        endif
!        write(iulog,*) 'Yi-Chi808:pfld(kbcon),pfld(kbcon1)',pfld(i,kbcon(i)),pfld(i,kbcon1(i))
!        write(iulog,*) 'Yi-Chi808:dthk=25,cnvflg',tem,cnvflg
      enddo
!!
      totflg = .true.
      do i = 1, im
        totflg = totflg .and. (.not. cnvflg(i))
      enddo
      if(totflg) return
!!
!
!  determine first guess cloud top as the level of zero buoyancy
!
      do i = 1, im
        flg(i) = cnvflg(i)
        ktcon(i) = 1
      enddo
      do k = 2, km1
      do i = 1, im
        if (flg(i).and.k .lt. kmax(i)) then
          if(k.gt.kbcon1(i).and.dbyo(i,k).lt.0.) then
             ktcon(i) = k
             flg(i)   = .false.
          endif
        endif
      enddo
      enddo
!  Yi-Chi: first guess cloud top should > cloud base for cthk (150hPa)
      do i = 1, im
        if(cnvflg(i)) then
          tem = pfld(i,kbcon(i))-pfld(i,ktcon(i))
          if(tem.lt.cthk) cnvflg(i) = .false.
             ! +++ Yi-Chi
             if(.not.cnvflg(i)) then
             flgmask(i) = 3. ! cloud depth should larger than 150
             ! --- Yi-Chi
             endif
!          write(iulog,*) 'Yi-Chi839:pfld(kbcon),pfld(ktcon),cnvflg',kbcon(i),ktcon(i),cnvflg
        endif
      enddo
!!
      totflg = .true.
      do i = 1, im
        totflg = totflg .and. (.not. cnvflg(i))
      enddo
      if(totflg) return
!!
!
!  search for downdraft originating level above theta-e minimum
!
      do i = 1, im
        if(cnvflg(i)) then
           hmin(i) = heo(i,kbcon1(i))
           lmin(i) = kbmax(i)
           jmin(i) = kbmax(i)
        endif
      enddo
      do k = 2, km1
        do i = 1, im
          if (cnvflg(i) .and. k .le. kbmax(i)) then
            if(k.gt.kbcon1(i).and.heo(i,k).lt.hmin(i)) then
               lmin(i) = k + 1
               hmin(i) = heo(i,k)
            endif
          endif
        enddo
      enddo
!
!  make sure that jmin(i) is within the cloud
!
      do i = 1, im
        if(cnvflg(i)) then
          jmin(i) = min(lmin(i),ktcon(i)-1)
          jmin(i) = max(jmin(i),kbcon1(i)+1)
          if(jmin(i).ge.ktcon(i)) cnvflg(i) = .false.
        endif
      enddo
!
!  specify upper limit of mass flux at cloud base
!
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
!
!  compute cloud moisture property and precipitation
!
      do i = 1, im
        if (cnvflg(i)) then
          aa1(i) = 0.
          qcko(i,kb(i)) = qo(i,kb(i))
          ! +++ Yi-Chi : June 22, 2013
          qrcko(i,kb(i)) = qo(i,kb(i))
          ! ---
!         rhbar(i) = 0.
        endif
      enddo
      do k = 2, km1
        do i = 1, im
          if (cnvflg(i)) then
            if(k.gt.kb(i).and.k.lt.ktcon(i)) then
              dz    = zi(i,k) - zi(i,k-1)
              gamma = el2orc * qeso(i,k) / (to(i,k)**2)
              ! +++ Yi-Chi : August 2013
              !qrch = qeso(i,k)  &
              !    + gamma * dbyo(i,k) / (hvap * (1. + gamma))
              qrch = qrch2nd(gamma, dbyo(i,k),to(i,k),qeso(i,k))
              ! ---
!j
              tem  = 0.5 * (xlamue(i,k)+xlamue(i,k-1)) * dz
              tem1 = 0.5 * xlamud(i) * dz
              factor = 1. + tem - tem1
              qcko(i,k) = ((1.-tem1)*qcko(i,k-1)+tem*0.5* &
                          (qo(i,k)+qo(i,k-1)))/factor
          ! +++ Yi-Chi : June 22, 2013
          qrcko(i,k) = qcko(i,k)
          ! ---
!j
              dq = eta(i,k) * (qcko(i,k) - qrch)
!
!             rhbar(i) = rhbar(i) + qo(i,k) / qeso(i,k)
!
!  check if there is excess moisture to release latent heat
!
              if(k.ge.kbcon(i).and.dq.gt.0.) then
                etah = .5 * (eta(i,k) + eta(i,k-1))
                if(ncloud.gt.0..and.k.gt.jmin(i)) then
                  dp = 1000. * del(i,k)
                  qlk = dq / (eta(i,k) + etah * (c0 + c1) * dz)
                  dellal(i,k) = etah * c1 * dz * qlk * g / dp
                else
                  qlk = dq / (eta(i,k) + etah * c0 * dz)
                endif
                ! +++ Yi-Chi
                !qlic(i,k)   = qlk*eta(i,k)*g/dp ! Yi-Chi
                qlic_tmp(i,k)   = qlk*eta(i,k)*g/dp ! Yi-Chi
                ! ---
                aa1(i) = aa1(i) - dz * g * qlk
                qcko(i,k) = qlk + qrch
                pwo(i,k) = etah * c0 * dz * qlk
                pwavo(i) = pwavo(i) + pwo(i,k)
              endif
            endif
          endif
        enddo
      enddo
!
!     do i = 1, im
!       if(cnvflg(i)) then
!         indx = ktcon(i) - kb(i) - 1
!         rhbar(i) = rhbar(i) / float(indx)
!       endif
!     enddo
!
!  calculate cloud work function
!
      do k = 2, km1
        do i = 1, im
          if (cnvflg(i)) then
            if(k.ge.kbcon(i).and.k.lt.ktcon(i)) then
              dz1 = zo(i,k+1) - zo(i,k)
              gamma = el2orc * qeso(i,k) / (to(i,k)**2)
              ! +++ Yi-Chi : Aug 19, 2013
              qrch = qrch2nd(gamma, dbyo(i,k),to(i,k),qeso(i,k))
              dqaa1= qrch - qeso(i,k)        ! q_i-qbar
              dsaa1= dbyo(i,k) - hvap*dqaa1  ! s_i-sbar
              aa1(i) = aa1(i) + &
                     dz1 * (g / (cp * to(i,k))) &
                     * (dsaa1 + cp *to(i,k) * delta * dqaa1)
!              write(iulog,*) '(Yi-Chi:sas_conv.F90:1112)aa1,dqaa1,dsaa1',aa1(i),dqaa1,dsaa1
              !rfact =  1. + delta * cp * gamma &
              !        * to(i,k) / hvap
              !aa1(i) = aa1(i) +   &
              !        dz1 * (g / (cp * to(i,k))) &
              !        * dbyo(i,k) / (1. + gamma) &
              !        * rfact
              ! -----
              val = 0.
              aa1(i)= aa1(i)+     &
                     dz1 * g * delta * &
                     max(val,(qeso(i,k) - qo(i,k))) 
            end if
          end if
        enddo
      enddo
      do i = 1, im
        if(cnvflg(i).and.aa1(i).le.0.) cnvflg(i) = .false.
!        write(iulog,*) 'Yi-Chi975:aa1,cnvflg',aa1(i),cnvflg
      enddo
!!
      totflg = .true.
      do i=1,im
        totflg = totflg .and. (.not. cnvflg(i))
      enddo
      if(totflg) return
!       write(iulog,*) 'Yi-Chi962:kbcon,ktcon,cnvflg',kbcon,ktcon,cnvflg
!!
!
!  estimate the onvective overshooting as the level 
!    where the [aafac * cloud work function] becomes zero,
!    which is the final cloud top
!
      do i = 1, im
        if (cnvflg(i)) then
          aa2(i) = aafac * aa1(i)
        endif
      enddo
!
      do i = 1, im
        flg(i) = cnvflg(i)
        ktcon1(i) = kmax(i) - 1
      enddo
      do k = 2, km1
        do i = 1, im
          if (flg(i)) then
            if(k.ge.ktcon(i).and.k.lt.kmax(i)) then
              dz1 = zo(i,k+1) - zo(i,k)
              gamma = el2orc * qeso(i,k) / (to(i,k)**2)
              !++++++++++
              ! Yi-Chi
              ! Augus 19, 2013: recalculate aa2 with 2nd-order approx
              !----------
              qrch = qrch2nd(gamma, dbyo(i,k),to(i,k),qeso(i,k))	
              dqaa1= qrch - qeso(i,k)        ! q_i-qbar
              dsaa1= dbyo(i,k) - hvap*dqaa1  ! s_i-sbar
              tem =  dz1 * (g / (cp * to(i,k))) &
                     * (dsaa1 + cp *to(i,k) * delta * dqaa1)
              !rfact =  1. + delta * cp * gamma &
              !        * to(i,k) / hvap
              !tem   = dz1 * (g / (cp * to(i,k))) &
              !        * dbyo(i,k) / (1. + gamma) &
              !        * rfact
              !++++++++++
              ! Yi-Chi
              ! March 2013: add a condition for bouyancy < 0
              !----------
              aa2(i) = aa2(i) + tem
              !----
              !aa2(i) = aa2(i) + &
              !        dz1 * (g / (cp * to(i,k))) &
              !        * dbyo(i,k) / (1. + gamma) &
              !        * rfact
              ! if(a2(i).lt.0.) then
              if((aa2(i).lt.0.).or.(tem.gt.0.)) then
              !------------
                ktcon1(i) = k
                flg(i) = .false.
              endif
            endif
          endif
        enddo
      enddo
!
!  compute cloud moisture property, detraining cloud water 
!    and precipitation in overshooting layers 
!
      do k = 2, km1
        do i = 1, im
          if (cnvflg(i)) then
            if(k.ge.ktcon(i).and.k.lt.ktcon1(i)) then
              dz    = zi(i,k) - zi(i,k-1)
              gamma = el2orc * qeso(i,k) / (to(i,k)**2)
              ! +++ Yi-Chi : August 2013
              qrch = qrch2nd(gamma, dbyo(i,k),to(i,k),qeso(i,k))
              !qrch = qeso(i,k)  &
              !    + gamma * dbyo(i,k) / (hvap * (1. + gamma))
              ! ---
!j
              tem  = 0.5 * (xlamue(i,k)+xlamue(i,k-1)) * dz
              tem1 = 0.5 * xlamud(i) * dz
              factor = 1. + tem - tem1
              qcko(i,k) = ((1.-tem1)*qcko(i,k-1)+tem*0.5* &
                          (qo(i,k)+qo(i,k-1)))/factor
              ! +++ Yi-Chi : June 22, 2013
              qrcko(i,k) = qcko(i,k)
              ! ---
!j
              dq = eta(i,k) * (qcko(i,k) - qrch)
!
!  check if there is excess moisture to release latent heat
!
              if(dq.gt.0.) then
                etah = .5 * (eta(i,k) + eta(i,k-1))
                if(ncloud.gt.0.) then
                  dp = 1000. * del(i,k)
                  qlk = dq / (eta(i,k) + etah * (c0 + c1) * dz)
                  dellal(i,k) = etah * c1 * dz * qlk * g / dp
                else
                  qlk = dq / (eta(i,k) + etah * c0 * dz)
                endif
                ! +++ Yi-Chi
                qlic_tmp(i,k)   = qlk*eta(i,k)*g/dp ! Yi-Chi
                !qlic(i,k)   = qlk*eta(i,k)*g/dp ! Yi-Chi
                ! ---
                qcko(i,k) = qlk + qrch
                pwo(i,k) = etah * c0 * dz * qlk
                pwavo(i) = pwavo(i) + pwo(i,k)
              endif
            endif
          endif
        enddo
      enddo
!
! exchange ktcon with ktcon1
!
      do i = 1, im
        if(cnvflg(i)) then
          kk = ktcon(i)
          ktcon(i) = ktcon1(i)
          ktcon1(i) = kk
        endif
      enddo
!       write(iulog,*) 'Yi-Chi1049:kbcon,ktcon,cnvflg',kbcon,ktcon,cnvflg
!
!  this section is ready for cloud water
!
      if(ncloud.gt.0) then
!
!  compute liquid and vapor separation at cloud top
!
      do i = 1, im
        if(cnvflg(i)) then
          k = ktcon(i) - 1
          gamma = el2orc * qeso(i,k) / (to(i,k)**2)
          ! +++ Yi-Chi (Jan 3, 2013):
          qrch = qrch2nd(gamma, dbyo(i,k),to(i,k),qeso(i,k))
          !qrch = qeso(i,k)  &
          !    + gamma * dbyo(i,k) / (hvap * (1. + gamma))
          ! ---
          dq = qcko(i,k) - qrch
!
!  check if there is excess moisture to release latent heat
!
          if(dq.gt.0.) then
            qlko_ktcon(i) = dq
            qcko(i,k) = qrch
          endif
        endif
      enddo
      endif
!
!cccc if(lat.eq.latd.and.lon.eq.lond.and.cnvflg(i)) then
!cccc   print *, ' aa1(i) before dwndrft =', aa1(i)
!cccc endif
!
!------- downdraft calculations
!
!--- compute precipitation efficiency in terms of windshear
!
      do i = 1, im
        if(cnvflg(i)) then
          vshear(i) = 0.
        endif
      enddo
      do k = 2, km
        do i = 1, im
          if (cnvflg(i)) then
            if(k.gt.kb(i).and.k.le.ktcon(i)) then
              shear= sqrt((uo(i,k)-uo(i,k-1)) ** 2  &
                       + (vo(i,k)-vo(i,k-1)) ** 2)
              vshear(i) = vshear(i) + shear
            endif
          endif
        enddo
      enddo
      do i = 1, im
        if(cnvflg(i)) then
          vshear(i) = 1.e3 * vshear(i) / (zi(i,ktcon(i))-zi(i,kb(i)))
          e1=1.591-.639*vshear(i)  &
            +.0953*(vshear(i)**2)-.00496*(vshear(i)**3)
          edt(i)=1.-e1
          val =         .9
          edt(i) = min(edt(i),val)
          val =         .0
          edt(i) = max(edt(i),val)
          edto(i)=edt(i)
          edtx(i)=edt(i)
        endif
      enddo
!
!  determine detrainment rate between 1 and kbcon
!
      do i = 1, im
        if(cnvflg(i)) then
          sumx(i) = 0.
        endif
      enddo
      do k = 1, km1
      do i = 1, im
        if(cnvflg(i).and.k.ge.1.and.k.lt.kbcon(i)) then
          dz = zi(i,k+1) - zi(i,k)
          sumx(i) = sumx(i) + dz
        endif
      enddo
      enddo
      do i = 1, im
        beta = betas
        if(slimsk(i).eq.1.) beta = betal
        if(cnvflg(i)) then
          dz  = (sumx(i)+zi(i,1))/float(kbcon(i))
          tem = 1./float(kbcon(i))
          xlamd(i) = (1.-beta**tem)/dz
        endif
      enddo
!
!  determine downdraft mass flux
!
      do k = km1, 1, -1
        do i = 1, im
          if (cnvflg(i) .and. k .le. kmax(i)-1) then
           if(k.lt.jmin(i).and.k.ge.kbcon(i)) then
              dz        = zi(i,k+1) - zi(i,k)
              ptem      = xlamdd - xlamde
              etad(i,k) = etad(i,k+1) * (1. - ptem * dz)
           else if(k.lt.kbcon(i)) then
              dz        = zi(i,k+1) - zi(i,k)
              ptem      = xlamd(i) + xlamdd - xlamde
              etad(i,k) = etad(i,k+1) * (1. - ptem * dz)
           endif
          endif
        enddo
      enddo
!
!--- downdraft moisture properties
!
      do i = 1, im
        if(cnvflg(i)) then
          jmn = jmin(i)
          hcdo(i,jmn) = heo(i,jmn)
          qcdo(i,jmn) = qo(i,jmn)
          ! +++ Yi-Chi : June 22, 2013
          qrcdo(i,jmn)= qo(i,jmn)
          !qrcdo(i,jmn)= qeso(i,jmn)
          ! ---
          ucdo(i,jmn) = uo(i,jmn)
          vcdo(i,jmn) = vo(i,jmn)
          pwevo(i) = 0.
        endif
      enddo
!j
      do k = km1, 1, -1
        do i = 1, im
          if (cnvflg(i) .and. k.lt.jmin(i)) then
              dz = zi(i,k+1) - zi(i,k)
              if(k.ge.kbcon(i)) then
                 tem  = xlamde * dz
                 tem1 = 0.5 * xlamdd * dz
              else
                 tem  = xlamde * dz
                 tem1 = 0.5 * (xlamd(i)+xlamdd) * dz
              endif
              factor = 1. + tem - tem1
              ptem = 0.5 * tem - pgcon
              ptem1= 0.5 * tem + pgcon
              hcdo(i,k) = ((1.-tem1)*hcdo(i,k+1)+tem*0.5* &
                          (heo(i,k)+heo(i,k+1)))/factor
              ucdo(i,k) = ((1.-tem1)*ucdo(i,k+1)+ptem*uo(i,k+1) &
                          +ptem1*uo(i,k))/factor
              vcdo(i,k) = ((1.-tem1)*vcdo(i,k+1)+ptem*vo(i,k+1) &
                          +ptem1*vo(i,k))/factor
              dbyo(i,k) = hcdo(i,k) - heso(i,k)
          endif
        enddo
      enddo
!
      do k = km1, 1, -1
        do i = 1, im
          if (cnvflg(i).and.k.lt.jmin(i)) then
              gamma      = el2orc * qeso(i,k) / (to(i,k)**2)
              ! +++ Yi-Chi (Aug 2013)
              qrcdo(i,k) = qrch2nd(gamma,dbyo(i,k),to(i,k),qeso(i,k))
              !qrcdo(i,k) = qeso(i,k)+ &
              !       (1./hvap)*(gamma/(1.+gamma))*dbyo(i,k)
              ! ---
!             detad      = etad(i,k+1) - etad(i,k)
!j
              dz = zi(i,k+1) - zi(i,k)
              if(k.ge.kbcon(i)) then
                 tem  = xlamde * dz
                 tem1 = 0.5 * xlamdd * dz
              else
                 tem  = xlamde * dz
                 tem1 = 0.5 * (xlamd(i)+xlamdd) * dz
              endif
              factor = 1. + tem - tem1
              ! +++ Yi-Chi : June 22, 2013
              qcdo(i,k) = ((1.-tem1)*qrcdo(i,k+1)+tem*0.5* &
                          (qo(i,k)+qo(i,k+1)))/factor
              !qcdo(i,k) = ((1.-tem1)*qcdo(i,k+1)+tem*0.5* &
              !            (qo(i,k)+qo(i,k+1)))/factor
              ! ---
!j
!             pwdo(i,k)  = etad(i,k+1) * qcdo(i,k+1) -
!    &                     etad(i,k) * qrcdo(i,k)
!             pwdo(i,k)  = pwdo(i,k) - detad *
!    &                    .5 * (qrcdo(i,k) + qrcdo(i,k+1))
!j
! +++ Yi-Chi : June 22, 2013
               pwdo(i,k)  = etad(i,k) * (qcdo(i,k) - qrcdo(i,k))
!              pwdo(i,k)  = etad(i,k+1) * (qcdo(i,k) - qrcdo(i,k))
!              qcdo(i,k)  = qrcdo(i,k)
! ---
              pwevo(i)   = pwevo(i) + pwdo(i,k)
          endif
        enddo
      enddo
!
!--- final downdraft strength dependent on precip
!--- efficiency (edt), normalized condensate (pwav), and
!--- evaporate (pwev)
!
      do i = 1, im
        edtmax = edtmaxl
        if(slimsk(i).eq.0.) edtmax = edtmaxs
        if(cnvflg(i)) then
          if(pwevo(i).lt.0.) then
            edto(i) = -edto(i) * pwavo(i) / pwevo(i)
            edto(i) = min(edto(i),edtmax)
          else
            edto(i) = 0.
          endif
        endif
      enddo
!
!--- downdraft cloudwork functions
!
      do k = km1, 1, -1
        do i = 1, im
          if (cnvflg(i) .and. k .lt. jmin(i)) then
              gamma = el2orc * qeso(i,k) / to(i,k)**2
              dhh=hcdo(i,k)
              dt=to(i,k)
              dg=gamma
              dh=heso(i,k)
              dz=-1.*(zo(i,k+1)-zo(i,k))
              ! +++ Yi-Chi : Aug 19, 2013
              ! qrch = qrch2nd(gamma, dbyo(i,k),to(i,k),qeso(i,k))
              qrch = qrch2nd(dg,dhh-dh,dt,qeso(i,k))
              dqaa1= qrch - qeso(i,k)        ! q_i-qbar
              dsaa1= dhh  - dh - hvap*dqaa1  ! s_i-sbar
              aa1(i) = aa1(i) + &
                     edto(i) * dz * (g / (cp * dt)) &
                     * (dsaa1 + cp *dt * delta * dqaa1)
!              write(iulog,*) '(Yi-Chi:sas_conv.F90:1482)downdraft aa1,dqaa1,dsaa1',aa1(i),dqaa1,dsaa1
              !aa1(i)=aa1(i)+edto(i)*dz*(g/(cp*dt))*((dhh-dh)/(1.+dg)) &
              !      *(1.+delta*cp*dg*dt/hvap)
              ! --- Yi-Chi
              val=0.
              aa1(i)=aa1(i)+edto(i)*  &
             dz*g*delta*max(val,(qeso(i,k)-qo(i,k)))
          endif
        enddo
      enddo
      do i = 1, im
        if(cnvflg(i).and.aa1(i).le.0.) then
           cnvflg(i) = .false.
        endif
!        write(iulog,*) '(Yi-Chi:sas_conv.F90:1482)downdraft aa1<0,cnvflg',cnvflg(i)
      enddo
!!
      totflg = .true.
      do i=1,im
        totflg = totflg .and. (.not. cnvflg(i))
      enddo
      if(totflg) return
!!
!
!--- what would the change be, that a cloud with unit mass
!--- will do to the environment?
!
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
      do i = 1, im
        if(cnvflg(i)) then
          dp = 1000. * del(i,1)
          dellah(i,1) = edto(i) * etad(i,1) * (hcdo(i,1) &
                        - heo(i,1)) * g / dp
          ! +++ Yi-Chi
          dellaq(i,1) = edto(i) * etad(i,1) * (qrcdo(i,1) &
          !dellaq(i,1) = edto(i) * etad(i,1) * (qcdo(i,1) &
          ! ----
                        - qo(i,1)) * g / dp
          dellau(i,1) = edto(i) * etad(i,1) * (ucdo(i,1) &
                        - uo(i,1)) * g / dp
          dellav(i,1) = edto(i) * etad(i,1) * (vcdo(i,1) &
                        - vo(i,1)) * g / dp
          ! Yi-Chi:
          ddout(i,1)  = edto(i) * etad(i,1) * g / dp
          ! Yi-Chi : Apr26 : use fraction of cloud water
          !ddout(i,1)  = edto(i) * etad(i,1) * (qcdo(i,1) &
          !              - qo(i,1))/qcdo(i,1) * g / dp
        endif
      enddo
!
!--- changed due to subsidence and entrainment
!
      do k = 2, km1
        do i = 1, im
          if (cnvflg(i).and.k.lt.ktcon(i)) then
              aup = 1.
              if(k.le.kb(i)) aup = 0.
              adw = 1.
              if(k.gt.jmin(i)) adw = 0.
              dp = 1000. * del(i,k)
              dz = zi(i,k) - zi(i,k-1)
!
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
!
              tem  = 0.5 * (xlamue(i,k)+xlamue(i,k-1))
              tem1 = xlamud(i)
!
              if(k.le.kbcon(i)) then
                ptem  = xlamde
                ptem1 = xlamd(i)+xlamdd
              else
                ptem  = xlamde
                ptem1 = xlamdd
              endif

              ! Yi-Chi (Apr 1, 2012) : 
              !  output entrainment/detrainment rate for convective transport
              !  unit : fraction of unit mass flux (1/m)
              !    transform the entrainment/detrainment rate to g/Pa
              !    dp(Pa)
              eutmp(i,k) = aup*tem*eta(i,k-1)*dz*(g/dp)   ! updraft entrainment rate per unit length (m?)
              dutmp(i,k) = aup*tem1*eta(i,k-1)*dz*(g/dp)  ! updraft detrainment
              edtmp(i,k) = adw*edto(i)*ptem*etad(i,k)*dz*(g/dp)   ! downdraft entrainment
              ddtmp(i,k) = adw*edto(i)*ptem1*etad(i,k)*dz*(g/dp)  ! downdraft detrainment
              duqtmp(i,k)=  aup*tem1*eta(i,k-1)*.5*(qcko(i,k)+qcko(i,k-1))*dz*g/dp 

!j
              dellah(i,k) = dellah(i,k) +  &
          ((aup*eta(i,k)-adw*edto(i)*etad(i,k))*dv1h &
         - (aup*eta(i,k-1)-adw*edto(i)*etad(i,k-1))*dv3h  &
         - (aup*tem*eta(i,k-1)+adw*edto(i)*ptem*etad(i,k))*dv2h*dz &
         +  aup*tem1*eta(i,k-1)*.5*(hcko(i,k)+hcko(i,k-1))*dz  &
         +  adw*edto(i)*ptem1*etad(i,k)*.5*(hcdo(i,k)+hcdo(i,k-1))*dz &
              ) *g/dp
!j
              dellaq(i,k) = dellaq(i,k) +  &
          ((aup*eta(i,k)-adw*edto(i)*etad(i,k))*dv1q  &
         - (aup*eta(i,k-1)-adw*edto(i)*etad(i,k-1))*dv3q   &
         - (aup*tem*eta(i,k-1)+adw*edto(i)*ptem*etad(i,k))*dv2q*dz   &
! +++ Yi-Chi : June 22, 2013
         +  aup*tem1*eta(i,k-1)*.5*(qrcko(i,k)+qcko(i,k-1))*dz &
         +  adw*edto(i)*ptem1*etad(i,k)*.5*(qrcdo(i,k)+qcdo(i,k-1))*dz &
!         +  aup*tem1*eta(i,k-1)*.5*(qcko(i,k)+qcko(i,k-1))*dz   &
!         +  adw*edto(i)*ptem1*etad(i,k)*.5*(qrcdo(i,k)+qrcdo(i,k-1))*dz   &
! ---
              ) *g/dp
!j
              dellau(i,k) = dellau(i,k) +  &
          ((aup*eta(i,k)-adw*edto(i)*etad(i,k))*dv1u   &
         - (aup*eta(i,k-1)-adw*edto(i)*etad(i,k-1))*dv3u   &
         - (aup*tem*eta(i,k-1)+adw*edto(i)*ptem*etad(i,k))*dv2u*dz   &
         +  aup*tem1*eta(i,k-1)*.5*(ucko(i,k)+ucko(i,k-1))*dz   &
         +  adw*edto(i)*ptem1*etad(i,k)*.5*(ucdo(i,k)+ucdo(i,k-1))*dz   &
         -  pgcon*(aup*eta(i,k-1)-adw*edto(i)*etad(i,k))*(dv1u-dv3u)   &
              ) *g/dp
!j
              dellav(i,k) = dellav(i,k) +  &
          ((aup*eta(i,k)-adw*edto(i)*etad(i,k))*dv1v   &
         - (aup*eta(i,k-1)-adw*edto(i)*etad(i,k-1))*dv3v   &
         - (aup*tem*eta(i,k-1)+adw*edto(i)*ptem*etad(i,k))*dv2v*dz   &
         +  aup*tem1*eta(i,k-1)*.5*(vcko(i,k)+vcko(i,k-1))*dz   &
         +  adw*edto(i)*ptem1*etad(i,k)*.5*(vcdo(i,k)+vcdo(i,k-1))*dz   &
         -  pgcon*(aup*eta(i,k-1)-adw*edto(i)*etad(i,k))*(dv1v-dv3v)   &
              ) *g/dp
!j
          endif
        enddo
      enddo
!
!------- cloud top
!
      do i = 1, im
        if(cnvflg(i)) then
          indx = ktcon(i)
          dp = 1000. * del(i,indx)
          dv1h = heo(i,indx-1)
          dellah(i,indx) = eta(i,indx-1) *   &
                          (hcko(i,indx-1) - dv1h) * g / dp
          dv1q = qo(i,indx-1)
          dellaq(i,indx) = eta(i,indx-1) *   &
                          (qcko(i,indx-1) - dv1q) * g / dp
          dv1u = uo(i,indx-1)
          dellau(i,indx) = eta(i,indx-1) *   &
                          (ucko(i,indx-1) - dv1u) * g / dp
          dv1v = vo(i,indx-1)
          dellav(i,indx) = eta(i,indx-1) *   &
                          (vcko(i,indx-1) - dv1v) * g / dp
!
!  cloud water
!
          dellal(i,indx) = eta(i,indx-1) *   &
                          qlko_ktcon(i) * g / dp

          ! Yi-Chi (Apr 1, 2012) :
          !  output entrainment/detrainment rate for convective transport
          !  At the cloud top, all the mass detrained in the updraft
          !  downdraft won't reach to this level : (g/Pa)
          ! (Apr 26, 2012)
          dutmp(i,indx) = eta(i,indx-1)*(g/dp)  ! updraft detrainment
          !duout(i,indx) = eta(i,indx-1)*qlko_ktcon(i)/ &
          !                qcko(i,indx-1)*(g/dp)
          duqtmp(i,indx) = eta(i,indx-1) *   &
                          (qcko(i,indx-1) - dv1q) * g / dp
        endif
      enddo
!
!------- final changed variable per unit mass flux
!
      do k = 1, km
        do i = 1, im
          if (cnvflg(i).and.k .le. kmax(i)) then
            if(k.gt.ktcon(i)) then
              qo(i,k) = q1(i,k)
              to(i,k) = t1(i,k)
            endif
            if(k.le.ktcon(i)) then
              qo(i,k) = dellaq(i,k) * mbdt + q1(i,k)
              dellat = (dellah(i,k) - hvap * dellaq(i,k)) / cp
              to(i,k) = dellat * mbdt + t1(i,k)
              val   =           1.e-10
              qo(i,k) = max(qo(i,k), val  )
            endif
          endif
        enddo
      enddo

! Yi-Chi : Aug 14, 2012 : output updraft properties
     do i = 1, im
        if (cnvflg(i)) then
        do k = 1, kmax(i)
        qckoout(i,k) = qcko(i,k)      
        enddo
        endif
      enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!--- the above changed environment is now used to calulate the
!--- effect the arbitrary cloud (with unit mass flux)
!--- would have on the stability,
!--- which then is used to calculate the real mass flux,
!--- necessary to keep this change in balance with the large-scale
!--- destabilization.
!
!--- environmental conditions again, first heights
!
      do k = 1, km
        do i = 1, im
          if(cnvflg(i) .and. k .le. kmax(i)) then
            ! Yi-Chi (July 2012)-----------
            !esat = c1e*exp(c2e*(to(i,k)-tfreez)/(c3e+to(i,k)-tfreez))       ! esat(T) in mb
            !qeso(i,k) = 0.01*esat*100 ! mb->Pa
            qeso(i,k) = esat(to(i,k),flgfvps)
            !qeso(i,k) = 0.01 * fpvs(to(i,k))      ! fpvs is in pa
            !------------------------------
            qeso(i,k) = eps * qeso(i,k) / (pfld(i,k)+epsm1*qeso(i,k))
            val       =             1.e-8
            qeso(i,k) = max(qeso(i,k), val )
!           tvo(i,k)  = to(i,k) + delta * to(i,k) * qo(i,k)
          endif
        enddo
      enddo
!
!--- moist static energy
!
      do k = 1, km1
        do i = 1, im
          if(cnvflg(i) .and. k .le. kmax(i)-1) then
            dz = .5 * (zo(i,k+1) - zo(i,k))
            dp = .5 * (pfld(i,k+1) - pfld(i,k))
            ! Yi-Chi (March 2013)
            es = esat(to(i,k+1),flgfvps)
            ! Yi-Chi (July 2012)-------------
            !esat = c1e*exp(c2e*(to(i,k+1)-tfreez)/(c3e+to(i,k+1)-tfreez))       ! esat(T) in mb
            !es = 0.01*esat*100 ! mb->Pa
            !es = esat(to(i,k+1))
            !es = 0.01 * fpvs(to(i,k+1))      ! fpvs is in pa
            !--------------------------------
            pprime = pfld(i,k+1) + epsm1 * es
            qs = eps * es / pprime
            dqsdp = - qs / pprime
            desdt = es * (fact1 / to(i,k+1) + fact2 / (to(i,k+1)**2))
            dqsdt = qs * pfld(i,k+1) * desdt / (es * pprime)
            gamma = el2orc * qeso(i,k+1) / (to(i,k+1)**2)
            dt = (g * dz + hvap * dqsdp * dp) / (cp * (1. + gamma))
            dq = dqsdt * dt + dqsdp * dp
            to(i,k) = to(i,k+1) + dt
            qo(i,k) = qo(i,k+1) + dq
            po(i,k) = .5 * (pfld(i,k) + pfld(i,k+1))
          endif
        enddo
      enddo
      do k = 1, km1
        do i = 1, im
          if(cnvflg(i) .and. k .le. kmax(i)-1) then
            ! Yi-Chi (March 2013)
            qeso(i,k) = esat(to(i,k),flgfvps)
            ! Yi-Chi (July 2012)-------------
            !esat = c1e*exp(c2e*(to(i,k)-tfreez)/(c3e+to(i,k)-tfreez))       ! esat(T) in mb
            !qeso(i,k) = 0.01*esat*100 ! mb->Pa
            !qeso(i,k) = 0.01 * fpvs(to(i,k))      ! fpvs is in pa
            !--------------------------------
            qeso(i,k) = eps * qeso(i,k) / (po(i,k) + epsm1 * qeso(i,k))
            val1      =             1.e-8
            qeso(i,k) = max(qeso(i,k), val1)
            val2      =           1.e-10
            qo(i,k)   = max(qo(i,k), val2 )
!           qo(i,k)   = min(qo(i,k),qeso(i,k))
            heo(i,k)   = .5 * g * (zo(i,k) + zo(i,k+1)) +   &
                         cp * to(i,k) + hvap * qo(i,k)
            heso(i,k) = .5 * g * (zo(i,k) + zo(i,k+1)) +   &
                       cp * to(i,k) + hvap * qeso(i,k)
          endif
        enddo
      enddo
!      write(iulog,*) 'Yi-Chi(sas_conv) hvap,qo: ',hvap,qo
!      write(iulog,*) 'Yi-Chi(sas_conv) cp,to: ',cp,to
!      write(iulog,*) 'Yi-Chi(sas_conv) g,zo: ',g,zo
      do i = 1, im
        if(cnvflg(i)) then
          k = kmax(i)
          heo(i,k) = g * zo(i,k) + cp * to(i,k) + hvap * qo(i,k)
          heso(i,k) = g * zo(i,k) + cp * to(i,k) + hvap * qeso(i,k)
!         heo(i,k) = min(heo(i,k),heso(i,k))
        endif
      enddo
!       write(iulog,*) 'Yi-Chi(sas_conv) heo,heso: ',heo,heso
!
!**************************** static control
!
!------- moisture and cloud work functions
!
      do i = 1, im
        if(cnvflg(i)) then
          xaa0(i) = 0.
          xpwav(i) = 0.
        endif
      enddo
!
      do i = 1, im
        if(cnvflg(i)) then
          indx = kb(i)
          hcko(i,indx) = heo(i,indx)
          qcko(i,indx) = qo(i,indx)
        endif
      enddo
      do k = 2, km1
        do i = 1, im
          if (cnvflg(i)) then
            if(k.gt.kb(i).and.k.le.ktcon(i)) then
              dz = zi(i,k) - zi(i,k-1)
              tem  = 0.5 * (xlamue(i,k)+xlamue(i,k-1)) * dz
              tem1 = 0.5 * xlamud(i) * dz
              factor = 1. + tem - tem1
              hcko(i,k) = ((1.-tem1)*hcko(i,k-1)+tem*0.5*   &
                          (heo(i,k)+heo(i,k-1)))/factor
            endif
          endif
        enddo
      enddo
      do k = 2, km1
        do i = 1, im
          if (cnvflg(i)) then
            if(k.gt.kb(i).and.k.lt.ktcon(i)) then
              dz = zi(i,k) - zi(i,k-1)
              gamma = el2orc * qeso(i,k) / (to(i,k)**2)
              xdby = hcko(i,k) - heso(i,k)
              ! +++ Yi-Chi : Aug 2013
              xqrch = qrch2nd(gamma, xdby,to(i,k),qeso(i,k))
              !xqrch = qeso(i,k)    &
              !     + gamma * xdby / (hvap * (1. + gamma))
              ! ---
!j
              tem  = 0.5 * (xlamue(i,k)+xlamue(i,k-1)) * dz
              tem1 = 0.5 * xlamud(i) * dz
              factor = 1. + tem - tem1
              qcko(i,k) = ((1.-tem1)*qcko(i,k-1)+tem*0.5*    &
                          (qo(i,k)+qo(i,k-1)))/factor
!j
              dq = eta(i,k) * (qcko(i,k) - xqrch)
!
              if(k.ge.kbcon(i).and.dq.gt.0.) then
                etah = .5 * (eta(i,k) + eta(i,k-1))
                if(ncloud.gt.0..and.k.gt.jmin(i)) then
                  qlk = dq / (eta(i,k) + etah * (c0 + c1) * dz)
                else
                  qlk = dq / (eta(i,k) + etah * c0 * dz)
                endif
                if(k.lt.ktcon1(i)) then
                  xaa0(i) = xaa0(i) - dz * g * qlk
                endif
                qcko(i,k) = qlk + xqrch
                xpw = etah * c0 * dz * qlk
                xpwav(i) = xpwav(i) + xpw
              endif
            endif
            if(k.ge.kbcon(i).and.k.lt.ktcon1(i)) then
              dz1 = zo(i,k+1) - zo(i,k)
              gamma = el2orc * qeso(i,k) / (to(i,k)**2)
              ! +++ Yi-Chi :
              qrch = qrch2nd(gamma, xdby,to(i,k),qeso(i,k))
              dqaa1= qrch - qeso(i,k)        ! q_i-qbar
              dsaa1= xdby - hvap*dqaa1  ! s_i-sbar
              xaa0(i) = xaa0(i) + &
                     dz1 * (g / (cp * to(i,k))) &
                     * (dsaa1 + cp *to(i,k) * delta * dqaa1) 
              !rfact =  1. + delta * cp * gamma &
              !        * to(i,k) / hvap
              !xaa0(i) = xaa0(i)       &
              !       + dz1 * (g / (cp * to(i,k)))  &
              !       * xdby / (1. + gamma)  &
              !       * rfact
              !----------------------
              val=0.
              xaa0(i)=xaa0(i)+    &
                      dz1 * g * delta *  &
                      max(val,(qeso(i,k) - qo(i,k)))
            endif
          endif
        enddo
      enddo
!
!------- downdraft calculations
!
!--- downdraft moisture properties
!
      do i = 1, im
        if(cnvflg(i)) then
          jmn = jmin(i)
          hcdo(i,jmn) = heo(i,jmn)
          qcdo(i,jmn) = qo(i,jmn)
          ! +++ Yi-Chi : June 22, 2013
          qrcd(i,jmn) = qo(i,jmn)
          !qrcd(i,jmn) = qeso(i,jmn)
          ! ---
          xpwev(i) = 0.
        endif
      enddo
!j
      do k = km1, 1, -1
        do i = 1, im
          if (cnvflg(i) .and. k.lt.jmin(i)) then
              dz = zi(i,k+1) - zi(i,k)
              if(k.ge.kbcon(i)) then
                 tem  = xlamde * dz
                 tem1 = 0.5 * xlamdd * dz
              else
                 tem  = xlamde * dz
                 tem1 = 0.5 * (xlamd(i)+xlamdd) * dz
              endif
              factor = 1. + tem - tem1
              hcdo(i,k) = ((1.-tem1)*hcdo(i,k+1)+tem*0.5*  &
                          (heo(i,k)+heo(i,k+1)))/factor
          endif
        enddo
      enddo
!j
      do k = km1, 1, -1
        do i = 1, im
          if (cnvflg(i) .and. k .lt. jmin(i)) then
              dq = qeso(i,k)
              dt = to(i,k)
              gamma    = el2orc * dq / dt**2
              dh       = hcdo(i,k) - heso(i,k)
              ! +++ Yi-Chi (Aug 2013)
              qrcd(i,k) = qrch2nd(gamma, dh,dt,dq)
              !qrcd(i,k)=dq+(1./hvap)*(gamma/(1.+gamma))*dh
              ! ---
!             detad    = etad(i,k+1) - etad(i,k)
!j
              dz = zi(i,k+1) - zi(i,k)
              if(k.ge.kbcon(i)) then
                 tem  = xlamde * dz
                 tem1 = 0.5 * xlamdd * dz
              else
                 tem  = xlamde * dz
                 tem1 = 0.5 * (xlamd(i)+xlamdd) * dz
              endif
              factor = 1. + tem - tem1
! +++ Yi-Chi : June 22, 2013
               qcdo(i,k) = ((1.-tem1)*qrcd(i,k+1)+tem*0.5* &
!              qcdo(i,k) = ((1.-tem1)*qcdo(i,k+1)+tem*0.5* &
! ---
                          (qo(i,k)+qo(i,k+1)))/factor
!j
!             xpwd     = etad(i,k+1) * qcdo(i,k+1) -
!    &                   etad(i,k) * qrcd(i,k)
!             xpwd     = xpwd - detad *
!    &                 .5 * (qrcd(i,k) + qrcd(i,k+1))
!j
! +++ Yi-Chi : June 22, 2013
              xpwd     = etad(i,k) * (qcdo(i,k) - qrcd(i,k))
!              xpwd     = etad(i,k+1) * (qcdo(i,k) - qrcd(i,k))
!              qcdo(i,k)= qrcd(i,k)
! ---
              xpwev(i) = xpwev(i) + xpwd
          endif
        enddo
      enddo
!
      do i = 1, im
        edtmax = edtmaxl
        if(slimsk(i).eq.0.) edtmax = edtmaxs
        if(cnvflg(i)) then
          if(xpwev(i).ge.0.) then
            edtx(i) = 0.
          else
            edtx(i) = -edtx(i) * xpwav(i) / xpwev(i)
            edtx(i) = min(edtx(i),edtmax)
          endif
        endif
      enddo
!
!
!--- downdraft cloudwork functions
!
!
      do k = km1, 1, -1
        do i = 1, im
          if (cnvflg(i) .and. k.lt.jmin(i)) then
              gamma = el2orc * qeso(i,k) / to(i,k)**2
              dhh=hcdo(i,k)
              dt= to(i,k)
              dg= gamma
              dh= heso(i,k)
              dz=-1.*(zo(i,k+1)-zo(i,k))
              ! +++ Yi-Chi : august 2013
              qrch = qrch2nd(dg,dhh-dh,dt,qeso(i,k))
              dqaa1= qrch - qeso(i,k)        ! q_i-qbar
              dsaa1= dhh  - dh - hvap*dqaa1  ! s_i-sbar
              xaa0(i) = xaa0(i) + &
                     edtx(i) * dz * (g / (cp * dt)) &
                     * (dsaa1 + cp *dt * delta * dqaa1)
              !xaa0(i)=xaa0(i)+edtx(i)*dz*(g/(cp*dt))*((dhh-dh)/(1.+dg)) &
              !       *(1.+delta*cp*dg*dt/hvap)
              ! ---
              val=0.
              xaa0(i)=xaa0(i)+edtx(i)*   &
             dz*g*delta*max(val,(qeso(i,k)-qo(i,k)))
          endif
        enddo
      enddo
!
!  calculate critical cloud work function
!
      do i = 1, im
        if(cnvflg(i)) then
          if(pfld(i,ktcon(i)).lt.pcrit(15))then
            acrt(i)=acrit(15)*(975.-pfld(i,ktcon(i)))  &
                   /(975.-pcrit(15))
          else if(pfld(i,ktcon(i)).gt.pcrit(1))then
            acrt(i)=acrit(1)
          else
            k =  int((850. - pfld(i,ktcon(i)))/50.) + 2
            k = min(k,15)
            k = max(k,2)
            acrt(i)=acrit(k)+(acrit(k-1)-acrit(k))*  &
                (pfld(i,ktcon(i))-pcrit(k))/(pcrit(k-1)-pcrit(k))
          endif
        endif
      enddo
      do i = 1, im
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
!
!  modify critical cloud workfunction by cloud base vertical velocity
!
          if(pdot(i).le.w4) then
            acrtfct(i) = (pdot(i) - w4) / (w3 - w4)
          elseif(pdot(i).ge.-w4) then
            acrtfct(i) = - (pdot(i) + w4) / (w4 - w3)
          else
            acrtfct(i) = 0.
          endif
          val1    =             -1.
          acrtfct(i) = max(acrtfct(i),val1)
          val2    =             1.
          acrtfct(i) = min(acrtfct(i),val2)
          acrtfct(i) = 1. - acrtfct(i)
!
!  modify acrtfct(i) by colume mean rh if rhbar(i) is greater than 80 percent
!
!         if(rhbar(i).ge..8) then
!           acrtfct(i) = acrtfct(i) * (.9 - min(rhbar(i),.9)) * 10.
!         endif
!
!  modify adjustment time scale by cloud base vertical velocity
!
          dtconv(i) = dt2 + max((1800. - dt2),0.) *  &
                     (pdot(i) - w2) / (w1 - w2)
!         dtconv(i) = max(dtconv(i), dt2)
!         dtconv(i) = 1800. * (pdot(i) - w2) / (w1 - w2)
          dtconv(i) = max(dtconv(i),dtmin)
          dtconv(i) = min(dtconv(i),dtmax)
!
        endif
      enddo
!
!--- large scale forcing
!
      do i= 1, im
        if(cnvflg(i)) then
          fld(i)=(aa1(i)-acrt(i)* acrtfct(i))/dtconv(i)
          if(fld(i).le.0.) cnvflg(i) = .false.
        endif
!       write(iulog,*) '(Yi-Chi:sas_conv.F90:2067)quasi-equilibrium xk>0,fld,cnvflg',fld(i),cnvflg(i)
        if(cnvflg(i)) then
!         xaa0(i) = max(xaa0(i),0.)
          xk(i) = (xaa0(i) - aa1(i)) / mbdt
          if(xk(i).ge.0.) cnvflg(i) = .false.
        endif
!        write(iulog,*) '(Yi-Chi:sas_conv.F90:2067)quasi-equilibrium xk>0,xaa0,cnvflg',xaa0,cnvflg(i)
!
!--- kernel, cloud base mass flux
!
        if(cnvflg(i)) then
          xmb(i) = -fld(i) / xk(i)
          xmb(i) = min(xmb(i),xmbmax(i))
        endif
      enddo
      ! +++ Yi-Chi : output closure-related variables
      do i= 1, im
      if(cnvflg(i)) then
      xmbout(i)   = xmb(i)
      xkout(i)    = xk(i)
      fldout(i)   = fld(i)
      dtconvout(i)= dtconv(i)
      aa1crit(i)  = acrt(i)* acrtfct(i)
      xomega(i)   = pdot(i)
      endif
      enddo
      ! ---
!!
      totflg = .true.
      do i=1,im
        totflg = totflg .and. (.not. cnvflg(i))
      enddo
      if(totflg) return
!!
!
!  restore to,qo,uo,vo to t1,q1,u1,v1 in case convection stops
!
      do k = 1, km
        do i = 1, im
          if (cnvflg(i) .and. k .le. kmax(i)) then
            to(i,k) = t1(i,k)
            qo(i,k) = q1(i,k)
            uo(i,k) = u1(i,k)
            vo(i,k) = v1(i,k)
            ! Yi-Chi (March 2013)
            qeso(i,k) = esat(t1(i,k),flgfvps)
            ! Yi-Chi (July 2012)-------------
            !esat = c1e*exp(c2e*(t1(i,k)-tfreez)/(c3e+t1(i,k)-tfreez))       ! esat(T) in mb
            !qeso(i,k) = 0.01*esat*100 ! mb->Pa
            !qeso(i,k) = 0.01 * fpvs(t1(i,k))      ! fpvs is in pa
            !---------------------------------
            qeso(i,k) = eps * qeso(i,k) / (pfld(i,k) + epsm1*qeso(i,k))
            val     =             1.e-8
            qeso(i,k) = max(qeso(i,k), val )
          endif
        enddo
      enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!--- feedback: simply the changes from the cloud with unit mass flux
!---           multiplied by  the mass flux necessary to keep the
!---           equilibrium with the larger-scale.
!
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
          if (cnvflg(i) .and. k .le. kmax(i)) then
            if(k.le.ktcon(i)) then
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
          if (cnvflg(i) .and. k .le. kmax(i)) then
            if(k.le.ktcon(i)) then
            ! Yi-Chi (March 2013)
            ! replace with a function
            qeso(i,k) = esat(t1(i,k),flgfvps)
            ! Yi-Chi (July 2012)-------------
            !esat = c1e*exp(c2e*(t1(i,k)-tfreez)/(c3e+t1(i,k)-tfreez))       ! esat(T) in mb
            !qeso(i,k) = 0.01*esat*100 ! mb->Pa
            !qeso(i,k) = 0.01 * fpvs(t1(i,k))      ! fpvs is in pa
            !--------------------------------
              qeso(i,k) = eps * qeso(i,k)/(pfld(i,k) + epsm1*qeso(i,k))
              val     =             1.e-8
              qeso(i,k) = max(qeso(i,k), val )
            endif
          endif
        enddo
      enddo
!
      do i = 1, im
        rntot(i) = 0.
        delqev(i) = 0.
        delq2(i) = 0.
        flg(i) = cnvflg(i)
      enddo
      do k = km, 1, -1
        do i = 1, im
          if (cnvflg(i) .and. k .le. kmax(i)) then
            if(k.lt.ktcon(i)) then
              aup = 1.
              if(k.le.kb(i)) aup = 0.
              adw = 1.
              if(k.ge.jmin(i)) adw = 0.
              rain =  aup * pwo(i,k) + adw * edto(i) * pwdo(i,k)
              rntot(i) = rntot(i) + rain * xmb(i) * .001 * dt2
            endif
          endif
        enddo
      enddo
      do k = km, 1, -1
        do i = 1, im
          if (k .le. kmax(i)) then
            deltv(i) = 0.
            delq(i) = 0.
            qevap(i) = 0.
  
            if(cnvflg(i).and.k.lt.ktcon(i)) then
              aup = 1.
              if(k.le.kb(i)) aup = 0.
              adw = 1.
              if(k.ge.jmin(i)) adw = 0.
              rain =  aup * pwo(i,k) + adw * edto(i) * pwdo(i,k)
              rn(i) = rn(i) + rain * xmb(i) * .001 * dt2
              ! Yi-Chi: add rain production at each level Mar 2012
              rainout(i,k) = rain * xmb(i)*g/(1000. * del(i,k))
            endif
            if(flg(i).and.k.lt.ktcon(i)) then
              evef = edt(i) * evfact
              if(slimsk(i).eq.1.) evef=edt(i) * evfactl
!             if(slimsk(i).eq.1.) evef=.07
!             if(slimsk(i).ne.1.) evef = 0.
              qcond(i) = evef * (q1(i,k) - qeso(i,k)) &
                      / (1. + el2orc * qeso(i,k) / t1(i,k)**2)
              dp = 1000. * del(i,k)
              if(rn(i).gt.0..and.qcond(i).lt.0.) then
                qevap(i) = -qcond(i) * (1.-exp(-.32*sqrt(dt2*rn(i))))
                qevap(i) = min(qevap(i), rn(i)*1000.*g/dp)
                delq2(i) = delqev(i) + .001 * qevap(i) * dp / g
              endif
              if(rn(i).gt.0..and.qcond(i).lt.0..and.   &
                delq2(i).gt.rntot(i)) then
                qevap(i) = 1000.* g * (rntot(i) - delqev(i)) / dp
                flg(i) = .false.
              endif
              if(rn(i).gt.0..and.qevap(i).gt.0.) then
                q1(i,k) = q1(i,k) + qevap(i)
                t1(i,k) = t1(i,k) - elocp * qevap(i)
                rn(i) = rn(i) - .001 * qevap(i) * dp / g
                deltv(i) = - elocp*qevap(i)/dt2
                delq(i) =  + qevap(i)/dt2
                delqev(i) = delqev(i) + .001*dp*qevap(i)/g
              endif
              ! Yi-Chi: add output for evaporation Mar 2012: kg/kg/s
              !         remove evaporation from rainout (July12)
                  rainout(i,k)  = rainout(i,k) - qevap(i)/dt2
                  qevapout(i,k) = qevap(i)/dt2
              !---------------------------------------
              dellaq(i,k) = dellaq(i,k) + delq(i) / xmb(i)
              delqbar(i) = delqbar(i) + delq(i)*dp/g
              deltbar(i) = deltbar(i) + deltv(i)*dp/g
            endif
          endif
        enddo
      enddo
!      write(iulog,*) '(Yi-Chi:sas_conv:rainout:',rainout
!      write(iulog,*) '(Yi-Chi:sas_conv:qevapout:',qevapout     
!j
!    do i = 1, im
!     if(me.eq.31.and.cnvflg(i)) then
!     if(cnvflg(i)) then
!       print *, ' deep delhbar, delqbar, deltbar = ',
!    &             delhbar(i),hvap*delqbar(i),cp*deltbar(i)
!       print *, ' deep delubar, delvbar = ',delubar(i),delvbar(i)
!       print *, ' precip =', hvap*rn(i)*1000./dt2
!       print*,'pdif= ',pfld(i,kbcon(i))-pfld(i,ktcon(i))
!     endif
!     enddo
!   Yi-Chi : output the diagnostics used to check the balance
!            of the scheme
!    do i = 1, im
!    if(cnvflg(i)) then
!      write(iulog,*) ' deep delhbar, delqbar, deltbar = ', &
!         delhbar(i),hvap*delqbar(i),cp*deltbar(i)
!      write(iulog,*) ' deep delubar, delvbar = ',delubar(i),delvbar(i)
!      write(iulog,*) ' precip =', hvap*rn(i)*1000./dt2
!      write(iulog,*) 'pdif= ',pfld(i,kbcon(i))-pfld(i,ktcon(i))
!    endif
!    enddo

!  precipitation rate converted to actual precip
!  in unit of m instead of kg
!
      do i = 1, im
        if(cnvflg(i)) then
!
!  in the event of upper level rain evaporation and lower level downdraft
!    moistening, rn can become negative, in this case, we back out of the
!    heating and the moistening
!
          if(rn(i).lt.0..and..not.flg(i)) rn(i) = 0.
          if(rn(i).le.0.) then
            rn(i) = 0.
          else
            ktop(i) = ktcon(i)
            kbot(i) = kbcon(i)
            kcnv(i) = 1
            cldwrk(i) = aa1(i)
            ! +++ Yi-Chi
            ! Yi-Chi : add up all
             lengath = lengath + 1 ! add up # of grids with active convection
             ideep(lengath)= i     ! save the index in ideep
            ! ---
          endif
        endif
      enddo

     delqlbar = 0.0
!
!  cloud water
!
      if (ncloud.gt.0) then
!
      do k = 1, km
        do i = 1, im
          if (cnvflg(i) .and. rn(i).gt.0.) then
            if (k.gt.kb(i).and.k.le.ktcon(i)) then
              tem  = dellal(i,k) * xmb(i) * dt2
              tem1 = max(0.0, min(1.0, (tcr-t1(i,k))*tcrf))
              if (ql(i,k,2) .gt. -999.0) then
                ql(i,k,1) = ql(i,k,1) + tem * tem1            ! ice
                ql(i,k,2) = ql(i,k,2) + tem *(1.0-tem1)       ! water
                dlfsas(i,k) = tem*(1.0-tem1)/dt2 ! Yi-Chi: add to output Mar 2012
              else
                ! +++ Yi-Chi
                dp = 1000. * del(i,k)
                delqlbar(i) = delqlbar(i)+tem*dp/g
                ! ---
                ql(i,k,1) = ql(i,k,1) + tem
                dlfsas(i,k) = tem/dt2 ! Yi-Chi: add to output Mar 2012
              endif
            endif
          endif
        enddo
      enddo
!
      endif

!    do i = 1, im
!    if(cnvflg(i)) then
!      write(iulog,*) 'dt2',dt2
!      write(iulog,*) ' deep delqlbar ', hvap*delqlbar(i)/dt2
!      write(iulog,*) ' precip =', hvap*rn(i)*1000./dt2
!      write(iulog,*) ' precip + delqlbar =', hvap*rn(i)*1000./dt2 + hvap*delqlbar(i)/dt2
!      write(iulog,*) ' difq =', hvap*rn(i)*1000./dt2 + hvap*delqlbar(i)/dt2 + hvap*delqbar(i)
!      !write(iulog,*) 'dlfsas =', dlfsas
!      !write(iulog,*) 'ktcon:',ktcon(i)
!      !write(iulog,*) 'qlko_ktcon:',qlko_ktcon(i)
!    endif
!    enddo

!
      do k = 1, km
        do i = 1, im
          if(cnvflg(i).and.rn(i).le.0.) then
            if (k .le. kmax(i)) then
              t1(i,k) = to(i,k)
              q1(i,k) = qo(i,k)
              u1(i,k) = uo(i,k)
              v1(i,k) = vo(i,k)
            endif
          endif
        enddo
      enddo
! Yi-Chi : add up all 
!      do i=1,im
!       if(cnvflg(i)) then
!             lengath = lengath + 1 ! add up # of grids with active convection 
!             ideep(lengath)= i     ! save the index in ideep
!       endif
!      enddo
! Yi-Chi : output detrain/entrain mass flux for updraft and downdraft
!          aug 14, 2012 : add ktop(=kconvt) level account for detrainment
      do k = 1, km
      do i = 1, im
!        if(cnvflg(i).and.rn(i).gt.0.) then
          if(cnvflg(i)) then
          if(k.ge.kb(i) .and. k.le.ktop(i)) then
          ! Yi-Chi +++
              qlic(i,k) = 0.5*(qlic_tmp(i,k)+qlic_tmp(i,k+1))
              qlic(i,k) = qlic(i,k) * xmb(i) * dt2
          !qlic(i,k)  = qlic(i,k) * xmb(i) * dt2 !kg/kg: (Yi-Chi, Sept 2013)
          ! Yi-Chi ---
!          if(k.ge.kb(i) .and. k.lt.ktop(i)) then
              euout(i,k) = eutmp(i,k) * xmb(i)
              duout(i,k) = dutmp(i,k) * xmb(i)
              edout(i,k) = edtmp(i,k) * xmb(i)
              ddout(i,k) = ddtmp(i,k) * xmb(i)
              duqout(i,k)= duqtmp(i,k)* xmb(i)
          else
          ! Yi-Chi (July, 2012) : 
          !    deep convection has no detrain/entrain at all other levels
              euout(i,k) = 0.0
              duout(i,k) = 0.0
              edout(i,k) = 0.0
              ddout(i,k) = 0.0
              duqout(i,k)= 0.0
              qlic(i,k)  = 0.0 ! Yi-Chi : Sept 2013
          endif
          endif
        enddo
      enddo
!      write(iulog,*) 'Yi-Chi(sas_conv)jmin',jmin
!      write(iulog,*) 'Yi-Chi(sas_conv)etad',etad
!      write(iulog,*) 'Yi-Chi(sas_conv)edto',edto
!      write(iulog,*) 'Yi-Chi(sas_conv)edout',edout
!      write(iulog,*) 'Yi-Chi(sas_conv)ddout',ddout

!
! hchuang code change
!
      do k = 1, km
        do i = 1, im
          if(cnvflg(i).and.rn(i).gt.0.) then
            if(k.ge.kb(i) .and. k.lt.ktop(i)) then
              ud_mf(i,k) = eta(i,k) * xmb(i) * dt2 ! Yi-Chi: (kg m-2 s-1)* s = kg m-2
              ! Yi-Chi : ZM asks for mb/s, Apr 2012
              !          currently Pa/g
              ! Thus the unit needs to be transformed
              ud_mf(i,k) = ud_mf(i,k)*g/(100.*dt2)
            endif
            !write(iulog,*) 'Yi-Chi(sas_conv)kb(i)',kb(i)
            !write(iulog,*) 'Yi-Chi(sas_conv)ud_mf(i,:)',ud_mf(i,:)
          endif
        enddo
      enddo
      do i = 1, im
        if(cnvflg(i).and.rn(i).gt.0.) then
           k = ktop(i)-1
           dt_mf(i,k) = ud_mf(i,k)
        endif
      enddo
      do k = 1, km
        do i = 1, im
          if(cnvflg(i).and.rn(i).gt.0.) then
          ! Yi-Chi: July 2012
          ! set downdraft mass fluxes as zero below cloud base
           if(k.ge.kb(i) .and. k.le.jmin(i)) then
          !  if(k.ge.1 .and. k.le.jmin(i)) then
              dd_mf(i,k) = edto(i) * etad(i,k) * xmb(i) * dt2
              ! Yi-Chi : ZM asks for mb/s : Apr2012
              !          currently Pa/g
              ! Thus the unit needs to be transformed
              dd_mf(i,k) = dd_mf(i,k)*g/(100.*dt2)
            endif
          endif
        enddo
      enddo
!!


     ! ---------------------------------------- !
     ! Yi-Chi: Oct 2014                         !
     ! Writing main diagnostic output variables !
     ! ---------------------------------------- !
     call outfld( 'SAS_qcko'         , qckoout,       im,    lchnk )
     call outfld( 'SAS_duq'          , duqout ,       im,    lchnk )
     call outfld( 'SAS_flgmask'      , flgmask,       im,    lchnk )
     call outfld( 'SAS_cinmask'      , cinmask,       im,    lchnk )
     call outfld( 'SAS_heo'          , heoout ,       im,    lchnk )
     call outfld( 'SAS_heso'         , flgmask,       im,    lchnk )
     
     call outfld( 'SAS_xmb'         , xmbout   ,       im,    lchnk )
     call outfld( 'SAS_xk'          , xkout    ,       im,    lchnk )
     call outfld( 'SAS_fld'         , fldout   ,       im,    lchnk )
     call outfld( 'SAS_dtconv'      , dtconvout,       im,    lchnk )
     call outfld( 'SAS_aa1crit'     , aa1crit  ,       im,    lchnk )
     call outfld( 'SAS_xomega'      , xomega   ,       im,    lchnk )

      return
end subroutine sascnvn

function esat(T,flgforfpvs)
   use module_gfs_funcphys, only : fpvs
   use shr_kind_mod, only: r8 => shr_kind_r8
   use ppgrid
   use physconst,       only: tmelt
   real(r8), intent(in)  :: T                   ! temperature K
   logical,  intent(in)  :: flgforfpvs          ! True : use fpvs
   real(r8)              :: esat                ! vapor pressure
   !real(r8), intent(out) :: esat                ! saturated vapor prs
   real(r8),parameter :: c1e = 6.112_r8
   real(r8),parameter :: c2e = 17.67_r8
   real(r8),parameter :: c3e = 243.5_r8

   !print*,'function esat:Temp',T
   if(flgforfpvs) then
     esat = 0.01 * fpvs(T)
   else
     esat = c1e*exp(c2e*(T-tmelt)/(c3e+T-tmelt))
   endif

!   fpvs
!   Input argument list:
!     t          Real(krealfp) temperature in Kelvin
!
!   Output argument list:
!     fpvs       Real(krealfp) saturation vapor pressure in Pascals
   return
end function esat

!-------------------------------------------------
!   April 2 : Yi-chi 
!   add convtran_sas based on convtran for ZM scheme
!-------------------------------------------------
! subroutine convtran
!end subroutine convtran

subroutine convtran_sas(lchnk   , &
                    doconvtran,q       ,ncnst   ,mu      ,md      , &
                    du      ,eu      ,ed    ,dd  ,  dp   , &
!                    du      ,eu      ,ed      ,dp      ,dsubcld , &
                    jt      ,mx      ,ideep   ,il1g    ,il2g    , &
                    nstep   ,fracis  ,dqdt    ,dpdry   )
!-----------------------------------------------------------------------
! Purpose:
! Convective transport of trace species
! Mixing ratios may be with respect to either dry or moist air
! Method:
! Author: Yi-Chi Wang
!   mainly based on subroutine convtran by P. Rasch
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use constituents,    only: cnst_get_type_byind
   use ppgrid
   use abortutils, only: endrun
   implicit none
!-----------------------------------------------------------------------
!
! Input arguments
!
   integer, intent(in) :: lchnk                 ! chunk identifier
   integer, intent(in) :: ncnst                 ! number of tracers to transport
   logical, intent(in) :: doconvtran(ncnst)     ! flag for doing convective transport
   real(r8), intent(in) :: q(pcols,pver,ncnst)  ! Tracer array including moisture
   real(r8), intent(in) :: mu(pcols,pver)       ! Mass flux up
   real(r8), intent(in) :: md(pcols,pver)       ! Mass flux down
   real(r8), intent(in) :: du(pcols,pver)       ! Mass detraining from updraft
   real(r8), intent(in) :: eu(pcols,pver)       ! Mass entraining from updraft
   real(r8), intent(in) :: ed(pcols,pver)       ! Mass entraining from downdraft
! ------- Yi-Chi: Apr 2012 : add downdraft detrainment for SAS
   real(r8), intent(in) :: dd(pcols,pver)       ! Mass detraining from downdraft
! -------
   real(r8), intent(in) :: dp(pcols,pver)       ! Delta pressure between interfaces
!   real(r8), intent(in) :: dsubcld(pcols)       ! Delta pressure from cloud base to sfc
   real(r8), intent(in) :: fracis(pcols,pver,ncnst) ! fraction of tracer that is insoluble

   integer, intent(in) :: jt(pcols)         ! Index of cloud top for each column
   integer, intent(in) :: mx(pcols)         ! Index of cloud top for each column
   integer, intent(in) :: ideep(pcols)      ! Gathering array
   integer, intent(in) :: il1g              ! Gathered min lon indices over which to operate
   integer, intent(in) :: il2g              ! Gathered max lon indices over which to operate
   integer, intent(in) :: nstep             ! Time step index
   real(r8), intent(in) :: dpdry(pcols,pver)       ! Delta pressure between interfaces


! input/output

   real(r8), intent(out) :: dqdt(pcols,pver,ncnst)  ! Tracer tendency array

!--------------------------Local Variables------------------------------

   integer i                 ! Work index
   integer k                 ! Work index
   integer kbm               ! Highest altitude index of cloud base
   integer kk                ! Work index
   integer kkp1              ! Work index
   integer km1               ! Work index
   integer kp1               ! Work index
   integer ktm               ! Highest altitude index of cloud top
   integer m                 ! Work index

   real(r8) cabv                 ! Mix ratio of constituent above
   real(r8) cbel                 ! Mix ratio of constituent below
   real(r8) cdifr                ! Normalized diff between cabv and cbel
   real(r8) chat(pcols,pver)     ! Mix ratio in env at interfaces
   real(r8) cond(pcols,pver)     ! Mix ratio in downdraft at interfaces
   real(r8) const(pcols,pver)    ! Gathered tracer array
   real(r8) fisg(pcols,pver)     ! gathered insoluble fraction of tracer
   real(r8) conu(pcols,pver)     ! Mix ratio in updraft at interfaces
   real(r8) dcondt(pcols,pver)   ! Gathered tend array
   real(r8) small                ! A small number
   real(r8) mbsth                ! Threshold for mass fluxes
   real(r8) mupdudp              ! A work variable
! Yi-chi: Apr 1, 2012
   real(r8) mdpdudp              ! A work variable for downdraft mass flux- detrained mass
   real(r8) minc                 ! A work variable
   real(r8) maxc                 ! A work variable
   real(r8) fluxin               ! A work variable
   real(r8) fluxout              ! A work variable
   real(r8) netflux              ! A work variable
! Yi-Chi : Apr 20, 2012 : for debugging
   real(r8) fluxin_check(pcols,pver)               ! A work variable
   real(r8) fluxout_check(pcols,pver)               ! A work variable

   real(r8) dutmp(pcols,pver)       ! Mass detraining from updraft
   real(r8) eutmp(pcols,pver)       ! Mass entraining from updraft
   real(r8) edtmp(pcols,pver)       ! Mass entraining from downdraft
! Yi-chi: Apr 1, 2012
   real(r8) ddtmp(pcols,pver)       ! Mass detraining from downdraft
   real(r8) dptmp(pcols,pver)    ! Delta pressure between interfaces
! Yi-Chi: July 12, 2012
   real(r8) compcon
   real(r8) compde
!-----------------------------------------------------------------------
!
   small = 1.e-30_r8
! mbsth is the threshold below which we treat the mass fluxes as zero (in mb/s)
   mbsth = 1.e-15_r8

! Find the highest level top and bottom levels of convection
   ktm = pver
   kbm = pver
   do i = il1g, il2g
      ktm = min(ktm,jt(i))
      kbm = min(kbm,mx(i))
   end do

! Loop ever each constituent
   do m = 2, ncnst
      if (doconvtran(m)) then

         if (cnst_get_type_byind(m).eq.'dry') then
            do k = 1,pver
               do i =il1g,il2g
                  dptmp(i,k) = dpdry(i,k)
                  dutmp(i,k) = du(i,k)*dp(i,k)/dpdry(i,k)
                  eutmp(i,k) = eu(i,k)*dp(i,k)/dpdry(i,k)
                  edtmp(i,k) = ed(i,k)*dp(i,k)/dpdry(i,k)
                  ! Yi-Chi : Apr 1, 2012
                  ddtmp(i,k) = dd(i,k)*dp(i,k)/dpdry(i,k)
               end do
            end do
         else
            do k = 1,pver
               do i =il1g,il2g
                  dptmp(i,k) = dp(i,k)
                  dutmp(i,k) = du(i,k)
                  eutmp(i,k) = eu(i,k)
                  edtmp(i,k) = ed(i,k)
                  ! Yi-Chi : Apr 1, 2012
                  ddtmp(i,k) = dd(i,k)
               end do
            end do
         endif
!        dptmp = dp

! Gather up the constituent and set tend to zero
         do k = 1,pver
            do i =il1g,il2g
               const(i,k) = q(ideep(i),k,m)
               fisg(i,k) = fracis(ideep(i),k,m)
            end do
         end do

! From now on work only with gathered data

! Interpolate environment tracer values to interfaces
         do k = 1,pver
            km1 = max(1,k-1)
            do i = il1g, il2g
               minc = min(const(i,km1),const(i,k))
               maxc = max(const(i,km1),const(i,k))
               if (minc < 0) then
                  cdifr = 0._r8
               else
                  cdifr = abs(const(i,k)-const(i,km1))/max(maxc,small)
               endif

! If the two layers differ significantly use a geometric averaging
! procedure
               if (cdifr > 1.E-6_r8) then
                  cabv = max(const(i,km1),maxc*1.e-12_r8)
                  cbel = max(const(i,k),maxc*1.e-12_r8)
                  chat(i,k) = log(cabv/cbel)/(cabv-cbel)*cabv*cbel

               else             ! Small diff, so just arithmetic mean
                  chat(i,k) = 0.5_r8* (const(i,k)+const(i,km1))
               end if

! Provisional up and down draft values
               conu(i,k) = chat(i,k)
               cond(i,k) = chat(i,k)

!              provisional tends
               dcondt(i,k) = 0._r8

            end do
         end do

! Do levels adjacent to top (for updraft) and bottom (for downdraft)
         k = 2
         km1 = 1
         kk = pver
         do i = il1g,il2g
            mupdudp = mu(i,kk) + dutmp(i,kk)*dptmp(i,kk)
            if (mupdudp > mbsth) then
               conu(i,kk) = (+eutmp(i,kk)*fisg(i,kk)*const(i,kk)*dptmp(i,kk))/mupdudp
            endif
            ! Yi-Chi : Apr 19, 2012
            !  include detrained mass from downdraft mass flux
            if (md(i,k) < -mbsth) then
            cond(i,k) =  (-edtmp(i,km1)*fisg(i,km1)*const(i,km1)*dptmp(i,km1) &
                          +ddtmp(i,km1)*cond(i,km1)*dptmp(i,km1))/md(i,k)
            endif
! CAM-ZM:
!            if (md(i,k) < -mbsth) then
!            cond(i,k) =  (-edtmp(i,km1)*fisg(i,km1)*const(i,km1)*dptmp(i,km1))/md(i,k)
!            endif
         end do

! Updraft from bottom to top
         do kk = pver-1,1,-1
            kkp1 = min(pver,kk+1)
            do i = il1g,il2g
               mupdudp = mu(i,kk) + dutmp(i,kk)*dptmp(i,kk)
               if (mupdudp > mbsth) then
                  conu(i,kk) = (  mu(i,kkp1)*conu(i,kkp1)+eutmp(i,kk)*fisg(i,kk)* &
                                  const(i,kk)*dptmp(i,kk) )/mupdudp
               endif
            end do
         end do

! Downdraft from top to bottom
         do k = 3,pver
            km1 = max(1,k-1)
            do i = il1g,il2g
               if (md(i,k) < -mbsth) then
               ! Yi-Chi : Apr 19 2012
               !  include detrained mass from downdraft mass flux
               cond(i,k) =  (  md(i,km1)*cond(i,km1) &
                               + (-1)*edtmp(i,km1)*fisg(i,km1)*const(i,km1)*dptmp(i,km1) &
                               - (-1)*ddtmp(i,km1)*cond(i,km1)*dptmp(i,km1))/md(i,k)
! CAM-ZM
!                  cond(i,k) =  (  md(i,km1)*cond(i,km1)-edtmp(i,km1)*fisg(i,km1)*const(i,km1) &
!                                  *dptmp(i,km1) )/md(i,k)
                !write(iulog,*) 'Yi-Chi(sas_conv)mx',mx(i)
                !write(iulog,*) 'Yi-Chi(sas_conv)k',k
                !write(iulog,*) 'Yi-Chi(sas_conv)cond1',md(i,km1)*cond(i,km1)
                !write(iulog,*) 'Yi-Chi(sas_conv)cond2',(-1)*edtmp(i,km1)*fisg(i,km1)*const(i,km1)*dptmp(i,km1)
                !write(iulog,*) 'Yi-Chi(sas_conv)cond3',(-1)*ddtmp(i,km1)*cond(i,km1)*dptmp(i,km1)

               endif
            end do
         end do
         !write(iulog,*) 'Yi-Chi(sas_conv)ddtmp',ddtmp
         !write(iulog,*) 'Yi-Chi(sas_conv)edtmp',edtmp

!         write(iulog,*) 'Yi-Chi(sas_conv:before ktm)dcondt',dcondt
         !write(iulog,*) 'Yi-Chi(sas_conv)cond',cond
         !write(iulog,*) 'Yi-Chi(sas_conv)conu',conu

         !write(iulog,*) 'Yi-Chi(sas_conv)num of species: m',m
         do k = ktm,pver
            km1 = max(1,k-1)
            kp1 = min(pver,k+1)
            do i = il1g,il2g
! version 1 hard to check for roundoff errors
!               dcondt(i,k) =
!     $                  +(+mu(i,kp1)* (conu(i,kp1)-chat(i,kp1))
!     $                    -mu(i,k)*   (conu(i,k)-chat(i,k))
!     $                    +md(i,kp1)* (cond(i,kp1)-chat(i,kp1))
!     $                    -md(i,k)*   (cond(i,k)-chat(i,k))
!     $                   )/dp(i,k)

! version 2 hard to limit fluxes
!               fluxin =  mu(i,kp1)*conu(i,kp1) + mu(i,k)*chat(i,k)
!     $                 -(md(i,k)  *cond(i,k)   + md(i,kp1)*chat(i,kp1))
!               fluxout = mu(i,k)*conu(i,k)     + mu(i,kp1)*chat(i,kp1)
!     $                 -(md(i,kp1)*cond(i,kp1) + md(i,k)*chat(i,k))

! version 3 limit fluxes outside convection to mass in appropriate layer
! these limiters are probably only safe for positive definite quantitities
! it assumes that mu and md already satify a courant number limit of 1
!               fluxin =  mu(i,kp1)*conu(i,kp1)+ mu(i,k)*min(chat(i,k),const(i,km1)) &
!                         -(md(i,k)  *cond(i,k) + md(i,kp1)*min(chat(i,kp1),const(i,kp1)))
!               fluxout = mu(i,k)*conu(i,k) + mu(i,kp1)*min(chat(i,kp1),const(i,k)) &
!                         -(md(i,kp1)*cond(i,kp1) + md(i,k)*min(chat(i,k),const(i,k)))
!               if (abs(netflux) < max(fluxin,fluxout)*1.e-12_r8) then
!                  netflux = 0._r8
!               endif
!               dcondt(i,k) = netflux/dptmp(i,k)

! Yi-Chi:
!        term of compensation
     compcon = mu(i,k)*min(chat(i,k),const(i,km1)) + md(i,kp1)*min(chat(i,kp1),const(i,kp1)) &
               - mu(i,kp1)*min(chat(i,kp1),const(i,k)) - md(i,k)*min(chat(i,k),const(i,k))
     compde  = dutmp(i,k)*conu(i,k)*dptmp(i,k) + ddtmp(i,k)*cond(i,k)*dptmp(i,k)
     netflux = compcon + compde
     dcondt(i,k) = netflux/dptmp(i,k)
     ! write(iulog,*) 'Yi-Chi(sas_conv)k',k 
     ! write(iulog,*) 'Yi-Chi(sas_conv)comp detrainment',compde
     ! write(iulog,*) 'Yi-Chi(sas_conv)dutmp',dutmp(i,k),ddtmp(i,k)
     ! write(iulog,*) 'Yi-Chi(sas_conv)conu,cond',conu(i,k),cond(i,k)
     ! write(iulog,*) 'Yi-Chi(sas_conv)comp convection',compcon
     ! write(iulog,*) 'Yi-Chi(sas_conv)mu,md(i,k)',mu(i,k),md(i,k)
     ! Yi-Chi : Apr 20, 2012 debugging
     !          fluxin_check(i,k) = fluxin
     !          fluxout_check(i,k) = fluxout
! Yi-Chi's version for convective transport : july 11, 2012
!    : account for residual mass and detrainment into 
!    : 
!               write(iulog,*) 'Yi-Chi(sas_conv)k',k
!               write(iulog,*) 'Yi-Chi(sas_conv)mt,mu,md at k',dptmp(i,k)/gravit,mu(i,k),md(i,k)
!               write(iulog,*) 'Yi-Chi(sas_conv)in:mu*conu',mu(i,kp1)*conu(i,kp1)
!               write(iulog,*) 'Yi-Chi(sas_conv)in:mu*chat',mu(i,k)*min(chat(i,k),const(i,km1))
!               write(iulog,*) 'Yi-Chi(sas_conv)in:-md*cond',-1*md(i,k)*cond(i,k)
!               write(iulog,*) 'Yi-Chi(sas_conv)in:-md*chat',-1*md(i,kp1)*min(chat(i,kp1),const(i,kp1))
!               write(iulog,*) 'Yi-Chi(sas_conv)out:mu*conu',mu(i,k)*conu(i,k)
!               write(iulog,*) 'Yi-Chi(sas_conv)out:mu*chat',mu(i,kp1)*min(chat(i,kp1),const(i,k))
!               write(iulog,*) 'Yi-Chi(sas_conv)out:-md*cond',-1*md(i,kp1)*cond(i,kp1)
!               write(iulog,*) 'Yi-Chi(sas_conv)out:-md*chat',-1*md(i,k)*min(chat(i,k),const(i,k))
!               write(iulog,*) 'Yi-Chi(sas_conv)fluxin,fluxout',fluxin,fluxout
!               netflux = fluxin - fluxout
! Yi-Chi : July 11, 2012
!          add residual mass fluxes

            end do
         end do
 
!         write(iulog,*) 'Yi-Chi(sas_conv:before kbm)dcondt',dcondt
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!DIR$ NOINTERCHANGE
         do k = kbm,pver
            km1 = max(1,k-1)
            do i = il1g,il2g
               if (k == mx(i)) then

! version 1
!                  dcondt(i,k) = (1./dsubcld(i))*
!     $              (-mu(i,k)*(conu(i,k)-chat(i,k))
!     $               -md(i,k)*(cond(i,k)-chat(i,k))
!     $              )

! version 2
!                  fluxin =  mu(i,k)*chat(i,k) - md(i,k)*cond(i,k)
!                  fluxout = mu(i,k)*conu(i,k) - md(i,k)*chat(i,k)
! version 3
!                  fluxin =  mu(i,k)*min(chat(i,k),const(i,km1)) - md(i,k)*cond(i,k)
!                  fluxout = mu(i,k)*conu(i,k) - md(i,k)*min(chat(i,k),const(i,k))
!                  ! Yi-Chi : Apr 20, 2012 debugging
!                  fluxin_check(i,k) = fluxin
!                  fluxout_check(i,k) = fluxout

!                  netflux = fluxin - fluxout
!                  if (abs(netflux) < max(fluxin,fluxout)*1.e-12_r8) then
!                     netflux = 0._r8
!                  endif
!!                  dcondt(i,k) = netflux/dsubcld(i)
          compcon = mu(i,k)*min(chat(i,k),const(i,km1)) - md(i,k)*min(chat(i,k),const(i,k))
          compde  = dutmp(i,k)*conu(i,k)*dptmp(i,k) + ddtmp(i,k)*cond(i,k)*dptmp(i,k)
          netflux = compcon + compde
                  
                  dcondt(i,k) = netflux/dptmp(i,k)
               else if (k > mx(i)) then
!                  dcondt(i,k) = dcondt(i,k-1)
                  dcondt(i,k) = 0._r8
               end if
            end do
         end do
        
         !write(iulog,*) 'Yi-Chi(sas_conv)const',const 
         !write(iulog,*) 'Yi-Chi(sas_conv)dcondt',dcondt 
         !write(iulog,*) 'Yi-Chi(sas_conv)fluxin_check',fluxin_check
         !write(iulog,*) 'Yi-Chi(sas_conv)fluxout_check',fluxout_check
! Initialize to zero everywhere, then scatter tendency back to full array
         dqdt(:,:,m) = 0._r8
         do k = 1,pver
            kp1 = min(pver,k+1)
!DIR$ CONCURRENT
            do i = il1g,il2g
               dqdt(ideep(i),k,m) = dcondt(i,k)
            end do
         end do

      end if      ! for doconvtran

   end do

   return

end subroutine convtran_sas

REAL FUNCTION qrch2nd(gamma, dbyotmp,totmp,qesotmp)
! +++
! Yi-Chi (Jan 3, 2013):
!  add 2nd derivative for qrch (q_sat of the plume)
! hvap*hvap/(rv*cp)
  use physconst,       only: grav=>gravit,cp=>cpair, hvap=> latvap, &
                                 rv=>rh2o, rair, t0c=>tmelt
!  use physconst,       only: cpair, epsilo, gravit, latvap, tmelt, rair, &
!                             cpwv, rh2o
!use module_gfs_physcons, grav => con_g, cp => con_cp, hvap => con_hvap, &
!                   rv => con_rv, fv => con_fvirt
use module_gfs_machine ,  only : kind_phys
!      parameter(grav=gravit,cp=cpair, hvap=latvap)
!      parameter(rv=rh2o, t0c=tmelt, cvap= cpwv, cliq=cpliq,fv=(rh2o/rair)-1)

real(kind=kind_phys) A, B, C, F, detm
real(kind=kind_phys) gamma, dbyotmp,totmp,qesotmp
!real(kind=kind_phys) fv

!    grav=gravit
!    cp=cpair
!    hvap=latvap
!    rv=rh2o
!    t0c=tmelt
!    cvap= cpwv
!    cliq=cpliq
!    fv=(rv/rair)-1

    qrch2nd = 0.0
    A = (cp*gamma/(hvap*totmp**2) - 2*qesotmp/(totmp**3)) &
         /((2*cp*cp*rv)/(hvap*hvap))
    B = gamma + 1
    C = -1*dbyotmp
    detm = B**2 - 4*A*C
    !write(iulog,*) 'Yi-Chi(sas_conv)detm',detm
    if( detm .gt. 0.0 ) then
      if(A.eq.0.0) then
      qrch2nd  = qesotmp &
                 + (gamma*dbyotmp) / (hvap*(1.+gamma))
      else
      F = (-B + sqrt(detm)) / (2*A)
      qrch2nd  = qesotmp &
                 + (dbyotmp - F) / hvap
      endif
    endif
    if( detm .eq. 0.0 ) then
      F = -B/(2*A)
      qrch2nd  = qesotmp &
                 + (dbyotmp - F) / hvap
    endif
    if( detm .lt. 0.0 ) then
      !print*,'(Yi-Chi:sas_conv.F90,line 2735) this point detm < 0'
      qrch2nd  = qesotmp &
                 + (gamma*dbyotmp) / (hvap*(1.+gamma))
    endif
              !qrch = qeso(i,k)  &
              !    + gamma * dbyo(i,k) / (hvap * (1. + gamma))
   RETURN
end function qrch2nd

end module
