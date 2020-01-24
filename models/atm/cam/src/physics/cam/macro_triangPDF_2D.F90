
  module cldfrc_triang_macro

   use shr_kind_mod,     only: r8 => shr_kind_r8
   use ppgrid,           only: pcols, pver, pverp
!!cjshiu cesm1.0.3   use wv_saturation,    only: aqsat,qsat_ice
!   use wv_saturation,    only: qsat_water,qsat_ice

   implicit none
   private
   save

   public :: macro_triangpdf_2d ! calculation 
   public ::  astG_PDF_caldelta_single
   private:: cubic_exact

! ycwang revise the subroutine to triangular PDF : Feb 2015
   contains

   ! ------------------ !
   ! macro_triangpdf_2d !
   ! ------------------ !

subroutine macro_triangpdf_2d(ncol,  play, tlay, qcwat, qclwat, iccwat, &
           deltalpdf, cldlfupdf, deltaipdf, cldifupdf, deltatpdf, cldtfupdf, cldmfupdf, qscwatupdf, qsciceupdf)

!   use shr_kind_mod,     only: r8 => shr_kind_r8
!   use ppgrid,           only: pcols, pver, pverp
!cjshiu cesm1.0.3   use wv_saturation,    only: aqsat,qsat_ice
   use wv_saturation,    only: qsat_water,qsat_ice
   
  implicit none
!cjshiu   save
   
!-----------------
! input argument
!-----------------  
  integer , intent(in) :: ncol
!cjshiu  integer , intent(in) :: icol
!cjshiu  integer , intent(in) :: nlay
  real(r8), intent(in) :: play(pcols,pver)     ! layer pressure (Pa)
  real(r8), intent(in) :: tlay(pcols,pver)     ! layer temperature (K)
  real(r8), intent(in) :: qcwat(pcols,pver)    ! grid mean water vapor (kg/kg)
  real(r8), intent(in) :: qclwat(pcols,pver)   ! cloud liquid water (kg/kg)
  real(r8), intent(in) :: iccwat(pcols,pver)   ! cloud ice water (kg/kg)
!cjshiu  real(r8), intent(in) :: play(ncol,pver)     ! layer pressure (Pa)
!cjshiu  real(r8), intent(in) :: tlay(ncol,pver)     ! layer temperature (K)
!cjshiu  real(r8), intent(in) :: qcwat(ncol,pver)    ! grid mean water vapor (kg/kg)
!cjshiu  real(r8), intent(in) :: qclwat(ncol,pver)   ! cloud liquid water (kg/kg)
!cjshiu  real(r8), intent(in) :: iccwat(ncol,pver)   ! cloud ice water (kg/kg)
      
!-----------------
! output argument
!-----------------  
  real(r8), intent(out) :: deltalpdf(pcols,pver)         ! liquid cloud delta value
  real(r8), intent(out) :: cldlfupdf(pcols,pver)         ! liquid cloud fraction
  real(r8), intent(out) :: deltaipdf(pcols,pver)         ! ice cloud delta value
  real(r8), intent(out) :: cldifupdf(pcols,pver)         ! ice cloud fraction
  real(r8), intent(out) :: deltatpdf(pcols,pver)         ! total cloud water (ice+liquid) delta value
  real(r8), intent(out) :: cldtfupdf(pcols,pver)         ! total cloud fraction computed by (ice+liquid) pdf
  real(r8), intent(out) :: cldmfupdf(pcols,pver)         ! max of icd & liquid cloud fraction 
  real(r8), intent(out) :: qscwatupdf(pcols,pver)        ! saturation specific humidity over water 
  real(r8), intent(out) :: qsciceupdf(pcols,pver)        ! saturation specific humidity over ice 
!cjshiu  real(r8), intent(out) :: deltalpdf(ncol,pver)         ! liquid cloud delta value
!cjshiu  real(r8), intent(out) :: cldlfupdf(ncol,pver)         ! liquid cloud fraction
!cjshiu  real(r8), intent(out) :: deltaipdf(ncol,pver)         ! ice cloud delta value
!cjshiu  real(r8), intent(out) :: cldifupdf(ncol,pver)         ! ice cloud fraction
!cjshiu  real(r8), intent(out) :: deltatpdf(ncol,pver)         ! total cloud water (ice+liquid) delta value
!cjshiu  real(r8), intent(out) :: cldtfupdf(ncol,pver)         ! total cloud fraction computed by (ice+liquid) pdf
!cjshiu  real(r8), intent(out) :: cldmfupdf(ncol,pver)         ! max of icd & liquid cloud fraction 
        
!-----------------  
! local variable
!-----------------  
!cjshiu  real(r8), parameter :: qlmin = 2.e-18_r8                 ! minimum conc. of cloud liquid mass (as qsmall of NCEP value of 1.e-10 or 2.e-12 is used)
!cjshiu  real(r8), parameter :: qlmin = 2.e-12_r8                 ! minimum conc. of cloud liquid mass (as qsmall of NCEP value of 1.e-10 or 2.e-12 is used)
  real(r8), parameter :: qlmin = 2.e-10_r8                 ! minimum conc. of cloud liquid mass (as qsmall of NCEP value of 1.e-10 or 2.e-12 is used)
!cjshiu  real(r8), parameter :: qimin = 2.e-18_r8                 ! minimum conc. of cloud ice mass 
!cjshiu  real(r8), parameter :: qimin = 2.e-12_r8                 ! minimum conc. of cloud ice mass 
  real(r8), parameter :: qimin = 2.e-10_r8                 ! minimum conc. of cloud ice mass 
!cjshiu  real(r8), parameter :: qimin = 1.e-8_r8                 ! minimum conc. of cloud ice mass 
!cjshiu  real(r8), parameter :: qimin = 1.e-7_r8                 ! minimum conc. of cloud ice mass 
        
  real(r8), parameter :: qtmin = qlmin + qimin             ! minimum conc. of cloud ice+liquid mass 

!cjshiu  real(r8), parameter :: sup = 1.1_r8
!cjshiu  real(r8), parameter :: sup = 1.2_r8
  real(r8), parameter :: sup = 1.0_r8
      
!cjshiu  real(r8), parameter :: rhcrit = 0.5_r8
!cjshiu  real(r8), parameter :: rhcrit = 0.8_r8
!cjshiu  real(r8), parameter :: rhcrit = 0.83_r8
  real(r8), parameter :: rhcrit = 0.85_r8
!cjshiu  real(r8), parameter :: rhcrit = 0.7_r8
!cjshiu  real(r8), parameter :: rhcrit = 0.9_r8

  real(r8) :: esat(pcols,pver)     ! saturation vapor pressure
  real(r8) :: qsat(pcols,pver)     ! saturation specific humidity
  real(r8) :: eisat(pcols,pver)     ! saturation vapor pressure over ice
  real(r8) :: qisat(pcols,pver)     ! saturation specific humidity over ice
  real(r8) :: qscwat(pcols,pver)   ! saturation specific humidity over water
  real(r8) :: qscice(pcols,pver)   ! saturation specific humidity over ice
!cjshiu  real(r8) :: esat(ncol,pver)     ! saturation vapor pressure
!cjshiu  real(r8) :: qsat(ncol,pver)     ! saturation specific humidity
!cjshiu  real(r8) :: qscwat(ncol,pver)   ! saturation specific humidity over water

  real(r8) :: qlbar                ! cloud liquid water (no ice included)
  real(r8) :: qibar                ! cloud ice water (cloud ice only)
  real(r8) :: qtbar                ! cloud ice+liquid water 
  real(r8) :: qvbar                ! grid mean water vapor
  real(r8) :: qsc                  ! saturation specific humidity over water
!cjshiu  real(r8) :: qsci                  ! saturation specific humidity over ice
  real(r8) :: deltasqrt            ! square root of delta of uniform PDF
  real(r8) :: fice                 ! ratio of cloud ice water to total cloud condensed water (ice/(ice+liquid))
  real(r8) :: rh,rhi

  integer :: i
  integer :: k

! +++ Yi-Chi +++ !
! add output for subroutine astG_PDF_caldelta
! --- Yi-Chi --- !
  real(r8) :: a_out            ! cloud fraction calculated
  real(r8) :: delta_out        ! PDF width calculated
  real(r8) :: insideacos       ! dummy argument for calculating 

!----------------------------------------------

!cjshiu   i = icol
!cjshiu   k = nlay

!cjshiu add initial values
   esat(1:ncol,1:pver) = 0._r8     ! saturation vapor pressure over water
   qsat(1:ncol,1:pver) = 0._r8    ! saturation specific humidity over water
   eisat(1:ncol,1:pver) = 0._r8     ! saturation vapor pressure over ice
   qisat(1:ncol,1:pver) = 0._r8    ! saturation specific humidity over ice

   qscwat(1:ncol,1:pver) = 0._r8   ! saturation specific humidity over water
   qscice(1:ncol,1:pver) = 0._r8   ! saturation specific humidity over ice
   !qclwat(1:ncol,1:nlay) = 0._r8  ! cloud liquid water (no ice included)

   qlbar = 0._r8               ! cloud liquid water (no ice included)
   qibar = 0._r8               ! cloud ice water (cloud ice only)
   qvbar = 0._r8               ! grid mean water vapor
   qsc   = 0._r8               ! saturation specific humidity over water
!cjshiu   qsci   = 0._r8               ! saturation specific humidity over ice
   deltasqrt = 0._r8            ! square root of delta of uniform PDF
!cjshiu add end

   deltatpdf(:,:) = 0._r8                !cjshiu need to check if 1:ncol and 1:pver should be used for global simulation
   cldtfupdf(:,:) = 0._r8
   deltaipdf(:,:) = 0._r8
   cldifupdf(:,:) = 0._r8
   deltalpdf(:,:) = 0._r8
   cldlfupdf(:,:) = 0._r8
   qscwatupdf(:,:) = 0._r8
   qsciceupdf(:,:) = 0._r8
 
!cjshiu    print*,'cjshiu in macro_uniformPDF.F90  i= k= ncol= pcol= pver= ',i,k,ncol,pcols,pver 
!cjshiu    print*,'cjshiu in macro_uniformPDF_2D.F90  play= tlay= qcwat= qclwat= iccwat= ',play,tlay,qcwat,qclwat,iccwat 
   do k=1,pver
    do i=1,ncol
      call qsat_water(tlay(i,k), play(i,k), esat(i,k), qsat(i,k))
      qscwat(i,k) = qsat(i,k)
      qscwatupdf(i,k) = qscwat(i,k)
    enddo 
   enddo
!cjshiu cesm1.0.3   call aqsat(tlay, play, esat, qsat, pcols, &
!cjshiu cesm1.0.3               ncol, pver, 1, pver)
!cjshiu    print*,'after calling aqsat in macro_uniform tlay= play= east= qsat= ',tlay, play, esat, qsat

!cjshiu cesm1.0.3    qscwat(1:ncol,1:pver) = qsat(1:ncol,1:pver)
!cjshiu add qscice for relax sup assumption
   do k=1,pver
    do i=1,ncol
     call qsat_ice(tlay(i,k),play(i,k),eisat(i,k),qisat(i,k))
     qscice(i,k) = qisat(i,k)
     qsciceupdf(i,k) = qscice(i,k)
    enddo 
   enddo
!cjshiu end add qscice

   !qscwatupdf(1:ncol,1:nlay) = qscwat(1:ncol,1:nlay)   !cjshiu put to physical buffer for astG_PDF uniform PDF

   !*** yihsuan cloud liquid+ice water ***
   ! initialize output 
!cjshiu move to upper 
!cjshiu   deltatpdf(:,:) = 0._r8                !cjshiu need to check if 1:ncol and 1:pver should be used for global simulation
!cjshiu   cldtfupdf(:,:) = 0._r8
!cjshiu   deltaipdf(:,:) = 0._r8
!cjshiu   cldifupdf(:,:) = 0._r8
!cjshiu   deltalpdf(:,:) = 0._r8
!cjshiu   cldlfupdf(:,:) = 0._r8
!cjshiu   cldmfupdf(:,:) = 0._r8

   do k=1,pver
    do i=1,ncol
      qtbar = qclwat(i,k) + iccwat(i,k)
      qvbar = qcwat(i,k)

      !print*,'macro_uniformPDF.F90 total cld: qtbar= qvbar=', qvbar,qlbar
      ! if cloud water larger than threshold value,
      ! compute PDF
      if (qtbar .GT. qtmin) then

        ! cloud ice ratio
        fice = iccwat(i,k)/qtbar

        ! saturation vapor pressure is linear combination of cloud ice and liquid
!cjshiu        qsc  = (1._r8-fice)*qscwat(i,k) + fice*qscwat(i,k)*sup
        qsc  = (1._r8-fice)*qscwat(i,k) + fice*qscice(i,k)*sup

        if (qsc .LE. qvbar) then
!cjshiu          qvbar = 0.9999_r8*qsc
          qvbar = qsc
        endif

        ! +++ Yi-Chi +++ !
        ! ycwang: add triang PDF
        !call astG_PDF_caldelta( U_in, p_in, ql_in, qv_in, landfrac_in, snowh_in, a_out, Ga_out, delta_out, ncol )
       call astG_PDF_caldelta_single( qsc, qtbar, qvbar, a_out, delta_out)
       deltatpdf(i,k) = delta_out
       cldtfupdf(i,k) = a_out
       !print*,'macro_uniformPDF.F90 total cld: Case: qtbar>qtmin, delta=, a=', delta_out,a_out
       ! <<uniform>>
       ! deltasqrt = ( sqrt(qtbar) + sqrt((qsc-qvbar)) )
       ! deltatpdf(i,k) = deltasqrt * deltasqrt
       ! cldtfupdf(i,k) = 1._r8/(2._r8*deltatpdf(i,k))*(qtbar+qvbar+deltatpdf(i,k)-qsc)
       ! --- Yi-Chi --- !

      ! if cloud water smaller than threshold value,
      ! use relative humidity to determine cloud fraction
      else

        ! grid-scale relative humidity
        rh = qvbar/qscwat(i,k)
        !print*,'macro_uniformPDF.F90 total cld: Case: qtbar<qtmin: use RH, rh= ', rh
!cjshiu          print*,'cjshiu in macro_uniformPDF.F90 ice phase qvbar=  qscwat= i k ', qvbar,qscwat(i,k),i,k
        delta_out    = 1._r8 - rhcrit
        ! liquid phase cloud
        if ( rh .GE. rhcrit .and. tlay(i,k) .GT. 273.15_r8) then
        ! +++ Yi-Chi +++ !
        ! << uniform >>
         ! cldtfupdf(i,k) = 1._r8 - ((1._r8-rh)/(1._r8-rhcrit))**0.5
        insideacos = 1._r8 + (2._r8*sqrt(2._r8)/3._r8-1._r8)*delta_out
        if(rh.ge.insideacos) then
        cldtfupdf(i,k) = 1._r8
        else
        cldtfupdf(i,k) = 4._r8*(cos((1._r8/3._r8)*(acos((3._r8/2._r8/sqrt(2._r8))* &
                      (1._r8+(rh-1.0_r8)/delta_out))-2._r8*3.141592_r8)))**2._r8
        endif
        !print*,'macro_uniformPDF.F90 total cld: Case: rh>rhcrit: liq cld=', cldlfupdf(i,k)
        ! --- Yi-Chi --- !
        ! ice phase cloud
        elseif ( rh .GE. rhcrit .and. tlay(i,k) .LE. 273.15_r8) then
!cjshiu          print*,'cjshiu in macro_uniformPDF.F90 ice phase rh= ',rh
     !cjshiu add the following to fix when RH > 1.1 for clouf ice e.g. rh = 1.278     
          rh = min(rh,1.1_r8) 
        ! +++ Yi-Chi +++ !
        insideacos = 1._r8 + (2._r8*sqrt(2._r8)/3._r8-1._r8)*delta_out
        if(rh.ge.insideacos) then
        cldtfupdf(i,k) = 1._r8
        else
        cldtfupdf(i,k) = 4._r8*(cos((1._r8/3._r8)*(acos((3._r8/2._r8/sqrt(2._r8))* &
                      (1._r8+(rh-1.0_r8)/delta_out))-2._r8*3.141592_r8)))**2._r8
        endif
        !print*,'macro_uniformPDF.F90 total cld: Case: rh>rhcrit: ice cld=', cldlfupdf(i,k)
        ! <<uniform>>
        !  cldtfupdf(i,k) = 1._r8 - ((1.1_r8-rh)/(1.1_r8-rhcrit))**0.5
!cjshiu check          cldtfupdf(i,k) =  ((1.1_r8-rh)/(1.1_r8-rhcrit))**2._r8
!cjshiu need to check the above formula
        ! --- Yi-Chi --- !

        ! rh < rh_critical, no cloud
        else
          cldtfupdf(i,k) = 0._r8

        endif
        !print*,'macro_uniformPDF.F90 total cld: qtbar<qtmin, delta=, a=', delta_out,cldtfupdf(i,k)

      endif
     enddo
   enddo
   !*** yihsuan cloud liquid+ice water end ***

   !cjshiu cloud liquid water phase
    do k=1,pver
     do i=1,ncol
         qlbar = qclwat(i,k)
         qvbar = qcwat(i,k)
         qsc   = qscwat(i,k)
         !print*,'macro_uniformPDF.F90 liq cld: qvbar= qlbar  qsc=', qvbar,qlbar,qsc
        !cjshiu there is cases that qsc < qvbar (set to qvbar slightly lower than qsc to avoid floating invalid)
          if (qsc .LE. qvbar) then
!cjshiu             qvbar = 0.9999_r8*qsc
!cjshiu            print*,'cjshiu i, k, qsc, qvbar, qlbar=',i,k,qsc,qvbar,qlbar
            qvbar = qsc
          endif
        !print*,'qlbar= qvbar= qsc= ',qlbar,qvbar,qsc
      
      rh = qvbar/qsc
      !print*,'macro_uniformPDF.F90 liq cld: rh=', rh

      if (qlbar .GT. qlmin) then
        ! +++ Yi-Chi +++ !
        ! ycwang: add triang PDF
        !      call astG_PDF_caldelta( U_in, p_in, ql_in, qv_in, landfrac_in, snowh_in, a_out, Ga_out, delta_out, ncol )
        call astG_PDF_caldelta_single( qsc, qlbar, qvbar, a_out, delta_out )
        deltalpdf(i,k) = delta_out
        cldlfupdf(i,k) = a_out
        !print*,'macro_uniformPDF.F90 liq cld: Case: qlbar>qlmin: delta_out=, a_out=', delta_out, a_out
        ! <<uniform>>
        !deltasqrt = ( sqrt(qlbar) + sqrt((qsc-qvbar)) )
        !!print*,'delatasqrt= ',deltasqrt
        !deltalpdf(i,k) = deltasqrt * deltasqrt
        !!cjshiu may need to add zero check for delta
        !cldlfupdf(i,k) = 1._r8/(2._r8*deltalpdf(i,k))*(qlbar+qvbar+deltalpdf(i,k)-qsc)
        !!print*,'i=  k= deltalpdf(i,k)= cldlfupdf(i,k) = ',i,k,deltalpdf(i,k),cldlfupdf(i,k)
        ! --- Yi-Chi --- !
      elseif (qlbar .LE. qlmin) then
        !print*,'macro_uniformPDF.F90 liq cld: Case: qlbar<qlmin'
        delta_out   = 1.0_r8 - rhcrit
       if (rh .GT. rhcrit) then    ! case : RH > RHcrit
        !dV    = 1.0_r8 - rhcrit
        ! Yi-Chi adds 
        deltalpdf(i,k) = delta_out
        insideacos = 1._r8 + (2._r8*sqrt(2._r8)/3._r8-1._r8)*delta_out
        if(rh.ge.insideacos) then
        cldlfupdf(i,k) = 1._r8
        else
        cldlfupdf(i,k) = 4._r8*(cos((1._r8/3._r8)*(acos((3._r8/2._r8/sqrt(2._r8))* &
                      (1._r8+(rh-1.0_r8)/delta_out))-2._r8*3.141592_r8)))**2._r8
        endif
        !print*,'macro_uniformPDF.F90 liq cld: Case: rh>rhcrit, delta_out=, a_out=', delta_out, cldlfupdf(i,k)
        !deltalpdf(i,k) = qsc*(1._r8-rhcrit)
        !cldlfupdf(i,k) = 1._r8 - ((1._r8-rh)/(1._r8-rhcrit))**0.5
        ! --- Yi-Chi --- !
       elseif (rh .LE. rhcrit) then! case : RH < RHcrit 
        deltalpdf(i,k) = delta_out
        !cjshiu should test for the rhminl and rhminh as those of Park's usage, first assume to be 0.8
        !deltalpdf(i,k) = qsc*(1._r8-rhcrit)      !cjshiu which value can be used??
        !cjshiu may need to consider assume a very small value for cldfraction, first set 0
        cldlfupdf(i,k) = 0._r8
        !print*,'macro_uniformPDF.F90 liq cld: Case: rh<rhcrit, delta_out=, a_out=', delta_out, cldlfupdf(i,k)
        else
        !print*,'need to check why this will happen i,k,rh= ',i,k,rh
       endif  
      endif
     enddo
    enddo

   !cjshiu cloud ice phase (need to check the possible mixed phase clouds??)
    do k=1,pver
     do i=1,ncol
         qibar = iccwat(i,k)
         qvbar = qcwat(i,k)
!cjshiu         qsc   = qscwat(i,k)*sup
         qsc   = qscice(i,k)*sup
         !print*,'macro_uniformPDF.F90 ice cld: qvbar= qibar=  qsc=', qvbar,qibar,qsc
        !cjshiu there is cases that qsc < qvbar (set to qvbar slightly lower than qsc to avoid floating invalid)
        !cjshiu might consider increasing value of sup to fix the above situation (what's the physical meaning? behind this)
!cjshiu debug          if (qsc .LE. qvbar) then
!cjshiu debug            qvbar = 0.98_r8*qsc
!cjshiu debug          endif
          if (qscice(i,k) .LE. qvbar) then
!cjshiu old           qvbar = 0.9999_r8*qsc
            qvbar = qscice(i,k)
          endif
!cjshiu        print*,'cjshiu qibar= qvbar= qsc= qscice=',qibar,qvbar,qsc,qscice(i,k)

!cjshiu      rhi = qvbar/qsc
        rhi = qvbar/qscice(i,k)
        !print*,'macro_uniformPDF.F90 ice cld: rhi=', rhi
!cjshiu        rhi = (qvbar+qibar)/qscice(i,k)
!cjshiu check
!cjshiu        print*,'cjshiu rhi in uniform pdf i, k, rhi=',i,k,rhi

      if (qibar .GT. qimin) then
        ! +++ Yi-Chi +++ !
        ! ycwang: add triang PDF
        !      call astG_PDF_caldelta( U_in, p_in, ql_in, qv_in, landfrac_in, snowh_in, a_out, Ga_out, delta_out, ncol )
        call astG_PDF_caldelta_single( qsc, qibar, qvbar, a_out, delta_out )
        deltaipdf(i,k) = delta_out
        cldifupdf(i,k) = a_out
        !print*,'macro_uniformPDF.F90 ice cld: Case: qibar>qimin, delta_out=, a_out=', delta_out, cldifupdf(i,k)
        ! <<uniform>>
        !deltasqrt = ( sqrt(qibar) + sqrt((qsc-qvbar)) )
        !!print*,'delatasqrt= ',deltasqrt
        !deltaipdf(i,k) = deltasqrt * deltasqrt
        !!cjshiu may need to add zero check for delta
        !cldifupdf(i,k) = 1._r8/(2._r8*deltaipdf(i,k))*(qibar+qvbar+deltaipdf(i,k)-qsc)
        !!print*,'i=  k= deltaipdf(i,k)= cldifupdf(i,k) = ',i,k,deltaipdf(i,k),cldifupdf(i,k)
        ! --- Yi-Chi --- !
      elseif (qibar .LE. qimin) then
        !print*,'macro_uniformPDF.F90 ice cld: Case: qibar<qimin'
        delta_out    = 1._r8 - rhcrit
       if (rhi .GT. rhcrit) then
        ! +++ Yi-Chi +++ !
        deltaipdf(i,k) = delta_out
        !print*,'macro_uniformPDF.F90 ice cld: acos()=',(3._r8/2._r8/sqrt(2._r8))* &
        !              (1._r8+(rhi-1.0_r8)/delta_out)
        ! add a condition to avoid acos >1
        ! ycwang : this may be another tunable part for ice clouds.
        ! need to think what is the right critical values for ice clouds.
        insideacos = 1._r8 + (2._r8*sqrt(2._r8)/3._r8-1._r8)*delta_out
        if(rhi.ge.insideacos) then
        cldifupdf(i,k) = 1._r8
        else
        cldifupdf(i,k) = 4._r8*(cos((1._r8/3._r8)*(acos((3._r8/2._r8/sqrt(2._r8))* &
                      (1._r8+(rhi-1.0_r8)/delta_out))-2._r8*3.141592_r8)))**2._r8
        endif
         !print*,'macro_uniformPDF.F90 ice cld: Case: rh>rhcrit, delta_out=, a_out=', delta_out, cldifupdf(i,k)
        ! <<uniform>>
        ! deltaipdf(i,k) = qsc*(1._r8-rhcrit)
        !!cjshiu        cldifupdf(i,k) = 1._r8 - ((1.1_r8-rhi)/(1.1_r8-rhcrit))**0.5
        !cldifupdf(i,k) = 1._r8 - ((1._r8-rhi)/(1._r8-rhcrit))**0.5 !cjshiu if add qibar to rhi, rhi > 1.0 can occur and cause?
        ! --- Yi-Chi --- !
        elseif (rhi .LE. rhcrit) then
        ! +++ Yi-Chi +++ !
        ! <<triangular>>
        deltaipdf(i,k) = delta_out
        ! <<uniform>>
        !cjshiu should test for the rhminl and rhminh as those of Park's uasage, first assume to be 0.8
        !deltaipdf(i,k) = qsc*(1._r8-rhcrit)  !cjshiu how about assuming to be zero
        !cjshiu may need to consider assume a very small value for cldfraction, first set 0
        ! --- Yi-Chi --- !
        cldifupdf(i,k) = 0._r8
        !print*,'macro_uniformPDF.F90 ice cld: Case: rh<rhcrit, delta_out=, a_out=', delta_out, cldifupdf(i,k)
        else
        !print*,'need to check why this will happen i,k,rhi= ',i,k,rhi
       endif
      endif
     enddo
    enddo

   !*** yihsuan cloud fraction  ***
    do k=1,pver
     do i=1,ncol
        cldmfupdf(i,k) = max(cldlfupdf(i,k),cldifupdf(i,k))         !cjshiu temporarily way for the possible case of mixed phase clouds
     enddo
    enddo
   !*** yihsuan cloud fraction ***

   !cjshiu let cldtfupdf = max(CLDLFUPDF,CLDIFUPDF) need figure out the suitable way when not using temperature as proportion of cloud liquid and ice as the way now used by default CAM5 physics
!    do k=1,nlay
!     do i=1,ncol
!        cldtfupdf(i,k) = max(cldlfupdf(i,k),cldifupdf(i,k))         !cjshiu temporarily way for the possible case of mixed phase clouds
!     enddo
!    enddo

  return

end subroutine macro_triangpdf_2d

!   subroutine astG_PDF_caldelta_single( U, p, ql, qv, landfrac, snowh, a_out, Ga_out, delta_out)
   subroutine astG_PDF_caldelta_single( qs, ql, qv, a_out, delta_out)

   ! created : Yi-Chi (Jan 2015)
   ! This program is modified from subroutine astG_PDF of Park
   ! (1) add ql_in to refer for the distribution
   ! (2) output calculate dV(i.e. delta) to refer for
   ! Instead of specified width dV of triangular PDF,
   !   we use ql and qv to calculate dV and cloud fraction
   !   from the triangular shape distribution of total 
   !   water substance.
   ! Modified : Yi-Chi (April 2015)
   ! - the choice of roots from the 3-order polynomial is 
   !   corrected. dV should be < 1 and > -1 when multiple
   !   real roots exist.
   ! --------------------------------------------------------- !
   ! Park's comments                                           !
   ! --------------------------------------------------------- !
   ! Compute 'stratus fraction(a)' and Gs=(dU/da) from the     !
   ! analytical formulation of triangular PDF.                 !
   ! Here, 'dV' is the ratio of 'half-width of PDF / qs(p,T)', !
   ! so using constant 'dV' assume that width is proportional  !
   ! to the saturation specific humidity.                      !
   !    dV ~ 0.1.
   !    cldrh : RH of in-stratus( = 1 if no supersaturation)   !
   ! Note that if U > 1, Ga = 1.e10 instead of Ga = 0, that is !
   ! G is discontinuous across U = 1.  In fact, it does not    !
   ! matter whether Ga = 1.e10 or 0 at a = 1: I derived that   !
   ! they will produce the same results.                       !
   ! --------------------------------------------------------- !

   implicit none
   ! +++ Yi-Chi : March 23, 2014
   real(r8), intent(in)  :: ql                 ! liquid water
   real(r8), intent(out) :: delta_out             ! PDF width
   ! --- Yi-Chi
   !integer,  intent(in)  :: ncol
   real(r8), intent(in)  :: qs                  ! Relative humidity
!   real(r8), intent(in)  :: p                  ! Pressure [Pa]
   real(r8), intent(in)  :: qv                 ! Grid-mean water vapor specific humidity [kg/kg]
!   real(r8), intent(in)  :: landfrac           ! Land fraction
!   real(r8), intent(in)  :: snowh              ! Snow depth (liquid water equivalent)

   real(r8), intent(out) :: a_out                 ! Stratus fraction
   !real(r8), intent(out) :: Ga_out                ! dU/da

   !real(r8)              :: U
   real(r8)              :: a                     ! Stratus fraction
   !real(r8)              :: Ga                    ! dU/da

   ! Local variables
   integer :: i                                   ! Loop indexes
   real(r8) dV                                    ! Width of triangular PDF
   real(r8) cldrh                                 ! RH of stratus cloud
   !real(r8) rhmin                                 ! Critical RH
   !real(r8) rhwght

   ! +++ Yi-Chi : local variables for triangular
   real(r8)              :: rhtot, rhql, RHT
   real(r8)              :: a1,a2,a3,a4  ! coefficients for polynomial
   real(r8)              :: det ! determinant for polynomial
   real(r8)              :: delta1, delta2, delta3 ! for 3 roots and final root
   real(r8)              :: cldfrc1, cldfrc2, cldfrc3 ! for 3 roots and final root
   real(r8)              :: rhsmall,small ! small number
   real(r8)              :: cldfrc, delta
   ! ---
   ! Statement functions
   !logical land
   !land   = nint(landfrac) == 1

   ! ---------- !
   ! Parameters !
   ! ---------- !

   cldrh  = 1.0_r8


   ! ---------------- !
   ! Main computation !
   ! ---------------- !

   a_out     = 0._r8
   !Ga_out    = 0._r8


   !+++Yi-Chi: March 23, 2014
   !ql       = ql       ! liquid water
   !qs       = qv/U        ! saturated mixing ratio
   !U        = qv/qs

   rhtot    = (ql+qv)/qs  ! total relative humidity
   rhql     = ql/qs
   RHT      = rhql/(1-rhtot) ! normalized total RH

   small    = 1.e-10_r8
   rhsmall  = 1.e-3_r8
   !---Yi-Chi

   ! -------------------- !
   ! deal with exceptions !
   ! -------------------- !
   ! exception 1 : If ql is very small, cldfrac = 0.
   ! GO TO 1000
   if(ql.gt.(-1._r8*small).and.ql.lt.(1._r8*small)) then
     !print*,'Case: ql<small, cldfrc=0'
     !print*,'rhtot,rhql:',rhtot,rhql
     !cldfrc = 1._r8
     !delta  = 1._r8
     a_out     = 0._r8
     !Ga_out    = 1.e10_r8
     delta_out = 1._r8
     go to 1000
   endif
   ! ------------------------
   ! exception 2 : If U > 1 => cldfrc = 1
   if((rhtot-rhql).gt.(1._r8-rhsmall)) then
     !print*,'Case: U>1, cldfrc=1'
     !print*,'rhtot,rhql:',rhtot,rhql
     !cldfrc = 1._r8
     !delta  = 1._r8
     a_out     = 1._r8
     !Ga_out    = 1.e10_r8
     delta_out = 1._r8
     go to 1000
   endif
   ! ------------------------
   ! exception 3 : If rhtot ~0 => cldfrc = 0
   if(rhtot.lt.(0._r8+rhsmall)) then
     !print*,'Case: rhtot<rhsmall, cldfrc=0'
     !cldfrc = 1._r8
     !delta  = 1._r8
     a_out     = 0._r8
     !Ga_out    = 1.e10_r8
     delta_out = 1._r8
     go to 1000
   endif
   ! ------------------------
   ! exception 4 : If U ~0 => cldfrc = 0
   if((rhtot-rhql).lt.(0._r8+rhsmall)) then
     !print*,'Case: rh<rhsmall, cldfrc=0'
     !cldfrc = 1._r8
     !delta  = 1._r8
     a_out     = 0._r8
     !Ga_out    = 1.e10_r8
     delta_out = 1._r8
     go to 1000
   endif

   ! --------------------!
   ! For qs > qt         !
   ! --------------------!
   if(rhtot .lt. (1._r8-small)) then
     !print*,'Case:qs>qt'
     !print*,'rhtot,rhql:',rhtot,rhql
     a1 = 1._r8/6._r8
     a2 = -3._r8*(1/6._r8)
     a3 = (1._r8/6._r8)*(3._r8+6._r8*RHT)
     a4 = (1._r8/6._r8)*(-1._r8)
     !det = 18._r8*a1*a2*a3*a4 - 4._r8*(a2**3)*a4 + (a2**2)*(a3**2) - 4._r8*a1*(a3**3) - 27._r8*(a1**2)*(a4**2)
     call cubic_exact(a1,a2,a3,a4,delta1,delta2,delta3,det)
     !call cubic_exact(a1,a2,a3,a4,delta1,delta2,delta3)
       cldfrc1 = 0.5*((1-delta1)**2)
       cldfrc2 = 0.5*((1-delta2)**2)
       cldfrc3 = 0.5*((1-delta3)**2)
       delta   = -999.999
       if(cldfrc1.ge.small .and. cldfrc1.le.1._r8) then
       if(delta1.le.1) then
       delta = delta1
       endif
       endif
       if(cldfrc2.ge.small .and. cldfrc2.le.1._r8) then
       if(delta2.le.1) then
       delta = delta2
       endif
       endif
       if(cldfrc3.ge.small .and. cldfrc3.le.1._r8) then
       if(delta3.le.1) then
       delta = delta3
       endif
       endif
       !print*,'cldfrc',cldfrc1,cldfrc2,cldfrc3
       if(delta.lt.-999.0)then
         !print*,'cldfrc',cldfrc1,cldfrc2,cldfrc3
         !print*,'No reasonable cloud fraction; set cldfrc = 0'
       !  stop
         cldfrc     = 0._r8
         delta      = 1._r8
       else
         cldfrc= 0.5*((1-delta)**2)
         delta = (1-rhtot)/delta
       endif

     ! ------------------------------------- !
     ! Case : determinant = 0 => 2 real roots
     !        choose the one larger than 0, so delta > 0
     !if(det .eq. 0.0) then
     !  if(delta1 .gt.0 .and. delta3.lt.0 ) then
     !  delta = delta1
     !  elseif(delta1.lt.0 .and. delta3 .gt. 0) then
     !  delta = delta3
     !  else
     !  delta = delta1
     !  endif

     !! Case : determinant != 0
     !else
     !  if(delta1 .eq. real(delta1)) then
     !    delta = delta1
     !  elseif(delta2 .eq. real(delta2)) then
     !    delta = delta2
     !  else
     !    if(delta3 .eq. real(delta3)) then
     !    delta = delta3
     !    else
     !    print*,'no real roots!'
     !    endif
     !  endif
     !endif  ! determinant


   ! ------------------- !
   ! For qs < qt         !
   ! ------------------- !
   elseif(rhtot.gt.(1._r8+small)) then
     !print*,'qs<qt'
     !print*,'rhtot,rhql:',rhtot,rhql
     a1 = 1._r8/6._r8
     a2 = 3._r8*(1._r8/6._r8)
     a3 = (1._r8/6._r8)*(-3._r8-6._r8*RHT)
     a4 = (1._r8/6._r8)
     !det = 18._r8*a1*a2*a3*a4 - 4._r8*(a2**3)*a4 + (a2**2)*(a3**2) - 4._r8*a1*(a3**3) - 27._r8*(a1**2)*(a4**2)
     !print*,'a1,a2,a3,a4:',a1,a2,a3,a4
     !print*,'det:',det

     call cubic_exact(a1,a2,a3,a4,delta1,delta2,delta3,det)

     ! case : det = 0
     !     multiple real roots with two duplicate roots
     !if(det .eq. 0.0) then
       cldfrc1 = 1-0.5*((1+delta1)**2)
       cldfrc2 = 1-0.5*((1+delta2)**2)
       cldfrc3 = 1-0.5*((1+delta3)**2)
       delta   = -999.999
       if(cldfrc1.ge.small .and. cldfrc1.le.1._r8) then
       if((delta1.le.1).and.(delta1.ge.-1)) then
       delta = delta1
       endif
       endif
       if(cldfrc2.ge.small .and. cldfrc2.le.1._r8) then
       if((delta2.le.1).and.(delta2.ge.-1)) then
       delta = delta2
       endif
       endif
       if(cldfrc3.ge.small .and. cldfrc3.le.1._r8) then
       if((delta3.le.1).and.(delta3.ge.-1)) then
       delta = delta3
       endif
       endif
       !print*,'cldfrc',cldfrc1,cldfrc2,cldfrc3
       !print*,'delta:',delta1,delta2,delta3
       !if(delta.lt.-999.0)then
       !  !print*,'cldfrc',cldfrc1,cldfrc2,cldfrc3
       !  print*,'No reasonable cloud fraction'
       !  stop
       !endif
       if(delta.lt.-999.0)then
         !print*,'cldfrc',cldfrc1,cldfrc2,cldfrc3
         !print*,'No reasonable cloud fraction; set cldfrc = 0'
       !  stop
         cldfrc     = 0._r8
         delta      = 1._r8
       else
         cldfrc= 1-0.5*((1+delta)**2)
         delta = (1-rhtot)/delta
       endif
       !print*,'qt>qs; cldfrc=;delta=',cldfrc,delta
     !! Case : determinant != 0
     !else
     !  if(delta1 .eq. real(delta1)) then
     !    delta = delta1
     !  elseif(delta2 .eq. real(delta2)) then
     !    delta = delta2
     !  else
     !    if(delta3 .eq. real(delta3)) then
     !    delta = delta3
     !    else
     !    print*,'no real roots!'
     !    endif
     !  endif
     !endif  ! determinant

     !cldfrc= 1-0.5*((1+delta)**2)
     !delta = (1-rhtot)/delta


   ! ------------------- !
   ! For rhtot = 1       !
   ! ------------------- !
   else
     cldfrc= 0.5_r8
     delta = 6._r8*rhql
   endif

   a_out     = cldfrc
   delta_out = delta

   !print*,'final choice for cldfrc,delta(single)',cldfrc,delta

   ! +++ Yi-Chi :
   ! one condition:
   !  when environment is dry and width is too large (i.e. condensate exists)
   !        call cloud fraction = 0
   if( (cldfrc .lt. 0.1_r8).and.((rhtot-delta).lt.0._r8) ) then
     a_out    = 0._r8
     cldfrc   = 0._r8
     delta_out= 1._r8 -0.8_r8 ! dV = 1-rhmin
   !  print*,'too dry. Set cldfrc,delta as',cldfrc,delta
   endif
   !-----------
!1000     print*,'final: a_out,delta_out:',a_out,delta_out
1000   return
   end subroutine astG_PDF_caldelta_single

   ! ----------------- !
   ! End of subroutine !
   ! ----------------- !


   subroutine cubic_exact( a, b, c, d, x1, x2, x3, det)
   ! ----------------- !
   ! This function deals with ax^3 + bx^2 + cx + d
   ! And output the solution based on the exact formula.
   ! ----------------- !

   real(r8), intent(in)  :: a,b,c,d     ! Land fraction

   real(r8), intent(out) :: x1, x2, x3          ! Stratus fraction
   !real(r8), intent(out) :: x1
   !real(r8), intent(out) :: x2,x3

   real(r8), intent(out) :: det                     ! determinant
   real(r8)              :: p                     ! Pressure [Pa]
   real(r8)              :: delta0, delta1,deltasqrt        !
   real(r8)              :: small
   !real(r8)              :: u1,u2,u3              ! determinant
   real(r8)              :: u1
   complex      :: u2, u3
   !complex      :: u1,u2, u3
   complex      :: coroot,corootin,deltasqrtcmplx                     ! Pressure [Pa]
   complex      :: x1tmp, x2tmp, x3tmp

   small = 1.0e-8_r8
   det   = 18._r8*a*b*c*d - 4._r8*(b**3)*d + (b**2)*(c**2) - 4._r8*a*(c**3) - 27._r8*(a**2)*(d**2)

   delta0 = b**2 -3._r8*a*c
   delta1 = 2._r8*(b**3) - 9._r8*a*b*c + 27._r8*(a**2)*d

   !print*,'a,b,c,d:',a,b,c,d
   !print*,'delta1,delta0',delta1,delta0
   !print*,'det(cubic_exact)',det
   !print*,'small',small
   if ((det.lt.small) .and. (det.gt.-1*small)) then
   !if ((det.lt.small) .and. (det.gt.-1*small) .and.  &
   !    (delta0.lt.small) .and. (delta0.gt.-1._r8*small)) then
   ! case : D=0
   ! delta0 = 0
   !  print*,'det~0'
   !  print*,'delta0',delta0
     if((delta0.lt.small) .and. (delta0.gt.-1._r8*small)) then
       x1 = -b/(3._r8*a)
       x2 = x1
       x3 = x1
     else ! delta0 !=0
       x1 = (9._r8*a*d-b*c)/(2._r8*delta0)
       x2 = x1
       x3 = (4._r8*a*b*c-9._r8*(a**2)*d-(b**3))/(a*delta0)
     endif
   det = 0._r8
   !  print*,'x1,x2,x3(det=0)',x1,x2,x3
   else
     if(det.gt.small) then
      ! print*,'det>0 : distinct real roots.'
      ! print*,'rhtot,rhql',rhtot,rhql
      ! print*,'qv,ql,qs',qv,ql,qs
     endif
     ! case : D>0 or D<0
     u1 = 1._r8
     u2 = cmplx(-1._r8/2._r8,sqrt(3._r8)/2._r8)
     u3 = cmplx(-1._r8/2._r8,-sqrt(3._r8)/2._r8)
     !u2 = (-1+sqrt(-1._r8)*sqrt(3._r8))/2
     !u3 = (-1-sqrt(-1._r8)*sqrt(3._r8))/2

     deltasqrt = delta1**2-4._r8*(delta0**3)
     !if(deltasqrt.lt.0._r8) then
       deltasqrtcmplx = cmplx(deltasqrt,0.0)
       !print*,'sqrt(deltasqrt)', sqrt(deltasqrt)
       call ccbrt(2,deltasqrtcmplx,corootin)
       corootin = 0.5_r8*(delta1 + corootin)
       !corootin  = 0.5*(delta1+sqrt(deltasqrt))
       call ccbrt(3,corootin,coroot)
     !else
     !  coroot  = ((1._r8/2._r8)*(delta1 + sqrt(deltasqrt)))**(1/3)
     !endif
     !print*,'coroot:',coroot
     !print*,'a:',a
     !print*,'u1:',u1
     !print*,'u2:',u2
     !print*,'u3:',u3
     x1tmp = (-1._r8/3._r8)*(1._r8/a)*(b+u1*coroot+delta0/(u1*coroot))
     x2tmp = (-1._r8/3._r8)*(1._r8/a)*(b+u2*coroot+delta0/(u2*coroot))
     x3tmp = (-1._r8/3._r8)*(1._r8/a)*(b+u3*coroot+delta0/(u3*coroot))
     !print*,'u1*coroot',u1*coroot
     !print*,'delta0/(u1*coroot)',delta0/(u1*coroot)
     !print*,'x1tmp,x2tmp,x3tmp:',x1tmp,x2tmp,x3tmp

     if((aimag(x1tmp).ge.small).and.(aimag(x1tmp).le.-1._r8*small))then
     x1 =  -999.999
     else
     x1 = real(x1tmp)
     endif
     if((aimag(x2tmp).ge.small).and.(aimag(x2tmp).le.-1._r8*small))then
     x2 =  -999.999
     else
     x2 = real(x2tmp)
     endif
     if((aimag(x3tmp).ge.small).and.(aimag(x3tmp).le.-1._r8*small))then
     x3 =  -999.999
     else
     x3 = real(x3tmp)
     endif

   endif


   return
   end subroutine cubic_exact
 
   ! ----------------- !
   ! End of subroutine !
   ! ----------------- !
  subroutine ccbrt(sqrtnum,datain,datacube)
   ! ----------------- !
   ! This function deals with cubic root for complex number
   ! And output the solution based on the exact formula.
   ! ----------------- !

   integer, intent(in) :: sqrtnum
   complex, intent(in) :: datain
   complex, intent(out) :: datacube
   real(r8) :: datareal, dataimg, radius          ! Stratus fraction
   real(r8) :: radius_cube, theta_cube
   !real(r8), intent(out) :: x1
   !real(r8), intent(out) :: x2,x3
   datareal = real(datain)
   dataimg  = imag(datain)
   radius   = sqrt(datareal**2 + dataimg**2)

   radius_cube = radius**(1._r8/real(sqrtnum))
   ! print*,'datareal,dataimg:',datareal,dataimg
   ! print*,'radius,radius_cube:',radius,radius_cube
   ! ATAN2(Y,X) computes the arctangent of the complex number X + i Y.
   theta_cube  = atan2(dataimg,datareal)/real(sqrtnum)
   ! print*,'theta_cube,radius_cube:',theta_cube,radius_cube
   datacube    = cmplx(cos(theta_cube), sin(theta_cube))
   datacube    = radius_cube*datacube
   ! print*,'datacube:',datacube
   end subroutine ccbrt

   ! ----------------- !
   ! End of Subroutine !
   ! ----------------- !
end module cldfrc_triang_macro

