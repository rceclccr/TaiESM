subroutine macro_uniformpdf_2d(ncol,  play, tlay, qcwat, qclwat, iccwat, &
           deltalpdf, cldlfupdf, deltaipdf, cldifupdf, deltatpdf, cldtfupdf, cldmfupdf, qscwatupdf, qsciceupdf)

   use shr_kind_mod,     only: r8 => shr_kind_r8
   use ppgrid,           only: pcols, pver, pverp
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
!cjshiu  real(r8), parameter :: qlmin = 2.e-8_r8                 ! minimum conc. of cloud liquid mass (as qsmall of NCEP value of 1.e-10 or 2.e-12 is used)
!cjshiu  real(r8), parameter :: qimin = 2.e-18_r8                 ! minimum conc. of cloud ice mass 
!cjshiu  real(r8), parameter :: qimin = 2.e-12_r8                 ! minimum conc. of cloud ice mass 
  real(r8), parameter :: qimin = 2.e-10_r8                 ! minimum conc. of cloud ice mass 
!cjshiu  real(r8), parameter :: qimin = 2.e-8_r8                 ! minimum conc. of cloud ice mass 
!cjshiu  real(r8), parameter :: qimin = 1.e-8_r8                 ! minimum conc. of cloud ice mass 
!cjshiu  real(r8), parameter :: qimin = 1.e-7_r8                 ! minimum conc. of cloud ice mass 
        
  real(r8), parameter :: qtmin = qlmin + qimin             ! minimum conc. of cloud ice+liquid mass 

!cjshiu  real(r8), parameter :: sup = 1.0_r8
!cjshiu  real(r8), parameter :: sup = 1.2_r8
!cjshiu  real(r8), parameter :: sup = 1.4_r8
!cjshiu  real(r8), parameter :: sup = 1.6_r8
  real(r8), parameter :: sup = 1.0_r8
      
!cjshiu  real(r8), parameter :: rhcrit = 0.5_r8
  real(r8), parameter :: rhcrit = 0.8_r8
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

        deltasqrt = ( sqrt(qtbar) + sqrt((qsc-qvbar)) )
        deltatpdf(i,k) = deltasqrt * deltasqrt
        cldtfupdf(i,k) = 1._r8/(2._r8*deltatpdf(i,k))*(qtbar+qvbar+deltatpdf(i,k)-qsc)

      ! if cloud water smaller than threshold value,
      ! use relative humidity to determine cloud fraction
      else

        ! grid-scale relative humidity
        rh = qvbar/qscwat(i,k)
!cjshiu          print*,'cjshiu in macro_uniformPDF.F90 ice phase qvbar=  qscwat= i k ', qvbar,qscwat(i,k),i,k

        ! liquid phase cloud
        if ( rh .GE. rhcrit .and. tlay(i,k) .GT. 273.15_r8) then
          cldtfupdf(i,k) = 1._r8 - ((1._r8-rh)/(1._r8-rhcrit))**0.5

        ! ice phase cloud
        elseif ( rh .GE. rhcrit .and. tlay(i,k) .LE. 273.15_r8) then
!cjshiu          print*,'cjshiu in macro_uniformPDF.F90 ice phase rh= ',rh
     !cjshiu add the following to fix when RH > 1.1 for clouf ice e.g. rh = 1.278     
          rh = min(rh,1.1_r8) 
!cjshiu          cldtfupdf(i,k) = 1._r8 - ((1.1_r8-rh)/(1.1_r8-rhcrit))**0.5
          cldtfupdf(i,k) = 1._r8 - ((1.1_r8-rh)/(1.1_r8-rhcrit))**0.5
!cjshiu check          cldtfupdf(i,k) =  ((1.1_r8-rh)/(1.1_r8-rhcrit))**2._r8
!cjshiu need to check the above formula

        ! rh < rh_critical, no cloud
        else
          cldtfupdf(i,k) = 0._r8

        endif

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
        !cjshiu there is cases that qsc < qvbar (set to qvbar slightly lower than qsc to avoid floating invalid)
          if (qsc .LE. qvbar) then
!cjshiu             qvbar = 0.9999_r8*qsc
!cjshiu            print*,'cjshiu i, k, qsc, qvbar, qlbar=',i,k,qsc,qvbar,qlbar
            qvbar = qsc
          endif
        !print*,'qlbar= qvbar= qsc= ',qlbar,qvbar,qsc
      
      rh = qvbar/qsc

      if (qlbar .GT. qlmin) then
        deltasqrt = ( sqrt(qlbar) + sqrt((qsc-qvbar)) )
        !print*,'delatasqrt= ',deltasqrt
        deltalpdf(i,k) = deltasqrt * deltasqrt
        !cjshiu may need to add zero check for delta
        cldlfupdf(i,k) = 1._r8/(2._r8*deltalpdf(i,k))*(qlbar+qvbar+deltalpdf(i,k)-qsc)
        !print*,'i=  k= deltalpdf(i,k)= cldlfupdf(i,k) = ',i,k,deltalpdf(i,k),cldlfupdf(i,k)
      elseif (qlbar .LE. qlmin) then
       if (rh .GT. rhcrit) then
        deltalpdf(i,k) = qsc*(1._r8-rhcrit)
        cldlfupdf(i,k) = 1._r8 - ((1._r8-rh)/(1._r8-rhcrit))**0.5
        elseif (rh .LE. rhcrit) then
        !cjshiu should test for the rhminl and rhminh as those of Park's usage, first assume to be 0.8
        deltalpdf(i,k) = qsc*(1._r8-rhcrit)      !cjshiu which value can be used??
        !cjshiu may need to consider assume a very small value for cldfraction, first set 0
        cldlfupdf(i,k) = 0._r8
        else
        print*,'need to check why this will happen i,k,rh= ',i,k,rh
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
!cjshiu        rhi = (qvbar+qibar)/qscice(i,k)
!cjshiu check
!cjshiu        print*,'cjshiu rhi in uniform pdf i, k, rhi=',i,k,rhi

     if (qibar .GT. qimin) then
!cjshiu      if (qibar .GT. qimin .AND. tlay(i,k) .LE. 273.15_r8) then
        deltasqrt = ( sqrt(qibar) + sqrt((qsc-qvbar)) )
        !print*,'delatasqrt= ',deltasqrt
        deltaipdf(i,k) = deltasqrt * deltasqrt
        !cjshiu may need to add zero check for delta
        cldifupdf(i,k) = 1._r8/(2._r8*deltaipdf(i,k))*(qibar+qvbar+deltaipdf(i,k)-qsc)
        !print*,'i=  k= deltaipdf(i,k)= cldifupdf(i,k) = ',i,k,deltaipdf(i,k),cldifupdf(i,k)
      elseif (qibar .LE. qimin) then
!cjshiu      elseif (qibar .LE. qimin .AND. tlay(i,k) .LE. 273.15_r8) then
       if (rhi .GT. rhcrit) then
        deltaipdf(i,k) = qsc*(1._r8-rhcrit)
!cjshiu        cldifupdf(i,k) = 1._r8 - ((1.1_r8-rhi)/(1.1_r8-rhcrit))**0.5
        cldifupdf(i,k) = 1._r8 - ((1._r8-rhi)/(1._r8-rhcrit))**0.5 !cjshiu if add qibar to rhi, rhi > 1.0 can occur and cause?
        elseif (rhi .LE. rhcrit) then
        !cjshiu should test for the rhminl and rhminh as those of Park's uasage, first assume to be 0.8
        deltaipdf(i,k) = qsc*(1._r8-rhcrit)  !cjshiu how about assuming to be zero
        !cjshiu may need to consider assume a very small value for cldfraction, first set 0
        cldifupdf(i,k) = 0._r8
        else
        print*,'need to check why this will happen i,k,rhi= ',i,k,rhi
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

end subroutine macro_uniformpdf_2d
