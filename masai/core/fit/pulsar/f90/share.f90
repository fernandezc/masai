!===============================================================================
! pulsar.share
!===============================================================================
! Copyright (C) 2006-2013 Christian Fernandez, Jean-Paul Amoureux
! JPA - Unite de Catalyse et Chimie du Solide, Lille, France.
! CF  - Laboratoire Catalyse et Spectrochimie, Caen, France.  
!       christian.fernandez@ensicaen.fr
! This software is governed by the CeCILL-B license under French law 
! and abiding by the rules of distribution of free software.  
! You can  use, modify and/ or redistribute the software under 
! the terms of the CeCILL-B license as circulated by CEA, CNRS and INRIA
! at the following URL "http://www.cecill.info".
! See Licence.txt in the main masai source directory
!===============================================================================


!------------------------
MODULE share
!------------------------
! other shared parameters
!------------------------
  use parameters

  integer*8 :: nt=-1
  integer*8 :: ng=-1
  logical*2 :: staticsample=.false.              ! (replace ispeed in the original version)
  real*8    :: dg=eps

  real*8    :: de(-4:4,-4:4)=eps                 ! wigner matrix
  real*8    :: dk(2,0:4,-2:2)=eps                !

  real*8    :: spins=eps                         ! spin i
  real*8    :: spini=eps                         ! spin s

  real*8    :: vls=eps                           ! larmor of spin s
  real*8    :: cqs=eps
  real*8    :: abundance=1.0d0
  real*8    :: etaqs=eps
  real*8    :: isos=eps
  real*8    :: csas=eps
  real*8    :: etacs=eps
  real*8    :: t2s=eps
  real*8    :: qq=eps                            ! effect of quality factor
  real*8    :: wm(4)=eps
  real*8    :: fst=eps
  integer*8 :: icoupled=-1                       ! if i exist
  real*8    :: vli=eps                           ! larmor of spin i
  real*8    :: cqi=eps
  real*8    :: etaqi=eps
  real*8    :: isoi=eps
  real*8    :: csai=eps
  real*8    :: etaci=eps
  real*8    :: t2i=eps
  real*8    :: jsi=eps                           ! j coupling between s and i
  real*8    :: dsi=eps                           ! dipolar coupling i-s
  integer*8 :: kcoupled=-1                       ! number of coupled k nucleus to s
  real*8    :: coefps=eps
  real*8    :: coefpi=eps
  real*8    :: coefss=eps
  real*8    :: coefsi=eps
  real*8    :: isost=eps

  complex*8 :: r2mdsi(0:2)                       ! r2m component for is dipolar tensor
  complex*8 :: r2mqi(0:2)
  complex*8 :: r2mcs(0:2)
  complex*8 :: r2mci(0:2)

  integer*8 :: ns=-1
  integer*8 :: ni=-1
  integer*8 :: nsi=-1
  real*8    :: fnexp=eps

  real*8 :: vvi(2)
  real*8 :: vvs(2)
  real*8 :: vvii(2)
  real*8 :: vvis(2)
  complex*8 :: ats(-4:4)                   ! arts and aits
  complex*8 :: aqs(-4:4)                  ! arqs(10,-4:4) and aiqs(10,-4:4)
  complex*8 :: acs(-4:4)                  ! arcs and aics
  complex*8 :: ads(-4:4)                     ! ads(4,:) correspondant  dipolaire i-s
  complex*8 :: aqi(-4:4)                     ! arqi and aiqi
  complex*8 :: aci(-4:4)                     ! arci and aici
  complex*8 :: aqpi(0:4)
  complex*8 :: aqps(0:4)
  complex*8 :: adcps(0:4)
  complex*8 :: aqsi(2,0:4)
  complex*8 :: aqss(2,0:4)
  complex*8 :: acspi(0:4)
  complex*8 :: adipolis(0:4)

  real*8 :: zl(0:4)=(/1.0d0,2.0d0,2.0d0,2.0d0,2.0d0/)

  real*8    :: qsangle(3)
  real*8    :: csangle(3)
  real*8    :: qiangle(3)
  real*8    :: ciangle(3)
  real*8    :: dsiangle(3)                   ! dipolar orientation

  real*8, dimension(:), allocatable :: spink              ! non -excited and non-observed nuclei
  real*8, dimension(:), allocatable :: vlk                ! larmor of spins k
  real*8, dimension(:), allocatable :: jsk                ! j coupling between s-k
  real*8, dimension(:), allocatable :: dsk                ! dipolar coupling s-k
  real*8, dimension(:,:), allocatable :: dskangle         ! dipolar orientation of s-k

  real*8, dimension(:), allocatable :: rfpowers 
  real*8, dimension(:), allocatable :: rfoffsets 
  real*8, dimension(:), allocatable :: rfphases 
  real*8, dimension(:), allocatable :: rfpoweri 
  real*8, dimension(:), allocatable :: rfoffseti 
  real*8, dimension(:), allocatable :: rfphasei 
  real*8, dimension(:), allocatable :: rflength 
  real*8, dimension(:), allocatable :: rfdelay 
  integer*8, dimension(:), allocatable :: decouple 

  integer*8, dimension(:), allocatable :: itour 
  integer*8, dimension(:), allocatable :: ntp 
  real*8, dimension(:), allocatable :: deltatp 
  integer*8, dimension(:), allocatable :: ntpadd 
  real*8, dimension(:), allocatable :: deltatpadd 
  real*8, dimension(:), allocatable :: t2ss 
  real*8, dimension(:), allocatable :: t2ii 
  real*8, dimension(:), allocatable :: t2si 

  complex*8, dimension(:,:), allocatable :: rott             ! rortt and roitt
  complex*8, dimension(:,:,:,:), allocatable :: qt           ! qrt and qit(0:40,0:80,11,-20:20)

  real*8, dimension(:), allocatable :: u1s 
  real*8, dimension(:), allocatable :: u2s 

  complex*8, dimension(:,:), allocatable :: r2mdsk           ! tensor components ks
  integer*8, dimension(:), allocatable :: nk 

  real*8, dimension(:), allocatable :: dipk 
  integer*8, dimension(:), allocatable :: mk 

  real*8, dimension(:,:,:), allocatable :: wr 
  complex*8, dimension(:,:,:,:), allocatable :: amp 
  real*8                 :: dnorm=eps


  complex*8, dimension(:,:), allocatable :: adks             ! ards and aids   correspondant au dipolaire k-s
  complex*8, dimension(:,:), allocatable :: a                ! ab(11,4),bb(11,4)
  complex*8, dimension(:,:), allocatable :: p                ! pr(11,-20:20) and pi(11,-20:20)
  complex*8, dimension(:,:), allocatable :: u                ! tr and ti
  complex*8, dimension(:,:), allocatable :: ui               ! tri(100,100),tii(100,100)
  complex*8, dimension(:,:), allocatable :: uf               ! trf(100,100),tif(100,100)
  complex*8, dimension(:,:), allocatable :: usi              ! trsi(100,100),tisi(100,100)
  complex*8, dimension(:,:), allocatable :: usf              ! trsf(100,100),tisf(100,100)
  complex*8, dimension(:,:), allocatable :: uint             ! trint(100,100),tiint(100,100)
  complex*8, dimension(:,:), allocatable :: ro               ! ror and roi
  complex*8, dimension(:,:), allocatable :: rot              ! rotr and roti
  complex*8, dimension(:,:), allocatable :: y                ! yr and yi
  complex*8, dimension(:,:), allocatable :: htst             ! htsrt(100,100) and htsit(100,100)
  complex*8, dimension(:,:), allocatable :: z                ! zr and zi
  real*8, dimension(:), allocatable :: vp 
  complex*8, dimension(:,:), allocatable :: expvp            ! exp(ivp)
  real*8, dimension(:,:), allocatable :: htsr             ! htsr(100,100)
  complex*8, dimension(:,:), allocatable :: sti              ! stri(11,11),stii(11,11)
  complex*8, dimension(:,:), allocatable :: rosto            !
  complex*8, dimension(:,:), allocatable :: q                ! qr(11,-20:20) and qi(11,-20:20)
  real*8, dimension(:), allocatable :: phasei 

  COMPLEX*8, dimension(:,:,:), allocatable :: roctp 
  integer*8 :: nctp=-1

  integer*8 :: STOPPED=-1

  complex*8 :: sumc
  complex*8 :: sums


END MODULE share

!-----------------------
subroutine reset_share()
!-----------------------

use share

      nt=-1
      ng=-1
      staticsample=.false.              ! (replace ispeed in the original version)
      dg=eps
      de=eps                            ! wigner matrix
      dk=eps                            !
      spins=eps                         ! spin i
      spini=eps                         ! spin s

      vls=eps                           ! larmor of spin s
      cqs=eps
      sigmaqs = eps
      abundance=1.0d0
      etaqs=eps
      isos=eps
      csas=eps
      etacs=eps
      t2s=eps
      qq=eps                            ! effect of quality factor
      wm=eps
      fst=eps
      icoupled=-1                       ! if i exist
      vli=eps                           ! larmor of spin i
      cqi=eps
      sigmaqi = eps
      etaqi=eps
      isoi=eps
      csai=eps
      etaci=eps
      t2i=eps
      jsi=eps                           ! j coupling between s and i
      dsi=eps                           ! dipolar coupling i-s
      kcoupled=-1                       ! number of coupled k nucleus to s
      coefps=eps
      coefpi=eps
      coefss=eps
      coefsi=eps
      isost=eps

      sumc = cmplx(0.0d0,0.0d0)
      sums = cmplx(0.0d0,0.0d0)



      r2mdsi=cmplx(0.0d0,0.0d0)                       ! r2m component for is dipolar tensor
      r2mqi=cmplx(0.0d0,0.0d0)
      r2mcs=cmplx(0.0d0,0.0d0)
      r2mci=cmplx(0.0d0,0.0d0)

      ns=-1
      ni=-1
      nsi=-1
      fnexp=eps

      vvi=0.0d0
      vvs=0.0d0
      vvii=0.0d0
      vvis=0.0D0
      ats=cmplx(0.0d0,0.0d0)                   ! arts and aits
      aqs=cmplx(0.0d0,0.0d0)                  ! arqs(10,-4:4) and aiqs(10,-4:4)
      acs=cmplx(0.0d0,0.0d0)                  ! arcs and aics
      ads=cmplx(0.0d0,0.0d0)                     ! ads(4,:) correspondant  dipolaire i-s
      aqi=cmplx(0.0d0,0.0d0)                     ! arqi and aiqi
      aci=cmplx(0.0d0,0.0d0)                     ! arci and aici
      aqpi=cmplx(0.0d0,0.0d0)
      aqps=cmplx(0.0d0,0.0d0)
      adcps=cmplx(0.0d0,0.0d0)
      aqsi=cmplx(0.0d0,0.0d0)
      aqss=cmplx(0.0d0,0.0d0)
      acspi=cmplx(0.0d0,0.0d0)
      adipolis=cmplx(0.0d0,0.0d0)

      zl(0:4)=(/1.0d0,2.0d0,2.0d0,2.0d0,2.0d0/)

      qsangle=0.0d0
      csangle=0.0d0
      qiangle=0.0d0
      ciangle=0.0d0
      dsiangle=0.0d0

      if (allocated(spink)) deallocate (spink)
      if (allocated(vlk)) deallocate (vlk)
      if (allocated(jsk)) deallocate (jsk)
      if (allocated(dsk)) deallocate (dsk)
      if (allocated(dskangle)) deallocate (dskangle)

      if (allocated(rfpowers)) deallocate (rfpowers)
      if (allocated(rfoffsets)) deallocate (rfoffsets)
      if (allocated(rfphases)) deallocate (rfphases)
      if (allocated(rfpoweri)) deallocate (rfpoweri)
      if (allocated(rfoffseti)) deallocate (rfoffseti)
      if (allocated(rfphasei)) deallocate (rfphasei)
      if (allocated(rflength)) deallocate (rflength)
      if (allocated(rfdelay)) deallocate (rfdelay)
      if (allocated(decouple)) deallocate (decouple)

      if (allocated(itour)) deallocate (itour)
      if (allocated(ntp)) deallocate (ntp)
      if (allocated(deltatp)) deallocate (deltatp)
      if (allocated(ntpadd)) deallocate (ntpadd)
      if (allocated(deltatpadd)) deallocate (deltatpadd)
      if (allocated(t2ss)) deallocate (t2ss)
      if (allocated(t2ii)) deallocate (t2ii)
      if (allocated(t2si)) deallocate (t2si)

      if (allocated(rott)) deallocate(rott)
      if (allocated(qt)) deallocate (qt)

      if (allocated(u1s)) deallocate (u1s)
      if (allocated(u2s)) deallocate (u2s)

      if (allocated(r2mdsk)) deallocate (r2mdsk)          ! tensor components ks
      if (allocated(nk)) deallocate (nk)

      if (allocated(dipk)) deallocate(dipk)
      if (allocated(mk)) deallocate(mk)

      if (allocated(wr)) deallocate (wr)
      if (allocated(amp)) deallocate (amp)
      dnorm=0.0d0


      if (allocated(adks)) deallocate (adks)
      if (allocated(a)) deallocate (a)
      if (allocated(p)) deallocate (p)
      if (allocated(u)) deallocate (u)
      if (allocated(ui)) deallocate (ui)
      if (allocated(uf)) deallocate (uf)
      if (allocated(usi)) deallocate (usi)
      if (allocated(usf)) deallocate (usf)
      if (allocated(uint)) deallocate (uint)
      if (allocated(ro)) deallocate (ro)
      if (allocated(rot)) deallocate (rot)
      if (allocated(y)) deallocate (y)
      if (allocated(htst)) deallocate (htst)
      if (allocated(z)) deallocate (z)
      if (allocated(vp)) deallocate (vp)
      if (allocated( expvp)) deallocate (expvp)
      if (allocated(htsr)) deallocate (htsr)
      if (allocated(sti)) deallocate (sti)
      if (allocated(rosto)) deallocate (rosto)
      if (allocated(q)) deallocate (q)
      if (allocated(phasei)) deallocate (phasei)
      if (allocated(roctp)) deallocate(roctp)
      nctp=-1

END subroutine

