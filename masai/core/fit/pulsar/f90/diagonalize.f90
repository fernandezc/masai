!===============================================================================
! pulsar.diagonalize
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


module diagonalize
!-----------------

  real*8, dimension(:), allocatable, save :: e 
  real*8, dimension(:,:), allocatable, save :: tau 
  real*8, dimension(:,:), allocatable, save :: hr 
  real*8, dimension(:,:), allocatable, save :: hi 
  real*8, dimension(:,:), allocatable, save :: zr 
  real*8, dimension(:,:), allocatable, save :: zi 


  contains

    subroutine dealloc()
    !-------------------

      if (allocated(e)) deallocate (e)
      if (allocated(tau)) deallocate (tau)
      if (allocated(hr)) deallocate (hr)
      if (allocated(hi)) deallocate (hi)
      if (allocated(zr)) deallocate (zr)
      if (allocated(zi)) deallocate (zi)

    end subroutine dealloc

    subroutine diag(n,hts,vp,z)
    !--------------------------
      integer*8 :: n
      real*8    :: vp(n)
      complex*8 :: z(n,n)
      complex*8 :: hts(n,n)

      integer*8 :: i

      if (.not.allocated(e)) allocate (e(n))
      if (.not.allocated(tau)) allocate (tau(2,n))
      if (.not.allocated(hr)) allocate (hr(n,n))
      if (.not.allocated(hi)) allocate (hi(n,n))
      if (.not.allocated(zr)) allocate (zr(n,n))
      if (.not.allocated(zi)) allocate (zi(n,n))

      vp=0.0d0
      zr=0.0d0
      zi=0.0d0
      do i=1,n
        zr(i,i)=1.0d0
      end do
      hr=real(hts)
      hi=aimag(hts)
                                                 ! todo: make the SUBROUTINEs handling complex arrays
      call htridi(n,hr,hi,vp,e,tau)
      call tql2(n,vp,e,zr)
      call htribk(n,hr,hi,tau,zr,zi)

      hts=cmplx(hr,hi)
      z=cmplx(zr,zi)

    end subroutine diag

    subroutine htridi(n,htsr,htsi,vp,e,tau)
    !--------------------------------------
      integer*8 :: n
      real*8    :: htsr(n,n)
      real*8    :: htsi(n,n)
      real*8    :: vp(n)
      real*8    :: e(n)
      real*8    :: tau(2,n)
      real*8 eps
      parameter(eps=1.0d-30)

      integer*8 :: i,j,ii,l,k,jp1
      real*8    :: h,g,scal,f,si,gi,hh,fi

      tau(1,n)=1.
      tau(2,n)=0.
      do 100 i=1,n
    100 vp(i)=htsr(i,i)
      do 300 ii=1,n
        i=n+1-ii
        l=i-1
        h=0.
        scal=0.
        if(l.lt.1) goto 130
        do 120 k=1,l
    120 scal=scal+abs(htsr(i,k))+abs(htsi(i,k))
        if(abs(scal).gt.eps) goto 140
        tau(1,l)=1.
        tau(2,l)=0.
    130 e(i)=0.
        goto 290
    140 do 150 k=1,l
            htsr(i,k)=htsr(i,k)/scal
            htsi(i,k)=htsi(i,k)/scal
    150 h=h+htsr(i,k)*htsr(i,k)+htsi(i,k)*htsi(i,k)
        g=sqrt(h)
        e(i)=scal*g
        f=cabs(cmplx(htsr(i,l),htsi(i,l)))
        if(abs(f).lt.eps) goto 160
        tau(1,l)=(htsi(i,l)*tau(2,i)-htsr(i,l)*tau(1,i))/f
        si=(htsr(i,l)*tau(2,i)+htsi(i,l)*tau(1,i))/f
        h=h+f*g
        g=1.+g/f
        htsr(i,l)=g*htsr(i,l)
        htsi(i,l)=g*htsi(i,l)
        if(l.eq.1) goto 270
        goto 170
    160 tau(1,l)=-tau(1,i)
        si=tau(2,i)
        htsr(i,l)=g
    170 f=0.
        do 240 j=1,l
        g=0.
        gi=0.
        do 180 k=1,j
        g=g+htsr(j,k)*htsr(i,k)+htsi(j,k)*htsi(i,k)
    180 gi=gi-htsr(j,k)*htsi(i,k)+htsi(j,k)*htsr(i,k)
        jp1=j+1
        if(l.lt.jp1) goto 220
        do 200 k=jp1,l
        g=g+htsr(k,j)*htsr(i,k)-htsi(k,j)*htsi(i,k)
    200 gi=gi-htsr(k,j)*htsi(i,k)-htsi(k,j)*htsr(i,k)
    220 e(j)=g/h
        tau(2,j)=gi/h
    240 f=f+e(j)*htsr(i,j)-tau(2,j)*htsi(i,j)
        hh=f/(h+h)
        do 260 j=1,l
        f=htsr(i,j)
        g=e(j)-hh*f
        e(j)=g
        fi=-htsi(i,j)
        gi=tau(2,j)-hh*fi
        tau(2,j)=-gi
        do 260 k=1,j
        htsr(j,k)=htsr(j,k)-f*e(k)-g*htsr(i,k)+fi*tau(2,k)+gi*htsi(i,k)
    260 htsi(j,k)=htsi(j,k)-f*tau(2,k)-g*htsi(i,k)-fi*e(k)-gi*htsr(i,k)
    270 do 280 k=1,l
        htsr(i,k)=scal*htsr(i,k)
    280 htsi(i,k)=scal*htsi(i,k)
        tau(2,l)=-si
    290 hh=vp(i)
        vp(i)=htsr(i,i)
        htsr(i,i)=hh
    300 htsi(i,i)=scal*sqrt(h)

    end subroutine htridi

    subroutine tql2(n,vp,e,z)
    !------------------------
      integer*8 :: n
      real*8    :: vp(n)
      real*8    :: e(n)
      real*8    :: z(n,n)

      integer*8 :: i,j,l,m, l1, mml, ii, k
      real*8    :: zachep,f,b,h,g,r,c,s,p

      zachep=2.**(-26)
      do 100 i=2,n
    100  e(i-1)=e(i)
      f=0.
      b=0.
      e(n)=0.
      do 240 l=1,n
       j=0
       h=zachep*(abs(vp(l))+abs(e(l)))
       if(b.lt.h) b=h
        do 110 m=l,n
    110     if(abs(e(m)).le.b) goto 120
    120     if(m.eq.l) goto 220
    130     if(j.eq.30) print *, "There is a big problem in 'htridi'"
            j=j+1
            l1=l+1
            g=vp(l)
            p=(vp(l1)-g)/(2.*e(l))
            r=sqrt(p*p+1.)
            vp(l)=e(l)/(p+sign(r,p))
            h=g-vp(l)
            do 140 i=l1,n
    140     vp(i)=vp(i)-h
            f=f+h
            p=vp(m)
            c=1.
            s=0.
            mml=m-l
            do 200 ii=1,mml
               i=m-ii
               g=c*e(i)
               h=c*p
               if(abs(p).lt.abs(e(i))) goto 150
               c=e(i)/p
               r=sqrt(c*c+1.)
               e(i+1)=s*p*r
               s=c/r
               c=1./r
               goto 160
    150        c=p/e(i)
               r=sqrt(c*c+1.)
               e(i+1)=s*e(i)*r
               s=1./r
               c=c*s
    160        p=c*vp(i)-s*g
               vp(i+1)=h+s*(c*g+s*vp(i))
               do 180 k=1,n
                  h=z(k,i+1)
                  z(k,i+1)=s*z(k,i)+c*h
    180           z(k,i)=c*z(k,i)-s*h
    200     continue
            e(l)=s*p
            vp(l)=c*p
            if(abs(e(l)).gt.b) goto 130
    220     vp(l)=vp(l)+f
    240  continue
         do 300 ii=2,n
            i=ii-1
            k=i
            p=vp(i)
            do 260 j=ii,n
               if(vp(j).ge.p) goto 260
               k=j
               p=vp(j)
    260     continue
            if(k.eq.i) goto 300
            vp(k)=vp(i)
            vp(i)=p
            do 280 j=1,n
               p=z(j,i)
               z(j,i)=z(j,k)
               z(j,k)=p
    280     continue
    300  continue

    end subroutine tql2

    subroutine htribk(n,htsr,htsi,tau,zr,zi)
    !---------------------------------------
      integer*8 :: n
      real*8    :: htsr(n,n)
      real*8    :: htsi(n,n)
      real*8    :: tau(2,n)
      real*8    :: zr(n,n)
      real*8    :: zi(n,n)
      real*8 eps
      parameter(eps=1.0d-30)

      integer*8 :: k,i,j,l
      real*8    :: h,s,si

      do 50 k=1,n
          do 50 j=1,n
          zi(k,j)=-zr(k,j)*tau(2,k)
  50      zr(k,j)=zr(k,j)*tau(1,k)
       do 140 i=2,n
          l=i-1
          h=htsi(i,i)
          if (abs(h).lt.eps) goto 140            ! (h.eq.0.0)
          do 130 j=1,n
             s=0.
             si=0.
             do 110 k=1,l
                s=s+htsr(i,k)*zr(k,j)-htsi(i,k)*zi(k,j)
  110           si=si+htsr(i,k)*zi(k,j)+htsi(i,k)*zr(k,j)
             s=s/h/h
             si=si/h/h
             do 120 k=1,l
                zr(k,j)=zr(k,j)-s*htsr(i,k)-si*htsi(i,k)
  120           zi(k,j)=zi(k,j)-si*htsr(i,k)+s*htsi(i,k)
  130     continue
  140  continue

    end subroutine htribk

end module diagonalize

