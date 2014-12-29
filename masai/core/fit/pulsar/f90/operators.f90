!===============================================================================
! pulsar.operators
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


module operators
!---------------
                                                 ! todo: make this module able to handle more than two nuclei
                                                 !
  implicit none

  real*8, allocatable    :: unity(:,:)           ! identity matrix

  real*8, allocatable    :: s2(:,:)              ! s^2
  real*8, allocatable    :: sx(:,:)              ! sx
  real*8, allocatable    :: sy(:,:)              ! sy
  real*8, allocatable    :: sz(:,:)              ! sz
  real*8, allocatable    :: sq(:,:)              ! 3sz^2-s^2
  real*8, allocatable    :: sq1(:,:)             ! sz(4s^2-8sz^2-1)/2
  real*8, allocatable    :: sq2(:,:)             ! sz(2s^2-2sz^2-1)/2
  real*8, allocatable    :: sp(:,:)              ! s+
  real*8, allocatable    :: sm(:,:)              ! s-

  real*8, allocatable    :: i2(:,:)              ! i^2
  real*8, allocatable    :: ix(:,:)              ! ix
  real*8, allocatable    :: iy(:,:)              ! iy
  real*8, allocatable    :: iz(:,:)              ! iz
  real*8, allocatable    :: iq(:,:)              ! 3iz^2-i^2
  real*8, allocatable    :: iq1(:,:)             ! iz(4i^2-8iz^2-1)/2
  real*8, allocatable    :: iq2(:,:)             ! iz(2i^2-2iz^2-1)/2
  real*8, allocatable    :: ip(:,:)              ! i+
  real*8, allocatable    :: im(:,:)              ! i-

  real*8, allocatable    :: sziz(:,:)            ! product iz*sz

END module

subroutine reset_operators()
!--------------------------

use share
use operators

implicit none

                                             !loop counters
  integer*8 :: i
  integer*8 :: j

  ! first deallocate the operators array if they are still allocated
  if (allocated(sz)) deallocate (sz)
  if (allocated(s2)) deallocate (s2)
  if (allocated(sx)) deallocate (sx)
  if (allocated(sy)) deallocate (sy)
  if (allocated(sq)) deallocate (sq)
  if (allocated(sq1)) deallocate (sq1)
  if (allocated(sq2)) deallocate (sq2)

  if (allocated(iz)) deallocate (iz)
  if (allocated(i2)) deallocate (i2)
  if (allocated(ix)) deallocate (ix)
  if (allocated(iy)) deallocate (iy)
  if (allocated(iq)) deallocate (iq)
  if (allocated(iq1)) deallocate (iq1)
  if (allocated(iq2)) deallocate (iq2)

  if (allocated(sziz)) deallocate (sziz)
  if (allocated(unity)) deallocate (unity)

  if (allocated(sp)) deallocate (sp)
  if (allocated(sp)) deallocate (sm)
  if (allocated(im)) deallocate (im)
  if (allocated(ip)) deallocate (ip)
  ns=nint(2.0d0*spins+1.0d0)
  ni=nint(2.0d0*spini+1.0d0)
  nsi=ns*ni

  allocate(unity(nsi,nsi))
  unity=0.0d0
  do j=1,ni
   do i=1,ns
     unity(i+(j-1)*ns,i+(j-1)*ns)=1.0d0
   end do
  end do

  allocate(sz(nsi,nsi))
  allocate(s2(nsi,nsi))
  allocate(iz(nsi,nsi))
  allocate(i2(nsi,nsi))
  sz=0.0d0
  s2=0.0d0
  iz=0.0d0
  i2=0.0d0
  do j=1,ni
    do i=1,ns
      s2(i+(j-1)*ns,i+(j-1)*ns)=spins*(spins+1.0d0)
      sz(i+(j-1)*ns,i+(j-1)*ns)=spins-dble(i-1)
      i2(i+(j-1)*ns,i+(j-1)*ns)=spini*(spini+1.0d0)
      iz(i+(j-1)*ns,i+(j-1)*ns)=spini-dble(j-1)
    end do
  end do
  allocate(sp(nsi,nsi))
  allocate(sm(nsi,nsi))
  allocate(sx(nsi,nsi))
  allocate(sy(nsi,nsi))
  allocate(sq(nsi,nsi))
  allocate(sq1(nsi,nsi))
  allocate(sq2(nsi,nsi))
  sp=0.0D0
  sm=0.0D0
  sp=cshift(sqrt(s2-matmul(sz,sz-unity)),-1,2) ! circular shift of unity row or one column to create s+ or s-
  sm=transpose(sp)
  sx=(sp+sm)/2.0d0
  sy=-(sp-sm)/2.0d0
  sq=3.0d0*matmul(sz,sz)-s2
  sq1=matmul(sz,4.0d0*s2-8.0d0*matmul(sz,sz)-unity)/2.0d0
  sq2=matmul(sz,2.0d0*s2-2.0d0*matmul(sz,sz)-unity)/2.0d0
  deallocate(sp)
  deallocate(sm)

  allocate(ip(nsi,nsi))
  allocate(im(nsi,nsi))
  allocate(ix(nsi,nsi))
  allocate(iy(nsi,nsi))
  allocate(iq(nsi,nsi))
  allocate(iq1(nsi,nsi))
  allocate(iq2(nsi,nsi))
  ip=0.0D0
  im=0.0D0
  ip=cshift(sqrt(i2-matmul(iz,iz-unity)),-ns,2)! circular shift of ns row or ns column to create i+ or i-
  im=transpose(ip)
  ix=(ip+im)/2.0d0
  iy=-(ip-im)/2.0d0
  iq=3.0d0*matmul(iz,iz)-i2
  iq1=matmul(iz,4.0d0*i2-8.0d0*matmul(iz,iz)-unity)/2.0d0
  iq2=matmul(iz,2.0d0*i2-2.0d0*matmul(iz,iz)-unity)/2.0d0
  deallocate(ip)
  deallocate(im)

  allocate (sziz(nsi,nsi))
  sziz=0.0D0
  sziz=matmul(sz,iz)

end subroutine reset_operators

