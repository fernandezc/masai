!===============================================================================
! pulsar.util
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


real*8 function radians(degrees)
!-------------------------------!
! conversion degrees to radians !
!-------------------------------!

  real*8,intent(in) :: degrees
  real*8 :: pi2
  parameter (pi2=0.017453292519943295769236907684886d0)
  radians=degrees*pi2
  return
end function radians

!logical*8 function isnull(angle)
!!-------------------------------
!  implicit none!
!
!  real*8 angle(:)!
!
!  real*8 eps
!  parameter(eps=1.0d-30)!
!
!  isnull=.false.
!  if (max(abs(maxval(angle)),abs(minval(angle))).le.eps) isnull=.true.
!
!end function isnull


subroutine wigner(angle,dw)
!-------------------------------------------!
! computation of the wigner matrix elements !
!-------------------------------------------!

  real*8, intent(in)  :: angle
  real*8, intent(inout) :: dw(-4:4,-4:4)

                                                 ! local variable
  real*8    :: c
  real*8    :: s

  real*8    :: sq38
  parameter ( sq38=0.61237243569579452d0 )

                                                 ! functions type declarations
  real*8    :: radians

  c=dcos(radians(angle))
  s=dsin(radians(angle))

  dw=0.0d0

  dw(-2,-2)=((1.0d0+c)**2)/4.0d0
  dw(2,2)=dw(-2,-2)
  dw(2,1)=-(1.0d0+c)*s/2.0d0
  dw(-2,-1)=-dw(2,1)
  dw(1,2)=-dw(2,1)
  dw(-1,-2)=dw(2,1)
  dw(2,0)=sq38*s*s
  dw(0,2)=dw(2,0)
  dw(-2,0)=dw(2,0)
  dw(0,-2)=dw(2,0)
  dw(2,-1)=-(1.0d0-c)*s/2.0d0
  dw(-1,2)=-dw(2,-1)
  dw(-2,1)=-dw(2,-1)
  dw(1,-2)=dw(2,-1)
  dw(2,-2)=((1.0d0-c)**2)/4.0d0
  dw(-2,2)=dw(2,-2)
  dw(1,1)=c*c-(1.0d0-c)/2.0d0
  dw(-1,-1)=dw(1,1)
  dw(1,-1)=(1.0d0+c)/2.0d0-c*c
  dw(-1,1)=dw(1,-1)
  dw(1,0)=-sq38*2.0d0*c*s
  dw(0,1)=-dw(1,0)
  dw(-1,0)=-dw(1,0)
  dw(0,-1)=dw(1,0)
  dw(0,0)=(3.0d0*c*c-1.0d0)/2.0d0

end subroutine wigner


subroutine orientate(euler,eta,r2m)
!------------------------------------------------------------!
! calculation of the r2m tensor components in a the m frame  !
!                                                            !
!      pas--(alpha,beta,gamma)-->mpas                        !
!                                                            !
!  mpas       2                             pas              !
! r       =  sum    d   (alpha,beta,gamma) r                 !
!  2m       m'=-2    m'm                    2m'              !
!------------------------------------------------------------!
  implicit none

                                                 ! args
  real*8, intent(in)  :: euler(3)
  real*8, intent(in)  :: eta
  complex*8, intent(out) :: r2m(0:2)
                                                 ! local
  real*8 :: dw(-4:4,-4:4)
  real*8 :: alpha, beta, gamma
  real*8 :: v2
  integer*2 :: m
  real*8 :: dm

  real*8 :: v0
  parameter (v0=0.81649658092772603273242802490196d0)! sqrt(2/3)

  real*8 :: radians

  alpha=radians(euler(1))
  beta=euler(2)                                  ! not yet converted to radians as the wigner function will do it
  gamma=radians(euler(3))

  v2=-eta/3.0d0                                  ! todo: why -eta/3 and not +eta

  call wigner(beta,dw)

  do m=1,2
    dm=dble(m)
    r2m(m)=cmplx(v0*dw(0,m)*dcos(dm*gamma)+ &
           v2*(dw(2,m)*dcos(dm*gamma+2.0d0*alpha)+dw(-2,m)*dcos(dm*gamma-2.0d0*alpha)), &
                 v0*dw(0,m)*dsin(dm*gamma)+ &
           v2*(dw(2,m)*dsin(dm*gamma+2.0d0*alpha)+dw(-2,m)*dsin(dm*gamma-2.0d0*alpha)))

  end do

  r2m(0)=cmplx(v0*dw(0,0)+2.0d0*v2*dw(2,0)*dcos(2.0d0*alpha))

end subroutine orientate

function num2str_i4n(i4n,fmat)
!-------------------------------------------
  implicit none
  integer*4, intent(in) :: i4n
  character(len=*), intent(in) ::  fmat
  character(len=32) ::  num2str_i4N
  integer*4 :: iochk
  num2str_i4n=''
  write(num2str_i4n,fmt=fmat,iostat=iochk) i4n
  if (iochk.ne.0) then
    num2str_i4n='# INVALID #'
  endif
  return
END function num2str_i4n

function num2str_i8n(i8n,fmat)
!-------------------------------------------
  implicit none
  integer*8, intent(in) :: i8n
  character(len=*), intent(in) ::  fmat
  character(len=32) ::  num2str_i8N
  integer*4 :: iochk
  num2str_i8n=''
  write(num2str_i8n,fmt=fmat,iostat=iochk) i8n
  if (iochk.ne.0) then
    num2str_i8n='# INVALID #'
  endif
  return
END function num2str_i8n

function num2str_r4n(r4n,fmat)
!-------------------------------------------
  implicit none
  real*4, intent(in) :: r4n
  character(len=*), intent(in) ::  fmat
  character(len=32) ::  num2str_r4N
  integer*4 :: iochk
  num2str_r4n=''
  write(num2str_r4n,fmt=fmat,iostat=iochk) r4n
  if (iochk.ne.0) then
    num2str_r4n='# INVALID #'
  endif
  return
END function num2str_r4n

function num2str_r8n(r8n,fmat)
!-------------------------------------------
  implicit none
  real*8 , intent(in) :: r8n
  character(len=*), intent(in) ::  fmat
  character(len=32) ::  num2str_r8N
  integer*4 :: iochk
!     erase any information in str.
  num2str_r8n=''
!     attempt to write i4n to str.
  write(num2str_r8n,fmt=fmat,iostat=iochk) r8n
  if (iochk.ne.0) then
!       invalid operation
    num2str_r8n='# INVALID #'
  endif
  return
END function num2str_r8n
!
