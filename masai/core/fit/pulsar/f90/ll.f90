!===============================================================================
! pulsar.ll
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


logical*8 function isnull(angle)
!-------------------------------
  implicit none

  real*8 angle(0)

  real*8 eps
  parameter(eps=1.0d-30)

  isnull=.false.
  if (max(abs(maxval(angle)),abs(minval(angle))).le.eps) isnull=.true.

end function isnull


subroutine ll(xx,yy,zz,ii,jj,r3)
!-------------------------------

  use share
  use operators
  use diagonalize

  implicit none

                                                 ! args
                                                 !-----
  real*8, intent(in) :: xx
  real*8, intent(in) :: yy
  real*8, intent(in) :: zz
  integer*8, intent(in) :: ii
  integer*8, intent(in) :: jj
  real*8, intent(in) :: r3

                                                 ! external function
                                                 !------------------
!  interface isnull
!    logical*8 function isnull(angle)
!      real*8 :: angle(:)
!      logical*8 :: isnull
!    end function
!  end interface

  external isnull
  logical*8 :: isnull

                                                 ! loop counters
  integer*8 :: i,j,k,l, kk
  integer*8 :: insb,ig,nc
  integer*8 :: itp, ib, ictp

                                                 ! other parameters
  real*8 :: cb
  real*8 :: sb
  real*8 :: ca
  real*8 :: sa
  real*8 :: c2a
  real*8 :: s2a
  real*8 :: rrrr
  integer*8 :: ins

  real*8    :: ddm
  real*8    :: gam
  real*8    :: arg
  !integer*8 :: nphi
  real*8    :: tat
  real*8    :: tp
  real*8    :: dtpp


  real*8    :: v0qi
  real*8    :: v0qs
  real*8    :: v0ci
  real*8    :: v0dcs
  real*8    :: v0dip
  real*8    :: tetaf
  real*8    :: ctetaf
  real*8    :: stetaf
  real*8    :: v0qis
  real*8    :: v0qii
  real*8    :: v0dcis
  real*8    :: v0cii
  real*8    :: tetai
  real*8    :: dtetai
  real*8    :: dtetaf
  real*8    :: cteta
  real*8    :: steta
  real*8    :: vpij
  real*8    :: dvpij
  real*8    :: cvp
  real*8    :: svp
  integer*8 :: indi
  integer*8 :: indj
  !real*8    :: divise
  real*8    :: gt
  real*8    :: arg2
  complex*8    :: args
  complex*8 :: signal
  real*8    :: dem
  complex*8 :: recph


  cb=zz
  sb=sqrt(abs(1.0d0-zz*zz))
  if(abs(sb)>eps) then
    ca=xx/sb
    sa=yy/sb
    c2a=(xx*xx-yy*yy)/sb/sb
    s2a=2.0d0*xx*yy/sb/sb
  else
    ca=1.0d0
    sa=0.0d0
    c2a=1.0d0
    s2a=0.0d0
  end if
  

                                                 ! calculation of the am
                                                 !-----------------------
  ats=cmplx(0.0d0,0.0D0)
  aqi=cmplx(0.0d0,0.0D0)
  aci=cmplx(0.0d0,0.0D0)
  aqs=cmplx(0.0d0,0.0D0)
  acs=cmplx(0.0d0,0.0D0)
  ads=cmplx(0.0d0,0.0D0)
  p=cmplx(0.0d0,0.0D0)
  q=cmplx(0.0d0,0.0D0)

  if(abs(coefps).gt.eps) call am(etaqs,aqs)

  if(abs(csas).gt.eps) then
    if(isnull(csangle)) then
      call am(etacs,acs)
    else
      call amdif(r2mcs,acs)
    end if
  end if
  
  
  if (kcoupled>0) adks=0.0d0
  do ins=1,kcoupled
    if(abs(dipk(ins)).ge.0.1) then
      if (isnull(dskangle(ins,:))) then
         call am(0.0d0, adks(ins,:))             ! we fill first the adks array
      else
         call amdif(r2mdsk(ins,:),adks(ins,:))
      end if
    end if
  end do

  do l=0,2
    ats(l)=sq6*csas*acs(l)                       ! then we create the ats array
    do ins=1,kcoupled
      ats(l)=ats(l)+sq6*dipk(ins)*adks(ins,l)
    end do
  end do

  if(abs(dsi).gt.eps.and.icoupled.gt.0) then
    if(isnull(dsiangle)) then
      call am(0.0d0,ads)
    else
      call amdif(r2mdsi,ads)
    end if
  end if

  if(abs(coefpi).gt.eps.and.icoupled.gt.0) then
    if (isnull(qiangle)) then
      call am(etaqi,aqi)
    else
      call amdif(r2mqi,aqi)
    end if
  end if

  if(abs(csai).gt.eps.and.icoupled.gt.0) then
    if (isnull(ciangle)) then
      call am(etaci, aci)
    else
      call amdif(r2mci,aci)
    end if
  end if
  
  aqss=cmplx(0.0d0,0.0D0)
  acspi=cmplx(0.0d0,0.0D0)
  aqsi=cmplx(0.0d0,0.0D0)
  aqpi=cmplx(0.0d0,0.0D0)
  aqps=cmplx(0.0d0,0.0D0)
  adcps=cmplx(0.0d0,0.0D0)
  adipolis=cmplx(0.0d0,0.0D0)

  do l=0,4
    if (staticsample.and.(l.ne.0)) cycle
    rrrr=zl(l)*de(l,0)
    aqpi(l)=rrrr*aqi(l)
    acspi(l)=rrrr*aci(l)*sq6*csai
    adipolis(l)=rrrr*ads(l)*sq6*dsi
    aqps(l)=rrrr*aqs(l)
    adcps(l)=rrrr*ats(l)
    do k=1,2
      do i=l-2,2
        aqsi(k,l)=aqsi(k,l)+zl(l)*dk(k,l,i)*aqi(i)*aqi(l-i)
        aqss(k,l)=aqss(k,l)+zl(l)*dk(k,l,i)*aqs(i)*aqs(l-i)
      end do
    end do
  end do

                                                 ! wr(ii,jj,k) and p(k,nsbs) are independant of gamma and of the phase cycling
                                                 ! so they are calculated once only for each crystal orientations.

  a=cmplx(0.0d0,0.0D0)
  do k=1,ns-1
    ddm=3.0d0*(2.0d0*dble(k)-2.0d0*spins-1.0d0)*coefps
    wr(ii,jj,k)=isost+real(u1s(k)*aqss(1,0)+u2s(k)*aqss(2,0)+adcps(0)-ddm*aqps(0))
    if(staticsample) cycle
    do l=1,4
      a(k,l)= (u1s(k)*aqss(1,l)+u2s(k)*aqss(2,l)+adcps(l)-ddm*aqps(l))/wm(l)! complex a correspond to the previous ab and bb
    end do
  end do

  gam=-dg
  do i=1,ng
    gam=gam+dg
    do k=1,ns-1
      arg=0.0d0
      do l=1,4
        arg=arg+real(a(k,l))*dsin(l*gam)-aimag(a(k,l))*dcos(l*gam)
      end do
      do insb=-nsb,nsb
        p(k,insb)=p(k,insb)+cmplx(dcos(arg-dble(insb)*gam),dsin(arg-dble(insb)*gam))
      end do
    end do
  end do

  gam=-dg
  tat=0.0d0
  
  ! loop on gamma angles
  ! --------------------
  loop_gamma: do ig=1,ng
    
    gam=gam+dg
    

    ! loop on the coherence selection (min of ntcp is one)
    !-----------------------------------------------------
    loop_ctp: do ictp=1,nctp

      ! initial reduced density matrix
      !-------------------------------
      ro=sz+iz*vli/vls

      tat=0.0d0

      ! loop on the preparation period
      !-------------------------------
      loop_ncycles: do nc=1,ncycles

        !TENTATIVE: CASE OF IDEAL PULSE ON S
        if (idealpulse.AND.nc.EQ.1) then
           ro=-cmplx(0.,sy)
        end if

        if ((.not.idealpulse)) then          !TODO add pulseangle parameter
                                             !and make possible to have more than one idealpulse in the sequence

          ! period with rf-fields
          !----------------------
          tp=-deltatp(nc)*twopim6/2.0d0

                                                   !   lorsque le pulse est court (<1 tour), tri et tii ne servent � rien : matrices identit�es
                                                   !   lorsque le pulse est long  (>1 tour), tri et tii ont �t� r�-initialis�es (280) et servent
                                                   !   de point de d�part pour le dernier petit pas de pulse.

          ui=cmplx(unity)                            ! initialize to the identity matrix

                                                   ! pulse-decomposition in small rf-steps
                                                   !--------------------------------------
          loop_ntp: do itp=1,ntp(nc)

            dtpp=deltatp(nc)*twopim6
            tp=tp+dtpp

            if(itour(nc).ge.0.and.itp.eq.ntp(nc)) then
                                                   ! dans le cas o� le pulse est plus long qu'un tour, on traite � la  fin le dernier petit pulse
              dtpp=deltatpadd(nc)*twopim6
              if(dtpp.lt.0.0000001d0) cycle loop_ntp
              tp=rflength(nc)*twopim6-dtpp/2.0d0
            end if

            v0qi=real(aqpi(0))
            v0qs=real(aqps(0))
            v0ci=real(acspi(0))
            v0dcs=real(adcps(0))
            v0dip=real(adipolis(0))
            vvi(:)=real(aqsi(:,0))
            vvs(:)=real(aqss(:,0))

            if (.not.staticsample) then
                                                   ! voir le probleme de la multiplication par *twopim6  (unite)
              do l=1,4
                tetaf=wm(l)*(tat+tp)+l*gam
                ctetaf=dcos(tetaf)
                stetaf=dsin(tetaf)
                v0qi=v0qi+real(aqpi(l))*ctetaf+aimag(aqpi(l))*stetaf
                v0qs=v0qs+real(aqps(l))*ctetaf+aimag(aqps(l))*stetaf
                v0ci=v0ci+real(acspi(l))*ctetaf+aimag(acspi(l))*stetaf
                v0dcs=v0dcs+real(adcps(l))*ctetaf+aimag(adcps(l))*stetaf
                v0dip=v0dip+real(adipolis(l))*ctetaf+aimag(adipolis(l))*stetaf
                vvi(:)=vvi(:)+real(aqsi(:,l))*ctetaf+aimag(aqsi(:,l))*stetaf
                vvs(:)=vvs(:)+real(aqss(:,l))*ctetaf+aimag(aqss(:,l))*stetaf
              end do

            end if

            htst= cmplx(   (isoi-rfoffseti(nc)+v0ci)*iz                       &! real part
                         + (isost-rfoffsets(nc)+v0dcs)*sz                     &
                         + (v0dip+jsi)*sziz                                   &
                         + coefpi*v0qi*iq                                     &
                         + coefps*v0qs*sq                                     &
                         + coefsi*(vvi(1)*iq1+vvi(2)*iq2)                     &
                         + coefss*(vvs(1)*sq1+vvs(2)*sq2)                     &
                         + rfpoweri(nc)*ix*dcos(rfphasei(nc))                 &
                         + rfpowers(nc)*sx*dcos(rfphases(nc))                 &
                       ,   rfpoweri(nc)*iy*dsin(rfphasei(nc))                 &! imaginary
                         + rfpowers(nc)*sy*dsin(rfphases(nc)) )


            ! diagonalization
            !----------------
            call diag(nsi,htst,vp,z)

            ! calculation of u+.ro(t).u
            !--------------------------
            expvp=0.0d0
            do i=1,nsi
              expvp(i,i)=cmplx(dcos(vp(i)*dtpp),dsin(vp(i)*dtpp))! exp(i*vp)
            end do
            y=matmul(z,expvp)                      ! y=z*exp(i*vp)
            u=matmul(y,transpose(conjg(z)))        ! u=y.z+

            uf=matmul(ui,u)                        ! uf=ui.u

            if(itp.eq.ntpadd(nc)) uint=uf
            ! storage of the intermediate evolution operator in the case where the "pulselength > rotorperiod"

            ! pulselength > rotorperiod
            !--------------------------
            if(itour(nc).ge.0.and.itp.eq.(ntp(nc)-1)) then

              usi=uf
              usf=uf

                                                   !   � la fin du 1� tour,on calcule la matrice densit� correspondant � itour+1 tours
                                                   !   entiers: c'est la puissance itour+1 de celle correspondant � 1 tour.
                                                   !   on doit donc remultiplier itour fois celle correspondant � 1 tour.
              do ib=1,itour(nc)
                uf = matmul(usf,usi)
                usf = uf
              end do

                                                   !   apr�s itour+1 il y a un reste, dont la 1� partie correspond � ntpint pas de dur�e dtp.
                                                   !   l'op�rateur d'�volution correspondant � ces ntpint pas, a �t� stoqu� -> uint
              if(ntpadd(nc).gt.0) uf = matmul(usf,uint)

            end if

                                                   ! here we are, except for the last small pulse fraction
                                                   !------------------------------------------------------
            ui=uf                                  ! on red�marre de l� pour traiter ce dernier petit pulse: on repart � loop_ntp

          end do loop_ntp

          rot=matmul(ro,uf)
          ro=matmul(transpose(conjg(uf)),rot)

        end if ! idealpulse

        ! period without rf-field and no homonuclear dipolar couplings
        ! (todo : make a test to replace by a pulse in this case)
        !-------------------------------------------------------------

        if (rfdelay(nc).gt.eps) then

          v0qis= real(aqps(0))*rfdelay(nc)*twopim6
          v0qii= real(aqpi(0))*rfdelay(nc)*twopim6
          v0dcis= real(adcps(0))*rfdelay(nc)*twopim6
          v0cii= real(acspi(0))*rfdelay(nc)*twopim6
          v0dip= real(adipolis(0))*rfdelay(nc)*twopim6
          vvis(:)= real(aqss(:,0))*rfdelay(nc)*twopim6
          vvii(:)= real(aqsi(:,0))*rfdelay(nc)*twopim6

          if(.not.staticsample) then

            do l=1,4
              tetai=wm(l)*(tat+rflength(nc)*twopim6)+l*gam
              tetaf=tetai+wm(l)*rfdelay(nc)*twopim6
              dtetai=mod(tetai,twopi)
              dtetaf=mod(tetaf,twopi)
              cteta=(dcos(dtetaf)-dcos(dtetai))/wm(l)
              steta=(dsin(dtetaf)-dsin(dtetai))/wm(l)
              v0qis=v0qis+real(aqps(l))*steta-aimag(aqps(l))*cteta
              v0qii=v0qii+real(aqpi(l))*steta-aimag(aqpi(l))*cteta
              v0dcis=v0dcis+real(adcps(l))*steta-aimag(adcps(l))*cteta
              v0cii=v0cii+real(acspi(l))*steta-aimag(acspi(l))*cteta
              v0dip=v0dip+real(adipolis(l))*steta-aimag(adipolis(l))*cteta
              vvis(:)=vvis(:)+real(aqss(:,l))*steta-aimag(aqss(:,l))*cteta
              vvii(:)=vvii(:)+real(aqsi(:,l))*steta-aimag(aqsi(:,l))*cteta
            end do

          end if                                 ! not static sample

          htsr= 0.0d0
          htsr=    ((isost-rfoffsets(nc))*rfdelay(nc)*twopim6+v0dcis)*sz  &
                 + (( isoi-rfoffseti(nc))*rfdelay(nc)*twopim6+v0cii)*iz   &
                 + decouple(nc)*(v0dip+jsi*rfdelay(nc)*twopim6)*sziz      &
                 + coefps*v0qis*sq                                        &
                 + coefpi*v0qii*iq                                        &
                 + coefss*(vvis(1)*sq1+vvis(2)*sq2)                       &
                 + coefsi*(vvii(1)*iq1+vvii(2)*iq2)


          do i=1,nsi
            do j=1,nsi
              vpij=htsr(i,i)-htsr(j,j)
              dvpij=mod(vpij,twopi)
              cvp=dcos(dvpij)
              svp=dsin(dvpij)
              ro(i,j)=ro(i,j)*cmplx(cvp,-svp)
            end do
          end do

        end if                                   ! rfdelays

        ! effect of t2 on each coherences
        !--------------------------------
        if (t2si(nc).le.0.99999999.or.t2ss(nc).le.0.99999999.or.t2ii(nc).le.0.99999999) then
          do i=1,nsi
            do j=1,nsi
              indi=int((i-1)/ns)
              indj=int((j-1)/ns)
              if(indi.eq.indj.and.i.ne.j) ro(i,j)=ro(i,j)*t2ss(nc)
              if(indi.ne.indj) then
                if((i-j).eq.(ns*((i-j)/ns))) then
                  ro(i,j)=ro(i,j)*t2ii(nc)
                else
                  ro(i,j)=ro(i,j)*t2si(nc)
                end if
              end if
            end do
          end do
        end if

        ! Coherence filtering  (up to now, only on the spin S)
        ! TODO : allow also a selection on spin I
        !-----------------------------------------------------
       
        if (allocated(ctp))then
          if (nc.ne.ncycles.AND.abs(nint(ctp(ictp,nc+1))).lt. 999) then
                                                                   ! pass if ctp doesn't exist or if it is the last pulse
                                                                   ! in the latter case we don't care about selecting something
                                                                   ! (the receiver select -1 anyway)

          ! we first store the roI density matrix for spin I, because it must remain unaffected by the selection on spin S
          !--------------------------------------------------------------------------------------------------------------
          sti=0.0d0
          if (ni>1) then
            do i=1,ni
              do j=1,ni
                do kk=1,ns
                  sti(i,j)=sti(i,j)+ro(kk+(i-1)*ns,kk+(j-1)*ns)/dble(ns)
                                                   !
                                                   ! TODO : check this (because the sum of the diagonal element is not zero
                                                   ! when there is no I spin?)
                end do
              end do
            end do
          end if

                                                   ! when we make a selection directly in the density matrix (for instance on the spin S),
                                                   ! we must keep the coherences ((mS+p,mI),(mS,mI)), so with mI kept constant. this is checked using
                                                   ! two parameters indi and indj

          do i=1,nsi
            do j=1,nsi
              indi=int((i-1)/ns)                    ! these two lines calculate the corresponding index on spin I (mI)
              indj=int((j-1)/ns)

              ! the selection must be done at mI constant, and thus when indi=indj
              ! the corresponding element is kept
              !--------------------------------------------------------------------

              ! we use nc+1 because the list ctp have its first element to 0 (before the first pulse)

                if((indi.eq.indj).AND.((j-i).eq.Nint(ctp(ictp,nc+1)))) then
                  cycle
                END if

             ! if ( abs(nint(ctp(ictp,nc+1))-10000) .lt. 999) then
             !   ! case of the pathway compression (symetrical pathways)
             !
             !   if(indi.eq.indj.and.abs(j-i).eq.abs(nint(ctp(ictp,nc+1))-10000)) cycle
             !
             ! end if

              ! if it is a non selected coherence, then it is set to 0
              !-------------------------------------------------------
              ro(i,j)=0.0d0

              ! now restore roI constant for all values of mS corresponding to identical values of mI.
              !---------------------------------------------------------------------------------------
              if((i-j).eq.(ns*((i-j)/ns))) ro(i,j)=sti(indi+1,indj+1)
            end do
          end do

         end if 
        end if ! test allocation of ctp
 
        tat=tat+(rflength(nc)+rfdelay(nc))*twopim6

      ! goto the calculation of the next pulse+delay
      !---------------------------------------------

      end do loop_ncycles

      !store the density matrix for the given CTP
      !------------------------------------------
      roctp(ictp,:,:)=ro(:,:)

    ! go to the calculation of the next CTP
    !-------------------------------------
    end do loop_ctp


    ! we compute the sum of the different CTP (sum along the dimension 1 of roctp)
    !-----------------------------------------------------------------------------
    ro=SUM(roctp,DIM=1)/DBLE(nctp)

    ! we compute the total density matrix weighted (R3) by the polar angle.
    !----------------------------------------------------------------------

    rott=rott+ro/r3

    ! End of the preparation period
    !------------------------------
    gt=spinningspeed*tat+gam
    do k=1,ns-1
      arg2=0.0d0
      do l=1,4
        arg2=arg2+aimag(a(k,l))*dcos(l*gt)-real(a(k,l))*dsin(l*gt)
      end do
      do insb=-nsb,nsb
        args=cmplx(dcos(arg2+insb*gt),dsin(arg2+insb*gt))
        signal=0.0d0
        do kk=0,ni-1
          signal=signal+ro(k+1+kk*ns,k+kk*ns)
        end do
        q(k,insb)=q(k,insb)+signal*args
      end do
    end do

  end do loop_gamma

  if (spectrum_compute) then

    ! observation period : sx et -sy
    !-------------------------------
    recph=exp(cmplx(0.d0,(rcph-90.)*pi/180.0d0))     ! apply receiver phase -90 when rcph=0.
    do k=1,ns-1
      do insb=-nsb,nsb
        dem=dnorm*sx(k,k+1)/(1.0d0+qq*(wr(ii,jj,k)+dble(insb)*spinningspeed)**2)/r3
        qt(ii,jj,k,insb)=qt(ii,jj,k,insb)+q(k,insb)
        amp(ii,jj,k,insb)=conjg(p(k,insb)*qt(ii,jj,k,insb))*dem
        amp(ii,jj,k,insb)=amp(ii,jj,k,insb)*recph
      end do
    end do

  END if

  contains
  !------------------------!
  ! ll internal procedures !
  !------------------------!
    subroutine am(eta,aa)
    !--------------------
      implicit none

      complex*8 :: aa(-4:4)
      real*8 :: eta

      aa(-2)= cmplx(sb*sb/2.0d0 -eta*(1.0d0+cb*cb)*c2a/6.0d0, -eta*cb*s2a/3.0d0)
      aa(-1)= cmplx(-sb*cb*(eta*c2a/3.0d0+1.0d0), -eta*sb*s2a/3.0d0)
      aa(0) = cmplx((3.0d0*cb*cb-1.0d0 -eta*sb*sb*c2a)/sq6)
      aa(1) =-conjg(aa(-1))
      aa(2) = conjg(aa(-2))

    end subroutine am

    subroutine amdif(r2m,aa)
    !-----------------------
      implicit none

      complex*8 :: aa(-4:4)
      complex*8 :: r2m(0:2)
      real*8 :: sq15
      parameter(sq15=1.2247448713915890490986420373529d0)
      real*8 :: fgm1,fgm2,fgp1,fgp2

      fgm1=ca*real(r2m(1))-sa*aimag(r2m(1))
      fgm2=c2a*real(r2m(2))-s2a*aimag(r2m(2))
      fgp1=ca*aimag(r2m(1))+sa*real(r2m(1))
      fgp2=c2a*aimag(r2m(2))+s2a*real(r2m(2))

      aa(0)=  cmplx((3.0d0*zz*zz-1.0d0)*real(r2m(0))/2.0d0+sq15*sb*(sb*fgm2-2.0d0*zz*fgm1))
      aa(1)=  cmplx(sb*zz*(sq15*real(r2m(0))-fgm2)+(2.0d0*zz*zz-1.0d0)*fgm1, sb*fgp2-zz*fgp1)
      aa(2)=  cmplx((sq15*sb*sb*real(r2m(0))+(1.0d0+zz*zz)*fgm2)/2.0d0+sb*zz*fgm1, -zz*fgp2-sb*fgp1)
      aa(-2)= conjg(aa(2))
      aa(-1)=-conjg(aa(1))

    end subroutine amdif

end subroutine ll
