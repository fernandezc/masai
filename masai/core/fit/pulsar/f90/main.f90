!===============================================================================
! pulsar.main
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

! TODO: Rewrite using less allocated arrays as this appear very time consuming

! reset ========================================================================
subroutine reset()
    ! Reset and deallocation of all arrays

    use diagonalize

    call dealloc()
    call reset_share()
    call reset_operators()
    call reset_parameters()

end subroutine reset
!===============================================================================

! compute ======================================================================
logical*1 function compute()
    ! main fonction called by python to compute a NMR spectra

    use operators
    use share

    implicit none

    ! function interfaces
    interface radians
        real*8 function radians(degrees)
            real*8,intent(in) :: degrees
        end function radians
    end interface

    interface wigner
        subroutine wigner(angle,dw)
        real*8, intent(in)  :: angle
        real*8, intent(inout) :: dw(-4:4,-4:4)
      end subroutine wigner
    end interface

    interface num2str
        function num2str_i4n(i4n,fmat)
        implicit none
        integer*4, intent(in) :: i4n
        character(len=*), intent(in) :: fmat
        character(len=32) num2str_i4n
      end function num2str_i4n
        function num2str_i8n(i8n,fmat)
        implicit none
        integer*8, intent(in) :: i8n
        character(len=*), intent(in) :: fmat
        character(len=32) num2str_i8n
      end function num2str_i8n
        function num2str_r4n(r4n,fmat)
        implicit none
        real*4, intent(in) :: r4n
        character(len=*), intent(in) :: fmat
        character(len=32) num2str_r4n
      end function num2str_r4n
        function num2str_r8n(r8n,fmat)
        implicit none
        real*8, intent(in) :: r8n
        character(len=*), intent(in) :: fmat
        character(len=32) num2str_r8n
      end function num2str_r8n
    end interface num2str

    ! loop index and flags
    integer*4 :: i=0
    integer*4 :: j=0
    integer*4 :: k=0
    integer*4 :: l=0
    integer*4 :: nc=0

    ! local parameters
    real*8    :: powder1=0.0d0
    real*8    :: powder2=0.0d0
    integer*8 :: inorm1(16)=(/2154752, 933774, 589820, 430020, &
                              338069, 278417, 236611, 205703, &
                              181927, 163068, 147749, 135059, &
                              124378, 115254, 107382, 100517/)
    integer*8 :: inorm2(16)=(/6122925, 3121593, 2088467, 1568265, &
                              1255355, 1046463, 897112, 785115, &
                              697918, 628170, 571061, 523524, &
                              483263, 448748, 418821, 392639/)
    real*8    :: rotorperiod=0.0d0 ! rotor period in usec (1e6/spinningspeed)
    real*8    :: rpf=0.0d0 ! rotor period fraction (due to rfstep)
    real*8    :: dms=0.0d0
    real*8    :: coefdip=0.0d0
    complex*8 :: tempr2m(0:2)=0.0d0
    real*8    :: tempangle(3)=0.0d0
    logical*8 :: wholesphere=.false.
    integer*8 :: nscans=1
    integer*8 :: ncoupled=0

    !----------------------------------------------!
    ! check parameters and prepare the calculation !
    !----------------------------------------------!

    ! check allocation of the parameters arrays
    ! in principle there is no need to do it here, as it is done by the python interface
    ! but I did it in the case where the program is used without python


    if (allocated(nucleus)) then
        if (size(nucleus,1).EQ.1) then
            if (allocated(indirect)) deallocate(indirect)
            if (allocated(dipole)) deallocate(dipole)
        end if
    else
        compute=.false.
        return
    end if

    !-------------------------------------------
    ! check if a sequence is defined...
    ! if it is not, there is no need to continue
    !-------------------------------------------
    if (allocated(pulse)) then
        ncycles=size(pulse,1)
    else
        compute=.false.
        return
    end if

    ! ------------------------------------------------
    ! Start evaluation of parameters and controls them
    ! ------------------------------------------------

    ! spectrometer
    ! TODO : Change the definition of NT and NG
    nt=5*accuracy
    powder1=dble(inorm1(nt/5))/1.0d7
    powder2=dble(inorm2(nt/5))/1.0d7

    ! spinning
    staticsample=.false.
    if(spinningangle.le.eps.or.spinningspeed.le.eps) then
        staticsample=.true.
        ng=1
        spinningangle=0.0d0
        spinningspeed=0.0d0
        nsb=0
    end if
    if (spinningspeed.gt.eps) then
        rotorperiod=1.0d6/spinningspeed ! in us
        ng=nint(sw/spinningspeed)+1
    else
        rotorperiod=1.0d30
    end if

    ng=max(nt,ng)
    dg=2.0d0*pi/dble(ng)
    ! cos and sin are not reconized by the linker when using f2py (on windows at least)
    ! (i don't know why as it works with g95: may be a problem with librairies)
    ! cos and sin are thus replaced by dcos and dsin (which means that all computations
    ! will be done in double-precision real*8)
    ! i decided to pass all arguments in double precision real*8 and integer*8 in the
    ! remaining of the program for sake of simplicity
    ! (may be this can be changed if there is a notable effect on speedness)
    de=0.0D0
    call wigner(spinningangle,de)

    do k=1,2
        do l=0,4
            do i=l-2,2
                dk(k,l,i)=de(i,k)*de(l-i,-k)
            end do
        end do
    end do

    ! spin system
    if (allocated(nucleus)) then
        spins=nucleus(1,2) ! spin of s
        vls=nucleus(1,3)   ! larmor of s (in hz)

        !-----------------------------------------------------------------------------
        !TODO: I have to CHECK if the sign of vls (gyromagnetic) may be usefull later!
        ! I didn't use it here due to the problems that arise in definition
        ! of chemical shift and so on
        !-----------------------------------------------------------------------------
        abundance=nucleus(1,4)
        isos=chemicalshift(1,2)+sr  ! isotropic chemical shift of s (in hz)
                                    ! modified by SR frequency
        etacs=chemicalshift(1,4)! csa assymetry of s
        csas=chemicalshift(1,3)/(3.0d0+etacs)! csa of s
        cqs=quadrupole(1,2)     ! cq of s
        etaqs=quadrupole(1,3)! etaq of s
        t2s=t2(1,2)! relaxation t2 of s
        do i=1,3
            csangle(i)=chemicalshift(1,4+i)
            qsangle(i)=quadrupole(1,3+i)
        end do
    else
        compute=.false.
        return
    endif

    ncoupled= size(nucleus,1)-1
    if (ncoupled.ge.1) then
        icoupled=1
        spini=nucleus(2,2)! spin of i
        vli=nucleus(2,3)! larmor of i
        if (spini.lt.eps.or.vli.lt.eps) then
            icoupled=0
            spini=0.0D0
            vli=0.0D0
        else
            isoi=chemicalshift(2,2)! isotropic chemical shift of i
            etaci=chemicalshift(2,4)! csa assymetry of i
            csai=chemicalshift(2,3)/(3.0d0+etaci)! csa of i
            cqi=quadrupole(2,2)! cq of i
            etaqi=quadrupole(2,3)! etaq of i
            t2i=t2(2,2)! relaxation t2 of i
            do i=1,3
                ciangle(i)=chemicalshift(2,4+i)
                qiangle(i)=quadrupole(2,3+i)
            end do
            ! the definition of j and d coupling is reported later
            jsi=0.0d0
            dsi=0.0d0
            dsiangle=0.0d0
        end if
    end if

    if (ncoupled.ge.2) then
        kcoupled=ncoupled-1
        allocate (spink(kcoupled))! spin of k
        allocate (vlk(kcoupled))! larmor of k
        allocate (jsk(kcoupled))! j s-k
        allocate (dsk(kcoupled))! dipolar s-k
        allocate (dskangle(kcoupled,3))! orientation of s-k dipole
        do i=1,kcoupled
            spink(i)=nucleus(i+2,2)
            vlk(i)=nucleus(i+2,3)
        end do
        jsk=0.0d0! these arrays will be set later
        dsk=0.0d0
        dskangle=0.0d0
    end if

    ! j coupling
    if (allocated(indirect)) then
        do i=1,size(indirect,1)
            ! remember or note that the first nucleus is s
            ! the second nucleus is i and the other are the k nuclei
            if ((nint(indirect(i,2)).eq.2).and.(icoupled.gt.0)) then
                ! then it is a j coupling between s and i
                jsi=indirect(i,3)
                ! the anisotropy of j is not yet implemented
                ! (even if the table can handle the definition of deltaj)
            end if
            if ((nint(indirect(i,2)).gt.2).and.(kcoupled.gt.0))  then
                ! in this case it is a j coupling between s and k
                jsk(nint(indirect(i,2))-2)=indirect(i,3)
            end if
        end do
    end if

    ! dipolar coupling
    if (allocated(dipole)) then
        do i=1,size(dipole,1)
            ! remember or note that the first nucleus is s
            ! the second nucleus is i and the other are the k nuclei
            if ((nint(dipole(i,2)).eq.2).and.(icoupled.gt.0)) then
                ! then it is a dipole coupling between s and i
                dsi=dipole(i,3)
                do j=1,3
                    dsiangle(j)=dipole(i,4+j)
                end do
            end if
            if ((dipole(i,2).gt.2).and.(kcoupled.gt.0))  then
                ! in this case it is a dipolar coupling between s and k
                dsk(nint(dipole(i,2))-2)=dipole(i,3)
                do j=1,3
                    dskangle(nint(dipole(i,2))-2,j)=dipole(i,4+j)
                end do
            end if
        end do
    end if

    ! spectrum limit
    fst=-sw/2.0d0! i think this can be changed later

    ! simulation parameters
    if (rfstep.lt.0.10d0) rfstep=0.1d0
    if (rfstep.gt.360.0d0) rfstep=360.0d0
    rpf=rotorperiod*rfstep/360.0d0! rpf: rotorperiod fraction in us
    qq=4.0d0*qfactor*qfactor/vls/vls

    do l=1,4
        wm(l)=dble(l)*spinningspeed
    end do

    ! detection: all=0, central=1, satellite=2
    if (nall.eq.999) nall=0
    if (spins<0.4) nall=1
    if (mod(spins,1.0d0).lt.0.001) nall=0 ! in the case of an integer spin s,
                                          !all transitions must be calculated

    ! define all operators
    call reset_operators()

    ! preparation pulses and delays
    allocate (rfpowers(ncycles))
    allocate (rfoffsets(ncycles))
    allocate (rfphases(ncycles))

    allocate (rfpoweri(ncycles))
    allocate (rfoffseti(ncycles))
    allocate (rfphasei(ncycles))

    allocate (rflength(ncycles))

    allocate (rfdelay(ncycles))
    allocate (decouple(ncycles))

    allocate (itour(ncycles))
    allocate (ntp(ncycles))
    allocate (deltatp(ncycles))
    allocate (ntpadd(ncycles))
    allocate (deltatpadd(ncycles))
    allocate (t2ss(ncycles))
    allocate (t2ii(ncycles))
    allocate (t2si(ncycles))

    do nc=1,ncycles
        rflength(nc)=pulse(nc,2)
        if (abs(rflength(nc)).le.eps) then
            idealpulse=.true.
        end if
        rfpowers(nc)=pulse(nc,3)
        rfoffsets(nc)=pulse(nc,4)
        rfphases(nc)=radians(pulse(nc,5))
        rfpoweri(nc)=pulse(nc,6)
        rfoffseti(nc)=pulse(nc,7)
        rfphasei(nc)=radians(pulse(nc,8))
        
        ! calculate the increment dtp for pulse effect integration
        itour(nc)=floor(rflength(nc)/rotorperiod)-1

        ! when the pulse length is less than one rotor period (itour=0),
        ! the pulse is divided in ntp intervals with duration dtp
        ! else each of the itour rotor period are divided in the same way
        ! in ntp intervals (duration deltatp)
        ! and the last period is divided in ntpadd intervals of the same duration,
        ! plus a small additional period dtpadd<deltatp

        ntp(nc)=min(rflength(nc),rotorperiod)/rpf+1
        deltatp(nc)=min(rflength(nc),rotorperiod)/dble(ntp(nc))
        !  i don't know why dtp doesnt work ... changed to deltatp
        !  probably a reserved keyword in f2py
        deltatpadd(nc)=0.0d0
        ntpadd(nc)=0
        if (itour(nc).ge.0) then
            ntpadd(nc)=floor( (rflength(nc)-dble(itour(nc)+1)*rotorperiod)/deltatp(nc) )
            deltatpadd(nc)=rflength(nc) &
                             -dble(itour(nc)+1)*rotorperiod-dble(ntpadd(nc))*deltatp(nc)
        end if

        rfdelay(nc)=delay(nc,2)
        t2ss(nc)=exp(-(rflength(nc)+rfdelay(nc))/t2s)
        t2ii(nc)=exp(-(rflength(nc)+rfdelay(nc))/t2i)
        t2si(nc)=(t2ss(nc)+t2ii(nc))/2.0d0
        
        ! decoupling !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! To check !!!!!!!
        if(delay(nc,3)>.1d0) then
            ! warning when decouple = 1 this means no decoupling...
            decouple(nc)=1
        else
            decouple(nc)=0
        end if
    end do

    ! coefficients for the quadrupole
    coefps=0.0d0
    coefpi=0.0d0
    if(spins.gt.0.6) coefps=-0.5d0*sq6*cqs/(4.0d0*spins*(2.0d0*spins-1.0d0))
    if(spini.gt.0.6) coefpi=-0.5d0*sq6*cqi/(4.0d0*spini*(2.0d0*spini-1.0d0))
    coefss=6.0d0*coefps**2/vls
    coefsi=0.0D0
    if (icoupled>0) coefsi=6.0d0*coefpi**2/vli

    allocate(u1s(ns),u2s(ns))

    do k=1,ns-1
        dms=spins+1.0d0 - dble(k)
        u1s(k)=-0.50d0*coefss*(24.0d0*dms*(dms-1.0d0)-4.0d0*spins*(spins+1.0d0)+9.0d0)
        u2s(k)=-0.25d0*coefss*(12.0d0*dms*(dms-1.0d0)-4.0d0*spins*(spins+1.0d0)+6.0d0)
    end do

    ! tensor component of the dipolar tensor / quadrupole s
    ! (todo: change this to have a reference with respect to the molecular frame)

    r2mdsi=cmplx(0.0d0,0.0D0)
    if (icoupled>0) call orientate(-dsiangle,0.0d0,r2mdsi)
    ! the sign - before dsiangle depend on the definition choosen (euler or polar angle)

    ! quadrupole i / quadrupole s
    r2mqi=cmplx(0.0d0,0.0D0)
    if (icoupled>0) call orientate(qiangle,etaqi,r2mqi)

    ! csa s / quadrupole s
    r2mcs=cmplx(0.0d0,0.0D0)
    call orientate(csangle,etacs,r2mcs)

    ! csa i / quadrupole s
    r2mci=cmplx(0.0d0,0.0D0)
    if (icoupled>0) call orientate(ciangle,etaci,r2mci)

    ! other dipolar coupling with non-excited nuclei (k)
    coefdip=1.0d0
    if (kcoupled>0) then
        allocate(r2mdsk(kcoupled,0:2))
        allocate(nk(kcoupled))
        allocate (dipk(kcoupled))
        allocate (mk(kcoupled))
        dipk=0.0d0
        mk=0
        r2mdsk=cmplx(0.0d0,0.0D0)
        do k=1,kcoupled
            nk(k)=nint(2.0d0*spink(k))
            coefdip=coefdip/(2.0d0*spink(k)+1.0d0)
            tempr2m=cmplx(0.0d0,0.0d0)
            do j=0,2
                tempangle(j)=dskangle(k,j)
            end do
            call orientate(-tempangle,0.0d0,tempr2m)
            do j=0,2
                r2mdsk(k,j)=tempr2m(j)
            end do
        end do
    end if

    ! frequencies and amplitude
    allocate (wr(0:nt,0:2*nt,ns))
    wr=0.0d0

    allocate (amp(0:nt,0:2*nt,ns,-nsb:nsb))
    amp=cmplx(0.0d0,0.0D0)

    ! ll arrays
    ! they are allocated here to avoid the multiple allocation-deallocation in ll
    if (kcoupled>0) allocate (adks(kcoupled,-4:4))
    allocate (a(ns,4))
    a=0.0d0

    allocate (p(ns,-nsb:nsb))
    p=0.0d0

    allocate (u(nsi,nsi))
    u=0.0d0

    allocate (ui(nsi,nsi))
    ui=0.0d0
    allocate (uf(nsi,nsi))
    uf=0.0d0
    allocate (usi(nsi,nsi))
    usi=0.0d0
    allocate (usf(nsi,nsi))
    usf=0.0d0
    allocate (uint(nsi,nsi))
    uint=0.0d0
    allocate (ro(nsi,nsi))
    ro=0.0d0
    allocate (rot(nsi,nsi))
    rot=0.0d0
    allocate (z(nsi,nsi))
    z=0.0d0
    allocate (vp(nsi))
    vp=0.0d0
    allocate (expvp(nsi,nsi))
    expvp=0.0d0
    allocate (y(nsi,nsi))
    y=0.0d0
    allocate (htst(nsi,nsi))
    htst=0.0d0
    allocate (htsr(nsi,nsi))
    htsr=0.0d0
    allocate (sti(ni,ni))
    sti=0.0d0
    allocate (rosto(nsi,nsi))
    rosto=0.0d0
    allocate (q(ns,-nsb:nsb))
    q=0.0d0
    ! and the density matrix
    allocate (rott(nsi,nsi))! rortt and roitt
    rott=0.0d0

    ! allocate some other important arrays
    allocate (qt(0:nt,0:2*nt,ns,-nsb:nsb))
    qt=cmplx(0.0d0,0.0d0)

    ! spec will contain the simulated spectrum
    ! (if zerospectrum is false then it is not erased for addition of data)
    ! this is not used in masai TOSUPPRESS
    zerospectrum = .true.
    if (zerospectrum.and.allocated(spec)) then
        deallocate (spec)
        allocate (spec(npts))
        spec=cmplx(0.0d0,0.0d0)
        sumc = cmplx(0.0d0, 0.0d0)
        sums = cmplx(0.0d0, 0.0d0)
    end if
    if (.not.allocated(spec)) then
    	! in all case if it is not allocated, we have to allocate it
        allocate (spec(npts))
        spec=cmplx(0.0d0,0.0d0)
        sumc = cmplx(0.0d0, 0.0d0)
        sums = cmplx(0.0d0, 0.0d0)
    end if

    if (allocated(ros)) deallocate (ros)
    allocate (ros(ns,ns))! contains the density matrix of S
    ros=cmplx(0.0d0,0.0d0)

    if (allocated(roi)) deallocate (roi)
    allocate (roi(ni,ni))! contains the density matrix of I
    roi=cmplx(0.0d0,0.0d0)

    ! preparation for phase cycling or coherence pathway (CTP) selection
    ! up to now, only the CTP selection is implemented

    if (allocated(ctp)) then
        nctp=size(ctp,1)
    else
        nctp=1
    end if

    allocate (roctp(nctp,nsi,nsi))
    roctp=0.0d0

    if (keepgoing) then
        !------------------------------------------------------------------
        ! make the recursive loop on the various k to calculate the spectra
        call loop_on_k(1_8)
        !------------------------------------------------------------------
        ! Last normalisation to unity area
        if (mod(ns,2) .eq. 0) then
            spec = spec /sumc    ! norm to 1 with respect to the central transition
        else
            spec = spec /sums
            ! or with respect to the whole transition for integer spin
        end if
    else
        exception=10
        exception_text="user cancellation"
    end if

    ! make the density matrix of S
    ros=0.0d0
    do i=1,ns
        do j=1,ns
            do k=0,ni-1
                ros(i,j)=ros(i,j)+rott(i+k*ns,j+k*ns)
            end do
        end do
    end do
    ros=ros/dble(ni)/fnexp

    ! make the density matrix of I
    if (icoupled>0.and.ni.GE.2) then
        roi=0.0d0
        do i=1,ni
            do j=1,ni
                do k=1,ns
                    roi(i,j)=roi(i,j)+rott(k+(i-1)*ns,k+(j-1)*ns)
                end do
            end do
        end do
        roi=roi*vls/vli/dble(ns)/fnexp
    end if

    ! deallocate (important if there is a new call to 'compute' subroutine by python)
    !--------------------------------------------------------------------------------
    deallocate (rfpowers)
    deallocate (rfoffsets)
    deallocate (rfphases)

    deallocate (rfpoweri)
    deallocate (rfoffseti)
    deallocate (rfphasei)

    deallocate (rflength)

    deallocate (rfdelay)
    deallocate (decouple)

    deallocate (itour)
    deallocate (ntp)
    deallocate (deltatp)
    deallocate (ntpadd)
    deallocate (deltatpadd)
    deallocate (t2ss)
    deallocate (t2ii)
    deallocate (t2si)
    deallocate (u1s,u2s)

    if (kcoupled>0) then
        deallocate (r2mdsk)
        deallocate (nk)
        deallocate (dipk)
        deallocate (mk)
        deallocate (adks)
        deallocate (spink)
        deallocate (vlk)
        deallocate (jsk)
        deallocate (dsk)
        deallocate (dskangle)
    end if

    deallocate (wr)
    deallocate (amp)
    deallocate (a)
    deallocate (p)
    deallocate (u)
    deallocate (ui)
    deallocate (uf)
    deallocate (usi)
    deallocate (usf)
    deallocate (uint)
    deallocate (ro)
    deallocate (rot)
    deallocate (z)
    deallocate (vp)
    deallocate (expvp)
    deallocate (y)
    deallocate (htst)
    deallocate (htsr)
    deallocate (sti)
    deallocate (rosto)
    deallocate (q)
    deallocate (rott)
    deallocate (qt)
    deallocate (roctp)

    !computation complete
    compute=.true.
    return

    !internal functions and subroutines===========================================
    contains

        ! loop on k ================================================================
        recursive subroutine loop_on_k(ik)
        !---------------------------------
            implicit none

            integer*8 :: ik ! index of the nucleus k
            integer*8 :: jk
            integer*8 :: k

            if (ik > kcoupled) then
                ! we are in the innermost loop
                ! execute the remaining of the program

                ! calculation of the chemical shift for each transition
                isost=isos
                do k=1,kcoupled
                    isost=isost+jsk(k)*dble(mk(k))/2.0d0
                end do

                wholesphere=.false.
                wholesphere=(wholesphere.or.&
                (abs(csas)*(abs(csangle(1))+abs(csangle(2))+abs(csangle(3))).gt.eps))
                wholesphere=(wholesphere.or.&
                (abs(cqs)*(abs(qsangle(1))+abs(qsangle(2))+abs(qsangle(3))).gt.eps))
                if (icoupled>0) wholesphere=(wholesphere.or.&
                    (abs(cqi)*(abs(qiangle(1))+abs(qiangle(2))+abs(qiangle(3))).gt.eps))
                if (icoupled>0) wholesphere=(wholesphere.or. &
                    (abs(dsiangle(1))+abs(dsiangle(2))+abs(dsiangle(3))).gt.eps)
                if (icoupled>0) wholesphere=(wholesphere.or. &
                    (abs(csai)*(abs(ciangle(1))+abs(ciangle(2))+abs(ciangle(3))).gt.eps))
                do k=1, kcoupled
                    wholesphere=(wholesphere.or.((abs(mk(k))*abs(dskangle(k,1))+ &
                                        abs(mk(k))*abs(dskangle(k,2))+ &
                                        abs(mk(k))*abs(dskangle(k,3))).gt.eps))
                end do


                ! reset to zero
                qt=0.0d0 ! is-it necessary here?      (in the loop?)

                ! todo: loop on pulse phases if there is a phase cycling
                ! -----------------------------------------------------
                ! todo: create a recursive function to make these loops

                if(.not.wholesphere) then
                    dnorm=abundance*dble(nt)*coefdip/dble(ng)**2/dble(ns)/dble(nscans)
                    call edgea()
                else
                    dnorm=abundance*dble(nt)*coefdip/dble(ng)**2 &
                                                           /dble(ns)/dble(nscans)/4.0d0
                    call edgeb()
                end if

                ! end of phases loops

                ! compute the final powder spectra with various integration
                ! depending on the wholesphere flag
                ! ----------------------------------------------------------

                ! if only the averaged density matrix values are needed (intensities),
                ! we can bypass the next step:

                if(.not.wholesphere) then
                    fnexp=dble(ng*nscans)*powder1/coefdip
                    if (spectrum_compute) call spectrea()
                else
                    fnexp=dble(ng*nscans)*powder2/coefdip
                    if (spectrum_compute)  call spectreb()
                end if
                return
            end if

            ! recursive call to the subroutine
            !---------------------------------
            ! replace the nested loops on the differents nucleus
            ! programming hint:
            ! do not use index instead of jk, ik that are defined outside this function

            do jk=-nk(ik),nk(ik),2
                mk(ik)=jk
                dipk(ik)=dsk(ik)*dble(jk)/2.0d0
                call loop_on_k(ik+1)
            end do

        end subroutine loop_on_k
        ! end loop on k ============================================================

        ! edgeA ====================================================================
        subroutine edgea()
            !-----------------
            implicit none
            integer*8 :: jj,ii,cr
            real*8    :: r,r3,xi,yi,zi

            cr=0
            do jj=0,nt
                do ii=0,nt-jj
                    r=sqrt(dble((nt-ii-jj)**2+ii*ii+jj*jj))
                    r3=r*r*r
                    xi=dble(nt-ii-jj)/r
                    yi=dble(ii)/r
                    zi=dble(jj)/r
                    call ll(xi,yi,zi,ii,jj,r3)
                    cr=cr+1
                end do
            end do
        end subroutine edgea
        !===========================================================================

        !===========================================================================
        subroutine edgeb()
        !   -----------------
            implicit none
            integer*8 :: jj,ii,cr
            real*8    :: r,r3,xi,yi,zi
            integer*8 :: nsbs
            !
            cr=0
            do jj=0,nt-1
                do ii=0,nt-jj
                    r=sqrt(dble((nt-ii-jj)**2+ii*ii+jj*jj))
                    r3=r*r*r
                    xi=dble(nt-ii-jj)/r
                    yi=dble(ii)/r
                    zi=dble(jj)/r
                    call ll(xi,yi,zi,ii,jj,r3)
                    cr=cr+1
                end do
                do ii=nt-jj+1,nt
                    r=sqrt(dble((nt-ii-jj)**2+(nt-jj)**2+(nt-ii)**2))
                    r3=r*r*r
                    xi=dble(nt-ii-jj)/r
                    yi=dble(nt-jj)/r
                    zi=dble(nt-ii)/r
                    call ll(xi,yi,zi,ii,jj,r3)
                    cr=cr+1
                end do
            end do
            !
            do jj=nt,2*nt-1
                do ii=jj-nt+1,nt-1
                    r=sqrt(dble((jj-nt-ii)**2+(nt-jj)**2+(nt-ii)**2))
                    r3=r*r*r
                    xi=dble(jj-nt-ii)/r
                    yi=dble(nt-jj)/r
                    zi=dble(nt-ii)/r
                    call ll(xi,yi,zi,ii,jj,r3)
                    cr=cr+1
                end do
                do ii=1,jj-nt
                    r=sqrt(dble((jj-nt-ii)**2+ii*ii+(2*nt-jj)**2))
                    r3=r*r*r
                    xi=dble(jj-nt-ii)/r
                    yi=-dble(ii)/r
                    zi=dble(2*nt-jj)/r
                    call ll(xi,yi,zi,ii,jj,r3)
                    cr=cr+1
                end do
            end do
            !
            do k=1,ns-1
                do jj=0,nt-1
                    wr(0,2*nt-jj,k)=wr(0,jj,k)
                    do nsbs=-nsb,nsb
                        amp(0,2*nt-jj,k,nsbs)=amp(0,jj,k,nsbs)
                    end do
                end do
                do ii=0,nt
                    wr(nt,nt+ii,k)=wr(ii,0,k)
                    do nsbs=-nsb,nsb
                        amp(nt,nt+ii,k,nsbs)=amp(ii,0,k,nsbs)
                    end do
                end do
                do jj=1,nt-1
                    wr(nt-jj,2*nt,k)=wr(nt,jj,k)
                    do nsbs=-nsb,nsb
                        amp(nt-jj,2*nt,k,nsbs)=amp(nt,jj,k,nsbs)
                    end do
                end do
            end do
            !
            r3=dble(nt*nt*nt)
            call ll(0.0d0,0.0d0,1.0d0,0_8,nt,r3)
            !
        end subroutine edgeb
        !===========================================================================

        !===========================================================================
        subroutine spectrea()
            !   --------------------
            implicit none
            integer*8 :: nsbs
            real*8 :: decal
            complex*8 :: ampt

            do k=1,ns-1
                !
                if(nall.eq.1.and.(2*k).ne.ns) cycle
                if(nall.eq.2.and.(2*k).eq.ns) cycle
                !
                do nsbs=-nsb,nsb
                    decal=nsbs*spinningspeed
                    do i=0,nt-2
                        do j=0,nt-2-i
                            ampt=amp(i+1,j,k,nsbs)+amp(i,j+1,k,nsbs)+amp(i,j,k,nsbs)
                            if ((2*k).eq.ns) then
                                sumc = sumc + ampt
                            else
                                sums = sums + ampt
                            end if
                            call tent(wr(i+1,j,k),wr(i,j+1,k),wr(i,j,k),ampt,decal)
                            ampt=amp(i+1,j,k,nsbs)+amp(i,j+1,k,nsbs)+amp(i+1,j+1,k,nsbs)
                            if ((2*k).eq.ns) then
                                sumc = sumc + ampt
                            else
                                sums = sums + ampt
                            end if
                            call tent(wr(i+1,j,k),wr(i,j+1,k),wr(i+1,j+1,k),ampt,decal)
                        end do
                    end do
                    do i=0,nt-1
                        j=nt-1-i
                        ampt=amp(i+1,j,k,nsbs)+amp(i,j+1,k,nsbs)+amp(i,j,k,nsbs)
                        if ((2*k).eq.ns) then
                            sumc = sumc + ampt
                        else
                            sums = sums + ampt
                        end if
                        call tent(wr(i+1,j,k),wr(i,j+1,k),wr(i,j,k),ampt,decal)
                    end do
                end do
            end do
        end subroutine spectrea
        !===========================================================================

        !===========================================================================
        subroutine spectreb()
        !   --------------------
            implicit none
            !
            integer*8 :: nsbs
            real*8 :: decal
            complex*8 :: ampt
            !
            do k=1,ns-1
                !
                if(nall.eq.1.and.(2*k).ne.ns) cycle
                if(nall.eq.2.and.(2*k).eq.ns) cycle
                !
                do nsbs=-nsb,nsb
                    decal=nsbs*spinningspeed
                    do i=0,nt-1
                        do j=0,nt-1
                            ampt=amp(i+1,j,k,nsbs)+amp(i,j+1,k,nsbs)+amp(i,j,k,nsbs)
                            if ((2*k).eq.ns) then
                                sumc = sumc + ampt
                            else
                                sums = sums + ampt
                            end if
                            call tent(wr(i+1,j,k),wr(i,j+1,k),wr(i,j,k),ampt,decal)
                            ampt=amp(i+1,j,k,nsbs)+amp(i,j+1,k,nsbs)+amp(i+1,j+1,k,nsbs)
                            if ((2*k).eq.ns) then
                                sumc = sumc + ampt
                            else
                                sums = sums + ampt
                            end if
                            call tent(wr(i+1,j,k),wr(i,j+1,k),wr(i+1,j+1,k),ampt,decal)
                        end do
                        do  j=nt,2*nt-1
                            ampt=amp(i,j,k,nsbs)+amp(i+1,j+1,k,nsbs)+amp(i+1,j,k,nsbs)
                            if ((2*k).eq.ns) then
                                sumc = sumc + ampt
                            else
                                sums = sums + ampt
                            end if
                            call tent(wr(i,j,k),wr(i+1,j+1,k),wr(i+1,j,k),ampt,decal)
                            ampt=amp(i,j,k,nsbs)+amp(i+1,j+1,k,nsbs)+amp(i,j+1,k,nsbs)
                            if ((2*k).eq.ns) then
                                sumc = sumc + ampt
                            else
                                sums = sums + ampt
                            end if
                            call tent(wr(i,j,k),wr(i+1,j+1,k),wr(i,j+1,k),ampt,decal)
                        end do
                    end do
                end do
            end do
        end subroutine spectreb
        !===========================================================================

        !===========================================================================
        subroutine tent(w1,w2,w3,amplitude,shift)
            !   ----------------------------------------
            !   todo: - rewrite to suppress gotos and labels
            implicit none
            real*8, intent(in) :: w1
            real*8, intent(in) :: w2
            real*8, intent(in) :: w3
            complex*8, intent(in) :: amplitude
            real*8, intent(in) :: shift
            real*8 :: finc
            real*8 :: ff1
            real*8 :: ff2
            real*8 :: ff3
            real*8 :: fmin
            real*8 :: fmid
            real*8 :: fmax
            real*8 :: fdn
            real*8 :: fxd
            real*8 :: fmn
            complex*8 :: top
            integer*8 :: np
            integer*8 :: npmid
            integer*8 :: npmax
            real*8 :: f1
            real*8 :: f2
            finc=sw/dble(npts-1)
            ff1=max(w1,w2)
            ff2=max(w2,w3)
            ff3=max(w3,w1)
            fmin=min(w1,w2,w3)+shift
            fmid=min(ff1,ff2,ff3)+shift
            fmax=max(w1,w2,w3)+shift
            fdn=fmid-fmin+0.0000001d0
            fxd=fmax-fmid+0.0000001d0
            fmn=fmax-fmin+0.0000001d0
            top=amplitude/fmn
            np=int((fmin-fst)/finc)+1
            npmid=int((fmid-fst)/finc)+1
            npmax=int((fmax-fst)/finc)+1
            if(npmax.gt.npts.or.np.lt.1) goto 80
            if(np.ne.npmid) goto 10
            spec(np)=spec(np)+fdn*top
            goto 40
        10  f2=finc*dble(np)+fst
            spec(np)=spec(np)+(f2-fmin)**2*top/fdn
        20  np=np+1
            f1=f2
            if(np.eq.npmid) goto 30
            f2=finc*dble(np)+fst
            spec(np)=spec(np)+finc*(f2+f1-2.0d0*fmin)*top/fdn
            goto 20
        30  spec(np)=spec(np)+(fmid-f1)*(fdn+f1-fmin)*top/fdn
        40  if(np.ne.npmax) goto 50
            spec(np)=spec(np)+fxd*top
            goto 80
        50  f2=finc*dble(npmid)+fst
            spec(np)=spec(np)+(f2-fmid)*(fmax-f2+fxd)*top/fxd
        60  np=np+1
            f1=f2
            if(np.eq.npmax) goto 70
            f2=finc*dble(np)+fst
            spec(np)=spec(np)+finc*(2.0d0*fmax-f1-f2)*top/fxd
            goto 60
        70  spec(np)=spec(np)+(fmax-f1)**2*top/fxd
        80  continue
        end subroutine tent
        !=================================================================================
    !end contains
end function compute
!=========================================================================================

