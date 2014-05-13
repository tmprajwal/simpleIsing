! Ising2d.f95 :Ising 2D model of magnetic dipole
!		Author: Prajwal Mohanmurthy
!		        Mississippi State University
!               04/2012

program ising2d
  implicit none

  integer ::  maini,mainj,maink,cout,inte,mag
  integer,parameter :: n=100
  integer(kind=8) :: maxn	

  real(kind=8) :: ediff,ran,tt,ienergy,energy,iech,energy2,cv,cor
  integer(kind=8),dimension(n,n) :: s

  ediff = 0.0
  ienergy = 0.0
  iech = 0.0
  energy2 = 0.0
  cv = 0.0
  mag = 0
  tt = 2.2
  maxn = 10**8

  write(*,*) "--------------------------------------"

  do maini = 1,11
    tt = 2.0 + 0.1*(maini-1)
    write(*,*) tt, "--------------------------------------"
    ediff = 0.0
    ienergy = 0.0
    iech = 0.0
    energy2 = 0.0
    cv = 0.0
    mag = 0
    call init()
    !call inien()
    call isingsim()
    do mainj = 1,n/2
        call corf((n/2) - mainj + 1)
        write(*,*) (n/2) - mainj + 1,abs(cor)
    enddo
  enddo



!--------------------------------------------------------------------------

contains

  subroutine isingsim()
    integer :: ii, ij,i,j
    real(kind=8) :: compran

    iech = 0.0
    do i = 1,maxn
        call random_number(ran)
        compran = ran
        call random_number(ran)
        ii = 1+floor(ran*n)
        call random_number(ran)
        ij = 1+floor(ran*n)
        !write(*,*) ii,ij

        call deltau(ii,ij)
        !write(*,*) ediff
        if(ediff<0) then
            inte = -s(ii,ij)
            s(ii,ij) = inte
            iech = iech + ediff
        else if(compran < exp(-ediff/tt)) then
            inte = -s(ii,ij)
            s(ii,ij) = inte
            iech = iech + ediff
        endif
        energy = energy + iech
        energy2 = energy2 + (iech*iech)
    enddo
    energy = energy/maxn
    energy2 = energy2/maxn
    cv = (energy2 - (energy**2))/(tt**2)

  end subroutine isingsim

!--------------------------------------------------------------------------

    subroutine corf(r)

        integer,intent(in) :: r
        integer :: ci,cj,ccou
        real(kind=8) :: sumsp, sp

        ccou = 0
        sumsp = 0.0
        sp = 0.0

        do ci = 1,n-r
            do cj = 1,n-r
                sumsp = sumsp + (s(ci,cj)*s(ci,cj+r))
                sp = sp + s(ci,cj)
                ccou = ccou + 1
                sumsp = sumsp + (s(ci,cj)*s(ci+r,cj))
                sp = sp + s(ci,cj)
                ccou = ccou + 1
            enddo
        enddo
        cor = ((sumsp/ccou) - ((2*sp/ccou)**2))

    end subroutine corf

!--------------------------------------------------------------------------

  subroutine init()
  integer :: i,j  
  do i = 1,n
        do j = 1,n
            ran = 0
            call random_number(ran)
            !write(*,*) ran
            if(ran >= 0.5) then
                s(i,j) = 1
                !write(*,*) s(i,j), ran
            else
                s(i,j) = -1
                !write(*,*) s(i,j), ran
            endif
        enddo
    enddo
  end subroutine init

!--------------------------------------------------------------------------

  subroutine inien()

    integer :: dt,db,dl,dr,ini,inj

    dt= 0
    db= 0
    dl= 0
    dr= 0
    ienergy = 0.0

    do ini = 1,n
        do inj = 1,n
            if(ini==1) then
                dt = s(n,inj)
            else
                dt = s(ini-1,inj)
            endif
            if(ini==n) then
                db = s(1,inj)
            else
                db = s(ini+1,inj)
            endif
            if(inj==1) then
                dl = s(ini,n)
            else
                dl = s(ini,inj-1)
            endif
            if(inj==n) then
                dr = s(ini,1)
            else
                dr = s(ini,inj+1)
            endif
            !write(*,*) dt, dr, db, dl
            ienergy = ienergy + (-s(ini,inj)*(dt+dr+db+dl))
            !write(*,*) ienergy
       enddo
    enddo
    energy = ienergy/2
    energy2 = energy**2

  end subroutine inien

!--------------------------------------------------------------------------

  subroutine deltau(di,dj)

    integer,intent(in) :: di,dj
    integer :: dt,db,dl,dr

    dt= 0
    db= 0
    dl= 0
    dr= 0
    ediff = 0.0

    if(di==1) then
        dt = s(n,dj)
    else
        dt = s(di-1,dj)
    endif
    if(di==n) then
        db = s(1,dj)
    else
        db = s(di+1,dj)
    endif
    if(dj==1) then
        dl = s(di,n)
    else
        dl = s(di,dj-1)
    endif
    if(dj==n) then
        dr = s(di,1)
    else
        dr = s(di,dj+1)
    endif
    !write(*,*) dt, dr, db, dl
    ediff = (2.0*s(di,dj)*(dt+dr+db+dl))

  end subroutine deltau

!--------------------------------------------------------------------------

  subroutine printmat(a,n,m)
  ! print out a n x m matrix
    integer(kind=8),intent(in) :: a(:,:)
    integer,intent(in) :: n,m

    integer i,j

    do i=1,n
       !write(*,"(I3,$)") i
       do j=1,m
          write(*,"(I3,$)") a(i,j)
       enddo
       write(*,*)
    enddo

  end subroutine printmat

!--------------------------------------------------------------------------

end program ising2d

