! Ising2d.f95 :Ising 2D model of magnetic dipole
!		Author: Prajwal Mohanmurthy
!		        Mississippi State University
!               04/2012

program ising2d
  implicit none

  integer ::  i,j,maink,cout,maxn,inte,mag
  integer,parameter :: n=20

  real(kind=8) :: ediff,ran,tt,ienergy,energy,iech,energy2,cv,cor
  integer(kind=8),dimension(n,n) :: s

 do maink = 1,101
  tt = 1.0 + 0.1*(maink-1)
  ediff = 0.0
  ienergy = 0.0
  iech = 0.0
  energy2 = 0.0
  cv = 0.0
  mag = 0
  tt = 4.0
  maxn = n**6

  call init()
  call inien()
  !write(*,*) "Energy, Energy^2, CV:"
  !write(*,*) energy, energy2
  !call printmat(s,n,n)

  call isingsim()
  call magnetization()
  !write(*,*) "--------------------------------------"
  !write(*,*) "Energy, Energy^2, CV, Magnetization:"
  write(*,*) tt, energy,energy2,cv,abs(mag)
  !call printmat(s,n,n)
 enddo

!--------------------------------------------------------------------------

contains

  subroutine isingsim()
    integer :: ii, ij
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

    subroutine magnetization()

        integer :: mi, mj
        mag = 0

        do mi = 1,n
            do mj = 1,n
                mag = mag + s(mi,mj)
            enddo
        enddo
    end subroutine magnetization

!--------------------------------------------------------------------------

  subroutine init()
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

    integer :: dt,db,dl,dr,dtr,dbl,ini,inj

    dt= 0
    db= 0
    dl= 0
    dr= 0
    dtr = 0
    dbl = 0
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
            if(inj==1 .and. ini==n) then
                dtr = s(1,n)
            else if(ini==n) then
                dtr = s(1,inj-1)
            else if(inj==1) then
                dtr = s(ini+1,n)
            else
                dtr = s(ini+1,inj-1)
            endif
            if(inj==n .and. ini==1) then
                dbl = s(n,1)
            else if(ini==1) then
                dbl = s(n,inj+1)
            else if(inj==n) then
                dbl = s(ini-1,1)
            else
                dbl = s(ini-1,inj+1)
            endif
            !write(*,*) dt, dr, db, dl
            ienergy = ienergy + (-s(ini,inj)*(dt+dr+db+dl+dbl+dtr))
            !write(*,*) ienergy
       enddo
    enddo
    energy = ienergy/2
    energy2 = energy**2

  end subroutine inien

!--------------------------------------------------------------------------

  subroutine deltau(di,dj)

    integer,intent(in) :: di,dj
    integer :: dt,db,dl,dr,dtr,dbl

    dt= 0
    db= 0
    dl= 0
    dr= 0
    dtr = 0
    dbl = 0
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
    if(dj==1 .and. di==n) then
        dtr = s(1,n)
    else if(di==n) then
        dtr = s(1,dj-1)
    else if(dj==1) then
        dtr = s(di+1,n)
    else
        dtr = s(di+1,dj-1)
    endif
    if(dj==n .and. di==1) then
        dbl = s(n,1)
    else if(di==1) then
        dbl = s(n,dj+1)
    else if(dj==n) then
        dbl = s(di-1,1)
    else
        dbl = s(di-1,dj+1)
    endif
    ediff = (2.0*s(di,dj)*(dt+dr+db+dl+dtr+dbl))

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
