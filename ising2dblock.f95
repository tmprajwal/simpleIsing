! Ising2d.f95 :Ising 2D model of magnetic dipole
!		Author: Prajwal Mohanmurthy
!		        Mississippi State University
!               04/2012

program ising2d
  implicit none

  integer ::  maini,mainj,maink,cout,inte,mag,mag2
  integer,parameter :: n=30
  integer,parameter :: n2=10
  integer(kind=8) :: maxn

  real(kind=8) :: ediff,ran,tt
  integer(kind=8),dimension(n,n) :: s
  integer(kind=8),dimension(n2,n2) :: s2

  maxn = n2**6
  maini=1

  !write(*,*) "--------------------------------------"

  do maini = 1,101
    tt = 1.0 + 0.1*(maini-1)
    !write(*,*) tt, "--------------------------------------"
    mag2 = 0
    mag = 0
    call init()
    call isingsim()
    !call printmat(s,n,n)
    !write(*,*) "--------------------------------------"
    call magnetization()
    call blocks()
    !call printmat(s2,n2,n2)
    call magnetization2()
    write(*,*) tt, abs(mag), abs(mag2)
  enddo

!--------------------------------------------------------------------------

contains

  subroutine isingsim()
    integer :: ii, ij,i,j
    real(kind=8) :: compran

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
        else if(compran < exp(-ediff/tt)) then
            inte = -s(ii,ij)
            s(ii,ij) = inte
        endif
    enddo

  end subroutine isingsim

!--------------------------------------------------------------------------

    subroutine blocks()

        integer :: si, sums,outi,outj,ini,inj
        sums = 0
        si = n/n2

        do outi=1,n2
            do outj=1,n2
                sums = 0
                do ini = ((outi-1)*si)+1,((outi-1)*si)+si
                    do inj = ((outj-1)*si)+1,((outj-1)*si)+si
                        sums = sums + s(ini,inj)
                     enddo
                enddo
                if(sums>0) then
                    s2(outi,outj) = 1
                else
                    s2(outi,outj) = -1
                endif
            enddo
        enddo

    end subroutine blocks

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

    subroutine magnetization2()

        integer :: mi, mj
        mag2 = 0

        do mi = 1,n2
            do mj = 1,n2
                mag2 = mag2 + s2(mi,mj)
            enddo
        enddo
    end subroutine magnetization2

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
