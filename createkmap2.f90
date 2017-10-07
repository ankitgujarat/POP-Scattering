subroutine createkmap(kp,xqq,nk1,nk2,nk3,nks,kmp)
!Subroutine to find the mapping of the k+q to the k point grid. For running this code make sure
!that the q vectors are near the gamma point which can be achieved by calculating phonons only 
!near the gamma point within 0.25 of the distance as stated in the literature
!**********************************************************************************************
implicit none
integer, intent(in) :: nk1,nk2,nk3,nks
real,intent(in) :: kp(:,:)
real,intent(in) :: xqq(:)
integer :: kmp(:)
!
!Local variables
!
real,allocatable ::xkq(:,:)
real :: xx,yy,zz,e,xk,yk,zk
integer :: ik ! index for the k vectors
integer :: n !Looping through all the k points to find the index of the k+q point
logical :: comm_q,in_grid,ind_fnd !index of the k+q point found
!allocate(kp(3,nks),kmap(nks))
!allocate(xqq(3))
allocate(xkq(3,nks))
e = 1.d-4
do ik=1,nks
    xx = kp(1,ik)*nk1
    yy = kp(2,ik)*nk2
    zz = kp(3,ik)*nk3
!The k point grid should also be positive
    if(xx .le. -e .or. yy .le. -e .or. zz .le. -e) then
        print *,'The k point grid is negative hence terminating calculation'
    end if
!Checking the commensuration of q point with respect to the k point grid because if the q point grid is not commensurate with the k point
!grid then the addition of the k and the q point will never lie in the k point grid thus making the calculation rather useless
!
!We are not converting any quantity from cartesian to crystal coordinate since the input file is in the crystal co-ordinate
    xx = xqq(1)*nk1
    yy = xqq(2)*nk2
    zz = xqq(3)*nk3
    comm_q = abs(xx-nint(xx)).le.e .and. &
             abs(yy-nint(yy)).le.e .and. &
             abs(zz-nint(zz)).le.e 
    if (.not.comm_q) then
        print *,'The q point grid is not commensurate with the k point grid hence terminating the calculation'
    end if
!Adding the q point to the k points and checking whether the k+q point lies on the k point grid
    xkq(:,ik) = kp(:,ik)+xqq(:)
    xx = xkq(1,ik)*nk1
    yy = xkq(2,ik)*nk2
    zz = xkq(3,ik)*nk3
    in_grid = abs(xx-nint(xx)).le.e .and. &
              abs(yy-nint(yy)).le.e .and. &
              abs(zz-nint(zz)).le.e
    if(.not.in_grid) then
        print *,'The k+q point doesnot lie on the k grid'
    end if
!If the k+q point is found to be on the k grid, then the next task is to find the index of the k+q point in the
!k grid. Since the q points are near the Gamma point, then I donot have to think about the G vector calculation
!since I donot expect any type of folding necessary 
!Here the choice of k and q points will be such that substracting 1 from the xx,yy and zz value will not be negative. It is like folding back to the
!Brillouin zone
    if (xx > (1- 1/nk1)*nk1) then
        xx = xx - 1
    end if
    if (yy > (1-1/nk2)*nk2) then
        yy = yy - 1
    end if
    if (zz > (1-1/nk3)*nk3) then
        zz = zz - 1
    end if
    if(nint(xx) == nk1) then
        xx = 0
    end if
    if(nint(yy) == nk2) then
        yy = 0
    end if
    if(nint(zz) == nk3) then
        zz = 0
    end if
!The idea implemented by the EPW team is to loop through all the k points in the grid and find the one that matches
!the newlt generated k+q point which is a very trivial way of doing it
    do n=1,nks
        xk = kp(1,n)*nk1
        yk = kp(2,n)*nk2
        zk = kp(3,n)*nk3
        ind_fnd = nint(xx) .eq. nint(xk) .and. &
                  nint(yy) .eq. nint(yk) .and. &
                  nint(zz) .eq. nint(zk) 
        if (ind_fnd) then
            kmp(ik) = n
        else
            print *,'Cannot find the index, something not working in the code above'
        end if
    end do
end do
end subroutine createkmap
!*************************************************************************************
!A subroutine can modify the arguments passed to it
!*************************************************************************************
!Check for the Brillouin Zone folding :- Why is it necessary???
!Got the reason for folding back to the Brillouin zone
!*************************************************************************************