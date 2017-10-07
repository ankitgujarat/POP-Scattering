program recip
implicit none
real,allocatable :: latvec(:,:),qp(:)
real,allocatable :: rlatvec(:,:)
real,allocatable :: a(:),b(:),c(:),mat(:,:)
real :: n,omega,vol
integer :: i,j,k,l,m,q
allocate (latvec(3,3),a(3),b(3),c(3),qp(3))
allocate (rlatvec(3,3),mat(3,3))
qp(1) = 0.75
qp(2) = 0.25
qp(3) = 0.50
do i=1,3
    do j=1,3
        read *,latvec(i,j)
    end do
end do
do i=1,3
    m=i
    q=i
    do j=1,3
        if (m+1 > 3) then
            m=3
            q=0
        end if
        a(j) = latvec(m,j)
        b(j) = latvec(q+1,j)
    end do
    do k=1,3
        do l=1,3
            n=a(k)*b(l)
            if (k .eq. l) then
                mat(k,l) = 0
            end if
            if (l-k==2 .or. l-k==-1) then
                mat(k,l) = (-1)*n
            else
                mat(k,l) = n
            end if
        end do
    end do
    c(1) = mat(2,3)+mat(3,2)
    c(2) = mat(1,3)+mat(3,1)
    c(3) = mat(1,2)+mat(2,1)
    rlatvec(i,1) = c(1)
    rlatvec(i,2) = c(2)
    rlatvec(i,3) = c(3)
end do
omega = 2*3.14/((b(1)*rlatvec(2,1)+b(2)*rlatvec(2,2)+b(3)*rlatvec(2,3)))
!print *,(1/omega)*2*3.14
vol = b(1)*rlatvec(2,1)+b(2)*rlatvec(2,2)+b(3)*rlatvec(2,3)
print *,vol
do i=1,3
  do j=1,3
    rlatvec(i,j)=rlatvec(i,j)*omega
    print *,rlatvec(i,j)
  end do
end do
qp(1) = qp(1)*(rlatvec(1,1)+rlatvec(2,1)+rlatvec(3,1))
qp(2) = qp(2)*(rlatvec(1,2)+rlatvec(2,2)+rlatvec(3,2))
qp(3) = qp(3)*(rlatvec(1,3)+rlatvec(2,3)+rlatvec(3,3))
print *,'********************************************************************'
print *,qp
deallocate (latvec,a,b,c,mat)
end program recip