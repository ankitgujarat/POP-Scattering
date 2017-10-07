program umatrices
!This program reads the U matrices from the input file and computes the matrix multiplication (overlap calculation) for the 
!Uk and Uk+q matrices using the map generated in the createkmap2.f90 file. The inspiration is largely drawn from the EPW
!teams code but it is just implemented for the q points near the Gamma point
implicit none
complex,allocatable :: uk(:,:)
integer :: i,j,nks,num_wann,ios
real :: a,b
nks = 64
num_wann = 4
allocate(uk(nks,num_wann*num_wann))
open(199,file='G:\Charge Transport Theory\Fortran Codes\input1_text.in',status='old',iostat=ios,action='read',form='formatted')
if (ios/=0) then
    print *,'Cannot read the u matrices file'
end if
do i=1,4
  read(199,*)
end do
do i=1,nks
  do j=1,num_wann*num_wann
    read(199,*) a,b
    uk(i,j) = cmplx(a,b)
  end do
  if(i<nks) then 
    read(199,*)
    read(199,*)
  end if
end do
print *,uk
close(199)
end program umatrices
!*****************************************************************************************
!Uk+q will be calculated on the fly from the Uk matrix read from the input file
!*****************************************************************************************