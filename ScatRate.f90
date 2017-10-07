program ScatRate
implicit none
!Declaring Variables
integer :: ios,nks,nmode,kp,mo,kf
real,allocatable :: kps(:,:),ener(:),mo_ener(:),scat_st(:),freq(:),q(:),scat_ra(:)
!Loop variables
integer :: i,j
real :: mo_ene,k1_ene,k2_ene,kf_ene,tolerance
!
!Reading the input kpoints and the energy
open(199,file='./kpt.in',status='old',iostat=ios,form='formatted')
if (ios/=0) then
  print *,'Error reading the k point input file'
end if
print *,'Reading k points: done'
read(199,*)
read(199,*)
read(199,*) nks
open(200,file='./ener.in',status='old',iostat=ios,form='formatted')
if (ios/=0) then
  print *,'Error reading energy input file'
end if
!Allocating the matrices
allocate(kps(3,nks),ener(nks),scat_ra(nks),q(3))
!Reading the k points and energy values that is generated from the wannier90 geninterp program
do i=1,nks
  read(199,*,iostat=ios) kps(:,i)
  if (ios/=0) then
    print *,'Error reading the input kpoints'
  end if
  read(200,*,iostat=ios) ener(i)
  if (ios/=0) then
    print *,'Error reading the input energy values'
  end if
  scat_ra(i) = 0
end do
print *,'K point and energy reading complete'
close(199)
close(200)
!Reading the inverse wavelength for each mode
nmode=30
allocate (scat_str(nmode))
open(199,file='./freq.in',status='old',iostat=ios,form='formatted')
if (ios/=0) then
  print *,'Error opening the mode energy file'
end if
do i=1,nmode
  read(199,*,iostat=ios) freq(i)
  if (ios/=0) then
    print *,'Error reading input energy per mode file'
  end if
  mo_ener(i) = 0
end do
read (199,*)
do i=1,nmode
  read(199,*,iostat=ios) scat_str(i)
  if (ios/=0) then
    print *,'Error reading the scattering strength/matrix element'
  end if
end do
close(199)
print *,'Frequency and matrix element reading complete'
!****************************************************************************************************
!Call a subroutine to convert the inverse of wavelength into energy (meV)
call conv_ener(freq,nmode,mo_ener)
!****************************************************************************************************
!In the computing strategy that I am implementing, we have to loop through the energy value twice to obtain the final state which is computationally too expensive
!Defining energy tolerance to be about 1meV (just a random guess)
tolerance = 1.0d-3
do mo=1,nmode
  mo_ene = mo_ener(mo)
  do kp=1,nks
    k1_ene = ener(kp)+mo_ene*10.0d-3
    k2_ene = ener(kp)-mo_ene*10.0d-3
    do ks=1,nks
      if ((abs(ener(ks)-k1_ene)) .lt. tolerance .or. (abs(ener(ks)-k2_ene)) .lt. tolerance)then
        print *,"Final state found"
        scat_ra(kp) = scat_ra(kp) + scat_str(mo)**2
      end if
    end do
  end do
end do 