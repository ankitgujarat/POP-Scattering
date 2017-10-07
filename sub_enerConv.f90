subroutine conv_ener(freq,nmode,mo_ener)
implicit none
!Local Variables
integer :: nmode,j
real :: mo_ener(nmode),freq(nmode)
!Conversion code assuming the frequency in the dynmat file is given in THz
do j = 1,nmode
    mo_ener(j) = (1240*freq(j)*10.0d15)/(3*10.0d17)
end do
end subroutine conv_ener