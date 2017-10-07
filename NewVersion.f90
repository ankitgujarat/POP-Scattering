program POPScat
implicit none
integer :: ios,nks,i,nk1,nk2,nk3,count,nqs,j,num_wann,nmodes,nat,v,l !Donot consider irreducible k points
character (len=30) :: coordinates,calculation,latvector
real,allocatable :: kps(:,:),qps(:,:),epsil(:,:),zeff(:,:,:),atmpos(:,:),disp(:,:,:),latvec(:,:),epmatl(:),freq(:),rlatvec(:,:)
integer,allocatable :: kmap(:)
complex,allocatable :: umat(:,:)
real :: alat,fac
!
namelist / inputpop / &
    calculation
!
namelist / kqpoints / &
    coordinates,nks,nk1,nk2,nk3,nqs,num_wann,nat,alat,latvector
!
namelist / phonon / &
	nmodes
! opening the file
open(199,file='G:\Charge Transport Theory\Calculations\Dir_res\i_tx.in',status='old',iostat=ios,action='read',form='formatted')
if (ios/=0) then
  print *,'Error reading input file'
end if
print *,'Read input text'
open(200,file='G:\Charge Transport Theory\Calculations\Dir_res\d_p.in',status='old',iostat=ios,action='read',form='formatted')
  if (ios/=0) then
    print *,'Error reading input file'
end if
open(202,file='G:\Charge Transport Theory\Calculations\Dir_res\qptx1.kpt',status='old',iostat=ios,action='read',form='formatted')
if (ios/=0) then
  print *,'Error reading the fine grid q point'
end if
!print *,'Read displacement pattern'
! reading input file
read(199,nml=inputpop,iostat=ios)
if (ios/=0) then
    print *,'Error reading input namelist1 ',ios
end if
print *,calculation
!
read(199,nml=kqpoints,iostat=ios)
if (ios/=0) then
    print *,'Error reading input namelist2 ',ios
end if
!
print *,'Read kpoint parameters'
read(199,nml=phonon,iostat=ios)
if (ios/=0) then
  print *,'Error reading input namelist'
end if
print *,'Read phonon parameters'
!print *,coordinates,nks,nk1,nk2,nk3,nqs
allocate(kmap(nks),epsil(3,3),zeff(3,3,nat),atmpos(3,nat),freq(3*nat),rlatvec(3,3))
allocate (kps(3,nks),qps(3,nqs),disp(3,nat,3*nat),latvec(3,3))
allocate (umat(nks/4,num_wann**2),epmatl(nmodes))
print *,'Matrix allocation complete'
do i=1,nmodes
  epmatl(i) = 0
end do
read(199,*)
do i=1,3
  read(199,*,iostat=ios) epsil(:,i)
  if (ios/=0) then
    print *,'Error reading the dielectric matrix'
  end if
end do
print *,'Read dielectric tensor'
!print *,epsil
!print *,'***********************************************************************************************'
read(199,*)
do i=1,nat
  do j=1,3
    read(199,*,iostat=ios) zeff(:,j,i)
    if (ios/=0) then
      print *,'Error reading the effective charge'
    end if
  end do
end do
print *,'Read effective charge tensor'
read(199,*)
do i=1,nat
  read(199,*,iostat=ios) atmpos(:,i)
  if (ios/=0) then
    print *,'Error reading atomic positions'
  end if
end do
print *,'Read atomic positions'
read(199,*)
do i=1,3
  read(199,*,iostat=ios) latvec(:,i)
  if (ios/=0) then
    print *,'Error reading lattice vectors'
  end if
end do
print *,'Read lattice vectors'
read(199,*)
do i=1,nks
    read (199,*,iostat=ios) kps(:,i)
    kmap(i) = 0
    if (ios/=0) then
        print *,'Error reading k points'
    end if
end do
print *,'Read input k points'
!read(199,*)
!do i=1,nqs
!    read(199,*,iostat=ios) qps(:,i)
!    if (ios/=0) then
!        print *,'Error reading q points'
!    end if
!end do
!print *,'Read input q points'
do i=1,3*nat
  do j=1,nat
    read(200,*,iostat=ios)disp(:,j,i)
    if (ios/=0) then
      print *,'Error reading displacement patterns'
    end if
  end do
end do
print *,'Read displacement pattern'
read(200,*)
do i=1,3*nat
  read(200,*,iostat=ios) freq(i)
  if (ios/=0) then
    print *,'Error reading the frequency components'
  end if
end do
print *,'Read frequency'
do i=1,nks/4
  do j=1,num_wann**2
    umat(i,j) = 0
  end do
end do
do j=1,nat
  print *,'Calling reciprocal lattice vector for atm',j
  call recip(atmpos(:,j),latvec,latvector,alat)
 end do
!call umatrices(nks,num_wann,umat)
read(202,*)
do i=1,nqs
  read(202,*,iostat=ios) qps(:,i)
  if (ios/=0) then
    print *,'Error reading q point'
  end if
end do
close(202)
print *,'Read fine grid q point complete'
!************************************************************************************************
!Calculation of the prefactor
!************************************************************************************************
!fac = (4*3.142)/711.31
do i=1,nqs 
  !convert q vector from crystal to cartesian coordinate
  !Avoiding singularity in q by translating it by reciprocal lattice vector 
  if (qps(1,i)==0.0 .and. qps(2,i)==0.0 .and. qps(3,i)==0.0) then
    qps(1,1) = 0.000
    qps(2,1) = 2.177
    qps(3,1) = 0.588
  end if
  print *,'Calling reciproal lattice vector for q',i
  call recip(qps(:,i),latvec,latvector,alat)
  print *,'Call to POP Scattering calculation'
  call popelements(qps(:,i),zeff,epsil,atmpos,disp,freq,nat,epmatl)
end do
close(199)
close(200)
open(201,file='G:\Charge Transport Theory\Calculations\Dir_res\El-Ph.out',action='write',iostat=ios,status='new',form='formatted')
if (ios/=0) then
  print *,'Error writing output file'
end if
!write(201,*) epmatl
!close(201)
fac = 4*3.142/678.776
do v = 1,3*nat
  write(201,*) epmatl(v)*fac
end do
print *,'Program completed successfully'
do v = 1,3
  print *,epsil(1,v),epsil(2,v),epsil(3,v)
end do
end program POPScat
!
!*****************************************************************************************************
!*****************************************************************************************************
subroutine recip(qp,latvec,latvectyp,alat)
implicit none
real :: latvec(3,3),qp(3)
real,allocatable :: rlatvec(:,:)
real,allocatable :: a(:),b(:),c(:),mat(:,:)
real :: n,omega,alat,lat
character (len=20) :: latvectyp
integer :: lt,in,k,l,m,qpos
allocate (a(3),b(3),c(3),rlatvec(3,3))
allocate (mat(3,3))
print *,'Inside subroutine recip'
print *,'co ordinates under consideration',qp
lat=alat
if (latvectyp == 'cartesian') then
  print *,'inside IF loop'
  lat = 1.0
end if
do lt=1,3
    m=lt
    qpos=lt
    do in=1,3
        if (m+1 > 3) then
            m=3
            qpos=0
        end if
        a(in) = latvec(m,in)
        b(in) = latvec(qpos+1,in)
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
    rlatvec(lt,1) = c(1)
    rlatvec(lt,2) = c(2)
    rlatvec(lt,3) = c(3)
end do
omega = 2*3.14/((b(1)*rlatvec(2,1)+b(2)*rlatvec(2,2)+b(3)*rlatvec(2,3)))
do lt=1,3
  do in=1,3
    rlatvec(lt,in)=rlatvec(lt,in)*omega/lat
    !print *,rlatvec(i,j)
  end do
end do
qp(1) = qp(1)*(rlatvec(1,1)+rlatvec(2,1)+rlatvec(3,1))
qp(2) = qp(2)*(rlatvec(1,2)+rlatvec(2,2)+rlatvec(3,2))
qp(3) = qp(3)*(rlatvec(1,3)+rlatvec(2,3)+rlatvec(3,3))
print *,'After',qp
deallocate (a,b,c,mat)
print *,'End Subroutine'
end subroutine recip
!
!Reciprocal lattice vectors are in 1/bohr so are the q co-ordinates
!**********************************************************************************************************
!**********************************************************************************************************
subroutine popelements(q,zeff,epsil,atmpos,disp,freq,nat,epmatl)
!q: the q point under consideration
!zeff: the born effective charge
!epsil : the dielectric tensor
!atmpos : the atomic positions
!disp: the atomic displacement pattern
!freq: the frequency of each mode
!nat: no of atoms per unit cell in the crystal system
!latvectype: the type of lattice vector as input ('cartesian' or 'crystal')
implicit none
integer :: nat
real :: q(3),zeff(3,3,nat),epsil(3,3),atmpos(3,nat),disp(3,nat,3*nat),freq(3*nat),epmatl(3*nat)
!
!define local variables
integer :: mode,natm,lt
real :: qeq,arg,ener
real,allocatable :: qz(:)
complex :: z,qze,qze1,qze2,qze3
!
print *,'Inside subroutine pop'
allocate(qz(3))
qeq = q(1)*(epsil(1,1)*q(1) + epsil(2,1)*q(2) + epsil(3,1)*q(3)) + &
      q(2)*(epsil(1,2)*q(1) + epsil(2,2)*q(2) + epsil(3,2)*q(3)) + &
      q(3)*(epsil(1,3)*q(1) + epsil(2,3)*q(2) + epsil(3,3)*q(3))
do mode = 1,3*nat
  qze=0
  qze1=0
  qze2=0
  qze3=0
  if (freq(mode)> 0.2) then
    do natm = 1,nat
        qz(1) = 0
        qz(2) = 0
        qz(3) = 0
        !qze = 0
        arg = -2*3.142*(q(1)*atmpos(1,natm) + q(2)*atmpos(2,natm) + q(3)*atmpos(3,natm))
        z = cmplx(cos(arg),sin(arg))
        do lt=1,3
            qz(lt)=(q(1)*zeff(lt,1,natm))+(q(2)*zeff(lt,2,natm))+(q(3)*zeff(lt,3,natm))
        end do
        qze1 = qz(1)*disp(1,natm,mode) + qz(2)*disp(2,natm,mode) + qz(3)*disp(3,natm,mode)
        qze2 = qze1
        !Have removed the z term from the above line to check the consistency of the code, but one should include 
        !the complex z term
        qze = qze+qze2
    end do
    !Cannot figure out the origin of the term exp(-qeq/4) 
    ener = ((2*(6.626*10d-34)*(freq(mode)*10d12)*(6.242*10d18))/27.211)**(-0.5) !in atomic units
    !ener = ener*4*3.142/678.776
    qze3 = (qze/qeq)*ener*cmplx(1)*exp(-qeq/4)
    epmatl(mode) = epmatl(mode)+abs(qze3)
  else
    epmatl(mode) = epmatl(mode)+0.0
  end if
end do
end subroutine popelements

!There is still some discrepency with the evaluation of constants terms which is a concern to get definite scattering rates 
