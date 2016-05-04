! samples the partition function of a square ice ising model H=J\sum_v Q_v, where Q_v=\sum_{i \in v} \sigma^{z}_i
! using single spin flips and short-loop updates 


! Lattice variables
module lattice
save 
integer(4) lx,ly,nh,ns ! system's dimensions
integer(4), allocatable :: neig(:,:)  ! table of neighbours
integer(4), allocatable :: vertex_to_spin(:,:) ! given a vertex, what spins belong to it
integer(4), allocatable :: conn(:,:) ! connections
integer(4), allocatable :: ordered(:),visited(:,:) ! loop update variables 
integer(4), allocatable :: sub(:) 
end module lattice

! spin configuration variables
module configuration
save
integer(4), allocatable :: spin(:),spint(:) !spin configuration
integer(4),allocatable :: Qc(:) ! charge Q at vertex c
real(8) beta ! temperature
real(8) En ! 
end module configuration

module measurements
save
integer(4), allocatable :: sq(:,:)
real(8), allocatable :: cost(:),sint(:)
real(8),allocatable :: qv(:,:)
end module measurements

program squareice
use lattice
use configuration
implicit none
integer(4) i,j,k,m,nbinsThermal,nbinsProduction,binl,went,countervisits,avvisit,cloop,nloops ! monte carlo: details of the sampling
integer(4)iseedd
real(8) :: drand1 ! random number
real(8) :: ave 
real(8), external :: energy

iseedd=0

!input
read(5,*)lx ! size of the system along x
read(5,*)ly ! size of the system along y
read(5,*)beta ! inverse temperature
read(5,*)iseedd ! seed for random number generetor
read(5,*)nbinsThermal ! number of bins in the thermalization
read(5,*)nbinsProduction ! number of bins during production run
read(5,*)binl ! samples per bin
read(5,*)nloops !number of loops per update


!initialization
call glattice() ! generating lattice
if (iseedd==0) call ransystem(iseedd) ! initializaitng random number generator
call rand_init(iseedd) ! initialize random number generator
call init() ! initial spin configuration
En=energy() ! initial energy
!write(6,*)'Energy init', En

!thermalization
do k=1,nbinsThermal
 do m=1,binl
 
  do i=1,nh 
   call spinflip() ! single spin flips
  end do
 
  cloop=0
  avvisit=0 
  
  do i=1,nloops
   call loopupdate(went,countervisits) ! loop updates
   if(went==1)then
      cloop=cloop+1
       avvisit=avvisit+countervisits ! counting visits
   end if
  end do
 
  ! adjusting the number of loops so that the updates are effective and decorrelate samples effectively 
  if(cloop>0)then
   ave=dble( avvisit)/dble(cloop)
   if(nloops<ns/(2*ave))then
      nloops=nloops+1
   elseif(nloops>ns/(2*ave))then
     nloops=nloops-1
     if(nloops<1)nloops=1
   end if  
  end if

 
 end do
 !write(6,*)'diffE',En-energy()  
end do

!write(6,*)'number of loops',nloops


! production
do k=1,nbinsProduction
 do m=1,binl

  do i=1,nh
   call spinflip()
  end do
 
  do i=1,nloops
   call loopupdate(went,countervisits)
  ! write(6,*)'E',energy() 
  end do
  
 end do
 write(6,*)spin 
 write(16,*)'diffE',energy()
 !open(10,file='1config.dat',status='replace') 
 !write(10,*)spin(:) 
 !close(10)
 !CALL SYSTEM('ipython activation.py')
 
  
end do


end program squareice


!-- loop update  ---!

subroutine loopupdate(went,countervisits)
use configuration
use lattice
implicit none
integer(4)i,itetra,ctetra,countervisits,spin2(2),outs,go,vstart,vfinal,tet,toflip,went
real(8) :: drand1
integer(4), external :: charge

! initializing variables tracking which vertices are visited and in what order.
ordered=-1
visited=-1
countervisits = 1  
go=0 
 
!! select a random vertex to start the loop
itetra=drand1()*ns+1
ctetra=itetra ! current vertex

!write(6,*)'loop construction'
loop: do 
     !write(6,*)'visiting tetrahedron:', ctetra 
     if(visited(ctetra,1)==-1)then !ctetra has not been visited

       visited(ctetra,1)=1 ! mark ctetra as visited
       ordered(countervisits)=ctetra ! ctetra vertex is visited at countervisits visit
       visited(ctetra,2)=countervisits ! ctetra is visited at countervisits visit
       countervisits=countervisits+1
       if(charge(ctetra)==0) then ! if the vertex is charge free
           call inoutspin(ctetra,spin2,countervisits)  ! finds which are the possible spins pointing out of the vertex
           !write(6,*)'spins in ctetra',vertex_to_spin(ctetra,:),spin(vertex_to_spin(ctetra,:))
           outs=drand1()*2+1 ! picks one at random
           visited(ctetra,3)=spin2(outs) ! marks what spin was used as output in vertex ctetra 
           if (conn(spin2(outs),1)==ctetra) then
              ctetra=conn(spin2(outs),2)  ! pick the next vertex to be visited
           elseif (conn(spin2(outs),2)==ctetra)then
              ctetra=conn(spin2(outs),1)
           end if
       else ! charge is present at ctetra vertex  
        go=-1 
        exit loop
       end if

     elseif(visited(ctetra,1)==1)then ! the vertex has been visited before
       !! loop found
       vstart=visited(ctetra,2)
       vfinal=countervisits-1
       go=1
       !write(6,*)'loop found',vstart,vfinal
       exit loop  
     end if 
end do loop

spint=spin
! flipping the spins along the loop if there was one successful loop
if (go==1)then
  do i=vstart,vfinal
    tet=ordered(i)
    toflip=visited(tet,3)
    !write(6,*)'flipping spin', toflip 
    spin(toflip)=-spin(toflip)  
  end do
  !write(6,*)'did spin change?',sum(abs(spin-spint)) 
end if  

went=go
!    visits=countervisits;
  

end subroutine loopupdate

subroutine inoutspin(ctetra,spin2,x)
use configuration
use lattice
implicit none
integer(4)ctetra,spin2(2),coun,i,x

spin2=0
coun=1
do i=1,4
if(spin(vertex_to_spin(ctetra,i))==(-1)**x)then
   spin2(coun)=vertex_to_spin(ctetra,i)
   coun=coun+1
end if

end do

!spin2=0
!coun=1
!do i=1,4
!if(i==1.or.i==2)then
!  if(spin(vertex_to_spin(ctetra,i))==-1)then
!   spin2(coun)=vertex_to_spin(ctetra,i)
!   coun=coun+1 
!  end if
!elseif(i==3.or.i==4)then
!  if(spin(vertex_to_spin(ctetra,i))==1)then
!   spin2(coun)=vertex_to_spin(ctetra,i)
!   coun=coun+1
!  end if
!end if

!end do
!!write(6,*)'possible output spin', ctetra, spin2


end subroutine inoutspin

! single spin flip
subroutine spinflip()
use configuration
use lattice
implicit none
integer(4)i,ri,t1,t2,q1,q2,q1p,q2p
real(8) :: drand1,ediff,boltz,prob
integer(4), external :: charge

!do i=1,1*nh
  ri=drand1()*nh+1
  q1=charge(conn(ri,1))
  q2=charge(conn(ri,2))
  spin(ri)=-spin(ri) 
  q1p=charge(conn(ri,1))
  q2p=charge(conn(ri,2))
  spin(ri)=-spin(ri)
  ediff=q1p**2+q2p**2-q1**2-q2**2
  if(ediff<=0)then
    spin(ri)=-spin(ri)
    Qc(conn(ri,1))=q1p
    Qc(conn(ri,2))=q2p 
    En=En+ediff  
  else
    prob=drand1()
    boltz=exp(-ediff*beta)
    if(boltz>prob)then
     spin(ri)=-spin(ri)
     Qc(conn(ri,1))=q1p
     Qc(conn(ri,2))=q2p
     En=En+ediff
    end if
  end if   
!end do

!write(6,*)'Energy',En/ns
end subroutine spinflip
! random number seed from system
subroutine ransystem(iseedd)
implicit none
integer(4)iseedd
CALL SYSTEM('echo $RANDOM>random.dat')
open(unit=100,file='random.dat',form='formatted',status='unknown')
read(100,*) iseedd
!write(6,*)'seed',iseedd 
close(100)
CALL SYSTEM('rm random.dat')
end subroutine ransystem

!-----Energy-------!
!
FUNCTION energy()
USE lattice
IMPLICIT NONE
REAL(8) :: energy
integer(4), external :: charge
integer(4)i,j 

energy=0.0d0
do i=1,ns
energy=energy+(charge(i))**2 
end do
!energy=energy/dble(ns)
return
END FUNCTION energy

!-----Charge-------!
function charge(v)
USE configuration
USE lattice
integer(4)i,v
integer(4)charge
charge=0
do i=1,4
  charge=charge+spin(vertex_to_spin(v,i))
end do
!do i=1,4
! if(i==1.or.i==2)then 
!  charge=charge-spin(vertex_to_spin(v,i))
! elseif(i==3.or.i==4)then
!  charge=charge+spin(vertex_to_spin(v,i))
! end if 
!end do
return
END FUNCTION charge

subroutine init()
use configuration
use lattice
implicit none
integer(4)i,j
real(8) :: drand1
integer(4), external :: charge

allocate(spin(nh),spint(nh),Qc(ns))
do j=1,10*nh
do i=1,nh
 ! random initialization
  spin(i)=2*int(2.*drand1())-1
end do
end do

do i=1,nh
  ! AF and charge free configuration ( if we have PBC, which we do)
  spin(i)=(-1)**(i+1)
end do

do i=1,ns
 Qc(i)=charge(i)
 !write(6,*)'charges',i,Qc(i)
end do

end subroutine init



subroutine glattice()
use lattice
use measurements
implicit none
integer(4)i,j,k,tet

ns=lx*ly
nh=2*lx*ly
allocate(neig(ns,4),vertex_to_spin(ns,4),conn(nh,2),ordered(ns), visited(ns,4),sub(ns))

! generating nearest neighbors
neig=0
k=1
do j=1,ly
 do i=1,lx

  if(i<lx)then
   neig(k,1)=k+1
  elseif(i==lx)then
   neig(k,1)=k-(lx-1)
  end if 
  if(j<ly)then  
   neig(k,2)=k+lx
  elseif(j==ly)then
   neig(k,2)=k-lx*(ly-1)
  end if
  if(i>1)then
   neig(k,3)=k-1
  elseif(i==1)then 
   neig(k,3)=k+(lx-1)  
  end if
  if(j>1)then
   neig(k,4)=k-lx
  elseif(j==1)then
   neig(k,4)=k+lx*(ly-1) 
  end if 
  k=k+1
 end do
end do

!write(6,*)"Neighbors"
!do i=1,ns
! write(6,*)i,neig(i,:)
!end do


!write(6,*)"vertices"
! generating vertex to spin table
do i=1,ns
 vertex_to_spin(i,1)=2*(i-1)+1
 vertex_to_spin(i,2)=2*(i-1)+2
 vertex_to_spin(i,3)=2*(neig(i,3)-1)+1
 vertex_to_spin(i,4)=2*(neig(i,4)-1)+2  
 !write(6,*)i-1,vertex_to_spin(i,:)-1
end do

!write(6,*)"tetrahedra sharing spin i"
! generating connection between spin and which vertex
do i=1,nh
   if (mod(i,2)==1) then
     tet=i/2+1 
     conn(i,1)=tet
     conn(i,2)=neig(tet,1)
   elseif(mod(i,2)==0)then
     tet=i/2  
     conn(i,1)=tet
     conn(i,2)=neig(tet,2) 
   end if
 !write(6,*)i,conn(i,:)
end do

tet=1
do j=1,ly
 do i=1,lx
  sub(tet)=(-1)**(i+j)
  tet=tet+1 
 end do
end do
!do i=1,ns
!write(6,*)i,sub(i)
!end do
!stop

allocate(sq(ns,2),qv(ns,2),cost(ns),sint(ns)) ! 2 sublattices

end subroutine glattice
