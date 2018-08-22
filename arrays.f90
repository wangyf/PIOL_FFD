! Allocate arrays
module m_arrays
implicit none
contains

subroutine arrays
use m_globals
integer :: j, k, l

! 3d
j = nm(1)
k = nm(2)
l = nm(3)


allocate(         &
    s1(j,k,l),    &
    s2(j,k,l)     )

allocate(          &
    w1(j,k,l,3),   &
    w2(j,k,l,3) )

s1 = 0.
s2 = 0.
w1 = 0.
w2 = 0. 


! user customized
!*****************************************************
!2d source plane
allocate(xx(nm(1),nm(2),nm(3)), &
         yy(nm(1),nm(2),nm(3)), &
         zz(nm(1),nm(2),nm(3)))

! cartesian coordinate of receiver
allocate(rx(nrec), &
         ry(nrec), &
         rz(nrec))
         
! moment rate function read from files
allocate(mrij(nm(1),nm(2),nm(3),6))
allocate(lsummrij(nt,6))
allocate(gsummrij(nt,6))

! depth phases time series window
! travelling time and amplitude with each pair of subfault and receiver
if(nt0>1) allocate(P(nt0,nrec,3),&
	             recP(nt0,nrec),&
	             tp(nm(1),nm(2),nm(3),nsourceoffset,nrec),&
	             ap(nm(1),nm(2),nm(3),nsourceoffset,nrec),&
	             rayp(nm(1),nm(2),nm(3),nsourceoffset,nrec))
	             
if(nt1>1) allocate(PP(nt1,nrec,3),&
			 recPP(nt1,nrec),&
			 tpp(nm(1),nm(2),nm(3),nsourceoffset,nrec),&
	             app(nm(1),nm(2),nm(3),nsourceoffset,nrec),&
	             raypp(nm(1),nm(2),nm(3),nsourceoffset,nrec))
	             
if(nt2>1) allocate(SP(nt2,nrec,3),&
			 recSP(nt2,nrec),&
			 tsp(nm(1),nm(2),nm(3),nsourceoffset,nrec),&
	             asp(nm(1),nm(2),nm(3),nsourceoffset,nrec),&
	             raysp(nm(1),nm(2),nm(3),nsourceoffset,nrec))
	             
if(nt3>1) allocate(S(nt3,nrec,3),SV(nt3,nrec,3),SH(nt3,nrec,3),&
			 recS(nt3,nrec),recSV(nt3,nrec),recSH(nt3,nrec),&
			 ts(nm(1),nm(2),nm(3),nsourceoffset,nrec),&
	             as(nm(1),nm(2),nm(3),nsourceoffset,nrec),&
	             rays(nm(1),nm(2),nm(3),nsourceoffset,nrec))


! assign initial value

xx = 0
yy = 0
zz = 0

x0 = 0
y0 = 0
z0 = 0

rx = 0
ry = 0
rz = 0

mrij = 0

if (nt0>1) then
    P=0
    recP=0
    tp=0
    ap=0
   rayp=0
end if
if (nt1>1) then
    PP=0
    recPP=0
    tpp=0
    app=0
    raypp=0
end if
if (nt2>1) then
    SP=0
    recSP=0
     tsp=0
     asp=0
     raysp=0
end if
if (nt3>1) then
    S=0
    recS=0
    ts=0
    as=0
    rays=0
end if

!any finite number can indicate success
!*****************************************************

end subroutine

end module
