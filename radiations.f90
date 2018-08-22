! compute seismic waves and spectra
module m_radiations
implicit none
contains

! P wave radiation pattern term in amplitude
subroutine Prdpt(mrij,x,y,z,rp)
use m_util
real,intent(in) :: mrij(6)
real,intent(in) :: x, y, z
real,intent(out) :: rp(3)
real :: r,itheta,iphi, &
        st,ct,sp,cp,e_r(3),e_t(3),e_p(3)

call cat2sph(x,y,z,r,itheta,iphi)

st = sin(itheta)
ct = cos(itheta)
sp = sin(iphi)
cp = cos(iphi)

e_r = (/ ct*cp, ct*sp, st/)  !P wave vector Aki and Richard P78
!e_t = (/ st*cp, st*sp, -ct/) !SV wave
!e_p = (/ -sp, cp, 0/) !SH wave

rp = e_r * vmv(mrij,e_r,e_r) !Aki and Richard P111

end subroutine


! S wave radiation pattern term in amplitude
subroutine Srdpt(mrij,x,y,z,rpv,rph)
use m_util
real,intent(in) :: mrij(6)
real,intent(in) :: x, y, z
real,intent(out) :: rpv(3),rph(3)
real :: r,itheta,iphi, &
        st,ct,sp,cp,e_r(3),e_t(3),e_p(3)

call cat2sph(x,y,z,r,itheta,iphi)

st = sin(itheta)
ct = cos(itheta)
sp = sin(iphi)
cp = cos(iphi)

e_r = (/ ct*cp, ct*sp, st/)  !P wave vector
e_t = (/ st*cp, st*sp, -ct/) !SV wave
e_p = (/ -sp, cp, 0./) !SH wave

rpv = e_t * vmv(mrij,e_t,e_r) !SV wave

rph = e_p * vmv(mrij,e_p,e_r) !SV wave

end subroutine


! pP wave radiation pattern term in amplitude
subroutine PPrdpt(mrij,sx,sy,sz,rx,ry,rz,rp)
use m_util
real,intent(in) :: mrij(6)
real,intent(in) :: sx,sy,sz,rx,ry,rz
real,intent(out) :: rp(3)
real :: e_r1(3),e_r(3) !ray vector of pP and P
real :: x, y, z, r, itheta, iphi, sp, cp, ii

x = rx - sx
y = ry - sy
z = rz - sz

call cat2sph(x,y,z,r,itheta,iphi)
sp=sin(iphi)
cp=cos(iphi)

ii = abs(atan(sqrt(x**2+y**2)/abs(sz+rz)))

e_r1 = (/sin(ii)*cp, sin(ii)*sp, cos(ii)/)
e_r  = (/sin(ii)*cp, sin(ii)*sp, -cos(ii)/)

rp = e_r * vmv(mrij,e_r1,e_r1)

end subroutine


! sP wave radiation pattern term in amplitude
subroutine SPrdpt(mrij,sx,sy,sz,rx,ry,rz,rp)
use m_util
real,intent(in) :: mrij(6)
real,intent(in) :: sx,sy,sz,rx,ry,rz
real,intent(out) :: rp(3)
real :: e_r1(3),e_t1(3),e_r(3)
real :: x, y, z, r, itheta, iphi, sp, cp, ii, jj, iflip

x = rx - sx
y = ry - sy
z = rz - sz

call cat2sph(x,y,z,r,itheta,iphi)
sp=sin(iphi)
cp=cos(iphi)

ii = abs(atan(sqrt(x**2+y**2)/abs(sz+rz)))

e_r1 = (/sin(jj)*cp, sin(jj)*sp, cos(jj)/) !polarization outgoing
e_t1 = (/cos(jj)*cp, cos(jj)*sp, -sin(jj)/) !outgoing 
e_r =  (/sin(ii)*cp, sin(ii)*sp, -cos(ii)/) !P incident to receiver


iflip = sum(e_t1(1:2)*e_r(1:2),1)
if (iflip > 0) then
  iflip = 1
elseif(iflip < 0) then
  iflip = -1
else
  iflip = 0
end if
rp = e_r * iflip * vmv(mrij,e_t1,e_r1)
end subroutine



! compute P wave seismogram
subroutine pwave
use m_globals
integer :: i,j,k,l, m, nit1,nit2
real :: x, y, z, mr(6), rp(3), time,tmpam(3),wei1,wei2

!$omp parallel do schedule(static) private(i,j,k,rp,x,y,z,tmpam,time,nit)
do k = 1,nrec
do m = 1, nsourceoffset
  do l = 1,nm(3)
  do j = 1, nm(2)
    do i = 1, nm(1)
  !do j = 35,35
    !do i = 1,1
    
      x = rx(k) - xx(i,j,l) - xyzoffset(m,1)
      y = ry(k) - yy(i,j,l) - xyzoffset(m,2)
      z = rz(k) - zz(i,j,l) - xyzoffset(m,3) !vector from source to receiver
      
      call Prdpt(mrij(i,j,l,1:6),x,y,z,rp)
      !if (ip == 0) print *,ap(i,j,k),rp
      !ap(i,j,l,k) = 1
      tmpam = ap(i,j,l,m,k)*rp
      !tmpam = mrij(i,j,l,5)
      !if (ip==1) print *,rp
      time = tm + tp(i,j,l,m,k) - t0s ! time after t0s
      !if (master) write(0,*) 'P ',time
      nit1 = floor(time/dt)+1
      nit2 = ceiling(time/dt)+1
      wei1 = 1-(time - (nit1-1)*dt)/dt
      wei2 = 1 - wei1
      if (nit2 > nt0 .or. nit1 < 1) then
          write(0,*) 't0s and t0e are not long enough for P'
          stop
      end if
      !$omp critical
      P(nit1,k,1:3) = P(nit1,k,1:3) + tmpam(1:3) * wei1
      P(nit2,k,1:3) = P(nit2,k,1:3) + tmpam(1:3) * wei2
      !$omp end critical
    end do
  end do
  end do
end do  
end do
!$omp end parallel do
end subroutine


! compute S wave seismogram
subroutine swave
use m_globals
integer :: i,j,l,k, m, nit1,nit2
real :: x, y, z, mr(6), rpv(3), rph(3), time,tmpam1(3),tmpam2(3),wei1,wei2

!$omp parallel do schedule(static) private(i,j,k,rp,x,y,z,tmpam,time,nit)
do k = 1,nrec
do m = 1, nsourceoffset
  do l = 1,nm(3)
  do j = 1, nm(2)
    do i = 1, nm(1)
    
      x = rx(k) - xx(i,j,l) - xyzoffset(m,1)
      y = ry(k) - yy(i,j,l) - xyzoffset(m,2)
      z = rz(k) - zz(i,j,l) - xyzoffset(m,3) !vector from source to receiver
      
      call Srdpt(mrij(i,j,l,1:6),x,y,z,rpv,rph)
      !if (ip == 0) print *,ap(i,j,k),rp
      tmpam1 = as(i,j,l,m,k)*rpv
      tmpam2 = as(i,j,l,m,k)*rph
      !if (ip==1) print *,rp
      time = tm + ts(i,j,l,m,k) - t3s ! time after t0s
      !if (master) write(0,*) 'P ',time
      nit1 = floor(time/dt)+1
      nit2 = ceiling(time/dt)+1
      wei1 = 1-(time - (nit1-1)*dt)/dt
      wei2 = 1 - wei1
      if (nit2 > nt3 .or. nit1 < 1) then
          write(0,*) 't3s and t3e are not long enough for S'
          stop
      end if
      !$omp critical
      SV(nit1,k,1:3) = SV(nit1,k,1:3) + tmpam1(1:3) * wei1
      SV(nit2,k,1:3) = SV(nit2,k,1:3) + tmpam1(1:3) * wei2
      
      SH(nit1,k,1:3) = SH(nit1,k,1:3) + tmpam2(1:3) * wei1
      SH(nit2,k,1:3) = SH(nit2,k,1:3) + tmpam2(1:3) * wei2
      
      S(nit1,k,1:3) = SV(nit1,k,1:3)+SH(nit1,k,1:3)
      S(nit2,k,1:3) = SV(nit2,k,1:3)+SH(nit2,k,1:3)
      !$omp end critical
    end do
  end do
  end do
end do  
end do
!$omp end parallel do
end subroutine




! compute pP wave seismogram
subroutine ppwave
use m_globals
integer :: i,j,l,k,m, nit1,nit2
real :: x, y, z, mr(6), rp(3), time, tmpam(3),wei1,wei2

!$omp parallel do schedule(static) private(i,j,k,rp,x,y,z,tmpam,time,nit)
do k = 1,nrec
do m = 1, nsourceoffset
  do l = 1, nm(3)
  do j = 1, nm(2)
    do i = 1, nm(1)
      
      call PPrdpt(mrij(i,j,l,1:6),xx(i,j,l)+xyzoffset(m,1),yy(i,j,l)+ &
                  xyzoffset(m,2),zz(i,j,l)+xyzoffset(m,3),rx(k), &
                 ry(k),rz(k),rp)
      tmpam = app(i,j,l,m,k)*rp
      
      time = tm + tpp(i,j,l,m,k) - t1s ! time after t1s

      nit1 = floor(time/dt)+1
      nit2 = ceiling(time/dt)+1
      wei1 = 1-(time - (nit1-1)*dt)/dt
      wei2 = 1 - wei1
      if (nit2 > nt1 .or. nit1 < 1) then
          write(0,*) 't1s and t1e are not long enough for pP'
          stop
      end if
      !$omp critical
      PP(nit1,k,1:3) = PP(nit1,k,1:3) + tmpam(1:3) * wei1
      PP(nit2,k,1:3) = PP(nit2,k,1:3) + tmpam(1:3) * wei2
      !$omp end critical
    end do
  end do
  end do
end do
end do
!$omp end parallel do
end subroutine

! compute sP wave seismogram
subroutine spwave
use m_globals
integer :: i,j,l,k, m, nit1,nit2
real :: x, y, z, mr(6), rp(3), time, tmpam(3),wei1,wei2

!$omp parallel do schedule(static) private(i,j,k,rp,x,y,z,tmpam,time,nit)
do k = 1,nrec
do m = 1, nsourceoffset
  do l = 1,nm(3)
  do j = 1, nm(2)
    do i = 1, nm(1)
      
      call SPrdpt(mrij(i,j,l,1:6),xx(i,j,l)+xyzoffset(m,1),yy(i,j,l)+ &
                  xyzoffset(m,2),zz(i,j,l)+xyzoffset(m,3),rx(k), &
                 ry(k),rz(k),rp)
      tmpam = asp(i,j,k,m,l)*rp
      
      time = tm + tsp(i,j,l,m,k) - t2s ! time after t2s
      
      nit1 = floor(time/dt)+1
      nit2 = ceiling(time/dt)+1
      wei1 = 1-(time - (nit1-1)*dt)/dt
      wei2 = 1 - wei1

      if (nit2 > nt2 .or. nit1 < 1) then
          write(0,*) 't2s and t2e are not long enough for sP'
          stop
      end if
      !$omp critical
      SP(nit1,k,1:3) = SP(nit1,k,1:3) + tmpam(1:3) * wei1
      SP(nit2,k,1:3) = SP(nit2,k,1:3) + tmpam(1:3) * wei2
      !$omp end critical
    end do
  end do
  end do
end do
end do  
!$omp end parallel do
end subroutine

real function vmv(mij,e1,e2)
real,intent(in) :: mij(6),e1(3),e2(3)

vmv = e1(1)*e2(1)*mij(1) + e1(2)*e2(1)*mij(6) + e1(3)*e2(1)*mij(5) + &
      e1(1)*e2(2)*mij(6) + e1(2)*e2(2)*mij(2) + e1(3)*e2(2)*mij(4) + &
      e1(1)*e2(3)*mij(5) + e1(2)*e2(3)*mij(4) + e1(3)*e2(3)*mij(3)

end function
end module