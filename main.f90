! Compute Far-field displacement from moment rate function from SORD
! Parallel IO Library from SORD code
program FarFieldDisplacement

! no need to change below
!------------------------
use m_collective
use m_globals
use m_parameters
use m_setup
use m_arrays
use m_fieldio
use m_util
use m_radiations
use m_assemble
implicit none
!------------------------


! general local variable definition
real :: start_time,stop_time

! user defined variable block **
! user-customized  
! **************************************  
real :: duration
real :: tmps(6),tmpr(6)
integer :: root(3) = 0
!***************************************

! first define thickness of ghost layers
nhalo = 0
call initialize( np0, ip ); master = ip == 0

if ( master ) start_time = MPI_wtime()
if ( master ) write( 0, * ) 
if ( master ) write( 0, * ) 'Depth Phases Generator from SORD'
if ( master ) write( 0, * ) 

call read_parameters
call setup   
call arrays 

if ( master ) write( 0, '(A)' ) 'Input Statistics'
if ( master ) write( 0, '(A)' ) 'Dimension: '
if ( master ) write( 0, '(4(A10,I6))') 'nx=',nn(1),'ny=',nn(2),'nz=',nn(3),'nt=',nt
if ( sync ) call barrier

!***************************************
! -- e.g. read in data by PIOL
if (debug > 1) write(0,*) ip3,nm
! Initial Values

! read in fault coordinate
  call readin_cor
  call receiver

! calculate ray from subfault to receiver

! based on t!s and t!e to switch which phase will be computed
duration = dt * (nt-1)
! P wave statistics
if (nt0 > 1 .and. t0s > 0.) then
  call Ptrvtm_coeff
  if (master) write(0,*)
  call traveltime('P','allmin')
  if(trvmin < t0s .and. master) then
	write(0,*) 'P: t0s < ',trvmin
	stop
  end if
  call traveltime('P','allmax')
  if(trvmax > t0e - duration .and. master) then
  	 write(0,*) 'P: t0e >',trvmax+duration
  	 stop
  end if
  call outtraveltime('P')
end if

! pP wave statistics
if (nt1 > 1 .and. t1s > 0.) then
  call PPtrvtm_coeff
  if (master) write(0,*)
  call traveltime('pP','allmin')
  if(trvmin < t1s .and. master) then
  	write(0,*) 'pP: t1s < ',trvmin
  	stop
  end if
  call traveltime('pP','allmax')
  if(trvmax > t1e - duration .and. master) then
  	write(0,*) 'pP: t1e >',trvmax+duration
  	stop
  end if
  call outtraveltime('pP')
end if

! sP wave statistics
if (nt2 > 1 .and. t2s > 0.) then
  call SPtrvtm_coeff
  if (master) write(0,*)
  call traveltime('sP','allmin')
  if(trvmin < t2s .and. master) then
  	write(0,*) 'sP: t2s < ',trvmin
  	stop
  end if
  call traveltime('sP','allmax')
  if(trvmax > t2e - duration .and. master) then
  	write(0,*) 'sP: t2e >',trvmax+duration
  	stop
  end if
  call outtraveltime('sP')
end if

! S wave statistics
if (nt3 > 1 .and. t3s > 0.) then
  call Strvtm_coeff
  if (master) write(0,*)
  call traveltime('S','allmin')
  if(trvmin < t3s .and. master) then
  	write(0,*) 'S: t2s < ',trvmin
  	stop
  end if
  call traveltime('S','allmax')
  if(trvmax > t3e - duration .and. master) then
  	write(0,*) 'S: t3e >',trvmax+duration
  	stop
  end if
  call outtraveltime('S')
end if

! before this checkpoint, program will judge if array setup can satisfy this scenario

if ( sync ) call barrier
! Parallel read-in data
timeloop: do it = its,ite,ditse

tm = tm0 + (it-1)*dt

! read in moment rate function
mrij = 0.
call fieldio('<','mr11',mrij(:,:,:,1))
call fieldio('<','mr22',mrij(:,:,:,2))
call fieldio('<','mr33',mrij(:,:,:,3))
call fieldio('<','mr23',mrij(:,:,:,4))
call fieldio('<','mr31',mrij(:,:,:,5))
call fieldio('<','mr12',mrij(:,:,:,6))

! convert from big-endian
if (bigendian) call SWAPR44D(mrij)
if (any(isnan(mrij)) .or. maxval(mrij) > huge(0.0d0)) then
	write( 0, * ) 'ERROR: NaN/Inf in READING mrij'
	stop
end if
! compute summed moment rate function
! sum this moment rate tensor into summr
lsummrij(it,:) = 0.
lsummrij(it,:) = sum(sum(sum(mrij,1),1),1)
gsummrij(it,:) = 0.
tmps = lsummrij(it,1:6)
call rreduce1(tmpr,tmps,'allsum',root)
gsummrij(it,1:6) = tmpr

! compute 
if (nt0>1 .and. t0s > 0) call pwave
if (nt1>1 .and. t1s > 0) call ppwave
if (nt2>1 .and. t2s > 0) call spwave
if (nt3>1 .and. t3s > 0) call swave

! output time stamp
if (master) then
	write( 0, '(a)', advance='no' ) '.'
	if ( modulo( it, 50 ) == 0 .or. it == nt ) write( 0, '(i6)' ) it
end if
end do timeloop


! reduce waveforms to master processor
  call assembly
!***************************************
!stop timing
if (master) then
  stop_time = MPI_wtime()
  write(0,*) 'Runing time = ',stop_time - start_time, ' sec'
  write(0,*) 'Finish'
end if
call finalize

end program FarFieldDisplacement

!==================== main program ends ==================================
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~




!==================== subroutine starts ==================================
!------------------------------------------------
!read in coordinate of sources (moment rate)
subroutine readin_cor
use m_globals
use m_collective
use m_fieldio
use m_util
implicit none

call fieldio('<','faultx',xx(:,:,:))
call fieldio('<','faulty',yy(:,:,:))
call fieldio('<','faultz',zz(:,:,:))

if (bigendian) then
	call SWAPR43D(xx(:,:,:))
	call SWAPR43D(yy(:,:,:))
	call SWAPR43D(zz(:,:,:))
end if

! rule out nan and inf after reading
if (any(isnan(xx)) .or. maxval(xx) > huge(0.0d0)) then
	write( 0, * ) 'ERROR: NaN/Inf in READING faultx'
	stop
end if
if (any(isnan(yy)) .or. maxval(yy) > huge(0.0d0)) then
	write( 0, * ) 'ERROR: NaN/Inf in READING faulty'
	stop
end if
if (any(isnan(zz)) .or. maxval(zz) > huge(0.0d0)) then
	write( 0, * ) 'ERROR: NaN/Inf in READING faultz'
	stop
end if

if (usrorig > 0) call neworigin

end subroutine readin_cor

!--------------------------------------------------------------------------------------
! compute average coordinate of faults to get approximate origin (used for receivers)
subroutine neworigin
use m_globals
use m_collective
implicit none

real :: xsum,ysum,zsum
integer :: rootip(3)

xsum = sum(xx)
ysum = sum(yy)
zsum = sum(zz)

rootip = 0

call rreduce0(x0,xsum,'allsum',rootip)
call rreduce0(y0,ysum,'allsum',rootip)
if (usrorig > 1) call rreduce0(z0,zsum,'allsum',rootip)

x0 = x0 / (nn(1)*nn(2)*nn(3))
y0 = y0 / (nn(1)*nn(2)*nn(3))
if (usrorig > 1) z0 = z0 / (nn(1)*nn(2)*nn(3))

if (master) write(0,*) 'New origin is ',x0,y0,z0

end subroutine


!------------------------------------------------
! calculate theta and phi for each station
subroutine receiver
use m_globals
use m_util
implicit none

integer :: i,j,k
real :: iphi, itheta

k = 0
do j = 1, ntheta
  do i = 1, nphi !theta is 0->pi/2; phi = -pi->pi
     k = k + 1
     itheta = otheta + dtheta * (j - 1)
     iphi =  ophi + dphi * (i - 1) ! phi and theta have an offset to avoid nodal plane
     rx(k) = x0 + radius * cos(itheta*deg2rad) * cos(iphi*deg2rad) !sin (degree)
     ry(k) = y0 + radius * cos(itheta*deg2rad) * sin(iphi*deg2rad)
     rz(k) = z0 + radius * sin(itheta*deg2rad)
     ! azimuthal angle is pi/2 - iphi
     ! take-off angle is pi/2 - itheta
  end do
end do
end subroutine receiver

!-----------------
! P wave amplitude coefficient
subroutine Ptrvtm_coeff
use m_globals
use m_util
implicit none

integer :: i,j,k,l,m
real :: x, y, z, r, itheta, iphi

do k = 1,nrec
do m = 1, nsourceoffset
  do l = 1,nm(3)
  do j = 1, nm(2)
    do i = 1, nm(1)

      x = rx(k) - xx(i,j,l) - xyzoffset(m,1)
      y = ry(k) - yy(i,j,l) - xyzoffset(m,2)
      z = rz(k) - zz(i,j,l) - xyzoffset(m,3) !vector from source to receiver

      call cat2sph(x,y,z,r,itheta,iphi)
     
      tp(i,j,l,m,k) = r/vp + xyzoffset(m,4) ! travel time + time offsets
      ap(i,j,l,m,k) = 1./(4*pi*rho*vp**3*r)
      rayp(i,j,l,m,k) = cos(itheta)/vp
      
    end do
  end do
  end do
end do
end do
end subroutine

!-----------------
! S wave amplitude coefficient
subroutine Strvtm_coeff
use m_globals
use m_util
implicit none

integer :: i,j,k,l,m
real :: x, y, z, r, itheta, iphi

do k = 1,nrec
do m = 1, nsourceoffset
  do l = 1, nm(3)
  do j = 1, nm(2)
    do i = 1, nm(1)

      x = rx(k) - xx(i,j,l) - xyzoffset(m,1)
      y = ry(k) - yy(i,j,l) - xyzoffset(m,2)
      z = rz(k) - zz(i,j,l) - xyzoffset(m,3) !vector from source to receiver

      call cat2sph(x,y,z,r,itheta,iphi)
     
      ts(i,j,l,m,k) = r/vs + xyzoffset(m,4)! travel time
      as(i,j,l,m,k) = 1./(4*pi*rho*vs**3*r)
      rays(i,j,l,m,k) = cos(itheta)/vs
      
    end do
  end do
  end do
end do
end do
end subroutine

!-----------------
! pP wave amplitude coefficient
subroutine PPtrvtm_coeff
use m_globals
use m_util
implicit none

integer :: i,j,k,l,m
real :: x, y, z, r, sinp, cosp, itheta, iphi, ii, jj, &
        r1,r2,rtot, pray, rcpp

do k = 1,nrec
do m = 1, nsourceoffset
  do l = 1, nm(3)
  do j = 1, nm(2)
    do i = 1, nm(1)
    
      x = rx(k) - xx(i,j,l) - xyzoffset(m,1)
      y = ry(k) - yy(i,j,l) - xyzoffset(m,2)
      z = rz(k) - zz(i,j,l) - xyzoffset(m,3) !vector from source to receiver
      
      call cat2sph(x,y,z,r,itheta,iphi)
      
      sinp = sin(iphi)
      cosp = cos(iphi)
      
      ii = abs(atan(sqrt(x**2+y**2)/abs(rz(k)+zz(i,j,l)))) ! P incidence
      jj = asin(vs/vp*sin(ii)) !equivalent S incidence
      
      r1=abs(zz(i,j,l)/cos(ii))
	r2=abs(rz(k)/cos(ii))
	rtot = r1+r2
	
	if (rtot < r) rtot=r
	if (rtot > 2*r) rtot = 2*r
	
	pray = sin(ii)/vp
	
	rcpp = (4*pray**2*cos(ii)*cos(jj)/vp/vs-(1/vs**2-2*pray**2)**2)/ &
		 (4*pray**2*cos(ii)*cos(jj)/vp/vs+(1/vs**2-2*pray**2)**2)
	
	tpp(i,j,l,m,k) = rtot/vp + xyzoffset(m,4)! travel time
      app(i,j,l,m,k) = rcpp/(4*pi*rho*vp**3*rtot)
      raypp(i,j,l,m,k) = pray
      
    end do
  end do
  end do
end do
end do
end subroutine

!-----------------
! sP wave amplitude coefficient
subroutine SPtrvtm_coeff
use m_globals
use m_collective
use m_util
implicit none

integer :: i,j,l,k,m
real :: x, y, z, r, sinp, cosp, itheta, iphi, ii, jj, &
        r1,r2,rtot, pray, rcsp, ratio, c1, c2, c3, c4

do k = 1,nrec
do m = 1, nsourceoffset
  do l = 1, nm(3)
  do j = 1, nm(2)
    do i = 1, nm(1)
    
      x = rx(k) - xx(i,j,l) - xyzoffset(m,1)
      y = ry(k) - yy(i,j,l) - xyzoffset(m,2)
      z = rz(k) - zz(i,j,l) - xyzoffset(m,3) !vector from source to receiver
      
      call cat2sph(x,y,z,r,itheta,iphi)
      
      sinp = sin(iphi)
      cosp = cos(iphi)
      
      !parameters used to estimate misfit function
      c1 = sqrt(x**2+y**2)
      c2 = rz(k)
      c3 = zz(i,j,l)
      c4 = vs/vp
      
	call bisection(c1,c2,c3,c4,ii)
	
      jj = asin(vs/vp*sin(ii)) !% up going S wave
      
      r1=abs(zz(i,j,l)/cos(jj))      ! S path
      r2=abs(rz(k)/cos(ii))        ! P path   
      rtot=r1+r2
      
      if (rtot < r) rtot=r
	if (rtot > 2*r) rtot = 2*r
	
	pray = sin(ii)/vp

	
	rcsp = 4/vp*pray*cos(jj)*(1/vs**2-2*pray**2)/&
             ((1/vs**2-2*pray**2)**2+4*pray**2*cos(ii)*cos(jj)/vp/vs)
             
      RATIO=((1-vp**2*pray**2)/(1-vs**2*pray**2))**(1/4)
      
	tsp(i,j,l,m,k) = r1/vs + r2/vp + xyzoffset(m,4)! travel time
      asp(i,j,l,m,k) = rcsp/(4*pi*rho*vs**2*vp*rtot)*ratio
      raysp(i,j,l,m,k) = pray
    end do
  end do
  end do
end do
end do
end subroutine


!-----------------------------------------------
! estimate travel time
subroutine traveltime(phase,op)
use m_globals
use m_collective
implicit none

character(*) :: phase,op !op can be 'sum,'min', 'max'
real :: ttime, avetime
integer :: root(3) = 0


select case(phase)
case('P')
	select case(op)
	case('ave');	ttime = sum(tp(:,:,:,:,:))
	case('allmin');	ttime = minval(tp(:,:,:,:,:))
	case('allmax');	ttime = maxval(tp(:,:,:,:,:))
	end select
	call rreduce0(avetime,ttime,op,root)
	if (op=='ave') avetime = avetime/(nn(1)*nn(2)*nn(3)*nrec*nsourceoffset)
case('S')
	select case(op)
	case('ave');	ttime = sum(ts(:,:,:,:,:))
	case('allmin');	ttime = minval(ts(:,:,:,:,:))
	case('allmax');	ttime = maxval(ts(:,:,:,:,:))
	end select
	call rreduce0(avetime,ttime,op,root)
	if (op=='ave') avetime = avetime/(nn(1)*nn(2)*nn(3)*nrec*nsourceoffset)	
case('pP')
	select case(op)
	case('ave');	ttime = sum(tpp(:,:,:,:,:))
	case('allmin');	ttime = minval(tpp(:,:,:,:,:))
	case('allmax');	ttime = maxval(tpp(:,:,:,:,:))
	end select
	call rreduce0(avetime,ttime,op,root)
	if (op=='ave') avetime = avetime/(nn(1)*nn(2)*nn(3)*nrec*nsourceoffset)
case('sP')
	select case(op)
	case('ave');	ttime = sum(tsp(:,:,:,:,:))
	case('allmin');	ttime = minval(tsp(:,:,:,:,:))
	case('allmax');	ttime = maxval(tsp(:,:,:,:,:))
	end select
	call rreduce0(avetime,ttime,op,root)
	if (op=='ave') avetime = avetime/(nn(1)*nn(2)*nn(3)*nrec*nsourceoffset)
end select

	if (master) write(0,*) phase, &
	           ' wave ',op,' is ',avetime, ' sec'
	if(op == 'allmin') trvmin = avetime
	if(op == 'allmax') trvmax = avetime

end subroutine

! ------------------------------------------------
! write out minimum travle time of each receiver
subroutine outtraveltime(phase)
use m_globals
use m_collective
use m_fieldio
implicit none

character(*) :: phase
real, allocatable,dimension(:) :: stime, rtime
integer :: root(3) = 0
integer :: i

allocate(stime(nrec),rtime(nrec))
select case(phase)
	case('P') 
		do i = 1, nrec 
			stime(i) = minval(tp(:,:,:,:,i))
		end do
	case('pP') 
		do i = 1, nrec 
			stime(i) = minval(tpp(:,:,:,:,i))
		end do
	case('sP') 
		do i = 1, nrec 
			stime(i) = minval(tsp(:,:,:,:,i))
		end do
	case('S') 
		do i = 1, nrec 
			stime(i) = minval(ts(:,:,:,:,i))
		end do
end select
		
call rreduce1(rtime,stime,'min',root)

if(master) call mstrwrite1d(rtime,trim(oufile)//'/'//trim(phase)//'-mintraveltime')
end subroutine

