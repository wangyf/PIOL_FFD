! Parallel IO Library from SORD code
program PIOL
use m_collective
use m_globals
use m_parameters
use m_setup
use m_arrays
use m_fieldio
use m_util

implicit none


interface 
  subroutine lsq(nt,dt,tms,t0,tp,su,tr,tmsfit)
    use m_util
    integer,intent(in) :: nt
    real,intent(in) :: dt,t0,tp,su
    real,intent(in) :: tms(:)
    real,intent(out) :: tr,tmsfit(:)
  end subroutine
end interface



! general variable defination
integer :: i,j,k,l

! user defined variable block
! user-customized   
real, dimension(:), allocatable :: time 
real, dimension(:,:,:),allocatable :: f1,f2,f3
real, dimension(:,:,:),allocatable :: psv,t0,tp,tr,te,su,ct !tr = te - t0
real, dimension(:,:,:,:), allocatable :: tms,tmsint,tmsfit !tms is slip velocity tmsint is slip
integer,dimension(:,:,:), allocatable :: t0_ind,tp_ind,tr_ind,te_ind
real :: sur_t0=0.01,sur_te=0.95 ! slip ratio corresponding to t0 and te
!!!!

! first define thickness of ghost layers
nhalo = 0

call initialize( np0, ip ); master = ip == 0
call read_parameters
call setup   
call arrays                          
!if ( master ) write( 0, * ) 'User defined a name here supported by PIOL'
if ( master ) write( 0, '(A)' ) 'Rupture Statistics'
if ( master ) write( 0, '(A)' ) 'Dimension: '
if ( master ) write( 0, '(4(A10,I6))') 'nx=',nn(1),'ny=',nn(2),'nz=',nn(3),'nt=',nt
if ( sync ) call barrier

! -- e.g. read in data by PIOL
if (debug > 1) write(0,*) ip3,nm
allocate(f1(nm(1),nm(2),nm(3)))
allocate(f2(nm(1),nm(2),nm(3)))
allocate(f3(nm(1),nm(2),nm(3)))
allocate(tms(nm(1),nm(2),nm(3),nt))
allocate(tmsfit(nm(1),nm(2),nm(3),nt))
allocate(tmsint(nm(1),nm(2),nm(3),nt))
allocate(time(nt))
allocate(psv(nm(1),nm(2),nm(3)))
allocate(t0(nm(1),nm(2),nm(3)))
allocate(tp(nm(1),nm(2),nm(3)))
allocate(tr(nm(1),nm(2),nm(3)))
allocate(te(nm(1),nm(2),nm(3)))
allocate(su(nm(1),nm(2),nm(3)))
allocate(ct(nm(1),nm(2),nm(3)))
allocate(t0_ind(nm(1),nm(2),nm(3)))
allocate(tp_ind(nm(1),nm(2),nm(3)))
allocate(tr_ind(nm(1),nm(2),nm(3)))
allocate(te_ind(nm(1),nm(2),nm(3)))

! Initial Values
f1 = 0.
f2 = 0.
f3 = 0.
tms = 0.
tmsfit = 0.
tmsint = 0.
time = 0.
psv = 0.
t0 = 0.
tp = 0.
tr = 0.
te = 0.
su = 0.
ct = 0.

t0_ind = 0
tp_ind = 0
tr_ind = 0
te_ind = 0


if ( sync ) call barrier
! Parallel read-in data
timeloop: do it = 1,nt

if (debug > 2) write(0,*) 'it=',it,io0%next%field
call fieldio('<','sv1',f1)
call fieldio('<','sv2',f2)
call fieldio('<','sv3',f3)

if( master .and. (mod(it,50)==0 .or. it == nt)) write(0,*) '..........',it

! store time series of each fault element in each node
time(it) = (it-1)*dt
tms(:,:,:,it) = sqrt(f1 * f1 + f2 * f2 + f3 * f3)
if (it > 1) tmsint(:,:,:,it) = tmsint(:,:,:,it-1) + dt * (tms(:,:,:,it-1) + tms(:,:,:,it) ) / 2.
end do timeloop
if (debug > 2) write(0,*) io0%next%field

! final slip
su(:,:,:) = tmsint(:,:,:,nt)
!if(master) write(0,*) nm

l = 0
! measure peak slip velocity
do k = 1, nm(3)
  do j = 1, nm(2)
    do i = 1, nm(1)
    
        l = l + 1
        if (master .and. mod(l,10)==0) write(0,'(F8.2 A)') 100.*l/(nm(1)*nm(2)*nm(3)), '% complete'
        ! peak slip velocity
        psv(i,j,k) = maxval(tms(i,j,k,:))
        
        if ( psv(i,j,k) > 1e-3 ) then
          
          ! peak time
          tp_ind(i,j,k) = maxloc(tms(i,j,k,:),1)
          tp(i,j,k) = time(tp_ind(i,j,k))
          
          ! slip ratio
          tmsint(i,j,k,:) = tmsint(i,j,k,:) / su(i,j,k)
        
          t0_ind(i,j,k) = find_first(tmsint(i,j,k,:),nt,'>=',sur_t0)
          te_ind(i,j,k) = find_last(tmsint(i,j,k,:),nt,'<=',sur_te)
          
          tp_ind(i,j,k) = tp_ind(i,j,k) - t0_ind(i,j,k)
          tp(i,j,k) = tp_ind(i,j,k) * dt
          
          if (t0_ind(i,j,k) > te_ind(i,j,k) ) then
            write(0,*) 'Error in picking t0 and te (t0 is larger than te here)'
            stop
          end if
          
          ! initial and ending time
          t0(i,j,k) = time(t0_ind(i,j,k))
          te(i,j,k) = time(te_ind(i,j,k))
          
          !if (debug > 0 .and. master) &
            !write(0,'(3I5,4F8.4)') i,j,k,t0(i,j,k),tp(i,j,k),te(i,j,k),psv(i,j,k)
          
          ! rise time
          !tr(i,j,k) = te(i,j,k) - t0(i,j,k)
          !tr_ind(i,j,k) = te_ind(i,j,k) - t0_ind(i,j,k)
          
          ! Yoffe function fitting
          !call lsq(nt,dt,tms(i,j,k,:),t0(i,j,k),tp(i,j,k),su(i,j,k),tr(i,j,k),tmsfit(i,j,k,:))

           
          ! New definition of rise time from moment
          ! centroid time of slip rate function
          ct(i,j,k) = trapz(tms(i,j,k,:)*time(:),nt,dt)/su(i,j,k)
          tr(i,j,k) = sqrt(trapz((time(:) - ct(i,j,k))*(time(:) - ct(i,j,k))*tms(i,j,k,:),nt,dt)/su(i,j,k))
        end if
        
    end do 
  end do
end do
call barrier
if ( master) write(0,*) 'Statistics process complete'

! static field when it == nt
it = nt
! output field
call fieldio('>','psv',psv)
!if(master) write(0,*) '1'
call fieldio('>','t0',t0)
!if(master) write(0,*) '2'
call fieldio('>','tp',tp)
!if(master) write(0,*) '3'
call fieldio('>','te',te)
!if(master) write(0,*) '4'
call fieldio('>','tr',tr)
!if(master) write(0,*) '5'
call fieldio('>','su',su)
!if(master) write(0,*) '6'
if ( master ) write(0,*) 'Static field write complete'

call barrier

do it = 1,nt
  if( master .and. (mod(it,50)==0 .or. it == nt)) write(0,*) '..........',it
  call fieldio('>','svm', tms(:,:,:,it))
!if (nancheck3d(tmsfit(:,:,:,it))==1) then
  call fieldio('>','tinti',tmsfit(:,:,:,it))
!else
!  if(master) write(0,*) 'Nan Found, Stop'
!  stop
!end if
end do

if ( master ) write(0,*) 'Time-variable field write complete'

if ( sync ) call barrier
call finalize
if (master) write(0,*) '--Done--'
 
deallocate(f1)
deallocate(f2)
deallocate(f3)
deallocate(tms)
deallocate(tmsfit)
deallocate(tmsint)
deallocate(time)
deallocate(psv)
deallocate(t0)
deallocate(tp)
deallocate(tr)
deallocate(te)
deallocate(su)
deallocate(t0_ind)
deallocate(tp_ind)
deallocate(tr_ind)
deallocate(te_ind)

end program

subroutine lsq(nt,dt,tms,t0,tp,su,tr,tmsfit)
use m_util
implicit none
integer,intent(in) :: nt
real,intent(in) :: dt,t0,tp,su
real,intent(in) :: tms(:)
real,intent(out) :: tr,tmsfit(:)

!local
real,dimension(nt) :: t,tmpt
integer :: i,np
real :: p,p1,p2,dp,fmin,misfit
!

do i = 1, nt
  t(i) = (i - 1 ) * dt
end do

! range for tr searching
p1 = 0
!p2 = t(nt)/2 ! or user define
p2 = 10
dp = 0.25
np = NINT((p2-p1)/dp+1)

fmin = 1e9
do i = 1, np
  p = p1 + dp * (i - 1)
  
  call tinti(nt,dt,t0,tp,p,su,tmpt)
  misfit = sqrt(sum((tmpt-tms)*(tmpt-tms),1)/nt)  

  if (misfit < fmin) then
    fmin = misfit
    tr = p
    tmsfit = tmpt
  end if
end do


end subroutine
