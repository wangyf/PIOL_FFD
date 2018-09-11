! Read model parameters
module m_parameters
implicit none
contains

subroutine read_parameters
use m_globals
use m_fieldio
use m_strings
integer :: ios, i, nargs, j
character(32) :: key
character(1) :: op,delims=' '
character(2560) :: line,args(3)

! I/O pointers
allocate( io0 )
io => io0
io%next => io0
io%field = 'head'

open( 1, file='parameters', status='old' )

doline: do

! Read line
read( 1, '(a)', iostat=ios ) line
if ( ios /= 0 ) exit doline

! Strip comments and punctuation
str = line
i = scan( str, '#' )
if ( i > 0 ) str(i:) = ' '
do
    i = scan( str, "()[]{}'," )
    if ( i == 0 ) exit
    str(i:i) = ' '
end do

! Read key val pair
if ( str == '' ) cycle doline
read( str, *, iostat=ios ) key
write(0,*) str
! Select input key
select case( key )
case( 'fieldio', '' )
!------------------- no need to change ------------------------
case( 'nn' );           read( str, *, iostat=ios ) key, op, nn
case( 'nt' );           read( str, *, iostat=ios ) key, op, nt
case( 'dx' );           read( str, *, iostat=ios ) key, op, dx
case( 'dt' );           read( str, *, iostat=ios ) key, op, dt
case( 'tm0' );          read( str, *, iostat=ios ) key, op, tm0
case( 'np3' );          read( str, *, iostat=ios ) key, op, np3
case( 'debug' );          read( str, *, iostat=ios ) key, op, debug
case( 'itio' );         read( str, *, iostat=ios ) key, op, itio
case( 'itstats' );      read( str, *, iostat=ios ) key, op, itstats
! # input:  0=separate files, 1=MPI-IO, -1=non-collective MPI-IO
! # output: 0=separate files, 1=MPI-IO, -1=non-collective MPI-IO
case( 'mpin' );         read( str, *, iostat=ios ) key, op, mpin
case( 'mpout' );        read( str, *, iostat=ios ) key, op, mpout
!-------------------------------------------------------------

! user customized parameter input in the following
!*************************************************************
case( 'rho' );    read( str, *, iostat=ios ) key,op, rho
case( 'vp' );    read( str, *, iostat=ios ) key,op, vp
case( 'vs' );    read( str, *, iostat=ios ) key,op, vs
case( 'dphi' )   
    read( str, *, iostat=ios ) key,op, dphi
    ophi = -180 + dphi/2!dphi (default)
case( 'ophi' );    read( str, *, iostat=ios ) key,op, ophi !offset phi
case( 'dtheta')
    read( str, *, iostat=ios ) key,op, dtheta
    otheta = 0 + dtheta/2 ! default dtheta 
case( 'otheta');   read( str, *, iostat=ios ) key,op, otheta !offset theta
case( 'radius');   read( str, *, iostat=ios ) key,op, radius !radius for focal sphere
case( 't0s' );    read( str, *, iostat=ios ) key,op, t0s  !time window starts for P
case( 't0e' );    read( str, *, iostat=ios ) key,op, t0e  !time window ends for P
case( 't1s' );    read( str, *, iostat=ios ) key,op, t1s  !pP
case( 't1e' );    read( str, *, iostat=ios ) key,op, t1e  !pP
case( 't2s' );    read( str, *, iostat=ios ) key,op, t2s  !sP
case( 't2e' );    read( str, *, iostat=ios ) key,op, t2e  !sP
case( 't3s' );    read( str, *, iostat=ios ) key,op, t3s  !S
case( 't3e' );    read( str, *, iostat=ios ) key,op, t3e  !S
case( 'comp' );    read( str, *, iostat=ios ) key,op, comp  !0 is amplitude 1 is all components
case( 'usrorig' );    read( str, *, iostat=ios ) key,op, usrorig  !use average fault as origin (not default)
case( 'endian' );    read( str, *, iostat=ios ) key,op, endian !1 is big-endian input 0 is little-endian
case( 'nsourceoffset' )
    read( str, *, iostat=ios ) key,op, nsourceoffset
    allocate(sourceoffset(4*nsourceoffset)) 
    allocate(xyzoffset(nsourceoffset,4))
case( 'sourceoffset')
	read( str, *, iostat=ios ) key,op, sourceoffset
	do j = 1, nsourceoffset
		xyzoffset(j,:) = sourceoffset((j-1)*4+1:(j-1)*4+4)
	end do
case( 'infile' ) !folder has '/'
   call parse(str, delims, args, nargs) 
   infile = args(3)
case( 'oufile' )
   call parse(str, delims, args, nargs) 
   oufile = args(3)
!*************************************************************
case default
    select case( key(1:1) )
    case( '=', '+' )
        call pappend
        io%ib = -1
        !XXXread( str, *, iostat=ios ) io%mode, io%nc
        if (debug > 1 .and. master) write(0,*) str
        read( str, *, iostat=ios ) io%mode, io%nc, io%tfunc, &
            io%period, io%x1, io%x2, io%nb, io%ii, io%filename, &
            io%val, io%field
        if(debug>1 .and. master ) write(0,*) io%mode,io%filename,io%val,io%field
    case default; ios = 1
    end select
end select

! Error check
if ( ios /= 0 ) then
    if ( master ) write( 0, * ) 'bad input: ', trim( line )
    stop
end if

end do doline

close( 1 )

!look for it1 and it2
! this is only used for uniformly time-dependent input
io => io0
loop: do while( io%next%field /= 'head' )
ioprev => io
io => io%next

if((io%mode(2:2)=='r' .or. io%mode(2:2)=='R') .and. &
   (io%ii(1,4)/=0 .and. io%ii(2,4)/=0 .and. io%ii(2,4)>=io%ii(1,4))) then
its = io%ii(1,4)
ite = io%ii(2,4)
ditse = io%ii(3,4)
	if (master .and. its < 1 .or. ite > nt) then
		write(0,*) 'reset nt (nt < ite)'
		stop
	end if
end if
end do loop


!output check
verb = master .and. debug > 0
if (master) then
write(0,*) 'nn = ',nn
write(0,*) 'np3 = ',np3
write(0,*) 'Src nt = ',nt
write(0,*) 'Sample it range: ',its,ite,ditse
write(0,*) 'dt = ',dt
write(0,*) 'P wave starts at ',t0s, ' s'
write(0,*) 'P wave ends at ',t0e, ' s'
end if

end subroutine

end module