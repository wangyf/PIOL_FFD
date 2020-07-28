! Setup model dimensions
module m_setup
implicit none
contains

subroutine setup
use m_globals
use m_collective
use m_util
integer :: nl(3),nmore(3),i

! user customized
!*************************
nt0 = nint((t0e-t0s)/dt)+1
nt1 = nint((t1e-t1s)/dt)+1
nt2 = nint((t2e-t2s)/dt)+1
nt3 = nint((t3e-t3s)/dt)+1 ! additional S wave not within P wave package

if (lattice == 'lon-lat') then
    ntheta = nint(90./dtheta)
    nphi = nint(360./dphi)
    nrec = nphi*ntheta
else if (lattice == 'Fibonacci') then
    nrec = nfibonacci
end if

if(master) write(0,*) 'nrec= ',nrec
if (nt0>1 .and. t0s > 0 .and. master) write(0,*) 'P wave nt0 = ',nt0
if (nt1>1 .and. t1s > 0 .and. master) write(0,*) 'pP wave nt1 = ',nt1
if (nt2>1 .and. t2s > 0 .and. master) write(0,*) 'sP wave nt2 = ',nt2
if (nt3>1 .and. t3s > 0 .and. master) write(0,*) 'S wave nt3 = ',nt3

if (master) then
    if (nt0>1 .and. t0s > 0) then
        write(0,*) 'P wave starts at ',t0s, ' s'
        write(0,*) 'P wave ends at ',t0e, ' s'
    end if
    if (nt1>1 .and. t1s > 0) then
        write(0,*) 'pP wave starts at ',t1s, ' s'
        write(0,*) 'pP wave ends at ',t1e, ' s'
    end if
    if (nt2>1 .and. t2s > 0) then
        write(0,*) 'sP wave starts at ',t2s, ' s'
        write(0,*) 'sP wave ends at ',t2e, ' s'
    end if
    if (nt3>1 .and. t3s > 0) then
        write(0,*) 'S wave starts at ',t3s, ' s'
        write(0,*) 'S wave ends at ',t3e, ' s'
    end if
end if

bigendian = endian > 0
if (master .and. bigendian) write(0,'(A)') 'Endianness of Input data is different with this machine'
!*************************




! Partition for parallelization
if ( np0 == 1 ) np3 = 1
nl3 = (nn - 1) / np3 + 1

nl3 = max( nl3, nhalo )

call rank( ip3, ipid, np3 )
nl = nl3
nnoff = nl3 * ip3 - nhalo

do i = 1, 3
        nmore(i) = nn(i) - np3(i) * (nl3(i) - 1)
        if(ip3(i) .gt. nmore(i) - 1) then
                nl(i) = nl3(i)-1
                nnoff(i) = nnoff(i) - (ip3(i) - nmore(i))
        end if
end do

! Size of arrays
nm = nl + 2 * nhalo

i1core = 1  + nhalo
i2core = nm - nhalo

! Debugging
verb = master .and. debug > 1
sync = debug > 0
if ( debug > 2 ) then
    write( str, "( a,i6.6,a )" ) 'debug/db', ipid, '.py'
    open( 1, file=str, status='replace' )
    write( 1, "( 'ip      = ', i8                                            )" ) ip
    write( 1, "( 'ipid    = ', i8                                            )" ) ipid
    write( 1, "( 'np3     = ', i8, 2(',', i8)                                )" ) np3
    write( 1, "( 'ip3     = ', i8, 2(',', i8)                                )" ) ip3
    write( 1, "( 'nn      = ', i8, 2(',', i8)                                )" ) nn
    write( 1, "( 'nm      = ', i8, 2(',', i8)                                )" ) nm
    write( 1, "( 'nhalo   = ', i8, 2(',', i8)                                )" ) nhalo
    write( 1, "( 'nnoff   = ', i8, 2(',', i8)                                )" ) nnoff
    close( 1 )
end if

end subroutine

end module

