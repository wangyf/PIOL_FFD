! Miscellaneous utilities
module m_util
implicit none
real, parameter :: pi = 4.*atan(1.)
real, parameter :: deg2rad = pi/180
real, parameter :: epsilon = 1e-3
real, parameter :: resolution = 1e-6
contains

! Array reciprocal
subroutine invert( f )
real, intent(inout) :: f(:,:,:)
integer :: n(3), j, k, l
n = (/ size(f,1), size(f,2), size(f,3) /)
do l = 1, n(3)
do k = 1, n(2)
do j = 1, n(1)
    if ( f(j,k,l) /= 0.0 ) f(j,k,l) = 1.0 / f(j,k,l)
end do
end do
end do
end subroutine

! Squared distance to x0
subroutine radiusR( r, x, x0, i1, i2 )
real, intent(out) :: r(:,:,:)
real, intent(in) :: x(:,:,:,:), x0(3)
integer, intent(in) :: i1(3), i2(3)
integer :: n(3), j, k, l
n = (/ size(r,1), size(r,2), size(r,3) /)
if ( any( i1 < 1 .or. i2 > n ) ) stop 'error in radius'
do l = i1(3), i2(3)
do k = i1(2), i2(2)
do j = i1(1), i2(1)
    r(j,k,l) = &
    ( x(j,k,l,1) - x0(1) ) * ( x(j,k,l,1) - x0(1) ) + &
    ( x(j,k,l,2) - x0(2) ) * ( x(j,k,l,2) - x0(2) ) + &
    ( x(j,k,l,3) - x0(3) ) * ( x(j,k,l,3) - x0(3) )
end do
end do
end do
end subroutine

! Average of local eight values
subroutine average( f2, f1, i1, i2, d )
real, intent(out) :: f2(:,:,:)
real, intent(in) :: f1(:,:,:)
integer, intent(in) :: i1(3), i2(3), d
integer :: n(3), j, k, l
n = (/ size(f1,1), size(f1,2), size(f1,3) /)
if ( any( i1 < 1 .or. i2 > n ) ) stop 'error in average'
do l = i1(3), i2(3)
do k = i1(2), i2(2)
do j = i1(1), i2(1)
    f2(j,k,l) = 0.125 * &
    ( f1(j,k,l) + f1(j+d,k+d,l+d) &
    + f1(j,k+d,l+d) + f1(j+d,k,l) &
    + f1(j+d,k,l+d) + f1(j,k+d,l) &
    + f1(j+d,k+d,l) + f1(j,k,l+d) )
end do
end do
end do
call set_halo( f2, 0.0, i1, i2 )
end subroutine

! Set array to real value outside specified region
subroutine set_halo( f, r, i1, i2 )
real, intent(inout) :: f(:,:,:)
real, intent(in) :: r
integer, intent(in) :: i1(3), i2(3)
integer :: n(3), i3(3), i4(3)
n = (/ size(f,1), size(f,2), size(f,3) /)
i3 = min( i1, n + 1 )
i4 = max( i2, 0 )
if ( n(1) > 1 ) f(:i3(1)-1,:,:) = r
if ( n(2) > 1 ) f(:,:i3(2)-1,:) = r
if ( n(3) > 1 ) f(:,:,:i3(3)-1) = r
if ( n(1) > 1 ) f(i4(1)+1:,:,:) = r
if ( n(2) > 1 ) f(:,i4(2)+1:,:) = r
if ( n(3) > 1 ) f(:,:,i4(3)+1:) = r
end subroutine

! L2 vector norm
subroutine vector_norm( f, w, i1, i2, di )
real, intent(out) :: f(:,:,:)
real, intent(in) :: w(:,:,:,:)
integer, intent(in) :: i1(3), i2(3), di(3)
integer :: n(3), j, k, l
n = (/ size(f,1), size(f,2), size(f,3) /)
if ( any( i1 < 1 .or. i2 > n ) ) stop 'error in vector_norm'
do l = i1(3), i2(3), di(3)
do k = i1(2), i2(2), di(2)
do j = i1(1), i2(1), di(1)
    f(j,k,l) = &
    w(j,k,l,1) * w(j,k,l,1) + &
    w(j,k,l,2) * w(j,k,l,2) + &
    w(j,k,l,3) * w(j,k,l,3)
end do
end do
end do
end subroutine

! Frobenius tensor norm - much faster than L2 norm for tensors
subroutine tensor_norm( f, w1, w2, i1, i2, di )
real, intent(out) :: f(:,:,:)
real, intent(in) :: w1(:,:,:,:), w2(:,:,:,:)
integer, intent(in) :: i1(3), i2(3), di(3)
integer :: n(3), j, k, l
n = (/ size(f,1), size(f,2), size(f,3) /)
if ( any( i1 < 1 .or. i2 > n ) ) stop 'error in tensor_norm'
do l = i1(3), i2(3), di(3)
do k = i1(2), i2(2), di(2)
do j = i1(1), i2(1), di(1)
    f(j,k,l) = &
    w1(j,k,l,1) * w1(j,k,l,1) + &
    w1(j,k,l,2) * w1(j,k,l,2) + &
    w1(j,k,l,3) * w1(j,k,l,3) + &
    ( w2(j,k,l,1) * w2(j,k,l,1) &
    + w2(j,k,l,2) * w2(j,k,l,2) &
    + w2(j,k,l,3) * w2(j,k,l,3) ) * 2.
end do
end do
end do
end subroutine

! In-place linear interpolation 
subroutine interpolate( f, i3, i4, di )
real, intent(inout) :: f(:,:,:)
integer, intent(in) :: i3(3), i4(3), di(3)
integer :: i1(3), i2(3), n(3), i, j, k, l, d
real :: h1, h2
n = (/ size(f,1), size(f,2), size(f,3) /)
i1 = i3
i2 = i4
where( i1 < 1 ) i1 = i1 + (-i1 / di + 1) * di
where( i2 > n ) i2 = i1 + (n - i1) / di * di
d = di(1)
do i = 1, d - 1
    h1 = 1.0 / d * i
    h2 = 1.0 / d * ( d - i )
    do l = i1(3), i2(3), di(3)
    do k = i1(2), i2(2), di(2)
    do j = i1(1), i2(1) - d, d
        f(j+i,k,l) = h1 * f(j,k,l) + h2 * f(j+d,k,l)
    end do
    end do
    end do
end do
d = di(2)
do i = 1, d - 1
    h1 = 1.0 / d * i
    h2 = 1.0 / d * ( d - i )
    do l = i1(3), i2(3), di(1)
    do k = i1(2), i2(2) - d, d
    do j = i1(1), i2(1)
        f(j,k+i,l) = h1 * f(j,k,l) + h2 * f(j,k+d,l)
    end do
    end do
    end do
end do
d = di(3)
do i = 1, d - 1
    h1 = 1.0 / d * i
    h2 = 1.0 / d * ( d - i )
    do l = i1(3), i2(3) - d, d
    do k = i1(2), i2(2)
    do j = i1(1), i2(1)
        f(j,k,l+i) = h1 * f(j,k,l) + h2 * f(j,k,l+d)
    end do
    end do
    end do
end do
end subroutine

! Time function
real function time_function( tfunc, tm, dt, period )
character(*), intent(in) :: tfunc
real, intent(in) :: tm, dt, period
real, parameter :: pi = 3.14159265
real :: t
time_function = 0.0
select case( tfunc )
case( 'ramp')
    if ( abs(tm) < period ) then
        time_function = 0.5*(1-cos(pi*tm/period))
    else
        time_function = 1.0
    end if
!    if ( tm < period ) time_function = 0.5*pi/period*sin(pi*tm/period) 
case( 'const'  )
    time_function = 1.0
case( 'delta'  )
    if ( abs( tm ) < 0.25 * dt ) time_function = 1.0
case( 'brune' )
    time_function = -exp( -tm / period ) / period * (tm + period) + 1.0
case( 'dbrune' )
    time_function =  exp( -tm / period ) / period ** 2.0 * tm
case( 'ddbrune' )
    time_function = -exp( -tm / period ) / period ** 3.0 * (tm - period) 
case( 'gaussian' )
    t = ( tm - 4.0 * period ) / period
    time_function = exp( -0.5 * t * t ) / ( period * sqrt( 2.0 * pi ) )
case( 'dgaussian', 'ricker1' )
    t = tm - period
    time_function = t * exp( -2.0 * (pi * t / period) ** 2.0 )
case( 'ddgaussian', 'ricker2' )
    t = ( pi * (tm - period) / period ) ** 2.0
    time_function = (1.0 - 2.0 * t) * exp( -t )
case default
    write( 0, * ) 'invalid time func: ', trim( tfunc )
    stop
end select
end function

! Timer
real function timer( i )
integer, intent(in) :: i
integer, save :: clock0, clockrate, clockmax
integer(8), save :: timers(8)
integer :: clock1
if ( i == 0 ) then
    call system_clock( clock0, clockrate, clockmax )
    timer = 0
    timers = 0
else
    call system_clock( clock1 )
    timers = timers + clock1 - clock0
    if ( clock0 > clock1 ) timers = timers + clockmax
    clock0 = clock1
    timer = real( timers(i) ) / real( clockrate )
    timers(:i) = 0
end if
end function

!user customized part
!****************************************************
! start


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!byteswap for converting between big-endian and little-endian
	SUBROUTINE SWAPR42D(F)
	IMPLICIT NONE
	REAL,INTENT(INOUT),DIMENSION(:,:) :: F
	INTEGER :: N(2), I, J
	
	N(1) = SIZE(F,1)
	N(2) = SIZE(F,2)
	
	DO J = 1, N(2) 
	  DO I = 1, N(1)
	    CALL SWAP_F4(F(I,J))
	  END DO
	END DO
	
	END SUBROUTINE


	SUBROUTINE SWAPR43D(F)
	IMPLICIT NONE
	REAL,INTENT(INOUT),DIMENSION(:,:,:) :: F
	INTEGER :: N(3), I, J, K
	
	N(1) = SIZE(F,1)
	N(2) = SIZE(F,2)
	N(3) = SIZE(F,3)
	
	DO K = 1, N(3)
	 DO J = 1, N(2) 
	  DO I = 1, N(1)
	    CALL SWAP_F4(F(I,J,K))
	  END DO
	 END DO
	END DO
	
	END SUBROUTINE
	
	SUBROUTINE SWAPR44D(F)
	IMPLICIT NONE
	REAL,INTENT(INOUT),DIMENSION(:,:,:,:) :: F
	INTEGER :: N(4), I, J, K, L
	
	N(1) = SIZE(F,1)
	N(2) = SIZE(F,2)
	N(3) = SIZE(F,3)
	N(4) = SIZE(F,4)
	
	DO L = 1, N(4)
	DO K = 1, N(3)
	 DO J = 1, N(2) 
	  DO I = 1, N(1)
	    CALL SWAP_F4(F(I,J,K,L))
	  END DO
	 END DO
	END DO
	END DO
	END SUBROUTINE
	
SUBROUTINE SWAP_F4 (FLOAT4)

      IMPLICIT NONE

      REAL, INTENT(IN OUT) :: FLOAT4

      INTEGER(KIND=1), DIMENSION(4) :: BYTE_ARR, BYTE_ARR_TMP
      INTEGER :: I


      BYTE_ARR = TRANSFER (FLOAT4, BYTE_ARR)

      BYTE_ARR_TMP = BYTE_ARR

      DO I = 1, 4
         BYTE_ARR(I) = BYTE_ARR_TMP(5-I)
      END DO

      FLOAT4 = TRANSFER (BYTE_ARR, FLOAT4)

      RETURN

END SUBROUTINE SWAP_F4

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!------------------
! cartesian coordinate into spherical coordinate
subroutine cat2sph(x,y,z,r,theta,phi)
real,intent(in) :: x, y, z
real,intent(out) :: r,phi,theta
real,parameter :: pi = 4.*atan(1.) 

r = sqrt(x**2 + y**2 + z**2)
theta = pi/2 - acos(z / r) !theta ranges from -pi/2 to pi/2
phi = atan2(y, x)          !phi ranges from -pi to pi
end subroutine



integer function find_first(buf,n,operatr,val)
real,intent(in) :: buf(:)
integer,intent(in) :: n
character(*),intent(in) :: operatr
real,intent(in) :: val

integer :: i
find_first = -9999
do i = 1,n
    select case( operatr )
    case('>=')
      if(buf(i)>= val) then 
        find_first = i
        exit
      end if
    case('<=')
      if(buf(i)<= val) then 
        find_first = i
        exit
      end if    
    case default
      write(0,*) 'Error: invalid type in find_first'
      stop
    end select
end do
end function

integer function find_last(buf,n,operatr,val)
real,intent(in) :: buf(:)
integer,intent(in) :: n
character(*),intent(in) :: operatr
real,intent(in) :: val

integer :: i
find_last = -9999
do i = 1,n
    select case( operatr )
    case('>=')
      if(buf(n-i+1)>= val) then 
        find_last = n-i+1
        exit
      end if
    case('<=')
      if(buf(n-i+1)<= val) then 
        find_last = n-i+1
        exit
      end if    
    case default
      write(0,*) 'Error: invalid type in find_last'
      stop
    end select
end do
end function

subroutine tinti(npts_stf,dt,T0,Tp,Tr,D,stfo)
implicit none
integer,intent(in) :: npts_stf
real,intent(in) :: dt,T0,Tp,Tr,D
real,intent(out) :: stfo(:)

integer :: i,it0,j

! tinti parameter
real :: c2,k,ts
real,dimension(npts_stf) :: t,c1,c3,c4,c5,c6,stf
!
stfo = 0.
do i = 1, npts_stf
  t(i) = (i - 1 ) * dt
end do

   ts=Tp
   c2=(3.0/8.0)*pi*(Tr**2.0)
   k=2.0/(pi*Tr*ts**2.0)

   c1=((0.5*t+0.25*Tr)*sqrt(t*(Tr-t))+(t*Tr-(Tr**2.0))*asin(sqrt(t/Tr))     &
      -0.75*(Tr**2.0)*atan(sqrt((Tr-t)/t)))

   c3=((ts-t-0.5*Tr)*sqrt((t-ts)*(Tr-t+ts))+Tr*(2.0*Tr-2.0*t+2.0*ts)        &
      *(asin(sqrt((t-ts)/Tr)))+(3.0/2.0)*Tr**2.0*(atan(sqrt((Tr-t+ts)       &
      /(t-ts))))) 

   c4=((-1.0*ts+0.5*t+0.25*Tr)*sqrt((t-2.0*ts)*(Tr-t+2.0*ts))+Tr*((-1.0*Tr) &
      +t-2.0*ts)*asin(sqrt((t-2.0*ts)/Tr))-0.75*Tr**2.0*atan(sqrt((Tr-t     &
      +2.0*ts)/(t-2.0*ts)))) 

   c5=((pi/2.0)*Tr*(t-Tr))

   c6=((pi/2.0)*Tr*(2.0*ts-t+Tr))
   
   ! here below create Yoffe's stf
   if (Tr .GT. ts*2.0) then
      do j=1,npts_stf
         if (t(j) .LE. 0.0) then
            stf(j)=0.0
         else if (0.0 .LT. t(j) .AND. t(j) .LT. ts) then
            stf(j)=D*k*(c1(j)+c2)
            else if (ts .LE. t(j) .AND. t(j) .LT. 2.0*ts) then
                stf(j)=D*k*(c1(j)-c2+c3(j))
            else if (2.0*ts .LE. t(j) .AND. t(j) .LT. Tr) then
                stf(j)=D*k*(c1(j)+c3(j)+c4(j))
            else if (Tr .LE. t(j) .AND. t(j) .LT. Tr+ts) then
                stf(j)=D*k*(c5(j)+c3(j)+c4(j))
            else if (Tr+ts .LE. t(j) .AND. t(j) .LT. Tr+2.0*ts) then
                stf(j)=D*k*(c4(j)+c6(j))
            else if (Tr+2.0*ts .LE. t(j)) then
                stf(j)=0.0
         end if
      enddo
   endif
   
   if ((ts .LT. Tr) .AND. (Tr .LT. 2.0*ts)) then
      do j=1,npts_stf
         if (t(j) .LE. 0.0) then
                stf(j)=0.0
             else if (0.0 .LT. t(j) .AND. t(j) .LT. ts) then
                stf(j)=D*k*(c1(j)+c2)
             else if (ts .LE. t(j) .AND. t(j) .LT. Tr) then
                stf(j)=D*k*(c1(j)-c2+c3(j))
             else if (Tr .LE. t(j) .AND. t(j) .LT. 2.0*ts) then
                stf(j)=D*k*(c5(j)+c3(j)-c2)
             else if (2.0*ts .LE. t(j) .AND. t(j) .LT. Tr*ts) then
                stf(j)=D*k*(c5(j)+c3(j)+c4(j))
             else if (Tr*ts .LE. t(j) .AND. t(j) .LT. Tr+2.0*ts) then
                stf(j)=D*k*(c4(j)+c6(j))
             else if (Tr+2.0*ts .LE. t(j)) then
                stf(j)=0.0
         endif
      enddo
   endif

! shift T0
it0 = nint(T0/dt)-1
do i = 1, npts_stf-it0
  stfo(i+it0) = stf(i) 
  !if(isnan(stf(i))) write(0,*) i,tp,tr,t0
end do
end subroutine

integer function nancheck4d(f)
implicit none
real,intent(in) :: f(:,:,:,:)
integer :: i,j,k,l,i1,i2,i3,i4
i1 = size(f,1)
i2 = size(f,2)
i3 = size(f,3)
i4 = size(f,4)
nancheck4d = 1

do i = 1, i1
  do j = 1, i2
    do k = 1, i3
      do l = 1, i4
        if(isnan(f(i,j,k,l))) then
          nancheck4d = 0
          write(0,*) 'Nan Found',i,j,k,l
          cycle
        end if
      end do
    end do
  end do
end do
end function

integer function nancheck3d(f)
implicit none
real,intent(in) :: f(:,:,:)
integer :: i,j,k,i1,i2,i3
i1 = size(f,1)
i2 = size(f,2)
i3 = size(f,3)
nancheck3d = 1

do i = 1, i1
  do j = 1, i2
    do k = 1, i3
        if(isnan(f(i,j,k))) then
          nancheck3d = 0
          write(0,*) 'Nan Found',i,j,k
          cycle
        end if
    end do
  end do
end do
end function

real function trapz(tms,nt,dt)
real,intent(in) :: tms(:),dt
integer,intent(in) :: nt

integer :: i

trapz = 0
do i = 2, nt
  trapz = trapz + ( tms(i) + tms(i-1) ) / 2. * dt
end do
end function 

subroutine bisection(c1,c2,c3,c4,ii)
real,intent(in) :: c1,c2,c3,c4
real,intent(out) :: ii
   REAL            :: Left,  fLeft
   REAL            :: Right, fRight
   REAL            :: Root
   REAL,parameter :: pi = 4.*atan(1.)

! Initial guess   
   Left = 0.
   Right = pi/2
   
   fLeft  = Funct(Left,c1,c2,c3,c4)                 ! compute their function values
   fRight = Funct(Right,c1,c2,c3,c4)
    
   IF (fLeft*fRight > 0.0)  THEN
      WRITE(*,*)  '*** ERROR: f(Left)*f(Right) must be negative ***'
   ELSE
      Root = Solve(Left, Right, resolution,c1,c2,c3,c4)
      !write(0,'(A,E18.10)') 'Root is ',root
   END IF
   ii = Root
end subroutine

REAL FUNCTION  Solve(Left, Right, Tolerance,c1,c2,c3,c4)
IMPLICIT  NONE
REAL, INTENT(IN) :: Left, Right, Tolerance,c1,c2,c3,c4 !c1-c4 parameter in function
REAL             :: a, Fa, b, Fb, c, Fc

      a = Left                          ! save Left and Right
      b = Right

      Fa = Funct(a,c1,c2,c3,c4)                     ! compute the function values
      Fb = Funct(b,c1,c2,c3,c4)

      IF (ABS(Fa) < Tolerance) THEN     ! if f(a) is already small
         Solve = a                      ! then a is a root
      ELSE IF (ABS(Fb) < Tolerance) THEN     ! is f(b) is small
         Solve = b                      ! then b is a root
      ELSE                              ! otherwise,
         DO                             ! iterate ....
            c  = (a + b)/2.0            !   compute the middle point
            Fc = Funct(c,c1,c2,c3,c4)               !   and its function value
            IF (ABS(Fc) < Tolerance) THEN    ! is it very small?
               Solve = c                ! yes, c is a root
               EXIT
            ELSE IF (Fa*Fc < 0.0) THEN  ! do f(a)*f(c) < 0 ?
               b  = c                   ! replace b with c
               Fb = Fc                  ! and f(b) with f(c)
            ELSE                        ! then f(c)*f(b) < 0 holds
               a  = c                   ! replace a with c
               Fa = Fc                  ! and f(a) with f(c)
            END IF

            if (abs(a-b) < tolerance) then
            ! interpolate to root 
              solve = (fa*b-fb*a)/(fa-fb) 
              exit
            end if 
         END DO                         ! go back and do it again
      END IF
END FUNCTION  Solve

real function Funct(kk,c1,c2,c3,c4)
    real, intent(in) :: kk,c1,c2,c3,c4
    
    Funct = c1-abs(c2*tan(kk))-abs(c3*tan(asin(c4*sin(kk))))
    
end function

!************************************************

end module

