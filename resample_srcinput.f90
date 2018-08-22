program resample_srcinput
! This little program is used to resample 
! source input (slip/moment rate and geometry)
!
! Yongfei Wang
! Jun 24, 2018
Interface
	SUBROUTINE SWAPR42D(F)	
		REAL,INTENT(INOUT),DIMENSION(:,:) :: F
	end subroutine
end interface

INTEGER(kind=8) :: FILE_SIZE,iostat, i,j,it,sumit,iswap,&
                   nx,ny,nt,ix1,ix2,dx,iy1,iy2,dy,it1,it2,dt,&
			       nxr,nyr,ntr
CHARACTER(LEN=20) :: FILENAME,filenameout
real,allocatable,dimension(:,:) :: fieldio

write(*,'(A)') 'Program: resample_src_input '
write(*,'(A)') ''
write(*,'(A)') 'Please input the grid file name: xx'
READ(*,*) FILENAME
write(*,*) FILENAME
write(*,'(A)') 'Please input original dimension nx ny nt'
READ(*,*) nx,ny,nt
write(*,*) nx,ny,nt

INQUIRE(FILE = FILENAME, SIZE=FILE_SIZE)
if (FILE_SIZE/4 /= nx*ny*nt) then
	write(*,'(A,4I10)') 'nx or ny is wrong in number',nx,ny,nt,FILE_SIZE/4
	stop
end if

! read in allocatable arrays
write(*,'(A)') 'Please input if swap byte 0 no 1 yes'
READ(*,*) iswap
if (iswap == 1) then
	write(*,'(A)') iswap,'Swap Bytes'
else
	write(*,'(A)') iswap,'No Swap Bytes'
end if

write(*,'(A)') 'Please input resample parameter: ix1 ix2 dx iy1 iy2 dy it1 it2 dt'
READ(*,*) ix1, ix2, dx, iy1, iy2, dy, it1, it2, dt
write(*,'(A,9I8)') 'Confirm: ',ix1, ix2, dx, iy1, iy2, dy, it1, it2, dt
! compute new dimension
nxr = floor((ix2 - ix1 )/ real(dx) ) + 1
nyr = floor((iy2 - iy1 )/ real(dx) ) + 1
ntr = floor((it2 - it1 )/ real(dt) ) + 1

open(1, file=filename, recl=4*nx*ny, form='unformatted', access='direct', &
        status='old' )
open(2, file=trim(filename)//'.resampled', recl=4*nxr*nyr, form='unformatted', &
		access='direct', status='replace' )
		
allocate(fieldio(nx,ny))
sumit = 0
do it = it1,it2,dt
	sumit = sumit + 1
	write( 0, '(a)', advance='no' ) '.'
	if ( modulo( it, 50 ) == 0 .or. it == nt ) write( 0, '(i6)' ) it
	read(1,rec=it,iostat=ios) ((fieldio(i,j),i=1,nx),j=1,ny)
	if (iswap == 1) then
		call SWAPR42D(fieldio)
	end if
	write(2,rec=sumit,iostat=ios) ((fieldio(i,j),i=ix1,ix2,dx),j=iy1,iy2,dy)
end do

close(1)
close(2)

! finish processing of grid file

end program resample_srcinput

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


!***********************************************************************************************************************************
!  SWAP_I2
!
!  Swap bytes for a two-byte INTEGER.
!  After calling this subroutine, the input integer will be replaced by the output integer.
!***********************************************************************************************************************************

      SUBROUTINE SWAP_I2 (BYTE2)

      IMPLICIT NONE

      INTEGER(KIND=2), INTENT(IN OUT) :: BYTE2

      INTEGER(KIND=1), DIMENSION(2) :: BYTE_ARR, BYTE_ARR_TMP
      INTEGER :: I


      BYTE_ARR = TRANSFER (BYTE2, BYTE_ARR)

      BYTE_ARR_TMP = BYTE_ARR

      DO I = 1, 2
         BYTE_ARR(I) = BYTE_ARR_TMP(3-I)
      END DO

      BYTE2 = TRANSFER (BYTE_ARR, BYTE2)

      RETURN

      END SUBROUTINE SWAP_I2





!***********************************************************************************************************************************
!  SWAP_I4
!
!  Swap bytes for a four-byte INTEGER.
!  After calling this subroutine, the input integer will be replaced by the output integer.
!***********************************************************************************************************************************

      SUBROUTINE SWAP_I4 (BYTE4)

      IMPLICIT NONE

      INTEGER(KIND=4), INTENT(IN OUT) :: BYTE4

      INTEGER(KIND=1), DIMENSION(4) :: BYTE_ARR, BYTE_ARR_TMP
      INTEGER :: I


      BYTE_ARR = TRANSFER (BYTE4, BYTE_ARR)

      BYTE_ARR_TMP = BYTE_ARR

      DO I = 1, 4
         BYTE_ARR(I) = BYTE_ARR_TMP(5-I)
      END DO

      BYTE4 = TRANSFER (BYTE_ARR, BYTE4)

      RETURN

      END SUBROUTINE SWAP_I4





!***********************************************************************************************************************************
!  SWAP_F4
!
!  Swap bytes for a four-byte REAL.
!  After calling this subroutine, the input number will be replaced by the output number.
!***********************************************************************************************************************************

      SUBROUTINE SWAP_F4 (FLOAT4)

      IMPLICIT NONE

      REAL(KIND=4), INTENT(IN OUT) :: FLOAT4

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





!***********************************************************************************************************************************
!  SWAP_F8
!
!  Swap bytes for an eight-byte REAL.
!  After calling this subroutine, the input number will be replaced by the output number.
!***********************************************************************************************************************************

      SUBROUTINE SWAP_F8 (FLOAT8)

      IMPLICIT NONE

      REAL(KIND=8), INTENT(IN OUT) :: FLOAT8

      INTEGER(KIND=1), DIMENSION(8) :: BYTE_ARR, BYTE_ARR_TMP
      INTEGER :: I


      BYTE_ARR = TRANSFER (FLOAT8, BYTE_ARR)

      BYTE_ARR_TMP = BYTE_ARR

      DO I = 1, 8
         BYTE_ARR(I) = BYTE_ARR_TMP(9-I)
      END DO

      FLOAT8 = TRANSFER (BYTE_ARR, FLOAT8)

      RETURN

      END SUBROUTINE SWAP_F8

