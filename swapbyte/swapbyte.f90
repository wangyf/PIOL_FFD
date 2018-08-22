PROGRAM BYTE_SWAP
! we first need to know file name and file size
INTEGER(kind=8) :: FILE_SIZE, i
CHARACTER(LEN=20) :: FILENAME,filenameout
!real,allocatable,dimension(:) :: fieldio
real :: fieldio
write(*,'(A30)') 'Program: byte_swap '
write(*,'(A30)') 'Swap byte of the file required to name'
write(*,'(A30)') 'Please input the file name: -- '
READ(*,*) FILENAME

INQUIRE(FILE = FILENAME, SIZE=FILE_SIZE)

write(*,'(A10,I10,A5)') 'File size is ', file_size,' byte'
!allocate(fieldio(file_size/4))

open(1, file=filename, recl=4, form='unformatted', access='direct', &
        status='old' )
open(2, file=trim(filename)//'.byteswap', recl=4, form='unformatted', access='direct', &
        status='replace' )
do i = 1, file_size/4   
read(1,rec=i) fieldio
call SWAP_F4(fieldio)
write(2,rec=i) fieldio
write(*,'(F6.2,A3)') real(i)/file_size*4*100,'%'
end do
close(1)
close(2)
END PROGRAM BYTE_SWAP




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

