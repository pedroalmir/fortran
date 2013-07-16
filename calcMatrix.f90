PROGRAM calcMatrix
!---------------------------------------------------------------------
!
!  This program ...
!
!  @author Flavio and Pedro Almir
!
!---------------------------------------------------------------------
IMPLICIT NONE
! Declare pedigree
REAL*8, DIMENSION(:,:), ALLOCATABLE :: pedigree
! Declare L (What is this?)
REAL*8, ALLOCATABLE, DIMENSION (:,:) :: L
! Declare A (What is this?)
REAL*8, ALLOCATABLE, DIMENSION (:,:) :: A
! Declare u (What is this?)
REAL*8, ALLOCATABLE, DIMENSION (:) :: vectorU
! Declare v (What is this?)
REAL*8, ALLOCATABLE, DIMENSION (:) :: vectorV
! Declare auxiliar variables
INTEGER :: i, j, k, row, col, animalID, fatherID, motherID
REAL*8 :: aux, pow, valueAux
INTEGER, DIMENSION (:,:), ALLOCATABLE :: matrix
! Declare the number of animals
INTEGER :: animals
! Declare the number of columns
INTEGER :: columns
! File Name
character (40) :: inputFileName
character (40) :: outputFileName
!---------------------------------------------------------------!
!							Commands							!
!---------------------------------------------------------------!
! Default file
inputFileName = 'input.txt'
outputFileName = 'output.txt'
! Default value = 3
columns = 3
! Get from file the number of animals. It's first line
open(unit = 11, file = inputFileName, status = 'old', action = 'read')
read(11, *) animals
! Initialization of variables
allocate(matrix(animals, columns))
! Read input file
write(*,'(/, 3X, A)') "Input:"
do row = 1, animals
	read(11, *) (matrix(row, col), col = 1, columns)
	write(*, "(100(3X, I3.1))") (matrix(row, col), col = 1, size(matrix, 2))
end do
! Initialization of variables
allocate (pedigree(animals + 1, columns + 1))
allocate (L(animals + 1, animals + 1))
allocate (A(animals + 1, animals + 1))
allocate (vectorU(animals + 1))
allocate (vectorV(animals + 1))
! Initilialization of pedigree
do row = 1, animals
	do col = 1, columns
		pedigree(row+1, col) = matrix(row, col)
	end do
end do
! Print Pedigree
write(*,'(/, 3X, A)') "Pedigree:"
do i=1, size(pedigree,1)
	write(*,"(100(3X,F12.10))") (pedigree(i,j), j=1, size(pedigree,2))
end do
            
! Write out the number of animals
! write(*,100) "The number of animais is ", animals, "."
!! 100 format (A, I6.2, A)

do i = 1, animals
	do j = 1, animals
		if(i <= j) then
			if(i == j) then
				if((pedigree(j,1) == 0.0) .AND. (pedigree(j,2) == 0.0)) then
					L(j,i) = 1.0
				else if((pedigree(j,1) > 0.0) .AND. (pedigree(j,2) > 0.0)) then
					aux = (1.0 - (0.25 * (vectorU(INT(pedigree(j, 1))) + vectorU(INT(pedigree(j, 1))))));
					pow = aux ** 0.5
					L(j,i) = pow
				else if((pedigree(j,1) > 0.0) .AND. (pedigree(j,2) == 0.0)) then
					aux = (1.0 - (0.25 * (vectorU(INT(pedigree(j, 1))))));
					pow = aux ** 0.5
					L(j,i) = pow
				else if((pedigree(j,1) == 0.0) .AND. (pedigree(j,2) > 0.0)) then
					aux = (1.0 - (0.25 * (vectorU(INT(pedigree(j, 2))))))
					pow = aux ** 0.5
					L(j,i) = pow
				end if
				vectorV(j) = L(j, i)
				vectorU(j) = vectorU(j) + (vectorV(j) * vectorV(j))
			else
				if(((i > pedigree(j,1)) .AND. (i > pedigree(j,2))) .OR. ((pedigree(j,1) == 0) .AND. (pedigree(j,2) == 0))) then
					L(j,i) = 0.0
				else if (i < pedigree(j,1) .AND. i < pedigree(j,2)) then
					L(j,i) = (0.5 * vectorV(INT(pedigree(j,1))) + (0.5 * vectorV(INT(pedigree(j,2)))))
				else if (i < pedigree(j,1) .AND. i == pedigree(j,2)) then
					L(j,i) = (0.5 * vectorV(INT(pedigree(j,1))) + (0.5 * vectorV(INT(pedigree(j,2)))))
				else if (i == pedigree(j,1) .AND. i < pedigree(j,2)) then
					L(j,i) = (0.5 * vectorV(INT(pedigree(j,1))) + (0.5 * vectorV(INT(pedigree(j,2)))))
				else if ((i > pedigree(j,2)) .AND. (i <= pedigree(j,1))) then
					L(j,i) = 0.5 * vectorV(INT(pedigree(j,1)))
				else if ((i > pedigree(j,1)) .AND. (i <= pedigree(j,2))) then
					L(j,i) = 0.5 * vectorV(INT(pedigree(j,2)))
				end if
				vectorV(j) = L(j, i);
				vectorU(j) = vectorU(j) + (vectorV(j) * vectorV(j))
			end if
		end if
	end do
end do

valueAux = 0.0;
do i = 1, animals
	do j = 1, animals
		do k = 1, animals
			valueAux = valueAux + (L(i, k) * L(j, k))
		end do
		A(i, j) = valueAux
		valueAux = 0.0
	end do
end do

! Print MatrixA
write(*,'(/, 3X, A)') "Matrix A:"
do i=1, size(A,1)-1
	write(*,"(100(3X, F12.10))") (A(i,j), j=1, size(A,2)-1)
end do

! Print Coeficiente de Parentesco
write(*, '(/, 3X, A, F12.10)') "Coeficiente de parentesco: ", A(1, animals)

open(unit=12, file=outputFileName, ACTION="write", STATUS="replace")
! Write Input
write(12,'(/, 3X, A)') "Input:"
do row = 1, animals
	write(12, "(100(3X, I3.1))") (matrix(row, col), col = 1, size(matrix, 2))
end do
! Write Pedigree
write(12,'(/, 3X, A)') "Pedigree:"
do i=1, size(pedigree,1)
	write(12,"(100(3X,F12.10))") (pedigree(i,j), j=1, size(pedigree,2))
end do
! Write MatrixA
write(12,'(/, 3X, A)') "Matrix A:"
do i=0, size(A,1)-1
	write(12,"(100(3X, F12.10))") (A(i,j), j=0, size(A,2)-1)
end do
! Write Parentesco
write(12, '(/, 3X, A, F12.10)') "Coeficiente de parentesco: ", A(1, animals)

deallocate(pedigree)
deallocate(L)
deallocate(A)
deallocate(vectorU)
deallocate(vectorV)

END PROGRAM calcMatrix
