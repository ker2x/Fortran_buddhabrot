PROGRAM MANDEL
USE omp_lib
use ifport
IMPLICIT NONE

CHARACTER ( len = 255 ), PARAMETER :: filename = 'buddhabrot.ppm'
INTEGER, PARAMETER :: file_out_unit = 10


! mix/max grid coordinate
REAL, PARAMETER :: xmin = -1.0, xmax = 2.0, ymin = -1.3, ymax =1.3 

! Arbitrary maximum # of iterations
INTEGER,   PARAMETER :: n_max=100
INTEGER,   PARAMETER :: grid_resolution = 1000
INTEGER*8, PARAMETER :: batchSize = 10000000

!Track pixel exposure by color
INTEGER :: exposureRMap(grid_resolution, grid_resolution)
INTEGER :: exposureGMap(grid_resolution, grid_resolution)
INTEGER :: exposureBMap(grid_resolution, grid_resolution)

INTEGER*8 :: i
REAL    :: x,y
COMPLEX :: z, c
INTEGER :: iter
INTEGER :: tempX, tempY

INTEGER :: ppm_i
INTEGER :: ppm_j
INTEGER :: ppm_jhi
INTEGER :: ppm_jlo

!initialize exposureMap to 1
exposureRMap = 1
exposureGMap = 1
exposureBMap = 1


x = RANDOM(1)
DO i=1, batchSize
  x = RANDOM(0) * 4. - 2.        !old code was : x = RANDOM(0) * (xmax-xmin) + xmin
  y = RANDOM(0) * 4. - 2.        !old code was : y = RANDOM(0) * (ymax-ymin) + ymin
  z = CMPLX(x,y)                 !choose a random point on complex plane
  IF (notInMSet(z, n_max)) THEN  !if it espace out of the mandelbrot set
    c = z                        !then
    DO iter=1, n_max             !iterate and plot orbit
      z = z*z + c                !mandelbrot formula : Z = Z²+C
      IF(CABS(z) < 4) THEN
        TempX = INT(grid_resolution * (REAL(z) + xmax) / (xmax - xmin)) 
        TempY = INT(grid_resolution * (AIMAG(z) + ymax) / (ymax - ymin))
        IF((TempX > 0) .AND. (TempX < grid_resolution) .AND. (TempY > 0) .AND. (TempY < grid_resolution)) THEN
          IF(iter < 20) THEN
              exposureRMap(TempX, TempY) = exposureRMap(TempX, TempY) + 1
            ELSE IF(iter < 50) THEN
              exposureBMap(TempX, TempY) = exposureBMap(TempX, TempY) + 1
            ELSE IF(iter < 100) THEN
              exposureGMap(TempX, TempY) = exposureGMap(TempX, TempY) + 1
          ENDIF
        END IF
      END IF
    END DO
  END IF
END DO


exposureRMap = exposureRMap - (MINVAL(exposureRMap) )
exposureGMap = exposureGMap - (MINVAL(exposureGMap) )
exposureBMap = exposureBMap - (MINVAL(exposureBMap) )
exposureRMap = exposureRMap / (MAXVAL(exposureRMap)/255. )
exposureGMap = exposureGMap / (MAXVAL(exposureGMap)/128. )
exposureBMap = exposureBMap / (MAXVAL(exposureBMap)/255. )
!write(*,*) MAXVAL(exposureRMap) , MAXVAL(exposureGMap) , MAXVAL(exposureBMap)
!write(*,*) MINVAL(exposureRMap) , MINVAL(exposureGMap) , MINVAL(exposureBMap)


open ( unit = file_out_unit, file = filename, status = 'replace', &
       form = 'formatted', access = 'sequential')
write ( file_out_unit, '(a2)' ) 'P3'
write ( file_out_unit, '(i5,2x,i5)' ) grid_resolution, grid_resolution
write ( file_out_unit, '(i5)' ) 255

do ppm_i = 1, grid_resolution
do ppm_jlo = 1, grid_resolution, 4
    ppm_jhi = min ( ppm_jlo + 3, grid_resolution )
    write ( file_out_unit, '(12i5)' ) ( exposureRMap(ppm_i,ppm_j), exposureGMap(ppm_i,ppm_j),exposureBMap(ppm_i,ppm_j), ppm_j = ppm_jlo,ppm_jhi )
  end do
end do


close ( unit = file_out_unit )


CONTAINS

	PURE FUNCTION notInMset(c, n_max)
	  COMPLEX, INTENT(IN) :: c
	  INTEGER, INTENT(IN) :: n_max
	  INTEGER             :: n 
	  COMPLEX             :: z
	  LOGICAL             :: notInMSet
	  z = c
	  n = 0

	  DO WHILE (ABS(z) < 4.0 .AND. (n < n_max))
	    z = z*z + c
	    n = n + 1
	  END DO

	  IF (n >= n_max) THEN
	    notInMset = .FALSE.
	  ELSE
	    notInMset = .TRUE.
	  END IF
	END FUNCTION notInMset


END
