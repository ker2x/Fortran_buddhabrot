PROGRAM MANDEL
USE omp_lib
IMPLICIT NONE

CHARACTER ( len = 255 ), PARAMETER :: filename = 'buddhabrot.ppm'
INTEGER, PARAMETER :: file_out_unit = 10


INTEGER,   PARAMETER :: n_max=100
INTEGER,   PARAMETER :: grid_resolution = 800
INTEGER,   PARAMETER :: zpower = 2 
INTEGER*8, PARAMETER :: batchSize = 50000000
REAL,      PARAMETER :: escapeOrbit = 4
REAL,      PARAMETER :: xmin = -1.3, xmax = 2.0, ymin = -1.3, ymax =1.3 

REAL,   PARAMETER :: intensityR = 255.
REAL,   PARAMETER :: intensityG = 255.
REAL,   PARAMETER :: intensityB = 255.

!Track pixel exposure by color
INTEGER :: exposureRMap(grid_resolution, grid_resolution)
INTEGER :: exposureGMap(grid_resolution, grid_resolution)
INTEGER :: exposureBMap(grid_resolution, grid_resolution)

INTEGER :: maxRExposure, maxGExposure, maxBExposure
INTEGER :: minRExposure, minGExposure, minBExposure
INTEGER :: maxExposure, minExposure

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

DO i=1, batchSize
  CALL RANDOM_NUMBER(x)
  CALL RANDOM_NUMBER(y)
  z = CMPLX(x*4. - 2. ,y*4. - 2.) !choose a random point on complex plane
  IF (notInMSet(z, n_max)) THEN   !if it espace out of the mandelbrot set
    c = z                        !then
    DO iter=1, n_max              !iterate and plot orbit
      z = z**zpower + c                !mandelbrot formula : Z = Z²+C
      IF(ABS(z) < escapeOrbit) THEN   !usefull when n_max > 1000
        TempX = INT(grid_resolution * (REAL(z) + xmax) / (xmax - xmin)) 
        TempY = INT(grid_resolution * (AIMAG(z) + ymax) / (ymax - ymin))
        IF((TempX > 0) .AND. (TempX < grid_resolution) .AND. (TempY > 0) .AND. (TempY < grid_resolution)) THEN
          IF((iter > 2)  .AND. (iter < 50)) THEN
            exposureRMap(TempX, TempY) = exposureRMap(TempX, TempY) + 1
          END IF
          IF((iter > 30) .AND. (iter < 70)) THEN
            exposureGMap(TempX, TempY) = exposureGMap(TempX, TempY) + 1
          END IF
          IF((iter > 50) .AND. (iter < 100)) THEN
            exposureBMap(TempX, TempY) = exposureBMap(TempX, TempY) + 1
          ENDIF
        END IF
      END IF  !(cabs(z)<4)
    END DO
  END IF
END DO

maxRExposure = MAXVAL(exposureRMap)
minRExposure = MINVAL(exposureRMap)
maxGExposure = MAXVAL(exposureGMap)
minGExposure = MINVAL(exposureGMap)
maxBExposure = MAXVAL(exposureBMap)
minBExposure = MINVAL(exposureBMap)
!write(*,*) maxRExposure, minRExposure, maxGExposure, minGExposure, maxBExposure, minBExposure

minExposure = MIN(minRExposure, minGExposure, minBExposure)
maxExposure = MAX(maxRExposure, maxGExposure, maxBExposure)

exposureRMap = exposureRMap - minExposure
exposureGMap = exposureGMap - minExposure
exposureBMap = exposureBMap - minExposure
exposureRMap = (exposureRMap / REAL(maxRExposure))*intensityR 
exposureGMap = (exposureGMap / REAL(maxGExposure))*intensityG 
exposureBMap = (exposureBMap / REAL(maxBExposure))*intensityB 

!maxRExposure = MAXVAL(exposureRMap)
!minRExposure = MINVAL(exposureRMap)
!maxGExposure = MAXVAL(exposureGMap)
!minGExposure = MINVAL(exposureGMap)
!maxBExposure = MAXVAL(exposureBMap)
!minBExposure = MINVAL(exposureBMap)
!write(*,*) maxRExposure, minRExposure, maxGExposure, minGExposure, maxBExposure, minBExposure

open ( unit = file_out_unit, file = filename, status = 'replace', &
       form = 'formatted', access = 'sequential')
write ( file_out_unit, '(a2)' ) 'P3'
write ( file_out_unit, '(i5,2x,i5)' ) grid_resolution, grid_resolution
write ( file_out_unit, '(i5)' ) INT(MAX(intensityR,intensityG,intensityB))

do ppm_i = 1, grid_resolution
do ppm_jlo = 1, grid_resolution, 4
    ppm_jhi = min ( ppm_jlo + 3, grid_resolution )
    write ( file_out_unit, '(12i5)' ) & 
    ( exposureRMap(ppm_i,ppm_j), exposureGMap(ppm_i,ppm_j),exposureBMap(ppm_i,ppm_j), ppm_j = ppm_jlo,ppm_jhi )
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
	  z = 0
	  n = 0

	  DO WHILE (ABS(z) < escapeOrbit .AND. (n < n_max))
	    z = z**zpower + c
	    n = n + 1
	  END DO

	  IF (n >= n_max) THEN
	    notInMset = .FALSE.
	  ELSE
	    notInMset = .TRUE.
	  END IF
	END FUNCTION notInMset


END
