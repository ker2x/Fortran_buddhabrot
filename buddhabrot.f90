PROGRAM MANDEL
USE omp_lib
IMPLICIT NONE
 
CHARACTER ( len = 255 ), PARAMETER :: filename = 'buddhabrot.ppm'
INTEGER, PARAMETER :: file_out_unit = 10
 
 
INTEGER, PARAMETER :: n_max=100
INTEGER, PARAMETER :: grid_resolution = 512
INTEGER, PARAMETER :: zpower = 2
INTEGER*8, PARAMETER :: batchSize = 10000000
REAL, PARAMETER :: escapeOrbit = 2
REAL, PARAMETER :: xmin = -1.0, xmax = 2.0, ymin = -1.3, ymax =1.3
 
REAL, PARAMETER :: intensityR = 2048.
REAL, PARAMETER :: intensityG = 2048.
REAL, PARAMETER :: intensityB = 2048.
 
!Track pixel exposure by color
INTEGER :: exposureRMap(grid_resolution, grid_resolution)
INTEGER :: exposureGMap(grid_resolution, grid_resolution)
INTEGER :: exposureBMap(grid_resolution, grid_resolution)
 
INTEGER :: maxRExposure, maxGExposure, maxBExposure
INTEGER :: minRExposure, minGExposure, minBExposure
INTEGER :: maxExposure, minExposure
 
INTEGER*8 :: i
REAL :: x,y
COMPLEX :: z, c
INTEGER :: iter
INTEGER :: tempX, tempY, tempYm
 
INTEGER :: ppm_i
INTEGER :: ppm_j
INTEGER :: ppm_jhi
INTEGER :: ppm_jlo
 
!initialize exposureMap to 1
exposureRMap = 0
exposureGMap = 0
exposureBMap = 0
 
!$OMP PARALLEL DO DEFAULT(PRIVATE),SHARED(exposureRMap, exposureGMap, exposureBMap)
DO i=1, batchSize
  CALL RANDOM_NUMBER(x)
  CALL RANDOM_NUMBER(y)
  z = CMPLX(0,0)
  !c = CMPLX(x*2.5 - 2 ,y*2.6 - 1.3) !choose a random point on complex plane
  c = CMPLX(x*2.5 - 2 ,y*1.3) !choose a random point on complex plane
  IF (notInMSet(c, n_max)) THEN !if it espace out of the mandelbrot set
    DO iter=1, n_max !iterate and plot orbit
      z = z**zpower + c !mandelbrot formula : Z = Z²+C
      !IF(ABS(z) < escapeOrbit) THEN
        TempX = INT(grid_resolution * (REAL(z) + xmax) / (xmax - xmin))
        TempY = INT(grid_resolution * (AIMAG(z) + ymax) / (ymax - ymin))
        TempYm = INT(grid_resolution/2 - (TempY - grid_resolution/2))
        IF((TempX > 0) .AND. (TempX < grid_resolution) .AND. (TempY > 0) .AND. (TempY < grid_resolution)) THEN
          IF((iter > 2) .AND. (iter < 10)) THEN
            exposureRMap(TempX, TempY)  = exposureRMap(TempX, TempY) + 1
            exposureRMap(TempX, TempYm) = exposureRMap(TempX, TempYm) + 1
          END IF
          IF((iter > 10) .AND. (iter < 50)) THEN
            exposureGMap(TempX, TempY)  = exposureGMap(TempX, TempY) + 1
            exposureGMap(TempX, TempYm) = exposureGMap(TempX, TempYm) + 1
          END IF
          IF((iter > 50) .AND. (iter < 100)) THEN
            exposureBMap(TempX, TempY)  = exposureBMap(TempX, TempY) + 1
            exposureBMap(TempX, TempYm) = exposureBMap(TempX, TempYm) + 1
          ENDIF
        END IF
      !END IF !(cabs(z)<4)
    END DO
  END IF
END DO
!$END PARALLEL
 


maxRExposure = MAXVAL(exposureRMap)
!minRExposure = MINVAL(exposureRMap)
maxGExposure = MAXVAL(exposureGMap)
!minGExposure = MINVAL(exposureGMap)
maxBExposure = MAXVAL(exposureBMap)
!minBExposure = MINVAL(exposureBMap)
!write(*,*) maxRExposure, minRExposure, maxGExposure, minGExposure, maxBExposure, minBExposure
 
!minExposure = MIN(minRExposure, minGExposure, minBExposure)
maxExposure = MAX(maxRExposure, maxGExposure, maxBExposure)
 
!exposureRMap = exposureRMap - minExposure
!exposureGMap = exposureGMap - minExposure
!exposureBMap = exposureBMap - minExposure
exposureRMap = (exposureRMap / REAL(maxRExposure))*intensityR
exposureGMap = (exposureGMap / REAL(maxGExposure))*intensityG
exposureBMap = (exposureBMap / REAL(maxBExposure))*intensityB
 
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
  INTEGER :: n
  COMPLEX :: z
  LOGICAL :: notInMSet
  z = CMPLX(0,0) 
  n = 0
  IF(((ABS(c - CMPLX(-1,0) )) < 0.25) .OR. (( ABS( 1.0 - SQRT(1-(4*c)) ) < 1.0 ) ) ) THEN
    notInMset = .FALSE.
  ELSE
    DO WHILE (ABS(z) < escapeOrbit .AND. (n < n_max))
      z = z**zpower + c
      n = n + 1
    END DO
 
    IF (n >= n_max) THEN
      notInMset = .FALSE.
    ELSE
      notInMset = .TRUE.
    END IF
  END IF
END FUNCTION notInMset
 
 
END
 

