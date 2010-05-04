PROGRAM MANDEL
USE omp_lib
IMPLICIT NONE
 
CHARACTER ( len = 255 ), PARAMETER :: filename = 'buddhabrot.ppm'
INTEGER, PARAMETER :: file_out_unit = 10
 
 
INTEGER, PARAMETER :: n_max=1000000
INTEGER, PARAMETER :: grid_resolution = 1024
INTEGER, PARAMETER :: grid_center = grid_resolution/2
INTEGER, PARAMETER :: zpower = 2
INTEGER, PARAMETER :: miniter = 10000
INTEGER*8, PARAMETER :: batchSize = 1000000
REAL, PARAMETER :: escapeOrbit = 4
REAL, PARAMETER :: xmin = -1.0, xmax = 2.0, ymin = -1.3, ymax =1.3
 
REAL, PARAMETER :: intensity = 255.
 
INTEGER :: exposureMap(grid_resolution, grid_resolution)
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
exposureMap = 0
 
!$OMP PARALLEL DO DEFAULT(PRIVATE),SHARED(exposureMap)
DO i=1, batchSize
  CALL RANDOM_NUMBER(x)
  CALL RANDOM_NUMBER(y)
  z = CMPLX(0,0)
  c = CMPLX((x*2.5 - 2) ,y*1.3 ) !choose a random point on complex plane
  IF (notInMSet(c, n_max)) THEN !if it espace out of the mandelbrot set
    DO iter=1, n_max !iterate and plot orbit
      z = z**zpower + c !mandelbrot formula : Z = Z²+C
      IF(iter .GE. miniter) THEN
        TempX = INT(grid_resolution * (REAL(z) + xmax) / (xmax - xmin))
        TempY = INT(grid_resolution * (AIMAG(z) + ymax) / (ymax - ymin))
        TempYm = INT(grid_center - (TempY - grid_resolution/2))
        IF((TempX > 0) .AND. (TempX < grid_resolution) .AND. (TempY > 0) .AND. (TempY < grid_resolution)) THEN
            exposureMap(TempX, TempY)  = exposureMap(TempX, TempY) + 1
            exposureMap(TempX, TempYm) = exposureMap(TempX, TempYm) + 1
        END IF
      END IF 
    END DO
  END IF
END DO
!$END PARALLEL
 
!minExposure = MINVAL(exposureMap)
maxExposure = MAXVAL(exposureMap)
!write(*,*) maxExposure, minExposure
 
exposureMap = SQRT(exposureMap / REAL(maxExposure))*intensity
 
open ( unit = file_out_unit, file = filename, status = 'replace', &
       form = 'formatted', access = 'sequential')
write ( file_out_unit, '(a2)' ) 'P3'
write ( file_out_unit, '(i5,2x,i5)' ) grid_resolution, grid_resolution
write ( file_out_unit, '(i5)' ) INT(intensity)
 
do ppm_i = 1, grid_resolution
do ppm_jlo = 1, grid_resolution, 4
    ppm_jhi = min ( ppm_jlo + 3, grid_resolution )
    write ( file_out_unit, '(12i5)' ) &
    ( exposureMap(ppm_i,ppm_j), exposureMap(ppm_i,ppm_j),exposureMap(ppm_i,ppm_j), ppm_j = ppm_jlo,ppm_jhi )
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
  IF(((ABS(c - CMPLX(-1,0) )) < 0.25) .OR. ((ABS( 1.0 - SQRT(1-(4*c)) ))  < 1.0 ) ) THEN
    notInMset = .FALSE.
  ELSE
    DO WHILE (ABS(z) < escapeOrbit .AND. (n < n_max))
      z = z**zpower + c
      n = n + 1
    END DO
 
    IF (n >= n_max) THEN
      notInMset = .FALSE.
    ELSE IF (n < miniter) THEN
      notInMset = .FALSE.  !Dirty trick to ignore point that exit before miniter
    ELSE
      notInMset = .TRUE.
    END IF
  END IF
END FUNCTION notInMset
 
 
END
 

