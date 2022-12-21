PROGRAM nemo_angle

 IMPLICIT NONE

include "netcdf.inc"

 INTEGER, PARAMETER :: stdout = 6
 INTEGER, PARAMETER :: stderr = 0

 INTEGER, PARAMETER :: sp  = SELECTED_REAL_KIND(6,37)
 INTEGER, PARAMETER :: dp  = SELECTED_REAL_KIND(12,307)
 INTEGER, PARAMETER :: wp  = dp

 INTEGER, PARAMETER :: i4  = SELECTED_INT_KIND(9)
 INTEGER, PARAMETER :: i8  = SELECTED_INT_KIND(18)

 ! Constants
 REAL(wp), PARAMETER :: zero    = 0.0_wp
 REAL(wp), PARAMETER :: one     = 1.0_wp
 REAL(wp), PARAMETER :: two     = 2.0_wp
 REAL(wp), PARAMETER :: four    = 4.0_wp
 !
 REAL(wp), PARAMETER :: half    = 0.5_wp
 REAL(wp), PARAMETER :: quarter = 0.25_wp
 !
 REAL(wp), PARAMETER :: dcircle     = 360.0_wp
 REAL(wp), PARAMETER :: dhalfcircle = 180.0_wp
 !
 REAL(wp), PARAMETER :: rpi     = 3.14159265358979323846_wp
 REAL(wp), PARAMETER :: deg2rad = rpi / dhalfcircle  !: conversion from degree into radian
 REAL(wp), PARAMETER :: rad2deg = dhalfcircle / rpi  !: conversion from radian into degree
 REAL(wp), PARAMETER :: miss_value = -1.0e30_wp
 !
 REAL(wp), PARAMETER :: eps1 = 1.e-8_wp
 REAL(wp), PARAMETER :: eps2 = 2.e-15_wp

 CHARACTER(LEN=16), PARAMETER :: conf = '1_BIZoo'

 INTEGER, PARAMETER :: ngrid = 4
 INTEGER, PARAMETER :: nvar  = 2
 CHARACTER(LEN=1), PARAMETER  :: grids(ngrid) = (/'T', 'U', 'V', 'F'/)
 CHARACTER(LEN=3), PARAMETER  :: vars(nvar)   = (/'sin', 'cos'/)

 REAL(wp), ALLOCATABLE, DIMENSION(:,:) ::   &
      glamt, gphit,   &  ! T point coordinates
      glamu, gphiu,   &  ! U point coordinates
      glamv, gphiv,   &  ! V point coordinates
      glamf, gphif       ! F point coordinates

 REAL(wp), ALLOCATABLE, DIMENSION(:,:) ::   &
      gsint, gcost,   &  ! cos/sin between model grid lines and NP direction at T point
      gsinu, gcosu,   &  ! cos/sin between model grid lines and NP direction at U point
      gsinv, gcosv,   &  ! cos/sin between model grid lines and NP direction at V point
      gsinf, gcosf       ! cos/sin between model grid lines and NP direction at F point

 REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: work

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !
 ! DOMAIN CONFIGURATION
 !
 INTEGER, PARAMETER :: jpiglo = 2567
 INTEGER, PARAMETER :: jpjglo = 2693

 INTEGER, PARAMETER :: nperio = 4
 INTEGER, PARAMETER :: npolj  = nperio
 !
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 INTEGER, PARAMETER :: jpi = jpiglo
 INTEGER, PARAMETER :: jpj = jpjglo
 INTEGER, PARAMETER :: jpim1 = jpi-1
 INTEGER, PARAMETER :: jpjm1 = jpj-1
 INTEGER, PARAMETER :: jpni  = 1
 INTEGER, PARAMETER :: nlci  = jpi
 INTEGER, PARAMETER :: nlcj  = jpj

 INTEGER :: ierr   = 0
 INTEGER :: iostat = 0
 INTEGER :: i,j,n,m
 INTEGER :: nv, ng

 ! NetCDF variables
 CHARACTER(LEN=NF_MAX_NAME) :: input_filename  = TRIM(conf)//'_coordinates.nc'
! CHARACTER(LEN=NF_MAX_NAME) :: input_filename = 'eNEAT36_TWIN_coordinates.nc'
 CHARACTER(LEN=NF_MAX_NAME) :: output_filename = TRIM(conf)//'_angle.nc'
 !
 CHARACTER(LEN=NF_MAX_NAME) :: str1 = ''
 CHARACTER(LEN=NF_MAX_NAME) :: str2 = ''
 !
 INTEGER :: ncid_in, varid_in(nvar*ngrid)
 INTEGER :: ncid_out, varid_out(4*nvar*ngrid)
 INTEGER :: x_id, y_id
 INTEGER :: lon_dim_id, lat_dim_id
 INTEGER :: lon_id, lat_id



!  WRITE(stdout,FMT='(2E30.16)') deg2rad*rad2deg, rad2deg*deg2rad
!  WRITE(stdout,FMT='(3E30.16)') rpi, rpi-two*ASIN(one), rpi-ACOS(-one)

   ! Open output file
   iostat = nf_open(TRIM(input_filename), NF_NOWRITE, ncid_in)
   if (iostat /= NF_NOERR) call nf_handle_err(iostat, __LINE__)

   iostat = nf_inq_dimid(ncid_in, 'x', x_id)
   if (iostat /= NF_NOERR) call nf_handle_err(iostat, __LINE__)

   iostat = nf_inq_dimlen(ncid_in, x_id, n)
   if (iostat /= NF_NOERR) call nf_handle_err(iostat, __LINE__)

   IF (n /= jpiglo) THEN
     WRITE(stderr,*) 'ERROR: Wrong grid size! Expected x=', jpiglo, ' found x=',n
     STOP 'X'
   END IF

   iostat = nf_inq_dimid(ncid_in, 'y', y_id)
   if (iostat /= NF_NOERR) call nf_handle_err(iostat, __LINE__)

   iostat = nf_inq_dimlen(ncid_in, y_id, n)
   if (iostat /= NF_NOERR) call nf_handle_err(iostat, __LINE__)

   IF (n /= jpjglo) THEN
     WRITE(stderr,*) 'ERROR: Wrong grid size! Expected y=', jpjglo, ' found y=',n
     STOP 'Y'
   END IF


   ALLOCATE( glamt(jpi,jpj), gphit(jpi,jpj),   & 
      &      glamu(jpi,jpj), gphiu(jpi,jpj),   & 
      &      glamv(jpi,jpj), gphiv(jpi,jpj),   &  
      &      glamf(jpi,jpj), gphif(jpi,jpj), STAT=ierr )
   IF( ierr /= 0 )   stop 'angle: unable to allocate coordinate arrays'


   n=1

   iostat = nf_inq_varid(ncid_in, 'glamt', varid_in(n))
   if (iostat /= NF_NOERR) call nf_handle_err(iostat, __LINE__)

   iostat = nf_get_var_double(ncid_in, varid_in(n), glamt)
   if (iostat /= NF_NOERR) call nf_handle_err(iostat, __LINE__)

   n=n+1

   iostat = nf_inq_varid(ncid_in, 'glamu', varid_in(n))
   if (iostat /= NF_NOERR) call nf_handle_err(iostat, __LINE__)

   iostat = nf_get_var_double(ncid_in, varid_in(n), glamu)
   if (iostat /= NF_NOERR) call nf_handle_err(iostat, __LINE__)

   n=n+1

   iostat = nf_inq_varid(ncid_in, 'glamv', varid_in(n))
   if (iostat /= NF_NOERR) call nf_handle_err(iostat, __LINE__)

   iostat = nf_get_var_double(ncid_in, varid_in(n), glamv)
   if (iostat /= NF_NOERR) call nf_handle_err(iostat, __LINE__)

   n=n+1

   iostat = nf_inq_varid(ncid_in, 'glamf', varid_in(n))
   if (iostat /= NF_NOERR) call nf_handle_err(iostat, __LINE__)

   iostat = nf_get_var_double(ncid_in, varid_in(n), glamf)
   if (iostat /= NF_NOERR) call nf_handle_err(iostat, __LINE__)

   n=n+1

   iostat = nf_inq_varid(ncid_in, 'gphit', varid_in(n))
   if (iostat /= NF_NOERR) call nf_handle_err(iostat, __LINE__)

   iostat = nf_get_var_double(ncid_in, varid_in(n), gphit)
   if (iostat /= NF_NOERR) call nf_handle_err(iostat, __LINE__)

   n=n+1

   iostat = nf_inq_varid(ncid_in, 'gphiu', varid_in(n))
   if (iostat /= NF_NOERR) call nf_handle_err(iostat, __LINE__)

   iostat = nf_get_var_double(ncid_in, varid_in(n), gphiu)
   if (iostat /= NF_NOERR) call nf_handle_err(iostat, __LINE__)

   n=n+1

   iostat = nf_inq_varid(ncid_in, 'gphiv', varid_in(n))
   if (iostat /= NF_NOERR) call nf_handle_err(iostat, __LINE__)

   iostat = nf_get_var_double(ncid_in, varid_in(n), gphiv)
   if (iostat /= NF_NOERR) call nf_handle_err(iostat, __LINE__)

   n=n+1

   iostat = nf_inq_varid(ncid_in, 'gphif', varid_in(n))
   if (iostat /= NF_NOERR) call nf_handle_err(iostat, __LINE__)

   iostat = nf_get_var_double(ncid_in, varid_in(n), gphif)
   if (iostat /= NF_NOERR) call nf_handle_err(iostat, __LINE__)

   iostat = nf_close(ncid_in)
   if (iostat /= NF_NOERR) call nf_handle_err(iostat, __LINE__)


   ALLOCATE( gsint(jpi,jpj), gcost(jpi,jpj),   & 
      &      gsinu(jpi,jpj), gcosu(jpi,jpj),   & 
      &      gsinv(jpi,jpj), gcosv(jpi,jpj),   &  
      &      gsinf(jpi,jpj), gcosf(jpi,jpj),   &
      &      work (jpi,jpj), STAT=ierr )
   IF( ierr /= 0 )   stop 'angle: unable to allocate angle arrays'

   CALL angle

   ! Open output file
   iostat = nf_create(TRIM(output_filename), IOR(NF_CLASSIC_MODEL, NF_NETCDF4), ncid_out)
   if (iostat /= NF_NOERR) call nf_handle_err(iostat, __LINE__)

   ! Define global attributes
   str1='CF-1.7'
   iostat = nf_put_att_text(ncid_out, NF_GLOBAL, 'Conventions', LEN_TRIM(str1), TRIM(str1))
   if (iostat /= NF_NOERR) call nf_handle_err(iostat, __LINE__)

   str1=TRIM(conf)//' grid angles related variables'
   iostat = nf_put_att_text(ncid_out, NF_GLOBAL, 'title', LEN_TRIM(str1), TRIM(str1))
   if (iostat /= NF_NOERR) call nf_handle_err(iostat, __LINE__)

   str1='CMCC'
   iostat = nf_put_att_text(ncid_out, NF_GLOBAL, 'institution', LEN_TRIM(str1), TRIM(str1))
   if (iostat /= NF_NOERR) call nf_handle_err(iostat, __LINE__)

   str1='NEMO '//TRIM(conf)//' coordinates.nc'
   iostat = nf_put_att_text(ncid_out, NF_GLOBAL, 'source', LEN_TRIM(str1), TRIM(str1))
   if (iostat /= NF_NOERR) call nf_handle_err(iostat, __LINE__)

   str1='Thu Mar 14 15:46:27 UTC 2019: nemo_angle.F90'
   iostat = nf_put_att_text(ncid_out, NF_GLOBAL, 'history', LEN_TRIM(str1), TRIM(str1))
   if (iostat /= NF_NOERR) call nf_handle_err(iostat, __LINE__)

   str1='Based on NEMO code: geo2ocean.F90, lbclnk.F90, lbcnfd.F90'
   iostat = nf_put_att_text(ncid_out, NF_GLOBAL, 'comment', LEN_TRIM(str1), TRIM(str1))
   if (iostat /= NF_NOERR) call nf_handle_err(iostat, __LINE__)


   ! Define dimensions
!   iostat = nf_def_dim(ncid_out, 'time', NF_UNLIMITED, time_dim_id)
!   if (iostat /= NF_NOERR) call nf_handle_err(iostat, __LINE__)

   iostat = nf_def_dim(ncid_out, 'x', jpiglo, lon_dim_id)
   if (iostat /= NF_NOERR) call nf_handle_err(iostat, __LINE__)

   iostat = nf_def_dim(ncid_out, 'y', jpjglo, lat_dim_id)
   if (iostat /= NF_NOERR) call nf_handle_err(iostat, __LINE__)

   ! Define time variable
!   iostat = nf_def_var(ncid_out, 'time', NF_DOUBLE, 1, (/ time_dim_id /), time_id)
!   if (iostat /= NF_NOERR) call nf_handle_err(iostat, __LINE__)

   ! Define attributes
!   str1='time'
!   iostat = nf_put_att_text(ncid_out, time_id, 'standard_name', LEN_TRIM(str1), TRIM(str1))
!   if (iostat /= NF_NOERR) call nf_handle_err(iostat, __LINE__)

!   str1='months since 1978-01-16 12:00:00'
!   iostat = nf_put_att_text(ncid_out, time_id, 'units', LEN_TRIM(str1), TRIM(str1))
!   if (iostat /= NF_NOERR) call nf_handle_err(iostat, __LINE__)

!   str1='standard'
!   iostat = nf_put_att_text(ncid_out, time_id, 'calendar', LEN_TRIM(str1), TRIM(str1))
!   if (iostat /= NF_NOERR) call nf_handle_err(iostat, __LINE__)

!   str1='time'
!   iostat = nf_put_att_text(ncid_out, time_id, 'long_name', LEN_TRIM(str1), TRIM(str1))
!   if (iostat /= NF_NOERR) call nf_handle_err(iostat, __LINE__)

!   str1='T'
!   iostat = nf_put_att_text(ncid_out, time_id, 'axis', LEN_TRIM(str1), TRIM(str1))
!   if (iostat /= NF_NOERR) call nf_handle_err(iostat, __LINE__)

   ! Define longitude coordinate variable
   iostat = nf_def_var(ncid_out, 'lon', NF_DOUBLE, 2, (/ lon_dim_id, lat_dim_id /), lon_id)
   if (iostat /= NF_NOERR) call nf_handle_err(iostat, __LINE__)

   iostat = nf_def_var_deflate(ncid_out, lon_id, 1, 1, 4)
   if (iostat /= NF_NOERR) call nf_handle_err(iostat, __LINE__)

   ! Define attributes
   str1='longitude'
   iostat = nf_put_att_text(ncid_out, lon_id, 'standard_name', LEN_TRIM(str1), TRIM(str1))
   if (iostat /= NF_NOERR) call nf_handle_err(iostat, __LINE__)

   str1='degrees_east'
   iostat = nf_put_att_text(ncid_out, lon_id, 'units', LEN_TRIM(str1), TRIM(str1))
   if (iostat /= NF_NOERR) call nf_handle_err(iostat, __LINE__)

   str1='longitude'
   iostat = nf_put_att_text(ncid_out, lon_id, 'long_name', LEN_TRIM(str1), TRIM(str1))
   if (iostat /= NF_NOERR) call nf_handle_err(iostat, __LINE__)

   str1='X'
   iostat = nf_put_att_text(ncid_out, lon_id, 'axis', LEN_TRIM(str1), TRIM(str1))
   if (iostat /= NF_NOERR) call nf_handle_err(iostat, __LINE__)

   ! Define latitude coordinate variable
   iostat = nf_def_var(ncid_out, 'lat', NF_DOUBLE, 2, (/ lon_dim_id, lat_dim_id /), lat_id)
   if (iostat /= NF_NOERR) call nf_handle_err(iostat, __LINE__)

   iostat = nf_def_var_deflate(ncid_out, lat_id, 1, 1, 4)
   if (iostat /= NF_NOERR) call nf_handle_err(iostat, __LINE__)

   ! Define attributes
   str1='latitude'
   iostat = nf_put_att_text(ncid_out, lat_id, 'standard_name', LEN_TRIM(str1), TRIM(str1))
   if (iostat /= NF_NOERR) call nf_handle_err(iostat, __LINE__)

   str1='degrees_north'
   iostat = nf_put_att_text(ncid_out, lat_id, 'units', LEN_TRIM(str1), TRIM(str1))
   if (iostat /= NF_NOERR) call nf_handle_err(iostat, __LINE__)

   str1='latitude'
   iostat = nf_put_att_text(ncid_out, lat_id, 'long_name', LEN_TRIM(str1), TRIM(str1))
   if (iostat /= NF_NOERR) call nf_handle_err(iostat, __LINE__)

   str1='Y'
   iostat = nf_put_att_text(ncid_out, lat_id, 'axis', LEN_TRIM(str1), TRIM(str1))
   if (iostat /= NF_NOERR) call nf_handle_err(iostat, __LINE__)


   n=1

 if (1==0) then
   DO nv=1, nvar

     IF (TRIM(vars(nv)) == 'sin') THEN
       str2 = 'Sine'
     ELSEIF (TRIM(vars(nv)) == 'cos') THEN
       str2 = 'Cosine'
     END IF

     DO ng=1, ngrid

       ! Define variable
       str1='g'//TRIM(vars(nv))//lower(TRIM(grids(ng)))
       iostat = nf_def_var(ncid_out, TRIM(str1), NF_DOUBLE, 2, (/ lon_dim_id, lat_dim_id /), varid_out(n))
       if (iostat /= NF_NOERR) call nf_handle_err(iostat, __LINE__)

       ! Define fill value
       iostat = nf_def_var_fill(ncid_out, varid_out(n), 0, miss_value)
       if (iostat /= NF_NOERR) call nf_handle_err(iostat, __LINE__)

       iostat = nf_def_var_deflate(ncid_out, varid_out(n), 1, 1, 4)
       if (iostat /= NF_NOERR) call nf_handle_err(iostat, __LINE__)

       ! Define attributes
!       str1=''
!       iostat = nf_put_att_text(ncid_out, varid_out(n), 'standard_name', LEN_TRIM(str1), TRIM(str1))
!       if (iostat /= NF_NOERR) call nf_handle_err(iostat, __LINE__)

       str1='1'
       iostat = nf_put_att_text(ncid_out, varid_out(n), 'units', LEN_TRIM(str1), TRIM(str1))
       if (iostat /= NF_NOERR) call nf_handle_err(iostat, __LINE__)

       str1=TRIM(str2)//' of the angle between model grid lines and North Pole direction at '//TRIM(grids(ng))//' points'
       iostat = nf_put_att_text(ncid_out, varid_out(n), 'long_name', LEN_TRIM(str1), TRIM(str1))
       if (iostat /= NF_NOERR) call nf_handle_err(iostat, __LINE__)

       iostat = nf_put_att_double(ncid_out, varid_out(n), 'missing_value', NF_DOUBLE, 1, miss_value)
       if (iostat /= NF_NOERR) call nf_handle_err(iostat, __LINE__)

       str1='lon lat'
       iostat = nf_put_att_text(ncid_out, varid_out(n), 'coordinates', LEN_TRIM(str1), TRIM(str1))
       if (iostat /= NF_NOERR) call nf_handle_err(iostat, __LINE__)

       n=n+1

     END DO
   END DO
 end if

!   DO nv=1, nvar

!     IF (TRIM(vars(nv)) == 'sin') THEN
       str2 = 'sine'
!     ELSEIF (TRIM(vars(nv)) == 'cos') THEN
!       str2 = 'sosine'
!     END IF

     DO ng=1, ngrid

       ! Define variable
       str1=lower(TRIM(grids(ng)))//'angle'
       iostat = nf_def_var(ncid_out, TRIM(str1), NF_DOUBLE, 2, (/ lon_dim_id, lat_dim_id /), varid_out(n))
       if (iostat /= NF_NOERR) call nf_handle_err(iostat, __LINE__)

       ! Define fill value
       iostat = nf_def_var_fill(ncid_out, varid_out(n), 0, miss_value)
       if (iostat /= NF_NOERR) call nf_handle_err(iostat, __LINE__)

       iostat = nf_def_var_deflate(ncid_out, varid_out(n), 1, 1, 4)
       if (iostat /= NF_NOERR) call nf_handle_err(iostat, __LINE__)

       ! Define attributes
!       str1=''
!       iostat = nf_put_att_text(ncid_out, varid_out(n), 'standard_name', LEN_TRIM(str1), TRIM(str1))
!       if (iostat /= NF_NOERR) call nf_handle_err(iostat, __LINE__)

       str1='degrees'
       iostat = nf_put_att_text(ncid_out, varid_out(n), 'units', LEN_TRIM(str1), TRIM(str1))
       if (iostat /= NF_NOERR) call nf_handle_err(iostat, __LINE__)

       str1='Angle between model curvilinear grid and geographic grid at '//TRIM(grids(ng))//' points'!//' (from '//TRIM(str2)//')'
       iostat = nf_put_att_text(ncid_out, varid_out(n), 'long_name', LEN_TRIM(str1), TRIM(str1))
       if (iostat /= NF_NOERR) call nf_handle_err(iostat, __LINE__)

       iostat = nf_put_att_double(ncid_out, varid_out(n), 'missing_value', NF_DOUBLE, 1, miss_value)
       if (iostat /= NF_NOERR) call nf_handle_err(iostat, __LINE__)

       str1='lon lat'
       iostat = nf_put_att_text(ncid_out, varid_out(n), 'coordinates', LEN_TRIM(str1), TRIM(str1))
       if (iostat /= NF_NOERR) call nf_handle_err(iostat, __LINE__)

       n=n+1

     END DO
!   END DO

 if (1==0) then
   DO nv=1, nvar

     IF (TRIM(vars(nv)) == 'sin') THEN
       str2 = 'Sine'
     ELSEIF (TRIM(vars(nv)) == 'cos') THEN
       str2 = 'Cosine'
     END IF

     DO ng=1, ngrid

       ! Define variable
       str1='g'//TRIM(vars(nv))//lower(TRIM(grids(ng)))//'2'
       iostat = nf_def_var(ncid_out, TRIM(str1), NF_DOUBLE, 2, (/ lon_dim_id, lat_dim_id /), varid_out(n))
       if (iostat /= NF_NOERR) call nf_handle_err(iostat, __LINE__)

       ! Define fill value
       iostat = nf_def_var_fill(ncid_out, varid_out(n), 0, miss_value)
       if (iostat /= NF_NOERR) call nf_handle_err(iostat, __LINE__)

       iostat = nf_def_var_deflate(ncid_out, varid_out(n), 1, 1, 4)
       if (iostat /= NF_NOERR) call nf_handle_err(iostat, __LINE__)

       ! Define attributes
!       str1=''
!       iostat = nf_put_att_text(ncid_out, varid_out(n), 'standard_name', LEN_TRIM(str1), TRIM(str1))
!       if (iostat /= NF_NOERR) call nf_handle_err(iostat, __LINE__)

       str1='1'
       iostat = nf_put_att_text(ncid_out, varid_out(n), 'units', LEN_TRIM(str1), TRIM(str1))
       if (iostat /= NF_NOERR) call nf_handle_err(iostat, __LINE__)

       str1=TRIM(str2)//' of the angle between model grid lines and North Pole direction at '//TRIM(grids(ng))//' points'
       iostat = nf_put_att_text(ncid_out, varid_out(n), 'long_name', LEN_TRIM(str1), TRIM(str1))
       if (iostat /= NF_NOERR) call nf_handle_err(iostat, __LINE__)

       iostat = nf_put_att_double(ncid_out, varid_out(n), 'missing_value', NF_DOUBLE, 1, miss_value)
       if (iostat /= NF_NOERR) call nf_handle_err(iostat, __LINE__)

       str1='lon lat'
       iostat = nf_put_att_text(ncid_out, varid_out(n), 'coordinates', LEN_TRIM(str1), TRIM(str1))
       if (iostat /= NF_NOERR) call nf_handle_err(iostat, __LINE__)

       n=n+1

     END DO
   END DO

   DO nv=1, nvar

     IF (TRIM(vars(nv)) == 'sin') THEN
       str2 = 'Sine'
     ELSEIF (TRIM(vars(nv)) == 'cos') THEN
       str2 = 'Cosine'
     END IF

     DO ng=1, ngrid

       ! Define variable
       str1='g'//TRIM(vars(nv))//lower(TRIM(grids(ng)))//'3'
       iostat = nf_def_var(ncid_out, TRIM(str1), NF_DOUBLE, 2, (/ lon_dim_id, lat_dim_id /), varid_out(n))
       if (iostat /= NF_NOERR) call nf_handle_err(iostat, __LINE__)

       ! Define fill value
       iostat = nf_def_var_fill(ncid_out, varid_out(n), 0, miss_value)
       if (iostat /= NF_NOERR) call nf_handle_err(iostat, __LINE__)

       iostat = nf_def_var_deflate(ncid_out, varid_out(n), 1, 1, 4)
       if (iostat /= NF_NOERR) call nf_handle_err(iostat, __LINE__)

       ! Define attributes
!       str1=''
!       iostat = nf_put_att_text(ncid_out, varid_out(n), 'standard_name', LEN_TRIM(str1), TRIM(str1))
!       if (iostat /= NF_NOERR) call nf_handle_err(iostat, __LINE__)

       str1='1'
       iostat = nf_put_att_text(ncid_out, varid_out(n), 'units', LEN_TRIM(str1), TRIM(str1))
       if (iostat /= NF_NOERR) call nf_handle_err(iostat, __LINE__)

       str1=TRIM(str2)//' of the angle between model grid lines and North Pole direction at '//TRIM(grids(ng))//' points'
       iostat = nf_put_att_text(ncid_out, varid_out(n), 'long_name', LEN_TRIM(str1), TRIM(str1))
       if (iostat /= NF_NOERR) call nf_handle_err(iostat, __LINE__)

       iostat = nf_put_att_double(ncid_out, varid_out(n), 'missing_value', NF_DOUBLE, 1, miss_value)
       if (iostat /= NF_NOERR) call nf_handle_err(iostat, __LINE__)

       str1='lon lat'
       iostat = nf_put_att_text(ncid_out, varid_out(n), 'coordinates', LEN_TRIM(str1), TRIM(str1))
       if (iostat /= NF_NOERR) call nf_handle_err(iostat, __LINE__)

       n=n+1

     END DO
   END DO
 end if

  ! End definition part
  iostat = nf_enddef(ncid_out)
  if (iostat /= NF_NOERR) call nf_handle_err(iostat, __LINE__)

  ! Fill coordinate variables using grid T
  iostat = nf_put_var_double(ncid_out, lon_id, glamt)
  if (iostat /= NF_NOERR) call nf_handle_err(iostat, __LINE__)

  iostat = nf_put_var_double(ncid_out, lat_id, gphit)
  if (iostat /= NF_NOERR) call nf_handle_err(iostat, __LINE__)

  ! Write sine/cosine of the angle

   n=1

 if (1==0) then
   iostat = nf_put_var_double(ncid_out, varid_out(n), gsint)
   if (iostat /= NF_NOERR) call nf_handle_err(iostat, __LINE__)
   n=n+1

   iostat = nf_put_var_double(ncid_out, varid_out(n), gsinu)
   if (iostat /= NF_NOERR) call nf_handle_err(iostat, __LINE__)
   n=n+1

   iostat = nf_put_var_double(ncid_out, varid_out(n), gsinv)
   if (iostat /= NF_NOERR) call nf_handle_err(iostat, __LINE__)
   n=n+1

   iostat = nf_put_var_double(ncid_out, varid_out(n), gsinf)
   if (iostat /= NF_NOERR) call nf_handle_err(iostat, __LINE__)
   n=n+1

   iostat = nf_put_var_double(ncid_out, varid_out(n), gcost)
   if (iostat /= NF_NOERR) call nf_handle_err(iostat, __LINE__)
   n=n+1

   iostat = nf_put_var_double(ncid_out, varid_out(n), gcosu)
   if (iostat /= NF_NOERR) call nf_handle_err(iostat, __LINE__)
   n=n+1

   iostat = nf_put_var_double(ncid_out, varid_out(n), gcosv)
   if (iostat /= NF_NOERR) call nf_handle_err(iostat, __LINE__)
   n=n+1

   iostat = nf_put_var_double(ncid_out, varid_out(n), gcosf)
   if (iostat /= NF_NOERR) call nf_handle_err(iostat, __LINE__)
   n=n+1
 end if

  ! Write the angle computed from sine/cosine

   work(:,:) = rad2deg*ASIN(gsint(:,:))!; CALL lbc_lnk( work, 'T', -one )
   iostat = nf_put_var_double(ncid_out, varid_out(n), work)
   if (iostat /= NF_NOERR) call nf_handle_err(iostat, __LINE__)
   n=n+1

   work(:,:) = rad2deg*ASIN(gsinu(:,:))!; CALL lbc_lnk( work, 'U', -one )
   iostat = nf_put_var_double(ncid_out, varid_out(n), work)
   if (iostat /= NF_NOERR) call nf_handle_err(iostat, __LINE__)
   n=n+1

   work(:,:) = rad2deg*ASIN(gsinv(:,:))!; CALL lbc_lnk( work, 'V', -one )
   iostat = nf_put_var_double(ncid_out, varid_out(n), work)
   if (iostat /= NF_NOERR) call nf_handle_err(iostat, __LINE__)
   n=n+1

   work(:,:) = rad2deg*ASIN(gsinf(:,:))!; CALL lbc_lnk( work, 'F', -one )
   iostat = nf_put_var_double(ncid_out, varid_out(n), work)
   if (iostat /= NF_NOERR) call nf_handle_err(iostat, __LINE__)
   n=n+1

!   work(:,:) = rad2deg*ACOS(gcost(:,:)); CALL lbc_lnk( work, 'T', -one )
!   iostat = nf_put_var_double(ncid_out, varid_out(n), work)
!   if (iostat /= NF_NOERR) call nf_handle_err(iostat, __LINE__)
!   n=n+1

!   work(:,:) = rad2deg*ACOS(gcosu(:,:)); CALL lbc_lnk( work, 'U', -one )
!   iostat = nf_put_var_double(ncid_out, varid_out(n), work)
!   if (iostat /= NF_NOERR) call nf_handle_err(iostat, __LINE__)
!   n=n+1

!   work(:,:) = rad2deg*ACOS(gcosv(:,:)); CALL lbc_lnk( work, 'V', -one )
!   iostat = nf_put_var_double(ncid_out, varid_out(n), work)
!   if (iostat /= NF_NOERR) call nf_handle_err(iostat, __LINE__)
!   n=n+1

!   work(:,:) = rad2deg*ACOS(gcosf(:,:)); CALL lbc_lnk( work, 'F', -one )
!   iostat = nf_put_var_double(ncid_out, varid_out(n), work)
!   if (iostat /= NF_NOERR) call nf_handle_err(iostat, __LINE__)
!   n=n+1

  ! Write re-computed sine/cosine of the angle

 if (1==0) then
   work(:,:) = SIN(deg2rad*(rad2deg*ASIN(gsint(:,:))))
   iostat = nf_put_var_double(ncid_out, varid_out(n), work)
   if (iostat /= NF_NOERR) call nf_handle_err(iostat, __LINE__)
   n=n+1

   work(:,:) = SIN(deg2rad*(rad2deg*ASIN(gsinu(:,:))))
   iostat = nf_put_var_double(ncid_out, varid_out(n), work)
   if (iostat /= NF_NOERR) call nf_handle_err(iostat, __LINE__)
   n=n+1

   work(:,:) = SIN(deg2rad*(rad2deg*ASIN(gsinv(:,:))))
   iostat = nf_put_var_double(ncid_out, varid_out(n), work)
   if (iostat /= NF_NOERR) call nf_handle_err(iostat, __LINE__)
   n=n+1

   work(:,:) = SIN(deg2rad*(rad2deg*ASIN(gsinf(:,:))))
   iostat = nf_put_var_double(ncid_out, varid_out(n), work)
   if (iostat /= NF_NOERR) call nf_handle_err(iostat, __LINE__)
   n=n+1

   work(:,:) = COS(deg2rad*(rad2deg*ACOS(gcost(:,:))))
   iostat = nf_put_var_double(ncid_out, varid_out(n), work)
   if (iostat /= NF_NOERR) call nf_handle_err(iostat, __LINE__)
   n=n+1

   work(:,:) = COS(deg2rad*(rad2deg*ACOS(gcosu(:,:))))
   iostat = nf_put_var_double(ncid_out, varid_out(n), work)
   if (iostat /= NF_NOERR) call nf_handle_err(iostat, __LINE__)
   n=n+1

   work(:,:) = COS(deg2rad*(rad2deg*ACOS(gcosv(:,:))))
   iostat = nf_put_var_double(ncid_out, varid_out(n), work)
   if (iostat /= NF_NOERR) call nf_handle_err(iostat, __LINE__)
   n=n+1

   work(:,:) = COS(deg2rad*(rad2deg*ACOS(gcosf(:,:))))
   iostat = nf_put_var_double(ncid_out, varid_out(n), work)
   if (iostat /= NF_NOERR) call nf_handle_err(iostat, __LINE__)
   n=n+1

  ! Write re-computed sine/cosine of the angle

   work(:,:) = SIN(deg2rad*(rad2deg*ACOS(gcost(:,:))))
   iostat = nf_put_var_double(ncid_out, varid_out(n), work)
   if (iostat /= NF_NOERR) call nf_handle_err(iostat, __LINE__)
   n=n+1

   work(:,:) = SIN(deg2rad*(rad2deg*ACOS(gcosu(:,:))))
   iostat = nf_put_var_double(ncid_out, varid_out(n), work)
   if (iostat /= NF_NOERR) call nf_handle_err(iostat, __LINE__)
   n=n+1

   work(:,:) = SIN(deg2rad*(rad2deg*ACOS(gcosv(:,:))))
   iostat = nf_put_var_double(ncid_out, varid_out(n), work)
   if (iostat /= NF_NOERR) call nf_handle_err(iostat, __LINE__)
   n=n+1

   work(:,:) = SIN(deg2rad*(rad2deg*ACOS(gcosf(:,:))))
   iostat = nf_put_var_double(ncid_out, varid_out(n), work)
   if (iostat /= NF_NOERR) call nf_handle_err(iostat, __LINE__)
   n=n+1

   work(:,:) = COS(deg2rad*(rad2deg*ASIN(gsint(:,:))))
   iostat = nf_put_var_double(ncid_out, varid_out(n), work)
   if (iostat /= NF_NOERR) call nf_handle_err(iostat, __LINE__)
   n=n+1

   work(:,:) = COS(deg2rad*(rad2deg*ASIN(gsinu(:,:))))
   iostat = nf_put_var_double(ncid_out, varid_out(n), work)
   if (iostat /= NF_NOERR) call nf_handle_err(iostat, __LINE__)
   n=n+1

   work(:,:) = COS(deg2rad*(rad2deg*ASIN(gsinv(:,:))))
   iostat = nf_put_var_double(ncid_out, varid_out(n), work)
   if (iostat /= NF_NOERR) call nf_handle_err(iostat, __LINE__)
   n=n+1

   work(:,:) = COS(deg2rad*(rad2deg*ASIN(gsinf(:,:))))
   iostat = nf_put_var_double(ncid_out, varid_out(n), work)
   if (iostat /= NF_NOERR) call nf_handle_err(iostat, __LINE__)
   n=n+1
 end if

   iostat = nf_close(ncid_out)
   if (iostat /= NF_NOERR) call nf_handle_err(iostat, __LINE__)

















 CONTAINS

   SUBROUTINE lbc_nfd( pt2d, cd_type, psgn, pr2dj )
      !!----------------------------------------------------------------------
      !!                  ***  routine lbc_nfd  ***
      !!
      !! ** Purpose :   2D lateral boundary condition : North fold treatment
      !!       without processor exchanges. 
      !!
      !! ** Method  :   
      !!
      !! ** Action  :   pt2d with updated values along the north fold
      !!----------------------------------------------------------------------
      CHARACTER(len=1)        , INTENT(in   ) ::   cd_type   ! define the nature of ptab array grid-points
      !                                                      ! = T , U , V , F , W points
      REAL(wp)                , INTENT(in   ) ::   psgn      ! control of the sign change
      !                                                      !   = -1. , the sign is changed if north fold boundary
      !                                                      !   =  1. , the sign is kept  if north fold boundary
      REAL(wp), DIMENSION(jpi,jpj), INTENT(inout) ::   pt2d      ! 2D array on which the boundary condition is applied
      INTEGER , OPTIONAL      , INTENT(in   ) ::   pr2dj     ! number of additional halos
      !
      INTEGER  ::   ji, jl, ipr2dj
      INTEGER  ::   ijt, iju, ijpj, ijpjm1
      !!----------------------------------------------------------------------

      SELECT CASE ( jpni )
      CASE ( 1 )     ;   ijpj = nlcj      ! 1 proc only  along the i-direction
      CASE DEFAULT   ;   ijpj = 4         ! several proc along the i-direction
      END SELECT
      !
      IF( PRESENT(pr2dj) ) THEN           ! use of additional halos
         ipr2dj = pr2dj
         IF( jpni > 1 )   ijpj = ijpj + ipr2dj
      ELSE
         ipr2dj = 0 
      ENDIF
      !
      ijpjm1 = ijpj-1


      SELECT CASE ( npolj )
      !
      CASE ( 3, 4 )                       ! *  North fold  T-point pivot
         !
         SELECT CASE ( cd_type )
         !
         CASE ( 'T' , 'W' )                               ! T- , W-points
            DO jl = 0, ipr2dj
               DO ji = 2, jpiglo
                  ijt=jpiglo-ji+2
                  pt2d(ji,ijpj+jl) = psgn * pt2d(ijt,ijpj-2-jl)
               END DO
            END DO
            pt2d(1,ijpj)   = psgn * pt2d(3,ijpj-2)
            DO ji = jpiglo/2+1, jpiglo
               ijt=jpiglo-ji+2
               pt2d(ji,ijpj-1) = psgn * pt2d(ijt,ijpj-1)
            END DO
         CASE ( 'U' )                                     ! U-point
            DO jl = 0, ipr2dj
               DO ji = 1, jpiglo-1
                  iju = jpiglo-ji+1
                  pt2d(ji,ijpj+jl) = psgn * pt2d(iju,ijpj-2-jl)
               END DO
            END DO
            pt2d(   1  ,ijpj  ) = psgn * pt2d(    2   ,ijpj-2)
            pt2d(jpiglo,ijpj  ) = psgn * pt2d(jpiglo-1,ijpj-2)
            pt2d(1     ,ijpj-1) = psgn * pt2d(jpiglo  ,ijpj-1)   
            DO ji = jpiglo/2, jpiglo-1
               iju = jpiglo-ji+1
               pt2d(ji,ijpjm1) = psgn * pt2d(iju,ijpjm1)
            END DO
         CASE ( 'V' )                                     ! V-point
            DO jl = -1, ipr2dj
               DO ji = 2, jpiglo
                  ijt = jpiglo-ji+2
                  pt2d(ji,ijpj+jl) = psgn * pt2d(ijt,ijpj-3-jl)
               END DO
            END DO
            pt2d( 1 ,ijpj)   = psgn * pt2d( 3 ,ijpj-3) 
         CASE ( 'F' )                                     ! F-point
            DO jl = -1, ipr2dj
               DO ji = 1, jpiglo-1
                  iju = jpiglo-ji+1
                  pt2d(ji,ijpj+jl) = psgn * pt2d(iju,ijpj-3-jl)
               END DO
            END DO
            pt2d(   1  ,ijpj)   = psgn * pt2d(    2   ,ijpj-3)
            pt2d(jpiglo,ijpj)   = psgn * pt2d(jpiglo-1,ijpj-3)
            pt2d(jpiglo,ijpj-1) = psgn * pt2d(jpiglo-1,ijpj-2)      
            pt2d(   1  ,ijpj-1) = psgn * pt2d(    2   ,ijpj-2)      
         CASE ( 'I' )                                     ! ice U-V point (I-point)
            DO jl = 0, ipr2dj
               pt2d(2,ijpj+jl) = psgn * pt2d(3,ijpj-1+jl)
               DO ji = 3, jpiglo
                  iju = jpiglo - ji + 3
                  pt2d(ji,ijpj+jl) = psgn * pt2d(iju,ijpj-1-jl)
               END DO
            END DO
         CASE ( 'J' )                                     ! first ice U-V point
            DO jl =0, ipr2dj
               pt2d(2,ijpj+jl) = psgn * pt2d(3,ijpj-1+jl)
               DO ji = 3, jpiglo
                  iju = jpiglo - ji + 3
                  pt2d(ji,ijpj+jl) = psgn * pt2d(iju,ijpj-1-jl)
               END DO
            END DO
         CASE ( 'K' )                                     ! second ice U-V point
            DO jl =0, ipr2dj
               pt2d(2,ijpj+jl) = psgn * pt2d(3,ijpj-1+jl)
               DO ji = 3, jpiglo
                  iju = jpiglo - ji + 3
                  pt2d(ji,ijpj+jl) = psgn * pt2d(iju,ijpj-1-jl)
               END DO
            END DO
         END SELECT
         !
      CASE ( 5, 6 )                        ! *  North fold  F-point pivot
         !
         SELECT CASE ( cd_type )
         CASE ( 'T' , 'W' )                               ! T-, W-point
            DO jl = 0, ipr2dj
               DO ji = 1, jpiglo
                  ijt = jpiglo-ji+1
                  pt2d(ji,ijpj+jl) = psgn * pt2d(ijt,ijpj-1-jl)
               END DO
            END DO
         CASE ( 'U' )                                     ! U-point
            DO jl = 0, ipr2dj
               DO ji = 1, jpiglo-1
                  iju = jpiglo-ji
                  pt2d(ji,ijpj+jl) = psgn * pt2d(iju,ijpj-1-jl)
               END DO
            END DO
            pt2d(jpiglo,ijpj) = psgn * pt2d(1,ijpj-1)
         CASE ( 'V' )                                     ! V-point
            DO jl = 0, ipr2dj
               DO ji = 1, jpiglo
                  ijt = jpiglo-ji+1
                  pt2d(ji,ijpj+jl) = psgn * pt2d(ijt,ijpj-2-jl)
               END DO
            END DO
            DO ji = jpiglo/2+1, jpiglo
               ijt = jpiglo-ji+1
               pt2d(ji,ijpjm1) = psgn * pt2d(ijt,ijpjm1)
            END DO
         CASE ( 'F' )                               ! F-point
            DO jl = 0, ipr2dj
               DO ji = 1, jpiglo-1
                  iju = jpiglo-ji
                  pt2d(ji,ijpj+jl) = psgn * pt2d(iju,ijpj-2-jl)
               END DO
            END DO
            pt2d(jpiglo,ijpj) = psgn * pt2d(1,ijpj-2)
            DO ji = jpiglo/2+1, jpiglo-1
               iju = jpiglo-ji
               pt2d(ji,ijpjm1) = psgn * pt2d(iju,ijpjm1)
            END DO
         CASE ( 'I' )                                  ! ice U-V point (I-point)
            pt2d( 2 ,ijpj:ijpj+ipr2dj) = zero
            DO jl = 0, ipr2dj
               DO ji = 2 , jpiglo-1
                  ijt = jpiglo - ji + 2
                  pt2d(ji,ijpj+jl)= half * ( pt2d(ji,ijpj-1-jl) + psgn * pt2d(ijt,ijpj-1-jl) )
               END DO
            END DO
         CASE ( 'J' )                                  ! first ice U-V point
            pt2d( 2 ,ijpj:ijpj+ipr2dj) = zero
            DO jl = 0, ipr2dj
               DO ji = 2 , jpiglo-1
                  ijt = jpiglo - ji + 2
                  pt2d(ji,ijpj+jl)= pt2d(ji,ijpj-1-jl)
               END DO
            END DO
         CASE ( 'K' )                                  ! second ice U-V point
            pt2d( 2 ,ijpj:ijpj+ipr2dj) = zero
            DO jl = 0, ipr2dj
               DO ji = 2 , jpiglo-1
                  ijt = jpiglo - ji + 2
                  pt2d(ji,ijpj+jl)= pt2d(ijt,ijpj-1-jl)
               END DO
            END DO
         END SELECT
         !
      CASE DEFAULT                           ! *  closed : the code probably never go through
         !
         SELECT CASE ( cd_type)
         CASE ( 'T' , 'U' , 'V' , 'W' )                 ! T-, U-, V-, W-points
            pt2d(:, 1:1-ipr2dj     ) = zero
            pt2d(:,ijpj:ijpj+ipr2dj) = zero
         CASE ( 'F' )                                   ! F-point
            pt2d(:,ijpj:ijpj+ipr2dj) = zero
         CASE ( 'I' )                                   ! ice U-V point
            pt2d(:, 1:1-ipr2dj     ) = zero
            pt2d(:,ijpj:ijpj+ipr2dj) = zero
         CASE ( 'J' )                                   ! first ice U-V point
            pt2d(:, 1:1-ipr2dj     ) = zero
            pt2d(:,ijpj:ijpj+ipr2dj) = zero
         CASE ( 'K' )                                   ! second ice U-V point
            pt2d(:, 1:1-ipr2dj     ) = zero
            pt2d(:,ijpj:ijpj+ipr2dj) = zero
         END SELECT
         !
      END SELECT
      !
   END SUBROUTINE lbc_nfd

   SUBROUTINE lbc_lnk( pt2d, cd_type, psgn, cd_mpp, pval )
      !!---------------------------------------------------------------------
      !!                 ***  ROUTINE lbc_lnk  ***
      !!
      !! ** Purpose :   set lateral boundary conditions on a 2D array (non mpp case)
      !!
      !! ** Method  :   psign = -1 :    change the sign across the north fold
      !!                      =  1 : no change of the sign across the north fold
      !!                      =  0 : no change of the sign across the north fold and
      !!                             strict positivity preserved: use inner row/column
      !!                             for closed boundaries.
      !!----------------------------------------------------------------------
      CHARACTER(len=1)            , INTENT(in   )           ::   cd_type   ! nature of pt3d grid-points
      REAL(wp), DIMENSION(jpi,jpj), INTENT(inout)           ::   pt2d      ! 2D array on which the lbc is applied
      REAL(wp)                    , INTENT(in   )           ::   psgn      ! control of the sign 
      CHARACTER(len=3)            , INTENT(in   ), OPTIONAL ::   cd_mpp    ! MPP only (here do nothing)
      REAL(wp)                    , INTENT(in   ), OPTIONAL ::   pval      ! background value (for closed boundaries)
      !!
      REAL(wp) ::   zland
      !!----------------------------------------------------------------------

      IF( PRESENT( pval ) ) THEN   ;   zland = pval      ! set land value (zero by default)
      ELSE                         ;   zland = zero
      ENDIF

      IF (PRESENT(cd_mpp)) THEN
         ! only fill the overlap area and extra allows 
         ! this is in mpp case. In this module, just do nothing
      ELSE      
         !
         !                                     ! East-West boundaries
         !                                     ! ====================
         SELECT CASE ( nperio )
         !
         CASE ( 1 , 4 , 6 )                       !** cyclic east-west
            pt2d( 1 ,:) = pt2d(jpim1,:)               ! all points
            pt2d(jpi,:) = pt2d(  2  ,:)
            !
         CASE DEFAULT                             !** East closed  --  West closed
            SELECT CASE ( cd_type )
            CASE ( 'T' , 'U' , 'V' , 'W' )            ! T-, U-, V-, W-points
               pt2d( 1 ,:) = zland
               pt2d(jpi,:) = zland
            CASE ( 'F' )                              ! F-point
               pt2d(jpi,:) = zland
            END SELECT
            !
         END SELECT
         !
         !                                     ! North-South boundaries
         !                                     ! ======================
         SELECT CASE ( nperio )
         !
         CASE ( 2 )                               !**  South symmetric  --  North closed
            SELECT CASE ( cd_type )
            CASE ( 'T' , 'U' , 'W' )                   ! T-, U-, W-points
               pt2d(:, 1 ) = pt2d(:,3)
               pt2d(:,jpj) = zland
            CASE ( 'V' , 'F' )                         ! V-, F-points
               pt2d(:, 1 ) = psgn * pt2d(:,2)
               pt2d(:,jpj) = zland
            END SELECT
            !
         CASE ( 3 , 4 , 5 , 6 )                   !**  North fold  T or F-point pivot  --  South closed
            SELECT CASE ( cd_type )                    ! South : closed
            CASE ( 'T' , 'U' , 'V' , 'W' , 'I' )             ! all points except F-point
               pt2d(:, 1 ) = zland
            END SELECT
            !                                          ! North fold
            CALL lbc_nfd( pt2d(:,:), cd_type, psgn )
            !
         CASE DEFAULT                             !**  North closed  --  South closed
            SELECT CASE ( cd_type )
            CASE ( 'T' , 'U' , 'V' , 'W' )             ! T-, U-, V-, W-points
               pt2d(:, 1 ) = zland
               pt2d(:,jpj) = zland
            CASE ( 'F' )                               ! F-point
               pt2d(:,jpj) = zland
            END SELECT
            !
         END SELECT
         !
      ENDIF
      !    
   END SUBROUTINE lbc_lnk
   
   SUBROUTINE angle
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE angle  ***
      !! 
      !! ** Purpose :   Compute angles between model grid lines and the North direction
      !!
      !! ** Method  :
      !!
      !! ** Action  :   Compute (gsint, gcost, gsinu, gcosu, gsinv, gcosv, gsinf, gcosf) arrays:
      !!      sinus and cosinus of the angle between the north-south axe and the 
      !!      j-direction at t, u, v and f-points
      !!
      !! History :
      !!   7.0  !  96-07  (O. Marti )  Original code
      !!   8.0  !  98-06  (G. Madec )
      !!   8.5  !  98-06  (G. Madec )  Free form, F90 + opt.
      !!   9.2  !  07-04  (S. Masson)  Add T, F points and bugfix in cos lateral boundary
      !!----------------------------------------------------------------------
      INTEGER ::   ji, jj   ! dummy loop indices
      INTEGER ::   ierr     ! local integer
      REAL(wp) ::   &
         zlam, zphi,            &  ! temporary scalars
         zlan, zphh,            &  !    "         "
         zxnpt, zynpt, znnpt,   &  ! x,y components and norm of the vector: T point to North Pole
         zxnpu, zynpu, znnpu,   &  ! x,y components and norm of the vector: U point to North Pole
         zxnpv, zynpv, znnpv,   &  ! x,y components and norm of the vector: V point to North Pole
         zxnpf, zynpf, znnpf,   &  ! x,y components and norm of the vector: F point to North Pole
         zxvvt, zyvvt, znvvt,   &  ! x,y components and norm of the vector: between V points below and above a T point
         zxffu, zyffu, znffu,   &  ! x,y components and norm of the vector: between F points below and above a U point
         zxffv, zyffv, znffv,   &  ! x,y components and norm of the vector: between F points left  and right a V point
         zxuuf, zyuuf, znuuf       ! x,y components and norm of the vector: between U points below and above a F point
      !!----------------------------------------------------------------------

      ! ============================= !
      ! Compute the cosinus and sinus !
      ! ============================= !
      ! (computation done on the north stereographic polar plane)

      DO jj = 2, jpjm1
!CDIR NOVERRCHK
         DO ji = 2, jpi   ! vector opt.

            ! north pole direction & modulous (at t-point)
            zlam = glamt(ji,jj)
            zphi = gphit(ji,jj)
            zxnpt = zero - two * COS( deg2rad*zlam ) * TAN( rpi*quarter - deg2rad*zphi*half )
            zynpt = zero - two * SIN( deg2rad*zlam ) * TAN( rpi*quarter - deg2rad*zphi*half )
            znnpt = zxnpt*zxnpt + zynpt*zynpt

            ! north pole direction & modulous (at u-point)
            zlam = glamu(ji,jj)
            zphi = gphiu(ji,jj)
            zxnpu = zero - two * COS( deg2rad*zlam ) * TAN( rpi*quarter - deg2rad*zphi*half )
            zynpu = zero - two * SIN( deg2rad*zlam ) * TAN( rpi*quarter - deg2rad*zphi*half )
            znnpu = zxnpu*zxnpu + zynpu*zynpu

            ! north pole direction & modulous (at v-point)
            zlam = glamv(ji,jj)
            zphi = gphiv(ji,jj)
            zxnpv = zero - two * COS( deg2rad*zlam ) * TAN( rpi*quarter - deg2rad*zphi*half )
            zynpv = zero - two * SIN( deg2rad*zlam ) * TAN( rpi*quarter - deg2rad*zphi*half )
            znnpv = zxnpv*zxnpv + zynpv*zynpv

            ! north pole direction & modulous (at f-point)
            zlam = glamf(ji,jj)
            zphi = gphif(ji,jj)
            zxnpf = zero - two * COS( deg2rad*zlam ) * TAN( rpi*quarter - deg2rad*zphi*half )
            zynpf = zero - two * SIN( deg2rad*zlam ) * TAN( rpi*quarter - deg2rad*zphi*half )
            znnpf = zxnpf*zxnpf + zynpf*zynpf

            ! j-direction: v-point segment direction (around t-point)
            zlam = glamv(ji,jj  )
            zphi = gphiv(ji,jj  )
            zlan = glamv(ji,jj-1)
            zphh = gphiv(ji,jj-1)
            zxvvt =  two * COS( deg2rad*zlam ) * TAN( rpi*quarter - deg2rad*zphi*half )   &
               &  -  two * COS( deg2rad*zlan ) * TAN( rpi*quarter - deg2rad*zphh*half )
            zyvvt =  two * SIN( deg2rad*zlam ) * TAN( rpi*quarter - deg2rad*zphi*half )   &
               &  -  two * SIN( deg2rad*zlan ) * TAN( rpi*quarter - deg2rad*zphh*half )
            znvvt = SQRT( znnpt * ( zxvvt*zxvvt + zyvvt*zyvvt )  )
            znvvt = MAX( znvvt, eps2 )

            ! j-direction: f-point segment direction (around u-point)
            zlam = glamf(ji,jj  )
            zphi = gphif(ji,jj  )
            zlan = glamf(ji,jj-1)
            zphh = gphif(ji,jj-1)
            zxffu =  two * COS( deg2rad*zlam ) * TAN( rpi*quarter - deg2rad*zphi*half )   &
               &  -  two * COS( deg2rad*zlan ) * TAN( rpi*quarter - deg2rad*zphh*half )
            zyffu =  two * SIN( deg2rad*zlam ) * TAN( rpi*quarter - deg2rad*zphi*half )   &
               &  -  two * SIN( deg2rad*zlan ) * TAN( rpi*quarter - deg2rad*zphh*half )
            znffu = SQRT( znnpu * ( zxffu*zxffu + zyffu*zyffu )  )
            znffu = MAX( znffu, eps2 )

            ! i-direction: f-point segment direction (around v-point)
            zlam = glamf(ji  ,jj)
            zphi = gphif(ji  ,jj)
            zlan = glamf(ji-1,jj)
            zphh = gphif(ji-1,jj)
            zxffv =  two * COS( deg2rad*zlam ) * TAN( rpi*quarter - deg2rad*zphi*half )   &
               &  -  two * COS( deg2rad*zlan ) * TAN( rpi*quarter - deg2rad*zphh*half )
            zyffv =  two * SIN( deg2rad*zlam ) * TAN( rpi*quarter - deg2rad*zphi*half )   &
               &  -  two * SIN( deg2rad*zlan ) * TAN( rpi*quarter - deg2rad*zphh*half )
            znffv = SQRT( znnpv * ( zxffv*zxffv + zyffv*zyffv )  )
            znffv = MAX( znffv, eps2 )

            ! j-direction: u-point segment direction (around f-point)
            zlam = glamu(ji,jj+1)
            zphi = gphiu(ji,jj+1)
            zlan = glamu(ji,jj  )
            zphh = gphiu(ji,jj  )
            zxuuf =  two * COS( deg2rad*zlam ) * TAN( rpi*quarter - deg2rad*zphi*half )   &
               &  -  two * COS( deg2rad*zlan ) * TAN( rpi*quarter - deg2rad*zphh*half )
            zyuuf =  two * SIN( deg2rad*zlam ) * TAN( rpi*quarter - deg2rad*zphi*half )   &
               &  -  two * SIN( deg2rad*zlan ) * TAN( rpi*quarter - deg2rad*zphh*half )
            znuuf = SQRT( znnpf * ( zxuuf*zxuuf + zyuuf*zyuuf )  )
            znuuf = MAX( znuuf, eps2 )

            ! cosinus and sinus using scalar and vectorial products
            gsint(ji,jj) = ( zxnpt*zyvvt - zynpt*zxvvt ) / znvvt
            gcost(ji,jj) = ( zxnpt*zxvvt + zynpt*zyvvt ) / znvvt

            gsinu(ji,jj) = ( zxnpu*zyffu - zynpu*zxffu ) / znffu
            gcosu(ji,jj) = ( zxnpu*zxffu + zynpu*zyffu ) / znffu

            gsinf(ji,jj) = ( zxnpf*zyuuf - zynpf*zxuuf ) / znuuf
            gcosf(ji,jj) = ( zxnpf*zxuuf + zynpf*zyuuf ) / znuuf

            ! (caution, rotation of 90 degres)
            gsinv(ji,jj) = ( zxnpv*zxffv + zynpv*zyffv ) / znffv
            gcosv(ji,jj) =-( zxnpv*zyffv - zynpv*zxffv ) / znffv

         END DO
      END DO

      ! =============== !
      ! Geographic mesh !
      ! =============== !

      DO jj = 2, jpjm1
         DO ji = 2, jpi   ! vector opt.
            IF( MOD( ABS( glamv(ji,jj) - glamv(ji,jj-1) ), dcircle ) < eps1 ) THEN
               gsint(ji,jj) = zero
               gcost(ji,jj) = one
            ENDIF
            IF( MOD( ABS( glamf(ji,jj) - glamf(ji,jj-1) ), dcircle ) < eps1 ) THEN
               gsinu(ji,jj) = zero
               gcosu(ji,jj) = one
            ENDIF
            IF(      ABS( gphif(ji,jj) - gphif(ji-1,jj) )            < eps1 ) THEN
               gsinv(ji,jj) = zero
               gcosv(ji,jj) = one
            ENDIF
            IF( MOD( ABS( glamu(ji,jj) - glamu(ji,jj+1) ), dcircle ) < eps1 ) THEN
               gsinf(ji,jj) = zero
               gcosf(ji,jj) = one
            ENDIF
         END DO
      END DO

      gsint(:,:) = MAX(MIN(gsint(:,:), one), -one)
      gsinu(:,:) = MAX(MIN(gsinu(:,:), one), -one)
      gsinv(:,:) = MAX(MIN(gsinv(:,:), one), -one)
      gsinf(:,:) = MAX(MIN(gsinf(:,:), one), -one)
      gcost(:,:) = MAX(MIN(gcost(:,:), one), -one)
      gcosu(:,:) = MAX(MIN(gcosu(:,:), one), -one)
      gcosv(:,:) = MAX(MIN(gcosv(:,:), one), -one)
      gcosf(:,:) = MAX(MIN(gcosf(:,:), one), -one)

      ! =========================== !
      ! Lateral boundary conditions !
      ! =========================== !

      ! lateral boundary cond.: T-, U-, V-, F-pts, sgn
      CALL lbc_lnk( gcost, 'T', -one )   ;   CALL lbc_lnk( gsint, 'T', -one )
      CALL lbc_lnk( gcosu, 'U', -one )   ;   CALL lbc_lnk( gsinu, 'U', -one )
      CALL lbc_lnk( gcosv, 'V', -one )   ;   CALL lbc_lnk( gsinv, 'V', -one )
      CALL lbc_lnk( gcosf, 'F', -one )   ;   CALL lbc_lnk( gsinf, 'F', -one )

   END SUBROUTINE angle

   SUBROUTINE nf_handle_err(stat, line)

    INTEGER, INTENT(IN) :: stat, line

    IF (stat /= NF_NOERR) THEN
      WRITE(stderr,FMT='(A)') nf_strerror(stat)
      WRITE(stderr,FMT='(A,I4)') 'at line: ', line
      STOP
    END IF

   end subroutine nf_handle_err
 
   function lower(in) RESULT (out)
     implicit none
     character (*), intent(in)  :: in
     character(:), allocatable  :: out
     integer                    :: i, j
     out = in                        !transfer whole array
     do i = 1, LEN_TRIM(out)         !each character
       j = iachar(out(i:i))         !get the ASCII position
       select case (j)
       case (65:90)             !The upper case characters
         out(i:i) = ACHAR(j+32) !Offset to the lower case
       end select
     end do
   end function lower

   function upper(in) RESULT (out)
     implicit none
     character (*), intent(in)  :: in
     character(:), allocatable  :: out
     integer                    :: i, j
     out = in                        !transfer whole array
     do i = 1, LEN_TRIM(out)         !each character
       j = iachar(out(i:i))         !get the ASCII position
       select case (j)
       case (97:122)             !The lower case characters
         out(i:i) = ACHAR(j-32) !Offset to the upper case
       end select
     end do
   end function upper

END PROGRAM nemo_angle
