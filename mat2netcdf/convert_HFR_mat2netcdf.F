      program conv_rad

c  This program will read HIOOS HFR data as supplied by P. Flament's group 
c  and produce ascii files to be ingested by the National HRF network.  The
c  specs were provided by Mark Otero at Scripps.  This program is based on
c  the Matlab supplied "matdemo2.F"
c  
c  Note:  The RADLAB group provides files every 15 minutes, and these are
c  archived at /export/lawelawe2/radlab/hioos.  The files contain the following:
c     tzone       3-character time zone (UT) 
c     t           date/time in decimal days since 0/0/0000
c     slat        latitude of radar (decimal deg)
c     slon        longitude of radar (decimal deg)
c     r(jm)       radial distance of cell from radar (km)
c     theta(im)   bearing of cell from radar (deg)
c     Acc(im,jm)  accuracy of the Bragg estimate in cell
c     Pow(im,jm)  relative power in dB of cell
c     Ur(im,jm)   radial current of cell in m/s, positive away from radar
c     Var(im,jm)  variance of grid interpolation in cell
c     X(im,jm)    longitude of cell (decimal deg)
c     Y(im,jm)    latitude of cell (decimal deg)
c
c
c
c  jimp 12-03-2009
c  jimp 06-14-2010: changes suggested by Mark Otero, including removing
c                   negative ranges
c  jimp 08-15-2010: looks like the grid changed on August 1; added ifdef's
c                   to accomodate smaller grids
c  jimp 09-10-2010: added support for KAL site
c  jimp 03-11-2013: added support for KNA site
c  jimp 03-12-2014: added support for KKH site
c  jimp 05-12-2015: added support for PPK site
c  jimp 09-15-2015: added support for KAP site
c  jimp 04-28-2016: major redo to add a buffer to range and bearing; these
c                   change without warning, so now we have a max setting,
c                   fill all with missing, then add valid data where appr-
c                   opriate, rather than manually setting these.

#include "netcdf.inc"
#include "fintrf.h"

      parameter ( imax = 245, jmax = 125 )
      parameter ( imaxpts = imax*jmax, nvars = 4 )

      mwpointer matOpen, matGetDir, matGetVariableInfo, mxGetM, mxGetN
      mwpointer mxGetNumberOfElements, mxGetNumberOfDimensions
      mwpointer matGetVariable, mxGetData
      mwpointer mfp, dir, adir(100), pa, pm
      mwsize ndir
      integer   i, stat, ndim, stime
      integer   rnge_dim, bear_dim, lev_dim, itim_dim
      integer   lon_id, lat_id, lev_id, tim_id, bear_id, rnge_id
      integer   count(4), start(4), dims(4), var_id(nvars)
      integer*8 ivmax
      integer*4 iindex
      real*8    data(imaxpts)
      character*32 name(100)
      character*3 location
      character*24 outfile
      character* 40 TITLE01, TITLE02, TITLE03, TITLE04, TITLE05, TITLE06
      character*240 TITLEA
      character*4 cvar(nvars)
      real X(imaxpts), Y(imaxpts), EACC(imaxpts), VELO(imaxpts)
      real varian(imaxpts), POW(imaxpts), r(imaxpts), theta(imaxpts)
      real bear(imax), rnge(jmax)
      real xlon(imax,jmax), ylat(imax,jmax)
      real var(nvars,imax,jmax), var2(imax,jmax,1,1)
       
c read time from direct input

      xmiss = -999.0
      pi = 4.0 * atan ( 1.0 )
      zlev = 0.0
      read ( *, 18 ) iyear, iday, ihr, imin, location
      write ( outfile, 19 ) 'RDL_', location, "_", iyear, "_",
     +                iday, "_", ihr, imin, ".nc"
      iy = iyear - 2008
      stime = iy * 365 * 24 * 60 + ( iyear - 2005 ) / 4 * 24 * 60 +
     +              ( iday - 1 )  * 24 * 60 + ihr * 60 + imin

      if ( location.eq.'kok' ) then
         write ( TITLE01, 29) "Results from Koko Head HFR, raw 15-minut"
         write ( TITLE02, 29) "e data provided by SOEST RADLAB group, P"
         write ( TITLE03, 29) ". Flament, PI.  Data are radial velociti"
         write ( TITLE04, 29) "es......................................"
         write ( TITLE05, 29) "x-origin = -157.70330..................."
         write ( TITLE06, 29) "y-origin =   21.26139..................."
         write (TITLEA, 39) TITLE01, TITLE02, TITLE03, TITLE04, TITLE05,
     +                    TITLE06
      endif
      if ( location.eq.'kak' ) then
         write ( TITLE01, 29) "Results from Kakaako HFR, raw 15-minute "
         write ( TITLE02, 29) "data provided by SOEST RADLAB group, P. "
         write ( TITLE03, 29) "Flament, PI.  Data are radial velocities"
         write ( TITLE04, 29) "x-origin = -157.85830..................."
         write ( TITLE05, 29) "y-origin =   21.29170..................."
         write ( TITLE06, 29) "........................................"
         write (TITLEA, 39) TITLE01, TITLE02, TITLE03, TITLE04, TITLE05,
     +                    TITLE06
      endif
      if ( location.eq.'kal' ) then
         write ( TITLE01, 29) "Results from Kalaeloa HFR, raw 15-minute"
         write ( TITLE02, 29) " data provided by SOEST RADLAB group, P."
         write ( TITLE03, 29) " Flament, PI.  Data are radial velocitie"
         write ( TITLE04, 29) "s......................................."
         write ( TITLE05, 29) "x-origin = -157.08360..................."
         write ( TITLE06, 29) "y-origin =   21.29750..................."
         write (TITLEA, 39) TITLE01, TITLE02, TITLE03, TITLE04, TITLE05,
     +                    TITLE06
      endif
      if ( location.eq.'kna' ) then
         write ( TITLE01, 29) "Results from Kaena Pt HFR, raw 15-minute"
         write ( TITLE02, 29) " data provided by SOEST RADLAB group, P."
         write ( TITLE03, 29) " Flament, PI.  Data are radial velocitie"
         write ( TITLE04, 29) "s......................................."
         write ( TITLE05, 29) "x-origin = -158.2897...................."
         write ( TITLE06, 29) "y-origin =   21.57500..................."
         write (TITLEA, 39) TITLE01, TITLE02, TITLE03, TITLE04, TITLE05,
     +                    TITLE06
      endif
      if ( location.eq.'kkh' ) then
         write ( TITLE01, 29) "Results from Keaukaha HFR, raw 15-minute"
         write ( TITLE02, 29) " data provided by SOEST RADLAB group, P."
         write ( TITLE03, 29) " Flament, PI.  Data are radial velocitie"
         write ( TITLE04, 29) "s......................................."
         write ( TITLE05, 29) "x-origin = -155.04830..................."
         write ( TITLE06, 29) "y-origin =   19.73200..................."
         write (TITLEA, 39) TITLE01, TITLE02, TITLE03, TITLE04, TITLE05,
     +                    TITLE06
      endif
      if ( location.eq.'ppk' ) then
         write ( TITLE01, 29) "Results from Pepeekeo HFR, raw 15-minute"
         write ( TITLE02, 29) " data provided by SOEST RADLAB group, P."
         write ( TITLE03, 29) " Flament, PI.  Data are radial velocitie"
         write ( TITLE04, 29) "s......................................."
         write ( TITLE05, 29) "x-origin = -155.08280..................."
         write ( TITLE06, 29) "y-origin =   19.84720..................."
         write (TITLEA, 39) TITLE01, TITLE02, TITLE03, TITLE04, TITLE05,
     +                    TITLE06
      endif
      if ( location.eq.'kap' ) then
         write ( TITLE01, 29) "Results from Kapolei HFR, raw 15-minute "
         write ( TITLE02, 29) "data provided by SOEST RADLAB group, P. "
         write ( TITLE03, 29) "Flament, PI.  Data are radial velocities"
         write ( TITLE04, 29) "........................................"
         write ( TITLE05, 29) "x-origin = -158.1169...................."
         write ( TITLE06, 29) "y-origin =   21.3114...................."
         write (TITLEA, 39) TITLE01, TITLE02, TITLE03, TITLE04, TITLE05,
     +                    TITLE06
      endif
18    format ( i4, 1x, i3.3, 1x, i2.2, 1x, i2.2, 1x, a3 )
19    format ( a4, a3, a1, i4, a1, 
     +         i3.3, a1, i2.2, i2.2, a3 )
29    format ( a40 )
39    format ( 6(a40) )

c set all values to missing

      do 20 i = 1, imax
      do 20 j = 1, jmax
         bear(i) = xmiss
         rnge(j) = xmiss
         xlon(i,j) = xmiss
         ylat(i,j) = xmiss
         var2(i,j,1,1) = xmiss
         do 30 k = 1, 4
            var(k,i,j) = xmiss
30       continue
20    continue

c open the file and directory of arrays, "ndir" is essentially the number
c of variables

      mfp = matOpen ( 'radar_in.mat', 'r' )
      if ( mfp .eq. 0 ) then
         write ( 6, * ) 'Can''t open input file'
         stop
      end if
      dir = matgetdir ( mfp, iindex )
      if ( dir .eq. 0 ) then
         write ( 6, * ) 'Can''t read directory.'
         stop
      endif

c  copy pointer into an array of pointers

      ndir = iindex
      call mxCopyPtrToPtrArray ( dir, adir, ndir )

c  copy pointer to character string

      do 40 i = 1, int ( ndir )
         call mxCopyPtrToCharacter ( adir(i), name(i), 32 )
40    continue

c  read data for each array; there are "ndir" arrays, and there are
c  "ivmax" points in each array

      do 50 j = 1, int ( ndir )
         pa = matGetVariableInfo ( mfp, name(j) )
         jm = mxGetM ( pa )
         im = mxGetN ( pa )
         ivmax = mxGetNumberOfElements ( pa )
         ndim = mxGetNumberOfDimensions ( pa )
         pm = matGetVariable ( mfp, name(j) )
         call mxCopyPtrToReal8 ( mxGetData(pm), data, ivmax )
         do 60 i = 1, int ( ivmax )
            if ( name(j).eq.'X' ) then
                X(i) = data(i)
                ibearmax = im
                jrngemax = jm
            elseif ( name(j).eq.'Y' ) then
                Y(i) = data(i)
            elseif ( name(j).eq.'Ur' ) then
               VELO(i) = data(i)
            elseif ( name(j).eq.'theta' ) then
               theta(i) = data(i)
            elseif ( name(j).eq.'r' ) then
               r(i) = data(i)
            elseif ( name(j).eq.'Acc' ) then
               EACC(i) = data(i)
            elseif ( name(j).eq.'Pow' ) then
               POW(i) = data(i)
            elseif ( name(j).eq.'Var' ) then
               varian(i) = data(i)
            elseif ( name(j).eq.'slat' ) then
               slat = data(i)
            elseif ( name(j).eq.'slon' ) then
               slon = data(i)
            endif
60       continue
         call mxDestroyArray(pa)
         call mxDestroyArray(pm)
50    continue

      do 70 j = 1, jrngemax
         if ( r(j).lt.0.0 ) then
            rnge(j) = 0.0
         else
            rnge(j) = r(j)
         endif
         if ( rnge(j).ne.rnge(j) ) rnge(j) = xmiss
70    continue
      do 80 i = 1, ibearmax
         bear(i) =  theta(i)
         if ( bear(i).ne.bear(i) ) bear(i) = xmiss
80    continue

      ij = 1
      do 90 i = 1, ibearmax
      do 90 j = 1, jrngemax
c check for NaN's (actually check that it's not a NaN)
         if ( X(ij).eq.X(ij) ) then
            xlon(i,j) = X(ij)
         else
            xlon(i,j) = xmiss
         endif
         if ( Y(ij).eq.Y(ij) ) then
            ylat(i,j) = Y(ij)
         else
            ylat(i,j) = xmiss
         endif
         if ( EACC(ij).eq.EACC(ij) ) then
            var(1,i,j) = EACC(ij)
         else
            var(1,i,j) = xmiss
         endif
         if ( POW(ij).eq.POW(ij) ) then
            var(2,i,j) = POW(ij)
         else
            var(2,i,j) = xmiss
         endif
         if ( VELO(ij).eq.VELO(ij) ) then
            var(3,i,j) = VELO(ij)
         else
            var(3,i,j) = xmiss
         endif
         if ( varian(ij).eq.varian(ij) ) then
            var(4,i,j) = varian(ij)
         else
            var(4,i,j) = xmiss
         endif
         ij = ij + 1
90    continue

c  close file

      stat = matClose ( mfp )
      if (stat .ne. 0) then
         write(6,*) 'Error closing ''matdemo.mat''.'
         stop
      end if

c  make output file

c Define variables

      write ( cvar(1), 49 ) "accu"
      write ( cvar(2), 49 ) "powr"
      write ( cvar(3), 49 ) "urad"
      write ( cvar(4), 49 ) "vari"
49    format ( a4 )

c Open and initialize the output netcdf files

      print *, "Opening file ", outfile
      iret = nf_create ( outfile, NF_CLOBBER, ncid )
      call check_err ( iret )

      iret = nf_put_att_text ( ncid, nf_global, 'title', 240, TITLEA )
      call check_err ( iret )

      iret = nf_put_att_text ( ncid, nf_global, 'Conventions',
     +                         6, 'CF-1.4' )


c Define dimensions and variables

      print *, "Defining dimensions"
      iret = nf_def_dim ( ncid, 'range', jmax, rnge_dim )
      iret = nf_def_dim ( ncid, 'bearing', imax, bear_dim )
      iret = nf_def_dim ( ncid, 'z', 1, lev_dim )
      iret = nf_def_dim ( ncid, 'time', NF_UNLIMITED, itim_dim )
      call check_err ( iret )

      print *, "Defining variables"
c 0-D variables
      iret = nf_def_var ( ncid, 'x-origin', NF_REAL, 0, 0, islon_id )
      iret = nf_def_var ( ncid, 'y-origin', NF_REAL, 0, 0, islat_id )
      call check_err ( iret )
c 1-D variables
      dims(1) = rnge_dim
      iret = nf_def_var ( ncid, 'range', NF_REAL, 1, dims, rnge_id )
      dims(1) = bear_dim
      iret = nf_def_var ( ncid, 'bearing', NF_REAL, 1, dims, bear_id )
      dims(1) = lev_dim
      iret = nf_def_var ( ncid, 'z', NF_REAL, 1, dims, lev_id )
      dims(1) = itim_dim
      iret = nf_def_var ( ncid, 'time', NF_INT, 1, dims, tim_id )
      call check_err ( iret )

c 2-D variables
      dims(1) = bear_dim
      dims(2) = rnge_dim
      iret = nf_def_var ( ncid, 'lat', NF_REAL, 2, dims, lat_id )
      dims(1) = bear_dim
      dims(2) = rnge_dim
      iret = nf_def_var ( ncid, 'lon', NF_REAL, 2, dims, lon_id )

c 4-D variables (skip depth...)
      dims(4) = itim_dim
      dims(3) = lev_dim
      dims(2) = rnge_dim
      dims(1) = bear_dim
      do 100 i = 1, nvars
         iret = nf_def_var ( ncid, cvar(i), NF_REAL, 4, dims, var_id(i))
100   continue
      call check_err ( iret )

c  Assign attributes

      print *, "Assigning attributes"
c x-origin
      iret = nf_put_att_text ( ncid, islon_id, 'long_name', 30,
     +                         'longitude of HFR at zero range' )
      iret = nf_put_att_text ( ncid, islon_id, 'short_name', 4,
     +                         'slon' )
      iret = nf_put_att_text ( ncid, islon_id, 'units', 12,
     +                         'degrees_east' )
c y-origin
      iret = nf_put_att_text ( ncid, islat_id, 'long_name', 29,
     +                         'latitude of HFR at zero range' )
      iret = nf_put_att_text ( ncid, islat_id, 'short_name', 4,
     +                         'slat' )
      iret = nf_put_att_text ( ncid, islat_id, 'units', 13,
     +                         'degrees_north' )
c range
      iret = nf_put_att_text ( ncid, rnge_id, 'long_name', 13,
     +                         'range to cell' )
c      iret = nf_put_att_text ( ncid, rnge_id, 'standard_name', 13,
c     +                         'range_to_cell' )
c      iret = nf_put_att_text ( ncid, rnge_id, 'short_name', 5,
c     +                         'range' )
      iret = nf_put_att_text ( ncid, rnge_id, 'units', 2,
     +                         'km' )
      iret = nf_put_att_real ( ncid, rnge_id, '_FillValue',
     +                         NF_REAL, 1, xmiss )
c bearing
      iret = nf_put_att_text ( ncid, bear_id, 'long_name', 36,
     +                         'bearing to cell clockwise from north' )
c      iret = nf_put_att_text ( ncid, bear_id, 'standard_name', 36,
c     +                         'bearing_to_cell_clockwise_from_north' )
c      iret = nf_put_att_text ( ncid, bear_id, 'short_name', 7,
c     +                         'bearing' )
      iret = nf_put_att_text ( ncid, bear_id, 'units', 7,
     +                         'degrees' )
      iret = nf_put_att_real ( ncid, bear_id, '_FillValue',
     +                         NF_REAL, 1, xmiss )
c time
      iret = nf_put_att_text ( ncid, tim_id, 'long_name', 4,
     +                         'time' )
      iret = nf_put_att_text ( ncid, tim_id, 'standard_name', 4,
     +                         'time' )
c      iret = nf_put_att_text ( ncid, tim_id, 'short_name', 4,
c     +                         'time' )
      iret = nf_put_att_text ( ncid, tim_id, 'axis', 1, 't' )
      iret = nf_put_att_text ( ncid, tim_id, 'units', 33,
     +                         'minutes since 2008-01-01 00:00:00' )
c depth
      iret = nf_put_att_text ( ncid, lev_id, 'long_name', 26,
     +                         'depth below mean sea level' )
      iret = nf_put_att_text ( ncid, lev_id, 'standard_name', 5,
     +                         'depth' )
c      iret = nf_put_att_text ( ncid, lev_id, 'short_name', 5,
c     +                         'depth' )
      iret = nf_put_att_text ( ncid, lev_id, 'axis', 1, 'z' )
      iret = nf_put_att_text ( ncid, lev_id, 'units', 6,
     +                         'meters' )
      call check_err ( iret )
c latitude
      iret = nf_put_att_text ( ncid, lat_id, 'long_name', 8,
     +                         'latitude' )
      iret = nf_put_att_text ( ncid, lat_id, 'standard_name', 8,
     +                         'latitude' )
c      iret = nf_put_att_text ( ncid, lat_id, 'short_name', 3,
c     +                         'lat' )
      iret = nf_put_att_text ( ncid, lat_id, 'axis', 1, 'y' )
      iret = nf_put_att_text ( ncid, lat_id, 'units', 13,
     +                         'degrees_north' )
      iret = nf_put_att_real ( ncid, lat_id, '_FillValue',
     +                         NF_REAL, 1, xmiss )
      call check_err ( iret )
c longitude
      iret = nf_put_att_text ( ncid, lon_id, 'long_name', 9,
     +                         'longitude' )
      iret = nf_put_att_text ( ncid, lon_id, 'standard_name', 9,
     +                         'longitude' )
c      iret = nf_put_att_text ( ncid, lon_id, 'short_name', 3,
c     +                         'lon' )
      iret = nf_put_att_text ( ncid, lon_id, 'axis', 1, 'x' )
      iret = nf_put_att_text ( ncid, lon_id, 'units', 12,
     +                         'degrees_east' )
      iret = nf_put_att_real ( ncid, lon_id, '_FillValue',
     +                         NF_REAL, 1, xmiss )
      call check_err ( iret )
c accuracy
      iret = nf_put_att_text ( ncid, var_id(1), 'long_name', 38,
     +                'accuracy of the Bragg estimate in cell' )
c      iret = nf_put_att_text ( ncid, var_id(1), 'standard_name', 23,
c     +                               'bragg_estimate_accuracy' )
c      iret = nf_put_att_text ( ncid, var_id(1), 'short_name', 4,
c     +                         'accu' )
      iret = nf_put_att_text ( ncid, var_id(1), 'units', 3,
     +                         'm/s')
c      range(1) = 0.0
c      range(2) = 20.0
c      iret = nf_put_att_real ( ncid, var_id(1), 'valid_range',
c     +                           NF_REAL, 2, range )
      iret = nf_put_att_real ( ncid, var_id(1), '_FillValue',
     +                         NF_REAL, 1, xmiss )
      call check_err ( iret )
c power
      iret = nf_put_att_text ( ncid, var_id(2), 'long_name', 22,
     +                         'relative power of cell' )
c      iret = nf_put_att_text ( ncid, var_id(2), 'standard_name', 10,
c     +                                     'cell_power' )
c      iret = nf_put_att_text ( ncid, var_id(2), 'short_name', 4,
c     +                         'powr' )
      iret = nf_put_att_text ( ncid, var_id(2), 'units', 2,
     +                         'dB')
c      range(1) = 0.0
c      range(2) = 30.0
c      iret = nf_put_att_real ( ncid, var_id(2), 'valid_range',
c     +                           NF_REAL, 2, range )
      iret = nf_put_att_real ( ncid, var_id(2), '_FillValue',
     +                         NF_REAL, 1, xmiss )
      call check_err ( iret )
c radial current speed
      iret = nf_put_att_text ( ncid, var_id(3), 'long_name', 22,
     +                         'radial current in cell' )
      iret = nf_put_att_text ( ncid, var_id(3), 'standard_name', 46,
     +                'radial_sea_water_velocity_away_from_instrument' )
c      iret = nf_put_att_text ( ncid, var_id(3), 'short_name', 4,
c     +                         'urad' )
      iret = nf_put_att_text ( ncid, var_id(3), 'units', 3,
     +                         'm/s')
c      range(1) = 0.0
c      range(2) = 360.0
c      iret = nf_put_att_real ( ncid, var_id(3), 'valid_range',
c     +                           NF_REAL, 2, range )
      iret = nf_put_att_real ( ncid, var_id(3), '_FillValue',
     +                         NF_REAL, 1, xmiss )
      call check_err ( iret )
c variance
      iret = nf_put_att_text ( ncid, var_id(4), 'long_name', 30,
     +                         'variance of grid interpolation' )
c      iret = nf_put_att_text ( ncid, var_id(4), 'standard_name', 27,
c     +                         'grid_interpolation_variance' )
c      iret = nf_put_att_text ( ncid, var_id(4), 'short_name', 4,
c     +                         'vari' )
      iret = nf_put_att_text ( ncid, var_id(4), 'units', 7,
     +                         'unknown')
c      range(1) = 0.0
c      range(2) = 30.0
c      iret = nf_put_att_real ( ncid, var_id(4), 'valid_range',
c     +                           NF_REAL, 2, range )
      iret = nf_put_att_real ( ncid, var_id(4), '_FillValue',
     +                         NF_REAL, 1, xmiss )
      call check_err ( iret )

      iret = nf_enddef ( ncid )
      call check_err ( iret )

c Write data

      print*,"Writing data"
      start(1) = 1
      count(1) = imax
      iret = nf_put_vara_real ( ncid, bear_id, start, count, bear )
      call check_err ( iret )

      start(1) = 1
      count(1) = jmax
      iret = nf_put_vara_real ( ncid, rnge_id, start, count, rnge )
      call check_err ( iret )

      start(1) = 1
      count(1) = 1
      iret = nf_put_vara_int ( ncid, tim_id, start, count, stime )
      call check_err ( iret )

      start(1) = 1
      count(1) = 1
      iret = nf_put_vara_real ( ncid, lev_id, start, count, zlev )
      call check_err ( iret )

      start(1) = 1
      count(1) = 1
      iret = nf_put_vara_real ( ncid, islon_id, start, count, slon )
      call check_err ( iret )

      start(1) = 1
      count(1) = 1
      iret = nf_put_vara_real ( ncid, islat_id, start, count, slat )
      call check_err ( iret )

      start(1) = 1
      count(1) = imax
      start(2) = 1
      count(2) = jmax
      iret = nf_put_vara_real ( ncid, lat_id, start, count, ylat )
      call check_err ( iret )

      start(1) = 1
      count(1) = imax
      start(2) = 1
      count(2) = jmax
      iret = nf_put_vara_real ( ncid, lon_id, start, count, xlon )
      call check_err ( iret )

      do 110 ivar = 1, nvars
         start(1) = 1
         count(1) = imax
         start(2) = 1
         count(2) = jmax
         start(3) = 1
         count(3) = 1
         start(4) = 1
         count(4) = 1
         do 120 i = 1, imax
         do 120 j = 1, jmax
            var2(i,j,1,1) = var(ivar,i,j)
120      continue
         id = var_id(ivar)
         iret = nf_put_vara_real ( ncid, id, start, count, var2 )
         call check_err ( iret )
110   continue

      iret = nf_close ( ncid )
      call check_err ( iret )

      stop
      end

c  ----------------------------------------------------------------------------
      subroutine check_err(iret)
c  ----------------------------------------------------------------------------
      integer iret
      include 'netcdf.inc'
      if ( iret .ne. NF_NOERR ) then
         print *, nf_strerror ( iret )
         stop
      endif
      end
