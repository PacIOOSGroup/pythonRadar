      program conv_rad

c  This program will read HIOOS HFR data as supplied by P. Flament's group 
c  and produce ascii files to be ingested by the National HRF network.  The
c  specs were provided by Mark Otero at Scripps.  This program is based on
c  the Matlab supplied "matdemo2.F"
c  
c  Note:  The RADLAB group provides files periodically (was 15 minutes, then
c  hourly), and these are archived at /export/lawelawe2/radlab/hioos.  The 
c  files contain the following:
c     tzone       3-character time zone (UT) 
c     t           date/time in decimal days since 0/0/0000
c     slat        latitude of radar (decimal deg)
c     slon        longitude of radar (decimal deg)
c     RangeResolutionKMeters (RRKM)
c     TransmitCenterFreqMHz (TCFM)
c     DopplerResolutionHzPerBin (DRBH)
c     CurrentVelocityLimit (CVLM)
c     TransmitSweepRateHz (TSRH)
c     TransmitBandwidthKHz (TBKH)
c     SpectraDopplerCells (SDCS)
c     r(im)       radial distance of cell from radar (km)
c     theta(jm)   bearing of cell from radar (deg)
c     Acc(im,jm)  accuracy of the Bragg estimate in cell
c     Pow(im,jm)  relative power in dB of cell
c     Ur(im,jm)   radial current of cell in cm/s, positive away from radar
c     Var(im,jm)  variance of grid interpolation in cell
c     X(im,jm)    longitude of cell (decimal deg)
c     Y(im,jm)    latitude of cell (decimal deg)
c
c  These are then converted to the following on output:
c     X(imax)     longitude of grid cell
c     Y(imax)     latitude of grid cell
c     VELU(imax)  u-component of velocity
c     VELV(imax)  v-component of velocity
c     EACC(imax)  
c     XDST(imax)  x-distance of cell from radar location
c     YDST(imax)  y-distance of cell from radar location
c     RNGE(imax)  range
c     BEAR(imax)  bearing
c     VEL(imax)   along-track velocity, positive toward radar (opposite from input)
c     HEAD(imax)  current direction
c
c
c  jimp 12-17-2013: this program was copied from the one used for PacIOOS
c                   HFR; it is essentially the same but modified for the
c                   two non-PacIOOS sites VIL and QTA
c  jimp 05-11-2015: added additional site SAU
c  jimp 05-27-2015: added additional site CON
c  jimp 05-30-2015: read metadata from input files instead of hard-coding

#include "fintrf.h"

      parameter ( imaxpts = 22000 )
      mwpointer matOpen, matGetDir, matGetVariableInfo, mxGetM, mxGetN
      mwpointer mxGetNumberOfElements, mxGetNumberOfDimensions
      mwpointer matGetVariable, mxGetData, mxGetClassID
      mwpointer mfp, dir, adir(100), pa, pm
      mwsize ndir
      integer*4 iindex
      integer*8 imax, im, jm, idclass

      integer   i, stat, ndim
      real*8    data(imaxpts), xl, yl, rad_head
      character*32 name(100)
      character*42  cdata
      character*3 location
      character*29 outfile
      character*42 SITE, MANU, TSRH, SDCS, RRKM, TCFM, DRHB, CVLM, TBKH
      real X(imaxpts), Y(imaxpts), VELU(imaxpts), VELV(imaxpts)
      real EACC(imaxpts), XDST(imaxpts), YDST(imaxpts), RNGE(imaxpts)
      real BEAR(imaxpts), VELO(imaxpts), HEAD(imaxpts), XLON(imaxpts)
      real YLAT(imaxpts)
      real r(imaxpts), theta(imaxpts), ACC(imaxpts), VEL(imaxpts)
      real rad_lat, rad_lon, rad_bear, rang
       
c read time from direct input

      pi = 4.0 * atan ( 1.0 )
      read ( *, 18 ) iyear, imon, iday, ihr, imin, isec, location
      write ( outfile, 19 ) 'RDL_', location, "_", iyear, "_",
     +                imon, "_", iday, "_", ihr, imin, isec, ".ruv"
18    format ( i4, 1x, i2.2, 1x, i2.2, 1x, i2.2, 1x, i2.2, 1x, i2.2, 
     +         1x, a3 )
19    format ( a4, a3, a1, i4, a1, 
     +         i2.2, a1, i2.2, a1, i2.2, i2.2, i2.2, a4 )


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

      do 10 i = 1, int ( ndir )
         call mxCopyPtrToCharacter ( adir(i), name(i), 32 )
10    continue

c  read data for each array; there are "ndir" arrays, and there are
c  "imax" points in each array

      do 20 j = 1, int ( ndir )
         pa = matGetVariableInfo ( mfp, name(j) )
         im = mxGetM ( pa )
         jm = mxGetN ( pa )
         imax = mxGetNumberOfElements ( pa )
         ndim = mxGetNumberOfDimensions ( pa )
         pm = matGetVariable ( mfp, name(j) )
         idclass = mxGetClassID ( pm )
         if ( idclass.eq.6 ) then
            call mxCopyPtrToReal8 ( mxGetData(pm), data, imax )
            if ( name(j).eq.'X' ) then
               i2 = im
               j2 = jm
            endif
            do 30 i = 1, int ( imax )
               if ( name(j).eq.'slon' ) slon = real ( data(i) )
               if ( name(j).eq.'slat' ) slat = real ( data(i) )
               if ( name(j).eq.'X' ) X(i) = real ( data(i) )
               if ( name(j).eq.'Y' ) Y(i) = real ( data(i) )
               if ( name(j).eq.'Ur' ) VELO(i) = real ( data(i) )
               if ( name(j).eq.'theta' ) theta(i) = real ( data(i) )
               if ( name(j).eq.'r' ) r(i) = real ( data(i) )
               if ( name(j).eq.'Acc' ) EACC(i) = real ( data(i) )
30          continue
         endif
         if ( idclass.eq.4 ) then
            call mxGetString ( pm, cdata, imax )
            if ( name(j).eq.'Manufacturer' ) MANU = cdata
            if ( name(j).eq.'Site' ) SITE = cdata
            if ( name(j).eq.'RangeResolutionKMeters' ) RRKM = cdata
            if ( name(j).eq.'TransmitCenterFreqMHz' )  TCFM = cdata
            if ( name(j).eq.'DopplerResolutionHzPerBin' ) DRHB = cdata
            if ( name(j).eq.'CurrentVelocityLimit' ) CVLM = cdata
            if ( name(j).eq.'TransmitSweepRateHz' ) TSRH = cdata
            if ( name(j).eq.'TransmitBandwidthKHz' ) TBKH = cdata
            if ( name(j).eq.'SpectraDopplerCells' ) SDCS = cdata
         endif
         call mxDestroyArray(pa)
         call mxDestroyArray(pm)
20    continue

      ij = 1
      iomax = 1
      do 40 j = 1, j2
      do 40 i = 1, i2
         xc1 = X(ij)
         xc2 = Y(ij)
         xc3 = theta(j)
         xc4 = r(i)
         xc5 = VELO(ij)
         xc6 = EACC(ij)
c check for NaN's (actually check that it's not a NaN)
         if ( xc1.eq.xc1.and.xc2.eq.xc2.and.
     +        xc3.eq.xc3.and.xc4.eq.xc4.and. 
     +        xc5.eq.xc5.and.xc6.eq.xc6 ) then
            BEAR(iomax) = theta(j) 
            RNGE(iomax) = r(i)
            XLON(iomax) = X(ij)
            YLAT(iomax) = Y(ij)
            if ( r(i).ge.0.0 ) then
               XDST(iomax) = r(i) * sin ( theta(j) * pi / 180.0 )
               YDST(iomax) = r(i) * cos ( theta(j) * pi / 180.0 )
               VELU(iomax) = VELO(ij) * sin ( theta(j) * pi / 180.0 )
               VELV(iomax) = VELO(ij) * cos ( theta(j) * pi / 180.0 )
               VEL(iomax) = -1.0 * VELO(ij)
c convert velocities to cm/s
               VELU(iomax) = 100.0 * VELU(iomax)
               VELV(iomax) = 100.0 * VELV(iomax)
               VEL(iomax) = 100.0 * VEL(iomax)
               ACC(iomax) = 100.0 * EACC(ij)
               rad_bear = theta(j) * pi / 180.0
               rad_lat = slat * pi / 180.0
               rad_lon = slon * pi / 180.0
               rang = 1000.0 * r(i)
               call dirct1 ( rad_lat, rad_lon, yl, xl, 
     +                       rad_bear, rad_head, rang )
               HEAD(iomax) = real ( rad_head * 180.0 / pi )
               iomax = iomax + 1
            endif
         endif
         ij = ij + 1
40    continue
      iomax = iomax - 1

c  close file

      stat = matClose ( mfp )
      if (stat .ne. 0) then
         write(6,*) 'Error closing ''matdemo.mat''.'
         stop
      end if

c  check to see if there are no valid data; if not, don't open the file
      if ( iomax.eq.0 ) go to 999

c  open output file

      open ( 16, file = outfile, form = 'formatted' )

      print*,"writing to ", outfile
      write ( 16, 29 ) '%CTF: 1.00'
      write ( 16, 29 ) '%FileType: LLUV rdls'
      write ( 16, 59 ) '%Manufacturer: ', MANU
      write ( 16, 29 ) '%LLUVSpec: 1.00 2007 12 06'
      write ( 16, 59 ) '%Site: ', SITE
      write ( 16, 39 ) '%TimeStamp: ',iyear, imon, iday, ihr, imin, isec
      write ( 16, 29 ) '%TimeZone: "UTC" +0.00 0'
      write ( 16, 49 ) '%Origin:  ', slat, slon
      write ( 16, 29 ) 
     +  '%GreatCircle: "WGS84" 6378137.000 298.257223562997'
      write ( 16, 29 ) '%GeodVersion: "Vincenty (1979)" 2.0 2002 10 01'
      write ( 16, 59 ) '%RangeResolutionKMeters: ', RRKM
      write ( 16, 59 ) '%TransmitCenterFreqMHz: ', TCFM
      write ( 16, 59 ) '%DopplerResolutionHzPerBin: ', DRHB
      write ( 16, 59 ) '%CurrentVelocityLimit: ', CVLM
      write ( 16, 29 ) '%MergedCount: 0001'
      write ( 16, 59 ) '%TransmitSweepRateHz: ', TSRH
      write ( 16, 59 ) '%TransmitBandwidthKHz: ', TBKH
      write ( 16, 59 ) '%SpectraDopplerCells: ', SDCS
      write ( 16, 29 ) '%TableType: LLUV RDL1'
      write ( 16, 29 ) '%TableColumns: 11'
      write ( 16, 59 ) '%TableColumnTypes: LOND LATD VELU VELV EACC ',
     +                                'XDST YDST RNGE BEAR VELO HEAD'
      write ( 16, 69 ) '%TableRows:    ', iomax
      write ( 16, 29 ) '%TableStart:'
      write ( 16, 79 ) '%%   Longitude   Latitude    U comp   V comp',  
     +                 '  Accuracy  X Distance   Y Distance',
     +                 '     Range  Bearing  Velocity  Direction'
      write ( 16, 79 ) '%%     (deg)       (deg)     (cm/s)   (cm/s)',  
     +                 '   (cm/s)      (km)         (km)   ',
     +                 '      (km) (deg NCW)  (cm/s)   (deg NCW)'

      do 60 i = 1, iomax
         write ( 16, 89 ) XLON(i), YLAT(i), VELU(i), VELV(i), ACC(i), 
     +       XDST(i), YDST(i), RNGE(i), BEAR(i), VEL(i), HEAD(i)
60    continue

      write ( 16, 29 ) '%TableEnd:'
      write ( 16, 29 ) '%End:'

29    format ( a ) 
39    format ( a12, i4, 1x, i2.2, 1x, i2.2, 1x, i2.2, 1x, 
     +                                          i2.2, 1x, i2.2 )
49    format ( a10, f12.7, 2x, f12.7 )
59    format ( a, a )
69    format ( a, i5.5 )
79    format ( a44, a35, a40 )
89    format ( 1x, f12.7, 1x, f12.7, 1x, f8.3, 1x, f8.3, 1x, f8.3, 1x
     +        f12.4, 1x, f12.4, 1x, f9.3, 1x, f8.4, 1x, f9.3, 1x, f8.4 )

999   continue
      stop
      end

c  -----------------------------------------------------------------------------
      subroutine dirct1 ( glat1, glon1, glat2, glon2, bear, head, rang )
c  -----------------------------------------------------------------------------
c
c *** solution of the geodetic direct problem after t.vincenty
c *** modified rainsford's method with helmert's elliptical terms
c *** effective in any azimuth and at any distance short of antipodal
c
c *** a is the semi-major axis of the reference ellipsoid
c *** f is the flattening of the reference ellipsoid
c *** latitudes and longitudes in radians positive north and east
c *** azimuths in radians clockwise from north
c *** geodesic distance s assumed in units of semi-major axis a
c
c *** programmed for cdc-6600 by lcdr l.pfeifer ngs rockville md 20feb75
c *** modified for system 360 by john g gergen ngs rockville md 750608
c
c      implicit real*8 (a-h,o-z)
      real*8 head, glat2, glon2 
      real*8 c, cf, cu, cy, cz, c2a
      real*8 d, eps, e, f, r, sa, su, sy, tu, x, y
      real glat1, glon1, bear, rang
      data eps / 0.5d-13 /
      pi = 4.0 * atan ( 1.0 )
      rad = 180.0 / pi
c these are specific to GRS80 / WGS84  (NAD83) 
      a = 6378137.d0
      f = 1.d0 / 298.25722210088d0

      r = 1.-f
      tu = r * sin ( glat1 ) / cos ( glat1 )
      sf = sin ( bear )
      cf = cos ( bear )
      head = 0.0
      if ( cf.ne.0.0 ) head = atan2( tu, cf ) * 2.0
      cu = 1.0 / sqrt ( tu * tu + 1.0 )
      su = tu * cu
      sa = cu * sf
      c2a = -sa * sa + 1.0
      x = sqrt ( ( 1.0 / r / r - 1.0 ) * c2a + 1.0 ) + 1.0
      x = ( x - 2.0 ) / x
      c = 1.0 - x
      c = ( x * x / 4.0 + 1 ) / c
      d = ( 0.375d0 * x * x - 1.0 ) * x
      tu = rang / r / a / c
      y = tu
      iflag = 1
  100 sy = sin ( y )
      cy = cos ( y )
      cz = cos ( head + y )
      e = cz * cz * 2.0 - 1.0
      c = y
      x = e * cy
      y = e + e - 1.0
      y = ( ( ( sy * sy * 4.0 - 3.0 ) * y * cz * d / 6.0 + x ) * d / 4.0
     +                                  - cz ) * sy * d + tu
      iflag = iflag + 1
      if ( abs ( y - c ).gt.eps.and.iflag.lt.20 ) go to 100
      head = cu * cy * cf - su * sy
      c = r * sqrt ( sa * sa + head * head )
      d = su * cy + cu * sy * cf
      glat2 = atan2 ( d, c )
      c = cu * cy - su * sy * cf
      x = atan2 ( sy * sf, c )
      c = ( ( -3.0 * c2a + 4.0 ) * f + 4.0 ) * c2a * f / 16.0
      d = ( ( e * cy * c + cz ) * sy * c + y ) * sa
      glon2 = glon1 + x - ( 1.0 - c ) * d * f
      head = atan2 ( sa, head ) + pi
      return
      end
