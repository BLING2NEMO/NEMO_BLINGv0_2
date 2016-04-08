MODULE trcini_bling
   !!======================================================================
   !!                         ***  MODULE trcini_bling  ***
   !! TOP :   initialisation of the BLINGv0 tracers
   !!======================================================================
   !! History :   2.0  !  2007-12  (C. Ethe, G. Madec) Original code
   !!----------------------------------------------------------------------
#if defined key_bling
   !!----------------------------------------------------------------------
   !!   'key_bling'                                               CFC tracers
   !!----------------------------------------------------------------------
   !! trc_ini_bling   : BLINGv0 model initialisation
   !!----------------------------------------------------------------------
   USE par_trc         ! TOP parameters
   USE oce_trc
   USE trc
   USE trcsms_blingv0
   USE trcopt_blingv0
   USE trcext_blingv0

   !: PISCES indices for PO4, DOC, O2, F to initialize BLING fields
   !USE par_pisces , ONLY : jppo4, jpdoc, jpoxy, jpfer 

   ! BLINGv0 modules
   USE vars_bling

   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_ini_bling   ! called by trcini.F90 module

   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: trcini_bling.F90 -1   $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE trc_ini_bling
      !!----------------------------------------------------------------------
      !!                     ***  trc_ini_bling  ***  
      !!
      !! ** Purpose :   initialization for BLINGv0 model
      !!
      !! ** Method  : - Read the namcfc namelist and check the parameter values
      !!----------------------------------------------------------------------

      !                       ! Allocate BLINGv0 arrays

      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) ' trc_ini_bling: initialisation of BLINGv0 model'
      IF(lwp) WRITE(numout,*) ' ~~~~~~~~~~~~~~'

      ! Allocate BLINGv0 arrays
      IF( bling_alloc() /= 0 )   CALL ctl_stop( 'STOP', 'trc_ini_bling: unable to allocate BLINGv0 arrays' )
      
      rfact   = rdttrc(1)                          

      IF(lwp) WRITE(numout,*) '    Passive Tracer time step    rfact  = ', rfact, ' rdt = ', rdttra(1)

      IF( .NOT. ln_rsttr ) THEN
        ! initialize prognostic tracers
        !trn(:,:,:,jpPO4_bling) = trn(:,:,:,jppo4) * tmask(:,:,:) / 106.d0 * 1.d3
        !trn(:,:,:,jpDOP_bling) = trn(:,:,:,jpdoc) * tmask(:,:,:) / 106.d0 * 1.d3
        !trn(:,:,:,jpFed_bling) = trn(:,:,:,jpfer) * tmask(:,:,:) * 1.d3
        !trn(:,:,:,jpOxy_bling) = trn(:,:,:,jpoxy) * tmask(:,:,:) * 1.d3
        !trn(:,:,:,jpine_bling) = trn(:,:,:,jpPO4_bling) + trn(:,:,:,jpDOP_bling)

        trn(:,:,:,jpPO4_bling) = 2.174e-6_wp * tmask(:,:,:) 
        trn(:,:,:,jpDOP_bling) = 1.000e-8_wp * tmask(:,:,:) 
        trn(:,:,:,jpFed_bling) = 0.6E-9      * tmask(:,:,:) 
        trn(:,:,:,jpOxy_bling) = 177.6e-6_wp * tmask(:,:,:) 
        trn(:,:,:,jpDIC_bling) = 0.E-9       * tmask(:,:,:)
        trn(:,:,:,jpalk_bling) = 0.E-9       * tmask(:,:,:)
      ENDIF
         
      ! initialize diagnostic tracers
      fpop_b   (:,:)   = 0.d0
      fpofe_b  (:,:)   = 0.d0
      co2_csurf(:,:)   = 0.d0
      co2_alpha(:,:)   = 0.d0
      chl_bling(:,:,:) = 0.d0
      irr_mem  (:,:,:) = 0.d0
      co3_ion  (:,:,:) = 0.d0
      htotal   (:,:,:) = 1.e-8

      ! initialize cumulative sums
      sumtffed = 0.d0
      sumtfoxy = 0.d0
      sumbfpo4 = 0.d0
      sumbffed = 0.d0
      sumbfoxy = 0.d0

      CALL trc_opt_bling_init
      CALL trc_ext_init_bling

      ! Data from Mauna Loa station
      pco2= (/ 315.97,316.91,317.64,318.45,318.99,319.62,320.04,321.38,322.16,323.04 &
              ,324.62,325.68,326.32,327.45,329.68,330.18,331.08,332.05,333.78,335.41 &
              ,336.78,338.68,340.10,341.44,343.03,344.58,346.04,347.39,349.16,351.56 &
              ,353.07,354.35,355.57,356.38,357.07,358.82,360.80,362.59,363.71,366.65 &
              ,368.33,369.52,371.13,373.22,375.77,377.49,379.80,381.90,383.76,385.59 &
              ,387.37,389.85,391.63,393.82,396.48,398.61,400.83 /)

      !INTEGER :: nypco2=(/ 1959 1960 1961 1962 1963 1964 1965 1966 1967 1968 &
      !                    ,1969 1970 1971 1972 1973 1974 1975 1976 1977 1978 &
      !                    ,1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 &
      !                    ,1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 &
      !                    ,1999 2000 2001 2002 2003 2004 2005 2006 2007 2008 &
      !                    ,2009 2010 2011 2012 2013 2014 2015 /)

   END SUBROUTINE trc_ini_bling

#else
   !!----------------------------------------------------------------------
   !!   Dummy module                                        No BLINGv0 model
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_ini_bling             ! Empty routine
   END SUBROUTINE trc_ini_bling
#endif
   !!======================================================================
END MODULE trcini_bling
