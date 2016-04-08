MODULE trcext_blingv0

#if defined key_bling

   !!======================================================================
   !!                     ***  MODULE trcext_bling  ***
   !! TOP :  Iron external inputs
   !!======================================================================
   !! History :  MClaret@McGill@04/2016. CO2 external fluxes added
   !! History :  MClaret@McGill@04-07/2014
   !!======================================================================

   USE oce_trc
   USE par_trc
   USE trc
   USE fldread 
   USE iom

   USE vars_bling

   IMPLICIT NONE
   PRIVATE

   PUBLIC trc_ext_bling
   PUBLIC trc_ext_init_bling

   INTEGER :: numphpext, numfedext, numoxyext

   TYPE(FLD), ALLOCATABLE, DIMENSION(:) :: sf_patm_bling ! structure of input fields (file informations, fields read)
   TYPE(FLD), ALLOCATABLE, DIMENSION(:) :: sf_dust_bling

# include "top_substitute.h90"

CONTAINS

   SUBROUTINE trc_ext_bling (kt) 

      !-----------------------------------------------------------------------
      ! The prefixes "b" refers to a bottom flux and "s" to a sea surface flux
      !-----------------------------------------------------------------------

      INTEGER, INTENT( IN ) ::   kt   ! ocean time step

      INTEGER  :: ji, jj, jk, ikb

      REAL(wp) :: zrfact
      REAL(wp) :: sss, sst, o2fact
      REAL(wp) :: tt, tk, ts, ts2, ts3, ts4, ts5

      REAL(wp) :: ztc, ztc2, ztc3, zws, zkgwan, xconv
      REAL(wp) :: zsch_co2, zsch_o2
      REAL(wp) :: co2_sat, co2_sur, tmp1, tmp2
      REAL(wp) :: o2_alpha, o2_alpha2, o2_sat, o2_sur, pco2_atm

      REAL(wp) :: foxy, fe_2_p_sed
      REAL(wp) :: sumbpo4_glob, sumtffed_glob, sumbfed_glob, sumtfoxy_glob, sumboxy_glob

      REAL(wp), POINTER, DIMENSION(:,:) :: tmask_ikb
      REAL(wp), POINTER, DIMENSION(:,:) :: sch_no_term_co2, co2_flx
      REAL(wp), POINTER, DIMENSION(:,:) :: sch_no_term_o2 , o2_flx
      REAL(wp), POINTER, DIMENSION(:,:) :: bpo4, bfed, boxy, bdic, balk

      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) ' trc_ext_bling:  BLINGv0 model'
      IF(lwp) WRITE(numout,*) ' ~~~~~~~~~~~~~~'

      CALL wrk_alloc( jpi, jpj, tmask_ikb )
      CALL wrk_alloc( jpi, jpj, sch_no_term_co2, co2_flx )
      CALL wrk_alloc( jpi, jpj, sch_no_term_o2 , o2_flx  )
      CALL wrk_alloc( jpi, jpj, bpo4, bfed, boxy, bdic, balk )

      !---------------------------------------------------------------------
      ! Calculate air-sea exhange for O2/CO2
      !---------------------------------------------------------------------

      ! Get atmospheric pressure if read from file
      IF( ln_patm_bling ) THEN
         CALL fld_read( kt, 1, sf_patm_bling )
         patm_bling(:,:) = sf_patm_bling(1)%fnow(:,:,1)
      ENDIF

      ! Get pco2 atmospheric value
      IF ( (nyear .GT. 1958) .AND. (nyear .LT. 2016) ) THEN
        pco2_atm=pco2(nyear-1958)
      ELSE
        pco2_atm=pco2_cte  
      ENDIF

      xconv  = 0.01_wp / 3600._wp ! convert from cm/h to m/s (piston velocity)
      o2fact = 1.d0 /22.3916d0 ! convert from ml/L to mol/m3

      !---------------------------------------------------------------------
      ! Compute gas exchange coefficients (wind and ice effects considered)
      ! Ported exactly from PISCES scheme p4zflx
!CDIR NOVERRCHK
      DO jj = 1, jpj
!CDIR NOVERRCHK
         DO ji = 1, jpi
            ztc  = MIN( 35., tsn(ji,jj,1,jp_tem) )
            ztc2 = ztc * ztc
            ztc3 = ztc * ztc2 

            ! Compute the schmidt Number both O2 and CO2
            zsch_co2 = 2073.1 - 125.62 * ztc + 3.6276 * ztc2 - 0.043126 * ztc3
            zsch_o2  = 1953.4 - 128.0  * ztc + 3.9918 * ztc2 - 0.050091 * ztc3

            !  wind speed 
            zws  = wndm(ji,jj) * wndm(ji,jj)
            zkgwan = 0.3 * zws  + 2.5 * ( 0.5246 + 0.016256 * ztc + 0.00049946  * ztc2 )
            zkgwan = zkgwan * xconv * ( 1.d0 - fr_i(ji,jj) ) * tmask(ji,jj,1)

            ! air-sea gas transfer velocity (piston velocity). Units are m/s
            sch_no_term_co2(ji,jj) = zkgwan * SQRT( 660.d0 / zsch_co2 )
            sch_no_term_o2 (ji,jj) = zkgwan * SQRT( 660.d0 / zsch_o2 )

         END DO
      END DO

      ! Gas air-sea fluxes
      DO jj=1,jpj
        DO ji=1,jpi

          zrfact = rfact / fse3t(ji,jj,1) ! [s/m]

          !-----------------
          ! CO2 air-sea flux
          ! saturation concentration (solubility from co2calc subroutine) * piston v [mol/m3 m/s]
          co2_sat=patm_bling(ji,jj)*pco2_atm*1.e-6*co2_alpha(ji,jj)*sch_no_term_co2(ji,jj) 

          ! surface concentration (computed in co2calc subroutine) * piston v [mol/m3 m/s]
          co2_sur=co2_csurf(ji,jj)*sch_no_term_co2(ji,jj)

          ! air-sea flux [mol/m3 m/s]
          co2_flx(ji,jj)=co2_sat-co2_sur

          ! add trend
          trn(ji,jj,1,jpDIC_bling)=trn(ji,jj,1,jpDIC_bling)+co2_flx(ji,jj)*zrfact

          !----------------
          ! O2 air-sea flux
          ! Compute O2 solubility (GARCIA and GORDON 1992, eq.8)
          !sst=MIN( 35.,tsn(ji,jj,1,jp_tem) )
          sst=tsn(ji,jj,1,jp_tem)
          sss=tsn(ji,jj,1,jp_sal)

          tt=298.15d0-sst
          tk=273.15d0+sst
          ts=LOG(tt/tk)
          ts2=ts *ts
          ts3=ts2*ts
          ts4=ts3*ts
          ts5=ts4*ts

          ! [mol/m3]
          tmp1=a_0+a_1*ts+a_2*ts2+a_3*ts3+a_4*ts4+a_5*ts5
          tmp2=sss*(b_0+b_1*ts+b_2*ts2+b_3*ts3)+c_0*sss*sss
          o2_alpha=o2fact*EXP(tmp1+tmp2)

          ! [mol/m3/atm]
          o2_alpha2=o2_alpha/po2_atm

          ! saturation concentration * piston v [mol/m2/s]
          !o2_sat=patm_bling(ji,jj)*po2_atm*o2_alpha2*sch_no_term_o2(ji,jj)
          o2_sat=po2_atm*o2_alpha2*sch_no_term_o2(ji,jj)

          ! surface concentration * piston v [mol/m2/s]
          o2_sur=trn(ji,jj,1,jpOxy_bling)*sch_no_term_o2(ji,jj)

          ! Flux of oxygen [mol/m2/s]
          o2_flx(ji,jj)=o2_sat-o2_sur

          ! add trend
          trn(ji,jj,1,jpOxy_bling)=trn(ji,jj,1,jpOxy_bling)+o2_flx(ji,jj)*zrfact

          !sumtfoxy=sumtfoxy+o2flx(ji,jj)*zrfact*cvol(ji,jj,1)*tmask(ji,jj,1)

        ENDDO
      ENDDO

      !---------------------------------------------------------------------
      ! Calculate external fluxes for iron. 
      !---------------------------------------------------------------------

      !Get dust field [mol Fe/m2/s]
      IF( ln_dust_bling ) THEN
         CALL fld_read( kt, 1, sf_dust_bling )
         dust_bling(:,:) = sf_dust_bling(1)%fnow(:,:,1)
      ENDIF

      fe_2_p_sed=106.0e-4

      DO jj = 1, jpj
         DO ji = 1, jpi
            ! [s/m] -> dust deposition in mol Fe/m3
            zrfact  = rfact / fse3t(ji,jj,1) 

            ! [mol Fe/m3]
            trn(ji,jj,1,jpFed_bling) =  trn(ji,jj,1,jpFed_bling) + dust_bling(ji,jj)*zrfact

            !sumtffed=sumtffed+dust_bling(ji,jj)/fse3t(ji,jj,1)*cvol(ji,jj,1)*tmask(ji,jj,1)

         ENDDO
      ENDDO

      ! -------------------------------------
      ! Exchange with sediments
      ! Calculate external bottom fluxes for tracer_vertdiff. Positive fluxes
      ! are from the water column into the seafloor. For P, the bottom flux  
      ! puts the sinking flux reaching the bottom cell into the water column 
      ! through diffusion. For iron, the sinking flux disappears into the 
      ! sediments if bottom waters are oxic (assumed adsorbed as oxides),
      ! while an efflux of dissolved iron occurs dependent on the supply of
      ! reducing organic matter (scaled by the org-P sedimentation rate).
      ! If bottom waters are anoxic, the sinking flux of Fe is returned to
      ! the water column. Note this is not appropriate for very long runs
      ! with an anoxic ocean (iron will keep accumulating forever).
      ! For oxygen, the consumption of oxidant required to respire  
      ! the settling flux of organic matter (in support of the
      ! PO4 bottom flux) diffuses from the bottom water into the sediment.
      ! Do not bury any C - all goes back to water column (for DIC).
      ! Return all CaCO3 to the water column to preserve the alkalinity inventory
      ! Note that the flux of NO3 out of the sediment, inferred from 
      ! the PO4 flux, causes a negative flux of alkalinity.
      ! -------------------------------------

      DO jj = 1, jpj
        DO ji = 1, jpi

            ! mbkt is a matrix containing the vertical index of the
            ! bottom layer at each horizontal point
            ikb    = mbkt(ji,jj)
	    tmask_ikb(ji,jj)=tmask(ji,jj,ikb)

            foxy  = trn(ji,jj,ikb,jpOxy_bling)
            zrfact = rfact / fse3t(ji,jj,ikb) 

            ! [mol P/m2/s]
            bpo4(ji,jj)=fpop_b(ji,jj)

            ! [mol O2/m2/s]
            boxy(ji,jj) = -oxy2p*fpop_b(ji,jj)

            ! [mol Fe/m2/s]
            IF (foxy>oxy_min) THEN
               bfed(ji,jj)= fe_2_p_sed*fpop_b(ji,jj)
            ELSE
               bfed(ji,jj)=(fe_2_p_sed*fpop_b(ji,jj)+fpofe_b(ji,jj))
            ENDIF

            ! [mol C/m2/s]
            bdic(ji,jj)=fpop_b(ji,jj)*c2p+fcaco3_b(ji,jj)

            ! [mol eq/m2/s]
            balk(ji,jj)=2.d0*fcaco3_b(ji,jj)-fpop_b(ji,jj)*n2p

            ! Add the bottom flux trend
            ! [mol/m3]
            trn(ji,jj,ikb,jpPO4_bling) = trn(ji,jj,ikb,jpPO4_bling) + bpo4(ji,jj)*zrfact
            trn(ji,jj,ikb,jpFed_bling) = trn(ji,jj,ikb,jpFed_bling) + bfed(ji,jj)*zrfact
            trn(ji,jj,ikb,jpOxy_bling) = trn(ji,jj,ikb,jpOxy_bling) + boxy(ji,jj)*zrfact
            trn(ji,jj,ikb,jpDIC_bling) = trn(ji,jj,ikb,jpDIC_bling) + bdic(ji,jj)*zrfact
            trn(ji,jj,ikb,jpalk_bling) = trn(ji,jj,ikb,jpalk_bling) + balk(ji,jj)*zrfact

            bpo4(ji,jj)=bpo4(ji,jj)*tmask(ji,jj,ikb)
            bfed(ji,jj)=bfed(ji,jj)*tmask(ji,jj,ikb)
            boxy(ji,jj)=boxy(ji,jj)*tmask(ji,jj,ikb)

         END DO
      END DO

      IF (ln_bling_ext) THEN

        !IF (narea==1) THEN

          IF( kt == nittrc000 ) THEN
            CALL ctl_opn( numphpext, 'php.extflx' , 'REPLACE', 'FORMATTED', 'SEQUENTIAL', -1, 6, .FALSE., narea )
            CALL ctl_opn( numfedext, 'fed.extflx' , 'REPLACE', 'FORMATTED', 'SEQUENTIAL', -1, 6, .FALSE., narea )
            CALL ctl_opn( numoxyext, 'oxy.extflx' , 'REPLACE', 'FORMATTED', 'SEQUENTIAL', -1, 6, .FALSE., narea )
          ENDIF

          WRITE(numphpext,9500) kt,  glob_sum( bpo4(:,:)*e1e2t(:,:) )
          WRITE(numfedext,9501) kt,  glob_sum( bfed(:,:)*e1e2t(:,:) ), glob_sum( dust_bling(:,:)*e1e2t(:,:)*tmask(:,:,1) )
          WRITE(numoxyext,9501) kt,  glob_sum( boxy(:,:)*e1e2t(:,:) ), glob_sum(     o2_flx(:,:)*e1e2t(:,:)*tmask(:,:,1) )

        !ENDIF

9500  FORMAT(i10,e18.10)
9501  FORMAT(i10,2(e18.10))

      ENDIF

      IF (ln_diatrc) THEN

        IF (lk_iomput) THEN
          CALL iom_put( "po4_btf", bpo4      (:,:)*tmask_ikb(:,:) )
          CALL iom_put( "oxy_stf", o2_flx    (:,:)*tmask_ikb(:,:) )
          CALL iom_put( "oxy_btf", boxy      (:,:)*tmask_ikb(:,:) )
          CALL iom_put( "fed_stf", dust_bling(:,:)*tmask_ikb(:,:) )
          CALL iom_put( "fed_btf", bfed      (:,:)*tmask_ikb(:,:) )
          CALL iom_put( "fed_bur", fpofe_b   (:,:)*tmask_ikb(:,:) )
          CALL iom_put( "dic_stf", co2_flx   (:,:)*tmask_ikb(:,:) )
          CALL iom_put( "dic_btf", bdic      (:,:)*tmask_ikb(:,:) )
          CALL iom_put( "alk_btf", balk      (:,:)*tmask_ikb(:,:) )
        ENDIF

        IF (  .NOT. lk_iomput ) THEN
          trc2d( :,:,jp_bling0_2d+1 ) = bpo4      (:,:)*tmask_ikb(:,:)
          trc2d( :,:,jp_bling0_2d+2 ) = o2_flx    (:,:)*tmask_ikb(:,:)
          trc2d( :,:,jp_bling0_2d+3 ) = boxy      (:,:)*tmask_ikb(:,:)
          trc2d( :,:,jp_bling0_2d+4 ) = dust_bling(:,:)*tmask_ikb(:,:)
          trc2d( :,:,jp_bling0_2d+5 ) = bfed      (:,:)*tmask_ikb(:,:)
          trc2d( :,:,jp_bling0_2d+6 ) = fpofe_b   (:,:)*tmask_ikb(:,:)
          trc2d( :,:,jp_bling0_2d+7 ) = co2_flx   (:,:)*tmask_ikb(:,:)
          trc2d( :,:,jp_bling0_2d+8 ) = bdic      (:,:)*tmask_ikb(:,:)
          trc2d( :,:,jp_bling0_2d+9 ) = balk      (:,:)*tmask_ikb(:,:)
        ENDIF

      ENDIF

      CALL wrk_dealloc( jpi, jpj, tmask_ikb )
      CALL wrk_dealloc( jpi, jpj, sch_no_term_co2, co2_flx )
      CALL wrk_dealloc( jpi, jpj, sch_no_term_o2 , o2_flx  )
      CALL wrk_dealloc( jpi, jpj, bpo4, bfed, boxy, bdic, balk )

   END SUBROUTINE trc_ext_bling

   SUBROUTINE trc_ext_init_bling

      INTEGER  :: jm, ierr
      INTEGER  :: numdust, ntimes_dust
      INTEGER  :: numpatm, ntimes_patm

      REAL(wp), DIMENSION(12) :: zsteps                 ! times records

      REAL(wp) , ALLOCATABLE, DIMENSION(:,:,:) :: zdust, zpatm
      !!----------------------------------------------------------------------

      ALLOCATE( sf_dust_bling(1), STAT=ierr )    
      IF( ierr > 0 )   CALL ctl_stop( 'STOP', 'trc_ext_init_bling: unable to allocate sf_apr structure' )

      CALL fld_fill( sf_dust_bling, (/ sn_dust_bling /), cn_dir_bling, 'trc_ext_init_bling', 'Iron from sediment ', 'namblingiron' )

      ALLOCATE( sf_dust_bling(1)%fnow(jpi,jpj,1) )
      IF( sn_dust_bling%ln_tint ) ALLOCATE( sf_dust_bling(1)%fdta(jpi,jpj,1,2) )

      CALL iom_open (  TRIM( sn_dust_bling%clname ) , numdust )

      CALL iom_gettime( numdust, zsteps, kntime=ntimes_dust)  ! get number of record in file

      ALLOCATE( zdust(jpi,jpj,ntimes_dust) )
      DO jm = 1, ntimes_dust
         CALL iom_get( numdust, jpdom_data, TRIM( sn_dust_bling%clvar ), zdust(:,:,jm), jm )
      ENDDO
       
      CALL iom_close( numdust )
      DEALLOCATE( zdust)

      !!-------------
      
      IF( ln_patm_bling ) THEN

        ALLOCATE( sf_patm_bling(1), STAT=ierr )       
        IF( ierr > 0 )   CALL ctl_stop( 'STOP', 'p4z_flx: unable to allocate sf_patm structure' )

        CALL fld_fill( sf_patm_bling, (/ sn_patm_bling /), cn_dir_patm_bling, 'trc_ext_init_bling', 'Atmospheric pressure ', 'namblingpatm' )

        ALLOCATE( sf_patm_bling(1)%fnow(jpi,jpj,1)   )
        IF( sn_patm_bling%ln_tint )  ALLOCATE( sf_patm_bling(1)%fdta(jpi,jpj,1,2) )

        CALL iom_open (  TRIM( sn_patm_bling%clname ) , numpatm )

        CALL iom_gettime( numpatm, zsteps, kntime=ntimes_patm)  ! get number of record in file

        ALLOCATE( zpatm(jpi,jpj,ntimes_patm) )

        DO jm = 1, ntimes_patm
          CALL iom_get( numpatm, jpdom_data, TRIM( sn_patm_bling%clvar ), zpatm(:,:,jm), jm )
        ENDDO
       
        CALL iom_close( numpatm )
        DEALLOCATE( zpatm)

      ELSE

        patm_bling(:,:)=1.e0;

      ENDIF

   END SUBROUTINE trc_ext_init_bling

#else
   !!----------------------------------------------------------------------
   !!   Dummy module                                        No BLINGv0 model
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_ext_bling 
      WRITE(*,*) 'trc_ext_bling: You should not have seen this print',  kt
   END SUBROUTINE trc_ext_bling

#endif

END MODULE trcext_blingv0
