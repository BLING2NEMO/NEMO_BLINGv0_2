MODULE trcsms_blingv0

#if defined key_bling
   !!======================================================================
   !!                     ***  MODULE trcsms_bling  ***
   !! TOP :  BLINGv0 model main routine. Computes the sources and sinks
   !!======================================================================
   !! History :  MClaret@McGill@04/2016. Carbon cycle included
   !! History :  MClaret@McGill@04-07/2014. Ported from BLINGv0 in GFDL
   !!-----------------------------------------------------------------------
   USE par_trc         ! TOP parameters
   USE oce_trc         ! Ocean variables
   USE trc             ! TOP variables
   USE trdmod_oce
   USE trdmod_trc
   USE iom

   ! BLINGv0 specific modules
   USE trcopt_blingv0
   USE trcext_blingv0
   USE vars_bling
   USE FMS_ocmip2_co2calc_mod

   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_sms_bling       ! called by trcsms.F90 module

   INTEGER :: numsms, numsms2, numphp, numfed, numoxy, numcar ! logical unit for phosphate budget

#  include "top_substitute.h90"

CONTAINS

   SUBROUTINE trc_sms_bling( kt )
      !!-------------------------------------------------------------------------------
      !!                     ***  trc_sms_bling  ***
      !!
      !! ** History : default leap-frog scheme(MY_TRC) changed to FORWARD scheme
      !!-------------------------------------------------------------------------------
      !! The prefixes "f" refers to a "field" and "j" to a volumetric rate
      !! j is followed by the currency unit (p for phosphate, ca for calcium carbonate)
      !!-------------------------------------------------------------------------------
      
      INTEGER, INTENT(in) ::   kt   ! ocean time-step index

      INTEGER  :: ji, jj, jk, jn

      REAL(wp) :: fpo4, fdop, ffed, foxy
      REAL(wp) :: ztra
      ! Irradiance k
      REAL(wp) :: po4_up, thetamax_fe, alpha_chl 
      ! Production
      REAL(wp) :: pc_tot, biomass_p_ts, theta, chl_dia, mulamb0expkT
      ! Phosphorous
      REAL(wp) :: frac_pop, jp_dop, fe2p_up, fpopkm1
      REAL(wp) :: zzz, wsink, oxy_up
      ! Iron
      REAL(wp) :: jfe_pofe, fpofekm1
      REAL(wp) :: dum5, dum2, dum3
!DO_CARBON
      ! Carbonate system
      REAL(wp) :: s_over_p, co3_solubility, fcaco3km1
      REAL(wp), POINTER, DIMENSION(:,:,:) :: jca_uptake, jca_reminp, fcaco3
      REAL(wp), POINTER, DIMENSION(:,:,:) :: jdic, jalk
      REAL(wp), POINTER, DIMENSION(:,:,:) :: zremin_caco3, frac_lg
      REAL(wp), POINTER, DIMENSION(:,:)   :: pco2_surf
!DO_CARBON

      REAL(wp), POINTER, DIMENSION(:,:,:) :: rho
      REAL(wp), POINTER, DIMENSION(:,:,:) :: expkT, irr_inst, irr_mix, irrk, pc_m, mu
      REAL(wp), POINTER, DIMENSION(:,:,:) :: def_fe, feprime, kfe_eq_lig, fpofe
      REAL(wp), POINTER, DIMENSION(:,:,:) :: jp_pop, fpop, zremin
      REAL(wp), POINTER, DIMENSION(:,:,:) :: jp_uptake, jp_remin, jp_recycle
      REAL(wp), POINTER, DIMENSION(:,:,:) :: jfe_uptake, jfe_remin, jfe_recycle
      REAL(wp), POINTER, DIMENSION(:,:,:) :: jfe_ads_inorg, jfe_ads_org
      REAL(wp), POINTER, DIMENSION(:,:,:) :: jpo4, jdop, jfed, joxy

      REAL(wp), POINTER, DIMENSION(:)     :: dum4
      REAL(wp), POINTER, DIMENSION(:,:,:) :: xnegtr

      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('trc_sms_bling')
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) ' trc_sms_bling:  BLINGv0 model'
      IF(lwp) WRITE(numout,*) ' ~~~~~~~~~~~~~~'

      ! allocate matrix variables
      CALL wrk_alloc( jpi, jpj, jpk, rho )
      CALL wrk_alloc( jpi, jpj, jpk, expkT, irr_inst, irr_mix, irrk, pc_m, mu )
      CALL wrk_alloc( jpi, jpj, jpk, def_fe, feprime, kfe_eq_lig, fpofe )
      CALL wrk_alloc( jpi, jpj, jpk, jp_pop, fpop, zremin )
      CALL wrk_alloc( jpi, jpj, jpk, jp_uptake,  jp_remin,  jp_recycle )
      CALL wrk_alloc( jpi, jpj, jpk, jfe_uptake, jfe_remin, jfe_recycle )
      CALL wrk_alloc( jpi, jpj, jpk, jfe_ads_inorg, jfe_ads_org )
      CALL wrk_alloc( jpi, jpj, jpk, jpo4, jdop, jfed, joxy )
      !DO_CARBON
      CALL wrk_alloc( jpi, jpj, jpk, jca_uptake, jca_reminp, fcaco3 )
      CALL wrk_alloc( jpi, jpj, jpk, jdic, jalk )
      CALL wrk_alloc( jpi, jpj, jpk, zremin_caco3, frac_lg )
      CALL wrk_alloc( jpi, jpj     , pco2_surf )

      CALL wrk_alloc( jpk, dum4 )
      CALL wrk_alloc( jpi, jpj, jpk, xnegtr )

      !IF( (narea==1) .and. ln_bling_mass ) THEN      !   Write values for phosphate budget
      IF( ln_bling_mass ) THEN      !   Write values for phosphate budget
        CALL trc_sms_bling_mass_conserv (kt)
      ENDIF

      rho(:,:,:)=rhd(:,:,:)*1035.e0+1035.e0

      jk=1
      CALL FMS_ocmip2_co2calc( 1,jpi,1,jpj,1,jpi,1,jpj,tmask(:,:,jk)-tmask(:,:,jk)+1.e0                &
                              ,tsn(:,:,jk,jp_tem),tsn(:,:,jk,jp_sal)                                   &
                              ,trn(:,:,jk,jpDIC_bling)/rho(:,:,jk),trn(:,:,jk,jpPO4_bling)/rho(:,:,jk) &
                              ,trn(:,:,jk,jpPO4_bling)/rho(:,:,jk),trn(:,:,jk,jpalk_bling)/rho(:,:,jk) &
                              ,htotal(:,:,jk)*htotal_scale_lo                                          &
                              ,htotal(:,:,jk)*htotal_scale_hi                                          &
                              ,htotal(:,:,jk)                                                          & 
                              ,co3_ion=co3_ion(:,:,jk)                                                 &
                              ,co2star=co2_csurf,alpha=co2_alpha                                       &
                              ,pCO2surf=pco2_surf                                                       )

      ! convert from mol/kg to mol/m3
      co2_alpha(:,:)  =co2_alpha(:,:)  *rho(:,:,1)
      co2_csurf(:,:)  =co2_csurf(:,:)  *rho(:,:,1)

      DO jk=2, jpk

          CALL FMS_ocmip2_co2calc( 1,jpi,1,jpj,1,jpi,1,jpj,tmask(:,:,jk)-tmask(:,:,jk)+1.e0                &
                                  ,tsn(:,:,jk,jp_tem),tsn(:,:,jk,jp_sal)                                   &
                                  ,trn(:,:,jk,jpDIC_bling)/rho(:,:,jk),trn(:,:,jk,jpPO4_bling)/rho(:,:,jk) &
                                  ,trn(:,:,jk,jpPO4_bling)/rho(:,:,jk),trn(:,:,jk,jpalk_bling)/rho(:,:,jk) &
                                  ,htotal(:,:,jk)*htotal_scale_lo                              &
                                  ,htotal(:,:,jk)*htotal_scale_hi                              &
                                  ,htotal(:,:,jk),co3_ion=co3_ion(:,:,jk)                       )

          ! MC note 04/02/16: I do not use the mask. NEMO advects fields w/o mask,
          ! and I'm not sure how masking co3_ion and htotal could affect DIC/ALK fields.

      ENDDO

      ! NOTE: htotal and co3_ion are left to mol/kg deliberately. 
      ! See zremin_caco3 computation for comment on co3_ion units
      !----------------------------------------------------------------------
      ! BLING model

      CALL trc_opt_bling_rb  (kt, irr_inst, irr_mix)  ! optical model (red and blue wavelengths)
      !CALL trc_opt_bling_rgb (kt, irr_inst, irr_mix)  ! optical model (red, blue, and green wavelengths)

      DO ji=1, jpi
          DO jj=1, jpj

              dum4(:) = 0.d0

              DO jk=1, jpk

               ! ----------------------------------------------------------
               ! negative trophic variables DO not contribute to the fluxes
               ! ----------------------------------------------------------

               ! [mol/m3]
               fpo4 = MAX( 0.e0, trn(ji,jj,jk,jpPO4_bling) )
               fdop = MAX( 0.e0, trn(ji,jj,jk,jpDOP_bling) )
               ffed = MAX( 0.e0, trn(ji,jj,jk,jpFed_bling) )
               foxy = trn(ji,jj,jk,jpOxy_bling)

               ! ----------------------------------------------------------
               ! TEMPERATURE DEPENDENCE
               ! NB: The temperature effect of Eppley (1972) is used instead 
               !     of that in Geider et al (1997) for both simplicity and 
               !     to incorporate combined effects on uptake, incorporation
               !     into organic matter and photorespiration.  Values of PCmax
               !     are normalized to 0C rather than 20C in Geider et al.(1997)
               ! ----------------------------------------------------------

               ! [no units]
               expkT(ji,jj,jk)=EXP(kappa_eppley*tsn(ji,jj,jk,jp_tem))

               ! ----------------------------------------------------------
               ! Phytoplankton are assumed to grow according to the general properties 
               ! described in Geider (1997). This formulation gives a biomass-specific 
               ! growthrate as a function of light, nutrient limitation, and 
               ! temperature. We modify this relationship slightly here, as described 
               ! below, and also use the assumption of steady state growth vs. loss to 
               ! derive a simple relationship between growth rate, biomass and uptake.
               ! ----------------------------------------------------------
               ! First, we calculate the limitation terms for PO4 and Fe, and the 
               ! Fe-limited Chl:C maximum.
               ! The light-saturated maximal photosynthesis rate term (pc_m) is simply 
               ! the product of a prescribed maximal photosynthesis rate (pc_0), the 
               ! Eppley temperature dependence, and a Liebig limitation (the minimum
               ! of Michaelis-Menton PO4-limitation, or iron-limitation).
               ! The iron limitation term has a lower limit of def_fe_min 
               ! and is scaled by (k_fe_2_p + fe_2_p_max) / fe_2_p_max
               ! so that it approaches 1 as fed approaches infinity. Thus, 
               ! it's of comparable magnitude to the PO4 limitation term.
               !
               ! Fe limitation acts in two additional mechanisms:
               ! 1. By reducing the maximum achievable Chl:C ratio 
               ! (theta) below a prescribed, Fe-replete maximum value (thetamax), to 
               ! approach a prescribed minimum Chl:C (thetamin) under extreme
               ! Fe-limitation.
               ! 2. By reducing alpha (the initial slope of the P-I curve) under Fe-
               ! limitation.
               ! ----------------------------------------------------------

               ! Iron uptake [no units]
               fe2p_up = fe2p_max * ffed / (kfe + ffed)  ![mol Fe/mol P]
               def_fe(ji,jj,jk)  = def_fe_min + (1.d0-def_fe_min) &
                                  *fe2p_up/(kfe2p_up+fe2p_up)*(kfe2p_up+fe2p_max)/fe2p_max 

               ! Phosphate uptake [no units]
               po4_up = fpo4 /( kpo4 + fpo4 )

               ! Maximum production (units of pc_0 [s-1])
               pc_m(ji,jj,jk) = pc_0 * expkT(ji,jj,jk) * MIN(po4_up,def_fe(ji,jj,jk)) 

               ! Iron limitation on photosyntesis machinery
               thetamax_fe=thetamax_lo + (thetamax_hi - thetamax_lo)*def_fe(ji,jj,jk)
               alpha_chl  =alpha_min   + (alpha_max   - alpha_min  )*def_fe(ji,jj,jk)

               !-----------------------------------------------------------------------
               ! Next, the nutrient-limited efficiency of algal photosystems, Irrk, is
               ! calculated. This requires a prescribed quantum yield, alpha.
               ! The iron deficiency term is included here as a multiplier of the 
               ! thetamax_fe to represent the importance of Fe in forming chlorophyll
               ! accessory antennae, which do not affect the Chl:C but still affect the
               ! phytoplankton ability to use light (eg Stzrepek & Harrison Nature 
               ! 2004).
               !-----------------------------------------------------------------------

               ! I_k [W/m2]
               irrk(ji,jj,jk) = pc_m(ji,jj,jk) / (epsln + alpha_chl*thetamax_fe) + irr_mem(ji,jj,jk)*0.5d0

               !-----------------------------------------------------------------------
               ! Now we can calculate the carbon-specific photosynthesis rate, pc_tot.
               !-----------------------------------------------------------------------
               
               ![s-1]
               pc_tot = pc_m(ji,jj,jk)*(1.d0-EXP(-irr_mix(ji,jj,jk)/(irrk(ji,jj,jk)+epsln)))
               !pc_tot = pc_m(ji,jj,jk)*(1.d0-EXP(-irr_inst(ji,jj,jk)/(irrk(ji,jj,jk)+epsln)))

               !-----------------------------------------------------------------------
               ! Next, we account for the maintenance effort that phytoplankton must 
               ! exert in order to combat decay. This is prescribed as a fraction of the
               ! light-saturated photosynthesis rate, resp_frac. The result of this is 
               ! to set a level of energy availability below which net growth (and 
               ! therefore nutrient uptake) is zero, given by resp_frac * pc_m.
               !-----------------------------------------------------------------------

               ! Net total production [s-1]
               mu(ji,jj,jk) = MAX (0.d0,pc_tot-resp_frac*pc_m(ji,jj,jk))

               !-----------------------------------------------------------------------
               ! We now must convert this net carbon-specific growth rate to nutrient 
               ! uptake rates, the quantities we are interested in. Since we have no 
               ! explicit biomass tracer, we use the result of Dunne et al. (GBC, 2005) 
               ! to calculate an implicit biomass from the uptake rate through the  
               ! application of a simple idealized grazing law. This has the effect of 
               ! reducing uptake in low growth-rate regimes and increasing uptake in 
               ! high growth-rate regimes - essentially a non-linear amplification of 
               ! the growth rate variability. The result is:
               !-----------------------------------------------------------------------

               ! Biomass (units of pstar [mol P/m3])
               mulamb0expkT = mu(ji,jj,jk)/(lambda0*expkT(ji,jj,jk))  ![no units]
               biomass_p_ts = p_star*mulamb0expkT*(1.d0+(mulamb0expkT)**2)

               IF (kt==nittrc000) biomass_p(ji,jj,jk)=biomass_p_ts

               biomass_p(ji,jj,jk) =   biomass_p(ji,jj,jk) &
                                    + (biomass_p_ts-biomass_p(ji,jj,jk))*MIN(1.d0,gam_biomass*rfact)!*tmask(ji,jj,jk)

               !if (ji==80 .and. jj==60 .and. jk==1) write(*,'(I3,5(1X,E11.4))') kt, &
               !pc_0, expkT(ji,jj,jk), po4_up, def_fe(ji,jj,jk), pc_tot,mu(ji,jj,jk),biomass_p(ji,jj,jk)

               !-----------------------------------------------------------------------
               ! We can now use the diagnostic biomass to calculate the chlorophyll
               ! concentration:
               !-----------------------------------------------------------------------

               ! Chl:C ration [g chl/g C]
               theta   = thetamax_fe / (1.d0 + (thetamax_fe*alpha_chl*irr_mem(ji,jj,jk))/(2.d0*pc_m(ji,jj,jk)+epsln) )

               ! Chl biomass [ug chl/m3]
               chl_dia = biomass_p(ji,jj,jk) * c2p * 12.011e+6 * theta !* tmask(ji,jj,jk) 
               chl_bling(ji,jj,jk) = MAX(chl_min, chl_dia)

               !--------------------------------------------------
               ! PHOSPHORUS CYCLE
               !--------------------------------------------------
               ! The uptake of nutrients is assumed to contribute to the growth of
               ! phytoplankton, which subsequently die and are consumed by heterotrophs.
               ! This can involve the transfer of nutrient elements between many
               ! organic pools, both particulate and dissolved, with complex histories.
               ! We take a simple approach here, partitioning the total uptake into two
               ! fractions - sinking and non-sinking - as a function of temperature, 
               ! following Dunne et al. (2005). 
               ! Then, the non-sinking fraction is further subdivided, such that the 
               ! majority is recycled instantaneously to the inorganic nutrient pool,
               ! representing the fast turnover of labile dissolved organic matter via
               ! the microbial loop, and the remainder is converted to semi-labile
               ! dissolved organic matter. Iron and phosphorus are treated identically 
               ! for the first step, but all iron is recycled instantaneously in the
               ! second step (i.e. there is no dissolved organic iron pool).
               !-----------------------------------------------------------------------

               ! Phosphorous uptake flux [mol P/m3/s]
               jp_uptake(ji,jj,jk)=biomass_p(ji,jj,jk)*mu(ji,jj,jk)

               ! [no units]
               frac_pop=(phi_sm+phi_lg*(mulamb0expkT)**2) / (1+(mulamb0expkT)**2) * EXP(kappa_remin*tsn(ji,jj,jk,jp_tem))

               ! [mol P/m3/s]
               jp_pop(ji,jj,jk)=frac_pop*jp_uptake(ji,jj,jk)

               !-----------------------------------------------------------------------
               ! Then the remainder is divided between instantaneously recycled and
               ! long-lived dissolved organic matter,
               !-----------------------------------------------------------------------
               ! [mol P/m3/s]
               jp_dop=phi_dop*(jp_uptake(ji,jj,jk)-jp_pop(ji,jj,jk))

               ! [mol P/m3/s]
               jp_recycle(ji,jj,jk)=jp_uptake(ji,jj,jk)-jp_pop(ji,jj,jk)-jp_dop

!DO_CARBON<
               !---------------------------------------------------------------------
               ! As a helpful diagnostic, the implied fraction of production by large 
               ! phytoplankton is calculated, also following Dunne et al. 2005. This
               ! could be done more simply, but is done here in a complicated
               ! way as a sanity check. Looks fine.
               ! Note the calculation is made in P units, rather than C.
               ! This is also used for the CaCO3 production.

               ! MC. Note: s_over_p = (mu/lambda)**2 for gamma_b=0. 
               ! To do the maths consider jp_uptake=mu*B and lambda=lambda0*expkT  
               ! [no units]
               s_over_p = ( ( 1.d0 + 4.d0 * jp_uptake(ji,jj,jk) / (lambda0*expkT(ji,jj,jk)*p_star) )**0.5 - 1.d0 ) / 2.d0

               ![no units]
               frac_lg(ji,jj,jk) = s_over_p / (1.d0 + s_over_p)

               !-----------------------------------------------------------------------
               ! Calcium carbonate production
               ! Alkalinity is consumed through the production of CaCO3. Here, this is
               ! simply a linear function of the implied growth rate of small
               ! phytoplankton, which gave a reasonably good fit to the global 
               ! observational synthesis of Dunne (in prep., 2009).
               ! This is consistent with the findings of Jin et al. (GBC,2006).
               !-----------------------------------------------------------------------

               ! [mol Ca/m3/s]
               jca_uptake(ji,jj,jk) = (1.d0-frac_lg(ji,jj,jk))*jp_uptake(ji,jj,jk)*ca2p         
!>DO_CARBON
               !---------------------------------------------------------------------
               ! IRON
               !---------------------------------------------------------------------
               ! Iron is then taken up as a function of PO4 uptake and iron limitation,
               ! with a maximum Fe:P uptake ratio of fe2p_max:
               !-----------------------------------------------------------------------

               ! Iron uptake/POFe/recycle fluxes (units of jp_uptake [mol Fe/m3/s])
               jfe_uptake(ji,jj,jk) =jp_uptake(ji,jj,jk)*fe2p_up
               jfe_pofe             =frac_pop*jfe_uptake(ji,jj,jk)
               jfe_recycle(ji,jj,jk)=jfe_uptake(ji,jj,jk)-jfe_pofe

               !-----------------------------------------------------------------------
               ! Calculate free and inorganically associated iron concentration for
               ! scavenging.
               ! We assume that there is a 
               ! spectrum of iron ligands present in seawater, with varying binding
               ! strengths and whose composition varies with light and iron 
               ! concentrations. For example, photodissocation of ligand complexes 
               ! occurs under bright light, weakening the binding strength 
               ! (e.g. Barbeau et al., Nature 2001), while at very low iron 
               ! concentrations (order kfe_eq_lig_femin), siderophores are thought
               ! to be produced as a response to extreme
               ! iron stress.
               ! In anoxic waters, iron should be reduced, and therefore mostly 
               ! immune to scavenging. Easiest way to do this is to skip the feprime
               ! calculation if oxygen is less than 0.
               !-----------------------------------------------------------------------

               ! Calculate Fe prime [mol Fe/m3]
               if (foxy > oxy_min ) then
                 dum5       = irr_inst(ji,jj,jk)**2/(kfe_eq_lig_irr**2+irr_inst(ji,jj,jk)**2)        ! [no units]
                 dum2       = max(  0.e0, min( 1.e0, 1.2e0*(ffed-kfe_eq_lig_femin)/(epsln+ffed) )  ) ! [no units]
                 kfe_eq_lig(ji,jj,jk) = kfe_eq_lig_max -(kfe_eq_lig_max-kfe_eq_lig_min)*dum5*dum2    ! [m3/mol Fe]
                 feprime(ji,jj,jk)= 1.e0 + kfe_eq_lig(ji,jj,jk)*(felig_bkg - ffed)                   ! [no units]
                 ! units of (kfe_eq_lig)-1[mol Fe/m3]
                 feprime(ji,jj,jk)= (-feprime(ji,jj,jk) + sqrt(feprime(ji,jj,jk)**2 + 4.e0*kfe_eq_lig(ji,jj,jk)*ffed) )/(2.e0*kfe_eq_lig(ji,jj,jk))
               else
                 feprime(ji,jj,jk) = 0.e0
               endif

               ! [mol Fe/m3/s]
               jfe_ads_inorg(ji,jj,jk) = min( 0.5d0/rfact, kfe_inorg*sqrt(feprime(ji,jj,jk)) )*feprime(ji,jj,jk)

               ! [mol Fe/m3/s]
               dum4(jk) = jfe_pofe + jfe_ads_inorg(ji,jj,jk)

               !--------------------------------------------------
               ! COMPUTE TRENDS w/o remineralization processes
               !--------------------------------------------------
               !
               ! [mol P/m3/s]
               jpo4(ji,jj,jk) =   jp_recycle(ji,jj,jk) + gamma_dop*fdop -jp_uptake(ji,jj,jk)
               jdop(ji,jj,jk) = - gamma_dop*fdop + phi_dop*(jp_uptake(ji,jj,jk)-jp_pop(ji,jj,jk))
               ! [mol Fe/m3/s]
               jfed(ji,jj,jk) =   jfe_recycle(ji,jj,jk)-jfe_uptake(ji,jj,jk)-jfe_ads_inorg(ji,jj,jk)

               !DO_CARBON
               ! [mol C/m3/s]
               jdic(ji,jj,jk) = - jca_uptake(ji,jj,jk) 
               jalk(ji,jj,jk) = - 2.d0*jca_uptake(ji,jj,jk)

            ENDDO


            !-----------------------------------------------------------------------
            ! SINKING AND REMINERALIZATION
            !-----------------------------------------------------------------------
            ! Calculate the remineralization lengthscale matrix, zremin, a function 
            ! of z. Sinking rate (wsink) is constant over the upper wsink0_z metres,
            ! then  increases linearly with depth.
            ! The remineralization rate is a function of oxygen concentrations,
            ! following a Holling type 2 dependence, decreasing to a minimum value
            ! of remin_min. This is ad hoc, following work by Bianchi, Sarmiento,
            ! Galbraith and Kwon (unpublished).
            !-----------------------------------------------------------------------
            ! In general, the flux at the bottom of a grid cell should equal
            ! Fb = (Ft + Prod*dz) / (1 + zremin*dz)
            ! where Ft is the flux at the top, and prod*dz is the integrated 
            ! production of new sinking particles within the layer.
            ! Since Ft=0 in the first layer,
            !---------------------------------------------------------------------
            ! Calculate co3_ion first. Used to compute caco3 remineralization lengthscale
            ! Also calculate co2 fluxes csurf and alpha for the next round of exchange
            ! Note a scaled value of the PO4, rather than SiOH3, is used for all 
            ! calculations since there is no prognostic silica cycle 
            ! GFDL subroutine used to compute H+ that includes the OCMIP2 protocol
            !---------------------------------------------------------------------

            ! k=1: surface layer
            jk=1
            ! [m]
            zzz  =fse3t(ji,jj,jk)

            ! Sinking rate [m/s]
            IF (zzz .lt. wsink0_z) THEN
               wsink=wsink0
            ELSE
               wsink=wsink0+wsink_acc*(zzz-wsink0_z)
            ENDIF

            ! Remineralization lengthscale (oxygen dependent process)
            ! [mol O2/m3/s]
            foxy = trn(ji,jj,jk,jpOxy_bling)       
            ! [no units]
            oxy_up =foxy**2 / (koxy**2 + foxy**2) 
            ! [m-1]
            zremin(ji,jj,jk) =gamma_pop*(oxy_up*(1.d0-remin_min)+remin_min)/(epsln+wsink)
 
            ! [mol P/m2/s]
            fpop(ji,jj,jk)    = jp_pop(ji,jj,jk)*fse3t(ji,jj,jk)/(1.d0+fse3t(ji,jj,jk)*zremin(ji,jj,jk)) 
            ! [mol P/m3/s]
            jp_remin(ji,jj,jk)=(jp_pop(ji,jj,jk)*fse3t(ji,jj,jk)-fpop(ji,jj,jk))/(epsln+fse3t(ji,jj,jk))

            !-----------------------------------------------------------------------
!DO_CARBON<-
            ! Using Sayles for calcite solubility (will change to Mucci later) [mol/kg]
            ! NOTE: co3_solubility units cancel with co3_ion units also in mol/kg.
            ! This algorithm is ported exactly from blingv0 in GFDL, which units are in mol/kg.
            co3_solubility=max(  4.95e-7    * exp( 0.05021/(tsn(ji,jj,jk,jp_tem)+273.15)*zzz )  &
                                *3.42031e-3 * rho0_co3sol * rho0_co3sol / max(epsln, tsn(ji,jj,jk,jp_sal))    &
                               , epsln )

            ! CaCO3 dissolution lengthscale is a function of the saturation
            ! state, CO3 / CO3solubility, such that the lengthscale decreases
            ! from infinity at CO3 > CO3solubility to zero at CO3 = 0.

            ! [m-1]
            zremin_caco3(ji,jj,jk)=  1.d0 / ca_remin_depth  &
                                  *(1.d0 - min(1.d0, co3_ion(ji,jj,jk)/(co3_solubility+epsln) ))

            ! Generate CaCO3 sinking flux, and dissolve it through the water column. 
            ! Same as for other elements - see below for more detailed explanation.

            ! MC note: fcaco3 is computed as fpop, and jca_remin as jp_remin
            ! [mol C/m2/s]
            fcaco3    (ji,jj,1)= jca_uptake(ji,jj,jk)*fse3t(ji,jj,jk)/(1.d0+fse3t(ji,jj,jk)*zremin_caco3(ji,jj,jk))
            jca_reminp(ji,jj,1)=(jca_uptake(ji,jj,jk)*fse3t(ji,jj,jk)-fcaco3(ji,jj,jk))/(epsln+fse3t(ji,jj,jk))

            ! Add remineralization terms to sources
            ! [mol C/m2/s]
            jdic(ji,jj,jk)=jdic(ji,jj,jk)+     jca_reminp(ji,jj,jk)
            ! [mol eq/m2/s]
            jalk(ji,jj,jk)=jalk(ji,jj,jk)+2.d0*jca_reminp(ji,jj,jk)
!DO_CARBON>
            !-----------------------------------------------------------------------
            ! Now, calculate the Fe adsorption using this fpop:
            ! The absolute first order rate constant is calculated from the 
            ! concentration of organic particles, after Parekh et al. (2005). Never
            !  allowed to be greater than 1/2dt for numerical stability.
            !-----------------------------------------------------------------------

            !Iron [gC/m3]
            dum3                 = (fpop(ji,jj,jk)*c2p*12.011d0/(epsln+wsink))**(0.58)
            ! [mol Fe/m3/s]
            jfe_ads_org(ji,jj,jk) = min (0.5d0/rfact, kfe_org*dum3)*feprime(ji,jj,jk)
            ! [mol Fe/m2/s]
            dum4(jk)              =( dum4(jk)+jfe_ads_org(ji,jj,jk) )*fse3t(ji,jj,jk)
            fpofe(ji,jj,jk)       =  dum4(jk)/(1.d0+fse3t(ji,jj,jk)*zremin(ji,jj,jk))
            ! [mol Fe/m3/s]
            jfe_remin(ji,jj,jk)   =( dum4(jk)-fpofe(ji,jj,jk))/(epsln+fse3t(ji,jj,jk) )

            ! Add remineralization terms to trends
            ! [mol P/m3/s]
            jpo4(ji,jj,1)=jpo4(ji,jj,1)+(1.d0-phi_dop)*jp_remin(ji,jj,1)
            jdop(ji,jj,1)=jdop(ji,jj,1)+       phi_dop*jp_remin(ji,jj,1)
            ! [mol Fe/m3/s]
            jfed(ji,jj,1)=jfed(ji,jj,1)+jfe_remin(ji,jj,1)-jfe_ads_org(ji,jj,1)   

            ! k=2:NK: rest of the water column
            DO jk=2, jpk

               fpopkm1   = fpop  (ji,jj,jk-1)
               fcaco3km1 = fcaco3(ji,jj,jk-1)
               fpofekm1  = fpofe (ji,jj,jk-1)
    
               zzz  =zzz+fse3t(ji,jj,jk)

               IF (zzz .lt. wsink0_z) THEN
                  wsink=wsink0
               ELSE
                  wsink=wsink0+wsink_acc*(zzz-wsink0_z)
               ENDIF

               ! Phosphorous
               foxy = trn(ji,jj,jk,jpOxy_bling)
               oxy_up =foxy**2 / (koxy**2 + foxy**2)
               ! [m-1]
               zremin(ji,jj,jk) =gamma_pop*(oxy_up*(1.d0-remin_min)+remin_min)/(epsln+wsink)

               ! [mol P/m2/s]
               fpop(ji,jj,jk)    =(fpopkm1+jp_pop(ji,jj,jk)*fse3t(ji,jj,jk))/(1.d0+fse3t(ji,jj,jk)*zremin(ji,jj,jk)) 
               ! [mol P/kg/s]
               jp_remin(ji,jj,jk)=(fpopkm1+jp_pop(ji,jj,jk)*fse3t(ji,jj,jk)-fpop(ji,jj,jk))/(epsln+fse3t(ji,jj,jk))

!DO_CARBON<
               co3_solubility=max(  4.95e-7    * exp( 0.05021/(tsn(ji,jj,jk,jp_tem)+273.15)*zzz )  &
                                   *3.42031e-3 * rho0_co3sol * rho0_co3sol / max(epsln, tsn(ji,jj,jk,jp_sal))    &
                                  , epsln )

               zremin_caco3(ji,jj,jk)=  1.d0 / ca_remin_depth  &
                                      *(1.d0 - min(1.d0, co3_ion(ji,jj,jk)/(co3_solubility+epsln) ))

               fcaco3(ji,jj,jk)=(fcaco3km1+jca_uptake(ji,jj,jk)*fse3t(ji,jj,jk))   &
                                / ( 1.d0+fse3t(ji,jj,jk) * zremin_caco3(ji,jj,jk))

               jca_reminp(ji,jj,jk)=  (fcaco3km1+jca_uptake(ji,jj,jk)*fse3t(ji,jj,jk)-fcaco3(ji,jj,jk)) &
                                    / (epsln+fse3t(ji,jj,jk))

               ! Add remineralization terms to sources
               ! [mol C/m3/s]
               jdic(ji,jj,jk)=jdic(ji,jj,jk)+     jca_reminp(ji,jj,jk)
               ! [mol eq/m3/s]
               jalk(ji,jj,jk)=jalk(ji,jj,jk)+2.d0*jca_reminp(ji,jj,jk)
!DO_CARBON>
               ! Iron
               ! [gC/m3]
               dum3            = (fpop(ji,jj,jk)*c2p*12.011d0/(epsln+wsink))**(0.58)
               ! [mol Fe/m3/s]
               jfe_ads_org(ji,jj,jk) = min (0.5d0/rfact, kfe_org*dum3)*feprime(ji,jj,jk)
               ! [mol Fe/m2/s]
               dum4(jk)        =( dum4(jk)+jfe_ads_org(ji,jj,jk) )*fse3t(ji,jj,jk)
               fpofe(ji,jj,jk) = (fpofekm1 + dum4(jk))/(1.d0+fse3t(ji,jj,jk)*zremin(ji,jj,jk))
               ! [mol Fe/m3/s]
               jfe_remin(ji,jj,jk) = (fpofekm1+dum4(jk)-fpofe(ji,jj,jk))/(epsln+fse3t(ji,jj,jk))

               ! Save fPOP and fPOFe at the bottom grid cell to compute bottom fluxes
               ! [mol/m2/s]
               IF (jk==mbkt(ji,jj)) THEN
                  fpop_b  (ji,jj) = fpop(ji,jj,jk)
                  fpofe_b (ji,jj) = fpofe(ji,jj,jk)
                  fcaco3_b(ji,jj) = fcaco3(ji,jj,jk)  !DO_CARBON
               ENDIF
           
               ! Add remineralization terms to trends
               ! [mol P/m3/s]
               jpo4(ji,jj,jk)=jpo4(ji,jj,jk)+(1.d0-phi_dop)*jp_remin(ji,jj,jk)
               jdop(ji,jj,jk)=jdop(ji,jj,jk)+      phi_dop *jp_remin(ji,jj,jk)
               ! [mol Fe/m3/s]
               jfed(ji,jj,jk)=jfed(ji,jj,jk)+jfe_remin(ji,jj,jk)-jfe_ads_org(ji,jj,jk)

            ENDDO 

            !-----------------------------------------------------------------------
            ! OXYGEN
            !-----------------------------------------------------------------------
            ! Assuming constant P:O ratio.
            ! Optional prevention of negative oxygen (does not conserve ocean 
            ! redox potential) or alternatively it can be allowed to go negative, 
            ! keeping track of an implicit nitrate deficit 
            ! plus sulfate reduction.
            !-----------------------------------------------------------------------

            DO jk=1, jpk
               ! [mol O2/m3]
               foxy = trn(ji,jj,jk,jpOxy_bling)
               IF ( (ln_prev_o2lt0) .and. (foxy<oxy_min) ) then
                   ! O2 below hypoxic threshold, no O2 consumption
                   joxy(ji,jj,jk)=0.d0
               ELSE
                   ! [mol O2/m3/s]
                   joxy(ji,jj,jk)=-oxy2p*jpo4(ji,jj,jk)
               ENDIF
            ENDDO

         ENDDO
      ENDDO

      ! Add dissolved inorganic carbon terms from po4 final fluxes
      jdic(:,:,:)=jdic(:,:,:)+jpo4(:,:,:)*c2p
      jalk(:,:,:)=jalk(:,:,:)-jpo4(:,:,:)*n2p

      tra(:,:,:,jpPO4_bling) = tra(:,:,:,jpPO4_bling) + jpo4(:,:,:)*rfact
      tra(:,:,:,jpDOP_bling) = tra(:,:,:,jpDOP_bling) + jdop(:,:,:)*rfact
      tra(:,:,:,jpFed_bling) = tra(:,:,:,jpFed_bling) + jfed(:,:,:)*rfact
      tra(:,:,:,jpOxy_bling) = tra(:,:,:,jpOxy_bling) + joxy(:,:,:)*rfact

      !jdic(:,:,:)=0.e0
      tra(:,:,:,jpDIC_bling) = tra(:,:,:,jpDIC_bling) + jdic(:,:,:)*rfact

      !jalk(:,:,:)=0.e0
      tra(:,:,:,jpalk_bling) = tra(:,:,:,jpalk_bling) + jalk(:,:,:)*rfact

      !test if concentrations fall below 0
      xnegtr(:,:,:) = 1.e0
      DO jn = jp_bling0, jp_bling1
         DO jk = 1, jpk
            DO jj = 1, jpj
               DO ji = 1, jpi
                  IF( ( trn(ji,jj,jk,jn) + tra(ji,jj,jk,jn) ) < 0.e0 ) THEN 
                     ztra             = ABS(  ( trn(ji,jj,jk,jn) - rtrn ) &
                                            / ( tra(ji,jj,jk,jn) + rtrn ) )
                     xnegtr(ji,jj,jk) = MIN( xnegtr(ji,jj,jk),  ztra )
                  ENDIF
              END DO
            END DO
         END DO
      END DO

      ! Prgonostic tracer fields
      DO jn = jp_bling0, jp_bling1
         !write(*,'(I3,2(1X,E11.4))') jp_bling0, trn(ji,jj,jk,jn),tra(ji,jj,jk,jn)
         trn(:,:,:,jn) = trn(:,:,:,jn) + xnegtr(:,:,:) * tra(:,:,:,jn)
      END DO

      !DO jn = jp_bling0, jp_bling1
      !   trn(:,:,:,jn) = trn(:,:,:,jn) + tra(:,:,:,jn)
      !   DO jk = 1, jpk
      !      DO jj = 1, jpj
      !         DO ji = 1, jpi
      !           trn(ji,jj,jk,jn)=MAX( 0.e0, trn(ji,jj,jk,jn) )
      !         END DO
      !      END DO
      !   END DO
      !END DO

      ! Copy new arrays to trb (tracer fields before) and set tra to zero
      ! to compute tracer gradients with tracer fields after ecological forcing
      tra(:,:,:,jp_bling0:jp_bling1) = 0.e0

      ! add external fluxes
      CALL trc_ext_bling (kt)

      ! Copy to trb to use BLING updated tracer values to compute transport
      ! trends
      DO jn=jp_bling0, jp_bling1
         trb(:,:,:,jn)=trn(:,:,:,jn)
      ENDDO

      DO jn=jp_bling0, jp_bling1
         CALL lbc_lnk( trn(:,:,:,jn), 'T', 1. )
         CALL lbc_lnk( trb(:,:,:,jn), 'T', 1. )
         CALL lbc_lnk( tra(:,:,:,jn), 'T', 1. )
      ENDDO

      !IF ( ln_bling_mass .and.  narea == 1 ) THEN
      IF ( ln_bling_mass  ) THEN

        IF( kt == nittrc000 ) THEN 
          CALL ctl_opn( numsms ,  'sms.budget' , 'REPLACE', 'FORMATTED', 'SEQUENTIAL', -1, 6, .FALSE., narea )
          CALL ctl_opn( numsms2, 'sms2.budget' , 'REPLACE', 'FORMATTED', 'SEQUENTIAL', -1, 6, .FALSE., narea )
        ENDIF


        WRITE(UNIT=numsms,FMT='(i10,6(3x,e18.10))')  kt  &
                           , glob_sum( jp_recycle(:,:,:)*cvol(:,:,:) )*rfact &
                           , glob_sum(   jp_remin(:,:,:)*cvol(:,:,:) )*(1.d0-phi_dop)*rfact &
                           , glob_sum ( -jp_uptake(:,:,:)*cvol(:,:,:) )*rfact &
     , glob_sum( phi_dop*(jp_uptake(:,:,:)-jp_pop(:,:,:))*cvol(:,:,:) )*rfact &
              , glob_sum ( jp_remin(:,:,:)*phi_dop*cvol(:,:,:) )*rfact

        WRITE(UNIT=numsms2,FMT='(i10,3(3x,e18.10))')  kt  &
                          , glob_sum (jpo4(:,:,:)*cvol(:,:,:))*rfact &
                          , glob_sum (jdop(:,:,:)*cvol(:,:,:))*rfact &
                          , glob_sum (jdic(:,:,:)*cvol(:,:,:))*rfact 
      ENDIF

      IF (ln_diatrc) THEN
         IF (lk_iomput) THEN
            CALL iom_put(      "expkT"  ,        expkT(:,:,:)*tmask(:,:,:) )
            CALL iom_put(   "irr_inst"  ,     irr_inst(:,:,:)*tmask(:,:,:) )
            CALL iom_put(    "irr_mix"  ,      irr_mix(:,:,:)*tmask(:,:,:) )
            CALL iom_put(       "irrk"  ,         irrk(:,:,:)*tmask(:,:,:) )
            CALL iom_put(       "pc_m"  ,         pc_m(:,:,:)*tmask(:,:,:) )
            CALL iom_put(         "mu"  ,           mu(:,:,:)*tmask(:,:,:) )
            CALL iom_put(  "biomass_p"  ,    biomass_p(:,:,:)*tmask(:,:,:) )
            CALL iom_put(       "fpop"  ,         fpop(:,:,:)*tmask(:,:,:) )
            CALL iom_put(     "zremin"  ,       zremin(:,:,:)*tmask(:,:,:) )
            CALL iom_put(     "def_fe"  ,       def_fe(:,:,:)*tmask(:,:,:) )
            CALL iom_put(    "feprime"  ,      feprime(:,:,:)*tmask(:,:,:) )
            CALL iom_put( "kfe_eq_lig"  ,   kfe_eq_lig(:,:,:)*tmask(:,:,:) )
            CALL iom_put(      "fpofe"  ,        fpofe(:,:,:)*tmask(:,:,:) )
            !CALL iom_put(   "frac_pop"  ,    frac_pop(:,:,:)*tmask(:,:,:) )

            CALL iom_put(     "jp_pop"  ,       jp_pop(:,:,:)*fse3t(:,:,:)*tmask(:,:,:) )
            CALL iom_put(  "jp_uptake"  ,    jp_uptake(:,:,:)*fse3t(:,:,:)*tmask(:,:,:) )
            CALL iom_put(   "jp_remin"  ,     jp_remin(:,:,:)*fse3t(:,:,:)*tmask(:,:,:) )
            CALL iom_put( "jp_recycle"  ,   jp_recycle(:,:,:)*fse3t(:,:,:)*tmask(:,:,:) )                           
            CALL iom_put( "jfe_uptake"  ,   jfe_uptake(:,:,:)*fse3t(:,:,:)*tmask(:,:,:) )                          
            CALL iom_put(  "jfe_remin"  ,    jfe_remin(:,:,:)*fse3t(:,:,:)*tmask(:,:,:) )                          
            CALL iom_put( "jfe_recycle" ,  jfe_recycle(:,:,:)*fse3t(:,:,:)*tmask(:,:,:) )                          
            CALL iom_put( "jfe_ads_org" ,  jfe_ads_org(:,:,:)*fse3t(:,:,:)*tmask(:,:,:) )                           
            CALL iom_put("jfe_ads_inorg",jfe_ads_inorg(:,:,:)*fse3t(:,:,:)*tmask(:,:,:) )                          
            CALL iom_put(        "jpo4" ,         jpo4(:,:,:)*fse3t(:,:,:)*tmask(:,:,:) )
            CALL iom_put(        "jdop" ,         jdop(:,:,:)*fse3t(:,:,:)*tmask(:,:,:) )
            CALL iom_put(        "jfed" ,         jfed(:,:,:)*fse3t(:,:,:)*tmask(:,:,:) )
            CALL iom_put(        "joxy" ,         joxy(:,:,:)*fse3t(:,:,:)*tmask(:,:,:) )

!DO CARBON<-
            CALL iom_put( "jca_uptake" ,  jca_uptake(:,:,:)*fse3t(:,:,:)*tmask(:,:,:) )
            CALL iom_put( "jca_reminp" ,  jca_reminp(:,:,:)*fse3t(:,:,:)*tmask(:,:,:) )
            CALL iom_put(      "jdic"  ,        jdic(:,:,:)*fse3t(:,:,:)*tmask(:,:,:) )
            CALL iom_put(      "jalk"  ,        jalk(:,:,:)*fse3t(:,:,:)*tmask(:,:,:) )

            CALL iom_put(    "frac_lg"  ,      frac_lg(:,:,:)*tmask(:,:,:) )
            CALL iom_put("zremin_caco3" , zremin_caco3(:,:,:)*tmask(:,:,:) )
            CALL iom_put(     "fcaco3"  ,       fcaco3(:,:,:)*tmask(:,:,:) )
            CALL iom_put(  "pco2_surf"  ,    pco2_surf(:,:)  *tmask(:,:,1) )
!DO CARBON>

         ENDIF

         IF ( .NOT. lk_iomput )  THEN
           trc2d( :,:,jp_bling0_2d ) = pco2_surf(:,:)*tmask(:,:,1)

           trc3d( :,:,:,jp_bling0_3d    ) =    chl_bling(:,:,:)*tmask(:,:,:)
           trc3d( :,:,:,jp_bling0_3d+1  ) =         jpo4(:,:,:)*fse3t(:,:,:)*tmask(:,:,:)
           trc3d( :,:,:,jp_bling0_3d+2  ) =         jdop(:,:,:)*fse3t(:,:,:)*tmask(:,:,:)
           trc3d( :,:,:,jp_bling0_3d+3  ) =       jp_pop(:,:,:)*fse3t(:,:,:)*tmask(:,:,:)
           trc3d( :,:,:,jp_bling0_3d+4  ) =    jp_uptake(:,:,:)*fse3t(:,:,:)*tmask(:,:,:)
           trc3d( :,:,:,jp_bling0_3d+5  ) =   jp_recycle(:,:,:)*fse3t(:,:,:)*tmask(:,:,:)
           trc3d( :,:,:,jp_bling0_3d+6  ) =     jp_remin(:,:,:)*fse3t(:,:,:)*tmask(:,:,:)
           trc3d( :,:,:,jp_bling0_3d+7  ) =         jfed(:,:,:)*fse3t(:,:,:)*tmask(:,:,:)
           trc3d( :,:,:,jp_bling0_3d+8  ) =   jfe_uptake(:,:,:)*fse3t(:,:,:)*tmask(:,:,:)
           trc3d( :,:,:,jp_bling0_3d+9  ) =  jfe_recycle(:,:,:)*fse3t(:,:,:)*tmask(:,:,:)
           trc3d( :,:,:,jp_bling0_3d+10 ) =    jfe_remin(:,:,:)*fse3t(:,:,:)*tmask(:,:,:)
           trc3d( :,:,:,jp_bling0_3d+11 ) =  jfe_ads_org(:,:,:)*fse3t(:,:,:)*tmask(:,:,:)
           trc3d( :,:,:,jp_bling0_3d+12 ) =jfe_ads_inorg(:,:,:)*fse3t(:,:,:)*tmask(:,:,:)
           trc3d( :,:,:,jp_bling0_3d+13 ) =         joxy(:,:,:)*fse3t(:,:,:)*tmask(:,:,:)
           trc3d( :,:,:,jp_bling0_3d+14 ) =         jdic(:,:,:)*fse3t(:,:,:)*tmask(:,:,:)
           trc3d( :,:,:,jp_bling0_3d+15 ) =   jca_uptake(:,:,:)*fse3t(:,:,:)*tmask(:,:,:)
           trc3d( :,:,:,jp_bling0_3d+16 ) =   jca_reminp(:,:,:)*fse3t(:,:,:)*tmask(:,:,:)
           trc3d( :,:,:,jp_bling0_3d+17 ) =         jalk(:,:,:)*fse3t(:,:,:)*tmask(:,:,:)

           trc3d( :,:,:,jp_bling0_3d+18 ) =        fpop(:,:,:)*tmask(:,:,:)
           trc3d( :,:,:,jp_bling0_3d+19 ) =       fpofe(:,:,:)*tmask(:,:,:)
           trc3d( :,:,:,jp_bling0_3d+20 ) =       expkT(:,:,:)*tmask(:,:,:)
           trc3d( :,:,:,jp_bling0_3d+21 ) =    irr_inst(:,:,:)*tmask(:,:,:)
           trc3d( :,:,:,jp_bling0_3d+22 ) =     irr_mix(:,:,:)*tmask(:,:,:)
           trc3d( :,:,:,jp_bling0_3d+23 ) =        irrk(:,:,:)*tmask(:,:,:)
           trc3d( :,:,:,jp_bling0_3d+24 ) =        pc_m(:,:,:)*tmask(:,:,:)
           trc3d( :,:,:,jp_bling0_3d+25 ) =          mu(:,:,:)*tmask(:,:,:)
           trc3d( :,:,:,jp_bling0_3d+26 ) =   biomass_p(:,:,:)*tmask(:,:,:)
           trc3d( :,:,:,jp_bling0_3d+27 ) =      zremin(:,:,:)*tmask(:,:,:)
           trc3d( :,:,:,jp_bling0_3d+28 ) =      def_fe(:,:,:)*tmask(:,:,:)
           trc3d( :,:,:,jp_bling0_3d+29 ) =     feprime(:,:,:)*tmask(:,:,:)
           trc3d( :,:,:,jp_bling0_3d+30 ) =  kfe_eq_lig(:,:,:)*tmask(:,:,:)
           trc3d( :,:,:,jp_bling0_3d+31 ) =     frac_lg(:,:,:)*tmask(:,:,:)
           trc3d( :,:,:,jp_bling0_3d+32 ) =zremin_caco3(:,:,:)*tmask(:,:,:)
           trc3d( :,:,:,jp_bling0_3d+33 ) =      fcaco3(:,:,:)*tmask(:,:,:)
           trc3d( :,:,:,jp_bling0_3d+34 ) =     co3_ion(:,:,:)*tmask(:,:,:)
           trc3d( :,:,:,jp_bling0_3d+35 ) =      htotal(:,:,:)*tmask(:,:,:)
         ENDIF
      ENDIF

      CALL wrk_dealloc( jpi, jpj, jpk, rho )
      CALL wrk_dealloc( jpi, jpj, jpk, expkT, irr_inst, irr_mix, irrk, pc_m, mu )
      CALL wrk_dealloc( jpi, jpj, jpk, def_fe, feprime, kfe_eq_lig, fpofe )
      CALL wrk_dealloc( jpi, jpj, jpk, jp_pop, fpop, zremin )
      CALL wrk_dealloc( jpi, jpj, jpk, jp_uptake, jp_remin, jp_recycle )
      CALL wrk_dealloc( jpi, jpj, jpk, jfe_uptake, jfe_remin, jfe_recycle )
      CALL wrk_dealloc( jpi, jpj, jpk, jfe_ads_inorg, jfe_ads_org )
      CALL wrk_dealloc( jpi, jpj, jpk, jpo4, jdop, jfed, joxy )
!DO_CARBON
      CALL wrk_dealloc( jpi, jpj, jpk, jdic, jalk )
      CALL wrk_dealloc( jpi, jpj, jpk, jca_uptake, jca_reminp, fcaco3 )
      CALL wrk_dealloc( jpi, jpj, jpk, zremin_caco3, frac_lg )
      CALL wrk_dealloc( jpi, jpj     , pco2_surf )

      CALL wrk_dealloc( jpk, dum4 )
      CALL wrk_dealloc( jpi, jpj, jpk, xnegtr  )

      IF( nn_timing == 1 )  CALL timing_stop('trc_sms_bling')
      !
   END SUBROUTINE trc_sms_bling

   SUBROUTINE trc_sms_bling_mass_conserv (kt)

      INTEGER, INTENT(in) ::   kt   ! ocean time-step index

      REAL(wp) :: globvol
      REAL(wp) :: sum_phosp, sum_fed, sum_oxy, sum_carbon
      
      !!-------------------------------------------------------------------

      globvol= glob_sum( cvol(:,:,:) )

      IF( kt == nittrc000 ) THEN 

        CALL ctl_opn( numphp, 'php.budget' , 'REPLACE', 'FORMATTED', 'SEQUENTIAL', -1, 6, .FALSE., narea )
        CALL ctl_opn( numfed, 'fed.budget' , 'REPLACE', 'FORMATTED', 'SEQUENTIAL', -1, 6, .FALSE., narea )
        CALL ctl_opn( numoxy, 'oxy.budget' , 'REPLACE', 'FORMATTED', 'SEQUENTIAL', -1, 6, .FALSE., narea )
        CALL ctl_opn( numcar, 'car.budget' , 'REPLACE', 'FORMATTED', 'SEQUENTIAL', -1, 6, .FALSE., narea )

      ENDIF

      ! total mass of BLING tracers
      ! PO4+DOP is conserved over the whole domain
      sum_phosp  = glob_sum ( ( trn(:,:,:,jpPO4_bling) + trn(:,:,:,jpDOP_bling)     )*cvol(:,:,:)  )

      ! DIC+DOP*106 is conserved over the whole domain w/o air-sea fluxes
      sum_carbon = glob_sum ( ( trn(:,:,:,jpDIC_bling) + trn(:,:,:,jpDOP_bling)*106 )*cvol(:,:,:)  )

      ! Non-conservative tracers
      sum_fed    = glob_sum (   trn(:,:,:,jpFed_bling)*cvol(:,:,:)  )
      sum_oxy    = glob_sum (   trn(:,:,:,jpOxy_bling)*cvol(:,:,:)  )

      IF( lwp ) THEN
        WRITE(numphp,9500) kt,  sum_phosp , globvol
        WRITE(numfed,9500) kt,  sum_fed   , globvol
        WRITE(numoxy,9500) kt,  sum_oxy   , globvol
        WRITE(numcar,9500) kt,  sum_carbon, globvol
      ENDIF

9500  FORMAT(i10,2(3x,e18.10))    

   END SUBROUTINE trc_sms_bling_mass_conserv

#else
   !!----------------------------------------------------------------------
   !!   Dummy module                                       No BLINGv0 model
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_sms_bling             ! Empty routine
      WRITE(*,*) 'trc_sms_bling: You should not have seen this print! error?', kt
   END SUBROUTINE trc_sms_bling
#endif
   !!======================================================================
END MODULE trcsms_blingv0
