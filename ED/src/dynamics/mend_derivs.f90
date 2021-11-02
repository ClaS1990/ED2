!======================================
! TO DO:
!======================================
! 1. Gangsheng questions: 
!    (a) Why are N vnup_amb so different compared to Zhu?  
!    (b) Why are Michaelis-Menten constants so different?
!    (c) Some important numbers (microbial uptake of DOC) differ dramatically
!        between 2013 and 2015 papers.
!    (d) Note: For microbial P uptake Vmax and kmm, applying same temperature 
!        sensitivities as for nitrogen.
!    (e) Note: For plant N and P uptake Vmax and kmm, applying same 
!        temperature sensitivities as for microbes.
!    (f) What are the POM pools? Are they appropriate?
!    (g) Added rk4 integrator
!    (h) Added P cycle, but need realistic phosphatase.
!    (i) Added vertical layers.
!    (j) Added plant litter inputs, plant uptake, external cycles.
! 2. Add plant cost for transporter enzyme production?
! 3. Implement transfers between different detritus pools.
! 4. Remove any hard-codings of soil density.  Check appropriateness of slden values.
! 5. Implement parameter values, including temperature sensitivities, for
!    phosphatase dynamics.
!======================================
! NOTES:
! 1. All detritus is being dumped into a generic SOM pool.  Fine for now,
!    but I will eventually want to have transfers.
! 2. plant2som_exchange:  
!    (a) plant inputs were getting zeroed in the CLM-CNP scheme.  
!        This was called after the monthly updates.  This is a problem 
!        because (i) in CLM scheme, disturbance loadings did not occur 
!        till the following time step; (ii) I'd be getting zeros when using 
!        MEND.  To address (i), I am now calling CLM-CNP scheme after 
!        disturbance.  To address (ii), I am zeroing plant inputs at the
!        beginning of vegetation_dynamics, which is called daily.  
!        This means that inputs would be applied to MEND on the following day.
!    (b) Inputs applied in mend_exchange.f90, plant2som.
! 3. BNF set to constant rate.
! 4. co2_lost as heterotrophic respiration.  Had to make substantial edits.
!    (a) Edited bgc_dyn.f90 so that rh would not be updated.  
!    (b) Had to make edits to rk4_derivs.f90 so that it did not use rh or 
!        cwd_rh to update the co2 of the canopy air space.  
!    (c) Instead, I incremented canopy air space CO2 derivatives in 
!        som_feedback_fluxes.  
!    (d) The checks of the CO2 budget used the rh calculated at the beginning
!        of the day.  However, I am now calculating rh dynamically.  
!        Therefore, I had to disable the check on the CO2 budget.  There may
!        be a safer way to implement this.
! 5. Zhu: what are the units of [E^plant_P], [E^plant_N]?  Are they gBiomass, 
!    gC, gN, gP, or what?  Similarly:  K^{mic,NH4}_M, K^{mic,NO3}_M,
!    K^{plant,NH4}_M, K^{plant,NO3}_M, K^{plant,P}_M, K^{mic,P}_M.
!    Finally, enz2biomass_plant conversion factor units?  Emailed on 7/6.
!    REPLY: QING CONFIRMED THAT THIS IS THE CASE.
! 6. plant_input_{C,N,P} edits: in growth_balive, I set the 
!    C leakage to 0.  I left the extra_nitrogen and extra_P.  
!======================================

Module mend_derivs
  use mend_consts_coms, only: decomp_consts
  use ed_max_dims, only: n_pft
  implicit none

Contains

  subroutine mend_derivs_layer(npom, consts, &
       pom_c, input_pom_c, d_pom_c, dom_c,   &
       input_dom_c, d_dom_c, enz_pom_c, d_enz_pom_c, mom_c, d_mom_c, &
       qom_c, d_qom_c, enz_mom_c, d_enz_mom_c, amb_c, d_amb_c, dmb_c,  &
       d_dmb_c, co2_lost, d_co2_lost, nmin, d_nmin, &
       nitr, d_nitr, &
       pom_n, input_pom_n, d_pom_n, dom_n,   &
       input_dom_n, d_dom_n, enz_pom_n, d_enz_pom_n, mom_n, d_mom_n, &
       qom_n, d_qom_n, enz_mom_n, d_enz_mom_n, amb_n, d_amb_n, dmb_n,  &
       d_dmb_n, nh4, input_nh4, d_nh4, no3, input_no3, d_no3, &
       pom_p, input_pom_p, d_pom_p, dom_p,   &
       input_dom_p, d_dom_p, enz_pom_p, d_enz_pom_p, mom_p, d_mom_p, &
       qom_p, d_qom_p, enz_mom_p, d_enz_mom_p, amb_p, d_amb_p, dmb_p,  &
       d_dmb_p, psol, input_psol, d_psol, plab, d_plab, ngas_lost,  &
       d_ngas_lost, enz_ptase_c, d_enz_ptase_c, enz_ptase_n, d_enz_ptase_n, &
       enz_ptase_p, d_enz_ptase_p, &
       d_nh4_plant, nh4_bnf, d_nh4_bnf, d_no3_plant, &
       c_leach, d_c_leach, n_leach, d_n_leach, p_leach, d_p_leach, &
       d_p_plant, pocc, d_pocc, ppar, d_ppar, input_ppar, &
       enz_plant_n, enz_plant_p, vnh4up_plant, vno3up_plant, vpup_plant, &
       water_drainage, bulk_den, pi1, wfp)
    implicit none

    integer, intent(in) :: npom
    type(decomp_consts) :: consts

    real, intent(in) :: bulk_den
    real, intent(in) :: water_drainage
    real, intent(in) :: dom_c
    real, intent(out) :: d_dom_c
    real, intent(in) :: amb_c
    real, intent(out) :: d_amb_c
    real, intent(in) :: dmb_c
    real, intent(out) :: d_dmb_c
    real, intent(in) :: mom_c
    real, intent(out) :: d_mom_c
    real, intent(in) :: qom_c
    real, intent(out) :: d_qom_c
    real, dimension(npom), intent(in) :: pom_c
    real, dimension(npom), intent(out) :: d_pom_c
    real, dimension(npom), intent(in) :: input_pom_c
    real, intent(in) :: input_dom_c
    real, dimension(npom), intent(in) :: enz_pom_c
    real, dimension(npom), intent(out) :: d_enz_pom_c
    real, intent(in) :: enz_mom_c
    real, intent(out) :: d_enz_mom_c
    real, intent(in) :: enz_ptase_c
    real, intent(out) :: d_enz_ptase_c

    real, intent(in) :: dom_n
    real, intent(out) :: d_dom_n
    real, intent(in) :: amb_n
    real, intent(out) :: d_amb_n
    real, intent(in) :: dmb_n
    real, intent(out) :: d_dmb_n
    real, intent(in) :: mom_n
    real, intent(out) :: d_mom_n
    real, intent(in) :: qom_n
    real, intent(out) :: d_qom_n
    real, dimension(npom), intent(in) :: pom_n
    real, dimension(npom), intent(out) :: d_pom_n
    real, dimension(npom), intent(in) :: input_pom_n
    real, intent(in) :: input_dom_n
    real, dimension(npom), intent(in) :: enz_pom_n
    real, dimension(npom), intent(out) :: d_enz_pom_n
    real, intent(in) :: enz_mom_n
    real, intent(out) :: d_enz_mom_n
    real, intent(in) :: enz_ptase_n
    real, intent(out) :: d_enz_ptase_n

    real, intent(in) :: dom_p
    real, intent(out) :: d_dom_p
    real, intent(in) :: amb_p
    real, intent(out) :: d_amb_p
    real, intent(in) :: dmb_p
    real, intent(out) :: d_dmb_p
    real, intent(in) :: mom_p
    real, intent(out) :: d_mom_p
    real, intent(in) :: qom_p
    real, intent(out) :: d_qom_p
    real, dimension(npom), intent(in) :: pom_p
    real, dimension(npom), intent(out) :: d_pom_p
    real, dimension(npom), intent(in) :: input_pom_p
    real, intent(in) :: input_dom_p
    real, dimension(npom), intent(in) :: enz_pom_p
    real, dimension(npom), intent(out) :: d_enz_pom_p
    real, intent(in) :: enz_mom_p
    real, intent(out) :: d_enz_mom_p
    real, intent(in) :: enz_ptase_p
    real, intent(out) :: d_enz_ptase_p

    real, intent(in) :: co2_lost
    real, intent(out) :: d_co2_lost
    real, intent(in) :: nmin
    real, intent(out) :: d_nmin
    real, intent(in) :: nitr
    real, intent(out) :: d_nitr
    real, intent(in) :: ngas_lost
    real, intent(out) :: d_ngas_lost
    real, intent(in) :: nh4
    real, intent(out) :: d_nh4
    real, intent(in) :: no3
    real, intent(out) :: d_no3
    real, intent(in) :: psol
    real, intent(out) :: d_psol
    real, intent(in) :: plab
    real, intent(out) :: d_plab
    real, intent(in) :: pocc
    real, intent(out) :: d_pocc
    real, intent(in) :: ppar
    real, intent(out) :: d_ppar
    real, intent(in) :: input_nh4
    real, intent(in) :: input_no3
    real, intent(in) :: input_psol
    real, intent(in) :: input_ppar
    real, intent(out), dimension(n_pft) :: d_nh4_plant
    real, intent(in) :: nh4_bnf
    real, intent(out) :: d_nh4_bnf
    real, intent(out), dimension(n_pft) :: d_no3_plant
    real, intent(in) :: c_leach
    real, intent(out) :: d_c_leach
    real, intent(in) :: n_leach
    real, intent(out) :: d_n_leach
    real, intent(in) :: p_leach
    real, intent(out) :: d_p_leach
    real, intent(out), dimension(n_pft) :: d_p_plant

    integer :: ipom

    real, dimension(npom) :: decomp_pom_c
    real, dimension(npom) :: decomp_pom_n
    real, dimension(npom) :: decomp_pom_p
    real :: decomp_mom_c
    real :: decomp_mom_n
    real :: decomp_mom_p
    real :: adsorp_c
    real :: adsorp_n
    real :: adsorp_p
    real :: desorp_c
    real :: desorp_n
    real :: desorp_p
    real :: micr_uptake_c
    real :: micr_uptake_n
    real :: micr_uptake_p
    real :: dormancy_c
    real :: reactiv_c
    real :: dormancy_n
    real :: reactiv_n
    real :: dormancy_p
    real :: reactiv_p
    real :: resp_amb
    real :: resp_dmb
    real :: turnover_amb_dom_c
    real, dimension(npom) :: turnover_amb_pom_c
    real :: turnover_amb_dom_n
    real, dimension(npom) :: turnover_amb_pom_n
    real :: turnover_amb_dom_p
    real, dimension(npom) :: turnover_amb_pom_p
    real, dimension(npom) :: enz_prod_pom_c
    real :: enz_prod_mom_c
    real, dimension(npom) :: enz_prod_pom_n
    real :: enz_prod_mom_n
    real, dimension(npom) :: enz_prod_pom_p
    real :: enz_prod_mom_p
    real, dimension(npom) :: turnover_enz_pom_c
    real :: turnover_enz_mom_c
    real, dimension(npom) :: turnover_enz_pom_n
    real :: turnover_enz_mom_n
    real, dimension(npom) :: turnover_enz_pom_p
    real :: turnover_enz_mom_p
    real :: min_amb_n
    real :: min_amb_p
    real :: immob_amb_nh4
    real :: immob_amb_no3
    real :: immob_amb_psol
    real :: min_dmb_n
    real :: min_dmb_p
    real :: ovflow_amb_c
    real :: ovflow_dmb_c
    real :: nitr_rate
    real :: denitr_rate
    real :: d_p_sum
    real :: turnover_enz_ptase_c    
    real :: turnover_enz_ptase_n    
    real :: turnover_enz_ptase_p    
    real :: enz_prod_ptase_c
    real :: enz_prod_ptase_n
    real :: enz_prod_ptase_p
    real :: f_lab_p
    real :: decomp_dom_p
    real, dimension(n_pft) :: f_plant_nh4
    real, dimension(n_pft) :: f_plant_no3
    real :: occlu_flux
    real, dimension(n_pft) :: f_plant_p
    real :: bnf_rate
    real :: c_dom_leach_rate
    real :: n_dom_leach_rate
    real :: p_dom_leach_rate
    real :: f_nit_nh4
    real :: f_den_no3
    real :: no3_leach_rate
    real :: psol_leach_rate
    real :: weath_flux
    real, intent(in), dimension(n_pft) :: enz_plant_n
    real, intent(in), dimension(n_pft) :: enz_plant_p
    real, intent(in), dimension(n_pft) :: vnh4up_plant
    real, intent(in), dimension(n_pft) :: vno3up_plant
    real, intent(in), dimension(n_pft) :: vpup_plant
    real, intent(in) :: pi1
    real, intent(in) :: wfp
    integer :: ipft

    d_pom_c(:) = 0.
    d_dom_c = 0.
    d_enz_pom_c(:) = 0.
    d_mom_c = 0.
    d_qom_c = 0.
    d_enz_mom_c = 0.
    d_amb_c = 0.
    d_dmb_c = 0.
    d_enz_ptase_c = 0.

    d_pom_n(:) = 0.
    d_dom_n = 0.
    d_enz_pom_n(:) = 0.
    d_mom_n = 0.
    d_qom_n = 0.
    d_enz_mom_n = 0.
    d_amb_n = 0.
    d_dmb_n = 0.
    d_enz_ptase_n = 0.

    d_pom_p(:) = 0.
    d_dom_p = 0.
    d_enz_pom_p(:) = 0.
    d_mom_p = 0.
    d_qom_p = 0.
    d_enz_mom_p = 0.
    d_amb_p = 0.
    d_dmb_p = 0.

    d_co2_lost = 0.
    d_nmin = 0.
    d_nitr = 0.
    d_ngas_lost = 0.
    d_nh4 = 0.
    d_no3 = 0.
    d_plab = 0.
    d_pocc = 0.
    d_ppar = 0.
    d_psol = 0.
    d_enz_ptase_p = 0.
    
    d_nh4_plant = 0.
    d_nh4_bnf = 0.
    d_no3_plant = 0.
    d_c_leach = 0.
    d_n_leach = 0.
    d_p_leach = 0.
    d_p_plant = 0.


!    if(trim(consts%type) /= 'som' .and. trim(consts%type) /= 'litter')return
!    if(trim(consts%type) == 'som')return
!    if(trim(consts%type) == 'som')return

    if(trim(consts%type) /= 'som')return

    ! Inputs
    do ipom = 1, npom
       d_pom_c(ipom) = d_pom_c(ipom) + input_pom_c(ipom)
       d_pom_n(ipom) = d_pom_n(ipom) + input_pom_n(ipom)
       d_pom_p(ipom) = d_pom_p(ipom) + input_pom_p(ipom)
!print*,1,d_pom_p(ipom), input_pom_p(ipom)
    enddo

    d_dom_c = d_dom_c + input_dom_c
    d_dom_n = d_dom_n + input_dom_n
    d_dom_p = d_dom_p + input_dom_p

    d_nh4 = d_nh4 + input_nh4
    d_no3 = d_no3 + input_no3
    d_psol = d_psol + input_psol
    d_ppar = d_ppar + input_ppar

    ! Decomp
    call decomp_fluxes(npom, decomp_pom_c, decomp_pom_n, decomp_pom_p, &
         decomp_mom_c, decomp_mom_n, decomp_mom_p, &
         enz_pom_c, pom_c, pom_n, pom_p, &
         enz_mom_c, mom_c, mom_n, mom_p, &
         consts)

    do ipom = 1, npom
       d_pom_c(ipom) = d_pom_c(ipom) - decomp_pom_c(ipom)
       d_pom_n(ipom) = d_pom_n(ipom) - decomp_pom_n(ipom)
       d_pom_p(ipom) = d_pom_p(ipom) - decomp_pom_p(ipom)
!print*,2,d_pom_p(ipom), decomp_pom_p(ipom)

       d_dom_c = d_dom_c + consts%decomp_dissolved_frac * decomp_pom_c(ipom)
       d_dom_n = d_dom_n + consts%decomp_dissolved_frac * decomp_pom_n(ipom)
       d_dom_p = d_dom_p + consts%decomp_dissolved_frac * decomp_pom_p(ipom)

       d_mom_c = d_mom_c + (1.0 - consts%decomp_dissolved_frac) *  &
            decomp_pom_c(ipom)
       d_mom_n = d_mom_n + (1.0 - consts%decomp_dissolved_frac) *  &
            decomp_pom_n(ipom)
       d_mom_p = d_mom_p + (1.0 - consts%decomp_dissolved_frac) *  &
            decomp_pom_p(ipom)
    enddo
    d_mom_c = d_mom_c - decomp_mom_c
    d_mom_n = d_mom_n - decomp_mom_n
    d_mom_p = d_mom_p - decomp_mom_p

    d_dom_c = d_dom_c + decomp_mom_c
    d_dom_n = d_dom_n + decomp_mom_n
    d_dom_p = d_dom_p + decomp_mom_p

    ! Adsorption/Desorption
    call sorp_fluxes(adsorp_c, adsorp_n, adsorp_p,  &
         desorp_c, desorp_n, desorp_p, &
         dom_c, qom_c, dom_n, qom_n, dom_p, qom_p, consts)
    d_qom_c = d_qom_c + adsorp_c - desorp_c
    d_qom_n = d_qom_n + adsorp_n - desorp_n
    d_qom_p = d_qom_p + adsorp_p - desorp_p

    d_dom_c = d_dom_c - adsorp_c + desorp_c
    d_dom_n = d_dom_n - adsorp_n + desorp_n
    d_dom_p = d_dom_p - adsorp_p + desorp_p

    ! Microbial dynamics
    call micr_dyn(npom, micr_uptake_c, micr_uptake_n, micr_uptake_p, &
         dom_c, dom_n, dom_p, &
         amb_c, amb_n, amb_p,  &
         dormancy_c, dormancy_n, dormancy_p, &
         reactiv_c, reactiv_n, reactiv_p, &
         dmb_c, dmb_n, dmb_p,  &
         resp_amb, resp_dmb,  &
         turnover_amb_dom_c, turnover_amb_dom_n, turnover_amb_dom_p, &
         turnover_amb_pom_c, turnover_amb_pom_n, turnover_amb_pom_p, &
         min_amb_n, min_amb_p, min_dmb_n, min_dmb_p, &
         immob_amb_nh4, immob_amb_no3, immob_amb_psol,  &
         nh4, no3, psol, &
         ovflow_dmb_c, ovflow_amb_c, f_plant_nh4, f_plant_no3, f_plant_p, &
         f_lab_p, enz_plant_n, enz_plant_p, plab, vnh4up_plant, &
         vno3up_plant, vpup_plant,  &
         consts, bulk_den, f_nit_nh4, f_den_no3)

    ! Microbial uptake
    d_amb_c = d_amb_c + micr_uptake_c
    d_amb_n = d_amb_n + micr_uptake_n
    d_amb_p = d_amb_p + micr_uptake_p

    d_dom_c = d_dom_c - micr_uptake_c
    d_dom_n = d_dom_n - micr_uptake_n
    d_dom_p = d_dom_p - micr_uptake_p

    ! Dormancy/Reactivation
    d_amb_c = d_amb_c - dormancy_c + reactiv_c
    d_amb_n = d_amb_n - dormancy_n + reactiv_n
    d_amb_p = d_amb_p - dormancy_p + reactiv_p

    d_dmb_c = d_dmb_c + dormancy_c - reactiv_c
    d_dmb_n = d_dmb_n + dormancy_n - reactiv_n
    d_dmb_p = d_dmb_p + dormancy_p - reactiv_p

    ! Respiration
    d_amb_c = d_amb_c - resp_amb
    d_amb_c = d_amb_c - ovflow_amb_c

    d_co2_lost = d_co2_lost + resp_amb
    d_co2_lost = d_co2_lost + ovflow_amb_c

    d_dmb_c = d_dmb_c - resp_dmb
    d_dmb_c = d_dmb_c - ovflow_dmb_c

    d_co2_lost = d_co2_lost + resp_dmb
    d_co2_lost = d_co2_lost + ovflow_dmb_c

    ! Turnover
    d_amb_c = d_amb_c - turnover_amb_dom_c
    d_amb_n = d_amb_n - turnover_amb_dom_n
    d_amb_p = d_amb_p - turnover_amb_dom_p

    d_dom_c = d_dom_c + turnover_amb_dom_c
    d_dom_n = d_dom_n + turnover_amb_dom_n
    d_dom_p = d_dom_p + turnover_amb_dom_p

    do ipom = 1, npom
       d_amb_c = d_amb_c - turnover_amb_pom_c(ipom)
       d_amb_n = d_amb_n - turnover_amb_pom_n(ipom)
       d_amb_p = d_amb_p - turnover_amb_pom_p(ipom)

       d_pom_c(ipom) = d_pom_c(ipom) + turnover_amb_pom_c(ipom)
       d_pom_n(ipom) = d_pom_n(ipom) + turnover_amb_pom_n(ipom)
       d_pom_p(ipom) = d_pom_p(ipom) + turnover_amb_pom_p(ipom)
    enddo

    ! Mineralization
    d_amb_n = d_amb_n - min_amb_n
    d_amb_p = d_amb_p - min_amb_p

    d_nh4 = d_nh4 + min_amb_n
    d_psol = d_psol + min_amb_p

    d_dmb_n = d_dmb_n - min_dmb_n
    d_dmb_p = d_dmb_p - min_dmb_p

    d_nh4 = d_nh4 + min_dmb_n
    d_psol = d_psol + min_dmb_p

    d_nmin = d_nmin + min_dmb_n + min_amb_n

    ! Immobilization
    d_amb_n = d_amb_n + immob_amb_nh4
    d_amb_n = d_amb_n + immob_amb_no3
    d_amb_p = d_amb_p + immob_amb_psol

    d_nh4 = d_nh4 - immob_amb_nh4
    d_no3 = d_no3 - immob_amb_no3
    d_psol = d_psol - immob_amb_psol

    d_nmin = d_nmin - immob_amb_nh4 - immob_amb_no3

    ! Plant uptake
    do ipft = 1, n_pft
       d_nh4_plant(ipft) = d_nh4_plant(ipft) + f_plant_nh4(ipft)
       d_no3_plant(ipft) = d_no3_plant(ipft) + f_plant_no3(ipft)
       d_p_plant(ipft) = d_p_plant(ipft) + f_plant_p(ipft)
       
       d_nh4 = d_nh4 - f_plant_nh4(ipft)
       d_no3 = d_no3 - f_plant_no3(ipft)
       d_psol = d_psol - f_plant_p(ipft)
    enddo

    ! Enzyme kinetics
    call enz_kin(npom, enz_prod_pom_c, enz_prod_pom_n, enz_prod_pom_p, &
         enz_prod_mom_c, enz_prod_mom_n, enz_prod_mom_p, &
         pom_c, amb_c, amb_n, amb_p, &
         turnover_enz_pom_c, turnover_enz_pom_n, turnover_enz_pom_p, &
         turnover_enz_mom_c, turnover_enz_mom_n, turnover_enz_mom_p, &
         enz_pom_c, enz_mom_c, &
         enz_prod_ptase_c, enz_prod_ptase_n, enz_prod_ptase_p, &
         turnover_enz_ptase_c, turnover_enz_ptase_n, turnover_enz_ptase_p, &
         enz_ptase_c, &
         consts)

    do ipom = 1, npom
       d_enz_pom_c(ipom) = d_enz_pom_c(ipom) + enz_prod_pom_c(ipom)

       d_amb_c = d_amb_c - enz_prod_pom_c(ipom)
       d_amb_n = d_amb_n - enz_prod_pom_n(ipom)
       d_amb_p = d_amb_p - enz_prod_pom_p(ipom)

       d_enz_pom_c(ipom) = d_enz_pom_c(ipom) - turnover_enz_pom_c(ipom)

       d_dom_c = d_dom_c + turnover_enz_pom_c(ipom)
       d_dom_n = d_dom_n + turnover_enz_pom_n(ipom)
       d_dom_p = d_dom_p + turnover_enz_pom_p(ipom)
    enddo

    d_enz_mom_c = d_enz_mom_c + enz_prod_mom_c

    d_amb_c = d_amb_c - enz_prod_mom_c
    d_amb_n = d_amb_n - enz_prod_mom_n
    d_amb_p = d_amb_p - enz_prod_mom_p

    d_enz_mom_c = d_enz_mom_c - turnover_enz_mom_c

    d_dom_c = d_dom_c + turnover_enz_mom_c
    d_dom_n = d_dom_n + turnover_enz_mom_n
    d_dom_p = d_dom_p + turnover_enz_mom_p

    d_enz_ptase_c = d_enz_ptase_c + enz_prod_ptase_c

    d_amb_c = d_amb_c - enz_prod_ptase_c
    d_amb_n = d_amb_n - enz_prod_ptase_n
    d_amb_p = d_amb_p - enz_prod_ptase_p

    d_enz_ptase_c = d_enz_ptase_c - turnover_enz_ptase_c

    d_dom_c = d_dom_c + turnover_enz_ptase_c
    d_dom_n = d_dom_n + turnover_enz_ptase_n
    d_dom_p = d_dom_p + turnover_enz_ptase_p

    ! Nitrification / Denitrification
    call nitr_denitr(nh4, no3, nitr_rate, denitr_rate, consts, d_co2_lost, &
         pi1, wfp, f_nit_nh4, f_den_no3)

    if(no3 < 1.0e-12)denitr_rate = 0.

    d_nh4 = d_nh4 - nitr_rate
    d_no3 = d_no3 + nitr_rate * 0.9 - denitr_rate

    d_ngas_lost = d_ngas_lost + denitr_rate + nitr_rate * 0.1

    d_nitr = d_nitr + nitr_rate

    call phosphatase_dyn(npom, pom_p, dom_p, consts, decomp_pom_p, &
       decomp_dom_p, enz_ptase_c)

    do ipom = 1, npom
       d_psol = d_psol + decomp_pom_p(ipom)
       d_pom_p(ipom) = d_pom_p(ipom) - decomp_pom_p(ipom)
    enddo
    d_psol = d_psol + decomp_dom_p
    d_dom_p = d_dom_p - decomp_dom_p

    bnf_rate = 1.2e-10 ! 5.e-11 ! mgN/gSOIL/s
    d_nh4_bnf = d_nh4_bnf + bnf_rate
    d_nh4 = d_nh4 + bnf_rate

    c_dom_leach_rate = consts%dom_leach_effic * dom_c * water_drainage
    n_dom_leach_rate = consts%dom_leach_effic * dom_n * water_drainage
    p_dom_leach_rate = consts%dom_leach_effic * dom_p * water_drainage
    no3_leach_rate = consts%no3_leach_effic * no3 * water_drainage
    psol_leach_rate = consts%psol_leach_effic * psol * water_drainage
    
    d_c_leach = d_c_leach + c_dom_leach_rate
    d_n_leach = d_n_leach + n_dom_leach_rate + no3_leach_rate
    d_p_leach = d_p_leach + p_dom_leach_rate + psol_leach_rate

    d_dom_c = d_dom_c - c_dom_leach_rate
    d_dom_n = d_dom_n - n_dom_leach_rate
    d_dom_p = d_dom_p - p_dom_leach_rate


    if(no3 < 1.0e-12)no3_leach_rate = 0.
    d_no3 = d_no3 - no3_leach_rate
    d_psol = d_psol - psol_leach_rate

    return
  end subroutine mend_derivs_layer

  subroutine phosphatase_dyn(npom, pom_p, dom_p, consts, decomp_pom_p, &
       decomp_dom_p, enz_ptase_c)
    implicit none
    type(decomp_consts) :: consts
    real :: denom
    integer, intent(in) :: npom
    real, intent(in), dimension(npom) :: pom_p
    real, intent(in) :: dom_p
    integer :: ipom
    real, intent(out) :: decomp_dom_p
    real, intent(out), dimension(npom) :: decomp_pom_p
    real, intent(in) :: enz_ptase_c


    !! TURNING OFF PTASE

!    denom = 1. + dom_p / consts%kmm_decomp_ptase_dom
    denom = 1.
    do ipom = 1, npom
       denom = denom + pom_p(ipom) / consts%kmm_decomp_ptase_pom(ipom)
    enddo

    do ipom = 1, npom
       decomp_pom_p(ipom) = consts%vmax_decomp_ptase_pom(ipom) *  &
            enz_ptase_c * pom_p(ipom) /   &
            (consts%kmm_decomp_ptase_pom(ipom) * denom)
    enddo

!    decomp_dom_p = consts%vmax_decomp_ptase_dom *  &
!         enz_ptase_c * dom_p /   &
!         (consts%kmm_decomp_ptase_dom * denom)
    decomp_dom_p = 0.

    return
  end subroutine phosphatase_dyn

  subroutine nitr_denitr(nh4, no3, nitr, denitr, consts, d_co2_lost, pi1, &
       wfp, f_nit_nh4, f_den_no3)
    implicit none
    
    real, intent(in) :: f_nit_nh4
    real, intent(in) :: f_den_no3
    type(decomp_consts) :: consts
    real, intent(in) :: nh4
    real, intent(in) :: no3
    real, intent(out) :: nitr
    real, intent(out) :: denitr
    real :: phi
    real :: delGrosso_M
    real :: delGrosso_A
    real, intent(in) :: d_co2_lost
    real, intent(in) :: pi1
    real, intent(in) :: wfp
    real :: water_function
    real :: no3_function
    real :: co2_function

    phi = 0.1  ! * 0. ! was just 0.1
    nitr = 2. * consts%v_nitr * nh4 !* f_nit_nh4
    denitr = consts%v_denitr * no3 + nitr * phi

!   Potter et al. 1996
!   P = slmsts - sfldcap
!   y = 0.477 * P**3 - 0.596 * P**2 + 0.437 * P + 0.564
!   D/D0 = P**y
    delGrosso_M = min(0.113,consts%v_denitr)*(-3.05)+0.36
    delGrosso_A = 9.0 - delGrosso_M * (d_co2_lost*86400.*1000.)
    water_function = 0.5 + atan(0.6*pi1*(0.1*(wfp*100.)-delGrosso_A))/pi1
    no3_function = 0.001 * (1.15 * (1000.*no3)**0.57) / 86400.
    co2_function = 0.001 * (0.1 * (d_co2_lost*86400.*1000.)**1.3 ) / 86400.
    denitr = min(no3_function,co2_function) * water_function !* f_den_no3

    return
  end subroutine nitr_denitr

  subroutine enz_kin(npom, enz_prod_pom_c, enz_prod_pom_n, enz_prod_pom_p, &
         enz_prod_mom_c, enz_prod_mom_n, enz_prod_mom_p, &
         pom_c, amb_c, amb_n, amb_p, &
         turnover_enz_pom_c, turnover_enz_pom_n, turnover_enz_pom_p, &
         turnover_enz_mom_c, turnover_enz_mom_n, turnover_enz_mom_p, &
         enz_pom_c, enz_mom_c, &
         enz_prod_ptase_c, enz_prod_ptase_n, enz_prod_ptase_p, &
         turnover_enz_ptase_c, turnover_enz_ptase_n, turnover_enz_ptase_p, &
         enz_ptase_c, &
         consts)
    implicit none

    type(decomp_consts) :: consts
    integer, intent(in) :: npom
    real, intent(in) :: amb_n
    real, intent(in) :: amb_p
    real, intent(out), dimension(npom) :: enz_prod_pom_c
    real, intent(out) :: enz_prod_mom_c
    real, intent(out), dimension(npom) :: enz_prod_pom_n
    real, intent(out) :: enz_prod_mom_n
    real, intent(out), dimension(npom) :: enz_prod_pom_p
    real, intent(out) :: enz_prod_mom_p
    integer :: ipom
    real, intent(in), dimension(npom) :: pom_c
    real, intent(in) :: amb_c
    real, intent(out), dimension(npom) :: turnover_enz_pom_c
    real, intent(out) :: turnover_enz_mom_c
    real, intent(out), dimension(npom) :: turnover_enz_pom_n
    real, intent(out) :: turnover_enz_mom_n
    real, intent(out), dimension(npom) :: turnover_enz_pom_p
    real, intent(out) :: turnover_enz_mom_p
    real, intent(in), dimension(npom) :: enz_pom_c
    real, intent(in) :: enz_mom_c

    real, intent(out) :: enz_prod_ptase_c
    real, intent(out) :: enz_prod_ptase_n
    real, intent(out) :: enz_prod_ptase_p
    real, intent(out) :: turnover_enz_ptase_c
    real, intent(out) :: turnover_enz_ptase_n
    real, intent(out) :: turnover_enz_ptase_p
    real, intent(in) :: enz_ptase_c


    do ipom = 1, npom
! CHECK THIS: ENZYME PRODUCTION COSTS A LOT OF N.  TURN IT OFF IF N-LIMITED?
       enz_prod_pom_c(ipom) = pom_c(ipom) / sum(pom_c) *  &
            consts%prod_frac_enz_pom * &
            consts%spec_maint_rate * amb_c
       if(amb_c/amb_n > consts%mb_cn_max)enz_prod_pom_c(ipom) = 0.
       enz_prod_pom_n(ipom) = enz_prod_pom_c(ipom) / consts%enz_pom_c2n(ipom)
       enz_prod_pom_p(ipom) = enz_prod_pom_c(ipom) / consts%enz_pom_c2p(ipom)

       turnover_enz_pom_c(ipom) = consts%enz_turnover_rate * enz_pom_c(ipom)
       turnover_enz_pom_n(ipom) = turnover_enz_pom_c(ipom) /  &
            consts%enz_pom_c2n(ipom)
       turnover_enz_pom_p(ipom) = turnover_enz_pom_c(ipom) /  &
            consts%enz_pom_c2p(ipom)
    enddo

    enz_prod_mom_c = consts%prod_frac_enz_mom * consts%spec_maint_rate * amb_c
    if(amb_c/amb_n > consts%mb_cn_max)enz_prod_mom_c = 0.
    enz_prod_mom_n = enz_prod_mom_c / consts%enz_mom_c2n
    enz_prod_mom_p = enz_prod_mom_c / consts%enz_mom_c2p

    turnover_enz_mom_c = consts%enz_turnover_rate * enz_mom_c
    turnover_enz_mom_n = turnover_enz_mom_c / consts%enz_mom_c2n
    turnover_enz_mom_p = turnover_enz_mom_c / consts%enz_mom_c2p
    
    enz_prod_ptase_c = consts%prod_frac_enz_ptase * &
         consts%spec_maint_rate * amb_c
    if(amb_c / amb_p < consts%mb_cp_max .and.   &
         amb_c / amb_n > consts%mb_cn_max)enz_prod_ptase_c = 0.
!    if(amb_c / amb_p > consts%mb_cp_max .and.   &
!         amb_c / amb_n < consts%mb_cn_max)enz_prod_ptase_c = 0.
    enz_prod_ptase_n = enz_prod_ptase_c / consts%enz_ptase_c2n
    enz_prod_ptase_p = enz_prod_ptase_c / consts%enz_ptase_c2p

    turnover_enz_ptase_c = consts%enz_turnover_rate * enz_ptase_c
    turnover_enz_ptase_n = turnover_enz_ptase_c / consts%enz_ptase_c2n
    turnover_enz_ptase_p = turnover_enz_ptase_c / consts%enz_ptase_c2p

    return
  end subroutine enz_kin

  subroutine micr_dyn(npom, micr_uptake_c, micr_uptake_n, micr_uptake_p, &
         dom_c, dom_n, dom_p, &
         amb_c, amb_n, amb_p,  &
         dormancy_c, dormancy_n, dormancy_p, &
         reactiv_c, reactiv_n, reactiv_p, &
         dmb_c, dmb_n, dmb_p,  &
         resp_amb, resp_dmb,  &
         turnover_amb_dom_c, turnover_amb_dom_n, turnover_amb_dom_p, &
         turnover_amb_pom_c, turnover_amb_pom_n, turnover_amb_pom_p, &
         min_amb_n, min_amb_p, min_dmb_n, min_dmb_p, &
         immob_amb_nh4, immob_amb_no3, immob_amb_psol,  &
         nh4, no3, psol, &
         ovflow_dmb_c, ovflow_amb_c, f_plant_nh4, f_plant_no3, f_plant_p, &
         f_lab_p, enz_plant_n, enz_plant_p, plab, vnh4up_plant,  &
         vno3up_plant, vpup_plant,  &
         consts, bulk_den, f_nit_nh4, f_den_no3)

    implicit none
    
    integer :: ipft
    real, intent(in) :: bulk_den
    type(decomp_consts) :: consts
    integer, intent(in) :: npom
    real, intent(out) :: micr_uptake_c
    real, intent(out) :: micr_uptake_n
    real, intent(out) :: micr_uptake_p
    real, intent(in) :: amb_c
    real, intent(in) :: amb_n
    real, intent(in) :: amb_p
    real, intent(in) :: dmb_c
    real, intent(in) :: dmb_n
    real, intent(in) :: dmb_p
    real, intent(in) :: dom_c
    real, intent(in) :: dom_n
    real, intent(in) :: dom_p
    real, intent(out) :: dormancy_c
    real, intent(out) :: reactiv_c
    real, intent(out) :: dormancy_n
    real, intent(out) :: reactiv_n
    real, intent(out) :: dormancy_p
    real, intent(out) :: reactiv_p
    real, intent(out) :: resp_amb
    real, intent(out) :: resp_dmb
    real, intent(out) :: turnover_amb_dom_c
    real, intent(in) :: plab
    real, dimension(npom), intent(out) :: turnover_amb_pom_c
    real, intent(out) :: turnover_amb_dom_n
    real, dimension(npom), intent(out) :: turnover_amb_pom_n
    real, intent(out) :: turnover_amb_dom_p
    real, dimension(npom), intent(out) :: turnover_amb_pom_p
    real, intent(out) :: f_lab_p
    real :: turnover_amb
    real, intent(out) :: min_amb_n
    real, intent(out) :: min_amb_p
    integer :: ipom
    real :: phi
    real :: ygn
    real :: vnup_amb
    real :: ygp
    real :: vpup_amb
    real :: rtp1
    real, intent(out) :: immob_amb_no3
    real, intent(out) :: immob_amb_nh4
    real, intent(out) :: immob_amb_psol
    real, intent(in) :: no3
    real, intent(in) :: nh4
    real, intent(in) :: psol
    real :: enz_lab_p
    real :: cn_amb
    real :: cn_dmb
    real :: cp_amb
    real :: cp_dmb
    real :: p_target
    real :: n_target
    real :: c_target
    real, intent(out) :: min_dmb_n
    real, intent(out) :: min_dmb_p
    real, intent(out) :: ovflow_dmb_c
    real, intent(out) :: ovflow_amb_c
    real :: f_immob_nh4, f_immob_no3
    real, intent(out) :: f_nit_nh4, f_den_no3
    real, intent(out), dimension(n_pft) :: f_plant_nh4
    real, intent(out), dimension(n_pft) :: f_plant_no3
    real, intent(out), dimension(n_pft) :: f_plant_p
    real, intent(in), dimension(n_pft) :: enz_plant_n
    real, intent(in), dimension(n_pft) :: enz_plant_p
    real, intent(in), dimension(n_pft) :: vnh4up_plant
    real, intent(in), dimension(n_pft) :: vno3up_plant
    real, intent(in), dimension(n_pft) :: vpup_plant
    real :: no3_eff, psol_eff
    real :: tot_enz_plant_n_nh4
    real :: tot_enz_plant_n_no3
    real :: tot_enz_plant_p

    tot_enz_plant_n_nh4 = sum(enz_plant_n/consts%ksnh4_plant)*bulk_den
    tot_enz_plant_n_no3 = sum(enz_plant_n/consts%ksno3_plant)*bulk_den
    tot_enz_plant_p = sum(enz_plant_p/consts%kspsol_plant)*bulk_den

    no3_eff = no3
    if(no3 < 1.0e-12)then
       no3_eff = 0.
    endif

    psol_eff = psol
    if(psol < 1.0e-12)then
       psol_eff = 0.
    endif

! NEED TO USE 2/3 POWER, CORRESPONDING TO SURFACE AREA.
! OTHERWISE, FUNDAMENTALLY UNSTABLE.  IT WOULD JUST BE EXPONENTIAL.
    micr_uptake_c = 1. / consts%growth_yield * (consts%Vmax_amb +  &
         consts%spec_maint_rate) * &
         dom_c * consts%amb_size_scale *  &
!         (amb_c / consts%amb_size_scale) /  &
         (amb_c / consts%amb_size_scale)**0.67 /  &
         (consts%kmm_amb + dom_c)
    micr_uptake_n = micr_uptake_c * dom_n / dom_c
    micr_uptake_p = micr_uptake_c * dom_p / dom_c

    dormancy_c = (1.0 - dom_c / (consts%kmm_amb + dom_c)) *   &
         consts%vmax_a2d * amb_c
    dormancy_n = dormancy_c * amb_n / amb_c
    dormancy_p = dormancy_c * amb_p / amb_c

    reactiv_c = dom_c / (consts%kmm_amb + dom_c) *   &
         consts%vmax_d2a * dmb_c
    reactiv_n = reactiv_c * dmb_n / dmb_c 
    reactiv_p = reactiv_c * dmb_p / dmb_c 

    resp_amb = (1. / consts%growth_yield - 1.) * (consts%Vmax_amb +   &
         consts%spec_maint_rate) * consts%amb_size_scale *  &
!         (amb_c / consts%amb_size_scale) * dom_c /  &
         (amb_c / consts%amb_size_scale)**0.67 * dom_c /  &
         (consts%kmm_amb + dom_c)

    resp_dmb = consts%dorm_resp_frac * consts%spec_maint_rate * dmb_c

    turnover_amb = consts%amb_turnover_rate * amb_c
    turnover_amb_dom_c = consts%turnover_frac_dom * turnover_amb
    turnover_amb_dom_n = turnover_amb_dom_c * amb_n / amb_c
    turnover_amb_dom_p = turnover_amb_dom_c * amb_p / amb_c
    
    do ipom = 1, npom
       turnover_amb_pom_c(ipom) = consts%turnover_frac_pom(ipom) * turnover_amb
       turnover_amb_pom_n(ipom) = turnover_amb_pom_c(ipom) * amb_n / amb_c
       turnover_amb_pom_p(ipom) = turnover_amb_pom_c(ipom) * amb_p / amb_c
    enddo

    ! N Mineralization
    phi = ((amb_c / amb_n - consts%mb_cn_min) / (consts%mb_cn_max -   &
         consts%mb_cn_min))**consts%wexp

    phi = min(max(0., phi), 1.)
    ygn = consts%ygn_base * phi

!    phi = max(0., phi)
!    ygn = consts%ygn_base * min(phi,1.)

    vnup_amb = consts%vnup_amb_base * phi
! CHECK PARAMETERS HERE.  THIS IS ADDED:
!    min_amb_n = (1. - ygn) * micr_uptake_c
!    min_amb_n = min((1. - ygn) * micr_uptake_c, micr_uptake_n)
    min_amb_n = (1. - ygn) * micr_uptake_n

    ! P mineralization
    phi = ((amb_c / amb_p - consts%mb_cp_min) / (consts%mb_cp_max -   &
         consts%mb_cp_min))**consts%wexp
    phi = min(max(0., phi), 1.)
    ygp = consts%ygp_base * phi
    vpup_amb = consts%vpup_amb_base * phi
! CHECK PARAMETERS HERE.  THIS IS ADDED:
!    min_amb_p = (1. - ygp) * micr_uptake_c
!    min_amb_p = min((1. - ygp) * micr_uptake_c, micr_uptake_p)
    min_amb_p = (1. - ygp) * micr_uptake_p

    ! N Immobilization:
! CHECK: SCALE WITH AMB_N OR AMB_C?
! SCALE WITH AMB TO SOME POWER?


    rtp1 =  1. + nh4 / (consts%ksnh4_amb / bulk_den) +   &
         no3_eff / (consts%ksno3_amb / bulk_den) +   &
         consts%micr_bio2enz * amb_c * 14. / (consts%ksnh4_amb / bulk_den) +   &
         tot_enz_plant_n_nh4 + &
         consts%enz_nit_n / consts%ksnh4_nit
    immob_amb_nh4 = vnup_amb * nh4 * consts%micr_bio2enz * amb_c * 14. /  &
         ((consts%ksnh4_amb / bulk_den) * rtp1) 

    rtp1 =  1. + nh4 / (consts%ksnh4_amb / bulk_den) +   &
         no3_eff / (consts%ksno3_amb / bulk_den) +   &
         consts%micr_bio2enz * amb_c * 14. / (consts%ksno3_amb / bulk_den) +  &
         tot_enz_plant_n_no3 + &
         consts%enz_den_n / consts%ksno3_den
    immob_amb_no3 = vnup_amb * no3_eff * consts%micr_bio2enz * amb_c * 14. /   &
         ((consts%ksno3_amb / bulk_den) * rtp1) 

    do ipft = 1, n_pft
       rtp1 =  1. + nh4 / (consts%ksnh4_plant(ipft) / bulk_den) +  &
            no3_eff / (consts%ksno3_plant(ipft) / bulk_den) +   &
            consts%micr_bio2enz * amb_c * 14. / (consts%ksnh4_amb / bulk_den) +  &
            tot_enz_plant_n_nh4 + &
            consts%enz_nit_n / consts%ksnh4_nit
       f_plant_nh4(ipft) = vnh4up_plant(ipft) * nh4 /   &
            ((consts%ksnh4_plant(ipft) / bulk_den) * rtp1)
!       f_plant_nh4(ipft) = vnh4up_plant(ipft) * nh4 * enz_plant_n(ipft) /   &
!            ((consts%ksnh4_plant(ipft) / bulk_den) * rtp1)
       
       rtp1 =  1. + nh4 / (consts%ksnh4_plant(ipft) / bulk_den) +  &
            no3_eff / (consts%ksno3_plant(ipft) / bulk_den) +   &
            consts%micr_bio2enz * amb_c * 14. / (consts%ksno3_amb / bulk_den) +  &
            tot_enz_plant_n_no3 + &
            consts%enz_den_n / consts%ksno3_den
       f_plant_no3(ipft) = vno3up_plant(ipft) * no3_eff /   &
            ((consts%ksno3_plant(ipft) / bulk_den) * rtp1)
!       f_plant_no3(ipft) = vno3up_plant(ipft) * no3_eff * enz_plant_n(ipft) /   &
!            ((consts%ksno3_plant(ipft) / bulk_den) * rtp1)
    enddo

    rtp1 =  1. +  &
         nh4 / (consts%ksnh4_nit / consts%eff_soil_depth / bulk_den) + &
         consts%micr_bio2enz * amb_c * 14. / (consts%ksnh4_amb / bulk_den) +   &
         tot_enz_plant_n_nh4 + &
         consts%enz_nit_n / consts%ksnh4_nit
    f_nit_nh4 = nh4 /  &
         (consts%ksnh4_nit / consts%eff_soil_depth / bulk_den * rtp1)

    rtp1 =  1. +  &
         no3_eff / (consts%ksno3_den / consts%eff_soil_depth / bulk_den) + &
         consts%micr_bio2enz * amb_c * 14.  / (consts%ksno3_amb / bulk_den) +   &
         tot_enz_plant_n_no3 + &
         consts%enz_den_n / consts%ksno3_den
    f_den_no3 = no3_eff /  &
         (consts%ksno3_den / consts%eff_soil_depth / bulk_den * rtp1)

    ! P immobilization
    ! Surfaces included
    enz_lab_p = max(0.,consts%vmax_surf / bulk_den - plab)

    do ipft = 1, n_pft
       rtp1 = 1.0 + (plab+psol_eff) / (consts%kspsol_plant(ipft) / bulk_den) +   &
            tot_enz_plant_p  +  &
            consts%micr_bio2enz * amb_c * 31.  / (consts%kspsol_amb / bulk_den) + &
            (consts%vmax_surf/bulk_den) / (consts%kspsol_lab / bulk_den)
       f_plant_p(ipft) = vpup_plant(ipft) * (plab+psol_eff) /  &
            ((consts%kspsol_plant(ipft) / bulk_den) * rtp1)
!       f_plant_p(ipft) = vpup_plant(ipft) * (plab+psol_eff) * enz_plant_p(ipft) /  &
!            ((consts%kspsol_plant(ipft) / bulk_den) * rtp1)
    enddo

    rtp1 = 1.0 + (plab+psol_eff) / (consts%kspsol_amb / bulk_den) +  &
         tot_enz_plant_p +   &
         consts%micr_bio2enz * amb_c * 31. / (consts%kspsol_amb / bulk_den) + &
         (consts%vmax_surf/bulk_den) / (consts%kspsol_lab / bulk_den)
    immob_amb_psol = vpup_amb * (plab+psol_eff) * consts%micr_bio2enz * amb_c * 31. /  &
         ((consts%kspsol_amb / bulk_den) * rtp1)

!    rtp1 = 1.0 + (plab+psol_eff) / (consts%kspsol_lab / bulk_den) +   &
!         enz_plant_p / (consts%kspsol_plant / bulk_den) +   &
!         amb_p / (consts%kspsol_amb / bulk_den) + &
!         (consts%vmax_surf/bulk_den) / (consts%kspsol_lab / bulk_den)
!    f_lab_p = psol_eff / (consts%kspsol_lab / bulk_den * rtp1) *  &
!         consts%vmax_surf / bulk_den * consts%kspsol_lab / bulk_den / &
!         (consts%kspsol_lab / bulk_den + psol_eff)**2
    rtp1 = consts%kspsol_lab / bulk_den * (1.0 +  &
         tot_enz_plant_p  +   &
         consts%micr_bio2enz * amb_c * 31. / (consts%kspsol_amb / bulk_den) + &
         (consts%vmax_surf/bulk_den) / (consts%kspsol_lab / bulk_den))
    f_lab_p = rtp1 * consts%vmax_surf/bulk_den /  &
         (rtp1 + plab + psol_eff)**2
!    if(plab < 1.0e-10)f_lab_p = 0.

    ovflow_dmb_c = 0.
    ovflow_amb_c = 0.
    min_dmb_p = 0.
    min_dmb_n = 0.
    ! Stoichiometric bounds: N
    cn_dmb = dmb_c / dmb_n
    if(cn_dmb < consts%mb_cn_min)then
       n_target = dmb_c / consts%mb_cn_min
       min_dmb_n = consts%stoich_adj_rate * (dmb_n - n_target)
    elseif(cn_dmb > consts%mb_cn_max)then
       c_target = dmb_n * consts%mb_cn_max
       ovflow_dmb_c = consts%stoich_adj_rate * (dmb_c - c_target)
    endif

    cn_amb = amb_c / amb_n
    if(cn_amb < consts%mb_cn_min)then
       n_target = amb_c / consts%mb_cn_min
       min_amb_n = min_amb_n + consts%stoich_adj_rate * (amb_n - n_target)
    elseif(cn_amb > consts%mb_cn_max)then
       c_target = amb_n * consts%mb_cn_max
       ovflow_amb_c = consts%stoich_adj_rate * (amb_c - c_target)
    endif

    ! Stoichiometric bounds: P
    cp_dmb = dmb_c / dmb_p
    if(cp_dmb < consts%mb_cp_min)then
       p_target = dmb_c / consts%mb_cp_min
       min_dmb_p = consts%stoich_adj_rate * (dmb_p - p_target)
    elseif(cp_dmb > consts%mb_cp_max)then
       c_target = dmb_p * consts%mb_cp_max
       ovflow_dmb_c = ovflow_dmb_c + consts%stoich_adj_rate *  &
            (dmb_c - c_target)
    endif

    cp_amb = amb_c / amb_p
    if(cp_amb < consts%mb_cp_min)then
       p_target = amb_c / consts%mb_cp_min
       min_amb_p = min_amb_p + consts%stoich_adj_rate * (amb_p - p_target)
    elseif(cp_amb > consts%mb_cp_max)then
       c_target = amb_p * consts%mb_cp_max
       ovflow_amb_c = ovflow_amb_c + consts%stoich_adj_rate *  &
            (amb_c - c_target)
    endif

    return
  end subroutine micr_dyn

  subroutine sorp_fluxes(adsorp_c, adsorp_n, adsorp_p, &
       desorp_c, desorp_n, desorp_p, &
       dom_c, qom_c, dom_n, qom_n, dom_p, qom_p, consts)
    implicit none

    type(decomp_consts) :: consts
    real, intent(out) :: adsorp_c
    real, intent(out) :: desorp_c
    real, intent(out) :: adsorp_n
    real, intent(out) :: desorp_n
    real, intent(out) :: adsorp_p
    real, intent(out) :: desorp_p
    real, intent(in) :: qom_c
    real, intent(in) :: dom_c
    real, intent(in) :: qom_n
    real, intent(in) :: dom_n
    real, intent(in) :: qom_p
    real, intent(in) :: dom_p

    adsorp_c = consts%sorp_rate *   &
         (1.0 - qom_c / consts%sorp_cap) * dom_c
    adsorp_n = adsorp_c * dom_n / dom_c
    adsorp_p = adsorp_c * dom_p / dom_c

    desorp_c = consts%desorp_rate * qom_c / consts%sorp_cap
    desorp_n = desorp_c * qom_n / qom_c
    desorp_p = desorp_c * qom_p / qom_c

    return
  end subroutine sorp_fluxes

  subroutine decomp_fluxes(npom, decomp_pom_c, decomp_pom_n, decomp_pom_p, &
       decomp_mom_c, decomp_mom_n, decomp_mom_p, &
       enz_pom_c, pom_c, pom_n, pom_p, &
       enz_mom_c, mom_c, mom_n, mom_p, &
       consts)

    implicit none

    type(decomp_consts) :: consts
    integer, intent(in) :: npom
    real, intent(out), dimension(npom) :: decomp_pom_c
    real, intent(out), dimension(npom) :: decomp_pom_n
    real, intent(out), dimension(npom) :: decomp_pom_p
    real, intent(in), dimension(npom) :: enz_pom_c
    real, intent(in), dimension(npom) :: pom_c
    real, intent(in), dimension(npom) :: pom_n
    real, intent(in), dimension(npom) :: pom_p
    real, intent(in) :: enz_mom_c
    integer :: ipom
    real, intent(out) :: decomp_mom_c
    real, intent(out) :: decomp_mom_n
    real, intent(out) :: decomp_mom_p
    real, intent(in) :: mom_c
    real, intent(in) :: mom_n
    real, intent(in) :: mom_p

    do ipom = 1, npom
       decomp_pom_c(ipom) = consts%vmax_decomp_pom(ipom) * enz_pom_c(ipom) *  &
            pom_c(ipom) / (consts%kmm_decomp_pom(ipom) + pom_c(ipom))
       decomp_pom_n(ipom) = decomp_pom_c(ipom) * pom_n(ipom) / pom_c(ipom)
       decomp_pom_p(ipom) = decomp_pom_c(ipom) * pom_p(ipom) / pom_c(ipom)
    enddo

    decomp_mom_c = consts%vmax_decomp_mom * enz_mom_c * mom_c / &
         (consts%kmm_decomp_mom + mom_c)
    decomp_mom_n = decomp_mom_c * mom_n / mom_c
    decomp_mom_p = decomp_mom_c * mom_p / mom_c

    return
  end subroutine decomp_fluxes
  
!   subroutine ncom(amb_n, amb_p, plant_n, plant_p, nh4, no3, pox,   &
!        f_immob_nh4, f_nit_nh4,  &
!        f_immob_no3, f_den_no3, f_plant_nh4, f_plant_no3, f_plant_p, &
!        f_immob_p, f_surf_p)
!     use ncom_consts
!     implicit none

!     real, intent(in) :: amb_n
!     real, intent(in) :: amb_p
!     real, intent(in) :: plant_n ! mgN/gSOIL
!     real, intent(in) :: plant_p ! mgP/gSOIL
!     real, intent(in) :: nh4 ! mgN/gSOIL
!     real, intent(in) :: no3 ! mgN/gSOIL
!     real, intent(in) :: pox ! mgP/gSOIL
!     real, intent(out) :: f_plant_nh4  !mgN/gSOIL
!     real, intent(out) :: f_plant_no3  !mgN/gSOIL
!     real, intent(out) :: f_immob_nh4  !mgN/gSOIL
!     real, intent(out) :: f_nit_nh4  !mgN/gSOIL
!     real, intent(out) :: f_immob_no3  !mgN/gSOIL
!     real, intent(out) :: f_den_no3  !mgN/gSOIL
!     real, intent(out) :: f_plant_p  !mgP/gSOIL
!     real, intent(out) :: f_immob_p  !mgP/gSOIL
!     real, intent(out) :: f_surf_p  !mgP/gSOIL

!     real :: enz_mic_n  !  mgN/gSOIL
!     real :: enz_mic_p  !  mgP/gSOIL
!     real :: enz_plant_n  !  mgN/gSOIL
!     real :: enz_plant_p  !  mgP/gSOIL
!     real :: denom ! unitless

!     enz_mic_n = enz_biomass_scale_mic * amb_n
!     enz_mic_p = enz_biomass_scale_mic * amb_p
!     enz_plant_n = enz_biomass_scale_plant * plant_n
!     enz_plant_p = enz_biomass_scale_plant * plant_p

!     denom = 1. + nh4 / kmm_plant_nh4 + no3 / kmm_plant_no3 + enz_plant_n / &
!          kmm_plant_nh4 + enz_mic_n / kmm_mic_nh4 + enz_nit_n / kmm_nit_nh4
!     f_plant_nh4 = k_plant_nh4 * nh4 * enz_plant_n / (kmm_plant_nh4 * denom)

!     denom = 1.0 + nh4 / kmm_mic_nh4 + no3 / kmm_mic_no3 + enz_plant_n /  &
!          kmm_plant_nh4 + enz_mic_n / kmm_mic_nh4 + enz_nit_n / kmm_nit_nh4
!     f_immob_nh4 = k_immob_nh4 * nh4 * enz_mic_n / (kmm_mic_nh4 * denom)

!     denom = 1.0 + nh4 / kmm_nit_nh4 + enz_plant_n / kmm_plant_nh4 + &
!          enz_mic_n / kmm_mic_nh4 + enz_nit_n / kmm_nit_nh4
!     f_nit_nh4 = k_nit_nh4 * nh4 * enz_nit_n / (kmm_nit_nh4 * denom)

! !    denom = 1. + nh4 / kmm_plant_nh4 + no3 / kmm_plant_no3 + enz_plant_n / &
! !         kmm_plant_no3 + enz_mic_n / kmm_mic_no3 + enz_den_n / kmm_den_no3
! !    f_plant_no3 = k_plant_no3 * no3 * enz_plant_n / (kmm_plant_no3 * denom)
!     f_plant_no3 = 0.

! !    denom = 1.0 + nh4 / kmm_mic_nh4 + no3 / kmm_mic_no3 + enz_plant_n /  &
! !         kmm_plant_no3 + enz_mic_n / kmm_mic_no3 + enz_den_n / kmm_den_no3
! !    f_immob_no3 = k_immob_no3 * no3 * enz_mic_n / (kmm_mic_no3 * denom)
!     f_immob_no3 = 0.

! !    denom = 1.0 + no3 / kmm_den_no3 + enz_plant_n / kmm_plant_no3 + &
! !         enz_mic_n / kmm_mic_no3 + enz_den_n / kmm_den_no3
! !    f_den_no3 = k_den_no3 * no3 * enz_den_n / (kmm_den_no3 * denom)
!     f_den_no3 = 0.

! !    denom = 1. + pox / kmm_plant_p + enz_plant_p / kmm_plant_p + enz_mic_p / &
! !         kmm_mic_p + enz_surf_p / kmm_surf_p
! !    f_plant_p = k_plant_p * pox * enz_plant_p / (kmm_plant_p * denom)
!     f_plant_p = 0.

! !    denom = 1. + pox / kmm_mic_p + enz_plant_p / kmm_plant_p + enz_mic_p / &
! !         kmm_mic_p + enz_surf_p / kmm_surf_p
! !    f_immob_p = k_mic_p * pox * enz_mic_p / (kmm_mic_p * denom)
!     f_immob_p = 0.

! !    denom = 1. + pox / kmm_surf_p + enz_plant_p / kmm_plant_p + enz_mic_p / &
! !         kmm_mic_p + enz_surf_p / kmm_surf_p
! !    f_surf_p = k_surf_p * pox * enz_mic_p / (kmm_surf_p * denom)
!     f_surf_p = 0.

!     return
!   end subroutine ncom

  subroutine mend_slow_P_layer(consts, plab, pocc, psec, ppar, psol, bulk_den)
    use mend_consts_coms, only: decomp_consts
    implicit none
    type(decomp_consts) :: consts
    real, intent(inout) :: ppar
    real, intent(inout) :: psol
    real, intent(inout) :: plab
    real, intent(inout) :: pocc
    real, intent(inout) :: psec
    real :: weath_flux
    real :: occlu_flux
    real :: adsorp_flux
    real :: desorp_flux
    real, intent(in) :: bulk_den

    ! Subroutine is called daily; rates are in 1/s.
    weath_flux = consts%weath_rate * ppar * 86400.

!!!!    weath_flux = consts%weath_rate / bulk_den
    psol = psol + weath_flux

!! In ed2_sensitivity, not updating ppar here.
    ppar = ppar - weath_flux

    occlu_flux = consts%occlu_rate * psec * 86400.
    psec = psec - occlu_flux
    pocc = pocc + occlu_flux

    adsorp_flux = consts%p_adsorp_rate * plab * 86400.
    plab = plab - adsorp_flux
    psec = psec + adsorp_flux

    desorp_flux = consts%p_desorp_rate * psec * 86400.
    psec = psec - desorp_flux
    plab = plab + desorp_flux

    return
  end subroutine mend_slow_P_layer

end Module mend_derivs
