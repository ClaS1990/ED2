Module mend_som
  implicit none

Contains

  subroutine mend_som_init(npom,  &
       pom_c, dom_c, enz_pom_c, mom_c, qom_c, enz_mom_c, amb_c, dmb_c,  &
       pom_n, dom_n, enz_pom_n, mom_n, qom_n, enz_mom_n, amb_n, dmb_n,  &
       pom_p, dom_p, enz_pom_p, mom_p, qom_p, enz_mom_p, amb_p, dmb_p,  &
       co2_lost, nmin, nitr, nh4, no3, psol, plab, ngas_lost, enz_ptase_c, enz_ptase_n, &
       enz_ptase_p, nh4_dep, no3_dep, ppar_dep, nh4_plant, nh4_bnf, &
       no3_plant, c_leach, n_leach, p_leach, p_plant, pocc, psec, ppar,  &
       enz_plant_n, enz_plant_p, vnh4up_plant, vno3up_plant, &
       vpup_plant, &
       plant_input_C_pom, plant_input_N_pom, plant_input_P_pom, &
       plant_input_C_dom, plant_input_N_dom, plant_input_P_dom, &
       consts, bulk_den, soil_cpct, soil_som_c2n, soil_totp, soil_extrp, isoilflg, &
       isi, ipy, eff_soil_depth)
    use ed_max_dims, only: n_pft
    use mend_consts_coms, only: decomp_consts
    use ed_work_vars, only: work_v
    implicit none

    integer, intent(in) :: isoilflg
    real, intent(in) :: soil_cpct
    real, intent(in) :: soil_som_c2n
    real, intent(in) :: soil_totp
    real, intent(in) :: soil_extrp
    type(decomp_consts) :: consts
    integer, intent(in) :: npom
    real, intent(out) :: dom_c
    real, intent(out) :: amb_c
    real, intent(out) :: dmb_c
    real, intent(out) :: mom_c
    real, intent(out) :: qom_c
    real, dimension(npom), intent(out) :: pom_c
    real, dimension(npom), intent(out) :: enz_pom_c
    real, intent(out) :: enz_mom_c
    real, intent(out) :: enz_ptase_c

    real, intent(out) :: dom_n
    real, intent(out) :: amb_n
    real, intent(out) :: dmb_n
    real, intent(out) :: mom_n
    real, intent(out) :: qom_n
    real, dimension(npom), intent(out) :: pom_n
    real, dimension(npom), intent(out) :: enz_pom_n
    real, intent(out) :: enz_mom_n
    real, intent(out) :: enz_ptase_n

    real, intent(out) :: dom_p
    real, intent(out) :: amb_p
    real, intent(out) :: dmb_p
    real, intent(out) :: mom_p
    real, intent(out) :: qom_p
    real, dimension(npom), intent(out) :: pom_p
    real, dimension(npom), intent(out) :: enz_pom_p
    real, intent(out) :: enz_mom_p
    real, intent(out) :: enz_ptase_p
    real, intent(out) :: no3_dep
    real, intent(out) :: nh4_dep

    real, intent(out) :: co2_lost
    real, intent(out) :: ngas_lost
    real, intent(out) :: nh4
    real, intent(out) :: no3
    real, intent(out) :: psol
    real, intent(out) :: plab
    real, intent(out) :: pocc
    real, intent(out) :: psec
    real, intent(out) :: ppar
    real, intent(out) :: ppar_dep

    real, intent(out), dimension(n_pft) :: nh4_plant
    real, intent(out) :: nh4_bnf
    real, intent(out), dimension(n_pft) :: no3_plant
    real, intent(out) :: nmin
    real, intent(out) :: nitr
    real, intent(out) :: c_leach
    real, intent(out) :: n_leach
    real, intent(out) :: p_leach
    real, intent(out), dimension(n_pft) :: p_plant
    real, intent(out), dimension(n_pft) :: enz_plant_n
    real, intent(out), dimension(n_pft) :: enz_plant_p
    real, intent(out), dimension(n_pft) :: vnh4up_plant
    real, intent(out), dimension(n_pft) :: vno3up_plant
    real, intent(out), dimension(n_pft) :: vpup_plant
    real, intent(out), dimension(npom) :: plant_input_C_pom
    real, intent(out), dimension(npom) :: plant_input_N_pom
    real, intent(out), dimension(npom) :: plant_input_P_pom
    real :: my_total, my_b, my_c
    real, intent(out) :: plant_input_C_dom
    real, intent(out) :: plant_input_N_dom
    real, intent(out) :: plant_input_P_dom
    real, intent(in) :: bulk_den
    real :: cpct
    real :: c2n
    real :: totp
    real :: extrp
    real :: c2p
    integer, intent(in) :: isi, ipy
    real :: Yang_conv_fac ! convert gP/m2 to mgP/gSoil
    real, intent(in) :: eff_soil_depth

    Yang_conv_fac = 1. / (eff_soil_depth * bulk_den)

    if(isoilflg == 3)then
       cpct = work_v(1)%orgc(isi,ipy) * 10.0 ! Divide by 100 to turn percentage into fraction.
       ! Then, multiply by 1000 because we want mg/g. Net is a factor of 10.
       c2n = work_v(1)%c2n(isi,ipy)
       extrp = work_v(1)%lab(isi,ipy) * Yang_conv_fac
       ! Tipping et al. (2016) Biogeochemistry
       !       c2p = (cpct*0.1) / (0.012*(cpct*0.1)**0.57)
       c2p = cpct / (work_v(1)%org(isi,ipy) * Yang_conv_fac)
    else
       cpct = soil_cpct
       c2n = soil_som_c2n
       totp = soil_totp
       extrp = soil_extrp
       ! Single number from Waring et al
!       c2p = cpct*0.1/(150.e-4)
       ! Use Tipping et al meta-analysis
       c2p = (cpct*0.1) / (0.012*(cpct*0.1)**0.57)
    endif

    pom_c(1) = 7./31. * cpct
    pom_n(1) = pom_c(1) / c2n
    pom_p(1) = pom_c(1) / c2p

    pom_c(2) = 2.5/31. * cpct
    pom_n(2) = pom_c(2) / c2n
    pom_p(2) = pom_c(2) / c2p

    dom_c = 0.15/31. * cpct
    dom_n = dom_c / c2n
    dom_p = dom_c / c2p

    enz_pom_c(1) = 0.014
    enz_pom_n(1) = enz_pom_c(1) / consts%enz_pom_c2n(1)
    enz_pom_p(1) = enz_pom_c(1) / consts%enz_pom_c2p(1)

    enz_pom_c(2) = 0.00044
    enz_pom_n(2) = enz_pom_c(2) / consts%enz_pom_c2n(2)
    enz_pom_p(2) = enz_pom_c(2) / consts%enz_pom_c2p(2)

    mom_c = 21./31. * cpct
    mom_n = mom_c / c2n
    mom_p = mom_c / c2p

    qom_c = 0.8 / 31. * cpct
    qom_n = qom_c / c2n
    qom_p = qom_c / c2p

    enz_mom_c = 0.013
    enz_mom_n = enz_mom_c / consts%enz_mom_c2n
    enz_mom_p = enz_mom_c / consts%enz_mom_c2p

    amb_c = 0.023
    amb_n = 0.0046
    amb_p = 0.00055

    dmb_c = 0.27
    dmb_n = 0.046
    dmb_p = 0.0055

    enz_ptase_c = 0.0075
    enz_ptase_n = enz_ptase_c / consts%enz_ptase_c2n
    enz_ptase_p = enz_ptase_c / consts%enz_ptase_c2p

    co2_lost = 0.
    ngas_lost = 0.
    nh4 = 1.0e-5
    no3 = 1.0e-5

    if(isoilflg /= 3)then
       ! See Yang et al. 2013, Biogeosciences, 10, 2525-2537
       !  the total is from measurements, but I am very rougly assuming
       !  a 3:1 ratio of occluded:secondary P, based on Yang et al. (2013)
       !  results for inceptisols.
       pocc = (totp - extrp) * 0.001 * 0.75
       psec = (totp - extrp) * 0.001 * 0.25
       ppar = 100. * 2. / bulk_den ! gP/m2 * 1000mg/g / 0.5 m / bulk_den / 1000kg/g
       my_total = extrp * 0.001
    else
       pocc = work_v(1)%occ(isi,ipy) * Yang_conv_fac
       psec = work_v(1)%sec(isi,ipy) * Yang_conv_fac
       ppar = work_v(1)%apa(isi,ipy) * Yang_conv_fac
       my_total = extrp
    endif

    my_b = consts%kmm_plang / bulk_den + consts%plab_max / bulk_den - my_total
    my_c = - consts%kmm_plang / bulk_den * my_total
    psol = (-my_b + sqrt(my_b**2-4.*my_c))*0.5
    plab = my_total - psol
    
    if(plab > consts%plab_max / bulk_den)then
       print*,'initialization problem.'
       print*,'initial plab > plab_max'
       stop
    endif

    nh4_plant = 0.
    nh4_bnf = 0.
    no3_plant = 0.
    nmin = 0.
    nitr = 0.
    c_leach = 0.
    n_leach = 0.
    p_leach = 0.
    p_plant = 0.

    nh4_dep = 0.
    no3_dep = 0.
    ppar_dep = 0.

    enz_plant_n = 0.
    enz_plant_p = 0.

    vnh4up_plant = 0.
    vno3up_plant = 0.
    vpup_plant = 0.

    plant_input_C_pom = 0.
    plant_input_C_dom = 0.
    plant_input_N_pom = 0.
    plant_input_N_dom = 0.
    plant_input_P_pom = 0.
    plant_input_P_dom = 0.

    return
  end subroutine mend_som_init

  subroutine mend_som_extern_forcing(ndep_rate, consts, slden, input_nh4, &
       input_no3, pdep_rate, input_ppar, year, ndep_appl, pdep_appl, ndep_met, pdep_met, &
       nh4dep_met, no3dep_met)

    use mend_consts_coms, only: decomp_consts
    use soil_coms, only: isoilflg
    use met_driver_coms, only: dep_scheme, has_nh4dep, has_no3dep

    implicit none

    type(decomp_consts) :: consts
    real, intent(in) :: slden
    real, intent(in) :: ndep_rate
    real, intent(out) :: input_nh4
    real, intent(out) :: input_no3
    real, intent(in) :: pdep_rate
    real, intent(out) :: input_ppar
    real, intent(in) :: ndep_appl
    real, intent(in) :: pdep_appl
    integer, intent(in) :: year
    real, intent(in) :: ndep_met, pdep_met, no3dep_met, nh4dep_met

    if(dep_scheme == 0)then
       if(has_nh4dep)then
          input_nh4 = nh4dep_met / (consts%eff_soil_depth * slden * 1000.)!*1e10
       else
          input_nh4 = ndep_rate / (consts%eff_soil_depth * slden * 1000.)
       endif
       if(has_no3dep)then
          input_no3 = no3dep_met / (consts%eff_soil_depth * slden * 1000.)
       else
          input_no3 = ndep_rate / (consts%eff_soil_depth * slden * 1000.)
       endif
       input_ppar = pdep_rate / (consts%eff_soil_depth * slden * 1000.)
    elseif(dep_scheme == 1)then
       input_nh4 = ndep_met  / (consts%eff_soil_depth * slden * 1000.)
       input_ppar = pdep_met  / (consts%eff_soil_depth * slden * 1000.)
       input_no3 = 0.
    endif

!    if(year >= 2015)then
!       input_nh4 = input_nh4 + ndep_appl / (consts%eff_soil_depth * slden * 1000.)
!       input_ppar = input_ppar + pdep_appl / (consts%eff_soil_depth * slden * 1000.)
!    endif

    return
  end subroutine mend_som_extern_forcing

end Module mend_som
