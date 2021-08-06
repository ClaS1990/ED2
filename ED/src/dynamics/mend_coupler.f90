Module mend_coupler
  implicit none

Contains

  subroutine mend_update_parameters_coupler(soil_tempk, soil_water,  &
       ntext_soil, pH)
    use mend_consts_coms, only: som_consts, som_consts_base
    use grid_coms, only: nzg
    use soil_coms, only: soil
    use consts_coms, only: wdns, grav
    use mend_diagnose, only: mend_update_parameters
    implicit none

    real, dimension(nzg), intent(in) :: soil_tempk
    real, dimension(nzg), intent(in) :: soil_water
    integer, dimension(nzg), intent(in) :: ntext_soil
    
    real :: tp
    real :: wp
    real, intent(in) :: pH
    real :: wfp
    integer :: ntxt

    ntxt = ntext_soil(nzg)
    tp = soil_tempk(nzg) - 273.15
    wfp = soil_water(nzg) / soil(ntxt)%slmsts
    wp = wdns * grav * soil(ntxt)%slpots / wfp**soil(ntxt)%slbs * 1.0e-6

    call mend_update_parameters(tp, wp, pH, wfp, som_consts, som_consts_base, &
         soil(ntxt)%slmsts, soil(ntxt)%sfldcap)

    return
  end subroutine mend_update_parameters_coupler

  subroutine mend_init(sens_params)
    use ed_state_vars, only: edgrid_g, edtype, polygontype, sitetype
    use nutrient_constants, only: soil_bulk_den, soil_ph
    use mend_consts_coms, only: mend_init_consts, som_consts
    use mend_som, only: mend_som_init
    use mend_state_vars, only: npom, nwood, mend_zero_vars, mend_mm_time
    use nutrient_constants, only: soil_cpct, soil_som_c2n, soil_totp,   &
         soil_extrp
    use soil_coms, only: soil, isoilflg
    implicit none

    type(edtype), pointer :: cgrid
    type(polygontype), pointer :: cpoly
    type(sitetype), pointer :: csite
    integer :: ipy
    integer :: isi
    integer :: ipa
    integer :: iwood
    integer, intent(in) :: sens_params

    write(*,*)'sens_params = ',sens_params
    call mend_init_consts(sens_params)

    mend_mm_time = 0.

    cgrid => edgrid_g(1)
    do ipy = 1, cgrid%npolygons
       cpoly => cgrid%polygon(ipy)
       do isi = 1, cpoly%nsites
          csite => cpoly%site(isi)
          do ipa = 1, csite%npatches

             call mend_zero_vars(csite%mend_mm, ipa, ipa)

             if(isoilflg(1) == 3)then
                csite%mend%bulk_den(ipa) = soil(cpoly%ntext_soil(1,isi))%slden
             else
                csite%mend%bulk_den(ipa) = soil_bulk_den
             endif

             csite%mend_mm%bulk_den(ipa) = csite%mend%bulk_den(ipa)
             csite%mend%pH(ipa) = soil_ph
             csite%mend_mm%pH(ipa) = csite%mend%pH(ipa)

             call mend_som_init(npom,  &
                  csite%mend%som%cvars%pom(:,ipa), csite%mend%som%cvars%dom(ipa), &
                  csite%mend%som%cvars%enz_pom(:,ipa), csite%mend%som%cvars%mom(ipa),  &
                  csite%mend%som%cvars%qom(ipa), csite%mend%som%cvars%enz_mom(ipa),  &
                  csite%mend%som%cvars%amb(ipa), csite%mend%som%cvars%dmb(ipa),   &
                  csite%mend%som%nvars%pom(:,ipa), csite%mend%som%nvars%dom(ipa),  &
                  csite%mend%som%nvars%enz_pom(:,ipa), csite%mend%som%nvars%mom(ipa),  &
                  csite%mend%som%nvars%qom(ipa), csite%mend%som%nvars%enz_mom(ipa),  &
                  csite%mend%som%nvars%amb(ipa), csite%mend%som%nvars%dmb(ipa),   &
                  csite%mend%som%pvars%pom(:,ipa), csite%mend%som%pvars%dom(ipa), &
                  csite%mend%som%pvars%enz_pom(:,ipa), csite%mend%som%pvars%mom(ipa),  &
                  csite%mend%som%pvars%qom(ipa), csite%mend%som%pvars%enz_mom(ipa),  &
                  csite%mend%som%pvars%amb(ipa), csite%mend%som%pvars%dmb(ipa),   &
                  csite%mend%som%fluxes%co2_lost(ipa),  &
                  csite%mend%som%fluxes%nmin(ipa), &
                  csite%mend%som%fluxes%nitr(ipa), &
                  csite%mend%som%invars%nh4(ipa), &
                  csite%mend%som%invars%no3(ipa), csite%mend%som%invars%psol(ipa), &
                  csite%mend%som%invars%plab(ipa), csite%mend%som%fluxes%ngas_lost(ipa), &
                  csite%mend%som%cvars%enz_ptase(ipa), csite%mend%som%nvars%enz_ptase(ipa), &
                  csite%mend%som%pvars%enz_ptase(ipa), csite%mend%som%fluxes%nh4_dep(ipa),  &
                  csite%mend%som%fluxes%no3_dep(ipa), csite%mend%som%fluxes%ppar_dep(ipa), &
                  csite%mend%som%fluxes%nh4_plant(:,ipa), csite%mend%som%fluxes%nh4_bnf(ipa), &
                  csite%mend%som%fluxes%no3_plant(:,ipa), csite%mend%som%fluxes%c_leach(ipa), &
                  csite%mend%som%fluxes%n_leach(ipa), csite%mend%som%fluxes%p_leach(ipa),  &
                  csite%mend%som%fluxes%p_plant(:,ipa), csite%mend%som%invars%pocc(ipa),  &
                  csite%mend%som%invars%psec(ipa), &
                  csite%mend%som%invars%ppar(ipa), csite%mend%som%plvars%enz_plant_n(:,ipa),  &
                  csite%mend%som%plvars%enz_plant_p(:,ipa), csite%mend%som%plvars%vnh4up_plant(:,ipa),  &
                  csite%mend%som%plvars%vno3up_plant(:,ipa),  &
                  csite%mend%som%plvars%vpup_plant(:,ipa),   &
                  csite%mend%som%plvars%plant_input_C_pom(:,ipa), &
                  csite%mend%som%plvars%plant_input_N_pom(:,ipa), &
                  csite%mend%som%plvars%plant_input_P_pom(:,ipa), &
                  csite%mend%som%plvars%plant_input_C_dom(ipa), &
                  csite%mend%som%plvars%plant_input_N_dom(ipa), &
                  csite%mend%som%plvars%plant_input_P_dom(ipa), &
                  som_consts, csite%mend%bulk_den(ipa),   &
                  soil_cpct, soil_som_c2n, soil_totp, soil_extrp, isoilflg(1), &
                  isi, ipy, som_consts%eff_soil_depth)
          enddo
       enddo
    enddo
    
    return
  end subroutine mend_init

  subroutine mend_extern_forcing(mend, ipa, ncohorts, broot, nplant, &
       pft, krdepth, slden, nstorage, pstorage,   &
       nstorage_max, pstorage_max, water_supply_nl, lai,enz_alloc_frac_p)
    use mend_state_vars, only: mend_model, nwood
    use mend_som, only: mend_som_extern_forcing
    use mend_consts_coms, only: som_consts
    use nutrient_constants, only: ndep_rate, pdep_rate, ndep_appl, pdep_appl, nlsl
    use ed_misc_coms, only: current_time
    use mend_plant, only: mend_som_plant_enzymes
    use soil_coms, only: nzg
    implicit none
    type(mend_model) :: mend
    integer, intent(in) :: ipa
    integer :: iwood
    integer, intent(in) :: ncohorts
    integer, intent(in), dimension(ncohorts) :: pft
    integer, intent(in), dimension(ncohorts) :: krdepth
    real, intent(in), dimension(ncohorts) :: broot
    real, intent(in), dimension(ncohorts) :: nplant
    real, intent(in), dimension(ncohorts) :: nstorage
    real, intent(in), dimension(ncohorts) :: pstorage
    real, intent(in), dimension(ncohorts) :: nstorage_max
    real, intent(in), dimension(ncohorts) :: pstorage_max
    real, intent(in), dimension(ncohorts) :: lai
    real, intent(in), dimension(ncohorts) :: water_supply_nl
    real, intent(in), dimension(ncohorts) :: enz_alloc_frac_p
    real :: broot_total
    integer :: ico
    real, intent(in) :: slden

    call mend_som_extern_forcing(ndep_rate, som_consts, slden, &
         mend%som%fluxes%nh4_dep(ipa), mend%som%fluxes%no3_dep(ipa), &
         pdep_rate, mend%som%fluxes%ppar_dep(ipa), current_time%year, &
         ndep_appl, pdep_appl)
    
    call mend_som_plant_enzymes(ncohorts, broot, nplant, pft,  &
         krdepth, slden, mend%som%plvars%enz_plant_n(:,ipa),  &
         mend%som%plvars%enz_plant_p(:,ipa), &
         mend%som%plvars%vnh4up_plant(:,ipa),  &
         mend%som%plvars%vno3up_plant(:,ipa), &
         mend%som%plvars%vpup_plant(:,ipa), som_consts, nstorage, pstorage, &
         nstorage_max, pstorage_max, water_supply_nl, lai, enz_alloc_frac_p)

    return
  end subroutine mend_extern_forcing

  subroutine mend_derivs_coupler(som, d_som, &
       csite, ipa, som_water_drainage, soil_water, d_can_co2, &
       d_co2budget_storage,ccapcani, ntext_soil)
    use mend_exchange, only: mend_plant2som_exchange, zero_exchange_vars, inc_exchange_vars, &
         mend_som2canopy_exchange
    use ed_state_vars, only: sitetype
    use mend_derivs, only: mend_derivs_layer
    use mend_consts_coms, only: som_consts
    use mend_state_vars, only: mend_vars, exchange_vars, npom
    use grid_coms, only: nzg
    use soil_coms, only: dslz, soil
    use consts_coms, only: wdns, pi1
    use nutrient_constants, only: nlsl

    implicit none

    real(kind=8), intent(in) :: som_water_drainage
    real(kind=8), dimension(nzg), intent(in) :: soil_water
    integer, intent(in) :: ipa
    type(sitetype), target :: csite
    type(mend_vars) :: som
    type(mend_vars) :: d_som
    real(kind=8), intent(in) :: ccapcani
    type(exchange_vars) :: plant2som

    real, dimension(npom) :: input_pom_c_net
    real, dimension(npom) :: input_pom_n_net
    real, dimension(npom) :: input_pom_p_net
    real :: input_dom_c_net
    real :: input_dom_n_net
    real :: input_dom_p_net
    real :: input_nh4_net
    real :: input_no3_net
    real :: input_psol_net
    real :: input_ppar_net
    real :: gm2_mgg
    real :: total_water
    integer :: k
    real :: som_water_drainage_ps  ! units: 1/s
    real(kind=8), intent(inout) :: d_can_co2
    real(kind=8), intent(inout) :: d_co2budget_storage
    real :: wfp
    integer, dimension(nzg) :: ntext_soil

    wfp = soil_water(nzg) / soil(ntext_soil(nzg))%slmsts
    gm2_mgg = 1. / (som_consts%eff_soil_depth * csite%mend%bulk_den(ipa))
    total_water = 0.
    do k = nlsl, nzg
       total_water = total_water + soil_water(k) * dslz(k)
    enddo
    ! No negative drainage.
    som_water_drainage_ps = max(0., som_water_drainage) / (total_water * wdns)

    call mend_plant2som_exchange(npom, &
         plant2som%pom_c, plant2som%pom_n, plant2som%pom_p, &
         plant2som%dom_c, plant2som%dom_n, plant2som%dom_p, &
         plant2som%nh4, plant2som%no3, plant2som%psol,   &
         csite%plant_input_C(:,ipa), csite%plant_input_N(:,ipa),  &
         csite%plant_input_P(:,ipa))

    d_som%plvars%plant_input_C_pom(:,ipa) = plant2som%pom_c * gm2_mgg * 1000./86400.
    d_som%plvars%plant_input_N_pom(:,ipa) = plant2som%pom_n * gm2_mgg * 1000./86400.
    d_som%plvars%plant_input_P_pom(:,ipa) = plant2som%pom_p * gm2_mgg * 1000./86400.
    d_som%plvars%plant_input_C_dom = plant2som%dom_c * gm2_mgg * 1000./86400.
    d_som%plvars%plant_input_N_dom = plant2som%dom_n * gm2_mgg * 1000./86400.
    d_som%plvars%plant_input_P_dom = plant2som%dom_p * gm2_mgg * 1000./86400.

    input_pom_c_net = (plant2som%pom_c) * &
         gm2_mgg * 1000./86400.
    input_dom_c_net = (plant2som%dom_c) * &
         gm2_mgg * 1000./86400.
    input_pom_n_net = (plant2som%pom_n) * &
         gm2_mgg * 1000./86400.
    input_dom_n_net = (plant2som%dom_n) * &
         gm2_mgg * 1000./86400.
    input_pom_p_net = (plant2som%pom_p) * &
         gm2_mgg * 1000./86400.
    input_dom_p_net = (plant2som%dom_p) * &
         gm2_mgg * 1000./86400.
    input_nh4_net = (plant2som%nh4) * &
         gm2_mgg  * 1000./86400. + som%fluxes%nh4_dep(ipa)
    input_no3_net = (plant2som%no3) * &
         gm2_mgg  * 1000./86400. + som%fluxes%no3_dep(ipa)
    input_psol_net = (plant2som%psol) * &
         gm2_mgg  * 1000./86400. + som%fluxes%ppar_dep(ipa)
    input_ppar_net = 0.

    call mend_derivs_layer(npom, som_consts,  &
         som%cvars%pom(:,ipa), input_pom_c_net, d_som%cvars%pom(:,ipa),  &
         som%cvars%dom(ipa), input_dom_c_net, d_som%cvars%dom(ipa),  &
         som%cvars%enz_pom(:,ipa), d_som%cvars%enz_pom(:,ipa),  &
         som%cvars%mom(ipa), d_som%cvars%mom(ipa),  &
         som%cvars%qom(ipa), d_som%cvars%qom(ipa),  &
         som%cvars%enz_mom(ipa), d_som%cvars%enz_mom(ipa),  &
         som%cvars%amb(ipa), d_som%cvars%amb(ipa),   &
         som%cvars%dmb(ipa), d_som%cvars%dmb(ipa),   &
         som%fluxes%co2_lost(ipa), d_som%fluxes%co2_lost(ipa),  &
         som%fluxes%nmin(ipa), d_som%fluxes%nmin(ipa),  &
         som%fluxes%nitr(ipa), d_som%fluxes%nitr(ipa),  &
         som%nvars%pom(:,ipa), input_pom_n_net, d_som%nvars%pom(:,ipa),  &
         som%nvars%dom(ipa), input_dom_n_net, d_som%nvars%dom(ipa),  &
         som%nvars%enz_pom(:,ipa), d_som%nvars%enz_pom(:,ipa),  &
         som%nvars%mom(ipa), d_som%nvars%mom(ipa),   &
         som%nvars%qom(ipa), d_som%nvars%qom(ipa),   &
         som%nvars%enz_mom(ipa), d_som%nvars%enz_mom(ipa),   &
         som%nvars%amb(ipa), d_som%nvars%amb(ipa),  &
         som%nvars%dmb(ipa), d_som%nvars%dmb(ipa),  &
         som%invars%nh4(ipa), input_nh4_net, d_som%invars%nh4(ipa), &
         som%invars%no3(ipa), input_no3_net, d_som%invars%no3(ipa),  &
         som%pvars%pom(:,ipa), input_pom_p_net, d_som%pvars%pom(:,ipa),  &
         som%pvars%dom(ipa), input_dom_p_net, d_som%pvars%dom(ipa),  &
         som%pvars%enz_pom(:,ipa), d_som%pvars%enz_pom(:,ipa),  &
         som%pvars%mom(ipa), d_som%pvars%mom(ipa),   &
         som%pvars%qom(ipa), d_som%pvars%qom(ipa),  &
         som%pvars%enz_mom(ipa), d_som%pvars%enz_mom(ipa),   &
         som%pvars%amb(ipa), d_som%pvars%amb(ipa),  &
         som%pvars%dmb(ipa), d_som%pvars%dmb(ipa),   &
         som%invars%psol(ipa), input_psol_net, d_som%invars%psol(ipa),  &
         som%invars%plab(ipa), d_som%invars%plab(ipa), &
         som%fluxes%ngas_lost(ipa), d_som%fluxes%ngas_lost(ipa), &
         som%cvars%enz_ptase(ipa), d_som%cvars%enz_ptase(ipa), &
         som%nvars%enz_ptase(ipa), d_som%nvars%enz_ptase(ipa), &
         som%pvars%enz_ptase(ipa), d_som%pvars%enz_ptase(ipa),  &
         d_som%fluxes%nh4_plant(:,ipa), &
         som%fluxes%nh4_bnf(ipa), d_som%fluxes%nh4_bnf(ipa), &
         d_som%fluxes%no3_plant(:,ipa), &
         som%fluxes%c_leach(ipa), d_som%fluxes%c_leach(ipa), &
         som%fluxes%n_leach(ipa), d_som%fluxes%n_leach(ipa),  &
         som%fluxes%p_leach(ipa), d_som%fluxes%p_leach(ipa), &
         d_som%fluxes%p_plant(:,ipa), &
         som%invars%pocc(ipa), d_som%invars%pocc(ipa), &
         som%invars%ppar(ipa), d_som%invars%ppar(ipa), input_ppar_net,  &
         som%plvars%enz_plant_n(:,ipa), som%plvars%enz_plant_p(:,ipa), &
         som%plvars%vnh4up_plant(:,ipa), som%plvars%vno3up_plant(:,ipa),  &
         som%plvars%vpup_plant(:,ipa), som_water_drainage_ps, &
         csite%mend%bulk_den(ipa), pi1, wfp)

    call mend_som2canopy_exchange(d_som%fluxes%co2_lost(ipa),  &
         csite%mend%bulk_den(ipa), som_consts, &
         d_can_co2, d_co2budget_storage, ccapcani)

    return
  end subroutine mend_derivs_coupler

  subroutine mend_update_diag(mend,ipa)
    use mend_state_vars, only: mend_model, nwood
    use mend_consts_coms, only: som_consts
    use mend_diagnose, only: mend_update_diag_layer
    implicit none
    type(mend_model) :: mend
    integer :: iwood
    integer, intent(in) :: ipa

    call mend_update_diag_layer(mend%som, som_consts, ipa, mend%bulk_den(ipa))

    return
  end subroutine mend_update_diag

  subroutine mend_slow_P(mend, ipa, ifert)
    use mend_state_vars, only: mend_model, nwood
    use mend_derivs, only: mend_slow_P_layer
    use mend_consts_coms, only: som_consts
    use mend_diagnose, only: mend_update_diag_layer
    use nutrient_constants, only: ndep_appl, pdep_appl
    implicit none
    type(mend_model) :: mend
    integer, intent(in) :: ipa
    integer :: iwood
    integer, intent(in) :: ifert

    call mend_slow_P_layer(som_consts, mend%som%invars%plab(ipa), &
         mend%som%invars%pocc(ipa), mend%som%invars%psec(ipa), &
         mend%som%invars%ppar(ipa), &
         mend%som%invars%psol(ipa), mend%bulk_den(ipa))
    call mend_update_diag_layer(mend%som, som_consts, ipa, mend%bulk_den(ipa))

    ! units were mg/m2/s
    if(ifert == 1)then
       mend%som%invars%nh4(ipa) = mend%som%invars%nh4(ipa) + ndep_appl*86400*365.25/3. / &
            (1000. * mend%bulk_den(ipa) * som_consts%eff_soil_depth)
       mend%som%invars%psol(ipa) = mend%som%invars%psol(ipa) + pdep_appl*86400*365.25/3. / &
            (1000. * mend%bulk_den(ipa) * som_consts%eff_soil_depth)
    endif

    return
  end subroutine mend_slow_P

end Module mend_coupler
