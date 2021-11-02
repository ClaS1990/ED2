Module mend_plant
  implicit none

  real, parameter :: water_supply_scale = 0.01

Contains

  subroutine mend_som_plant_enzymes(ncohorts, broot, nplant, pft, &
       krdepth, slden, enz_plant_n, &
       enz_plant_p, vnh4up_plant, vno3up_plant, vpup_plant, consts, &
       nstorage, pstorage, nstorage_max, pstorage_max, water_supply_nl, lai, &
       enz_alloc_frac_p)
    use mend_consts_coms, only: decomp_consts
    use pft_coms, only: root_beta
    use soil_coms, only: slz
    use nutrient_constants, only: nlsl
    use ed_max_dims, only: n_pft
    implicit none

    integer, intent(in) :: ncohorts
    integer, dimension(ncohorts), intent(in) :: pft
    integer, dimension(ncohorts), intent(in) :: krdepth
    real, dimension(ncohorts), intent(in) :: broot
    real, dimension(ncohorts), intent(in) :: nplant
    real, dimension(ncohorts), intent(in) :: water_supply_nl
    real, dimension(ncohorts), intent(in) :: lai
    real, intent(in) :: slden
    real, intent(inout), dimension(n_pft) :: enz_plant_n
    real, intent(inout), dimension(n_pft) :: enz_plant_p
    real, intent(inout), dimension(n_pft) :: vnh4up_plant
    real, intent(inout), dimension(n_pft) :: vno3up_plant
    real, intent(inout), dimension(n_pft) :: vpup_plant
    real, intent(in), dimension(ncohorts) :: nstorage
    real, intent(in), dimension(ncohorts) :: pstorage
    real, intent(in), dimension(ncohorts) :: nstorage_max
    real, intent(in), dimension(ncohorts) :: pstorage_max
    real, intent(in), dimension(ncohorts) :: enz_alloc_frac_p
    type(decomp_consts) :: consts
    integer :: ico
    real :: broot_nl
    real :: n_limit_factor
    real :: p_limit_factor
    real :: nstorage_total
    real :: nstorage_max_total
    real :: pstorage_total
    real :: pstorage_max_total
    real :: nplant_total
    integer, parameter :: nutr_limit_scheme = 1

    ! Plant enzyme carrier abundance.

    ! (1) In a typical tropical forest, Zhu et al. assume that the number of plant
    ! carrier enzymes is equal to the number of microbial carrier
    ! enzymes. But this statement is ambiguous. Does it apply to NH4
    ! carriers, PO4 carriers, sum of all carriers? Does it mean that the number of 
    ! enzymes is the same (which they do not simulate), or the amount of C in enzymes 
    ! (which they do simulate) is the same? 
    ! Here I assume that it applies to the number of enzymes, in moles. I further assume
    ! that it applies to each enzyme type individually.

    ! (2) Zhu et al. also say that the scaling factor between
    ! microbial biomass and enzyme abundance is 0.05. Does this mean 0.05 g C enzyme
    ! per g C microbe? Does it apply to all types of transport enzymes - is it a sum
    ! over types? Are transport enzymes generic? Further, Tang and Riley suggest a huge
    ! possible range for this number. Should transport enzymes rather scale with surface
    ! area? I will assume a generic transport enzyme. 
    ! "Abundance" is 0.05 amb_c [mg enz C / g soil].
    ! They assume amb_c = 0.1 gC/m2 or, say, 7e-4 mgC/gSoil.
    ! Enzyme amount is then 4e-5 mgC/gSoil.
    ! They estimate plant root biomass at 400 gC/m2, or 3 mgC/gSoil.
    ! Then, enzyme C per fine root C should be: 1e-5. This agrees with their 0.0000125.

    ! ANOTHER APPROACH.
    ! Say that it applies to number counts (moles).
    ! Number of enzymes is 0.05 * amb_c * (g soil) * (# microbes/ mg micr biomass C)
    ! Plants have same number.
    ! 0.05 * amb_c * (g soil) * (# microbes / mg micr biomass C) = (# plant enz / mg C fine root) (mg C fine root)
    ! 0.05 * amb_c * (g soil) * (# microbes / mg micr biomass C) = (# plant enz / mg C fine root) (3 mgC/gSoil) (g soil)
    ! (# plant enz / mg C fine root) = 1e-5 * (# microbes / mg micr biomass C) = K1
    ! This, I guess, would represent the sum of all transporter enzymes.
    ! 5e7 cells/gram soil, 5e10 cell/kgSoil, 4e13 cell/m3 soil, 8e12 cell/m2, 8e13 cell/g micr C, 8e10 cell/mg micr C

    ! (3) The MM constants are listed in grams/m2. Presumably this is in gN/m2, but I am
    ! not 100% sure. If so, I convert to gN/kgSoil. So conversion of the plant enzymes is:
    ! [E_plant_sum] = (8e5) * (1e6 mg C fine root / kg C fine root) * broot * nplant / soildepth / bulkden / AvoNum * (14 gN/mol)  
    
    ! (4) Going back to the first alternative,
    ! (g plant enz C /g plant fine root C) = 1e-5
    ! (mol plant enz C /g plant fine root C) = 1e-5 / 12 [g C / mol C]
    ! (mol plant enz N /g plant fine root C) = 1e-5 / 12 / C:N_plant_enz
    ! (g plant enz N /g plant fine root C) = 1e-5 / 12 / C:N_plant_enz * 14
    ! (g plant enz N /kg plant fine root C) = 1e-5 / 12 / C:N_plant_enz * 14 * 1000
    ! (g plant enz N /plant) = 1e-5 / 12 / C:N_plant_enz * 14 * 1000 * broot
    ! (g plant enz N /m2) = 1e-5 / 12 / C:N_plant_enz * 14 * 1000 * broot * nplant
    ! (g plant enz N /m3) = 1e-5 / 12 / C:N_plant_enz * 14 * 1000 * broot * nplant / soildepth
    ! (g plant enz N /kg soil) = 1e-5 / 12 / C:N_plant_enz * 14 * 1000 * broot * nplant / soildepth / bulkden
    ! Still, this doesn't make sense. We need the number of moles of transporter enzyme, and then
    ! multiply that by the MM. 
    ! moles enzyme = (# enzyme) / AvoNum = (total mg C in enzyme) / (mg C / enzyme) / AvoNum
    ! moles enzyme = (0.05 * amb_c) * (g soil) / (mgC/enzyme) / AvoNum
    ! (moles enzyme / gSoil) = (0.05 * amb_c) / (mgC/enzyme) / AvoNum
    ! (gN / gSoil) = (0.05 * amb_c) / (mgC/enzyme) / AvoNum * 14
    ! (gN / kgSoil) = (0.05 * amb_c) / (mgC/enzyme) / AvoNum * 14 * 1000
    ! (gN / m3) = (0.05 * amb_c) / (mgC/enzyme) / AvoNum * 14 * 1000 * 750
    ! (gN / m2) = (0.05 * amb_c) / (mgC/enzyme) / AvoNum * 14 * 1000 * 750 * 0.2
    ! amb_c = 0.1 gC/m2 / (0.2 m) / (750 kg/m3) = 7e-4 mgC/gSoil
    ! (gN / kgSoil) = 1e-22 / (mgC/enzyme)

    ! Scaling factor (gCenz / gCrootbiomass) * 
    ! fine root biomass (kgCfineroot/m2) *  ! kgCenz / m2
    ! 1000gCfineroot / kgCfineroot  /       ! gCenz / m2
    ! characteristic soil depth (m) /       ! gCenz / m3
    ! soil bulk density (g/m3) *            ! gCenz / gsoil
    ! 1000mgCenz / 1 gCenz     /            ! mgCenz / gsoil
    ! enzC:Nratio                           ! mgNenz / gsoil

    do ico = 1, ncohorts
       broot_nl = broot(ico) * (1. - root_beta(pft(ico))**  &
            (-slz(nlsl)/(-slz(krdepth(ico)))))

       if(nstorage(ico) < 0.001 * nstorage_max(ico))then
          n_limit_factor = 1.
       else
          n_limit_factor = 1.0 - (nstorage(ico) / &
               nstorage_max(ico))**consts%wexp
       endif
       n_limit_factor = max(min(1., n_limit_factor), 0.)
       p_limit_factor = 1.0 - (pstorage(ico) / &
            pstorage_max(ico))**consts%wexp
       p_limit_factor = max(min(1., p_limit_factor), 0.)
          
       enz_plant_n(pft(ico)) = enz_plant_n(pft(ico)) + &
            consts%enz2biomass_plant *   &
            nplant(ico) * broot_nl /   &
            consts%eff_soil_depth / slden * 14.
       enz_plant_p(pft(ico)) = enz_plant_p(pft(ico)) + &
            consts%enz2biomass_plant * enz_alloc_frac_p(ico) *   &
            nplant(ico) * broot_nl /   & 
            consts%eff_soil_depth / slden * 31.
       
       vnh4up_plant(pft(ico)) = vnh4up_plant(pft(ico)) + consts%vnh4up_plant_base *  &
            consts%enz2biomass_plant *   &
            nplant(ico) * broot_nl /   & 
            consts%eff_soil_depth / slden * 14. * &
            n_limit_factor
       vno3up_plant(pft(ico)) = vno3up_plant(pft(ico)) + consts%vno3up_plant_base *  &
            consts%enz2biomass_plant *   &
            nplant(ico) * broot_nl /   & 
            consts%eff_soil_depth / slden * 14. * &
            n_limit_factor
       vpup_plant(pft(ico)) = vpup_plant(pft(ico)) + consts%vpup_plant_base *  &
            consts%enz2biomass_plant * enz_alloc_frac_p(ico) *  &
            nplant(ico) * broot_nl /   & 
            consts%eff_soil_depth / slden * 31. * &
            p_limit_factor
    enddo

    return
  end subroutine mend_som_plant_enzymes

  subroutine mend_som_plant_feedback(nh4_plant, no3_plant, p_plant, slden,  &
       consts, ncohorts, nstorage, pstorage, nstorage_max, pstorage_max, &
       nplant, broot, rh, co2_lost, pft, krdepth, water_supply_nl, lai, &
       enz_alloc_frac_p)
    use ed_misc_coms, only: dtlsm
    use mend_consts_coms, only: decomp_consts
    use nutrient_constants, only: nlsl
    use ed_max_dims, only: n_pft
    use pft_coms, only: root_beta
    use soil_coms, only: slz, nzg
    implicit none

    real :: broot_nl
    integer, dimension(ncohorts), intent(in) :: pft
    integer, dimension(ncohorts), intent(in) :: krdepth
    integer, intent(in) :: ncohorts
    real, intent(inout), dimension(n_pft) :: nh4_plant
    real, intent(inout), dimension(n_pft) :: no3_plant
    real, intent(inout), dimension(n_pft) :: p_plant
    real, intent(in) :: slden
    type(decomp_consts) :: consts
    real, dimension(n_pft) :: plant_n_uptake
    real, dimension(n_pft) :: plant_p_uptake
    real, dimension(n_pft) :: total_n_activity
    real, dimension(n_pft) :: total_n_activity_nolimit
    real, dimension(n_pft) :: total_p_activity
    real, dimension(n_pft) :: total_p_activity_nolimit
    integer :: ico
    real, dimension(ncohorts) :: n_limit_factor
    real, dimension(ncohorts) :: p_limit_factor
    real, intent(inout), dimension(ncohorts) :: nstorage
    real, intent(inout), dimension(ncohorts) :: pstorage
    real, intent(in), dimension(ncohorts) :: lai
    real, intent(in), dimension(ncohorts) :: water_supply_nl
    real, intent(in), dimension(ncohorts) :: nstorage_max
    real, intent(in), dimension(ncohorts) :: pstorage_max
    real, intent(in), dimension(ncohorts) :: nplant
    real, intent(in), dimension(ncohorts) :: broot
    real, intent(in), dimension(ncohorts) :: enz_alloc_frac_p
    real, intent(inout) :: rh
    real :: co2_lost_units
    real, intent(in) :: co2_lost
    real :: nstorage_total
    real :: nstorage_max_total
    real :: pstorage_total
    real :: pstorage_max_total
    real :: nplant_total
    real :: total_nlim
    real :: total_plim
    real :: transp_fact

    ! kgN/m2
    plant_n_uptake = (nh4_plant + no3_plant) * 1.0e-6  * &
         1000. * slden * consts%eff_soil_depth  
    ! kgP/m2
    plant_p_uptake = p_plant * 1.0e-6  * &
         1000. * slden * consts%eff_soil_depth  

    total_n_activity = 0.
    total_p_activity = 0.
    do ico = 1, ncohorts
       broot_nl = broot(ico) * (1. - root_beta(pft(ico))**  &
            (-slz(nlsl)/(-slz(krdepth(ico)))))
       

       if(nstorage(ico) < 0.001 * nstorage_max(ico))then
          n_limit_factor(ico) = 1.
       else
          n_limit_factor(ico) = 1.0 - (nstorage(ico) / &
               nstorage_max(ico))**consts%wexp
       endif
       n_limit_factor(ico) = max(min(1., n_limit_factor(ico)), 0.)
       p_limit_factor(ico) = 1.0 - (pstorage(ico) / &
            pstorage_max(ico))**consts%wexp
       p_limit_factor(ico) = max(min(1., p_limit_factor(ico)), 0.)

       total_n_activity(pft(ico)) = total_n_activity(pft(ico)) + nplant(ico) *   &
            broot_nl * n_limit_factor(ico)
       total_p_activity(pft(ico)) = total_p_activity(pft(ico)) + nplant(ico) *   &
            broot_nl * p_limit_factor(ico) * enz_alloc_frac_p(ico)
    enddo

    do ico = 1, ncohorts
       broot_nl = broot(ico) * (1. - root_beta(pft(ico))**  &
            (-slz(nlsl)/(-slz(krdepth(ico)))))

       if(total_n_activity(pft(ico)) > 1.e-30)then
          nstorage(ico) = nstorage(ico) + plant_n_uptake(pft(ico)) * broot_nl /  &
               total_n_activity(pft(ico)) * n_limit_factor(ico)
       endif

       if(total_p_activity(pft(ico)) > 1.e-30)then
          pstorage(ico) = pstorage(ico) + plant_p_uptake(pft(ico)) * broot_nl /  &
               total_p_activity(pft(ico)) * p_limit_factor(ico) * enz_alloc_frac_p(ico)
       endif
    enddo

    ! gC/kgSoil
!    co2_lost_units = co2_lost
    ! gC/m3Soil
!    co2_lost_units = co2_lost_units * slden
    ! gC/m2Soil
!    co2_lost_units = co2_lost_units * consts%eff_soil_depth
    ! molC/m2Soil
!    co2_lost_units = co2_lost_units / 12.
    ! umolC/m2Soil
!    co2_lost_units = co2_lost_units * 1.0e6
    ! umolC/m2Soil/s
!    co2_lost_units = co2_lost_units / dtlsm

    ! Averaging this over the day.
!    rh = rh + co2_lost_units * dtlsm / 86400.

    return
  end subroutine mend_som_plant_feedback

end Module mend_plant
