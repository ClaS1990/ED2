Module mend_state_vars
  implicit none

  integer, parameter :: npom=2
  integer, parameter :: nwood=0
  real :: mend_mm_time, mend_dm_time

  ! Can be used for SOM or litter
  type org_vars
     real, allocatable, dimension(:,:) :: pom
     real, allocatable, dimension(:,:) :: enz_pom
     real, allocatable, dimension(:) :: dom
     real, allocatable, dimension(:) :: mom
     real, allocatable, dimension(:) :: qom
     real, allocatable, dimension(:) :: enz_mom
     real, allocatable, dimension(:) :: amb
     real, allocatable, dimension(:) :: dmb
     real, allocatable, dimension(:) :: enz_ptase
  end type org_vars

  type inorg_vars
     real, allocatable, dimension(:) :: nh4
     real, allocatable, dimension(:) :: no3
     real, allocatable, dimension(:) :: plab
     real, allocatable, dimension(:) :: psol
     real, allocatable, dimension(:) :: pocc
     real, allocatable, dimension(:) :: psec
     real, allocatable, dimension(:) :: ppar
  end type inorg_vars

  ! External fluxes (atmosphere, plants, runoff)
  type ext_fluxes
     real, allocatable, dimension(:) :: co2_lost
     real, allocatable, dimension(:) :: ngas_lost
     real, allocatable, dimension(:) :: nh4_dep
     real, allocatable, dimension(:) :: no3_dep
     real, allocatable, dimension(:) :: ppar_dep
     real, allocatable, dimension(:) :: c_leach
     real, allocatable, dimension(:) :: n_leach
     real, allocatable, dimension(:) :: p_leach
     real, allocatable, dimension(:) :: nh4_bnf
     real, allocatable, dimension(:,:) :: nh4_plant
     real, allocatable, dimension(:,:) :: no3_plant
     real, allocatable, dimension(:,:) :: p_plant
     real, allocatable, dimension(:) :: nh4_plant_sum
     real, allocatable, dimension(:) :: no3_plant_sum
     real, allocatable, dimension(:) :: p_plant_sum
     real, allocatable, dimension(:) :: nmin
     real, allocatable, dimension(:) :: nitr
  end type ext_fluxes

  type plant_vars
     real, allocatable, dimension(:,:) :: enz_plant_n
     real, allocatable, dimension(:,:) :: enz_plant_p
     real, allocatable, dimension(:,:) :: vnh4up_plant
     real, allocatable, dimension(:,:) :: vno3up_plant
     real, allocatable, dimension(:,:) :: vpup_plant
     real, allocatable, dimension(:,:) :: plant_input_C_pom
     real, allocatable, dimension(:) :: plant_input_C_dom
     real, allocatable, dimension(:,:) :: plant_input_N_pom
     real, allocatable, dimension(:) :: plant_input_N_dom
     real, allocatable, dimension(:,:) :: plant_input_P_pom
     real, allocatable, dimension(:) :: plant_input_P_dom
  end type plant_vars

  type mend_vars
     type(org_vars) :: cvars
     type(org_vars) :: nvars
     type(org_vars) :: pvars
     type(inorg_vars) :: invars
     type(ext_fluxes) :: fluxes
     type(plant_vars) :: plvars
  end type mend_vars

  ! Exchanges between SOM and litter
  type exchange_vars
     real, dimension(npom) :: pom_c
     real, dimension(npom) :: pom_n
     real, dimension(npom) :: pom_p
     real :: dom_c
     real :: dom_n
     real :: dom_p
     real :: nh4
     real :: no3
     real :: psol
  end type exchange_vars

  type mend_model
     type(mend_vars) :: som
!     type(mend_vars) :: litt
!     type(mend_vars), dimension(nwood) :: wood
     type(exchange_vars) :: exchange
     real, allocatable, dimension(:) :: bulk_den  ! kg/m3
     real, allocatable, dimension(:) :: pH
  end type mend_model

Contains

  subroutine mend_allocate(mend, npatches)
    implicit none
    integer, intent(in) :: npatches
    type(mend_model) :: mend
!    integer :: iwood

    call mend_allocate_type(mend%som, npatches)
!    call mend_allocate_type(mend%litt, npatches)
!    do iwood = 1, nwood
!       call mend_allocate_type(mend%wood(iwood), npatches)
!    enddo
!    call mend_allocate_exchange(mend%exchange, npatches)

    allocate(mend%bulk_den(npatches))
    allocate(mend%pH(npatches))

    return
  end subroutine mend_allocate

  subroutine mend_allocate_type(mvars, npatches)
    use ed_max_dims, only: n_pft
    implicit none
    integer, intent(in) :: npatches
    type(mend_vars) :: mvars

    call mend_allocate_type_org(mvars%cvars, npatches)
    call mend_allocate_type_org(mvars%nvars, npatches)
    call mend_allocate_type_org(mvars%pvars, npatches)

    allocate(mvars%invars%nh4(npatches))
    allocate(mvars%invars%no3(npatches))
    allocate(mvars%invars%plab(npatches))
    allocate(mvars%invars%psol(npatches))
    allocate(mvars%invars%pocc(npatches))
    allocate(mvars%invars%psec(npatches))
    allocate(mvars%invars%ppar(npatches))

    allocate(mvars%fluxes%co2_lost(npatches))
    allocate(mvars%fluxes%ngas_lost(npatches))
    allocate(mvars%fluxes%nh4_dep(npatches))
    allocate(mvars%fluxes%no3_dep(npatches))
    allocate(mvars%fluxes%ppar_dep(npatches))
    allocate(mvars%fluxes%c_leach(npatches))
    allocate(mvars%fluxes%n_leach(npatches))
    allocate(mvars%fluxes%p_leach(npatches))
    allocate(mvars%fluxes%nh4_bnf(npatches))
    allocate(mvars%fluxes%nh4_plant(n_pft,npatches))
    allocate(mvars%fluxes%no3_plant(n_pft,npatches))
    allocate(mvars%fluxes%p_plant(n_pft,npatches))
    allocate(mvars%fluxes%nh4_plant_sum(npatches))
    allocate(mvars%fluxes%no3_plant_sum(npatches))
    allocate(mvars%fluxes%p_plant_sum(npatches))
    allocate(mvars%fluxes%nmin(npatches))
    allocate(mvars%fluxes%nitr(npatches))

    allocate(mvars%plvars%enz_plant_n(n_pft,npatches))
    allocate(mvars%plvars%enz_plant_p(n_pft,npatches))
    allocate(mvars%plvars%vnh4up_plant(n_pft,npatches))
    allocate(mvars%plvars%vno3up_plant(n_pft,npatches))
    allocate(mvars%plvars%vpup_plant(n_pft,npatches))

    allocate(mvars%plvars%plant_input_C_pom(npom,npatches))
    allocate(mvars%plvars%plant_input_C_dom(npatches))
    allocate(mvars%plvars%plant_input_N_pom(npom,npatches))
    allocate(mvars%plvars%plant_input_N_dom(npatches))
    allocate(mvars%plvars%plant_input_P_pom(npom,npatches))
    allocate(mvars%plvars%plant_input_P_dom(npatches))

    return
  end subroutine mend_allocate_type

  subroutine mend_allocate_type_org(ovars, npatches)
    implicit none
    integer, intent(in) :: npatches
    type(org_vars) :: ovars

    allocate(ovars%pom(npom,npatches))
    allocate(ovars%enz_pom(npom,npatches))

    allocate(ovars%dom(npatches))
    allocate(ovars%mom(npatches))
    allocate(ovars%qom(npatches))
    allocate(ovars%enz_mom(npatches))
    allocate(ovars%amb(npatches))
    allocate(ovars%dmb(npatches))
    allocate(ovars%enz_ptase(npatches))

    return
  end subroutine mend_allocate_type_org

!  subroutine mend_allocate_exchange(exchange, npatches)
!    implicit none
!    type(exchange_vars) :: exchange
!    integer, intent(in) :: npatches

!    allocate(exchange%pom_c(npom,npatches))
!    allocate(exchange%pom_n(npom,npatches))
!    allocate(exchange%pom_p(npom,npatches))

!    allocate(exchange%dom_c(npatches))
!    allocate(exchange%dom_n(npatches))
!    allocate(exchange%dom_p(npatches))

!    allocate(exchange%nh4(npatches))
!    allocate(exchange%no3(npatches))
!    allocate(exchange%psol(npatches))

!    return
!  end subroutine mend_allocate_exchange

  subroutine mend_deallocate(mend)
    implicit none
    type(mend_model) :: mend

    call mend_deallocate_type(mend%som)

    deallocate(mend%bulk_den)
    deallocate(mend%pH)

    return
  end subroutine mend_deallocate

  subroutine mend_deallocate_type(mvars)
    implicit none
    type(mend_vars) :: mvars

    call mend_deallocate_type_org(mvars%cvars)
    call mend_deallocate_type_org(mvars%nvars)
    call mend_deallocate_type_org(mvars%pvars)

    deallocate(mvars%invars%nh4)
    deallocate(mvars%invars%no3)
    deallocate(mvars%invars%plab)
    deallocate(mvars%invars%psol)
    deallocate(mvars%invars%pocc)
    deallocate(mvars%invars%psec)
    deallocate(mvars%invars%ppar)

    deallocate(mvars%fluxes%co2_lost)
    deallocate(mvars%fluxes%ngas_lost)
    deallocate(mvars%fluxes%nh4_dep)
    deallocate(mvars%fluxes%no3_dep)
    deallocate(mvars%fluxes%ppar_dep)
    deallocate(mvars%fluxes%c_leach)
    deallocate(mvars%fluxes%n_leach)
    deallocate(mvars%fluxes%p_leach)
    deallocate(mvars%fluxes%nh4_bnf)
    deallocate(mvars%fluxes%nh4_plant)
    deallocate(mvars%fluxes%no3_plant)
    deallocate(mvars%fluxes%p_plant)
    deallocate(mvars%fluxes%nh4_plant_sum)
    deallocate(mvars%fluxes%no3_plant_sum)
    deallocate(mvars%fluxes%p_plant_sum)
    deallocate(mvars%fluxes%nmin)
    deallocate(mvars%fluxes%nitr)

    deallocate(mvars%plvars%enz_plant_n)
    deallocate(mvars%plvars%enz_plant_p)
    deallocate(mvars%plvars%vnh4up_plant)
    deallocate(mvars%plvars%vno3up_plant)
    deallocate(mvars%plvars%vpup_plant)

    deallocate(mvars%plvars%plant_input_C_pom)
    deallocate(mvars%plvars%plant_input_C_dom)
    deallocate(mvars%plvars%plant_input_N_pom)
    deallocate(mvars%plvars%plant_input_N_dom)
    deallocate(mvars%plvars%plant_input_P_pom)
    deallocate(mvars%plvars%plant_input_P_dom)

    return
  end subroutine mend_deallocate_type

  subroutine mend_deallocate_type_org(ovars)
    implicit none
    type(org_vars) :: ovars

    deallocate(ovars%pom)
    deallocate(ovars%enz_pom)

    deallocate(ovars%dom)
    deallocate(ovars%mom)
    deallocate(ovars%qom)
    deallocate(ovars%enz_mom)
    deallocate(ovars%amb)
    deallocate(ovars%dmb)
    deallocate(ovars%enz_ptase)

    return
  end subroutine mend_deallocate_type_org

!  subroutine mend_deallocate_exchange(exchange)
!    implicit none
!    type(exchange_vars) :: exchange

!    deallocate(exchange%pom_c)
!    deallocate(exchange%pom_n)
!    deallocate(exchange%pom_p)

!    deallocate(exchange%dom_c)
!    deallocate(exchange%dom_n)
!    deallocate(exchange%dom_p)

!    deallocate(exchange%nh4)
!    deallocate(exchange%no3)
!    deallocate(exchange%psol)

!    return
!  end subroutine mend_deallocate_exchange

  subroutine copy_mendtype(imend, omend, iipa, oipa)
    implicit none
    type(mend_model) :: imend
    type(mend_model) :: omend
    integer, intent(in) :: iipa
    integer, intent(in) :: oipa
!    integer :: iwood

    call copy_mendtype_type(imend%som,omend%som,iipa,oipa)
!    call copy_mendtype_type(imend%litt,omend%litt,iipa,oipa)
!    do iwood = 1, nwood
!       call copy_mendtype_type(imend%wood(iwood),omend%wood(iwood),iipa,oipa)
!    enddo
!    call copy_mendtype_exchange(imend%exchange, omend%exchange, ipa)

    omend%bulk_den(oipa) = imend%bulk_den(iipa)
    omend%pH(oipa) = imend%pH(iipa)

    return
  end subroutine copy_mendtype

  subroutine copy_mendtype_type(imvars, omvars, iipa, oipa)
    use ed_max_dims, only: n_pft
    implicit none
    type(mend_vars) :: imvars
    type(mend_vars) :: omvars
    integer, intent(in) :: iipa
    integer, intent(in) :: oipa
    integer :: ipom
    integer :: ipft

    call copy_mendtype_type_org(imvars%cvars, omvars%cvars, iipa, oipa)
    call copy_mendtype_type_org(imvars%nvars, omvars%nvars, iipa, oipa)
    call copy_mendtype_type_org(imvars%pvars, omvars%pvars, iipa, oipa)

    omvars%invars%nh4(oipa) = imvars%invars%nh4(iipa)
    omvars%invars%no3(oipa) = imvars%invars%no3(iipa)
    omvars%invars%plab(oipa) = imvars%invars%plab(iipa)
    omvars%invars%psol(oipa) = imvars%invars%psol(iipa)
    omvars%invars%pocc(oipa) = imvars%invars%pocc(iipa)
    omvars%invars%psec(oipa) = imvars%invars%psec(iipa)
    omvars%invars%ppar(oipa) = imvars%invars%ppar(iipa)

    omvars%fluxes%co2_lost(oipa) = imvars%fluxes%co2_lost(iipa)
    omvars%fluxes%ngas_lost(oipa) = imvars%fluxes%ngas_lost(iipa)
    omvars%fluxes%nh4_dep(oipa) = imvars%fluxes%nh4_dep(iipa)
    omvars%fluxes%no3_dep(oipa) = imvars%fluxes%no3_dep(iipa)
    omvars%fluxes%ppar_dep(oipa) = imvars%fluxes%ppar_dep(iipa)
    omvars%fluxes%c_leach(oipa) = imvars%fluxes%c_leach(iipa)
    omvars%fluxes%n_leach(oipa) = imvars%fluxes%n_leach(iipa)
    omvars%fluxes%p_leach(oipa) = imvars%fluxes%p_leach(iipa)
    omvars%fluxes%nh4_bnf(oipa) = imvars%fluxes%nh4_bnf(iipa)
    omvars%fluxes%nmin(oipa) = imvars%fluxes%nmin(iipa)
    omvars%fluxes%nitr(oipa) = imvars%fluxes%nitr(iipa)

    do ipft = 1, n_pft
       omvars%plvars%enz_plant_n(ipft,oipa) = imvars%plvars%enz_plant_n(ipft,iipa)
       omvars%plvars%enz_plant_p(ipft,oipa) = imvars%plvars%enz_plant_p(ipft,iipa)
       omvars%fluxes%nh4_plant(ipft,oipa) = imvars%fluxes%nh4_plant(ipft,iipa)
       omvars%fluxes%no3_plant(ipft,oipa) = imvars%fluxes%no3_plant(ipft,iipa)
       omvars%fluxes%p_plant(ipft,oipa) = imvars%fluxes%p_plant(ipft,iipa)
       omvars%plvars%vnh4up_plant(ipft,oipa) = imvars%plvars%vnh4up_plant(ipft,iipa)
       omvars%plvars%vno3up_plant(ipft,oipa) = imvars%plvars%vno3up_plant(ipft,iipa)
       omvars%plvars%vpup_plant(ipft,oipa) = imvars%plvars%vpup_plant(ipft,iipa)
    enddo

    omvars%fluxes%nh4_plant_sum(oipa) = imvars%fluxes%nh4_plant_sum(iipa)
    omvars%fluxes%no3_plant_sum(oipa) = imvars%fluxes%no3_plant_sum(iipa)
    omvars%fluxes%p_plant_sum(oipa) = imvars%fluxes%p_plant_sum(iipa)

    do ipom = 1, npom
       omvars%plvars%plant_input_C_pom(ipom,oipa) =   &
            imvars%plvars%plant_input_C_pom(ipom,iipa)
       omvars%plvars%plant_input_N_pom(ipom,oipa) =   &
            imvars%plvars%plant_input_N_pom(ipom,iipa)
       omvars%plvars%plant_input_P_pom(ipom,oipa) =   &
            imvars%plvars%plant_input_P_pom(ipom,iipa)
    enddo
    omvars%plvars%plant_input_C_dom(oipa) =   &
         imvars%plvars%plant_input_C_dom(iipa)
    omvars%plvars%plant_input_N_dom(oipa) =   &
         imvars%plvars%plant_input_N_dom(iipa)
    omvars%plvars%plant_input_P_dom(oipa) =   &
         imvars%plvars%plant_input_P_dom(iipa)

    return
  end subroutine copy_mendtype_type

  subroutine copy_mendtype_type_org(iovars, oovars, iipa, oipa)
    implicit none
    type(org_vars) :: iovars
    type(org_vars) :: oovars
    integer, intent(in) :: iipa
    integer, intent(in) :: oipa
    integer :: ipom

    do ipom = 1, npom
       oovars%pom(ipom,oipa) = iovars%pom(ipom,iipa)
       oovars%enz_pom(ipom,oipa) = iovars%enz_pom(ipom,iipa)
    enddo

    oovars%dom(oipa) = iovars%dom(iipa)
    oovars%mom(oipa) = iovars%mom(iipa)
    oovars%qom(oipa) = iovars%qom(iipa)
    oovars%enz_mom(oipa) = iovars%enz_mom(iipa)
    oovars%amb(oipa) = iovars%amb(iipa)
    oovars%dmb(oipa) = iovars%dmb(iipa)
    oovars%enz_ptase(oipa) = iovars%enz_ptase(iipa)
    
    return
  end subroutine copy_mendtype_type_org

!  subroutine copy_mendtype_exchange(iexchange, oexchange, ipa)
!    implicit none
!    type(exchange_vars) :: iexchange
!    type(exchange_vars) :: oexchange
!    integer, intent(in) :: ipa
!    integer :: ipom

!    do ipom = 1, npom
!       oexchange%pom_c(ipom,ipa) = iexchange%pom_c(ipom,ipa)
!       oexchange%pom_n(ipom,ipa) = iexchange%pom_n(ipom,ipa)
!       oexchange%pom_p(ipom,ipa) = iexchange%pom_p(ipom,ipa)
!    enddo

!    oexchange%dom_c(ipa) = iexchange%dom_c(ipa)
!    oexchange%dom_n(ipa) = iexchange%dom_n(ipa)
!    oexchange%dom_p(ipa) = iexchange%dom_p(ipa)

!    oexchange%nh4(ipa) = iexchange%nh4(ipa)
!    oexchange%no3(ipa) = iexchange%no3(ipa)
!    oexchange%psol(ipa) = iexchange%psol(ipa)

!    return
!  end subroutine copy_mendtype_exchange

  subroutine copy_mendtype_mask(imend, masksz, logmask, omend, inc)
    implicit none
    
    type(mend_model) :: imend
    integer, intent(in) :: inc
    integer, intent(in) :: masksz
    logical, dimension(inc), intent(in) :: logmask
    type(mend_model) :: omend
!    integer :: iwood

    call copy_mendtype_mask_type(imend%som, omend%som, masksz, logmask, inc)
!    call copy_mendtype_mask_type(imend%litt, omend%litt, masksz, logmask, inc)
!    do iwood = 1, nwood
!       call copy_mendtype_mask_type(imend%wood(iwood), omend%wood(iwood), masksz, logmask, inc)
!    enddo
!    call copy_mendtype_mask_exchange(imend%exchange, omend%exchange, masksz, logmask, inc)

    omend%bulk_den(1:masksz) = pack(imend%bulk_den, logmask)
    omend%pH(1:masksz) = pack(imend%pH, logmask)

    return
  end subroutine copy_mendtype_mask

  subroutine copy_mendtype_mask_type(imvars, omvars, masksz, logmask, inc)
    use ed_max_dims, only: n_pft
    implicit none
    type(mend_vars) :: imvars
    type(mend_vars) :: omvars
    integer, intent(in) :: inc
    integer, intent(in) :: masksz
    logical, dimension(inc), intent(in) :: logmask
    integer :: ipom
    integer :: ipft

    call copy_mendtype_mask_type_org(imvars%cvars, omvars%cvars, masksz, logmask, inc)
    call copy_mendtype_mask_type_org(imvars%nvars, omvars%nvars, masksz, logmask, inc)
    call copy_mendtype_mask_type_org(imvars%pvars, omvars%pvars, masksz, logmask, inc)

    omvars%invars%nh4(1:masksz) = pack(imvars%invars%nh4,logmask)
    omvars%invars%no3(1:masksz) = pack(imvars%invars%no3,logmask)
    omvars%invars%plab(1:masksz) = pack(imvars%invars%plab,logmask)
    omvars%invars%psol(1:masksz) = pack(imvars%invars%psol,logmask)
    omvars%invars%pocc(1:masksz) = pack(imvars%invars%pocc,logmask)
    omvars%invars%psec(1:masksz) = pack(imvars%invars%psec,logmask)
    omvars%invars%ppar(1:masksz) = pack(imvars%invars%ppar,logmask)

    omvars%fluxes%co2_lost(1:masksz) = pack(imvars%fluxes%co2_lost,logmask)
    omvars%fluxes%ngas_lost(1:masksz) = pack(imvars%fluxes%ngas_lost,logmask)
    omvars%fluxes%nh4_dep(1:masksz) = pack(imvars%fluxes%nh4_dep,logmask)
    omvars%fluxes%no3_dep(1:masksz) = pack(imvars%fluxes%no3_dep,logmask)
    omvars%fluxes%ppar_dep(1:masksz) = pack(imvars%fluxes%ppar_dep,logmask)
    omvars%fluxes%c_leach(1:masksz) = pack(imvars%fluxes%c_leach,logmask)
    omvars%fluxes%n_leach(1:masksz) = pack(imvars%fluxes%n_leach,logmask)
    omvars%fluxes%p_leach(1:masksz) = pack(imvars%fluxes%p_leach,logmask)
    omvars%fluxes%nh4_bnf(1:masksz) = pack(imvars%fluxes%nh4_bnf,logmask)
    omvars%fluxes%nmin(1:masksz) = pack(imvars%fluxes%nmin,logmask)
    omvars%fluxes%nitr(1:masksz) = pack(imvars%fluxes%nitr,logmask)

    do ipft = 1, n_pft
       omvars%plvars%enz_plant_n(ipft,1:masksz) = pack(imvars%plvars%enz_plant_n(ipft,:),logmask)
       omvars%plvars%enz_plant_p(ipft,1:masksz) = pack(imvars%plvars%enz_plant_p(ipft,:),logmask)
       omvars%fluxes%nh4_plant(ipft,1:masksz) = pack(imvars%fluxes%nh4_plant(ipft,:),logmask)
       omvars%fluxes%no3_plant(ipft,1:masksz) = pack(imvars%fluxes%no3_plant(ipft,:),logmask)
       omvars%fluxes%p_plant(ipft,1:masksz) = pack(imvars%fluxes%p_plant(ipft,:),logmask)
       omvars%plvars%vnh4up_plant(ipft,1:masksz) = pack(imvars%plvars%vnh4up_plant(ipft,:),logmask)
       omvars%plvars%vno3up_plant(ipft,1:masksz) = pack(imvars%plvars%vno3up_plant(ipft,:),logmask)
       omvars%plvars%vpup_plant(ipft,1:masksz) = pack(imvars%plvars%vpup_plant(ipft,:),logmask)
    enddo

    omvars%fluxes%nh4_plant_sum(1:masksz) = pack(imvars%fluxes%nh4_plant_sum,logmask)
    omvars%fluxes%no3_plant_sum(1:masksz) = pack(imvars%fluxes%no3_plant_sum,logmask)
    omvars%fluxes%p_plant_sum(1:masksz) = pack(imvars%fluxes%p_plant_sum,logmask)

    do ipom = 1, npom
       omvars%plvars%plant_input_C_pom(ipom,1:masksz) =   &
            pack(imvars%plvars%plant_input_C_pom(ipom,:),logmask)
       omvars%plvars%plant_input_N_pom(ipom,1:masksz) =   &
            pack(imvars%plvars%plant_input_N_pom(ipom,:),logmask)
       omvars%plvars%plant_input_P_pom(ipom,1:masksz) =   &
            pack(imvars%plvars%plant_input_P_pom(ipom,:),logmask)
    enddo
    omvars%plvars%plant_input_C_dom(1:masksz) =   &
         pack(imvars%plvars%plant_input_C_dom,logmask)
    omvars%plvars%plant_input_N_dom(1:masksz) =   &
         pack(imvars%plvars%plant_input_N_dom,logmask)
    omvars%plvars%plant_input_P_dom(1:masksz) =   &
         pack(imvars%plvars%plant_input_P_dom,logmask)

    return
  end subroutine copy_mendtype_mask_type

  subroutine copy_mendtype_mask_type_org(iovars, oovars, masksz, logmask, inc)
    implicit none
    type(org_vars) :: iovars
    type(org_vars) :: oovars
    integer, intent(in) :: inc
    integer, intent(in) :: masksz
    logical, dimension(inc), intent(in) :: logmask
    integer :: ipom
 
    do ipom = 1, npom
       oovars%pom(ipom,1:masksz) = pack(iovars%pom(ipom,:),logmask)
       oovars%enz_pom(ipom,1:masksz) = pack(iovars%enz_pom(ipom,:),logmask)
    enddo

    oovars%dom(1:masksz) = pack(iovars%dom,logmask)
    oovars%mom(1:masksz) = pack(iovars%mom,logmask)
    oovars%qom(1:masksz) = pack(iovars%qom,logmask)
    oovars%enz_mom(1:masksz) = pack(iovars%enz_mom,logmask)
    oovars%amb(1:masksz) = pack(iovars%amb,logmask)
    oovars%dmb(1:masksz) = pack(iovars%dmb,logmask)
    oovars%enz_ptase(1:masksz) = pack(iovars%enz_ptase,logmask)
   
    return
  end subroutine copy_mendtype_mask_type_org

!  subroutine copy_mendtype_mask_exchange(iexchange, oexchange, masksz, logmask, inc)
!    implicit none
!    type(exchange_vars) :: iexchange
!    type(exchange_vars) :: oexchange
!    integer, intent(in) :: masksz
!    logical, dimension(masksz), intent(in) :: logmask
!    integer, intent(in) :: inc
!    integer :: ipom

!    do ipom = 1, npom
!       oexchange%pom_c(ipom,1:inc) = pack(iexchange%pom_c(ipom,:),logmask)
!       oexchange%pom_n(ipom,1:inc) = pack(iexchange%pom_n(ipom,:),logmask)
!       oexchange%pom_p(ipom,1:inc) = pack(iexchange%pom_p(ipom,:),logmask)
!    enddo

!    oexchange%dom_c(1:inc) = pack(iexchange%dom_c,logmask)
!    oexchange%dom_n(1:inc) = pack(iexchange%dom_n,logmask)
!    oexchange%dom_p(1:inc) = pack(iexchange%dom_p,logmask)

!    oexchange%nh4(1:inc) = pack(iexchange%nh4,logmask)
!    oexchange%no3(1:inc) = pack(iexchange%no3,logmask)
!    oexchange%psol(1:inc) = pack(iexchange%psol,logmask)

!    return
!  end subroutine copy_mendtype_mask_exchange

  subroutine filltab_mendtype(nvar, npts, mend, igr, init, paglob_id, var_len, &
              var_len_global, max_ptrs, time_code)
    use ed_var_tables, only: vtable_edio_r, metadata_edio
    implicit none
    integer, intent(inout) :: nvar
    integer, intent(in) :: npts
    type(mend_model) :: mend
    integer, intent(in) :: igr
    integer, intent(in) :: init
    integer, intent(in) :: paglob_id
    integer, intent(in) :: var_len
    integer, intent(in) :: var_len_global
    integer, intent(in) :: max_ptrs
    integer :: iwood
    integer :: npatches
    character(len=64) :: type_string
    character(len=4) :: time_code

    npatches = npts

    if(time_code == 'mont')then
       type_string = 'MMEAN_MEND_'
    else
       type_string = 'MEND_'
    endif
    nvar=nvar+1
    call vtable_edio_r(npatches,mend%bulk_den,nvar,igr,init,paglob_id, &
         var_len,var_len_global,max_ptrs,trim(type_string)//'BULK_DEN :31:hist:'//trim(time_code)) 
    call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    nvar=nvar+1
    call vtable_edio_r(npatches,mend%pH,nvar,igr,init,paglob_id, &
         var_len,var_len_global,max_ptrs,trim(type_string)//'pH :31:hist:'//trim(time_code))
    call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 

    if(time_code == 'mont')then
       type_string = 'MMEAN_MEND_SOM_'
    else
       type_string = 'MEND_SOM_'
    endif
    call filltab_mendtype_type(nvar, npatches, mend%som, igr, init, paglob_id, var_len, &
              var_len_global, max_ptrs, type_string, time_code)

!    type_string = 'MEND_LITT_'
!    call filltab_mendtype_type(nvar, npatches, mend%litt, igr, init, paglob_id, var_len, &
!              var_len_global, max_ptrs, type_string)

!    do iwood = 1, nwood
!       write(type_string,'(a,i1.1,a)')'MEND_WOOD',iwood,'_'
!       call filltab_mendtype_type(nvar, npatches, mend%wood(iwood), igr, init, paglob_id,  &
!            var_len, var_len_global, max_ptrs, type_string)
!    enddo

!    type_string = 'MEND_EX_'
!    call filltab_mendtype_exchange(nvar, npatches, mend%exchange, igr, init, paglob_id,   &
!         var_len, var_len_global, max_ptrs, type_string)

    return
  end subroutine filltab_mendtype
  
  subroutine filltab_mendtype_type(nvar, npatches, mvars, igr, init, paglob_id, var_len, &
              var_len_global, max_ptrs, type_string, time_code)
    use ed_var_tables, only: vtable_edio_r, metadata_edio
    use ed_max_dims, only: n_pft
    implicit none
    integer, intent(inout) :: nvar
    integer, intent(in) :: npatches
    type(mend_vars) :: mvars
    integer, intent(in) :: igr
    integer, intent(in) :: init
    integer, intent(in) :: paglob_id
    integer, intent(in) :: var_len
    integer, intent(in) :: var_len_global
    integer, intent(in) :: max_ptrs
    character(len=*), intent(in) :: type_string
    character(len=4) :: time_code
    real, dimension(npatches) :: uptake_sum
    integer :: ipa, ipft

    
    call filltab_mendtype_type_org(nvar, npatches, mvars%cvars, igr, init, paglob_id, var_len, &
              var_len_global, max_ptrs, trim(type_string)//'C_', time_code)
    call filltab_mendtype_type_org(nvar, npatches, mvars%nvars, igr, init, paglob_id, var_len, &
              var_len_global, max_ptrs, trim(type_string)//'N_', time_code)
    call filltab_mendtype_type_org(nvar, npatches, mvars%pvars, igr, init, paglob_id, var_len, &
              var_len_global, max_ptrs, trim(type_string)//'P_', time_code)

    nvar=nvar+1
    call vtable_edio_r(npatches,mvars%invars%nh4,nvar,igr,init,paglob_id, &
         var_len,var_len_global,max_ptrs,trim(type_string)//'IN_NH4 :31:hist:'//trim(time_code))
    call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 

    nvar=nvar+1
    call vtable_edio_r(npatches,mvars%invars%no3,nvar,igr,init,paglob_id, &
         var_len,var_len_global,max_ptrs,trim(type_string)//'IN_NO3 :31:hist:'//trim(time_code))
    call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 

    nvar=nvar+1
    call vtable_edio_r(npatches,mvars%invars%plab,nvar,igr,init,paglob_id, &
         var_len,var_len_global,max_ptrs,trim(type_string)//'IN_PLAB :31:hist:'//trim(time_code))
    call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 

    nvar=nvar+1
    call vtable_edio_r(npatches,mvars%invars%psol,nvar,igr,init,paglob_id, &
         var_len,var_len_global,max_ptrs,trim(type_string)//'IN_PSOL :31:hist:'//trim(time_code))
    call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 

    nvar=nvar+1
    call vtable_edio_r(npatches,mvars%invars%pocc,nvar,igr,init,paglob_id, &
         var_len,var_len_global,max_ptrs,trim(type_string)//'IN_POCC :31:hist:'//trim(time_code))
    call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 

    nvar=nvar+1
    call vtable_edio_r(npatches,mvars%invars%psec,nvar,igr,init,paglob_id, &
         var_len,var_len_global,max_ptrs,trim(type_string)//'IN_PSEC :31:hist:'//trim(time_code))
    call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 

    nvar=nvar+1
    call vtable_edio_r(npatches,mvars%invars%ppar,nvar,igr,init,paglob_id, &
         var_len,var_len_global,max_ptrs,trim(type_string)//'IN_PPAR :31:hist:'//trim(time_code))
    call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 

    nvar=nvar+1
    call vtable_edio_r(npatches,mvars%fluxes%co2_lost,nvar,igr,init,paglob_id, &
         var_len,var_len_global,max_ptrs,trim(type_string)//'FL_CO2LOSS :31:hist:'//trim(time_code))
    call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 

    nvar=nvar+1
    call vtable_edio_r(npatches,mvars%fluxes%ngas_lost,nvar,igr,init,paglob_id, &
         var_len,var_len_global,max_ptrs,trim(type_string)//'FL_NGASLOSS :31:hist:'//trim(time_code))
    call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 

    nvar=nvar+1
    call vtable_edio_r(npatches,mvars%fluxes%nh4_dep,nvar,igr,init,paglob_id, &
         var_len,var_len_global,max_ptrs,trim(type_string)//'FL_NH4DEP :31:hist:'//trim(time_code))
    call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 

    nvar=nvar+1
    call vtable_edio_r(npatches,mvars%fluxes%no3_dep,nvar,igr,init,paglob_id, &
         var_len,var_len_global,max_ptrs,trim(type_string)//'FL_NO3DEP :31:hist:'//trim(time_code))
    call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 

    nvar=nvar+1
    call vtable_edio_r(npatches,mvars%fluxes%ppar_dep,nvar,igr,init,paglob_id, &
         var_len,var_len_global,max_ptrs,trim(type_string)//'FL_PPARDEP :31:hist:'//trim(time_code))
    call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 

    nvar=nvar+1
    call vtable_edio_r(npatches,mvars%fluxes%c_leach,nvar,igr,init,paglob_id, &
         var_len,var_len_global,max_ptrs,trim(type_string)//'FL_CLEACH :31:hist:'//trim(time_code))
    call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 

    nvar=nvar+1
    call vtable_edio_r(npatches,mvars%fluxes%n_leach,nvar,igr,init,paglob_id, &
         var_len,var_len_global,max_ptrs,trim(type_string)//'FL_NLEACH :31:hist:'//trim(time_code))
    call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 

    nvar=nvar+1
    call vtable_edio_r(npatches,mvars%fluxes%p_leach,nvar,igr,init,paglob_id, &
         var_len,var_len_global,max_ptrs,trim(type_string)//'FL_PLEACH :31:hist:'//trim(time_code))
    call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 

    nvar=nvar+1
    call vtable_edio_r(npatches,mvars%fluxes%nh4_bnf,nvar,igr,init,paglob_id, &
         var_len,var_len_global,max_ptrs,trim(type_string)//'FL_NH4BNF :31:hist:'//trim(time_code))
    call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 

    nvar=nvar+1
    call vtable_edio_r(npatches,mvars%fluxes%nh4_plant_sum,nvar,igr,init,paglob_id, &
         var_len,var_len_global,max_ptrs,trim(type_string)//'FL_NH4PL :31:hist:'//trim(time_code))
    call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 

    nvar=nvar+1
    call vtable_edio_r(npatches,mvars%fluxes%no3_plant_sum,nvar,igr,init,paglob_id, &
         var_len,var_len_global,max_ptrs,trim(type_string)//'FL_NO3PL :31:hist:'//trim(time_code))
    call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 

    nvar=nvar+1
    call vtable_edio_r(npatches,mvars%fluxes%nmin,nvar,igr,init,paglob_id, &
         var_len,var_len_global,max_ptrs,trim(type_string)//'FL_NMIN :31:hist:'//trim(time_code))
    call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 

    nvar=nvar+1
    call vtable_edio_r(npatches,mvars%fluxes%nitr,nvar,igr,init,paglob_id, &
         var_len,var_len_global,max_ptrs,trim(type_string)//'FL_NITR :31:hist:'//trim(time_code))
    call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 

    nvar=nvar+1
    call vtable_edio_r(npatches,mvars%fluxes%p_plant_sum,nvar,igr,init,paglob_id, &
         var_len,var_len_global,max_ptrs,trim(type_string)//'FL_PPL :31:hist:'//trim(time_code))
    call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 

!    nvar=nvar+1
!    call vtable_edio_r(npatches,mvars%plvars%enz_plant_n,nvar,igr,init,paglob_id, &
!         var_len,var_len_global,max_ptrs,trim(type_string)//'PL_ENZN :31:hist:mont:year') 
!    call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 

!    nvar=nvar+1
!    call vtable_edio_r(npatches,mvars%plvars%enz_plant_p,nvar,igr,init,paglob_id, &
!         var_len,var_len_global,max_ptrs,trim(type_string)//'PL_ENZP :31:hist:mont:year') 
!    call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 

!    nvar=nvar+1
!    call vtable_edio_r(npatches,mvars%plvars%vnh4up_plant,nvar,igr,init,paglob_id, &
!         var_len,var_len_global,max_ptrs,trim(type_string)//'PL_VNH4UP :31:hist:mont:year') 
!    call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 

!    nvar=nvar+1
!    call vtable_edio_r(npatches,mvars%plvars%vno3up_plant,nvar,igr,init,paglob_id, &
!         var_len,var_len_global,max_ptrs,trim(type_string)//'PL_VNO3UP :31:hist:mont:year') 
!    call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 

!    nvar=nvar+1
!    call vtable_edio_r(npatches,mvars%plvars%vpup_plant,nvar,igr,init,paglob_id, &
!         var_len,var_len_global,max_ptrs,trim(type_string)//'PL_VPUP :31:hist:mont:year') 
!    call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 

    nvar=nvar+1
    call vtable_edio_r(npatches*npom,mvars%plvars%plant_input_C_pom,nvar,igr,init,paglob_id, &
         var_len,var_len_global,max_ptrs,trim(type_string)//'PL_POMINPC :311:hist:'//trim(time_code))
    call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 

    nvar=nvar+1
    call vtable_edio_r(npatches*npom,mvars%plvars%plant_input_N_pom,nvar,igr,init,paglob_id, &
         var_len,var_len_global,max_ptrs,trim(type_string)//'PL_POMINPN :311:hist:'//trim(time_code))
    call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 

    nvar=nvar+1
    call vtable_edio_r(npatches*npom,mvars%plvars%plant_input_P_pom,nvar,igr,init,paglob_id, &
         var_len,var_len_global,max_ptrs,trim(type_string)//'PL_POMINPP :311:hist:'//trim(time_code))
    call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    
    nvar=nvar+1
    call vtable_edio_r(npatches,mvars%plvars%plant_input_C_dom,nvar,igr,init,paglob_id, &
         var_len,var_len_global,max_ptrs,trim(type_string)//'PL_DOMINPC :31:hist:'//trim(time_code))
    call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 

    nvar=nvar+1
    call vtable_edio_r(npatches,mvars%plvars%plant_input_N_dom,nvar,igr,init,paglob_id, &
         var_len,var_len_global,max_ptrs,trim(type_string)//'PL_DOMINPN :31:hist:'//trim(time_code))
    call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 

    nvar=nvar+1
    call vtable_edio_r(npatches,mvars%plvars%plant_input_P_dom,nvar,igr,init,paglob_id, &
         var_len,var_len_global,max_ptrs,trim(type_string)//'PL_DOMINPP :31:hist:'//trim(time_code))
    call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 

    return
  end subroutine filltab_mendtype_type

  subroutine filltab_mendtype_type_org(nvar, npatches, ovars, igr, init, paglob_id, var_len, &
              var_len_global, max_ptrs, type_string, time_code)
    use ed_var_tables, only: vtable_edio_r, metadata_edio
    implicit none
    integer, intent(inout) :: nvar
    integer, intent(in) :: npatches
    type(org_vars) :: ovars
    integer, intent(in) :: igr
    integer, intent(in) :: init
    integer, intent(in) :: paglob_id
    integer, intent(in) :: var_len
    integer, intent(in) :: var_len_global
    integer, intent(in) :: max_ptrs
    character(len=*), intent(in) :: type_string
    integer :: npts
    character(len=4) :: time_code

    npts = npatches * npom

    nvar=nvar+1
    call vtable_edio_r(npts,ovars%pom,nvar,igr,init,paglob_id, &
         var_len,var_len_global,max_ptrs,trim(type_string)//'POM :311:hist:'//trim(time_code))
    call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 
    
    nvar=nvar+1
    call vtable_edio_r(npts,ovars%enz_pom,nvar,igr,init,paglob_id, &
         var_len,var_len_global,max_ptrs,trim(type_string)//'ENZPOM :311:hist:'//trim(time_code))
    call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 

    npts = npatches

    nvar=nvar+1
    call vtable_edio_r(npts,ovars%dom,nvar,igr,init,paglob_id, &
         var_len,var_len_global,max_ptrs,trim(type_string)//'DOM :31:hist:'//trim(time_code))
    call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 

    nvar=nvar+1
    call vtable_edio_r(npts,ovars%mom,nvar,igr,init,paglob_id, &
         var_len,var_len_global,max_ptrs,trim(type_string)//'MOM :31:hist:'//trim(time_code))
    call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 

    nvar=nvar+1
    call vtable_edio_r(npts,ovars%qom,nvar,igr,init,paglob_id, &
         var_len,var_len_global,max_ptrs,trim(type_string)//'QOM :31:hist:'//trim(time_code))
    call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 

    nvar=nvar+1
    call vtable_edio_r(npts,ovars%enz_mom,nvar,igr,init,paglob_id, &
         var_len,var_len_global,max_ptrs,trim(type_string)//'ENZMOM :31:hist:'//trim(time_code))
    call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 

    nvar=nvar+1
    call vtable_edio_r(npts,ovars%amb,nvar,igr,init,paglob_id, &
         var_len,var_len_global,max_ptrs,trim(type_string)//'AMB :31:hist:'//trim(time_code))
    call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 

    nvar=nvar+1
    call vtable_edio_r(npts,ovars%dmb,nvar,igr,init,paglob_id, &
         var_len,var_len_global,max_ptrs,trim(type_string)//'DMB :31:hist:'//trim(time_code))
    call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 

    nvar=nvar+1
    call vtable_edio_r(npts,ovars%enz_ptase,nvar,igr,init,paglob_id, &
         var_len,var_len_global,max_ptrs,trim(type_string)//'ENZPTASE :31:hist:'//trim(time_code))
    call metadata_edio(nvar,igr,'No metadata available','[NA]','NA') 

    return
  end subroutine filltab_mendtype_type_org

  subroutine mend_zero_vars(mend, ip1, ip2)
    implicit none

    type(mend_model) :: mend
    integer :: iwood
    integer, intent(in) :: ip1
    integer, intent(in) :: ip2

    call mend_zero_rk4_type(mend%som, ip1, ip2)
!    call mend_zero_rk4_type(mend%litt, ip1, ip2)
!    do iwood = 1, nwood
!       call mend_zero_rk4_type(mend%wood(iwood), ip1, ip2)
!    enddo
    
!    mend%bulk_den = 0.
!    mend%pH = 0.
    
    return
  end subroutine mend_zero_vars

  subroutine mend_zero_rk4_type(vars, ip1, ip2)
    implicit none
    type(mend_vars) :: vars
    integer, intent(in) :: ip1
    integer, intent(in) :: ip2

    call mend_zero_rk4_type_org(vars%cvars, ip1, ip2)
    call mend_zero_rk4_type_org(vars%nvars, ip1, ip2)
    call mend_zero_rk4_type_org(vars%pvars, ip1, ip2)

    vars%fluxes%co2_lost(ip1:ip2) = 0.
    vars%fluxes%ngas_lost(ip1:ip2) = 0.
    vars%invars%nh4(ip1:ip2) = 0.
    vars%invars%no3(ip1:ip2) = 0.
    vars%invars%psol(ip1:ip2) = 0.
    vars%invars%plab(ip1:ip2) = 0.
    vars%fluxes%nh4_plant(:,ip1:ip2) = 0.
    vars%fluxes%nh4_bnf(ip1:ip2) = 0.
    vars%fluxes%no3_plant(:,ip1:ip2) = 0.
    vars%fluxes%nmin(ip1:ip2) = 0.
    vars%fluxes%nitr(ip1:ip2) = 0.
    vars%fluxes%c_leach(ip1:ip2) = 0.
    vars%fluxes%n_leach(ip1:ip2) = 0.
    vars%fluxes%p_leach(ip1:ip2) = 0.
    vars%fluxes%p_plant(:,ip1:ip2) = 0.
    vars%invars%pocc(ip1:ip2) = 0.
    vars%invars%psec(ip1:ip2) = 0.
    vars%invars%ppar(ip1:ip2) = 0.
    vars%plvars%enz_plant_n(:,ip1:ip2) = 0.
    vars%plvars%enz_plant_p(:,ip1:ip2) = 0.
    vars%plvars%vnh4up_plant(:,ip1:ip2) = 0.
    vars%plvars%vno3up_plant(:,ip1:ip2) = 0.
    vars%plvars%vpup_plant(:,ip1:ip2) = 0.

    vars%plvars%plant_input_C_pom(:,ip1:ip2) = 0.
    vars%plvars%plant_input_C_dom(ip1:ip2) = 0.
    vars%plvars%plant_input_N_pom(:,ip1:ip2) = 0.
    vars%plvars%plant_input_N_dom(ip1:ip2) = 0.
    vars%plvars%plant_input_P_pom(:,ip1:ip2) = 0.
    vars%plvars%plant_input_P_dom(ip1:ip2) = 0.

    vars%fluxes%nh4_dep(ip1:ip2) = 0.
    vars%fluxes%no3_dep(ip1:ip2) = 0.
    vars%fluxes%ppar_dep(ip1:ip2) = 0.

    return
  end subroutine mend_zero_rk4_type

  subroutine mend_zero_rk4_type_org(vars, ip1, ip2)
    implicit none

    type(org_vars) :: vars
    integer :: ipom
    integer, intent(in) :: ip1
    integer, intent(in) :: ip2

    vars%pom(:,ip1:ip2) = 0.
    vars%enz_pom(:,ip1:ip2) = 0.

    vars%dom(ip1:ip2) = 0.
    vars%mom(ip1:ip2) = 0.
    vars%enz_mom(ip1:ip2) = 0.
    vars%qom(ip1:ip2) = 0.
    vars%amb(ip1:ip2) = 0.
    vars%dmb(ip1:ip2) = 0.
    vars%enz_ptase(ip1:ip2) = 0.

    return
  end subroutine mend_zero_rk4_type_org

end Module mend_state_vars
