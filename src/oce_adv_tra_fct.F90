module oce_adv_tra_fct_interfaces
  interface
    subroutine oce_adv_tra_fct_init(mesh)
      use MOD_MESH
      use g_PARSUP
      type(t_mesh), intent(in), target  :: mesh
    end subroutine

    subroutine oce_tra_adv_fct(dttf_h, dttf_v, ttf, lo, adf_h, adf_v, mesh)
      use MOD_MESH
      use g_PARSUP
      type(t_mesh),  intent(in), target :: mesh
      real(kind=WP), intent(inout)      :: dttf_h(mesh%nl-1, myDim_nod2D+eDim_nod2D)
      real(kind=WP), intent(inout)      :: dttf_v(mesh%nl-1, myDim_nod2D+eDim_nod2D)
      real(kind=WP), intent(in)         :: ttf(mesh%nl-1, myDim_nod2D+eDim_nod2D)
      real(kind=WP), intent(in)         :: lo (mesh%nl-1, myDim_nod2D+eDim_nod2D)
      real(kind=WP), intent(inout)      :: adf_h(mesh%nl-1, myDim_edge2D)
      real(kind=WP), intent(inout)      :: adf_v(mesh%nl,  myDim_nod2D)
    end subroutine
 end interface
end module
#ifdef FESOMCUDA
module MOD_TRA_FCT_GPU
    USE, intrinsic :: ISO_C_BINDING
    type(c_ptr) :: nlevs_nod2D_gpu, nlevs_elem2D_gpu, nod_elem2D_gpu, nod_num_elem2D_gpu, elem2D_nodes_gpu,&
                   nod2D_edges_gpu, elem2D_edges_gpu, area_inv_gpu
    type(c_ptr) :: fct_lo_gpu, fct_ttf_gpu, fct_adf_v_gpu, fct_adf_h_gpu,  UV_rhs_gpu, fct_ttf_min_gpu,&
                   fct_ttf_max_gpu, fct_plus_gpu, fct_minus_gpu
    contains
    subroutine allocate_pinned_memory(arr, size_vert, size_hor)
        USE, intrinsic :: ISO_C_BINDING
        real, intent(out), pointer :: arr(:,:)
        integer, intent(in)        :: size_hor, size_vert
        type(c_ptr) :: cptr
        call allocate_pinned_doubles(cptr, size_hor * size_vert)
        call c_f_pointer(cptr, arr, (/size_vert, size_hor/))
    end subroutine
    end module
#endif
!
!
!===============================================================================
subroutine oce_adv_tra_fct_init(mesh)
    use MOD_MESH
#ifdef FESOMCUDA
    use MOD_TRA_FCT_GPU
#endif
    use O_MESH
    use o_ARRAYS
    use o_PARAM
    use g_PARSUP
    implicit none
    integer                  :: my_size, istat
    type(t_mesh), intent(in) , target :: mesh

#include "associate_mesh.h"

    my_size=myDim_nod2D+eDim_nod2D
    allocate(fct_LO(nl-1, my_size))        ! Low-order solution 
#ifdef FESOMCUDA
    call allocate_pinned_memory(adv_flux_hor, (nl-1), myDim_edge2D) ! antidiffusive hor. contributions / from edges
    call allocate_pinned_memory(adv_flux_ver, nl, myDim_nod2D)      ! antidiffusive ver. fluxes / from nodes
#else
    allocate(adv_flux_hor(nl-1,myDim_edge2D)) ! antidiffusive hor. contributions / from edges
    allocate(adv_flux_ver(nl, myDim_nod2D))   ! antidiffusive ver. fluxes / from nodes
#endif
    allocate(fct_ttf_max(nl-1, my_size),fct_ttf_min(nl-1, my_size))
    allocate(fct_plus(nl-1, my_size),fct_minus(nl-1, my_size))
    ! Initialize with zeros: 
    fct_LO=0.0_WP
    adv_flux_hor=0.0_WP
    adv_flux_ver=0.0_WP
    fct_ttf_max=0.0_WP
    fct_ttf_min=0.0_WP
    fct_plus=0.0_WP
    fct_minus=0.0_WP
#ifdef FESOMCUDA
    istat = 0
    call transfer_mesh(nlevs_nod2D_gpu, nlevels_nod2D, my_size, istat)
    if (istat /= 0) then
        write(0, *) "Error in transfer nlevels_nod2D to GPU"
    endif
    istat = 0
    call transfer_mesh(nlevs_elem2D_gpu, nlevels, myDim_elem2D, istat)
    if (istat /= 0) then
        write(0, *) "Error in transfer nlevels_elem2D to GPU"
    endif
    istat = 0
    call transfer_mesh(nod_elem2D_gpu, nod_in_elem2D, my_size * size(nod_in_elem2D, 1), istat)
    if (istat /= 0) then
        write(0, *) "Error in transfer nod_in_elem2D to GPU"
    endif
    istat = 0
    call transfer_mesh(nod_num_elem2D_gpu, nod_in_elem2D_num, my_size, istat)
    if (istat /= 0) then
        write(0, *) "Error in transfer nod_in_elem2D_num to GPU"
    endif
    istat = 0
    call transfer_mesh(elem2D_nodes_gpu, elem2D_nodes, myDim_elem2D * 3, istat)
    if (istat /= 0) then
        write(0, *) "Error in transfer elem2D_nodes to GPU"
    endif
    istat = 0
    call transfer_mesh(nod2D_edges_gpu, edges, myDim_edge2D * 2, istat)
    if (istat /= 0) then
        write(0, *) "Error in transfer edges to GPU"
    endif
    istat = 0
    call transfer_mesh(elem2D_edges_gpu, edge_tri, myDim_edge2D * 2, istat)
    if (istat /= 0) then
        write(0, *) "Error in transfer edge_tri to GPU"
    endif
    istat = 0
    call alloc_var(area_inv_gpu, area_inv, my_size * nl, istat)
    call transfer_var(area_inv_gpu, area_inv)
    if (istat /= 0) then
        write(0, *) "Error in alloc/transfer area_inv to GPU"
    endif
    istat = 0
    call alloc_var(fct_lo_gpu, fct_lo, my_size * (nl - 1), istat)
    if (istat /= 0) then
        write(0, *) "Error in alloc fct_lo to GPU"
    endif
    istat = 0
    call reserve_var(fct_ttf_gpu, my_size * (nl - 1), istat)
    if (istat /= 0) then
        write(0, *) "Error in alloc fct_ttf to GPU"
    endif
    istat = 0
    call alloc_var(fct_adf_v_gpu, adv_flux_ver, myDim_nod2D * nl, istat)
    if (istat /= 0) then
        write(0, *) "Error in alloc fct_adf_v to GPU"
    endif
    istat = 0
    call alloc_var(fct_adf_h_gpu, adv_flux_hor, myDim_edge2D * (nl - 1), istat)
    if (istat /= 0) then
        write(0, *) "Error in alloc adv_flux_hor to GPU"
    endif
    istat = 0
    call alloc_var(UV_rhs_gpu, UV_rhs, 2 * myDim_elem2D * (nl - 1), istat)
    if (istat /= 0) then
        write(0, *) "Error in alloc UV_rhs to GPU"
    endif
    istat = 0
    call alloc_var(fct_ttf_max_gpu, fct_ttf_max, my_size * (nl - 1), istat)
    if (istat /= 0) then
        write(0, *) "Error in alloc fct_ttf_max to GPU"
    endif
    istat = 0
    call alloc_var(fct_ttf_min_gpu, fct_ttf_min, my_size * (nl - 1), istat)
    if (istat /= 0) then
        write(0, *) "Error in alloc fct_ttf_min to GPU"
    endif
    istat = 0
    call alloc_var(fct_plus_gpu, fct_plus, my_size * (nl - 1), istat)
    if (istat /= 0) then
        write(0, *) "Error in alloc fct_plus to GPU"
    endif
    istat = 0
    call alloc_var(fct_minus_gpu, fct_minus, my_size * (nl - 1), istat)
    if (istat /= 0) then
        write(0, *) "Error in alloc fct_minus to GPU"
    endif
#endif    
    if (mype==0) write(*,*) 'FCT is initialized'
end subroutine oce_adv_tra_fct_init
!===============================================================================
subroutine oce_tra_adv_fct(dttf_h, dttf_v, ttf, lo, adf_h, adf_v, mesh)
    !
    ! 3D Flux Corrected Transport scheme
    ! Limits antidiffusive fluxes==the difference in flux HO-LO
    ! LO ==Low-order  (first-order upwind)
    ! HO ==High-order (3rd/4th order gradient reconstruction method)
    ! Adds limited fluxes to the LO solution   
    use MOD_MESH
#ifdef FESOMCUDA
    use MOD_TRA_FCT_GPU
#endif
    use O_MESH
    use o_ARRAYS
    use o_PARAM
    use g_PARSUP
    use g_CONFIG
    use g_comm_auto
    implicit none
    type(t_mesh),  intent(in), target :: mesh
    real(kind=WP), intent(inout)      :: dttf_h(mesh%nl-1, myDim_nod2D+eDim_nod2D)
    real(kind=WP), intent(inout)      :: dttf_v(mesh%nl-1, myDim_nod2D+eDim_nod2D)
    real(kind=WP), intent(in)         :: ttf(mesh%nl-1, myDim_nod2D+eDim_nod2D)
    real(kind=WP), intent(in)         :: lo (mesh%nl-1, myDim_nod2D+eDim_nod2D)
    real(kind=WP), intent(inout)      :: adf_h(mesh%nl-1, myDim_edge2D)
    real(kind=WP), intent(inout)      :: adf_v(mesh%nl,  myDim_nod2D)
    integer                           :: n, nz, k, elem, enodes(3), num, el(2), nl1, nl2, edge
    real(kind=WP)                     :: flux, ae,tvert_max(mesh%nl-1),tvert_min(mesh%nl-1) 
    real(kind=WP)                     :: flux_eps=1e-16
    real(kind=WP)                     :: bignumber=1e3
    integer                           :: vlimit=1
    integer                           :: alg_state !state of the algorithm, useful when combining C/cuda with Fortran

#include "associate_mesh.h"
    alg_state = 0
#ifdef FESOMCUDA
    call fct_ale_pre_comm_acc(  alg_state, fct_ttf_max_gpu, fct_ttf_min_gpu, fct_plus_gpu, fct_minus_gpu,&
                                fct_ttf_gpu, fct_LO_gpu, fct_adf_v_gpu, fct_adf_h_gpu, UV_rhs_gpu, area_inv_gpu,& 
                                myDim_nod2D, eDim_nod2D, myDim_elem2D, myDim_edge2D, mesh%nl, nlevs_nod2D_gpu,& 
                                nlevs_elem2D_gpu, elem2D_nodes_gpu, nod_num_elem2D_gpu, nod_elem2D_gpu,&
                                size(nod_in_elem2D, 1), nod2D_edges_gpu, elem2D_edges_gpu, vlimit, flux_eps,&
                                bignumber, dt)
#endif
    if (alg_state < 1) then
        ! --------------------------------------------------------------------------
        ! ttf is the tracer field on step n
        ! del_ttf is the increment 
        ! vlimit sets the version of limiting, see below
        ! --------------------------------------------------------------------------
        !___________________________________________________________________________
        ! a1. max, min between old solution and updated low-order solution per node
        do n=1,myDim_nod2D + edim_nod2d
            do nz=1, nlevels_nod2D(n)-1 
                fct_ttf_max(nz,n)=max(LO(nz,n), ttf(nz,n))
                fct_ttf_min(nz,n)=min(LO(nz,n), ttf(nz,n))
            end do
        end do
    end if
    if (alg_state < 2) then
        !___________________________________________________________________________
        ! a2. Admissible increments on elements
        !     (only layers below the first and above the last layer)
        !     look for max, min bounds for each element --> UV_rhs here auxilary array
        do elem=1, myDim_elem2D
            enodes=elem2D_nodes(:,elem)
            do nz=1, nlevels(elem)-1
                UV_rhs(1,nz,elem)=maxval(fct_ttf_max(nz,enodes))
                UV_rhs(2,nz,elem)=minval(fct_ttf_min(nz,enodes))
            end do
            if (nlevels(elem)<=nl-1) then
                do nz=nlevels(elem),nl-1
                    UV_rhs(1,nz,elem)=-bignumber
                    UV_rhs(2,nz,elem)= bignumber
                end do
            endif
        end do ! --> do elem=1, myDim_elem2D
    end if    
    if (alg_state < 3) then
        !___________________________________________________________________________
        ! a3. Bounds on clusters and admissible increments
        ! Vertical1: In this version we look at the bounds on the clusters
        !            above and below, which leaves wide bounds because typically 
        !            vertical gradients are larger.  
        if(vlimit==1) then
            !Horizontal
            do n=1, myDim_nod2D
                !___________________________________________________________________
                do nz=1,nlevels_nod2D(n)-1
                    ! max,min horizontal bound in cluster around node n in every 
                    ! vertical layer
                    ! nod_in_elem2D     --> elem indices of which node n is surrounded
                    ! nod_in_elem2D_num --> max number of surrounded elem 
                    tvert_max(nz)= maxval(UV_rhs(1,nz,nod_in_elem2D(1:nod_in_elem2D_num(n),n)))
                    tvert_min(nz)= minval(UV_rhs(2,nz,nod_in_elem2D(1:nod_in_elem2D_num(n),n)))
                end do
                
                !___________________________________________________________________
                ! calc max,min increment of surface layer with respect to low order 
                ! solution 
                fct_ttf_max(1,n)=tvert_max(1)-LO(1,n)
                fct_ttf_min(1,n)=tvert_min(1)-LO(1,n)
                
                ! calc max,min increment from nz-1:nz+1 with respect to low order 
                ! solution at layer nz
                do nz=2,nlevels_nod2D(n)-2  
                    fct_ttf_max(nz,n)=maxval(tvert_max(nz-1:nz+1))-LO(nz,n)
                    fct_ttf_min(nz,n)=minval(tvert_min(nz-1:nz+1))-LO(nz,n)
                end do
                ! calc max,min increment of bottom layer -1 with respect to low order 
                ! solution 
                nz=nlevels_nod2D(n)-1
                fct_ttf_max(nz,n)=tvert_max(nz)-LO(nz,n)
                fct_ttf_min(nz,n)=tvert_min(nz)-LO(nz,n)  
            end do
        end if

        !___________________________________________________________________________
        ! Vertical2: Similar to the version above, but the vertical bounds are more 
        ! local  
        if(vlimit==2) then
            do n=1, myDim_nod2D
                do nz=1,nlevels_nod2D(n)-1
                    tvert_max(nz)= maxval(UV_rhs(1,nz,nod_in_elem2D(1:nod_in_elem2D_num(n),n)))
                    tvert_min(nz)= minval(UV_rhs(2,nz,nod_in_elem2D(1:nod_in_elem2D_num(n),n)))
                end do
                do nz=2, nlevels_nod2D(n)-2
                    tvert_max(nz)=max(tvert_max(nz),maxval(fct_ttf_max(nz-1:nz+1,n)))
                    tvert_min(nz)=min(tvert_min(nz),minval(fct_ttf_max(nz-1:nz+1,n)))
                end do
                do nz=1,nlevels_nod2D(n)-1
                    fct_ttf_max(nz,n)=tvert_max(nz)-LO(nz,n)
                    fct_ttf_min(nz,n)=tvert_min(nz)-LO(nz,n)  
                end do
            end do
        end if
        
        !___________________________________________________________________________
        ! Vertical3: Vertical bounds are taken into account only if they are narrower than the
        !            horizontal ones  
        if(vlimit==3) then
            do n=1, myDim_nod2D
                do nz=1,nlevels_nod2D(n)-1
                    tvert_max(nz)= maxval(UV_rhs(1,nz,nod_in_elem2D(1:nod_in_elem2D_num(n),n)))
                    tvert_min(nz)= minval(UV_rhs(2,nz,nod_in_elem2D(1:nod_in_elem2D_num(n),n)))
                end do
                do nz=2, nlevels_nod2D(n)-2
                    tvert_max(nz)=min(tvert_max(nz),maxval(fct_ttf_max(nz-1:nz+1,n)))
                    tvert_min(nz)=max(tvert_min(nz),minval(fct_ttf_max(nz-1:nz+1,n)))
                end do
                do nz=1,nlevels_nod2D(n)-1
                    fct_ttf_max(nz,n)=tvert_max(nz)-LO(nz,n)
                    fct_ttf_min(nz,n)=tvert_min(nz)-LO(nz,n)  
                end do
            end do
        end if
    end if
    if (alg_state < 4) then 
        !___________________________________________________________________________
        ! b1. Split positive and negative antidiffusive contributions
        ! --> sum all positive (fct_plus), negative (fct_minus) antidiffusive 
        !     horizontal element and vertical node contribution to node n and layer nz
        !     see. R. Löhner et al. "finite element flux corrected transport (FEM-FCT)
        !     for the euler and navier stoke equation
        do n=1, myDim_nod2D
            do nz=1,nlevels_nod2D(n)-1
                fct_plus(nz,n)=0._WP
                fct_minus(nz,n)=0._WP
            end do
        end do
        
        !Vertical
        do n=1, myDim_nod2D
            do nz=1,nlevels_nod2D(n)-1
    !             fct_plus(nz,n)=fct_plus(nz,n)+ &
    !                             (max(0.0_WP,adf_v(nz,n))+max(0.0_WP,-adf_v(nz+1,n))) &
    !                             /hnode(nz,n)
    !             fct_minus(nz,n)=fct_minus(nz,n)+ &
    !                             (min(0.0_WP,adf_v(nz,n))+min(0.0_WP,-adf_v(nz+1,n))) &
    !                             /hnode(nz,n)
                fct_plus(nz,n) =fct_plus(nz,n) +(max(0.0_WP,adf_v(nz,n))+max(0.0_WP,-adf_v(nz+1,n)))
                fct_minus(nz,n)=fct_minus(nz,n)+(min(0.0_WP,adf_v(nz,n))+min(0.0_WP,-adf_v(nz+1,n)))
            end do
        end do
    end if

    if (alg_state < 5) then
        !Horizontal
        do edge=1, myDim_edge2D
            enodes(1:2)=edges(:,edge)   
            el=edge_tri(:,edge)
            nl1=nlevels(el(1))-1
            nl2=0
            if(el(2)>0) then
                nl2=nlevels(el(2))-1
            end if   
            do nz=1, max(nl1,nl2)
                fct_plus (nz,enodes(1))=fct_plus (nz,enodes(1)) + max(0.0_WP, adf_h(nz,edge))
                fct_minus(nz,enodes(1))=fct_minus(nz,enodes(1)) + min(0.0_WP, adf_h(nz,edge))  
                fct_plus (nz,enodes(2))=fct_plus (nz,enodes(2)) + max(0.0_WP,-adf_h(nz,edge))
                fct_minus(nz,enodes(2))=fct_minus(nz,enodes(2)) + min(0.0_WP,-adf_h(nz,edge)) 
            end do
        end do 
    end if

    if (alg_state < 6) then
        !___________________________________________________________________________
        ! b2. Limiting factors
        do n=1,myDim_nod2D
            do nz=1,nlevels_nod2D(n)-1
                flux=fct_plus(nz,n)*dt/area(nz,n)+flux_eps
                fct_plus(nz,n)=min(1.0_WP,fct_ttf_max(nz,n)/flux)
                flux=fct_minus(nz,n)*dt/area(nz,n)-flux_eps
                fct_minus(nz,n)=min(1.0_WP,fct_ttf_min(nz,n)/flux)
            end do
        end do
    end if
    
    ! fct_minus and fct_plus must be known to neighbouring PE
    call exchange_nod(fct_plus, fct_minus)
    
    !___________________________________________________________________________
    ! b3. Limiting   
    !Vertical
    do n=1, myDim_nod2D
        nz=1
        ae=1.0_WP
        flux=adf_v(nz,n)
        if(flux>=0.0_WP) then 
            ae=min(ae,fct_plus(nz,n))
        else
            ae=min(ae,fct_minus(nz,n))
        end if
        adf_v(nz,n)=ae*adf_v(nz,n) 
        
        do nz=2,nlevels_nod2D(n)-1
            ae=1.0_WP
            flux=adf_v(nz,n)
            if(flux>=0._WP) then 
                ae=min(ae,fct_minus(nz-1,n))
                ae=min(ae,fct_plus(nz,n))
            else
                ae=min(ae,fct_plus(nz-1,n))
                ae=min(ae,fct_minus(nz,n))
            end if            
            adf_v(nz,n)=ae*adf_v(nz,n)
        end do
    ! the bottom flux is always zero 
    end do

        call exchange_nod_end  ! fct_plus, fct_minus
    !Horizontal
    do edge=1, myDim_edge2D
        enodes(1:2)=edges(:,edge)
        el=edge_tri(:,edge)
        nl1=nlevels(el(1))-1
        nl2=0
        if(el(2)>0) then
            nl2=nlevels(el(2))-1
        end if  
        do nz=1, max(nl1,nl2)
            ae=1.0_WP
            flux=adf_h(nz,edge)
            
            if(flux>=0._WP) then
                ae=min(ae,fct_plus(nz,enodes(1)))
                ae=min(ae,fct_minus(nz,enodes(2)))
            else
                ae=min(ae,fct_minus(nz,enodes(1)))
                ae=min(ae,fct_plus(nz,enodes(2)))
            endif
            
            adf_h(nz,edge)=ae*adf_h(nz,edge)
        end do
    end do
end subroutine oce_tra_adv_fct
