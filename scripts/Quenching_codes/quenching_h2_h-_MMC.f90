module quenching_module

  Implicit None
  Double Precision, Parameter :: bohr=0.52917721092
  Double Precision, Parameter :: cmtoau=4.5563352527D-6 
  Double Precision, Parameter :: autocm=2.194746313D5
  Double Precision, Parameter :: autokcalmol=627.51
  Double Precision, Parameter :: autoK=315775.13
  Double Precision, Parameter :: melectron=1822.88839
  Double Precision, Parameter :: Hmass=1.00782503223*melectron
  Double Precision, Parameter :: Dmass=2.01410177812*melectron
  Double Precision, Parameter :: freq_cutoff=2*cmtoau 
  
  Double Precision :: max_H_displacement=0.1
  Double Precision :: max_mol_displacement=0.3
  Double Precision :: max_rot=0.2 ! radians  
  Character(len=2), Allocatable :: atom_type(:)
  Integer :: Dim,Natoms,Nparticles,NH2,Npairs,Npairs_H2
  
contains
  
  subroutine pot_energy(energy,x,Dim)
    Implicit None

    ! Variables used to compute the potential energy
    Integer, Intent(In) :: Dim
    Double Precision, Intent(In) :: x(Dim)
    Double Precision, Intent(Out) :: energy ! energy in atomic units
    Double Precision :: func
    
    energy=func(Natoms,x)
  end subroutine pot_energy
  
  subroutine gradient(grad,x,Dim)
    Implicit None
   
    ! Variables used to compute the gradient
    Integer, Intent(In) :: Dim
    Double Precision, Intent(In) :: x(Dim)
    Double Precision, Intent(Out) :: grad(Dim) ! atomic units
    
    call dfunc(Natoms,x,grad)
  end subroutine gradient
  
  function atom_mass(atom)
    Implicit None
    Double Precision :: atom_mass
    Character(len=2), Intent(In) :: atom
    
    if (atom=='H') then
       atom_mass=Hmass
    else if (atom=='D') then
       atom_mass=Dmass
    else
       write(*,*) 'atom ', atom, ' is not recognized'
       stop 
    endif
  end function atom_mass

  function distance(N,a,b)
    Implicit None
    Integer :: N
    Double Precision :: a(N),b(N),distance
    distance=sum((a-b)**2)
    distance=sqrt(distance)
  end function distance

  subroutine pair_distances(x,pair_dist)
    Implicit None
    Double Precision, Intent(In) :: x(Dim)
    Integer :: i,j,l
    Double Precision, Intent(Out) :: pair_dist(Npairs)
    
    ! Given a configuration of NH2 hydrogen molecules compute the H2-H^- and H2-H2 pair distances using the H2 centers of mass
    pair_dist=0d0
    l=0
    do i=1,NH2
       l=l+1
       pair_dist(l)=sum(((x(6*i-2:6*i)+x(6*i+1:6*i+3))/2-x(1:3))**2) ! H2-H^- pair distance
       do j=1,i-1
          l=l+1
          pair_dist(l)=sum(((x(6*i-2:6*i)+x(6*i+1:6*i+3))/2-(x(6*j-2:6*j)+x(6*j+1:6*j+3))/2)**2) ! H2-H2 pair distance
       enddo
    enddo
    pair_dist=dsqrt(pair_dist)
    if (Npairs>1) call hpsort(Npairs,pair_dist)
  end subroutine pair_distances

  subroutine cosines_angles(x,cosines)
    Implicit None
    Double Precision, Intent(In) :: x(Dim)
    Integer :: i,j,l
    Double Precision, Intent(Out) :: cosines(Npairs_H2)
    Double Precision :: vec(3),vecn(3,NH2)
    l=0
    ! Compute the distance between each pair of adjacent H atoms in the xyz file and obtain the cosines of the pairs of bond angles
    do i=1,NH2
       vec=x(6*i-2:6*i)-x(6*i+1:6*i+3)
       vecn(:,i)=vec/sqrt(dot_product(vec,vec))
       do j=1,i-1
          l=l+1
          cosines(l)=abs(dot_product(vecn(:,i),vecn(:,j)))
       enddo
    enddo
    if (Npairs_H2>1) call hpsort(Npairs_H2,cosines)
  end subroutine cosines_angles

  subroutine CM(Natoms,x)
    Implicit None
    Integer, Intent(In) :: Natoms
    Integer :: i
    Double Precision :: sum_mass,vec(3),x(Dim)
    
    vec=0d0
    sum_mass=0d0
    do i=1,Natoms
       vec(:)=vec(:)+x(3*i-2:3*i)*atom_mass(atom_type(i))
       sum_mass=sum_mass+atom_mass(atom_type(i))
    enddo
    vec=vec/sum_mass
    do i=1,Natoms
       x(3*i-2:3*i)=x(3*i-2:3*i)-vec(:)
    enddo
  end subroutine CM

  subroutine moment_of_inertia(x,MI_eigenvals)   
    Implicit None
    Integer :: i,IERR
    Double Precision, Intent(In) :: x(3,Natoms)
    Double Precision, Intent(Out) :: MI_eigenvals(3)
    Double Precision :: IT(3,3),dummy,FV1(3),FV2(3)

    ! Compute the center of mass of the configuration
    call CM(Natoms,x)
    
    ! Compute the elements of the inertia tensor
    IT=0d0
    do i=1,Natoms
       ! Diagonal elements
       IT(1,1)=IT(1,1)+atom_mass(atom_type(i))*(x(2,i)**2+x(3,i)**2) ! I_xx
       IT(2,2)=IT(2,2)+atom_mass(atom_type(i))*(x(3,i)**2+x(1,i)**2) ! I_yy
       IT(3,3)=IT(3,3)+atom_mass(atom_type(i))*(x(1,i)**2+x(2,i)**2) ! I_zz
       
       ! Off-diagonal elements, Symmetric matrix
       IT(1,2)=IT(1,2)-atom_mass(atom_type(i))*x(1,i)*x(2,i) ! I_xy
       IT(2,3)=IT(2,3)-atom_mass(atom_type(i))*x(2,i)*x(3,i) ! I_yz
       IT(3,1)=IT(3,1)-atom_mass(atom_type(i))*x(3,i)*x(1,i) ! I_zx
    enddo
    
    IT(2,1)=IT(1,2)
    IT(3,2)=IT(2,3)
    IT(1,3)=IT(3,1)
    
    ! Diagonalize the moment of inertia tensor and find the eigenvalues
    MI_eigenvals=0d0
    call RS(3,3,IT,MI_eigenvals,0,dummy,FV1,FV2,IERR)
    MI_eigenvals=MI_eigenvals*bohr**2/melectron
  end subroutine moment_of_inertia

  ! Sorts arrays in ascending order
  subroutine hpsort(N,RA)
    Implicit None
    Integer, Intent(In) :: N
    Double Precision, Intent(Inout) :: RA(N)
    Integer :: I,IR,J,L
    Double Precision :: RRA
    L=N/2+1
    IR=N
    !The index L will be decremented from its initial value during the
    !"hiring" (heap creation) phase. Once it reaches 1, the index IR 
    !will be decremented from its initial value down to 1 during the
    !"retirement-and-promotion" (heap selection) phase.
10  continue
    if (L > 1) then
       L=L-1
       RRA=RA(L)
    else
       RRA=RA(IR)
       RA(IR)=RA(1)
       IR=IR-1
       if (IR.eq.1) then
          RA(1)=RRA
          return
       endif
    endif
    I=L
    J=L+L
20  if (J.le.IR) then
       if (J < IR) then
          if (RA(J) < RA(J+1)) J=J+1
       endif
       if (RRA < RA(J)) then
          RA(I)=RA(J)
          I=J; J=J+J
       else
          J=IR+1
       endif
       goto 20
    endif
    RA(I)=RRA
    goto 10
  end subroutine hpsort

  ! Computes the average energy at temperature Temp
  ! Three Metropolis Monte Carlo moves for all particles:
  !   1) Attempt to displace a randomly chosen hydrogen atom.
  !   2) Attempt to displace a randomly chosen hydrogen molecule.
  !   3) Attempt to rotate a randomly chosen hydrogen molecule.  
  subroutine Metropolis_MC(NMC,x,beta,Rc,energy0,average_energy)
    Implicit None
    Integer, Intent(In) :: NMC ! Number of MC moves
    Double Precision, Intent(In) :: beta,Rc ! beta in (atomic units)**(-1), Rc in bohr
    Double Precision, Intent(Inout) :: x(Dim),energy0 ! x in bohr, energy0 in atomic units
    Double Precision, Intent(Out) :: average_energy ! average_energy in atomic units
    Integer :: i,j,k,Nrestart,accept1,accept2,accept3,all1,all2,all3
    Double Precision :: x0(6),energy,x_CM(3),CM_H2(3)
    Double Precision :: r1,t,s(3),u ! harvest for random_number
    Double Precision :: Rot(3,3) ! rotation matrix

    average_energy=0d0
    Nrestart=min(1000,NMC)
    accept1=0
    accept2=0
    accept3=0
    all1=0
    all2=0
    all3=0

    ! Main loop for the number of Monte Carlo iterations
    do j=1,NMC
       call random_number(r1) ! used to decide what to do
       i=floor(r1*NH2)+1 ! do something with H2 molecule i (we have H^-(HH)(HH)(HH)...) in the xyz file where H^- is always first
       if (i>NH2) i=NH2
       call random_number(t) ! used to accept the trial move 
       call random_number(s) ! generate a random vector on [0, 1)**3 used for a trial displacement/rotation
       s=2*(s-0.5d0)
       
       call random_number(r1)
       x0(1:6)=x(6*i-2:6*i+3) ! hold the original coordinates of the H2 molecule
       if (r1<0.33) then ! Attempt to displace the hydrogen molecule 
          all2=all2+1
          x(6*i-2:6*i)=x(6*i-2:6*i)+max_mol_displacement*s
          x(6*i+1:6*i+3)=x(6*i+1:6*i+3)+max_mol_displacement*s
          if (sum(((x(6*i-2:6*i)+x(6*i+1:6*i+3))/2)**2)<Rc**2) then ! accept or reject the translation move of the H2 molecule
             call pot_energy(energy,x,Dim)
             if (exp(beta*(energy0-energy))>=t) then
                ! accept the new configuration
                energy0=energy
                accept2=accept2+1
                goto 111
             endif
          endif
          ! reject the new configuration
          x(6*i-2:6*i+3)=x0(1:6)
111       continue
          
       else if (r1>0.33.and.r1<0.66) then ! perform rotation of the hydrogen molecule around its center of mass
          all3=all3+1
          s=max_rot*s
          Rot=reshape((/ 1d0,0d0,0d0, 0d0,cos(s(1)),-sin(s(1)), 0d0,sin(s(1)),cos(s(1)) /), (/3, 3/))
          Rot=matmul(Rot,reshape((/ cos(s(2)),0d0,-sin(s(2)), 0d0,1d0,0d0, sin(s(2)),0d0,cos(s(2)) /), (/3, 3/)))
          Rot=matmul(Rot,reshape((/ cos(s(3)),-sin(s(3)),0d0, sin(s(3)),cos(s(3)),0d0, 0d0,0d0,1d0 /), (/3, 3/)))
          CM_H2(:)=(x(6*i-2:6*i)+x(6*i+1:6*i+3))/2 ! center of mass of the hydrogen molecule with the n-th
          x(6*i-2:6*i)=CM_H2(:)+matmul(Rot,x(6*i-2:6*i)-CM_H2(:)) ! rotate hydrogen 1 with respect to the center of mass
          x(6*i+1:6*i+3)=CM_H2(:)+matmul(Rot,x(6*i+1:6*i+3)-CM_H2(:)) ! rotate hydrogen 2 with respect to the center of mass
          call pot_energy(energy,x,Dim)
          if (exp(beta*(energy0-energy))>=t) then
             ! accept the new configuration
             energy0=energy
             accept3=accept3+1
             goto 112
          endif
          ! reject the new configuration
          x(6*i-2:6*i+3)=x0(1:6)
112       continue
          
       else if (r1>0.66) then ! attempt to symmetrically translate two adjacent hydrogen atoms
          all1=all1+1
          call random_number(u)
          s(:)=max_H_displacement*(u-0.5)*(x(6*i-2:6*i)-x(6*i+1:6*i+3))
          x(6*i-2:6*i)=x(6*i-2:6*i)+s
          x(6*i+1:6*i+3)=x(6*i+1:6*i+3)-s
          call pot_energy(energy,x,Dim)
          if (exp(beta*(energy0-energy))>=t) then
             ! accept the new configuration
             energy0=energy
             accept1=accept1+1
             goto 113
          endif
          ! reject the new configuration
          x(6*i-2:6*i+3)=x0(1:6)
113       continue
       endif
       ! Compute the running sum of the output energies
       average_energy=average_energy+energy0
       
       if (mod(j,Nrestart)==0) then
          ! Adjust the displacement parameters every 1000 moves for ~50% acceptance rates
          if (dble(accept1)/all1<0.5) then
             max_H_displacement=max_H_displacement*0.9
          else
             max_H_displacement=max_H_displacement*1.1
          endif
          
          if (dble(accept2)/all2<0.5) then
             max_mol_displacement=max_mol_displacement*0.9
          else
             max_mol_displacement=max_mol_displacement*1.1 
          endif
          
          if (dble(accept3)/all3<0.5) then
             max_rot=max_rot*0.9
          else
             max_rot=max_rot*1.1
          endif

          accept1=0
          accept2=0
          accept3=0
          all1=0
          all2=0
          all3=0
       endif
    enddo
    ! Compute the average energy
    average_energy=average_energy/NMC
  end subroutine Metropolis_MC

  subroutine harm_approx_gs(x,energy,omega,IERR)
    Implicit None
    Double Precision, Intent(In) :: x(Dim)
    Integer, Intent(Out) :: IERR
    Double Precision, Intent(Out) :: energy,omega(Dim)
    Integer :: i,j
    Double Precision :: r(Dim),H(Dim,Dim),FV1(Dim),FV2(Dim),dummy
    Double Precision :: grad0(Dim),grad(Dim),sqrt_mass(Dim)
    Double Precision, Parameter :: s=1d-6 
    
    ! Compute the Hessian
    r=x
    do i=1,Natoms
       sqrt_mass(3*i-2:3*i)=sqrt(atom_mass(atom_type(i)))
    enddo
    call gradient(grad0,r,Dim)       
    do i=1,Dim
       r(i)=x(i)+s
       call gradient(grad,r,Dim)
       r(i)=x(i)
       do j=1,Dim
          H(i,j)=(grad(j)-grad0(j))/s
       enddo
    enddo
    
    ! Symmetrize the Hessian
    do i=1,Dim
       do j=1,i
          if (i.ne.j) H(i,j)=(H(i,j)+H(j,i))/2
          H(i,j)=H(i,j)/(sqrt_mass(i)*sqrt_mass(j)) ! Mass-Scaled Hessian
          if (i.ne.j) H(j,i)=H(i,j)
       enddo
    enddo
    
    ! Diagonalize the mass-scaled Hessian
    call RS(Dim,Dim,H,omega,0,dummy,FV1,FV2,IERR) 
    do i=Dim,1,-1
       omega(i)=sign(dsqrt(dabs(omega(i))),omega(i)) ! Frequencies
    enddo
    call pot_energy(energy,x,Dim)
    IERR=0
    do i=1,Dim
       if (omega(i)>freq_cutoff) then 
          energy=energy+omega(i)/2
       else
          IERR=IERR+1
       endif
    enddo
  end subroutine harm_approx_gs
end module quenching_module

Program water_quench

  use quenching_module
  
  Implicit None
  Integer :: i,j,k,n,IERR_H2,IERR_D2
  Integer :: Nquench,NMC,Nisomers,Nisomers_max
  Double Precision :: energy0,average_energy,energy_opt,Temp,Rc ! All energies and Temp in atomic units
  Double Precision, Allocatable :: x(:),x_quench(:),pair_dist(:),cosines(:) ! All distances in bohr
  Double Precision, Allocatable :: lib_energy_opt(:),lib_pair_dist(:,:),lib_cosines(:,:),omega_H2(:),omega_D2(:)
  Double Precision :: grad_tol,r_energy_thresh,r_pairs_thresh,r_cos_thresh,energy_thresh,x_CM(3)
  Double Precision :: energy_HA_H2,energy_HA_D2,epsil1,epsil2,epsil3,x1,x2,x3,MI_eigenvals_H2(3),MI_eigenvals_D2(3)
  Integer :: grad_maxfac,info,iter,ngrad,nfunc,count_unphysical,incr_write
  Character(len=50) :: coord_config0,file_num,coord_isomer 
  
  open(1,file='input.dat')
  read(1,*) Natoms
  read(1,*) coord_config0
  read(1,*) Nisomers_max
  read(1,*) Temp,Rc ! Temp in Kelvin, Rc in angstroms
  read(1,*) grad_tol,grad_maxfac,energy_thresh ! energy_thresh in kcal/mol
  read(1,*) Nquench,NMC ! Number of times quenching is done and number of MC moves
  read(1,*) r_energy_thresh,r_pairs_thresh,r_cos_thresh
  read(1,*) incr_write
  close(1)

  open(3,file='isomers.dat')
  write(3,*) '# Column 1: Isomer Index'
  write(3,*) '# Column 2: Quenched Isomer Energy (au)'
  write(3,*) '# Column 3: Conjugate Gradient Info'
  write(3,*) '# Column 4: Number of Grad Evals'
  
  open(4,file='count_unphysical.dat') 
  write(4,*) '# Column 1: Quenching Index'
  write(4,*) '# Column 2: Current Unphysical Config Total'
  write(4,*) '# Column 3: Conjugate Gradient Info'
  write(4,*) '# Column 4: Number of Grad Evals'

  open(9,file='energy_unquenched.dat')
  write(9,*) '# Column 1: Quenching Index'
  write(9,*) '# Column 2: Unquenched Config Energy (au)'
  write(9,*) '# Column 3: Average Energy from Metropolis MC (au)'

  open(11,file='ha_energies.dat')
  write(11,*) '# Column 1: Isomer Index'
  write(11,*) '# Column 2: Minimum/Classical Energy (au)'
  write(11,*) '# Colunn 3: HA Energy (H2) (au)'
  write(11,*) '# Colunn 4: HA Energy (D2) (au)'

  open(12,file='moments_of_inertia_isomers_h2.dat')
  write(12,*) '# Column 1: Isomer Number'
  write(12,*) '# Columns 2-4: Moments of Inertia (amu*A**2)'
  
  open(13,file='moments_of_inertia_isomers_d2.dat')
  write(13,*) '# Column 1: Isomer Number'
  write(13,*) '# Columns 2-4: Moments of Inertia (amu*A**2)'

  open(14,file='frequencies_h2.dat')
  write(14,*) '# Column 1: Isomer Number'
  write(14,*) '# Column 2: Normal Mode Indicies'
  write(14,*) '# Column 3: Normal Mode Frequncies (cm^-1)'

  open(15,file='frequencies_d2.dat')
  write(15,*) '# Column 1: Isomer Number'
  write(15,*) '# Column 2: Normal Mode Indicies'
  write(15,*) '# Column 3: Normal Mode Frequncies (cm^-1)'
  
  open(8,file='epsil_thresh_values.dat')

  call flush(3)
  call flush(4)
  call flush(9)
  call flush(11)
  call flush(12)
  call flush(13)
  call flush(14)
  call flush(15)

  ! Set the relevant parameters and allocate the desired memory
  Dim=3*Natoms
  Nparticles=(Natoms-1)/2+1
  NH2=(Natoms-1)/2
  Npairs=(Nparticles*(Nparticles-1))/2
  Npairs_H2=(NH2*(NH2-1))/2
  Temp=Temp/autoK
  epsil1=100
  epsil2=100 
  epsil3=100
  allocate(x(Dim),x_quench(Dim),atom_type(Natoms),pair_dist(Npairs),cosines(Npairs_H2),omega_H2(Dim),omega_D2(Dim))
  allocate(lib_energy_opt(Nisomers_max),lib_pair_dist(Npairs,Nisomers_max),lib_cosines(Npairs_H2,Nisomers_max))

  ! Read in the initial optimized configuration
  open(2,file=coord_config0)
  read(2,*) 
  read(2,*)
  do i=1,Natoms
     read(2,*) atom_type(i),x(3*i-2:3*i)
  enddo
  close(2)
  
  ! Convert from angstroms to bohr as H2-H^- PEF takes configurations only in bohr
  x=x/bohr
  Rc=Rc/bohr

  ! Translate each atom in the cluster about the fixed H^- ion at the origin
  do k=2,Natoms
     x(3*k-2:3*k)=x(3*k-2:3*k)-x(1:3)
  enddo
  x(1:3)=0

  ! Initialize the necessary arrays
  lib_energy_opt=0d0
  lib_pair_dist=0d0
  lib_cosines=0d0
  Nisomers=0
  count_unphysical=0

  ! Determine the energy of the initial configuration
  call pot_energy(energy0,x,Dim)

  do n=1,Nquench
     ! Use Metropolis MC to displace or rotate the configuration for NMC steps
     call Metropolis_MC(NMC,x,1/Temp,Rc,energy0,average_energy)
     
     ! Write out the energy of the unquenched configs and the average energy from Metropolis MC
     if (mod(n,incr_write)==0) then
        write(9,*) n,energy0,average_energy
        call flush(9)
     endif

     ! Quench the configuration output from Metropolis MC
     x_quench=x
     call cg_descent(grad_tol,grad_maxfac*Dim,x_quench,Dim,pot_energy,gradient,&
          info,energy_opt,iter,nfunc,ngrad,energy_thresh)

     ! If the configuration is physical, compute its center of mass pair distances and bond angle cosines
     if (info==0.or.info==2) then
        ! Obtain the ground state energy of the quenched config from HA using the H mass
        call harm_approx_gs(x_quench,energy_HA_H2,omega_H2,IERR_H2) ! energy_HA_H2 in atomic units
        call moment_of_inertia(x_quench,MI_eigenvals_H2)
        
        ! Switch the atom type to 'D' and obtain the ground state energy from HA using the D mass
        atom_type='D'
        call harm_approx_gs(x_quench,energy_HA_D2,omega_D2,IERR_D2) ! energy_HA_D2 in atomic units
        call moment_of_inertia(x_quench,MI_eigenvals_D2)
        ! Switch the atom type back to 'H' 
        atom_type='H'
        
        if (IERR_H2==6.and.IERR_D2==6) then ! the configuration is a minimum and try to assign it to one of the configs already in the library
           call pair_distances(x_quench,pair_dist)
           call cosines_angles(x_quench,cosines)
          
           ! Loop over all the configurations stored in the library
           do k=1,Nisomers
              if (dabs(energy_opt-lib_energy_opt(k))<r_energy_thresh &
                   .and.distance(Npairs,pair_dist,lib_pair_dist(:,k))<r_pairs_thresh &
                   .and.distance(Npairs_H2,cosines,lib_cosines(:,k))<r_cos_thresh) goto 100 ! configuration is already present in the library           
           enddo
           
           ! If the configuration is not found, add it to the library
           Nisomers=Nisomers+1
           if (Nisomers>Nisomers_max) then
              write(10,*) 'Stop: Nisomers=',Nisomers,'>',Nisomers_max ! Nconfigs exceeds Nconfigs_max
              stop
           else
              ! Add the current config's parameters to the library
              lib_energy_opt(Nisomers)=energy_opt
              lib_pair_dist(:,Nisomers)=pair_dist
              lib_cosines(:,Nisomers)=cosines
              write(12,*) Nisomers,(MI_eigenvals_H2(j),j=1,3)
              write(13,*) Nisomers,(MI_eigenvals_D2(j),j=1,3)
              do j=1,Dim
                 if (omega_H2(j)>freq_cutoff) write(14,*) Nisomers,j-6,omega_H2(j)*autocm
                 if (omega_D2(j)>freq_cutoff) write(15,*) Nisomers,j-6,omega_D2(j)*autocm
              enddo
              write(14,*)
              write(15,*)
              call flush(12)
              call flush(13)
              call flush(14)
              call flush(15)
              
              ! Compute the minimum values for epsilon (thresholds for quenching/assigning isomers)
              do j=1,Nisomers-1
                 x1=dabs(lib_energy_opt(j)-lib_energy_opt(Nisomers))
                 x2=distance(Npairs,lib_pair_dist(:,j),lib_pair_dist(:,Nisomers))
                 x3=distance(Npairs_H2,lib_cosines(:,j),lib_cosines(:,Nisomers))
                 if (epsil1>x1) epsil1=x1
                 if (epsil2>x2) epsil2=x2
                 if (epsil3>x3) epsil3=x3
              enddo
              
              ! Write out the HA energies to a file
              write(11,*) Nisomers,energy_opt,energy_HA_H2,energy_HA_D2 ! HA ground states in atomic units
              call flush(11)
                 
              ! Write out the new config to an xyz file in the library and some of its important parameters to another file
              write(file_num,'(i7)') Nisomers
              coord_isomer='library_opt_configs/'//trim(adjustl(file_num))//'.xyz'
              open(7,file=coord_isomer)
              write(7,*) Natoms
              write(7,*) energy_opt,'|au'
              do i=1,Natoms
                 write(7,*) atom_type(i),x_quench(3*i-2:3*i)*bohr
              enddo
              close(7)
              write(3,*) Nisomers,energy_opt,info,ngrad
              call flush(3)
           endif
        endif
     else if (info==1.or.info.ge.3) then ! don't bother with the assignment
        count_unphysical=count_unphysical+1
        write(4,*) n,count_unphysical,info,ngrad
        call flush(4)
        write(file_num,'(i7)') count_unphysical
        coord_isomer='unphysical_configs/'//trim(adjustl(file_num))//'_quenched.xyz'
        open(7,file=coord_isomer)
        write(7,*) Natoms
        write(7,*) energy_opt,'|au'
        do i=1,Natoms
           write(7,*) atom_type(i),x_quench(3*i-2:3*i)*bohr
        enddo
        close(7)
        write(file_num,'(i7)') count_unphysical
        coord_isomer='unphysical_configs/'//trim(adjustl(file_num))//'_unquenched.xyz'
        open(16,file=coord_isomer)
        write(16,*) Natoms
        write(16,*) energy0,'|au'
        do i=1,Natoms
           write(16,*) atom_type(i),x(3*i-2:3*i)*bohr
        enddo
        close(16)
     endif
100  continue
  enddo
  
  ! Write out the epsilon parameters at the very end
  write(8,*) 'The minimum energy difference between different configurations (au):',epsil1
  write(8,*) 'The minimum distance between different configurations using the center of mass H2-H2 distances (bohr):',epsil2
  write(8,*) 'The minimum distance between different configurations using the bond angle cosines:',epsil3
  
  close(3)
  close(4)
  close(8)
  close(9)
  close(11)
  close(12)
  close(13)
  close(14)
  close(15)
end program water_quench
