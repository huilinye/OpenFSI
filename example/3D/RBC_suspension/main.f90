program main
    implicit none
    integer:: i,j,k, num_rbc, num, n_atom, n_coe, n_order
    integer:: single_point, single_bond, single_angle, single_dihedral
    integer:: t_point, t_bond, t_angle, t_dihedral
    real, dimension(:), allocatable :: single_atom(:,:), single_bond_coef(:,:), single_angle_coef(:,:), single_dihedral_coef(:,:)
    real, dimension(:), allocatable :: atom_xyz(:,:), bond_coef(:,:), angle_coef(:,:), dihedral_coef(:,:)
    real, dimension(:), allocatable :: rnd(:)
    integer, dimension(:), allocatable :: single_bond_order(:,:), single_angle_order(:,:), single_dihedral_order(:,:)
    integer, dimension(:), allocatable :: bond(:,:), angle(:,:), dihedral(:,:)
    real(8) :: box(3,2), Dx, Dy, Dz, L_x, L_y, L_z, x_temp, y_temp, z_temp
    real(8) :: rotatex(3,3),rotatey(3,3),rotatez(3,3),rotate(3,3)
    integer :: x_tot, y_tot, z_tot
    character(len=10) :: single_rbc_name
    character(len = 20) :: ofilename
    character(len=10) :: temp
    real(8) :: t1,t2
    
    
    CALL CPU_TIME(t1)
    !============================================read data from single rbc file
    write(single_rbc_name,*) '1rbc.data'
    Open(10 , file=trim(adjustl(single_rbc_name)), status = 'old')
    read(10, *) 
    read(10, *)
    
    read(10, '(i5)') single_point
    read(10, '(i5)') single_bond
    read(10, '(i5)') single_angle
    read(10, '(i5)') single_dihedral
    
    !allocate the array size
    allocate(single_bond_coef(single_bond,8), single_angle_coef(single_angle,9), single_dihedral_coef(single_dihedral,3))
    allocate(single_atom(single_point,6),single_bond_order(single_bond,4), single_angle_order(single_angle,5), single_dihedral_order(single_dihedral,6))
    
    do i=1,18
        read(10, *)
    end do
    
    do i=1,single_bond
        read(10,*) single_bond_coef(i,1:8)
    end do
    
    do i=1,3
        read(10, *)
    end do
    
    do i=1,single_angle
        read(10,*) single_angle_coef(i,1:9)
    end do
    
    do i=1,3
        read(10, *)
    end do
    
    do i=1,single_dihedral
        read(10,*) single_dihedral_coef(i,1:3)
    end do
    
    do i=1,3
        read(10, *)
    end do
    
    do i=1,single_point
        read(10,*) single_atom(i,1:6)
    end do
    
    do i=1,3
        read(10, *)
    end do
    
    do i=1,single_bond
        read(10,*) single_bond_order(i,1:4)
    end do
    
    do i=1,3
        read(10, *)
    end do
    
    do i=1,single_angle
        read(10,*) single_angle_order(i,1:5)
    end do
    
    do i=1,3
        read(10, *)
    end do
    
    do i=1,single_dihedral
        read(10,*) single_dihedral_order(i,1:6)
    end do
!================================end read data    
    
!================================define the number of RBCs and distribute them in space
    open(11,file = 'parameters.dat',status='old') 
    !define the domain box size
    read(11,*) 
    read(11,*) Dx
    read(11,*) Dy
    read(11,*) Dz
    
    
    read(11,*)
    read(11,*)
    read(11,*) box(1,1), box(1,2)
    read(11,*) box(2,1), box(2,2)
    read(11,*) box(3,1), box(3,2)
    
    read(11,*)
    read(11,*)
    read(11,*) num_rbc !total number of RBCs
    
    !========================allocate the size of the array
    allocate(atom_xyz(num_rbc*single_point,6), bond_coef(num_rbc*single_bond,8), angle_coef(num_rbc*single_angle,9), dihedral_coef(num_rbc*single_dihedral,3))
    allocate(bond(num_rbc*single_bond,4), angle(num_rbc*single_angle,5), dihedral(num_rbc*single_dihedral,6))
    read(11,*)
    read(11,*)
    read(11,*) x_tot 
    read(11,*)
    read(11,*)
    read(11,*) y_tot 
    read(11,*)
    read(11,*)
    read(11,*) z_tot 
    
    close(11)
    L_x = Dx/x_tot
    L_y = Dy/y_tot
    L_z = Dz/z_tot
    
    !recenter the single RBC position
        x_temp = 0.d0
        y_temp = 0.d0
        z_temp = 0.d0
        
    do i=1,single_point
        x_temp = x_temp+single_atom(i,4)
        y_temp = y_temp+single_atom(i,5)
        z_temp = z_temp+single_atom(i,6)
    end do
    
    do i=1,single_point
    single_atom(i,4) = single_atom(i,4)-x_temp/single_point
    single_atom(i,5) = single_atom(i,5)-y_temp/single_point
    single_atom(i,6) = single_atom(i,6)-z_temp/single_point
    end do
    
    !========================= rotate RBC with random angle
    allocate(rnd(num_rbc))
    call random_seed()
    call random_number(rnd) 
    
    !========================= distribute the coordinates of all the atom
    
    do i=1,x_tot
        do j=1,y_tot
            do k=1,z_tot
                num = i+x_tot*(j-1)+x_tot*y_tot*(k-1)
                rotatex =  reshape((/1.0,0,0,0,cos(rnd(num)),sin(rnd(num)),0,-sin(rnd(num)),cos(rnd(num))/),(/3,3/))
                rotatey =  reshape((/cos(rnd(num)),0,-sin(rnd(num)),0,1.0,0,sin(rnd(num)),0,cos(rnd(num))/),(/3,3/))
                rotatez =  reshape((/cos(rnd(num)),sin(rnd(num)),0,-sin(rnd(num)),cos(rnd(num)),0,0,0,1.0/),(/3,3/))
                rotate = MATMUL(rotatex,rotatey)
                rotate = MATMUL(rotate,rotatez)
                do n_atom = 1,single_point
atom_xyz((num-1)*single_point+n_atom,1)=n_atom+(num-1)*single_point; !node number
atom_xyz((num-1)*single_point+n_atom,2)=num;                             !molecule type
atom_xyz((num-1)*single_point+n_atom,3)=1;                             !atom type
atom_xyz((num-1)*single_point+n_atom,4)=rotate(1,1)*single_atom(n_atom,4)+rotate(1,2)*single_atom(n_atom,5)+rotate(1,3)*single_atom(n_atom,6)+(i-1)*L_x+L_x/2.0+box(1,1)     !node position, 4->6: x,y,z
atom_xyz((num-1)*single_point+n_atom,5)=rotate(2,1)*single_atom(n_atom,4)+rotate(2,2)*single_atom(n_atom,5)+rotate(2,3)*single_atom(n_atom,6)+(j-1)*L_y+L_y/2.0+box(2,1)
atom_xyz((num-1)*single_point+n_atom,6)=rotate(3,1)*single_atom(n_atom,4)+rotate(3,2)*single_atom(n_atom,5)+rotate(3,3)*single_atom(n_atom,6)+(k-1)*L_z+L_z/2.0+box(3,1)
                end do
            end do
        end do
    end do
    
    !write bond information
    do i=1,x_tot
        do j=1,y_tot
            do k=1,z_tot
                num = i+x_tot*(j-1)+x_tot*y_tot*(k-1)
                
                do n_atom = 1,single_bond
                    
                    bond_coef((num-1)*single_bond+n_atom,1) = n_atom+(num-1)*single_bond
                    bond((num-1)*single_bond+n_atom,1) = n_atom+(num-1)*single_bond
                    bond((num-1)*single_bond+n_atom,2) = 1!n_atom+(num-1)*single_bond
                    
                    do n_coe = 2,8
                    bond_coef((num-1)*single_bond+n_atom,n_coe) = single_bond_coef(n_atom,n_coe)
                    end do
                    
                    do n_order = 3,4
                        bond((num-1)*single_bond+n_atom,n_order) = single_bond_order(n_atom,n_order) + (num-1)*single_point
                    end do
                end do
                
                end do
            end do
    end do

    !write angle information
    do i=1,x_tot
        do j=1,y_tot
            do k=1,z_tot
                num = i+x_tot*(j-1)+x_tot*y_tot*(k-1)
                
                do n_atom = 1,single_angle
                    
                    angle_coef((num-1)*single_angle+n_atom,1) = n_atom+(num-1)*single_angle
                    angle((num-1)*single_angle+n_atom,1) = n_atom+(num-1)*single_angle
                    angle((num-1)*single_angle+n_atom,2) = 1!n_atom+(num-1)*single_angle
                    
                    do n_coe = 2,9
                    angle_coef((num-1)*single_angle+n_atom,n_coe) = single_angle_coef(n_atom,n_coe)
                    end do
                    
                    do n_order = 3,5
                        angle((num-1)*single_angle+n_atom,n_order) = single_angle_order(n_atom,n_order) + (num-1)*single_point
                    end do
                end do
                
                end do
            end do
    end do
    
    !write dihedral information
    do i=1,x_tot
        do j=1,y_tot
            do k=1,z_tot
                num = i+x_tot*(j-1)+x_tot*y_tot*(k-1)
                
                do n_atom = 1,single_dihedral
                    
                    dihedral_coef((num-1)*single_dihedral+n_atom,1) = n_atom+(num-1)*single_dihedral
                    dihedral((num-1)*single_dihedral+n_atom,1) = n_atom+(num-1)*single_dihedral
                    dihedral((num-1)*single_dihedral+n_atom,2) = 1!n_atom+(num-1)*single_dihedral
                    
                    do n_coe = 2,3
                    dihedral_coef((num-1)*single_dihedral+n_atom,n_coe) = single_dihedral_coef(n_atom,n_coe)
                    end do
                    
                    do n_order = 3,6
                        dihedral((num-1)*single_dihedral+n_atom,n_order) = single_dihedral_order(n_atom,n_order) + (num-1)*single_point
                    end do
                end do
                
                end do
            end do
    end do
    
    !========================================= write data file
    write (ofilename, '(I,A8)') num_rbc, 'rbc.data'
    open(13,file = adjustl(trim(ofilename)),status='replace')
    write(13,*) 'Using fortran to create large system with huge numer of RBCs'
    write(13,*)
    write(13,*) size(atom_xyz,1), 'atoms'
    write(13,*) size(bond,1), 'bonds'
    write(13,*) size(angle,1), 'angles'
    write(13,*) size(dihedral,1), 'dihedrals'
    write(13,*) 0, 'impropers'
    write(13,*)
    write(13,*) 1, 'atom types'
    write(13,*) 1, 'bond types'
    write(13,*) 1, 'angle types'
    write(13,*) 1, 'dihedral types'
    write(13,*) 0, 'impropers'
    write(13,*) 
    write(13,*) box(1,1), box(1,2), 'xlo',' ', 'xhi'
    write(13,*) box(2,1), box(2,2), 'ylo', ' ','yhi'
    write(13,*) box(3,1), box(3,2), 'zlo', ' ','zhi'
    write(13,*)
    write(13,*) 'Masses'
    write(13,*)
    write(13,*) 1, 1.0
    write(13,*) 
    write(13,*) 'Bond Coeffs #wlc'
    write(13,*) 
    
    !do i=1,num_rbc*single_bond  !write bond information
    !write(13,'(I,7f20.10)') int(bond_coef(i,1)),bond_coef(i,2:8)
    !end do
    write(13,'(I,7f20.10)') 1,bond_coef(1,2:8)
    
    write(13,*)
    write(13,*) 'Angle Coeffs #rbc'
    write(13,*)
    
    !do i=1,num_rbc*single_angle  !write angle information
    !write(13,'(I,8f20.10)') int(angle_coef(i,1)),angle_coef(i,2:9)
    !end do
    write(13,'(I,8f20.10)') int(angle_coef(1,1)),angle_coef(1,2:9)
    
    write(13,*)
    write(13,*) 'Dihedral Coeffs #bend'
    write(13,*) 
    
    !do i=1,num_rbc*single_dihedral  !write dihedral information
    !write(13,'(I,2f20.10)') int(dihedral(i,1)),dihedral_coef(i,2:3)
    !end do
    write(13,'(I,2f20.10)') int(dihedral(1,1)),dihedral_coef(1,2:3)
    
    write(13,*) 
    write(13,*) 'Atoms'
    write(13,*)
    
    do i=1,num_rbc*single_point  !write atom information
    write(13,'(3I,3f20.10)') int(atom_xyz(i,1:3)),atom_xyz(i,4:6)
    end do
    
    write(13,*) 
    write(13,*) 'Bonds'
    write(13,*)
    
    do i=1,num_rbc*single_bond  
    write(13,'(4I)') bond(i,1:4)
    end do
    
    write(13,*) 
    write(13,*) 'Angles'
    write(13,*)
    
    do i=1,num_rbc*single_angle  
    write(13,'(5I)') angle(i,1:5)
    end do
    
    write(13,*) 
    write(13,*) 'Dihedrals'
    write(13,*)
    
    do i=1,num_rbc*single_dihedral  
    write(13,'(6I)') dihedral(i,1:6)
    end do
    
    
    close(13)
    
    
    CALL CPU_TIME(t2)
    
    write(*,*) t2-t1
    end program