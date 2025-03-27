   program iri_path_integral
!Programer: ABN & GGBA
	implicit none
	integer i,j,k,l,l2,nx,ny,nz,dx,dy,dz,nato,zato,i1,j1,itemp,jtemp,ktemp,ng,itest,itype,jtype,ktype,ans,atom,rotx,roty,rotz,npairs,t
	integer ii, itest2,dt
	doubleprecision,dimension(:,:,:),allocatable:: mat1,mat2,mat3,mat4,f1,f2,f3
    integer,dimension(:,:),allocatable::pair
    doubleprecision, dimension(:),allocatable::atoms_x,atoms_y,atoms_z,coord_x,coord_y,coord_z,values
	doubleprecision xtemp,ytemp,ztemp,result,dr,result_old,step_x,step_y,step_z
	doubleprecision p(0:8,1:3)
	doubleprecision x0,y0,z0,x,y,z,ntemp,azato,test,x1,x2,y1,y2,z1,z2,dist,x_0,y_0,z_0,da,step
	character*50 arq0,arq1,arq2,arq3,arq4,arq5,arq6,arq7,arq8,arq9,arq10,arq11,arq12
	character*3 ato_label
	character*100 comments
!!!CLEANING CACHE

    call execute_command_line('rm sign_grid.txt')
    call execute_command_line('rm molecule.xyz molecule.ini')
    call execute_command_line('rm path-interpolated.txt')
    call execute_command_line('rm line-and-sign.txt')
!!! #######################################################
!!! Reading sign(lamb2)*rho   file                        #
!!! #######################################################
    call execute_command_line('tail -n +16 func1.cub > sign.txt')
    print *, 'Please check the following information about the system:'
    open(unit=1,file='func1.cub',status='old')
    !open(unit=1,file='func2.cub',status='old') !To look at IRI file
    read(1,*) comments
    read(1,*) comments
    read(1,*) nato, x_0, y_0, z_0 !!Number of atoms and origin of coordinates
    read(1,*) nx, step_x !step is the size of the step in the grid
    read(1,*) ny, comments, step_y
    read(1,*) nz, comments, comments, step_z
    step_x = step_x * 0.529177
    step_y = step_y * 0.529177
    step_z = step_z * 0.529177
    x_0 = x_0 * 0.529177
    y_0 = y_0 * 0.529177
    z_0 = z_0 * 0.529177
    allocate(atoms_x(nato))
    allocate(atoms_y(nato))
    allocate(atoms_z(nato))
    open (unit=2,file='molecule.xyz',status='new')
    write(2,*) nato
    write(2,*)
    do i = 1, nato
        read(1,*) atom, comments, atoms_x(i), atoms_y(i), atoms_z(i)
        if (atom==6) then
            ato_label = 'C'
        elseif (atom==16) then
            ato_label = 'S'
        elseif (atom==7) then
            ato_label = 'N'
        elseif (atom==1) then
            ato_label = 'H'
        elseif (atom==8) then
            ato_label = 'O'
        elseif (atom==9) then
            ato_label = 'F'
        elseif (atom==15) then
            ato_label = 'P'
        elseif (atom==17) then
            ato_label = 'Cl'
        endif
        
        write(2,*) ato_label, atoms_x(i)*0.529177, atoms_y(i)*0.529177, atoms_z(i)*0.529177
    enddo
    atoms_x = atoms_x*0.529177
    atoms_y = atoms_y*0.529177
    atoms_z = atoms_z*0.529177
    close (2)
    close (1)
    print *, 'The system has', nato,'atoms',' coordinate system starts at', x_0, y_0, z_0
    print *, 'The grid has',nx, 'points in X,',ny,'points in Y, and',nz, 'points in Z, ','total of',nx*ny*nz,'points!'
    !print *, 'Is it correct? Yes=1, No=0'
    !read *, ans
    !if (ans /= 1) then
    !    goto 100
    !endif
!!! #######################################################
!!! Allocating the 3D array                               #
!!! #######################################################
    allocate(f1(nx, ny, nz))
    open(unit=2,file='sign.txt',status="old")
    do i = 1, nx
        do j = 1, ny
            read(2, *) (f1(i, j, k), k = 1, nz)
        end do
    end do
    close(2)
    allocate(coord_x(nx),coord_y(ny),coord_z(nz))
    do  i = 1, nx
        coord_x(i)=(x_0 + (i-1)*step_x)
    enddo
    do j = 1, ny
        coord_y(j)=(y_0 + (j-1)*step_y)
    enddo
    do k = 1, nz
        coord_z(k)=(z_0 + (k-1)*step_z)
    enddo
    !print *, 'The coordinate', coord_x(nx), coord_y(ny), coord_z(nz),'has a signlamb*rho of', f1(nx,ny,nz)

!!! #######################################################
!!! Writing coordinates and sign in one file (optional)   #
!!! #######################################################
    open (unit=3,file='sign_grid.txt',status='new')
    do i = 1, nx
        do j = 1, ny
            do k = 1, nz
                write(3,'(F10.5, 1X, F10.5, 1X, F10.5, 1X, E13.5)') coord_x(i), coord_y(j), coord_z(k), f1(i, j, k)
            enddo
        enddo
    enddo

!!! #######################################################
!!! Finding and writing the pairs for study               #
!!! #######################################################
    print *, 'Choose the pair of atoms you want to calculate the sign of the path:'
    rotx=0  !Molecule default rotation in X
    roty=0  !Molecule default rotation in Y
    rotz=0  !Molecule default rotation in Z
    200 call execute_command_line('rm jmol.script')
    open(unit=1,file='jmol.script',status='new')
    write(1,*) 'load molecule.xyz'
    write(1,*) 'zoom out'
    write(1,*) 'rotate x', rotx
    write(1,*) 'rotate y', roty
    write(1,*) 'rotate z', rotz
    write(1,*) 'select all'
    write(1,*) 'label'
    write(1,*) 'font label 30'
    write(1,*) 'write POVRAY "molecule"'
    close (1)
    call execute_command_line('jmol -s jmol.script')
    call execute_command_line('transp-bg.sh molecule.ini')
    call execute_command_line('povray molecule.ini')
    print *, 'were the atoms all visible? Yes=1, No=0'
    read *, ans
    if (ans /= 1) then
        print *, 'type the values for the rotation in X, Y, and Z axis, repectively'
        read *, rotx, roty, rotz
        goto 200
    endif
    !Reading the pairs chosen by user
    print *, 'How many pair do you wish to calculate?'
    read *, npairs
    allocate(pair(npairs,2),values(npairs))
    do i = 1, npairs
        print *, 'Type the atoms that make pair', i
        do j = 1, 2
            read *, pair(i,j)
        enddo
    enddo

!!! #######################################################
!!! Linking the coordinates to the pairs (not necessary)  #
!!! #######################################################
    !do i = 1, npairs
    !    do j = 1, 2
    !        k=pair(i,j)
    !        write (1,*) coord_x(k), coord_y(k), coord_z(k)
    !    enddo
    !enddo

!!! #######################################################
!!! Interpolating the line between the pairs              #
!!! #######################################################
    values=0.0
    dt=50000.        ! Line 'resolution', how many points in the line between two atoms
    do i = 1, npairs
    	write (arq0,'("saida_",i1,".txt")') i
	!print *, arq0       
	open(unit=10,file=arq0,status='unknown')
        j = pair(i,1)
        k = pair(i,2)
        result_old=-10

        do t = 0, dt                        
            
            dr = (1.0) / dt                            
            
            x = atoms_x(j) + t*dr* (atoms_x(k) - atoms_x(j))
            y = atoms_y(j) + t*dr* (atoms_y(k) - atoms_y(j))
            z = atoms_z(j) + t*dr* (atoms_z(k) - atoms_z(j))
            
            itemp = aint((x-x_0)/step_x)
            xtemp = 1 + ((x-x_0)/step_x)
            jtemp = aint((y-y_0)/step_y)
            ytemp = 1 + ((y-y_0)/step_y)
            ktemp = aint((z-z_0)/step_z)
            ztemp = 1 + ((z-z_0)/step_z)

            x1=x_0+(itemp*step_x)
            x2=x_0+((itemp+1)*step_x)
            y1=y_0+(jtemp*step_y)
            y2=y_0+((jtemp+1)*step_y)
            z1=z_0+(ktemp*step_z)
            z2=z_0+((ktemp+1)*step_z)

            p(0,1)=x
            p(0,2)=y
            p(0,3)=z

            p(1,1)=x1
            p(1,2)=y1
            p(1,3)=z1

            p(2,1)=x2
            p(2,2)=y1
            p(2,3)=z1

            p(3,1)=x1
            p(3,2)=y2
            p(3,3)=z1

            p(4,1)=x2
            p(4,2)=y2
            p(4,3)=z1

            p(5,1)=x1
            p(5,2)=y1
            p(5,3)=z2

            p(6,1)=x2
            p(6,2)=y1
            p(6,3)=z2

            p(7,1)=x1
            p(7,2)=y2
            p(7,3)=z2

            p(8,1)=x2
            p(8,2)=y2
            p(8,3)=z2

            test=10000
            do ii=1,8,1
            dist=sqrt((p(0,1)-p(ii,1))**2+(p(0,2)-p(ii,2))**2+(p(0,3)-p(ii,3))**2)   
            if (dist.lt.test) then          
                test=dist
                itest=ii
            end if   
            !print*, "Ponto ",i,": (",p(i,1)," ",p(i,2)," ",p(i,3)," ) -----> D(P0,Pi): ",dist
            end do

            call trilinear_interpolation (xtemp,ytemp,ztemp,f1,nx,ny,nz,result)
            write(10,*) x, y, z, result
                       
            !Integral of values bigger then zero to check the most repulsive path among the chosen paths
            da=((result_old+result)*(dr)*sqrt((atoms_x(k)-atoms_x(j))**2&
                &+(atoms_y(k)-atoms_y(j))**2+(atoms_z(k)-atoms_z(j))**2))/2
            if (da > 0) then
               ! values(i)=values(i)+abs(result-result_old)*(dr)*sqrt((atoms_x(k)-atoms_x(j))**2&
                !&+(atoms_y(k)-atoms_y(j))**2+(atoms_z(k)-atoms_z(j))**2)
                !print *, result, (dr**3), (atoms_x(k) - atoms_x(j)), (atoms_y(k) - atoms_y(j)), (atoms_z(k) - atoms_z(j))
            	values(i)=values(i)+da
            endif
            result_old=result
            !print*, "Ponto P0          : (",p(0,1)," ",p(0,2)," ",p(0,3)," )  -----> Value: ",result
            !print*, "Ponto mais proximo: (",p(itest,1)," ",p(itest,2)," ",p(itest,3)," )"          
        enddo
        print *, 'Linear Approximatation Approach'
        print *,'Pair', i,'atoms', pair(i,1),'and', pair(i,2),'Total repulsion:', values(i)
    	close (10)
    enddo
    close(1)
    close(2)
    call execute_command_line('rm -f sign.txt')
    100 stop
    end

subroutine trilinear_interpolation(x, y, z, values, nx, ny, nz, result)
  real(8), intent(in) :: x, y, z            ! Coordinates of the point to interpolate
  real(8), intent(in) :: values(nx, ny, nz) ! 3D array of values
  integer, intent(in) :: nx, ny, nz         ! Dimensions of the array
  real(8), intent(out) :: result            ! Interpolated value

  integer :: i, j, k
  real(8) :: xd, yd, zd
  real(8) :: c00, c01, c10, c11, c0, c1

  i = int(x)
  j = int(y)
  k = int(z)

  xd = x - real(i)
  yd = y - real(j)
  zd = z - real(k)

  c00 = values(i, j, k) * (1.0 - xd) + values(i + 1, j, k) * xd
  c01 = values(i, j + 1, k) * (1.0 - xd) + values(i + 1, j + 1, k) * xd
  c10 = values(i, j, k + 1) * (1.0 - xd) + values(i + 1, j, k + 1) * xd
  c11 = values(i, j + 1, k + 1) * (1.0 - xd) + values(i + 1, j + 1, k + 1) * xd

  c0 = c00 * (1.0 - yd) + c01 * yd
  c1 = c10 * (1.0 - yd) + c11 * yd

  result = c0 * (1.0 - zd) + c1 * zd

end subroutine trilinear_interpolation