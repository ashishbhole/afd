subroutine properties_convection_t2
use variables
implicit none
integer :: Nx, Ny, i, j, l, fid = 10
complex :: sumC, Aj, G
real :: dkh, dbeta_dkh
character(len=30) :: filename
real, allocatable :: kh(:), Nc(:)
real, allocatable :: modG(:,:), beta(:,:), phase_err(:,:), disp_err(:,:)

Ny = 157; Nx = 100

allocate(kh(1:Ny), Nc(1:Nx))
allocate(modG(Nx, Ny), beta(Nx, Ny), phase_err(Nx, Ny), disp_err(Nx, Ny))

do j = 1, Ny
   kh(j) = 0.02*j
   do i = 1, Nx
      Nc(i) = 0.02*i
      sumC = (0.0, 0.0)
      do l = 1, total_nodes
         sumC = sumC + Nc(i)*CM(node,l)*exp(iota*kh(j)*(l-node))
      enddo
      Aj = sumC

      call amplification_factor_two(Aj, G)
      modG(i, j) = sqrt( real(G)*real(G) + aimag(G)*aimag(G) )
      beta(i, j) = atan(-aimag(G)/real(G) )

      ! c_N/c = beta/(kh*Nc) ! need to correct by adding phase
      phase_err(i, j) = beta(i, j)/(kh(j)*Nc(i))
   enddo
enddo

dkh = kh(2) - kh(1)
do j = 1, Ny
   kh(j) = 0.02*j
   do i = 1, Nx
      Nc(i) = 0.02*i
      if(j==1)then
         dbeta_dkh = (beta(i, j+1)-beta(i, j))/dkh              
      elseif(j==Ny)then
         dbeta_dkh = (beta(i, j)-beta(i, j-1))/dkh              
      else
         dbeta_dkh = 0.5*(beta(i, j+1)-beta(i, j-1))/dkh
      endif
      disp_err(i, j) = (dbeta_dkh - (beta(i,j)/kh(j)))/(Nc(i)*kh(j))
   enddo
enddo

write(filename,'(a, a, a, a)') trim(time_integration), &
     & '_', trim(first_derivative), '.vtk'
open(fid, file=trim(filename))
write(fid,'("# vtk DataFile Version 3.0")')
write(fid,'("convection properties")')
write(fid,'("ASCII")')
write(fid,'("DATASET RECTILINEAR_GRID")')
write(fid,'("FIELD FieldData 2")')
write(fid,'("TIME 1 1 double")')
write(fid,'(e15.8)') 0.0
write(fid,'("CYCLE 1 1 int")')
write(fid,'(i8)') 0

write(fid,'("DIMENSIONS",i8, i8, i8)') size(Nc), size(kh), 1

write(fid,'("X_COORDINATES",i8," float")') size(Nc)
do i = 1, size(Nc)
   write(fid,*) Nc(i)
enddo

write(fid,'("Y_COORDINATES",i8," float")') size(kh)
do i = 1, size(kh)
   write(fid,*) kh(i)
enddo

write(fid,'("Z_COORDINATES",i8," float")') 1
write(fid,*) 0.0

write(fid,'("POINT_DATA",i16)') size(kh)*size(Nc)

write(fid,'("SCALARS modG float 1")')
write(fid,'("LOOKUP_TABLE default")')
do j = 1, size(kh)
do i = 1, size(Nc)
   write(fid, *) modG(i, j)
enddo
enddo

write(fid,'("SCALARS diff_err float 1")')
write(fid,'("LOOKUP_TABLE default")')
do j = 1, size(kh)
do i = 1, size(Nc)
   write(fid, *) phase_err(i, j)
enddo
enddo

write(fid,'("SCALARS disp_err float 1")')
write(fid,'("LOOKUP_TABLE default")')
do j = 1, size(kh)
do i = 1, size(Nc)
   write(fid, *) disp_err(i, j)
enddo
enddo

close(fid)

deallocate(kh, Nc, modG, beta, phase_err, disp_err)

end subroutine properties_convection_t2
