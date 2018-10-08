subroutine properties_convection_diffusion_t2
use variables
implicit none
integer :: Nx, Ny, Nz, i, j, k, l, m, fid = 10
complex :: sumC, sumD, Aj, G
real :: dkh, dbeta_dkh, dlogG_dkh
character(len=50) :: filename
real, allocatable :: kh(:), Nc(:), Pe(:)
real, allocatable :: modG(:,:,:), beta(:,:,:), phase_err(:,:,:)
real, allocatable :: diff_err(:,:,:), dispc_err(:,:,:), dispd_err(:,:,:)

Nz = 157 ; Nx = 50 ; Ny = 50

allocate(kh(1:Nz), Pe(1:Ny), Nc(1:Nx))
allocate(modG(Nx, Ny, Nz), beta(Nx, Ny, Nz), phase_err(Nx, Ny, Nz))
allocate(diff_err(Nx, Ny, Nz), dispc_err(Nx, Ny, Nz), dispd_err(Nx, Ny, Nz))

kh = 0.0 ; Pe = 0.0 ; Nc = 0.0
modG = 0.0 ; beta = 0.0 ; diff_err = 0.0 
phase_err = 0.0 ; dispc_err = 0.0 ; dispd_err = 0.0 

do k = 1, Nz
   kh(k) = 0.02*k
   do j = 1, Ny
      Pe(j) = 0.04*j
      sumD = (0.0, 0.0)
      do l = 1, total_nodes
         sumD = sumD + Pe(j)*DM(node,l)*exp(iota*kh(k)*(l-node))
      enddo
      do i = 1, Nx
         Nc(i) = 0.04*i
         sumC = (0.0, 0.0)
         do l = 1, total_nodes
            sumC = sumC + Nc(i)*CM(node,l)*exp(iota*kh(k)*(l-node))
         enddo
         Aj = sumC + sumD

         call amplification_factor_two(Aj, G)
         modG(i, j, k) = sqrt( real(G)*real(G) + aimag(G)*aimag(G) )
         beta(i, j, k) = atan(-aimag(G)/real(G) )

         ! c_N/c = beta/(kh*Nc) ! need to correct by adding phase
         phase_err(i, j, k) = beta(i, j, k)/(kh(k)*Nc(i))
         
         ! nu_n/nu = -log(mod(G))/(kh*kh*Pe)
         diff_err(i, j, k) = - log(modG(i, j, k)) / (kh(k)*kh(k)*Pe(j))
      enddo
   enddo
enddo

dkh = kh(2) - kh(1)
do k = 1, Nz
   kh(k) = 0.02*k
   do j = 1, Ny
      Pe(j) = 0.04*j
      do i = 1, Nx
         Nc(i) = 0.04*i
         if(i==1)then
            dbeta_dkh = (beta(i+1, j, k)-beta(i, j, k))/dkh
         elseif(i==Nx)then
            dbeta_dkh = (beta(i, j, k)-beta(i-1, j, k))/dkh
         else
            dbeta_dkh = 0.5*(beta(i+1, j, k)-beta(i-1, j, k))/dkh
         endif
         dispc_err(i, j, k) = (dbeta_dkh - (beta(i,j,k)/kh(k)))/(Nc(i)*kh(k))

         if(j==1)then
            dlogG_dkh = (modG(i, j+1, k)-modG(i, j, k))/dkh
         elseif(j==Ny)then
            dlogG_dkh = (modG(i, j, k)-modG(i, j-1, k))/dkh
         else
            dlogG_dkh = 0.5*(modG(i, j+1, k)-modG(i, j-1, k))/dkh
         endif
         dispd_err(i, j, k) = -(dlogG_dkh - (2.0*modG(i,j,k)/kh(k)))/(Pe(j)*Pe(j)*kh(k))
      enddo
   enddo
enddo

call system('rm -rf cd_properties')
call system('mkdir cd_properties')

do j = 1, Ny

write(filename,'(a, a, a, a, a, a, a, i3.3, a)') 'cd_properties/', trim(time_integration), '_',&
      & trim(first_derivative), '_', trim(second_derivative), '_Pe_', j,'.vtk'
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
write(fid,*) Pe(j)

write(fid,'("POINT_DATA",i16)') size(kh)*size(Nc)

write(fid,'("SCALARS G_by_Ge float 1")')
write(fid,'("LOOKUP_TABLE default")')
do k = 1, size(kh)
do i = 1, size(Nc)
   write(fid, *) modG(i, j, k)/exp(-(kh(k)*kh(k)*Pe(j)))
enddo
enddo

write(fid,'("SCALARS phase_err float 1")')
write(fid,'("LOOKUP_TABLE default")')
do k = 1, size(kh)
do i = 1, size(Nc)
   write(fid, *) phase_err(i, j, k)
enddo
enddo

write(fid,'("SCALARS diff_err float 1")')
write(fid,'("LOOKUP_TABLE default")')
do k = 1, size(kh)
do i = 1, size(Nc)
   write(fid, *) diff_err(i, j, k)
enddo
enddo

write(fid,'("SCALARS dispc_err float 1")')
write(fid,'("LOOKUP_TABLE default")')
do k = 1, size(kh)
do i = 1, size(Nc)
   write(fid, *) dispc_err(i, j, k)
enddo
enddo

write(fid,'("SCALARS dispd_err float 1")')
write(fid,'("LOOKUP_TABLE default")')
do k = 1, size(kh)
do i = 1, size(Nc)
   write(fid, *) dispd_err(i, j, k)
enddo
enddo

close(fid)

enddo   ! enddo over j for Pe

do m = 1, Nx

write(filename,'(a, a, a, a, a, a, a, i3.3, a)') 'cd_properties/', trim(time_integration), '_',&
      & trim(first_derivative), '_', trim(second_derivative), '_Nc_', m,'.vtk'
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

write(fid,'("DIMENSIONS",i8, i8, i8)') size(Pe), size(kh), 1

write(fid,'("X_COORDINATES",i8," float")') size(Pe)
do i = 1, size(Pe)
   write(fid,*) Pe(i)
enddo

write(fid,'("Y_COORDINATES",i8," float")') size(kh)
do i = 1, size(kh)
   write(fid,*) kh(i)
enddo

write(fid,'("Z_COORDINATES",i8," float")') 1
write(fid,*) Nc(m)

write(fid,'("POINT_DATA",i16)') size(kh)*size(Pe)

write(fid,'("SCALARS G_by_Ge float 1")')
write(fid,'("LOOKUP_TABLE default")')
do k = 1, size(kh)
do j = 1, size(Nc)
   write(fid, *) modG(m, j, k)/exp(-(kh(k)*kh(k)*Pe(j)))
enddo
enddo

write(fid,'("SCALARS phase_err float 1")')
write(fid,'("LOOKUP_TABLE default")')
do k = 1, size(kh)
do j = 1, size(Pe)
   write(fid, *) phase_err(m, j, k)
enddo
enddo

write(fid,'("SCALARS diff_err float 1")')
write(fid,'("LOOKUP_TABLE default")')
do k = 1, size(kh)
do j = 1, size(Pe)
   write(fid, *) diff_err(m, j, k)
enddo
enddo

write(fid,'("SCALARS dispc_err float 1")')
write(fid,'("LOOKUP_TABLE default")')
do k = 1, size(kh)
do j = 1, size(Pe)
   write(fid, *) dispc_err(m, j, k)
enddo
enddo

write(fid,'("SCALARS dispd_err float 1")')
write(fid,'("LOOKUP_TABLE default")')
do k = 1, size(kh)
do j = 1, size(Pe)
   write(fid, *) dispd_err(m, j, k)
enddo
enddo

close(fid)

enddo   ! enddo over m for Nc

deallocate(kh, Pe, Nc, modG, beta, diff_err, dispc_err, dispd_err)

end subroutine properties_convection_diffusion_t2
