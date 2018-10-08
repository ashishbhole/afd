subroutine properties_diffusion_t2
use variables
implicit none
integer :: Nx, Ny, i, j, l, fid = 10
complex :: sumD, Aj, G
real :: dkh, dlogG_dkh
character(len=30) :: filename
real, allocatable :: kh(:), Pe(:)
real, allocatable :: modG(:,:), beta(:,:), diff_err(:,:), disp_err(:,:)

Ny = 157; Nx = 100

allocate(kh(1:Ny), Pe(1:Nx))
allocate(modG(Nx, Ny), beta(Nx, Ny), diff_err(Nx, Ny), disp_err(Nx, Ny))

do j = 1, Ny
   kh(j) = 0.02*j
   do i = 1, Nx
      Pe(i) = 0.02*i
      sumD = (0.0, 0.0)
      do l = 1, total_nodes
         sumD = sumD + Pe(i)*DM(node,l)*exp(iota*kh(j)*(l-node))
      enddo
      Aj = sumD

      call amplification_factor_two(Aj, G)
      modG(i, j) = sqrt( real(G)*real(G) + aimag(G)*aimag(G) )
      beta(i, j) = atan(-aimag(G)/real(G) )

      ! nu_n/nu = -log(mod(G))/(kh*kh*Pe)
      diff_err(i, j) = - log(modG(i, j)) / (kh(j)*kh(j)*Pe(i))
   enddo
enddo

dkh = kh(2) - kh(1)
do j = 1, Ny
   kh(j) = 0.02*j
   do i = 1, Nx
      Pe(i) = 0.02*i
      if(j==1)then
         dlogG_dkh = (modG(i, j+1)-modG(i, j))/dkh
      elseif(j==Ny)then
         dlogG_dkh = (modG(i, j)-modG(i, j-1))/dkh
      else
         dlogG_dkh = 0.5*(modG(i, j+1)-modG(i, j-1))/dkh
      endif
      disp_err(i, j) = -(dlogG_dkh - (2.0*modG(i,j)/kh(j)))/(Pe(i)*Pe(i)*kh(j))
   enddo
enddo

write(filename,'(a, a, a, a)') trim(time_integration), &
     & '_', trim(second_derivative), '.vtk'
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
write(fid,*) 0.0

write(fid,'("POINT_DATA",i16)') size(kh)*size(Pe)

write(fid,'("SCALARS G_by_Ge float 1")')
write(fid,'("LOOKUP_TABLE default")')
do j = 1, size(kh)
do i = 1, size(Pe)
   write(fid, *) modG(i, j)/exp(-(kh(j)*kh(j)*Pe(i)))
enddo
enddo

write(fid,'("SCALARS phase float 1")')
write(fid,'("LOOKUP_TABLE default")')
do j = 1, size(kh)
do i = 1, size(Pe)
   write(fid, *) beta(i, j)
enddo
enddo

write(fid,'("SCALARS diff_err float 1")')
write(fid,'("LOOKUP_TABLE default")')
do j = 1, size(kh)
do i = 1, size(Pe)
   write(fid, *) diff_err(i, j)
enddo
enddo

write(fid,'("SCALARS disp_err float 1")')
write(fid,'("LOOKUP_TABLE default")')
do j = 1, size(kh)
do i = 1, size(Pe)
   write(fid, *) disp_err(i, j)
enddo
enddo

close(fid)

end subroutine properties_diffusion_t2
