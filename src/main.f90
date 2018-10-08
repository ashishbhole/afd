program main
use variables
implicit none
integer :: i
character(len=30) :: arg, input_filename

namelist /input_list/ model, first_derivative, second_derivative, time_integration_levels, &
                      time_integration, total_nodes, node

if(iargc() .ne. 1) then
   print*, "Please type : ./afd_1d input.in"
   call abort
endif

do i = 1, iargc()
  call getarg(i, arg)
  if(i == 1) input_filename = trim(arg)
enddo

write(*,*) " reading input file"
open(1, file=trim(input_filename))
  read(1, input_list)
close(1)
write(*,*) " reading input file complete"

print*, 'Model selected                                       : ', trim(model)
print*, 'Selected finite difference method for 1st derivative : ', trim(first_derivative)
print*, 'Selected finite difference method for 2nd derivative : ', trim(second_derivative)
print*, 'Selected time integration method                     : ', trim(time_integration)
print*, 'Levels of time-integration method                    : ', time_integration_levels

if(trim(model) == 'convection')then
   call first_derivative_matrix
   if(time_integration_levels == 2)then
      call properties_convection_t2
   elseif(time_integration_levels == 3)then
      !call properties_convection_t3     
   else
      print*, 'Specify valid time integration level: 2 or 3'
      call abort           
   endif
elseif(trim(model) == 'diffusion')then
   call second_derivative_matrix
   if(time_integration_levels == 2)then
      call properties_diffusion_t2
   elseif(time_integration_levels == 3)then
      !call properties_diffusion_t3
   else
      print*, 'Specify valid time integration level: 2 or 3'
      call abort
   endif   
elseif(trim(model) == 'convection-diffusion')then
   call first_derivative_matrix
   call second_derivative_matrix
   if(time_integration_levels == 2)then
      call properties_convection_diffusion_t2
   elseif(time_integration_levels == 3)then
      !call properties_convection_diffusion_t3
   else
      print*, 'Specify valid time integration level: 2 or 3'
      call abort
   endif   
else
   print*, 'Specify valid PDE model: convection, diffusion or convection-diffusion'
   call abort   
endif        

end program
