module variables

   character (len=30) :: model, first_derivative, second_derivative
   character (len=30) :: time_integration       
   integer :: time_integration_levels 
   integer :: total_nodes, node
   real, allocatable :: CM(:, :), DM(:, :)
   complex, parameter :: iota = (0.0, 1.0)

   integer :: nn, mm ! no of points in the stencil of FD formula at LHS and RHS
   real, allocatable :: ac(:), bc(:) 
   real, allocatable :: ac1(:), ac2(:), ac3(:), bc1(:), bc2(:), bc3(:)
   real, allocatable :: acN(:), acN1(:), acN2(:), bcN(:), bcN1(:), bcN2(:)   
   real, allocatable :: ad(:), bd(:)
   real, allocatable :: ad1(:), ad2(:), ad3(:), bd1(:), bd2(:), bd3(:)
   real, allocatable :: adN(:), adN1(:), adN2(:), bdN(:), bdN1(:), bdN2(:)   
   integer :: Nkh, NPe, NNc
   real :: kh_max, Pe_max, Nc_max
   real, parameter :: pi = 4.0*atan(1.0)
end module variables
