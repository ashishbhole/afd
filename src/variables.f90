module variables

   character (len=30) :: model, first_derivative, second_derivative
   character (len=30) :: time_integration       
   integer :: time_integration_levels 
   integer :: total_nodes, node
   real, allocatable :: CM(:, :), DM(:, :)
   complex, parameter :: iota = (0.0, 1.0)

end module variables
