subroutine first_derivative_matrix
use variables
integer:: i, N
real :: A(total_nodes,total_nodes), B(total_nodes,total_nodes)
real :: AI(total_nodes,total_nodes)
real :: alpha, ac, bc

allocate(CM(total_nodes, total_nodes))
A = 0.0 ; B = 0.0 ; CM = 0.0
alpha = 0.0 ; ac = 0.0 ; bc = 0.0

N = total_nodes

! For explicit methods, A is Identity matrix
do i = 1, N
   A(i,i) = 1.0
enddo

if(trim(second_derivative) == 'CD2')then
   B(1, 1) = -1.0 ; B(1, 2)  = 1.0
   B(N, N-1) = -1.0 ; B(N, N) = 1.0
   do i = 2, N-1
      B(i, i-1) = -0.5
      B(i, i+1) =  0.5
   enddo
elseif(trim(second_derivative) == 'CD4')then
   ! first and last node     
   B(1, 1) = -1.0 ; B(1, 2)  = 1.0
   B(N, N-1) = -1.0 ; B(N, N) = 1.0   
   ! second and second last node
   B(2, 1) =  -0.5 ; B(2, 3)   = 0.5
   B(N-1, N-2) = 0.5 ; B(N-1, N) = 0.5
   !interior nodes
   do i = 3, N-2
      B(i, i-2) = 1.0/12.0
      B(i, i-1) = -2.0/3.0
      B(i, i+1) = 2.0/3.0
      B(i, i+2) = -1.0/12.0
   enddo        
elseif(trim(second_derivative) == 'CD6')then
   ! first and last node     
   B(1, 1) = -1.0 ; B(1, 2)  = 1.0
   B(N, N-1) = -1.0 ; B(N, N) = 1.0   
   ! second and second last node
   B(2, 1) =  -0.5 ; B(2, 3)   = 0.5
   B(N-1, N-2) = 0.5 ; B(N-1, N) = 0.5
   ! third and third last node
   B(3, 1) = 1.0/12.0 ; B(3, 2) = -2.0/3.0 ; B(3, 4) = 2.0/3.0 ; B(3, 5) = -1.0/12.0
   B(N-2, N-4) = 1.0/12.0 ; B(N-2, N-3) = -2.0/3.0 ; B(N-2, N-1) = 2.0/3.0 ; B(N-2, N) = -1.0/12.0
   !interior nodes
   do i = 4, N-3
      B(i, i-3) = -1.0/60.0
      B(i, i-2) = 3.0/20.0
      B(i, i-1) = -3.0/4.0
      B(i, i+1) = 3.0/4.0
      B(i, i+2) = -3.0/20.0
      B(i, i+3) = 1.0/60.0
   enddo        
elseif(trim(second_derivative) == 'LELE')then
   alpha = 1.0/4.0
   ac = 14.0/9.0
   bc = 1.0/9.0
   ! LHS matrix
   ! interior node
   do i = 3, N-2
      A(i, i-1) = alpha
      A(i,i)    = 1.0
      A(i, i+1) = alpha      
   enddo   
   ! RHS matrix
   ! first and last node     
   B(1, 1) = -1.0 ; B(1, 2)  = 1.0
   B(N, N-1) = -1.0 ; B(N, N) = 1.0
   ! second and second last node
   B(2, 1) =  -0.5 ; B(2, 3)   = 0.5
   B(N-1, N-2) = 0.5 ; B(N-1, N) = 0.5   
   ! interior nodes
   do i = 3, N-2
      B(i, i-2) = -0.25*bc
      B(i, i-1) = -0.5*ac
      B(i, i+1) = 0.5*ac
      B(i, i+2) = 0.25*bc
   enddo        
else
   print*, 'Select available difference formula for first derivative'
   print*, 'Or add the stencil for the required formula in the code'
   call abort
endif

call inverse(A, AI)
CM = matmul(AI, B)

End Subroutine first_derivative_matrix
