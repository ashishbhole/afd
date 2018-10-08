subroutine second_derivative_matrix
use variables
integer:: i, N
real :: A(total_nodes,total_nodes), B(total_nodes,total_nodes)
real :: AI(total_nodes,total_nodes)
real :: alpha, ac, bc

allocate(DM(total_nodes, total_nodes))
A = 0.0 ; B = 0.0 ; DM = 0.0
alpha = 0.0 ; ac = 0.0 ; bc = 0.0

N = total_nodes

! For explicit methods, A is Identity matrix
do i = 1, N
   A(i,i) = 1.0
enddo

if(trim(second_derivative) == 'CD2')then
   ! first and last node     
   B(1, 1) =  1.0 ; B(1, 2)   = -2.0 ; B(1, 3)   = 1.0
   B(N, N) =  1.0 ; B(N, N-1) = -2.0 ; B(N, N-2) = 1.0
   !interior nodes
   do i = 2, N-1
      B(i, i-1) =  1.0
      B(i, i)   = -2.0
      B(i, i+1) =  1.0
   enddo
elseif(trim(second_derivative) == 'CD4')then
   ! first and last node     
   B(1, 1) =  1.0 ; B(1, 2)   = -2.0 ; B(1, 3)   = 1.0
   B(N, N-2) =  1.0 ; B(N, N-1) = -2.0 ; B(N, N) = 1.0
   ! second and second last node
   B(2, 1) =  1.0 ; B(2, 2) = -2.0 ; B(2, 3) =  1.0
   B(N-1, N-2) =  1.0 ; B(N-1, N-1) = -2.0 ; B(N-1, N) =  1.0
   !interior nodes
   do i = 3, N-2
      B(i, i-2) = -1.0/12.0
      B(i, i-1) =  4.0/3.0
      B(i, i)   = -5.0/2.0
      B(i, i+1) =  4.0/3.0
      B(i, i+2) = -1.0/12.0
   enddo        
elseif(trim(second_derivative) == 'LELE')then   
   alpha = 2.0/11.0
   ac = 12.0/11.0
   bc = 3.0/11.0
   ! LHS matrix
   ! interior node
   do i = 3, N-2
      A(i, i-1) = alpha
      A(i,i)    = 1.0
      A(i, i+1) = alpha
   enddo
   ! RHS matrix
   ! first and last node     
   B(1, 1) =  1.0 ; B(1, 2)   = -2.0 ; B(1, 3)   = 1.0
   B(N, N-2) =  1.0 ; B(N, N-1) = -2.0 ; B(N, N) = 1.0
   ! second and second last node
   B(2, 1) =  1.0 ; B(2, 2) = -2.0 ; B(2, 3) =  1.0
   B(N-1, N-2) =  1.0 ; B(N-1, N-1) = -2.0 ; B(N-1, N) =  1.0   
   ! interior nodes
   do i = 3, N-2
      B(i, i-2) = 0.25*bc
      B(i, i-1) = ac
      B(i, i)   = -0.5*bc-2.0*ac
      B(i, i+1) = ac
      B(i, i+2) = 0.25*bc
   enddo
else
   print*, 'Select a numercial method for second derivative'
   call abort
endif

call inverse(A, AI)
DM = matmul(AI, B)

end Subroutine second_derivative_matrix
