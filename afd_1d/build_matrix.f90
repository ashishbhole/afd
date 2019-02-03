module build_matrix
use variables
use tools

contains 

subroutine convection_matrix()
use variables
integer:: i, N
real :: A(total_nodes,total_nodes), B(total_nodes,total_nodes)
real :: AI(total_nodes,total_nodes)

N = total_nodes
if(.not. allocated(CM)) allocate(CM(N,N))
A = 0.0 ; B = 0.0 ; CM = 0.0

! interior nodes
do i = mm+1, N-mm
   A(i,i-nn:i+nn) = ac(-nn:nn)
   B(i,i-mm:i+mm) = bc(-mm:mm)
enddo

if(mm==1)then
   A(1,:) = ac1(:)
   B(1,:) = bc1(:) 

   A(N,:) = acN(:)
   B(N,:) = bcN(:)   
elseif(mm==2)then
   A(1,:) = ac1(:)
   B(1,:) = bc1(:)
   A(2,:) = ac2(:)
   B(2,:) = bc2(:)

   A(N,:) = acN(:)
   B(N,:) = bcN(:)
   A(N-1,:) = acN1(:)
   B(N-1,:) = bcN1(:)   
elseif(mm==3)then
   A(1,:) = ac1(:)
   B(1,:) = bc1(:)
   A(2,:) = ac2(:)
   B(2,:) = bc2(:)
   A(3,:) = ac3(:)
   B(3,:) = bc3(:)

   A(N,:) = acN(:)
   B(N,:) = bcN(:)
   A(N-1,:) = acN1(:)
   B(N-1,:) = bcN1(:)
   A(N-2,:) = acN2(:)
   B(N-2,:) = bcN2(:)
else
   print*, "something wrong in boundary conditions"
   print*, "check value of mm. It is = ", mm
   call abort   
endif

call inverse(A, AI)
CM = matmul(AI, B)

End Subroutine convection_matrix

subroutine diffusion_matrix()
use variables
integer:: i, N
real :: A(total_nodes,total_nodes), B(total_nodes,total_nodes)
real :: AI(total_nodes,total_nodes)

N = total_nodes
if(.not. allocated(DM)) allocate(DM(N,N))
A = 0.0 ; B = 0.0 ; DM = 0.0

! interior nodes
do i = mm+1, N-mm
   A(i,i-nn:i+nn) = ad(-nn:nn)
   B(i,i-mm:i+mm) = bd(-mm:mm)
enddo

if(mm==1)then
   A(1,:) = ad1(:)
   B(1,:) = bd1(:) 

   A(N,:) = adN(:)
   B(N,:) = bdN(:)   
elseif(mm==2)then
   A(1,:) = ad1(:)
   B(1,:) = bd1(:)
   A(2,:) = ad2(:)
   B(2,:) = bd2(:)

   A(N,:) = adN(:)
   B(N,:) = bdN(:)
   A(N-1,:) = adN1(:)
   B(N-1,:) = bdN1(:)   
elseif(mm==3)then
   A(1,:) = ad1(:)
   B(1,:) = bd1(:)
   A(2,:) = ad2(:)
   B(2,:) = bd2(:)
   A(3,:) = ad3(:)
   B(3,:) = bd3(:)

   A(N,:) = adN(:)
   B(N,:) = bdN(:)
   A(N-1,:) = adN1(:)
   B(N-1,:) = bdN1(:)
   A(N-2,:) = adN2(:)
   B(N-2,:) = bdN2(:)

else
   print*, "something wrong in boundary conditions"
   print*, "check value of mm. It is = ", mm
   call abort   
endif

call inverse(A, AI)
DM = matmul(AI, B)
End Subroutine diffusion_matrix

end module build_matrix
