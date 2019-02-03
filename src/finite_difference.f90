module finite_difference
use variables
contains

subroutine first_derivative_stencil
use variables
integer:: N

N = total_nodes
allocate(ac1(1:N), ac2(1:N), ac3(1:N))
allocate(bc1(1:N), bc2(1:N), bc3(1:N))
allocate(acN(1:N), acN1(1:N), acN2(1:N))
allocate(bcN(1:N), bcN1(1:N), bcN2(1:N))
ac1 = 0.0; ac2 = 0.0; ac3 = 0.0
bc1 = 0.0; bc2 = 0.0; bc3 = 0.0
acN = 0.0; acN1 = 0.0; acN2 = 0.0
bcN = 0.0; bcN1 = 0.0; bcN2 = 0.0

if(trim(first_derivative) == 'CD2')then
   ! no of coefficients at lhs and rhs in the stencil     
   nn = 0; mm = 1
   allocate(ac(-nn:nn), bc(-mm:mm))
   ac  = 0.0; bc  = 0.0

   !interior nodes
   ac(0) = 1.0
   bc(-1:1) = (/-0.5, 0.0, 0.5/)

   ! left boundaries
   ac1(1) = 1.0
   bc1(1:2) = (/-1.0, 1.0/)

   ! right boundaries
   acN(N) = 1.0
   bcN(N-1:N) = (/-1.0, 1.0/)
elseif(trim(first_derivative) == 'CD4')then
   ! no of coefficients at lhs and rhs in the stencil    
   nn = 0; mm = 2
   allocate(ac(-nn:nn), bc(-mm:mm))
   ac  = 0.0; bc  = 0.0

   !interior nodes
   ac(0) = 1.0
   bc(-2:2) = (/1.0/12.0, -2.0/3.0, 0.0, 2.0/3.0, -1.0/12.0/)

   ! left boundaries
   ac1(1) = 1.0
   bc1(1:2) = (/-1.0, 1.0/)
   ac2(2) = 1.0
   bc2(1:3) = (/-0.5, 0.0, 0.5/)

   !right boundaris
   acN(N) = 1.0
   bcN(1:2) = (/-1.0, 1.0/)
   acN1(N-1) = 1.0
   bcN1(N-2:N) = (/-0.5, 0.0, 0.5/)
elseif(trim(first_derivative) == 'CD6')then
   ! no of coefficients at lhs and rhs in the stencil    
   nn = 0; mm = 3
   allocate(ac(-nn:nn), bc(-mm:mm))
   ac  = 0.0; bc  = 0.0

   !interior nodes
   ac(0) = 1.0
   bc(-3:3) = (/-1.0/60.0, 3.0/20.0, -3.0/4.0, 0.0, 3.0/4.0, -3.0/20.0, 1.0/60.0/)

   ! left boundaries   
   ac1(1) = 1.0
   bc1(1:2) = (/-1.0, 1.0/)
   ac2(2) = 1.0
   bc2(1:3) = (/-0.5, 0.0, 0.5/)
   ac3(3) = 1.0
   bc3(1:4) = (/1.0/12.0, -2.0/3.0, 2.0/3.0, -1.0/12.0/)

  ! right boundaries
   acN(N) = 1.0
   bcN(N-1:N) = (/-1.0, 1.0/)
   acN1(N-1) = 1.0
   bcN1(N-2:N) = (/-0.5, 0.0, 0.5/)
   acN2(N-2) = 1.0
   bcN2(N-3:N) = (/1.0/12.0, -2.0/3.0, 2.0/3.0, -1.0/12.0/)
elseif(trim(first_derivative) == 'LELE')then
   ! no of coefficients at lhs and rhs in the stencil    
   nn = 1; mm = 2
   allocate(ac(-nn:nn), bc(-mm:mm))
   ac  = 0.0; bc  = 0.0
   !interior nodes
   ac(-1:1) = (/1.0/3.0, 1.0, 1.0/3.0/)
   bc(-2:2) = (/-1.0/36.0, -14.0/18.0, 0.0, 14.0/18.0, 1.0/36.0/)

   ! left boundaries
   ac1(1) = 1.0
   bc1(1:2) = (/-1.0, 1.0/)
   ac2(2) = 1.0
   bc2(1:3) = (/-0.5, 0.0, 0.5/)

   !right boundaris
   acN(N) = 1.0
   bcN(N-1:N) = (/-1.0, 1.0/)
   acN1(N-1) = 1.0
   bcN1(N-2:N) = (/-0.5, 0.0, 0.5/)

else
   print*, 'Select available difference formula for first derivative'
   print*, 'Or add the stencil for the required formula in the code'
   call abort
endif

End Subroutine first_derivative_stencil

subroutine second_derivative_stencil
use variables
real :: ak, bk

N = total_nodes
allocate(ad1(1:N), ad2(1:N), ad3(1:N))
allocate(bd1(1:N), bd2(1:N), bd3(1:N))
allocate(adN(1:N), adN1(1:N), adN2(1:N))
allocate(bdN(1:N), bdN1(1:N), bdN2(1:N))
ad1 = 0.0; ad2 = 0.0; ad3 = 0.0
bd1 = 0.0; bd2 = 0.0; bd3 = 0.0
adN = 0.0; adNm1 = 0.0; adNm2 = 0.0
bdN = 0.0; bdNm1 = 0.0; bdNm2 = 0.0

if(trim(second_derivative) == 'CD2')then
   ! no of coefficients at lhs and rhs in the stencil
   nn = 0; mm = 1
   allocate(ad(-nn:nn), bd(-mm:mm))
   ad  = 0.0; bd  = 0.0

   !interior nodes
   ad(0) = 1.0
   bd(-1:1) = (/1.0, -2.0, 1.0/)

   ! left boundaries
   ad1(1) = 1.0
   bd1(1:3) = (/1.0, -2.0, 1.0/)

   ! right boundaries
   adN(N) = 1.0
   bdN(N-2:N) = (/1.0, -2.0, 1.0/)
elseif(trim(second_derivative) == 'CD4')then
   ! no of coefficients at lhs and rhs in the stencil
   nn = 0; mm = 2
   allocate(ad(-nn:nn), bd(-mm:mm))
   ad  = 0.0; bd  = 0.0

   !interior nodes
   ad(0) = 1.0
   bd(-2:2) = (/-1.0/12.0, 4.0/3.0, -5.0/2.0, 4.0/3.0, -1.0/12.0/)

   ! left boundaries
   ad1(1) = 1.0
   bd1(1:3) = (/1.0, -2.0, 1.0/)
   ad2(2) = 1.0
   bd2(1:3) = (/1.0, -2.0, 1.0/)

   ! right boundaries
   adN(N) = 1.0
   bdN(N-2:N) = (/1.0, -2.0, 1.0/)
   adN1(N-1) = 1.0
   bdN1(N-2:N) = (/1.0, -2.0, 1.0/)
elseif(trim(second_derivative) == 'LELE')then
   ! no of coefficients at lhs and rhs in the stencil
   nn = 1; mm = 2
   allocate(ad(-nn:nn), bd(-mm:mm))
   ad  = 0.0; bd  = 0.0
   ak = 12.0/11.0; bk = 3.0/11.0

   !interior nodes
   ad(-1:1) = (/2.0/11.0, 1.0, 2.0/11.0/)
   bd(-2:2) = (/0.25*bk, ak, -0.5*bk-2.0*ak, ak, 0.25*bk/)

   ! left boundaries
   ad1(1) = 1.0
   bd1(1:3) = (/1.0, -2.0, 1.0/)
   ad2(2) = 1.0
   bd2(1:3) = (/1.0, -2.0, 1.0/)

   ! right boundaries
   adN(N) = 1.0
   bdN(N-2:N) = (/1.0, -2.0, 1.0/)
   adN1(N-1) = 1.0
   bdN1(N-2:N) = (/1.0, -2.0, 1.0/)
else
   print*, 'Select a numercial method for second derivative'
   call abort
endif

end Subroutine second_derivative_stencil

end module finite_difference
