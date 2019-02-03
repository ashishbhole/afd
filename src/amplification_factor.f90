module amplification_factor
use variables
contains
subroutine amplification_factor_two(Aj, G)
use variables
complex, intent(in)  :: Aj
complex, intent(out) :: G
real, parameter :: one_by_six = 1.0/6.0
real, parameter :: one_by_twenty_four = 1.0/24.0

if(trim(time_integration) == 'Explicit-Euler')then
   G = 1 - Aj
elseif(trim(time_integration) == 'Implicit_Euler')then
  G = 1.0/(1.0 - Aj)
elseif(trim(time_integration) == 'RK4')then
   G = 1 - Aj + 0.5*Aj*Aj + one_by_six*Aj*Aj*Aj - one_by_twenty_four*Aj*Aj*Aj*Aj
elseif(trim(time_integration) == 'Cranck_Nicholson')then
  G = (1.0 + 0.5*Aj)/(1.0 - 0.5*Aj)
else
   print*, 'Select a numercial method for time integration'
   call abort
endif
end Subroutine amplification_factor_two

subroutine amplification_factor_three(Aj, G1, G2)
use variables
complex, intent(in)  :: Aj
complex, intent(out) :: G1, G2

if(trim(time_integration) == 'Richardson')then
   G1 = Aj + csqrt(Aj*Aj + 1)
   G2 = Aj - csqrt(Aj*Aj + 1)
!elseif(trim(time_integration) == 'Gears')then
else
   print*, 'Select a numercial method for time integration'
   call abort
endif

end subroutine amplification_factor_three
  
end module amplification_factor
