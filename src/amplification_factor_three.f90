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
