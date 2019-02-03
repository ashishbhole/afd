module tools
contains

subroutine inverse(A, AI)
use variables
implicit none
real, intent(in)  :: A(total_nodes, total_nodes)
real, intent(out) :: AI(total_nodes, total_nodes)
real :: AA(total_nodes, 2*total_nodes)
integer:: I, J, K, L3, LI, LK, KP1, JJ, LJ, M2, M, N
real :: P1, T, AB, BIG

  N = total_nodes

   M = N + N
   M2 = N + 1
!      GENERATING THE AUGMENTED MATRIX AA
   DO I = 1, N            !ENTERING THE COEFFICIENT MATRIX A INTO  
   DO J = 1, N           !THE AUGMENTED MATRIX AA       
      AA(I, J) = A(I, J)
   END DO
   END DO

   DO  I = 1, N
   DO  J = M2, M
       AA(I, J) = 0.0
   END DO
   END DO

   DO I = 1, N       !ENTERING THE IDENTITY MATRIX I INTO
      J = I + N         !THE AUGMENTED MATRIX AA
      AA(I, J) = 1.0
   END DO
!      STARTING ROW TRANSFORMATION
   DO LJ = 1, N
        K = LJ
        IF(K .LT. N)THEN
           JJ = K
           BIG = ABS(AA(K, K))
           KP1 = K + 1
           DO I = KP1, N
              AB = ABS(AA(I, K))
              IF((BIG - AB) .LT. 0.0)THEN
                 BIG = AB               !PERFORMING PIVOTING AND
                 JJ = I                 !ROW TRANSFORMATION OPERATIONS
              ENDIF
          END DO
          IF ((JJ - K) .GT. 0.0)THEN
              DO J = K, M
                 T = AA(JJ, J)
                 AA(JJ, J) = AA(K, J)
                 AA(K, J) = T
            END DO
            ENDIF
         ENDIF
          P1 = AA(LJ, LJ)
          DO I = LJ, M
             AA(LJ, I) = AA(LJ, I)/P1
          END DO
          DO LK = 1, N
            T = AA(LK, LJ)
            DO LI = LJ, M
              IF((LK - LJ) .NE. 0)THEN
                AA(LK, LI) = AA(LK, LI) - AA(LJ, LI)*T
              ENDIF
      END DO
      END DO
   END DO

   DO I = 1, N
   DO J = M2, M
      L3 = J - N
      AI(I, L3) = AA(I, J) !GENERATING THE INVERTED MATRIX
   END DO
   END DO

End subroutine inverse

end module tools
