

   SUBROUTINE EULER(TRAJECTORY, X11, X12, X21, X22, H, N)
      REAL :: MSUN, G
      REAL, INTENT(IN) :: H, X11, X12, X21, X22
      INTEGER, INTENT(IN) :: N
      REAL, DIMENSION(2) :: X, V
      REAL, DIMENSION(3,N), INTENT(INOUT) :: TRAJECTORY
      PRINT*, "H=", H, "N=", N
      !f2py intent(in) N, H, X11, X12, X21, X22
      !f2py intent(in, out) TRAJECTORY
      MSUN = 1.98847E30
      G = 6.67430E-11

      X(1) = X11
      X(2) = X12
      V(1) = X21
      V(2) = X22
      DO I = 1, N
         X = X+H*V
         V = V-H*MSUN*G*X/(DOT_PRODUCT(X,X)**(3.0/2.0))
         TRAJECTORY(1, I) = H*I
         TRAJECTORY(2, I) = X(1)
         TRAJECTORY(3, I) = X(2)

      ENDDO

   END SUBROUTINE EULER
















