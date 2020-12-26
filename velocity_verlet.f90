

   SUBROUTINE VELOCITY_VERLET(TRAJECTORY, X11, X12, X21, X22, H, N)
      REAL :: MSUN, G
      REAL, INTENT(IN) :: H, X11, X12, X21, X22
      INTEGER, INTENT(IN) :: N
      REAL, DIMENSION(2) :: X, V, FXI, XP1, FXP1, VP1
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

      FXI = -MSUN*G*X/(DOT_PRODUCT(X,X)**(3.0/2.0))
      DO I = 1, N
         XP1 = X + H*V + H*H*FXI/2.0
         FXP1 = -MSUN*G*XP1/(DOT_PRODUCT(XP1,XP1)**(3.0/2.0))
         VP1 = V + H*(FXI+FXP1)/2.0
         FXI = FXP1
         X = XP1
         V = VP1
         TRAJECTORY(1, I) = H*I
         TRAJECTORY(2, I) = X(1)
         TRAJECTORY(3, I) = X(2)

      ENDDO

   END SUBROUTINE VELOCITY_VERLET
















