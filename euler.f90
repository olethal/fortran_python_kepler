

   SUBROUTINE EULER(TRAJECTORY, X11, X12, X21, X22, H, N)
      REAL :: MSUN, G, MEARTH
      REAL, INTENT(IN) :: H, X11, X12, X21, X22
      INTEGER, INTENT(IN) :: N
      REAL, DIMENSION(2) :: X, V, XNEW, VNEW
      REAL, DIMENSION(4,N), INTENT(INOUT) :: TRAJECTORY
      PRINT*, "H=", H, "N=", N
      !f2py intent(in) H, N, X11, X12, X21, X22
      !f2py intent(in, out) TRAJECTORY
      MSUN = 1.98847E30
      MEARTH = 1!2.99E-6*MSUN
      G = 6.67430E-11
      X(1) = X11
      X(2) = X12
      V(1) = X21
      V(2) = X22


      !X(1) = 0.0
      !X(2) = 150E9
      !V(1) = 30E3
      !V(2) = 0.0
      DO I = 1, N
         XNEW = X+H*V
         VNEW = V-H*MSUN*G*X/(DOT_PRODUCT(X,X)**(3.0/2.0))
         X = XNEW
         V = VNEW
         TRAJECTORY(4, I) = 0.5*MEARTH*DOT_PRODUCT(V, V)-G*MSUN*MEARTH/(DOT_PRODUCT(X, X))**0.5
         TRAJECTORY(1, I) = MEARTH*(X(1)*V(2)-X(2)*V(1))
         TRAJECTORY(2, I) = X(1)
         TRAJECTORY(3, I) = X(2)


      ENDDO

   END SUBROUTINE EULER
















