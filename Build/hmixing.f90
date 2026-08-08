      MODULE hmixing_mod
!
!svn $Id: hmixing.F 830 2017-01-24 21:21:11Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2017 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine computes time-dependent 3D horizontal mixing           !
!  coefficients.                                                       !
!                                                                      !
!  References:                                                         !
!                                                                      !
!    Smagorinsky, J, 1963: General circulation experiments with        !
!      the primitive equations: I. The basic experiment, Mon.          !
!      Wea. Rev., 91, 99-164.                                          !
!                                                                      !
!    Holland, W.R., J.C. Chow, and F.O. Bryan, 1998: Application       !
!      of a Third-Order Upwind Scheme in the NCAR Ocean Model, J.      !
!      Climate, 11, 1487-1493.                                         !
!                                                                      !
!    Webb, D.J., B.A. De Cuevas, and C.S. Richmond, 1998: Improved     !
!      Advection Schemes for Ocean Models, J. Atmos. Oceanic           !
!      Technol., 15, 1171-1187.                                        !
!                                                                      !
!    Griffies, S.M. and R.W. Hallberg, 2000: Biharmonic Friction       !
!      with a Smagorinsky-like Viscosity for Use in Large-Scale        !
!      Eddy-Permitting Ocean Models, Monthly Weather Review, 128,      !
!      2935-2946.                                                      !
!                                                                      !
!    Marchesiello, P., L. Debreu, and Xavien Couvelard, 2008:          !
!      Spurious diapycnal mixing in terrain-following coordinate       !
!      models" advection problem and solutions, DRAFT.                 !
!                                                                      !
!  This routine was adapted from a routine provided by Patrick         !
!  Marchiesello (April 2008).                                          !
!                                                                      !
!=======================================================================
!
      implicit none
      PRIVATE
      PUBLIC  :: hmixing
      CONTAINS
!
!***********************************************************************
      SUBROUTINE hmixing (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
      USE mod_mixing
      USE mod_ocean
      USE mod_stepping
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
      integer :: IminS, ImaxS, JminS, JmaxS
      integer :: LBi, UBi, LBj, UBj, LBij, UBij
!
!  Set horizontal starting and ending indices for automatic private
!  storage arrays.
!
      IminS=BOUNDS(ng)%Istr(tile)-4
      ImaxS=BOUNDS(ng)%Iend(tile)+3
      JminS=BOUNDS(ng)%Jstr(tile)-4
      JmaxS=BOUNDS(ng)%Jend(tile)+3
!
!  Determine array lower and upper bounds in the I- and J-directions.
!
      LBi=BOUNDS(ng)%LBi(tile)
      UBi=BOUNDS(ng)%UBi(tile)
      LBj=BOUNDS(ng)%LBj(tile)
      UBj=BOUNDS(ng)%UBj(tile)
!
!  Set array lower and upper bounds for MIN(I,J) directions and
!  MAX(I,J) directions.
!
      LBij=BOUNDS(ng)%LBij
      UBij=BOUNDS(ng)%UBij
!
      CALL wclock_on (ng, iNLM, 28)
      CALL hmixing_tile (ng, tile,                                      &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   IminS, ImaxS, JminS, JmaxS,                    &
     &                   nrhs(ng),                                      &
     &                   GRID(ng) % rmask,                              &
     &                   GRID(ng) % umask,                              &
     &                   GRID(ng) % vmask,                              &
     &                   GRID(ng) % pm,                                 &
     &                   GRID(ng) % pn,                                 &
     &                   GRID(ng) % omn,                                &
     &                   GRID(ng) % om_u,                               &
     &                   GRID(ng) % on_v,                               &
     &                   GRID(ng) % Hz,                                 &
     &                   GRID(ng) % z_r,                                &
     &                   MIXING(ng) % Hdiffusion,                       &
     &                   MIXING(ng) % diff3d_u,                         &
     &                   MIXING(ng) % diff3d_v,                         &
     &                   OCEAN(ng) % u,                                 &
     &                   OCEAN(ng) % v)
      CALL wclock_off (ng, iNLM, 28)
      RETURN
      END SUBROUTINE hmixing
!
!***********************************************************************
      SUBROUTINE hmixing_tile (ng, tile,                                &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         IminS, ImaxS, JminS, JmaxS,              &
     &                         nrhs,                                    &
     &                         rmask, umask, vmask,                     &
     &                         pm, pn, omn, om_u, on_v,                 &
     &                         Hz, z_r,                                 &
     &                         Hdiffusion,                              &
     &                         diff3d_u, diff3d_v,                      &
     &                         u, v)
!***********************************************************************
!
      USE mod_param
      USE mod_ncparam
      USE mod_scalars
!
      USE exchange_3d_mod
      USE mp_exchange_mod, ONLY : mp_exchange3d
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nrhs
!
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
      real(r8), intent(in) :: Hdiffusion(LBi:,LBj:)
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: pm(LBi:,LBj:)
      real(r8), intent(in) :: pn(LBi:,LBj:)
      real(r8), intent(in) :: omn(LBi:,LBj:)
      real(r8), intent(in) :: om_u(LBi:,LBj:)
      real(r8), intent(in) :: on_v(LBi:,LBj:)
      real(r8), intent(in) :: u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: v(LBi:,LBj:,:,:)
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(out) :: diff3d_u(LBi:,LBj:,:)
      real(r8), intent(out) :: diff3d_v(LBi:,LBj:,:)
!
!  Local variable declarations.
!
      integer :: i, j, k
      real(r8), parameter :: SmagorCoef = 0.1_r8
      real(r8), parameter :: PecletCoef = 1.0_r8 / 12.0_r8
      real(r8) :: DefRate, cff, clip_diff, clip_scale
!
!-----------------------------------------------------------------------
!  Set lower and upper tile bounds and staggered variables bounds for
!  this horizontal domain partition.  Notice that if tile=-1, it will
!  set the values for the global grid.
!-----------------------------------------------------------------------
!
      integer :: Istr, IstrB, IstrP, IstrR, IstrT, IstrM, IstrU
      integer :: Iend, IendB, IendP, IendR, IendT
      integer :: Jstr, JstrB, JstrP, JstrR, JstrT, JstrM, JstrV
      integer :: Jend, JendB, JendP, JendR, JendT
      integer :: Istrm3, Istrm2, Istrm1, IstrUm2, IstrUm1
      integer :: Iendp1, Iendp2, Iendp2i, Iendp3
      integer :: Jstrm3, Jstrm2, Jstrm1, JstrVm2, JstrVm1
      integer :: Jendp1, Jendp2, Jendp2i, Jendp3
!
      Istr   =BOUNDS(ng) % Istr   (tile)
      IstrB  =BOUNDS(ng) % IstrB  (tile)
      IstrM  =BOUNDS(ng) % IstrM  (tile)
      IstrP  =BOUNDS(ng) % IstrP  (tile)
      IstrR  =BOUNDS(ng) % IstrR  (tile)
      IstrT  =BOUNDS(ng) % IstrT  (tile)
      IstrU  =BOUNDS(ng) % IstrU  (tile)
      Iend   =BOUNDS(ng) % Iend   (tile)
      IendB  =BOUNDS(ng) % IendB  (tile)
      IendP  =BOUNDS(ng) % IendP  (tile)
      IendR  =BOUNDS(ng) % IendR  (tile)
      IendT  =BOUNDS(ng) % IendT  (tile)
      Jstr   =BOUNDS(ng) % Jstr   (tile)
      JstrB  =BOUNDS(ng) % JstrB  (tile)
      JstrM  =BOUNDS(ng) % JstrM  (tile)
      JstrP  =BOUNDS(ng) % JstrP  (tile)
      JstrR  =BOUNDS(ng) % JstrR  (tile)
      JstrT  =BOUNDS(ng) % JstrT  (tile)
      JstrV  =BOUNDS(ng) % JstrV  (tile)
      Jend   =BOUNDS(ng) % Jend   (tile)
      JendB  =BOUNDS(ng) % JendB  (tile)
      JendP  =BOUNDS(ng) % JendP  (tile)
      JendR  =BOUNDS(ng) % JendR  (tile)
      JendT  =BOUNDS(ng) % JendT  (tile)
!
      Istrm3 =BOUNDS(ng) % Istrm3 (tile)            ! Istr-3
      Istrm2 =BOUNDS(ng) % Istrm2 (tile)            ! Istr-2
      Istrm1 =BOUNDS(ng) % Istrm1 (tile)            ! Istr-1
      IstrUm2=BOUNDS(ng) % IstrUm2(tile)            ! IstrU-2
      IstrUm1=BOUNDS(ng) % IstrUm1(tile)            ! IstrU-1
      Iendp1 =BOUNDS(ng) % Iendp1 (tile)            ! Iend+1
      Iendp2 =BOUNDS(ng) % Iendp2 (tile)            ! Iend+2
      Iendp2i=BOUNDS(ng) % Iendp2i(tile)            ! Iend+2 interior
      Iendp3 =BOUNDS(ng) % Iendp3 (tile)            ! Iend+3
      Jstrm3 =BOUNDS(ng) % Jstrm3 (tile)            ! Jstr-3
      Jstrm2 =BOUNDS(ng) % Jstrm2 (tile)            ! Jstr-2
      Jstrm1 =BOUNDS(ng) % Jstrm1 (tile)            ! Jstr-1
      JstrVm2=BOUNDS(ng) % JstrVm2(tile)            ! JstrV-2
      JstrVm1=BOUNDS(ng) % JstrVm1(tile)            ! JstrV-1
      Jendp1 =BOUNDS(ng) % Jendp1 (tile)            ! Jend+1
      Jendp2 =BOUNDS(ng) % Jendp2 (tile)            ! Jend+2
      Jendp2i=BOUNDS(ng) % Jendp2i(tile)            ! Jend+2 interior
      Jendp3 =BOUNDS(ng) % Jendp3 (tile)            ! Jend+3
!
!-----------------------------------------------------------------------
!  Compute tracer horizontal diffusion coefficient as the sum of the
!  biharmonic diffusion term and a diapycnal hyperdiffusion term
!  resulting from the splitting into advective and diffusive components
!  the 3rd-order upstream bias horizontal advection.
!-----------------------------------------------------------------------
!
!  Following Holland et al. (1998) and Webb et al. (1998), the 3rd-order
!  upstream bias horizontal advection can be splitted into advective and
!  diffusive terms. The advective term is just the standard 4th-order
!  cenrtered differences operator.  The diffusive term is a 3rd-order
!  diapycnal hyperdiffusive operator which is proportional to the
!  absolute local velocity (Marchiesello et al., 2008):
!
!         B = 1/12 * ABS(u) * (GridSpace ** 3)
!
!  The resulting diapycnak diffusion, B, is scaled and clipped to achive
!  better stability. The squared-root of the total (diapycnal and
!  biharmonic) diffusion coefficient taken since the harmonic operator
!  is applied twice to achieve biharmonic mixing.
!
!  Squared-root biharmonic diffusion coefficients [SQRT(m4/s)] at
!  U-points.
!
      DO k=1,N(ng)
        DO j=Jstr-1,Jend+1
          DO i=IstrU-1,Iend+1
            diff3d_u(i,j,k)=0.5_r8*(Hdiffusion(i-1,j)+                  &
     &                              Hdiffusion(i  ,j))+                 &
     &                      PecletCoef*(om_u(i,j)**3)*                  &
     &                      ABS(u(i,j,k,nrhs))
            clip_scale=0.5_r8*(Hz(i-1,j,k)+Hz(i,j,k))/                  &
     &                 (z_r(i,j,k)-z_r(i-1,j,k))
            clip_diff=diff3d_u(i,j,k)*                                  &
     &                MIN(1.0_r8, clip_scale*clip_scale)
            diff3d_u(i,j,k)=MIN(clip_diff, diff3d_u(i,j,k))
            diff3d_u(i,j,k)=SQRT(diff3d_u(i,j,k))
            diff3d_u(i,j,k)=diff3d_u(i,j,k)*umask(i,j)
          END DO
        END DO
!
!  Squared-root biharmonic diffusion coefficients [SQRT(m4/s)] at
!  V-points.
!
        DO j=JstrV-1,Jend+1
          DO i=Istr-1,Iend+1
            diff3d_v(i,j,k)=0.5_r8*(Hdiffusion(i,j-1)+                  &
     &                              Hdiffusion(i,j  ))+                 &
     &                      PecletCoef*(on_v(i,j)**3)*                  &
     &                      ABS(v(i,j,k,nrhs))
            clip_scale=0.5_r8*(Hz(i,j,k)+Hz(i,j-1,k))/                  &
     &                 (z_r(i,j,k)-z_r(i,j-1,k))
            clip_diff=diff3d_v(i,j,k)*                                  &
     &                MIN(1.0_r8, clip_scale*clip_scale)
            diff3d_v(i,j,k)=MIN(clip_diff, diff3d_v(i,j,k))
            diff3d_v(i,j,k)=SQRT(diff3d_v(i,j,k))
            diff3d_v(i,j,k)=diff3d_v(i,j,k)*vmask(i,j)
          END DO
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Apply boundary conditions
!-----------------------------------------------------------------------
!
!  East-West gradient boundary conditions.
!
      IF (.not.EWperiodic(ng)) THEN
        IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
          DO k=1,N(ng)
            DO j=JstrV-1,Jend
              diff3d_u(Iend+1,j,k)=diff3d_u(Iend,j,k)
              diff3d_v(Iend+1,j,k)=diff3d_v(Iend,j,k)
            END DO
          END DO
        END IF
        IF (DOMAIN(ng)%Western_Edge(tile)) THEN
          DO k=1,N(ng)
            DO j=JstrV-1,Jend
              diff3d_u(Istr-1,j,k)=diff3d_u(Istr,j,k)
              diff3d_v(Istr-1,j,k)=diff3d_v(Istr,j,k)
            END DO
          END DO
        END IF
      END IF
!
!  North-South gradient boundary conditions.
!
      IF (.not.NSperiodic(ng)) THEN
        IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
          DO k=1,N(ng)
            DO i=IstrU-1,Iend
              diff3d_u(i,Jend+1,k)=diff3d_u(i,Jend,k)
              diff3d_v(i,Jend+1,k)=diff3d_v(i,Jend,k)
            END DO
          END DO
        END IF
        IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
          DO k=1,N(ng)
            DO i=IstrU-1,Iend
              diff3d_u(i,Jstr-1,k)=diff3d_u(i,Jstr,k)
              diff3d_v(i,Jstr-1,k)=diff3d_v(i,Jstr,k)
            END DO
          END DO
        END IF
      END IF
!
!  Boundary corners.
!
      IF (.not.(EWperiodic(ng).or.NSperiodic(ng))) THEN
        IF (DOMAIN(ng)%SouthWest_Corner(tile)) THEN
          DO k=1,N(ng)
            diff3d_u(Istr-1,Jstr-1,k)=0.5_r8*                           &
     &                                (diff3d_u(Istr  ,Jstr-1,k)+       &
     &                                 diff3d_u(Istr-1,Jstr  ,k))
            diff3d_v(Istr-1,Jstr-1,k)=0.5_r8*                           &
     &                                (diff3d_v(Istr  ,Jstr-1,k)+       &
     &                                 diff3d_v(Istr-1,Jstr  ,k))
          END DO
        END IF
        IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
          DO k=1,N(ng)
            diff3d_u(Iend+1,Jstr-1,k)=0.5_r8*                           &
     &                                (diff3d_u(Iend  ,Jstr-1,k)+       &
     &                                 diff3d_u(Iend+1,Jstr  ,k))
            diff3d_v(Iend+1,Jstr-1,k)=0.5_r8*                           &
     &                                (diff3d_v(Iend  ,Jstr-1,k)+       &
     &                                 diff3d_v(Iend+1,Jstr  ,k))
          END DO
        END IF
        IF (DOMAIN(ng)%Northern_Edge(tile).and.                         &
     &      DOMAIN(ng)%Western_Edge(tile)) THEN
          DO k=1,N(ng)
            diff3d_u(Istr-1,Jend+1,k)=0.5_r8*                           &
     &                                (diff3d_u(Istr  ,Jend+1,k)+       &
     &                                 diff3d_u(Istr-1,Jend  ,k))
            diff3d_v(Istr-1,Jend+1,k)=0.5_r8*                           &
     &                                (diff3d_v(Istr  ,Jend+1,k)+       &
     &                                 diff3d_v(Istr-1,Jend  ,k))
          END DO
        END IF
        IF (DOMAIN(ng)%Northern_Edge(tile).and.                         &
     &      DOMAIN(ng)%Eastern_Edge(tile)) THEN
          DO k=1,N(ng)
            diff3d_u(Iend+1,Jend+1,k)=0.5_r8*                           &
     &                                (diff3d_u(Iend  ,Jend+1,k)+       &
     &                                 diff3d_u(Iend+1,Jend  ,k))
            diff3d_v(Iend+1,Jend+1,k)=0.5_r8*                           &
     &                                (diff3d_v(Iend  ,Jend+1,k)+       &
     &                                 diff3d_v(Iend+1,Jend  ,k))
          END DO
        END IF
      END IF
!
!  Periodic boundary conditions.
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_u3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          diff3d_u)
        CALL exchange_v3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          diff3d_v)
      END IF
!
!  Exhange boundary data.
!
      CALL mp_exchange3d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    diff3d_u, diff3d_v)
      RETURN
      END SUBROUTINE hmixing_tile
      END MODULE hmixing_mod
