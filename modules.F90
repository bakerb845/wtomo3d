      MODULE INIT_MODULE
         IMPLICIT NONE
         REAL dx, dy, dz, xorig, yorig, zorig
         INTEGER nx, ny, nz
         LOGICAL freesurf(6), usemin
         SAVE
      END MODULE

      MODULE MODEL_MODULE
         IMPLICIT NONE
         REAL, ALLOCATABLE, DIMENSION(:,:,:) :: rho, rhop, vpr, vsr,  &
                                                qp, qs, hess3p, hess3s
         COMPLEX, ALLOCATABLE, DIMENSION(:,:,:) :: mu, da, mup, dap
         COMPLEX omega, keiy
         REAL kylen
         INTEGER iom, iky!, mx, my, mz
         LOGICAL qpex, qsex
         SAVE
      END MODULE MODEL_MODULE

      MODULE SRCPRM_MODULE
         IMPLICIT NONE
         COMPLEX, ALLOCATABLE, DIMENSION(:,:) :: source
         REAL, ALLOCATABLE, DIMENSION(:,:) :: dsd
         INTEGER, ALLOCATABLE, DIMENSION(:) :: isg
         INTEGER ipadsrc, nsg
         SAVE
      END MODULE SRCPRM_MODULE

      MODULE SOURCE_MODULE
         IMPLICIT NONE
         CHARACTER(2), ALLOCATABLE :: srctyp(:)
         REAL, ALLOCATABLE :: ain(:)
         INTEGER, ALLOCATABLE :: modenum(:)
         SAVE
      END MODULE SOURCE_MODULE

      MODULE RECPRM_MODULE
         IMPLICIT NONE
         COMPLEX, ALLOCATABLE, DIMENSION(:,:,:) :: receiver 
         SAVE
      END MODULE RECPRM_MODULE

      MODULE AGC_MODULE
         IMPLICIT NONE
         REAL, ALLOCATABLE, DIMENSION(:,:,:) :: agcpv, agcsv 
         SAVE
      END MODULE AGC_MODULE

      MODULE DEMUX_MODULE
         IMPLICIT NONE
         REAL, ALLOCATABLE, DIMENSION(:,:) :: prof(:,:)
         SAVE
      END MODULE

      MODULE PML_MODULE
         IMPLICIT NONE
         REAL pmld, pmlr, pmlf
         INTEGER ipml
         LOGICAL pml
         SAVE
      END MODULE

      MODULE COEFFS_MODULE
          IMPLICIT NONE

          COMPLEX lmu1, lmu2, lmu3, lmu4, lmu5, lmu6, lmu7, lmu8
          COMPLEX lpu1, lpu2, lpu3, lpu4, lpu5, lpu6, lpu7, lpu8
          COMPLEX lp2u1, lp2u2, lp2u3, lp2u4, lp2u5, lp2u6, lp2u7, lp2u8
          COMPLEX u1, u2, u3, u4, u5, u6, u7, u8
          COMPLEX da1, da2, da3, da4, da5, da6, da7, da8

          COMPLEX  daloc(3, 3, 3)
          COMPLEX  muloc(3, 3, 3)
          REAL     roloc(3, 3, 3)
          REAL     w1, w2, w3, wm1, wm2, wm3, wm4
      END MODULE COEFFS_MODULE 

      MODULE MAT_MODULE
         IMPLICIT NONE
         COMPLEX a1mUU, a1mVV, a1mWW, a1mUV, a1mVU, a1mUW, a1mWU, a1mVW, a1mWV, &
                 a2mUU, a2mVV, a2mWW, a2mUV, a2mVU, a2mUW, a2mWU, a2mVW, a2mWV, &
                 a3mUU, a3mVV, a3mWW, a3mUV, a3mVU, a3mUW, a3mWU, a3mVW, a3mWV, &
                 a4mUU, a4mVV, a4mWW, a4mUV, a4mVU, a4mUW, a4mWU, a4mVW, a4mWV, &
                 a5mUU, a5mVV, a5mWW, a5mUV, a5mVU, a5mUW, a5mWU, a5mVW, a5mWV, &
                 a6mUU, a6mVV, a6mWW, a6mUV, a6mVU, a6mUW, a6mWU, a6mVW, a6mWV, &
                 a7mUU, a7mVV, a7mWW, a7mUV, a7mVU, a7mUW, a7mWU, a7mVW, a7mWV, &
                 a8mUU, a8mVV, a8mWW, a8mUV, a8mVU, a8mUW, a8mWU, a8mVW, a8mWV, &
                 a9mUU, a9mVV, a9mWW, a9mUV, a9mVU, a9mUW, a9mWU, a9mVW, a9mWV, &
                 a1nUU, a1nVV, a1nWW, a1nUV, a1nVU, a1nUW, a1nWU, a1nVW, a1nWV, &
                 a2nUU, a2nVV, a2nWW, a2nUV, a2nVU, a2nUW, a2nWU, a2nVW, a2nWV, &
                 a3nUU, a3nVV, a3nWW, a3nUV, a3nVU, a3nUW, a3nWU, a3nVW, a3nWV, &
                 a4nUU, a4nVV, a4nWW, a4nUV, a4nVU, a4nUW, a4nWU, a4nVW, a4nWV, &
                 a5nUU, a5nVV, a5nWW, a5nUV, a5nVU, a5nUW, a5nWU, a5nVW, a5nWV, &
                 a6nUU, a6nVV, a6nWW, a6nUV, a6nVU, a6nUW, a6nWU, a6nVW, a6nWV, &
                 a7nUU, a7nVV, a7nWW, a7nUV, a7nVU, a7nUW, a7nWU, a7nVW, a7nWV, &
                 a8nUU, a8nVV, a8nWW, a8nUV, a8nVU, a8nUW, a8nWU, a8nVW, a8nWV, &
                 a9nUU, a9nVV, a9nWW, a9nUV, a9nVU, a9nUW, a9nWU, a9nVW, a9nWV, &
                 a1pUU, a1pVV, a1pWW, a1pUV, a1pVU, a1pUW, a1pWU, a1pVW, a1pWV, &
                 a2pUU, a2pVV, a2pWW, a2pUV, a2pVU, a2pUW, a2pWU, a2pVW, a2pWV, &
                 a3pUU, a3pVV, a3pWW, a3pUV, a3pVU, a3pUW, a3pWU, a3pVW, a3pWV, &
                 a4pUU, a4pVV, a4pWW, a4pUV, a4pVU, a4pUW, a4pWU, a4pVW, a4pWV, &
                 a5pUU, a5pVV, a5pWW, a5pUV, a5pVU, a5pUW, a5pWU, a5pVW, a5pWV, &
                 a6pUU, a6pVV, a6pWW, a6pUV, a6pVU, a6pUW, a6pWU, a6pVW, a6pWV, &
                 a7pUU, a7pVV, a7pWW, a7pUV, a7pVU, a7pUW, a7pWU, a7pVW, a7pWV, &
                 a8pUU, a8pVV, a8pWW, a8pUV, a8pVU, a8pUW, a8pWU, a8pVW, a8pWV, &
                 a9pUU, a9pVV, a9pWW, a9pUV, a9pVU, a9pUW, a9pWU, a9pVW, a9pWV
      END MODULE MAT_MODULE

      MODULE MATBG_MODULE
         COMPLEX a1mUU1, a1mVV1, a1mWW1, a1mUV1, a1mVU1, a1mUW1, a1mWU1, a1mVW1, a1mWV1, &
                 a2mUU1, a2mVV1, a2mWW1, a2mUV1, a2mVU1, a2mUW1, a2mWU1, a2mVW1, a2mWV1, &
                 a3mUU1, a3mVV1, a3mWW1, a3mUV1, a3mVU1, a3mUW1, a3mWU1, a3mVW1, a3mWV1, &
                 a4mUU1, a4mVV1, a4mWW1, a4mUV1, a4mVU1, a4mUW1, a4mWU1, a4mVW1, a4mWV1, &
                 a5mUU1, a5mVV1, a5mWW1, a5mUV1, a5mVU1, a5mUW1, a5mWU1, a5mVW1, a5mWV1, &
                 a6mUU1, a6mVV1, a6mWW1, a6mUV1, a6mVU1, a6mUW1, a6mWU1, a6mVW1, a6mWV1, &
                 a7mUU1, a7mVV1, a7mWW1, a7mUV1, a7mVU1, a7mUW1, a7mWU1, a7mVW1, a7mWV1, &
                 a8mUU1, a8mVV1, a8mWW1, a8mUV1, a8mVU1, a8mUW1, a8mWU1, a8mVW1, a8mWV1, &
                 a9mUU1, a9mVV1, a9mWW1, a9mUV1, a9mVU1, a9mUW1, a9mWU1, a9mVW1, a9mWV1, &
                 a1nUU1, a1nVV1, a1nWW1, a1nUV1, a1nVU1, a1nUW1, a1nWU1, a1nVW1, a1nWV1, &
                 a2nUU1, a2nVV1, a2nWW1, a2nUV1, a2nVU1, a2nUW1, a2nWU1, a2nVW1, a2nWV1, &
                 a3nUU1, a3nVV1, a3nWW1, a3nUV1, a3nVU1, a3nUW1, a3nWU1, a3nVW1, a3nWV1, &
                 a4nUU1, a4nVV1, a4nWW1, a4nUV1, a4nVU1, a4nUW1, a4nWU1, a4nVW1, a4nWV1, &
                 a5nUU1, a5nVV1, a5nWW1, a5nUV1, a5nVU1, a5nUW1, a5nWU1, a5nVW1, a5nWV1, &
                 a6nUU1, a6nVV1, a6nWW1, a6nUV1, a6nVU1, a6nUW1, a6nWU1, a6nVW1, a6nWV1, &
                 a7nUU1, a7nVV1, a7nWW1, a7nUV1, a7nVU1, a7nUW1, a7nWU1, a7nVW1, a7nWV1, &
                 a8nUU1, a8nVV1, a8nWW1, a8nUV1, a8nVU1, a8nUW1, a8nWU1, a8nVW1, a8nWV1, &
                 a9nUU1, a9nVV1, a9nWW1, a9nUV1, a9nVU1, a9nUW1, a9nWU1, a9nVW1, a9nWV1, &
                 a1pUU1, a1pVV1, a1pWW1, a1pUV1, a1pVU1, a1pUW1, a1pWU1, a1pVW1, a1pWV1, &
                 a2pUU1, a2pVV1, a2pWW1, a2pUV1, a2pVU1, a2pUW1, a2pWU1, a2pVW1, a2pWV1, &
                 a3pUU1, a3pVV1, a3pWW1, a3pUV1, a3pVU1, a3pUW1, a3pWU1, a3pVW1, a3pWV1, &
                 a4pUU1, a4pVV1, a4pWW1, a4pUV1, a4pVU1, a4pUW1, a4pWU1, a4pVW1, a4pWV1, &
                 a5pUU1, a5pVV1, a5pWW1, a5pUV1, a5pVU1, a5pUW1, a5pWU1, a5pVW1, a5pWV1, &
                 a6pUU1, a6pVV1, a6pWW1, a6pUV1, a6pVU1, a6pUW1, a6pWU1, a6pVW1, a6pWV1, &
                 a7pUU1, a7pVV1, a7pWW1, a7pUV1, a7pVU1, a7pUW1, a7pWU1, a7pVW1, a7pWV1, &
                 a8pUU1, a8pVV1, a8pWW1, a8pUV1, a8pVU1, a8pUW1, a8pWU1, a8pVW1, a8pWV1, &
                 a9pUU1, a9pVV1, a9pWW1, a9pUV1, a9pVU1, a9pUW1, a9pWU1, a9pVW1, a9pWV1
      END MODULE MATBG_MODULE

      MODULE MATF_MODULE
         IMPLICIT NONE
         COMPLEX e2uu, e2uv, e2vu, e2vv, &
                 e3uu, e3uv, e3vu, e3vv, &
                 e8uu, e8uv, e8vu, e8vv, &
                 e9uu, e9uv, e9vu, e9vv
      END MODULE

      MODULE GRADIENTS_MODULE
         IMPLICIT NONE
         COMPLEX, ALLOCATABLE, DIMENSION(:,:,:) :: gp1, gs1, gp2, gs2
         SAVE
      END MODULE GRADIENTS_MODULE

      MODULE PRECON_MODULE
         IMPLICIT NONE
         REAL, ALLOCATABLE, DIMENSION(:,:,:) :: preconp, precons
         SAVE
      END MODULE

      MODULE HASKELL_MODULE
         IMPLICIT NONE
         COMPLEX, ALLOCATABLE, DIMENSION(:,:) :: ugrn, vgrn, wgrn
         COMPLEX, ALLOCATABLE, DIMENSION(:) :: ugrn1f, vgrn1f, wgrn1f
         REAL, ALLOCATABLE, DIMENSION(:) :: vp1d, vs1d, vp1di, vs1di, &
                                            rh1d, h1d, vp1disc, vs1disc
         INTEGER nl
         SAVE
      END MODULE HASKELL_MODULE

      MODULE PYTABS_MODULE
         IMPLICIT NONE
         REAL, ALLOCATABLE, DIMENSION(:,:,:) :: pytab
         REAL, ALLOCATABLE, DIMENSION(:) :: pyz(:)
         SAVE
      END MODULE PYTABS_MODULE
