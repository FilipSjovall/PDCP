        !COMPILER-GENERATED INTERFACE MODULE: Thu Feb 27 16:53:58 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE PLANI8E_AXI__genmod
          INTERFACE 
            SUBROUTINE PLANI8E_AXI(KE,COORD,D_EL,ED,ES,NGP)
              INTEGER(KIND=4) :: NGP
              REAL(KIND=8) :: KE(16,16)
              REAL(KIND=8) :: COORD(8,2)
              REAL(KIND=8) :: D_EL(NGP,4,4)
              REAL(KIND=8) :: ED(16)
              REAL(KIND=8) :: ES(NGP,3,3)
            END SUBROUTINE PLANI8E_AXI
          END INTERFACE 
        END MODULE PLANI8E_AXI__genmod
