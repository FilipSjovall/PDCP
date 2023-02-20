        !COMPILER-GENERATED INTERFACE MODULE: Thu Feb 27 16:53:58 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE PLANI8F_AXI__genmod
          INTERFACE 
            SUBROUTINE PLANI8F_AXI(FE,COORD,ED,ES,NGP)
              INTEGER(KIND=4) :: NGP
              REAL(KIND=8) :: FE(1,16)
              REAL(KIND=8) :: COORD(8,2)
              REAL(KIND=8) :: ED(16)
              REAL(KIND=8) :: ES(NGP,3,3)
            END SUBROUTINE PLANI8F_AXI
          END INTERFACE 
        END MODULE PLANI8F_AXI__genmod
