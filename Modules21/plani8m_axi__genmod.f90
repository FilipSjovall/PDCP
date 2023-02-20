        !COMPILER-GENERATED INTERFACE MODULE: Thu Feb 27 16:53:58 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE PLANI8M_AXI__genmod
          INTERFACE 
            SUBROUTINE PLANI8M_AXI(ME,COORD,RHO,NGP)
              INTEGER(KIND=4) :: NGP
              REAL(KIND=8) :: ME(16,16)
              REAL(KIND=8) :: COORD(8,2)
              REAL(KIND=8) :: RHO
            END SUBROUTINE PLANI8M_AXI
          END INTERFACE 
        END MODULE PLANI8M_AXI__genmod
