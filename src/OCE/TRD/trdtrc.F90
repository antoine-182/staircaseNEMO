MODULE trdtrc
   !!======================================================================
   !!                       ***  MODULE trdtrc  ***
   !!  Dummy module
   !!======================================================================
   !!----------------------------------------------------------------------
   !!   Dummy module                                             NO TOP use
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE trd_trc( ptrtrd, kjn, ktrd, kt, Kmm )
      INTEGER ::   kt, kjn, ktrd   
      INTEGER ::   Kmm            ! time level index
      REAL    ::   ptrtrd(:,:,:)  
      WRITE(*,*) 'trd_trc : You should not have seen this print! error?', ptrtrd(1,1,1)
      WRITE(*,*) '  "      "      : You should not have seen this print! error?', kjn, ktrd, kt
   END SUBROUTINE trd_trc

   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: trdtrc.F90 12377 2020-02-12 14:39:06Z acc $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!======================================================================
END MODULE trdtrc
