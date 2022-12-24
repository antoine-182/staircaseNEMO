#if defined key_agrif
   SUBROUTINE Agrif_connection
   use dom_oce
   use agrif_parameters
      !!----------------------------------------------------------------------
      !!                 *** ROUTINE  Agrif_Sponge ***
      !!----------------------------------------------------------------------
      INTEGER  ::   ji, jj, ind1, ind2
      INTEGER  ::   ispongearea
      REAL(wp) ::   z1_spongearea
      !!----------------------------------------------------------------------
      !
         ! Define ramp from boundaries towards domain interior at T-points
         ! Store it in ztabramp

         ALLOCATE(ztabramp(jpi,jpj))
         ispongearea  = 1 + npt_connect * Agrif_irhox()
         z1_spongearea = 1._wp / REAL( ispongearea )
         
         ztabramp(:,:) = 0._wp

         ! --- West --- !
         IF( (nbondi == -1) .OR. (nbondi == 2) ) THEN
            ind1 = 1+nbghostcells
            ind2 = 1+nbghostcells + ispongearea 
            DO jj = 1, jpj
               DO ji = ind1, ind2                
                  ztabramp(ji,jj) = REAL( ind2 - ji ) * z1_spongearea * umask(ind1,jj,1)
               END DO
            ENDDO
         ENDIF

         ! --- East --- !
         IF( (nbondi == 1) .OR. (nbondi == 2) ) THEN
            ind1 = nlci - nbghostcells - ispongearea
            ind2 = nlci - nbghostcells
            DO jj = 1, jpj
               DO ji = ind1, ind2
                  ztabramp(ji,jj) = MAX( ztabramp(ji,jj), REAL( ji - ind1 ) * z1_spongearea * umask(ind2-1,jj,1) )
               ENDDO
            ENDDO
         ENDIF

         ! --- South --- !
         IF( (nbondj == -1) .OR. (nbondj == 2) ) THEN
            ind1 = 1+nbghostcells
            ind2 = 1+nbghostcells + ispongearea
            DO jj = ind1, ind2 
               DO ji = 1, jpi
                  ztabramp(ji,jj) = MAX( ztabramp(ji,jj), REAL( ind2 - jj ) * z1_spongearea * vmask(ji,ind1,1) )
               END DO
            ENDDO
         ENDIF

         ! --- North --- !
         IF( (nbondj == 1) .OR. (nbondj == 2) ) THEN
            ind1 = nlcj - nbghostcells - ispongearea
            ind2 = nlcj - nbghostcells
            DO jj = ind1, ind2
               DO ji = 1, jpi
                  ztabramp(ji,jj) = MAX( ztabramp(ji,jj), REAL( jj - ind1 ) * z1_spongearea * vmask(ji,ind2-1,1) )
               END DO
            ENDDO
         ENDIF
      !
      !
   END SUBROUTINE Agrif_connection
   
#else
subroutine agrif_connection_empty
end subroutine agrif_connection_empty
#endif