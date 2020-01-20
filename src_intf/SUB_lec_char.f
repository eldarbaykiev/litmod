! SUB_lec_char.f
       
       SUBROUTINE lec_char(U_lec,D_lec,TEXT)
       REAL(4)U_lec,D_lec,M_lec
       INTEGER i_lec
       CHARACTER(80) TEXT
       CHARACTER(1) dmch
       ILET=0
         TEXT=''
         CALL PGSCI(15)
         CALL PGRECT(0.,0.3,U_lec,D_lec)
         M_lec=(U_lec+D_lec)/2
         DO
          call pgband(0,1,0.,0.,xdm,ydm,dmch)
          i_lec=ichar(dmch)
          IF (i_lec==13) THEN
             IF(ILET>0) EXIT
             CYCLE
          ELSEIF (i_lec>47.AND.i_lec<58.OR.i_lec==45) THEN
             ILET=ILET+1
             TEXT(ILET:ILET)=dmch
          ELSEIF (ICHAR(dmch)==127.OR.ICHAR(dmch)==8) THEN
             IF(ILET<1) CYCLE
             TEXT(ILET:ILET)=''
             ILET=ILET-1
             IF (ILET==0) THEN
                ILET=1
                TEXT=''
             ENDIF
          ELSE
             IF((ICHAR(dmch)==67.OR.ICHAR(dmch)==99))THEN
              TEXT='-1'
              EXIT
             ELSE
              write(*,*)'INVALID NUMBER !'
             ENDIF
          ENDIF
         CALL PGSCI(15)
         CALL PGRECT(0.,0.3,U_lec,D_lec)
         CALL PGSCI(1)
         CALL PGPTEXT(.0,M_lec,0.,0.0,TRIM(TEXT))
         ENDDO
         CALL PGSCI(1)

         END SUBROUTINE
