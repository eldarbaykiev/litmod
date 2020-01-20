! xyz2vtk.f


      CHARACTER(100)name_f,name_nodes,name_vtk,fields(5)
      REAL(8),ALLOCATABLE:: COOR(:,:,:,:),FLD(:,:,:,:) 
      INTEGER,ALLOCATABLE::MAT(:,:,:)
      REAL(8) L_y_max,d_x,d_y,d_z,O_x,O_y,O_z
      INTEGER N_x,N_y,N_z,N_cells,n_node,CELL(8) 
      fields(1)='Temperature'  
      fields(2)='Density'
      fields(3)='Vp'
      fields(4)='Vs'
      fields(5)='Material'
      WRITE(*,*)'N_x,N_y,N_z...?'
      READ(*,*)N_x,N_y,N_z
      WRITE(*,*)'NAME OF FILE (WITHOUT .xyz)...?'
      READ(*,*)name_f
      WRITE(*,*)'NODE FILE NAME...?'
      READ(*,*)name_nodes
      ALLOCATE (COOR(0:N_x-1,0:N_y-1,0:N_z-1,3),
     *          FLD(0:N_x-1,0:N_y-1,0:N_z-1,4),
     *          MAT(0:N_x-1,0:N_y-1,0:N_z-1))
      name_vtk=TRIM(name_f)//'.vtk'
      name_f=TRIM(name_f)//'.xyz'
      name_nodes=TRIM(name_nodes)
      OPEN(1,file=TRIM(name_f),status='OLD')
      OPEN(2,file=TRIM(name_vtk),status='UNKNOWN')
      OPEN(3,file=TRIM(name_nodes),status='OLD')
      DO I_y=0,N_y-1
       DO I_x=0,N_x-1
        DO I_z=0,N_z-1
           READ(1,*)COOR(I_x,I_y,I_z,1),COOR(I_x,I_y,I_z,2),
     *               COOR(I_x,I_y,I_z,3),FLD(I_x,I_y,I_z,1),
     *               FLD(I_x,I_y,I_z,2),dm,FLD(I_x,I_y,I_z,3),
     *               FLD(I_x,I_y,I_z,4)           
           READ(3,*) dm,dm,dm,MAT(I_x,I_y,I_z) 
         ENDDO
        ENDDO
       ENDDO
        
       d_z=COOR(0,0,1,3)-COOR(0,0,0,3)
       d_x=COOR(0,1,0,2)-COOR(0,0,0,2)
       d_y=COOR(1,0,0,1)-COOR(0,0,0,1)
       O_x=COOR(0,0,0,1)
       O_y=COOR(0,0,0,2)
       O_z=COOR(0,0,0,3)
       L_y_max=MAXVAL(COOR(:,:,:,2))   
    
       WRITE(2,'(A)')'# vtk DataFile Version 3.0'
       WRITE(2,'(A)')TRIM(name_f)
       WRITE(2,'(A)')'ASCII'
       WRITE(2,'(A)')'DATASET STRUCTURED_POINTS'
       WRITE(2,'(A,1X,3I4)')'DIMENSIONS',N_x,N_y,N_z
       WRITE(2,'(A,1X,3F10.2)')'ORIGIN',O_x,O_y,O_z  
       WRITE(2,'(A,1X,3F10.2)')'SPACING',d_x,d_y,d_z
       write(*,*)'d_x,d_y,d_z:', d_x,d_y,d_z
       write(*,*)'A description of the material code can be found
     * in file: layers.info'
       write(*,*)'-100-->air; -200-->water; 0-->sub-lithospheric mantle'
!       DO I_x=0,N_x-1
!        DO I_y=0,N_y-1
!         DO I_z=0,N_z-1
!          WRITE(2,*)COOR(I_x,I_y,I_z,1),L_y_max-COOR(I_x,I_y,I_z,2),
!     *              COOR(I_x,I_y,I_z,3)
!         ENDDO
!        ENDDO
!       ENDDO
!       N_cells=(N_x-1)*(N_y-1)*(N_z-1)
!       WRITE(2,'(A,1X,I8,1X,I8)')'CELLS',N_cells,N_cells*9

!       DO I_x=0,N_x-2
!        DO I_y=0,N_y-2
!         DO I_z=0,N_z-2
!          ind=1
!          DO k=I_z,I_z+1
!           DO j=I_y,I_y+1
!            DO i=I_x,I_x+1
!             n_node=i*N_z*N_y+j*N_z+k
!             CELL(ind)=n_node
!             ind=ind+1
!            ENDDO
!           ENDDO
!          ENDDO  
!          WRITE(2,'(A,1X,I8,1X,I8,1X,I8,1X,I8,1X,I8,1X,I8,1X,I8,1X,I8)')
!     *         '8',CELL(:)
!         ENDDO
!        ENDDO
!       ENDDO 
       
!       WRITE(2,'(A,1X,I8)')'CELL_TYPES',N_cells
!       DO i=1,N_cells
!        WRITE(2,*)'11'
!       ENDDO

       WRITE(2,'(A,1X,I8)')'POINT_DATA',N_x*N_y*N_z
       DO i=1,5
        IF(i==5) THEN
         WRITE(2,'(A,1X,A,1X,A)')'SCALARS',
     *                                TRIM(fields(i)),'integer 1'
        ELSE
         WRITE(2,'(A,1X,A,1X,A)')'SCALARS',
     *                                TRIM(fields(i)),'double 1'
        ENDIF
        WRITE(2,'(A)')'LOOKUP_TABLE default'
        DO I_z=0,N_z-1
         DO I_y=0,N_y-1
          DO I_x=0,N_x-1
           IF(i==5) THEN
            WRITE(2,'(I4)')MAT(I_x,I_y,I_z)
           ELSE
            WRITE(2,*)FLD(I_x,I_y,I_z,i)
           ENDIF
          ENDDO
         ENDDO
        ENDDO

       ENDDO
       CLOSE(1)      
       CLOSE(2)
       CLOSE(3)

      END 


