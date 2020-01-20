! xyz2vtk.f


      CHARACTER(100)name_f,name_vtk,fields(4)
      REAL(8),ALLOCATABLE:: COOR(:,:,:,:),FLD(:,:,:,:) 
      REAL(8) L_y_max
      INTEGER N_x,N_y,N_z,N_cells,n_node,CELL(8) 
      fields(1)='Temperature'  
      fields(2)='Density'
      fields(3)='Vp'
      fields(4)='Vs'
      WRITE(*,*)'N_x,N_y,N_z...?'
      READ(*,*)N_x,N_y,N_z
      WRITE(*,*)'NAME OF FILE (WITHOUT .xyz)...?'
      READ(*,*)name_f
      ALLOCATE (COOR(0:N_x-1,0:N_y-1,0:N_z-1,3),
     *          FLD(0:N_x-1,0:N_y-1,0:N_z-1,4))
      name_vtk=TRIM(name_f)//'.vtk'
      name_f=TRIM(name_f)//'.xyz'
      OPEN(1,file=TRIM(name_f),status='OLD')
      OPEN(2,file=TRIM(name_vtk),status='UNKNOWN')
      DO I_y=0,N_y-1
       DO I_x=0,N_x-1
        DO I_z=0,N_z-1
           READ(1,*)COOR(I_x,I_y,I_z,1),COOR(I_x,I_y,I_z,2),
     *               COOR(I_x,I_y,I_z,3),FLD(I_x,I_y,I_z,1),
     *               FLD(I_x,I_y,I_z,2),dm,FLD(I_x,I_y,I_z,3),
     *               FLD(I_x,I_y,I_z,4)           
         ENDDO
        ENDDO
       ENDDO

       L_y_max=MAXVAL(COOR(:,:,:,2))   
    
       WRITE(2,'(A)')'# vtk DataFile Version 3.0'
       WRITE(2,'(A)')TRIM(name_f)
       WRITE(2,'(A)')'ASCII'
       WRITE(2,'(A)')'DATASET UNSTRUCTURED_GRID'
       WRITE(2,'(A,1X,I8,1X,A)')'POINTS',N_x*N_y*N_z,'float'
       DO I_x=0,N_x-1
        DO I_y=0,N_y-1
         DO I_z=0,N_z-1
          WRITE(2,*)COOR(I_x,I_y,I_z,1),L_y_max-COOR(I_x,I_y,I_z,2),
     *              COOR(I_x,I_y,I_z,3)
         ENDDO
        ENDDO
       ENDDO
       N_cells=(N_x-1)*(N_y-1)*(N_z-1)
       WRITE(2,'(A,1X,I8,1X,I8)')'CELLS',N_cells,N_cells*9

       DO I_x=0,N_x-2
        DO I_y=0,N_y-2
         DO I_z=0,N_z-2
          ind=1
          DO k=I_z,I_z+1
           DO j=I_y,I_y+1
            DO i=I_x,I_x+1
             n_node=i*N_z*N_y+j*N_z+k
             CELL(ind)=n_node
             ind=ind+1
            ENDDO
           ENDDO
          ENDDO  
          WRITE(2,'(A,1X,I8,1X,I8,1X,I8,1X,I8,1X,I8,1X,I8,1X,I8,1X,I8)')
     *         '8',CELL(:)
         ENDDO
        ENDDO
       ENDDO 
       
       WRITE(2,'(A,1X,I8)')'CELL_TYPES',N_cells
       DO i=1,N_cells
        WRITE(2,*)'11'
       ENDDO

       WRITE(2,'(A,1X,I8)')'POINT_DATA',N_x*N_y*N_z
       DO i=1,4

        WRITE(2,'(A,1X,A,1X,A)')'SCALARS',
     *                                TRIM(fields(i)),'float 1'
        WRITE(2,'(A)')'LOOKUP_TABLE default'
        DO I_x=0,N_x-1
         DO I_y=0,N_y-1
          DO I_z=0,N_z-1
           WRITE(2,*)FLD(I_x,I_y,I_z,i)
          ENDDO
         ENDDO
        ENDDO

       ENDDO
       CLOSE(1)      
       CLOSE(2)

      END 


