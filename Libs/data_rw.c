/*CREATED BY PIERRE-YVES AQUILANTI 2011*/
#include "data_rw.h"

PetscErrorCode writeBinaryScalarArray(const char * name, int nb, PetscScalar * array){
  int file_descriptor;
  PetscErrorCode ierr;


  PetscPrintf(PETSC_COMM_WORLD,"Outputing ScalarArray of size %d to file %s\n",nb,name);
  ierr=PetscBinaryOpen(name,FILE_MODE_WRITE,&file_descriptor);CHKERRQ(ierr);
  ierr=PetscBinarySynchronizedWrite(PETSC_COMM_WORLD,file_descriptor,array,nb,PETSC_SCALAR,PETSC_FALSE);CHKERRQ(ierr);
  ierr=PetscBinaryClose(file_descriptor);CHKERRQ(ierr);

  return ierr;
}

PetscErrorCode writeBinaryVecArray(const char * name, int nb, Vec * array){
  PetscInt size,i;
  PetscScalar * tmp_array;
  PetscErrorCode ierr;


  PetscPrintf(PETSC_COMM_WORLD,"Outputing VecArray of size %d to file %s\n",nb,name);
  ierr=VecGetSize(array[0],&size);CHKERRQ(ierr);

  for(i=0;i<nb;i++)
    ierr=VecGetArray(array[i],&tmp_array);CHKERRQ(ierr);
  ierr=writeBinaryScalarArray(name,nb*size,tmp_array);CHKERRQ(ierr);


  for(i=0;i<nb;i++)
    ierr=VecRestoreArray(array[i],&tmp_array);CHKERRQ(ierr);

  return ierr;
}

PetscErrorCode readBinaryScalarArray(const char * name, int * nb, PetscScalar * array){
  int file_descriptor;
  PetscErrorCode ierr;
  long size;

  getFileSize(name,&size);

  if(*nb<=0) *nb=(int)size/((int)sizeof(PetscScalar));
  if(size/sizeof(PetscScalar)!=*nb) {
    PetscPrintf(PETSC_COMM_WORLD,"Wrong ScalarArray size %d from file %s (supposed to be %ld (%ld))\n",*nb,name,size/sizeof(PetscScalar),size);
    return 1;
  }


  PetscPrintf(PETSC_COMM_WORLD,"Reading ScalarArray of size %ld from file %s\n",*nb,name);
  ierr=PetscBinaryOpen(name,FILE_MODE_READ,&file_descriptor);CHKERRQ(ierr);
  ierr=PetscBinarySynchronizedRead(PETSC_COMM_WORLD,file_descriptor,array,*nb,PETSC_SCALAR);CHKERRQ(ierr);
  ierr=PetscBinaryClose(file_descriptor);CHKERRQ(ierr);

  return ierr;
}

PetscErrorCode readBinaryVecArray(const char * name, int * nb, Vec * array){
  PetscInt size,i;
  PetscScalar * tmp_array;
  PetscErrorCode ierr;
    PetscPrintf(PETSC_COMM_WORLD,"Reading VecArray of size %d from file %s\n",*nb,name);
  ierr=VecGetSize(array[0],&size);CHKERRQ(ierr);

  for(i=0;i<*nb;i++)
    ierr=VecGetArray(array[i],&tmp_array);CHKERRQ(ierr);

  size*=(*nb);
  ierr=readBinaryScalarArray(name,&size,tmp_array);CHKERRQ(ierr);

  for(i=0;i<*nb;i++)
    ierr=VecRestoreArray(array[i],&tmp_array);CHKERRQ(ierr);

  return ierr;
}


PetscErrorCode getFileSize(const char * name, long * size){
  FILE * fptr;
  *size = 0L;

#ifdef LINUX
  struct stat fs;

  if(stat(name,&fs)!=0){
    perror("Cannot state file\n");
  }
  *size=fs.st_size;

#else
fptr=fopen(name,"rb");
  if(fptr!=NULL){
    fseek(fptr,0L,SEEK_END);
    *size = ftell(fptr);
    fclose(fptr);
  }
#endif

  printf("<<<<<<<<<<<<<<<<<<<<<<<<<<< size=%ld\n",*size);

  return 0;
}
