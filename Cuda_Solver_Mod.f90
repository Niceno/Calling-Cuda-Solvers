!==============================================================================!
!   Interface to cusolverDn and CUDA C functions                               !
!                                                                              !
!   Compile with:                                                              !
!   nvfortran -i8 -L /usr/local/cuda/lib64 another.f90 \                       !
!             -lcudart -lcublas -lcusparse -lcusolver                          !
!------------------------------------------------------------------------------!
  module Cuda_Solver_Mod

  integer*8, parameter :: CUDA_MEM_CPY_HOST_TO_DEVICE = 1
  integer*8, parameter :: CUDA_MEM_CPY_DEVICE_TO_HOST = 2

  ! Parameters to specify operation to be performed with the dense matrix
  ! (Found in: /usr/local/cuda/targets/x86_64-linux/include/cublas_api.h)
  integer*4, parameter :: CUBLAS_OP_N = 0
  integer*4, parameter :: CUBLAS_OP_T = 1

  ! More parameters for solvers can be found here:
  ! /usr/local/cuda/targets/x86_64-linux/include/cusolver_common.h

  interface

    !-----------------!
    !   Cuda_Malloc   !
    !-----------------!
    integer(c_int) function Cuda_Malloc  &
     (buffer,                            &
      size)                              &
      bind (C, name="cudaMalloc")

      use iso_c_binding
      implicit none
      type(c_ptr)              :: buffer
      integer(c_size_t), value :: size
    end function

    !------------------!
    !   Cuda_Mem_Cpy   !
    !------------------!
    integer(c_int) function Cuda_Mem_Cpy  &
     (dst,                                &
      src,                                &
      count,                              &
      kind)                               &
      bind (C, name="cudaMemcpy" )
      ! note for kind: CUDA_MEM_CPY_HOST_TO_DEVICE = 1
      !                CUDA_MEM_CPY_DEVICE_TO_HOST = 2

      use iso_c_binding
      implicit none
      type(c_ptr),       value :: dst
      type(c_ptr),       value :: src
      integer(c_size_t), value :: count
      integer(c_size_t), value :: kind
    end function

    !---------------!
    !   Cuda_Free   !
    !---------------!
    integer(c_int) function Cuda_Free  &
      (buffer)                         &
      bind(C, name="cudaFree")

      use iso_c_binding
      implicit none
      type(c_ptr), value :: buffer
    end function

    !-----------------------!
    !   Cuda_Mem_Get_Info   !
    !-----------------------!
    integer(c_int) function Cuda_Mem_Get_Info  &
      (fre, tot)                               &
      bind(C, name="cudaMemGetInfo")

      use iso_c_binding
      implicit none

      type(c_ptr), value :: fre
      type(c_ptr), value :: tot
    end function

    !-------------------------!
    !   Cu_Solver_Dn_Create   !
    !-------------------------!
    integer(c_int) function  Cu_Solver_Dn_Create  &
      (cusolver_hndl)                             &
      bind(C, name = "cusolverDnCreate")

      use iso_c_binding
      implicit none

      type(c_ptr) :: cusolver_hndl
    end function

    !--------------------------!
    !   Cu_Solver_Dn_Destroy   !
    !--------------------------!
    integer(c_int) function Cu_Solver_Dn_Destroy  &
      (cusolver_hndl)                             &
      bind(C, name = "cusolverDnDestroy")

      use iso_c_binding
      implicit none

      type(c_ptr), value :: cusolver_hndl
    end function

    !-------------------------------------!
    !   Cu_Solver_Dn_Dgetrf_Buffer_Size   !
    !-------------------------------------!
    integer(c_int) function Cu_Solver_Dn_Dgetrf_Buffer_Size  &
     (cusolver_hndl,                                         &
      m,                                                     &
      n,                                                     &
      pnt_a_gpu,                                             &
      lda,                                                   &
      lwork)                                                 &
      bind(C,name="cusolverDnDgetrf_bufferSize")

      use iso_c_binding
      implicit none

      type   (c_ptr), value :: cusolver_hndl   ! CUDA solver handle
      integer(c_int), value :: m
      integer(c_int), value :: n
      type   (c_ptr), value :: pnt_a_gpu
      integer(c_int), value :: lda
      type   (c_ptr), value :: lwork           ! size of workspace
    end function

    !-------------------------!
    !   Cu_Solver_Dn_Dgetrf   !
    !-------------------------!
    integer(c_int) function Cu_Solver_Dn_Dgetrf  &
     (cusolver_hndl,                             &
      m,                                         &
      n,                                         &
      pnt_a_gpu,                                 &
      lda,                                       &
      pnt_ws_gpu,                                &
      pnt_piv_gpu,                               &
      pnt_devinfo_gpu)                           &
      bind(C, name="cusolverDnDgetrf")

      use iso_c_binding
      implicit none

      type   (c_ptr), value :: cusolver_hndl   ! CUDA solver handle
      integer(c_int), value :: m
      integer(c_int), value :: n
      type   (c_ptr), value :: pnt_a_gpu
      integer(c_int), value :: lda
      type   (c_ptr), value :: pnt_ws_gpu      ! pointer to workspace on device
      type   (c_ptr), value :: pnt_piv_gpu     ! pivoting sequence on device
      type   (c_ptr), value :: pnt_devinfo_gpu
    end function

    !-------------------------!
    !   Cu_Solver_Dn_Dgetrs   !
    !-------------------------!
    integer(c_int) function Cu_Solver_Dn_Dgetrs  &
     (cusolver_hndl,                             &
      trans,                                     &
      n,                                         &
      nrhs,                                      &
      pnt_a_gpu,                                 &
      lda,                                       &
      pnt_piv_gpu,                               &
      pnt_b_gpu,                                 &
      ldb,                                       &
      pnt_devinfo_gpu)                           &
      bind(C, name="cusolverDnDgetrs")

      use iso_c_binding
      implicit none

      type   (c_ptr), value :: cusolver_hndl   ! CUDA solver handle
      integer(c_int), value :: trans
      integer(c_int), value :: n
      integer(c_int), value :: nrhs
      type   (c_ptr), value :: pnt_a_gpu
      integer(c_int), value :: lda
      type   (c_ptr), value :: pnt_piv_gpu     ! pivoting sequence on device
      type   (c_ptr), value :: pnt_b_gpu
      integer(c_int), value :: ldb
      type   (c_ptr), value :: pnt_devinfo_gpu
    end function

  end interface

  end module

