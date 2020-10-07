!==============================================================================!
!   Interface to cusolverDn and CUDA C functions                               !
!                                                                              !
!   Compile with:                                                              !
!   nvfortran -i8 -L /usr/local/cuda/lib64 another.f90 \                       !
!             -lcudart -lcublas -lcusparse -lcusolver                          !
!------------------------------------------------------------------------------!
  module Cuda_Solver_Mod

  interface

    !-----------------!
    !   Cuda_Malloc   !
    !-----------------!
    integer(c_int) function Cuda_Malloc  &
      (buffer, size)                      &
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
      (dst, src, count, kind)              &
      bind (C, name="cudaMemcpy" )
      ! note: CUDA_MEM_CPY_HOST_TO_DEVICE = 1
      ! note: CUDA_MEM_CPY_DEVICE_TO_HOST = 2

      use iso_c_binding
      implicit none
      type(c_ptr),       value :: dst, src
      integer(c_size_t), value :: count, kind
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
      (cusolver_Hndl)                             &
      bind(C, name = "cusolverDnCreate")

      use iso_c_binding
      implicit none

      type(c_ptr) :: cusolver_Hndl
    end function

    !--------------------------!
    !   Cu_Solver_Dn_Destroy   !
    !--------------------------!
    integer(c_int) function Cu_Solver_Dn_Destroy  &
      (cusolver_Hndl)                             &
      bind(C, name = "cusolverDnDestroy")

      use iso_c_binding
      implicit none

      type(c_ptr), value :: cusolver_Hndl
    end function

    !-------------------------------------!
    !   Cu_Solver_Dn_Sgetrf_Buffer_Size   !
    !-------------------------------------!
    integer(c_int) function Cu_Solver_Dn_Sgetrf_Buffer_Size  &
      (cusolver_Hndl, m, n, pnt_a_gpu, lda, Lwork)                 &
      bind(C,name="cusolverDnSgetrf_bufferSize") 

      use iso_c_binding
      implicit none

      type   (c_ptr), value :: cusolver_Hndl
      integer(c_int), value :: m
      integer(c_int), value :: n
      type   (c_ptr), value :: pnt_a_gpu
      integer(c_int), value :: lda
      type   (c_ptr), value :: Lwork
    end function

    !-------------------------!
    !   Cu_Solver_Dn_Sgetrf   !
    !-------------------------!
    integer(c_int) function Cu_Solver_Dn_Sgetrf                    &
      (cusolver_Hndl, m, n, pnt_a_gpu, lda, d_WS, d_Ipiv, d_devInfo)  &
      bind(C, name="cusolverDnSgetrf")

      use iso_c_binding
      implicit none

      type   (c_ptr), value :: cusolver_Hndl
      integer(c_int), value :: m
      integer(c_int), value :: n
      type   (c_ptr), value :: pnt_a_gpu
      integer(c_int), value :: lda
      type   (c_ptr), value :: d_WS
      type   (c_ptr), value :: d_Ipiv
      type   (c_ptr), value :: d_devInfo
    end function

    !-------------------------!
    !   Cu_Solver_Dn_Sgetrs   !
    !-------------------------!
    integer(c_int) function Cu_Solver_Dn_Sgetrs                               &
      (cusolver_Hndl, trans, n, nrhs, pnt_a_gpu, lda, d_Ipiv, pnt_b_gpu, ldb, d_devInfo)  &
      bind(C, name="cusolverDnSgetrs")

      use iso_c_binding
      implicit none

      type   (c_ptr), value :: cusolver_Hndl
      integer(c_int), value :: trans
      integer(c_int), value :: n
      integer(c_int), value :: nrhs
      type   (c_ptr), value :: pnt_a_gpu
      integer(c_int), value :: lda
      type   (c_ptr), value :: d_Ipiv
      type   (c_ptr), value :: pnt_b_gpu
      integer(c_int), value :: ldb
      type   (c_ptr), value :: d_devInfo
    end function

    !-------------------------------------!
    !   Cu_Solver_Dn_Dgetrf_Buffer_Size   !
    !-------------------------------------!
    integer(c_int) function Cu_Solver_Dn_Dgetrf_Buffer_Size  &
      (cusolver_Hndl, m, n, pnt_a_gpu, lda, Lwork)                 &
      bind(C,name="cusolverDnDgetrf_bufferSize") 

      use iso_c_binding
      implicit none

      type   (c_ptr), value :: cusolver_Hndl
      integer(c_int), value :: m
      integer(c_int), value :: n
      type   (c_ptr), value :: pnt_a_gpu
      integer(c_int), value :: lda
      type   (c_ptr), value :: Lwork
    end function

    !-------------------------!
    !   Cu_Solver_Dn_Dgetrf   !
    !-------------------------!
    integer(c_int) function Cu_Solver_Dn_Dgetrf                    &
      (cusolver_Hndl, m, n, pnt_a_gpu, lda, d_WS, d_Ipiv, d_devInfo)  &
      bind(C, name="cusolverDnDgetrf")

      use iso_c_binding
      implicit none

      type   (c_ptr), value :: cusolver_Hndl
      integer(c_int), value :: m
      integer(c_int), value :: n
      type   (c_ptr), value :: pnt_a_gpu
      integer(c_int), value :: lda
      type   (c_ptr), value :: d_WS
      type   (c_ptr), value :: d_Ipiv
      type   (c_ptr), value :: d_devInfo
    end function

    !-------------------------!
    !   Cu_Solver_Dn_Dgetrs   !
    !-------------------------!
    integer(c_int) function Cu_Solver_Dn_Dgetrs                               &
      (cusolver_Hndl, trans, n, nrhs, pnt_a_gpu, lda, d_Ipiv, pnt_b_gpu, ldb, d_devInfo)  &
      bind(C, name="cusolverDnDgetrs")

      use iso_c_binding
      implicit none

      type   (c_ptr), value :: cusolver_Hndl
      integer(c_int), value :: trans
      integer(c_int), value :: n
      integer(c_int), value :: nrhs
      type   (c_ptr), value :: pnt_a_gpu
      integer(c_int), value :: lda
      type   (c_ptr), value :: d_Ipiv
      type   (c_ptr), value :: pnt_b_gpu
      integer(c_int), value :: ldb
      type   (c_ptr), value :: d_devInfo
    end function

  end interface

  end module

