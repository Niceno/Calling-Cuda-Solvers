!==============================================================================!
!   Interface to cusolverDn and CUDA C functions                               !
!                                                                              !
!   Compile with:                                                              !
!   nvfortran -i8 -L /usr/local/cuda/lib64 another.f90 \                       !
!             -lcudart -lcublas -lcusparse -lcusolver                          !
!------------------------------------------------------------------------------!
  include 'Cuda_Solver_Mod.f90'

!==============================================================================!
  subroutine Error_Check(call_name, err)
!------------------------------------------------------------------------------!
  implicit none
!------------------------------------------------------------------------------!
  character(len=*) :: call_name
  integer*8        :: err
!------------------------------------------------------------------------------!

  if(err .ne. 0) then
    write (*,*)
    write (*, '(a,a,i2)')  call_name, " error: ", err
    write (*,*)
    stop
  end if

  end subroutine

!==============================================================================!
  program main
!------------------------------------------------------------------------------!
!   The purpose of this routine is to provide for GPU based solution of        !
!   large dense systems of equations in legacy FORTRAN Applications.           !
!                                                                              !
!   This is the development version of the routine which does not yet work ... !
!------------------------------------------------------------------------------!
  use iso_c_binding
  use Cuda_Solver_Mod
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  ! Matrix definition & host CPU storage variables
  integer(c_int)       :: m, n, lda, ldb
  real, allocatable    :: A(:,:)
  real, allocatable    :: B(:)
  real, allocatable    :: X(:)
  integer, allocatable :: Ipiv(:)
  integer, target      :: Lwork

  ! CPU equivalents of device variables
  integer devInfo
  integer(c_int) nrhs

  integer*8 Ipiv_size, devInfo_size, Lwork_size
  integer*8 Workspace

  ! Handle to device
  type(c_ptr) :: cusolver_Hndl

  ! Pointers to device memory
  type(c_ptr) :: pnt_a_gpu
  type(c_ptr) :: pnt_b_gpu
  type(c_ptr) :: d_Lwork
  type(c_ptr) :: d_WS
  type(c_ptr) :: d_Ipiv
  type(c_ptr) :: d_devInfo

  ! Function result variables
  integer*8 :: error

  ! Pointers to host CPU memory
  type(c_ptr) :: pnt_a_cpu
  type(c_ptr) :: pnt_b_cpu
  type(c_ptr) :: pnt_x_cpu
  type(c_ptr) :: CPU_Lwork_ptr

  target :: A, B, X, Ipiv

  type(c_ptr)          :: cpfre,cptot
  integer*8, target    :: free,total
  integer*8, parameter :: CUDA_MEM_CPY_DEVICE_TO_HOST = 2
  integer*8, parameter :: CUDA_MEM_CPY_HOST_TO_DEVICE = 1
  integer*4, parameter :: CUBLAS_OP_N = 0
  integer*4, parameter :: CUBLAS_OP_T = 1

  ! ================================================
  free  = 0
  total = 0
  cpfre = c_loc(free)
  cptot = c_loc(total)

  error = Cuda_Mem_Get_Info(cpfre,cptot)
  call Error_Check("Cuda_Mem_Get_Info", error)
  write (*, '(A, I12)') "  free mem: ", free
  write (*, '(A, I12)') " total mem: ", total

  lda  = 3
  ldb  = 3
  m    = lda
  n    = lda
  nrhs = 1

  allocate(A(lda,n))
  allocate(B(n))
  allocate(X(n))
  allocate(Ipiv(lda))
!  allocate(Lwork)

  Ipiv_size    = sizeof(Ipiv)
  devInfo_size = sizeof(devInfo)
  Lwork_size   = sizeof(Lwork)

  !--------------------------------------------!
  !   Example from cusolver_library.pdf        !
  !   "QR factorization dense linear solver"   !
  !                                            !
  !   DeFINE [A] AND [B]                       !
  !                                            !
  !           [A]         [x]  =    [B]        !
  !    | 1.0  2.0  3.0 | |1.0|    | 6.0|       !
  !    | 4.0  5.0  6.0 | |1.0| =  |15.0|       !
  !    | 2.0  1.0  1.0 | |1.0|    | 4.0|       !
  !--------------------------------------------!
  A(1,1) = 1.0;   A(1,2) = 2.0;  A(1,3) = 3.0
  A(2,1) = 4.0;   A(2,2) = 5.0;  A(2,3) = 6.0
  A(3,1) = 2.0;   A(3,2) = 1.0;  A(3,3) = 1.0

  B(1) = 6.0
  B(2) = 15.0
  B(3) = 4.0

  pnt_a_cpu     = c_loc(A)
  pnt_b_cpu     = c_loc(B)
  pnt_x_cpu     = c_loc(X)
  CPU_Lwork_ptr = c_loc(Lwork)

  !-----------------------------------!
  !   Step 1: Create cudense handle   !
  !-----------------------------------!
  error = Cu_Solver_Dn_Create(cusolver_Hndl)
  call Error_Check("Cu_Solver_Dn_Create", error)

  !------------------------------------!
  !   Step 2: copy A and B to Device   !
  !------------------------------------!
  error = Cuda_Malloc(pnt_a_gpu, sizeof(A))
  call Error_Check("Cuda_Malloc 1", error)

  error = Cuda_Malloc(pnt_b_gpu, sizeof(B))
  call Error_Check("Cuda_Malloc 2", error)

  ! Also allocate space for other device based variables
  error = Cuda_Malloc(d_Ipiv, Ipiv_size)
  call Error_Check("Cuda_Malloc 3", error)

  error = Cuda_Malloc(d_devInfo, devInfo_size)
  call Error_Check("Cuda_Malloc 4", error)

  error = Cuda_Malloc(d_Lwork, Lwork_size)
  call Error_Check("Cuda_Malloc 5", error)

  ! Copy A and B to device
  error = Cuda_Mem_Cpy(pnt_a_gpu,  &                 ! target
                       pnt_a_cpu,  &                 ! source
                       sizeof(A),  &                 ! size
                       CUDA_MEM_CPY_HOST_TO_DEVICE)  ! operation
  call Error_Check("Cuda_Mem_Cpy 1", error)

  error = Cuda_Mem_Cpy(pnt_b_gpu,  &                 ! target
                       pnt_b_cpu,  &                 ! source
                       sizeof(B),  &                 ! size
                       CUDA_MEM_CPY_HOST_TO_DEVICE)  ! operation
  call Error_Check("Cuda_Mem_Cpy 2", error)

  !---------------------------------------------------------------------------!
  !   Step 3: query working space of Sgetrf (and allocate memory on device)   !
  !---------------------------------------------------------------------------!
  Lwork = 5
  error = Cu_Solver_Dn_Dgetrf_Buffer_Size(cusolver_Hndl, m, n, pnt_a_gpu, lda, CPU_Lwork_ptr)
  call Error_Check("Cu_Solver_Dn_Dgetrf_Buffer_Size", error)

  write (*,*)
  write (*, '(A, I12)') " Lwork: ", Lwork
  write (*,*)

  Workspace = 4 * Lwork
  error = Cuda_Malloc(d_WS, Workspace)
  call Error_Check("Cuda_Malloc 6", error)

  !---------------------------------------------!
  !   Step 4: compute LU factorization of [A]   !
  !---------------------------------------------!
  error = Cu_Solver_Dn_Dgetrf(cusolver_Hndl, m, n, pnt_a_gpu, lda, d_WS, d_Ipiv, d_devInfo) 
  call Error_Check("Cu_Solver_Dn_Dgetrf", error)

  !-----------------------------------------------------------------!
  !   Step 5: compute solution vector [X] for right hand side [B]   !
  !-----------------------------------------------------------------!
  error = Cu_Solver_Dn_Dgetrs(cusolver_Hndl, CUBLAS_OP_N, n, nrhs, pnt_a_gpu, lda, d_Ipiv, pnt_b_gpu, ldb, d_devInfo)
  call Error_Check("Cu_Solver_Dn_Dgetrs", error)

  !------------------------------------------------!
  !   Step 6: copy solution vector stored in [B]   !
  !           on device into [X] vector on host    !
  !------------------------------------------------!
  error = Cuda_Mem_Cpy(pnt_x_cpu,  &                 ! target
                       pnt_b_gpu,  &                 ! source
                       sizeof(B),  &                 ! size
                       CUDA_MEM_CPY_DEVICE_TO_HOST)  ! operation
  call Error_Check("Cuda_Mem_Cpy 4", error)

  print *, 'Solution vector:'
  print *, x(1)
  print *, x(2)
  print *, x(3)

  !------------------------------------------------------------------!
  !   Step 7: free memory on device and release CPU-side resources   !
  !------------------------------------------------------------------!
  error = Cuda_Free(pnt_a_gpu)
  error = Cuda_Free(pnt_b_gpu)
  error = Cuda_Free(d_Ipiv)
  error = Cuda_Free(d_WS)
  error = Cuda_Free(d_Lwork)

  error = Cu_Solver_Dn_Destroy(cusolver_Hndl)

  !---------------------------------------------------!
  !   Step 8: deallocate memory on host before exit   !
  !---------------------------------------------------!
  deallocate(A)
  deallocate(B)
  deallocate(X)
  deallocate(Ipiv)

  end program
