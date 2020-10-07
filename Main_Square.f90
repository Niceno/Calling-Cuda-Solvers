!==============================================================================!
  program main
!------------------------------------------------------------------------------!
!   The purpose of this routine is to provide for GPU based solution of        !
!   large dense systems of equations in legacy FORTRAN Applications.           !
!------------------------------------------------------------------------------!
  use iso_c_binding
  use Cuda_Solver_Mod
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  ! Matrix definition & host CPU storage variables
  integer(c_int)       :: m, n, lda, ldb
  real,    allocatable :: a(:,:)   ! square matrix
  real,    allocatable :: b(:)     ! right hand side vector
  real,    allocatable :: x(:)     ! solution vector
  integer, allocatable :: piv(:)   ! host copy of pivoting sequence
  integer, target      :: lwork    ! size of workspace

  ! CPU equivalents of device variables
  integer        :: devinfo
  integer(c_int) :: nrhs           ! number of right hand sides?

  ! Handle to device
  type(c_ptr) :: cusolver_hndl

  ! Pointers to device memory
  type(c_ptr) :: pnt_a_gpu
  type(c_ptr) :: pnt_b_gpu
  type(c_ptr) :: pnt_lwork_gpu
  type(c_ptr) :: pnt_ws_gpu
  type(c_ptr) :: pnt_piv_gpu
  type(c_ptr) :: pnt_devinfo_gpu

  ! Function result variables
  integer :: error

  ! Parts of the linear system
  target :: a, b, x, piv

  ! Variables holding free and total memory
  integer, target :: free, total

  !-----------------------------------------------!
  !   Check the amount of free and total memory   !
  !-----------------------------------------------!
  free  = 0
  total = 0

  error = Cuda_Mem_Get_Info(c_loc(free), c_loc(total))
  call Error_Check("Cuda_Mem_Get_Info", error)
  write (*, '(a, i12)') "  free mem: ", free
  write (*, '(a, i12)') " total mem: ", total

  lda  = 3
  ldb  = 3
  m    = lda
  n    = lda
  nrhs = 1

  allocate(a(lda,n))
  allocate(b(n))
  allocate(x(n))
  allocate(piv(lda))
! allocate(lwork)

  !--------------------------------------------!
  !   Example from cusolver_library.pdf        !
  !   "QR factorization dense linear solver"   !
  !                                            !
  !   Define [A] and [B]                       !
  !                                            !
  !           [A]         [x]  =    [B]        !
  !    | 1.0  2.0  3.0 | |1.0|    | 6.0|       !
  !    | 4.0  5.0  6.0 | |1.0| =  |15.0|       !
  !    | 2.0  1.0  1.0 | |1.0|    | 4.0|       !
  !--------------------------------------------!
  a(1,1) = 1.0;   a(1,2) = 2.0;   a(1,3) = 3.0
  a(2,1) = 4.0;   a(2,2) = 5.0;   a(2,3) = 6.0
  a(3,1) = 2.0;   a(3,2) = 1.0;   a(3,3) = 1.0

  b(1) =  6.0
  b(2) = 15.0
  b(3) =  4.0

  !-----------------------------------!
  !   Step 1: Create cudense handle   !
  !-----------------------------------!
  error = Cu_Solver_Dn_Create(cusolver_hndl)
  call Error_Check("Cu_Solver_Dn_Create", error)

  !------------------------------------!
  !   Step 2: copy a and b to device   !
  !------------------------------------!
  error = Cuda_Malloc(pnt_a_gpu, sizeof(a))
  call Error_Check("Cuda_Malloc 1", error)

  error = Cuda_Malloc(pnt_b_gpu, sizeof(b))
  call Error_Check("Cuda_Malloc 2", error)

  ! Also allocate space for other device based variables
  error = Cuda_Malloc(pnt_piv_gpu, sizeof(piv))
  call Error_Check("Cuda_Malloc 3", error)

  error = Cuda_Malloc(pnt_devinfo_gpu, sizeof(devinfo))
  call Error_Check("Cuda_Malloc 4", error)

  error = Cuda_Malloc(pnt_lwork_gpu, sizeof(lwork))
  call Error_Check("Cuda_Malloc 5", error)

  ! Copy a and b to device
  error = Cuda_Mem_Cpy(pnt_a_gpu,  &                 ! target
                       c_loc(a),   &                 ! source
                       sizeof(a),  &                 ! size
                       CUDA_MEM_CPY_HOST_TO_DEVICE)  ! operation
  call Error_Check("Cuda_Mem_Cpy 1", error)

  error = Cuda_Mem_Cpy(pnt_b_gpu,  &                 ! target
                       c_loc(b),   &                 ! source
                       sizeof(b),  &                 ! size
                       CUDA_MEM_CPY_HOST_TO_DEVICE)  ! operation
  call Error_Check("Cuda_Mem_Cpy 2", error)

  !---------------------------------------------------------------------------!
  !   Step 3: query working space of Dgetrf (and allocate memory on device)   !
  !---------------------------------------------------------------------------!
  lwork = 5
  error = Cu_Solver_Dn_Dgetrf_Buffer_Size(cusolver_hndl,  &
                                          m,              &
                                          n,              &
                                          pnt_a_gpu,      &
                                          lda,            &
                                          c_loc(lwork))
  call Error_Check("Cu_Solver_Dn_Dgetrf_Buffer_Size", error)

  write (*,*)
  write (*, '(a, i12)') " lwork: ", lwork
  write (*,*)

  error = Cuda_Malloc(pnt_ws_gpu, 4 * lwork)  ! why 4 * lwork?
  call Error_Check("Cuda_Malloc 6", error)

  !---------------------------------------------!
  !   Step 4: compute LU factorization of [a]   !
  !---------------------------------------------!
  error = Cu_Solver_Dn_Dgetrf(cusolver_hndl,    &
                              m,                &
                              n,                &
                              pnt_a_gpu,        &
                              lda,              &
                              pnt_ws_gpu,       &
                              pnt_piv_gpu,      &
                              pnt_devinfo_gpu)
  call Error_Check("Cu_Solver_Dn_Dgetrf", error)

  !-----------------------------------------------------------------!
  !   Step 5: compute solution vector [x] for right hand side [b]   !
  !-----------------------------------------------------------------!
  error = Cu_Solver_Dn_Dgetrs(cusolver_hndl,    &
                              CUBLAS_OP_N,      &
                              n,                &
                              nrhs,             &
                              pnt_a_gpu,        &
                              lda,              &
                              pnt_piv_gpu,      &
                              pnt_b_gpu,        &
                              ldb,              &
                              pnt_devinfo_gpu)
  call Error_Check("Cu_Solver_Dn_Dgetrs", error)

  !------------------------------------------------!
  !   Step 6: copy solution vector stored in [b]   !
  !           on device into [x] vector on host    !
  !------------------------------------------------!
  error = Cuda_Mem_Cpy(c_loc(x),  &                  ! target
                       pnt_b_gpu,  &                 ! source
                       sizeof(b),  &                 ! size
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
  error = Cuda_Free(pnt_piv_gpu)
  error = Cuda_Free(pnt_ws_gpu)
  error = Cuda_Free(pnt_lwork_gpu)
  error = Cu_Solver_Dn_Destroy(cusolver_hndl)

  !---------------------------------------------------!
  !   Step 8: deallocate memory on host before exit   !
  !---------------------------------------------------!
  deallocate(a)
  deallocate(b)
  deallocate(x)
  deallocate(piv)

  end program
