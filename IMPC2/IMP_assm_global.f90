module imp_assm_global
    use imp_assembly_header
    use imp_timer_header
	use constants
    implicit none
    type(sys_assembly),allocatable::assm1(:)
	real(KREAL),allocatable::imp_pow(:,:)
    !type(sys_timer)::timer1
end module imp_assm_global