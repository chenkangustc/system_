module imp_re_input_global
    !use sys_reactor_header !与全局变量冲突
    use imp_re_input_header
    use imp_assm_global
    implicit none
    type(sys_re_input)::reInputdata
    !type(sys_reactor)::reactor
end module imp_re_input_global