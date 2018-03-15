module imp_property
    use constants
    contains
    function get_density(tin) result(density)
     real(KREAL),intent(in)::tin
     real(KREAL)::density
     density=11096-1.3326*tin
    end function get_density
end module imp_property
