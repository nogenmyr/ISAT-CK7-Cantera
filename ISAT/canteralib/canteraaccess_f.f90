module canteraAccess
interface

    integer function nasadata(n,thermo_ns)
        integer, intent(in) :: n
        double precision, intent(inout) :: thermo_ns(*)
    end function nasadata

end interface
end module canteraAccess
