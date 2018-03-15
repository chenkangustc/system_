!$
!===================================================================================================
!
!   driver for input parameter from xml format file
!---------------------------------------------------------------------------------------------------
!   Public subroutine lists:    Drving_xml_scan
!                               Drving_xml_read
!
!   Public type lists:          No
!
!===================================================================================================
module input_xml
    
    use constants
    use global_state
    use, intrinsic  :: ISO_FORTRAN_ENV
    
    use global
    use th_global
    
    use input_keyword
    use xml_interface
  
    implicit none 
    private
    public  :: Driving_xml_scan, Driving_xml_read
    
contains
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Driving_xml_scan (filename)
        
        character(len=*), intent(in)  :: filename
        
!!!        ! local variables
!!!        character(len=MAX_WORD_LEN)  :: a_line
!!!        character(len=MAX_WORD_LEN)  :: section_name
!!!        character(len=MAX_WORD_LEN)  :: keyword
!!!        
!!!        type(Node), pointer  :: mydoc => NULL()
!!!        type(Node), pointer  :: node_section => NuLL ()
!!!        type(Node), pointer  :: node_keyword => NULL ()
!!!        type(Node), pointer  :: node_data => NULL ()
!!!        
!!!        ! open file for scan
!!!        call open_xmldoc (mydoc, filename)
!!!        
!!!        if (check_for_node(mydoc, 'CASENAME:'))  then
!!!            call get_node_ptr(mydoc, 'CASENAME:', node_section)
!!!        end if
!!!        
!!!        if (check_for_node(mydoc, 'CONTROL:'))  then
!!!            call get_node_ptr(mydoc, 'CONTROL:', node_section)
!!!        end if
!!!        
!!!        if (check_for_node(mydoc, 'METHOD:'))  then
!!!            call get_node_ptr(mydoc, 'METHOD:', node_section)
!!!        end if
!!!        
!!!        if (check_for_node(mydoc, 'LINK:'))  then
!!!            call get_node_ptr(mydoc, 'LINK:', node_section)
!!!        end if
!!!        
!!!        if (check_for_node(mydoc, 'MATERIAL:'))  then
!!!            call get_node_ptr(mydoc, 'MATERIAL:', node_section)
!!!        end if
!!!        
!!!        if (check_for_node(mydoc, 'GEOMETRY:'))  then
!!!            call get_node_ptr(mydoc, 'GEOMETRY:', node_section)
!!!        end if
!!!        
!!!        if (check_for_node(mydoc, 'DEPLETION:'))  then
!!!            call get_node_ptr(mydoc, 'DEPLETION:', node_section)
!!!        end if
!!!        
!!!        if (check_for_node(mydoc, 'FEEDBACK:'))  then
!!!            call get_node_ptr(mydoc, 'FEEDBACK:', node_section)
!!!        end if
!!!        
!!!        if (check_for_node(mydoc, 'PERTURBATION:'))  then
!!!            call get_node_ptr(mydoc, 'PERTURBATION:', node_section)
!!!        end if
!!!        
!!!        if (check_for_node(mydoc, 'END:'))  then
!!!            call get_node_ptr(mydoc, 'END:', node_section)
!!!        end if
!!!        
!!!        ! close xml file
!!!        call close_xmldoc (mydoc)
        
    end subroutine Driving_xml_scan
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Driving_xml_read (filename)
        
        character(len=*), intent(in)  :: filename
        
!!!        ! local variables
!!!        character(len=MAX_WORD_LEN)  :: a_line
!!!        character(len=MAX_WORD_LEN)  :: section_name
!!!        character(len=MAX_WORD_LEN)  :: keyword
!!!        
!!!        type(Node), pointer  :: mydoc => NULL()
!!!        type(Node), pointer  :: node_section => NuLL ()
!!!        type(Node), pointer  :: node_keyword => NULL ()
!!!        type(Node), pointer  :: node_data => NULL ()
!!!        
!!!        ! open file for scan
!!!        call open_xmldoc (mydoc, filename)
!!!        
!!!        if (check_for_node(mydoc, 'CASENAME:'))  then
!!!            call get_node_ptr(mydoc, 'CASENAME:', node_section)
!!!        end if
!!!        
!!!        if (check_for_node(mydoc, 'CONTROL:'))  then
!!!            call get_node_ptr(mydoc, 'CONTROL:', node_section)
!!!        end if
!!!        
!!!        if (check_for_node(mydoc, 'METHOD:'))  then
!!!            call get_node_ptr(mydoc, 'METHOD:', node_section)
!!!        end if
!!!        
!!!        if (check_for_node(mydoc, 'LINK:'))  then
!!!            call get_node_ptr(mydoc, 'LINK:', node_section)
!!!        end if
!!!        
!!!        if (check_for_node(mydoc, 'MATERIAL:'))  then
!!!            call get_node_ptr(mydoc, 'MATERIAL:', node_section)
!!!        end if
!!!        
!!!        if (check_for_node(mydoc, 'GEOMETRY:'))  then
!!!            call get_node_ptr(mydoc, 'GEOMETRY:', node_section)
!!!        end if
!!!        
!!!        if (check_for_node(mydoc, 'DEPLETION:'))  then
!!!            call get_node_ptr(mydoc, 'DEPLETION:', node_section)
!!!        end if
!!!        
!!!        if (check_for_node(mydoc, 'FEEDBACK:'))  then
!!!            call get_node_ptr(mydoc, 'FEEDBACK:', node_section)
!!!        end if
!!!        
!!!        if (check_for_node(mydoc, 'PERTURBATION:'))  then
!!!            call get_node_ptr(mydoc, 'PERTURBATION:', node_section)
!!!        end if
!!!        
!!!        if (check_for_node(mydoc, 'END:'))  then
!!!            call get_node_ptr(mydoc, 'END:', node_section)
!!!        end if
!!!        
!!!        ! close xml file
!!!        call close_xmldoc (mydoc)
!!!        
    end subroutine Driving_xml_read
    
end module input_xml
