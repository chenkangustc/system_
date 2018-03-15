!$
!===================================================================================================
!
!   interface for FoX xml lib parsing by DoM method
!---------------------------------------------------------------------------------------------------
!   Public subroutine lists:    getChildrenByTagName
!
!   Public type lists:          xxx
!
!===================================================================================================
module xml_wrapper
    
!!!    use fox_m_fsys_array_str,   only: str_vs, vs_str, vs_str_alloc
!!!    use fox_m_fsys_format,      only: operator(//)
!!!    use fox_m_fsys_string,      only: toLower
!!!    use fox_m_utils_uri,        only: URI, parseURI, destroyURI, isAbsoluteURI, &
!!!                                &   rebaseURI, expressURI
!!!    use m_common_charset,       only: checkChars, XML1_0, XML1_1
!!!    use m_common_element,       only: element_t, get_element, attribute_t, &
!!!                                &   attribute_has_default, get_attribute_declaration, get_attlist_size
!!!    use m_common_namecheck,     only: checkQName, prefixOfQName, localPartOfQName, &
!!!                                &   checkName, checkPublicId, checkNCName
!!!    use m_common_struct,        only: xml_doc_state, init_xml_doc_state, destroy_xml_doc_state
!!!
!!!    use m_dom_error,            only: DOMException, throw_exception, inException, getExceptionCode, &
!!!                                &   NO_MODIFICATION_ALLOWED_ERR, NOT_FOUND_ERR, HIERARCHY_REQUEST_ERR, &
!!!                                &   WRONG_DOCUMENT_ERR, FoX_INTERNAL_ERROR, FoX_NODE_IS_NULL, FoX_LIST_IS_NULL, &
!!!                                &   INUSE_ATTRIBUTE_ERR, FoX_MAP_IS_NULL, INVALID_CHARACTER_ERR, NAMESPACE_ERR, &
!!!                                &   FoX_INVALID_PUBLIC_ID, FoX_INVALID_SYSTEM_ID, FoX_IMPL_IS_NULL, FoX_INVALID_NODE, &
!!!                                &   FoX_INVALID_CHARACTER, FoX_INVALID_COMMENT, FoX_INVALID_CDATA_SECTION, &
!!!                                &   FoX_INVALID_PI_DATA, NOT_SUPPORTED_ERR, FoX_INVALID_ENTITY, &
!!!                                &   INDEX_SIZE_ERR, FoX_NO_SUCH_ENTITY, FoX_HIERARCHY_REQUEST_ERR, &
!!!                                &   FoX_INVALID_URI
!!!
!!!    use m_dom_dom
!!!    use fox_dom
!!!    use fox_m_fsys_count_parse_input, only: countrts
!!!
!!!    implicit none
!!!
!!!contains
!!!    !$
!!!    !===============================================================================================
!!!    ! 
!!!    !===============================================================================================
!!!    function getChildrenByTagName (doc, tagName, name, ex) result(list) 
!!!    
!!!        type(Node), pointer                         :: doc
!!!        character(len=*), intent(in), optional      :: tagName, name
!!!        type(DOMException), intent(out), optional   :: ex
!!!        type(NodeList), pointer :: list
!!!        
!!!        type(NodeListPtr), pointer :: nll(:), temp_nll(:)
!!!        type(Node), pointer :: arg, this, treeroot
!!!        logical :: doneChildren, doneAttributes, allElements
!!!        integer :: i, i_tree
!!!        
!!!        list => null()
!!!        
!!!        if (.NOT. associated(doc)) then
!!!            if (getFoX_checks() .OR. FoX_NODE_IS_NULL<200) then
!!!                call throw_exception(FoX_NODE_IS_NULL, "getElementsByTagName", ex)
!!!                if (PRESENT(ex)) then
!!!                    if (inException(ex)) then
!!!                        return
!!!                    end if 
!!!                end if 
!!!            end if 
!!!        end if 
!!!        
!!!        if (doc%nodeType==DOCUMENT_NODE) then
!!!            if (PRESENT(name) .OR. .NOT. PRESENT(tagName)) then
!!!                if (getFoX_checks() .OR. FoX_INVALID_NODE<200) then
!!!                    call throw_exception(FoX_INVALID_NODE, "getElementsByTagName", ex)
!!!                    if (PRESENT(ex)) then
!!!                        if (inException(ex)) then
!!!                            return
!!!                        end if 
!!!                    end if 
!!!                end if 
!!!            end if 
!!!        else if  (doc%nodeType==ELEMENT_NODE) then
!!!            if (PRESENT(name) .OR. .NOT. PRESENT(tagName)) then
!!!                if (getFoX_checks() .OR. FoX_INVALID_NODE<200) then
!!!                    call throw_exception(FoX_INVALID_NODE, "getElementsByTagName", ex)
!!!                    if (PRESENT(ex)) then
!!!                        if (inException(ex)) then
!!!                            return
!!!                        end if 
!!!                    end if 
!!!                end if 
!!!            end if 
!!!        else      
!!!            if (getFoX_checks() .OR. FoX_INVALID_NODE<200) then
!!!                call throw_exception(FoX_INVALID_NODE, "getElementsByTagName", ex)
!!!                if (PRESENT(ex)) then
!!!                    if (inException(ex)) then
!!!                        return
!!!                    end if 
!!!                end if 
!!!            end if 
!!!        end if 
!!!        
!!!        if (doc%nodeType==DOCUMENT_NODE) then
!!!            arg => getDocumentElement(doc)
!!!        else
!!!            arg => doc
!!!        end if 
!!!        
!!!        allocate(list)
!!!        allocate(list%nodes(0))
!!!        list%element => doc
!!!        if (PRESENT(name)) list%nodeName => vs_str_alloc(name)
!!!        if (PRESENT(tagName)) list%nodeName => vs_str_alloc(tagName)
!!!        
!!!        allElements = (str_vs(list%nodeName)=="*")
!!!        
!!!        if (doc%nodeType==DOCUMENT_NODE) then
!!!            nll => doc%docExtras%nodelists
!!!        else if  (doc%nodeType==ELEMENT_NODE) then
!!!            nll => doc%ownerDocument%docExtras%nodelists
!!!        end if 
!!!        allocate(temp_nll(size(nll)+1))
!!!        do i = 1, size(nll)
!!!            temp_nll(i)%this => nll(i)%this
!!!        end do 
!!!        temp_nll(i)%this => list
!!!        deallocate(nll)
!!!        if (doc%nodeType==DOCUMENT_NODE) then
!!!            doc%docExtras%nodelists => temp_nll
!!!        else if  (doc%nodeType==ELEMENT_NODE) then
!!!            doc%ownerDocument%docExtras%nodelists => temp_nll
!!!        end if 
!!!        
!!!        treeroot => arg
!!!        i_tree = 0
!!!        doneChildren = .FALSE.
!!!        doneAttributes = .FALSE.
!!!        this => treeroot
!!!        do
!!!            if (.NOT. doneChildren .AND. .NOT. (getNodeType(this)==ELEMENT_NODE .AND. doneAttributes)) then
!!!                if (this%nodeType==ELEMENT_NODE) then
!!!                    if ((allElements .OR. str_vs(this%nodeName)==tagName) &
!!!                    .AND. .NOT.(getNodeType(doc)==ELEMENT_NODE .AND. associated(this, arg))) &
!!!                    call append(list, this)
!!!                    doneAttributes = .TRUE.
!!!                end if 
!!!            else
!!!                if (getNodeType(this)==ELEMENT_NODE .AND. .NOT. doneChildren) then
!!!                    doneAttributes = .TRUE.
!!!                else
!!!                end if 
!!!            end if 
!!!        
!!!        
!!!            if (.NOT. doneChildren) then
!!!                if (getNodeType(this)==ELEMENT_NODE .AND. .NOT. doneAttributes) then
!!!                    if (getLength(getAttributes(this))>0) then
!!!                        this => item(getAttributes(this), 0)
!!!                    else
!!!                        doneAttributes = .TRUE.
!!!                    end if 
!!!                else if  (hasChildNodes(this) .AND. .NOT. associated(getParentNode(this), treeroot)) then
!!!                    this => getFirstChild(this)
!!!                    doneChildren = .FALSE.
!!!                    doneAttributes = .FALSE.
!!!                else
!!!                    doneChildren = .TRUE.
!!!                    doneAttributes = .FALSE.
!!!                end if 
!!!                
!!!            else ! if doneChildren
!!!                if (associated(this, treeroot)) exit
!!!                if (getNodeType(this)==ATTRIBUTE_NODE) then
!!!                    if (i_tree<getLength(getAttributes(getOwnerElement(this)))-1) then
!!!                        i_tree= i_tree+ 1
!!!                        this => item(getAttributes(getOwnerElement(this)), i_tree)
!!!                        doneChildren = .FALSE.
!!!                    else
!!!                        i_tree= 0
!!!                        this => getOwnerElement(this)
!!!                        doneAttributes = .TRUE.
!!!                        doneChildren = .FALSE.
!!!                    end if 
!!!                else if  (associated(getNextSibling(this))) then
!!!                    this => getNextSibling(this)
!!!                    doneChildren = .FALSE.
!!!                    doneAttributes = .FALSE.
!!!                else
!!!                    this => getParentNode(this)
!!!                end if 
!!!            end if 
!!!        end do 
!!!
!!!    end function getChildrenByTagName
    
end module xml_wrapper
