module m_sax_tokenizer

#ifndef DUMMYLIB
  use fox_m_fsys_array_str, only: vs_str, str_vs, vs_str_alloc
  use m_common_charset, only: XML_WHITESPACE, &
    upperCase, isInitialNameChar
  use m_common_error, only: add_error, in_error
  use m_common_entities, only: entity_t, existing_entity, &
    expand_entity_text, expand_char_entity, &
    add_internal_entity, pop_entity_list, getEntityByName
  use m_common_namecheck, only: checkName, checkCharacterEntityReference

  use m_sax_reader, only: file_buffer_t, open_new_file, pop_buffer_stack, &
    push_chars, get_character, get_all_characters, &
    parse_text_declaration
  use m_sax_types ! everything, really
#ifdef PGF90
  use fox_m_utils_uri, only: URI
#endif

  implicit none
  private

  public :: sax_tokenize
  public :: normalize_attribute_text
  public :: expand_pe_text

contains

  subroutine sax_tokenize(fx, fb, eof)
    type(sax_parser_t), intent(inout) :: fx
    type(file_buffer_t), intent(inout) :: fb
    logical, intent(out) :: eof

    character :: c, q
    integer :: xv, phrase
    logical :: firstChar, ws_discard
    character, pointer :: tempString(:)

    xv = fx%xds%xml_version

    if (associated(fx%token)) deallocate(fx%token)
    fx%token => vs_str_alloc("")
    if (fx%nextTokenType/=TOK_NULL) then
      eof = .false.
      fx%tokenType = fx%nextTokenType
      fx%nextTokenType = TOK_NULL
      return
    endif
    fx%tokentype = TOK_NULL

    q = " "
    phrase = 0
    firstChar = .true.
    ws_discard = .false.
    do
      c = get_character(fb, eof, fx%error_stack)
      if (eof) then
        if (fx%state==ST_CHAR_IN_CONTENT) then
          if (phrase==1) then
            tempString => fx%token
            fx%token => vs_str_alloc(str_vs(fx%token)//"]")
            deallocate(tempString)
          elseif (phrase==2) then
            tempString => fx%token
            fx%token => vs_str_alloc(str_vs(fx%token)//"]]")
            deallocate(tempString)
          endif
          fx%tokenType = TOK_CHAR
        endif
        if (fx%tokenType/=TOK_NULL) then
          ! make sure we pass back this token before eof'ing
          eof = .false.
          return
        endif
      endif
      if (eof.or.in_error(fx%error_stack)) return
      if (fx%inIntSubset) then
        tempString => fx%xds%intSubset
        fx%xds%intSubset => vs_str_alloc(str_vs(tempString)//c)
        deallocate(tempString)
      endif

      select case (fx%state)
      case (ST_MISC)
        if (firstChar) ws_discard = .true.
        if (ws_discard) then
          if (verify(c, XML_WHITESPACE)>0) then
            if (c=="<") then
              ws_discard = .false.
            else
              call add_error(fx%error_stack, "Unexpected character found outside content")
            endif
          endif
        else
          if (c=="?") then
            fx%tokenType = TOK_PI_TAG
          elseif (c=="!") then
            fx%tokenType = TOK_BANG_TAG
          elseif (isInitialNameChar(c, xv)) then
            call push_chars(fb, c)
            fx%tokenType = TOK_OPEN_TAG
          else
            call add_error(fx%error_stack, "Unexpected character after <")
          endif
        endif

      case (ST_BANG_TAG)
        if (firstChar) then
          if (c=="-") then
            phrase = 1
          elseif (c=="[") then
            fx%tokenType = TOK_OPEN_SB
          elseif (verify(c,upperCase)==0) then
            deallocate(fx%token)
            fx%token => vs_str_alloc(c)
          else
            call add_error(fx%error_stack, "Unexpected character after <!")
          endif
        elseif (phrase==1) then
          if (c=="-") then
            fx%tokenType = TOK_OPEN_COMMENT
          else
            call add_error(fx%error_stack, "Unexpected character after <!-")
          endif
        elseif (verify(c,XML_WHITESPACE)>0) then
          tempString => fx%token
          fx%token => vs_str_alloc(str_vs(tempString)//c)
          deallocate(tempString)
        else
          call push_chars(fb, c)
          fx%tokenType = TOK_NAME
        endif

      case (ST_START_PI)
        ! grab until whitespace or ?
        if (verify(c, XML_WHITESPACE//"?")>0) then
          tempString => fx%token
          fx%token => vs_str_alloc(str_vs(tempString)//c)
          deallocate(tempString)
        else
          fx%tokenType = TOK_NAME
          if (c=="?") call push_chars(fb, c)
        endif

      case (ST_PI_CONTENTS)
        if (firstChar) ws_discard = .true.
        if (ws_discard) then
          if (verify(c, XML_WHITESPACE)/=0) then
            ws_discard = .false.
          else
            cycle
          endif
        endif
        if (phrase==1) then
          if (c==">") then
            fx%tokenType = TOK_CHAR
            fx%nextTokenType = TOK_PI_END
          elseif (c=="?") then
            ! The last ? didn't mean anything, but this one might.
            tempString => fx%token
            fx%token => vs_str_alloc(str_vs(tempString)//"?")
            deallocate(tempString)
          else
            phrase = 0
            tempString => fx%token
            fx%token => vs_str_alloc(str_vs(tempString)//"?"//c)
            deallocate(tempString)
          endif
        elseif (c=="?") then
          phrase = 1
        else
          tempString => fx%token
          fx%token => vs_str_alloc(str_vs(tempString)//c)
          deallocate(tempString)
        endif

      case (ST_START_COMMENT)
        select case(phrase)
        case (0)
          if (c=="-") then
            phrase = 1
          else
          tempString => fx%token
          fx%token => vs_str_alloc(str_vs(tempString)//c)
          deallocate(tempString)
          endif
        case (1)
          if (c=="-") then
            phrase = 2
          else
          tempString => fx%token
          fx%token => vs_str_alloc(str_vs(tempString)//"-"//c)
          deallocate(tempString)
            phrase = 0
          endif
        case (2)
          if (c==">") then
            fx%tokenType = TOK_CHAR
            fx%nextTokenType = TOK_COMMENT_END
            exit
          else
            call add_error(fx%error_stack, &
              "Expecting > after -- inside a comment.")
          endif
        end select

      case (ST_START_TAG)
        ! grab until whitespace or /, >
        if (verify(c, XML_WHITESPACE//"/>")>0) then
          tempString => fx%token
          fx%token => vs_str_alloc(str_vs(tempString)//c)
          deallocate(tempString)
        else
          fx%tokenType = TOK_NAME
          if (c==">") then
            fx%nextTokenType = TOK_END_TAG
          else
            call push_chars(fb, c)
          endif
        endif

      case (ST_START_CDATA_DECLARATION)
        if (firstChar) then
          if (verify(c, XML_WHITESPACE)==0) then
            call add_error(fx%error_stack, &
              "Whitespace not allowed around CDATA in section declaration")
          else
            deallocate(fx%token)
            fx%token => vs_str_alloc(c)
            ws_discard = .false.
          endif
        else
          if (verify(c, XML_WHITESPACE)==0) then
            call add_error(fx%error_stack, &
              "Whitespace not allowed around CDATA in section declaration")
          elseif (c=="[") then
            fx%tokenType = TOK_NAME
            if (c=="[") fx%nextTokenType = TOK_OPEN_SB
          else
            tempString => fx%token
            fx%token => vs_str_alloc(str_vs(tempString)//c)
            deallocate(tempString)
          endif
        endif


      case (ST_CDATA_CONTENTS)
        select case(phrase)
        case (0)
          if (c=="]") then
            phrase = 1
          else
            tempString => fx%token
            fx%token => vs_str_alloc(str_vs(tempString)//c)
            deallocate(tempString)
          endif
        case (1)
          if (c=="]") then
            phrase = 2
          else
            tempString => fx%token
            fx%token => vs_str_alloc(str_vs(tempString)//"]"//c)
            deallocate(tempString)
            phrase = 0
          endif
        case (2)
          if (c==">") then
            fx%tokenType = TOK_CHAR
            fx%nextTokenType = TOK_SECTION_END
          elseif (c=="]") then
            tempString => fx%token
            fx%token => vs_str_alloc(str_vs(tempString)//"]")
            deallocate(tempString)
          else
            tempString => fx%token
            fx%token => vs_str_alloc(str_vs(tempString)//"]]"//c)
            deallocate(tempString)
            phrase = 0
          endif
        end select

      case (ST_IN_TAG)
        if (firstChar) then
          if (verify(c,XML_WHITESPACE//"/>")>0) then
            call add_error(fx%error_stack, &
              "Whitespace required inside tag")
            exit
          endif
          ws_discard = .true.
        endif
        if (ws_discard) then
          if (verify(c,XML_WHITESPACE)>0) then
            if (c==">") then
              fx%tokenType = TOK_END_TAG
            elseif (c=="/") then
              phrase = 1
              ws_discard = .false.
            else
              deallocate(fx%token)
              fx%token => vs_str_alloc(c)
              ws_discard = .false.
            endif
          endif
        else
          if (phrase==1) then
            if (c==">") then
              fx%tokenType = TOK_END_TAG_CLOSE
            else
              call add_error(fx%error_stack, &
                "Unexpected character after / in element tag")
              exit
            endif
          else
            if (verify(c,XML_WHITESPACE//"=/>")==0) then
              fx%tokenType = TOK_NAME
              if (c=="=") then
                fx%nextTokenType = TOK_EQUALS
              elseif (c==">") then
                fx%nextTokenType = TOK_END_TAG
              else
                call push_chars(fb, c)
              endif
            else
              tempString => fx%token
              fx%token => vs_str_alloc(str_vs(tempString)//c)
              deallocate(tempString)
            endif
          endif
        endif

      case (ST_ATT_NAME)
        if (firstChar) ws_discard = .true.
        if (ws_discard) then
          if (verify(c,XML_WHITESPACE)>0) then
            if (c=="=") then
              fx%tokenType = TOK_EQUALS
            else
              call add_error(fx%error_stack, &
                "Unexpected character in element tag, expected =")
            endif
          endif
        endif

      case (ST_ATT_EQUALS)
        if (firstChar) ws_discard = .true.
        if (ws_discard) then
          if (verify(c,XML_WHITESPACE)>0) then
            if (verify(c,"'""")==0) then
              q = c
              ws_discard = .false.
            else
              call add_error(fx%error_stack, "Expecting "" or '")
            endif
          endif
        else
          if (c==q) then
            fx%tokenType = TOK_CHAR
          else
            tempString => fx%token
            fx%token => vs_str_alloc(str_vs(tempString)//c)
            deallocate(tempString)
          endif
        endif

      case (ST_CHAR_IN_CONTENT)
        if (c=="<".or.c=="&") then
          if (phrase==1) then
            tempString => fx%token
            fx%token => vs_str_alloc(str_vs(tempString)//"]")
            deallocate(tempString)
          elseif (phrase==2) then
            tempString => fx%token
            fx%token => vs_str_alloc(str_vs(tempString)//"]]")
            deallocate(tempString)
          endif
          fx%tokenType = TOK_CHAR
          if (c=="<") then
            call push_chars(fb, c)
          elseif (c=="&") then
            fx%nextTokenType = TOK_ENTITY
          endif
        elseif (c=="]") then
          if (phrase==0) then
            phrase = 1
          elseif (phrase==1) then
            phrase = 2
          else
            tempString => fx%token
            fx%token => vs_str_alloc(str_vs(tempString)//"]")
            deallocate(tempString)
          endif
        elseif (c==">") then
          if (phrase==1) then
            phrase = 0
            tempString => fx%token
            fx%token => vs_str_alloc(str_vs(tempString)//"]>")
            deallocate(tempString)
          elseif (phrase==2) then
            call add_error(fx%error_stack, "]]> forbidden in character context")
          else
            tempString => fx%token
            fx%token => vs_str_alloc(str_vs(tempString)//">")
            deallocate(tempString)
          endif
        elseif (phrase==1) then
          tempString => fx%token
          fx%token => vs_str_alloc(str_vs(tempString)//"]"//c)
          deallocate(tempString)
        elseif (phrase==2) then
          tempString => fx%token
          fx%token => vs_str_alloc(str_vs(tempString)//"]]"//c)
          deallocate(tempString)
        else
          tempString => fx%token
          fx%token => vs_str_alloc(str_vs(tempString)//c)
          deallocate(tempString)
        endif

      case (ST_TAG_IN_CONTENT)
        if (phrase==0) then
          if (c=="<") then
            phrase = 1
            ws_discard = .false.
          elseif (c=="&") then
            fx%tokenType = TOK_ENTITY
          else
            call add_error(fx%error_stack, "Unexpected character found in content")
          endif
        elseif (phrase==1) then
          if (c=="?") then
            fx%tokenType = TOK_PI_TAG
          elseif (c=="!") then
            fx%tokenType = TOK_BANG_TAG
          elseif (c=="/") then
            fx%tokenType = TOK_CLOSE_TAG
          elseif (isInitialNameChar(c, xv)) then
            call push_chars(fb, c)
            fx%tokenType = TOK_OPEN_TAG
          else
            call add_error(fx%error_stack, "Unexpected character after <")
          endif
        endif

      case (ST_START_ENTITY)
        if (verify(c,XML_WHITESPACE//";")>0) then
          tempString => fx%token
          fx%token => vs_str_alloc(str_vs(tempString)//c)
          deallocate(tempString)
        elseif (c==";") then
          fx%tokenType = TOK_NAME
        else
          call add_error(fx%error_stack, "Entity reference must be terminated with a ;")
        endif

      case (ST_CLOSING_TAG)
        if (verify(c,XML_WHITESPACE//">")>0) then
          tempString => fx%token
          fx%token => vs_str_alloc(str_vs(tempString)//c)
          deallocate(tempString)
        else
          fx%tokenType = TOK_NAME
          if (c==">") fx%nextTokenType = TOK_END_TAG
        endif

      case (ST_IN_CLOSING_TAG)
        if (verify(c, XML_WHITESPACE)>0) then
          if (c==">") then
            fx%tokenType = TOK_END_TAG
          else
            call add_error(fx%error_stack, "Unexpected character - expecting >")
          endif
        endif

      case (ST_IN_DOCTYPE, ST_DOC_NAME, ST_DOC_SYSTEM, ST_DOC_PUBLIC, &
        ST_DOC_DECL, ST_CLOSE_DOCTYPE)
        if (firstChar) then
          if (verify(c, XML_WHITESPACE)>0) then
            if (c==">") then
              fx%tokenType = TOK_END_TAG
            elseif (c=="[") then
              fx%tokenType = TOK_OPEN_SB
            else
              call add_error(fx%error_stack, &
                "Missing whitespace in doctype delcaration.")
            endif
          else
            ws_discard = .true.
          endif
        endif
        if (ws_discard) then
          if (verify(c, XML_WHITESPACE)>0) then
            if (verify(c, "'""")==0) then
              q = c
              deallocate(fx%token)
              fx%token => vs_str_alloc("")
              ws_discard = .false.
            elseif (c=="[") then
              fx%tokenType = TOK_OPEN_SB
            elseif (c==">") then
              fx%tokenType = TOK_END_TAG
            else
              deallocate(fx%token)
              fx%token => vs_str_alloc(c)
              ws_discard = .false.
            endif
          endif
        else
          if (q/=" ".and.c==q) then
            fx%tokenType = TOK_CHAR
          elseif (q==" ".and.verify(c, ">")==0.and.(fx%state==ST_IN_DOCTYPE)) then
            fx%nextTokenType = TOK_END_TAG
            fx%tokenType = TOK_NAME
          elseif (q==" ".and.verify(c, XML_WHITESPACE)==0) then
            call push_chars(fb, c)
            fx%tokenType = TOK_NAME
          else
            tempString => fx%token
            fx%token => vs_str_alloc(str_vs(tempString)//c)
            deallocate(tempString)
          endif
        endif

      case (ST_IN_SUBSET)
        call tokenizeDTD

      case (ST_START_PE)
        if (verify(c,XML_WHITESPACE//";")>0) then
          tempString => fx%token
          fx%token => vs_str_alloc(str_vs(tempString)//c)
          deallocate(tempString)
        elseif (c==";") then
          fx%tokenType = TOK_NAME
        else
          call add_error(fx%error_stack, "Entity reference must be terminated with a ;")
        endif

      end select
      
      firstChar = .false.
      if (fx%tokenType/=TOK_NULL) exit
    enddo

  contains

    subroutine tokenizeDTD

      if (c=="%") then
        if (fx%state_dtd==ST_DTD_START_COMMENT &
          .or.fx%state_dtd==ST_DTD_START_PI &
          .or.fx%state_dtd==ST_DTD_PI_CONTENTS &
          .or.fx%state_dtd==ST_DTD_ENTITY &
          .or.fx%state_dtd==ST_DTD_NOTATION_SYSTEM &
          .or.fx%state_dtd==ST_DTD_NOTATION_PUBLIC &
          .or.fx%state_dtd==ST_DTD_NOTATION_PUBLIC_2 &
          .or.fx%state_dtd==ST_DTD_ENTITY_PUBLIC &
          .or.fx%state_dtd==ST_DTD_ENTITY_SYSTEM) then
          ! % is perfectly legitimate
          continue
        elseif (fx%state_dtd==ST_DTD_SUBSET) then
          if (.not.fx%spaceBeforeEntity) then
            fx%spaceBeforeEntity = .true.
            call push_chars(fb, c)
            c = " "
          else
            fx%spaceBeforeEntity = .false.
            fx%tokenType = TOK_ENTITY
            return
          endif
        elseif (fx%state_dtd==ST_DTD_ATTLIST_CONTENTS &
          .and. q/="") then
            ! We are inside a ATTLIST attvalue, so apparent PErefs arent.
            continue
        elseif (fx%inIntSubset) then
          call add_error(fx%error_stack, &
            "Parameter entity reference not permitted inside markup for internal subset")
          return
        elseif (fx%state_dtd==ST_DTD_ATTLIST_CONTENTS &
          .or.fx%state_dtd==ST_DTD_ELEMENT_CONTENTS) then
          if (.not.associated(fx%content)) then
            ! content will not always be empty here;
            ! if we have two PErefs bang next to each other.
            fx%content => fx%token
            fx%token => vs_str_alloc("")
          endif
          fx%tokenType = TOK_ENTITY
          return
        elseif (fx%state_dtd==ST_DTD_ENTITY_ID) then
          ! % is ok if we are in the external subset
          continue
        else
          if (.not.fx%spaceBeforeEntity) then
            fx%spaceBeforeEntity = .true.
            call push_chars(fb, c)
            c = " "
          else
            fx%spaceBeforeEntity = .false.
            fx%tokenType = TOK_ENTITY
            return
          endif
        endif
      endif

      select case(fx%state_dtd)

      case (ST_DTD_SUBSET)
        if (firstChar) ws_discard = .true.
        if (ws_discard) then
          if (verify(c, XML_WHITESPACE)>0) then
            if (c=="]") then
              phrase = 1
              q = c
              ws_discard = .false.
            elseif (c=="<") then
              phrase = 1
              q = c
              ws_discard = .false.
            else
              call add_error(fx%error_stack, "Unexpected character found in document subset")
            endif
          endif
        elseif (phrase==1) then
          if (q=="<") then
            if (c=="?") then
              fx%tokenType = TOK_PI_TAG
            elseif (c=="!") then
              fx%tokenType = TOK_BANG_TAG
            else
              call add_error(fx%error_stack, "Unexpected character, expecting ! or ?")
            endif
          elseif (q=="]") then
            if (c=="]") then
              phrase = 2
            else
              fx%tokenType = TOK_CLOSE_SB
              call push_chars(fb, c)
              if (fx%inIntSubset) then
                tempString => fx%xds%intSubset
                allocate(fx%xds%intSubset(size(tempString)-2))
                fx%xds%intSubset = tempString(:size(tempString)-2)
                deallocate(tempString)
              endif
            endif
          endif
        elseif (phrase==2) then
          if (c==">") then
            fx%tokenType = TOK_SECTION_END
          else
            call add_error(fx%error_stack, "Unexpected character, expecting >")
          endif
        endif


      case (ST_DTD_BANG_TAG)
        if (firstChar) then
          if (c=="-") then
            phrase = 1
          elseif (c=="[") then
            fx%tokenType = TOK_OPEN_SB
          elseif (verify(c,upperCase)==0) then
            deallocate(fx%token)
            fx%token => vs_str_alloc(c)
          else
            call add_error(fx%error_stack, "Unexpected character after <!")
          endif
        elseif (phrase==1) then
          if (c=="-") then
            fx%tokenType = TOK_OPEN_COMMENT
          else
            call add_error(fx%error_stack, "Unexpected character after <!-")
          endif
        elseif (verify(c,XML_WHITESPACE)>0) then
          tempString => fx%token
          fx%token => vs_str_alloc(str_vs(tempString)//c)
          deallocate(tempString)
        else
          call push_chars(fb, c)
          fx%tokenType = TOK_NAME
        endif


      case (ST_DTD_START_SECTION_DECL)
        if (firstChar) then
          ws_discard = .true.
        endif
        if (ws_discard) then
          if (verify(c, XML_WHITESPACE)/=0) then
            deallocate(fx%token)
            fx%token => vs_str_alloc(c)
            ws_discard = .false.
          endif
        else
          if (verify(c, XML_WHITESPACE//"[")>0) then
            tempString => fx%token
            fx%token => vs_str_alloc(str_vs(tempString)//c)
            deallocate(tempString)
          else
            fx%tokenType = TOK_NAME
            if (c=="[") fx%nextTokenType = TOK_OPEN_SB
          endif
        endif

      case (ST_DTD_FINISH_SECTION_DECL)
        if (verify(c, XML_WHITESPACE)>0) then
          if (c/="[") then
            call add_error(fx%error_stack, &
              "Unexpected token found, expecting [")
          else
            fx%tokenType = TOK_OPEN_SB
          endif
        endif

      case (ST_DTD_IN_IGNORE_SECTION)
        select case(phrase)
        case (0)
          if (c=="<".or.c=="]") then
            phrase = 1
            q = c
          endif
        case (1)
          if ((q=="<".and.c=="!").or.(q=="]".and.c=="]")) then
            phrase = 2
          else
            phrase = 0
          endif
        case (2)
          if (q=="<".and.c=="[") then
            fx%tokenType = TOK_SECTION_START
          elseif (q=="]".and.c==">") then
            fx%tokenType = TOK_SECTION_END
          else
            phrase = 0
          endif
        end select


      case (ST_DTD_START_PI)
        ! grab until whitespace or ?
        if (verify(c, XML_WHITESPACE//"?")>0) then
          tempString => fx%token
          fx%token => vs_str_alloc(str_vs(tempString)//c)
          deallocate(tempString)
        else
          fx%tokenType = TOK_NAME
          if (c=="?") call push_chars(fb, c)
        endif

      case (ST_DTD_PI_CONTENTS)
        if (firstChar) ws_discard = .true.
        if (ws_discard) then
          if (verify(c, XML_WHITESPACE)/=0) then
            ws_discard = .false.
          else
            return
          endif
        endif
        if (phrase==1) then
          if (c==">") then
            fx%tokenType = TOK_CHAR
            fx%nextTokenType = TOK_PI_END
          elseif (c=="?") then
            ! The last ? didn't mean anything, but this one might.
            tempString => fx%token
            fx%token => vs_str_alloc(str_vs(tempString)//"?")
            deallocate(tempString)
          else
            phrase = 0
            tempString => fx%token
            fx%token => vs_str_alloc(str_vs(tempString)//"?"//c)
            deallocate(tempString)
          endif
        elseif (c=="?") then
          phrase = 1
        else
          tempString => fx%token
          fx%token => vs_str_alloc(str_vs(tempString)//c)
          deallocate(tempString)
        endif

      case (ST_DTD_START_COMMENT)
        select case(phrase)
        case (0)
          if (c=="-") then
            phrase = 1
          else
          tempString => fx%token
          fx%token => vs_str_alloc(str_vs(tempString)//c)
          deallocate(tempString)
          endif
        case (1)
          if (c=="-") then
            phrase = 2
          else
          tempString => fx%token
          fx%token => vs_str_alloc(str_vs(tempString)//"-"//c)
          deallocate(tempString)
            phrase = 0
          endif
        case (2)
          if (c==">") then
            fx%tokenType = TOK_CHAR
            tempString => fx%token
            fx%token => vs_str_alloc(str_vs(tempString)//c)
            deallocate(tempString)
            fx%nextTokenType = TOK_COMMENT_END
          else
            call add_error(fx%error_stack, &
              "Expecting > after -- inside a comment.")
          endif
        end select

      case (ST_DTD_ATTLIST, ST_DTD_ELEMENT, ST_DTD_ENTITY, &
        ST_DTD_ENTITY_PE, ST_DTD_NOTATION)
        if (firstChar) ws_discard = .true.
        if (ws_discard) then
          if (verify(c,XML_WHITESPACE)>0) then
            deallocate(fx%token)
            fx%token => vs_str_alloc(c)
            ws_discard = .false.
          endif
        elseif (verify(c,XML_WHITESPACE//">")>0) then
          tempString => fx%token
          fx%token => vs_str_alloc(str_vs(tempString)//c)
          deallocate(tempString)
        else
          if (c==">") then
            fx%nextTokenType = TOK_END_TAG
          else
            call push_chars(fb, c)
          endif
          fx%tokenType = TOK_NAME
        endif
        
      case (ST_DTD_ELEMENT_CONTENTS)
        if (c==">") then
          if (associated(fx%content)) then
            deallocate(fx%token)
            fx%token => fx%content
            fx%content => null()
          endif
          fx%tokenType = TOK_DTD_CONTENTS
          fx%nextTokenType = TOK_END_TAG
        else
          if (associated(fx%content)) then
            deallocate(fx%token)
            fx%token => vs_str_alloc(str_vs(fx%content)//c)
            deallocate(fx%content)
          else
            tempString => fx%token
            fx%token => vs_str_alloc(str_vs(tempString)//c)
            deallocate(tempString)
          endif
          if (c=="(") then
            fx%tokenType = TOK_OPEN_PAR
            fx%content => fx%token
            fx%token => vs_str_alloc("")
          elseif (c==")") then
            fx%tokenType = TOK_CLOSE_PAR
            fx%content => fx%token
            fx%token => vs_str_alloc("")
          endif
        endif

      case (ST_DTD_ATTLIST_CONTENTS)
        if (c==">") then
          fx%tokenType = TOK_DTD_CONTENTS
          fx%nextTokenType = TOK_END_TAG
        elseif (associated(fx%content)) then
          deallocate(fx%token)
          fx%token => vs_str_alloc(str_vs(fx%content)//c)
          deallocate(fx%content)
        else
          tempString => fx%token
          fx%token => vs_str_alloc(str_vs(tempString)//c)
          deallocate(tempString)
        endif
        if (c=="'".or.c=="""") then
          if (q==c) then
            q = ""
          else
            q = c
          endif
        endif

      case (ST_DTD_ENTITY_ID, ST_DTD_ENTITY_PUBLIC, ST_DTD_ENTITY_SYSTEM, &
        ST_DTD_ENTITY_NDATA, ST_DTD_ENTITY_END, ST_DTD_ENTITY_NDATA_VALUE, &
        ST_DTD_NOTATION_ID, ST_DTD_NOTATION_SYSTEM, ST_DTD_NOTATION_PUBLIC, &
        ST_DTD_NOTATION_PUBLIC_2, ST_DTD_NOTATION_END)
        if (firstChar) then
          if (verify(c, XML_WHITESPACE)>0) then
            if (c==">") then
              fx%tokenType = TOK_END_TAG
            else
              call add_error(fx%error_stack, "Missing whitespace in DTD.")
            endif
          else
            ws_discard = .true.
          endif
        elseif (ws_discard) then
          if (verify(c, XML_WHITESPACE)>0) then
            if (verify(c, "'""")==0) then
              q = c
              deallocate(fx%token)
              fx%token => vs_str_alloc("")
              ws_discard = .false.
            elseif (c==">") then
              fx%tokenType = TOK_END_TAG
            else
              deallocate(fx%token)
              fx%token => vs_str_alloc(c)
              ws_discard = .false.
            endif
          endif
        else
          if (q/=" ".and.c==q) then
            fx%tokenType = TOK_CHAR
          elseif (q==" ".and.verify(c, XML_WHITESPACE//">")==0) then
            fx%tokenType = TOK_NAME
            if (c==">") then
              fx%nextTokenType = TOK_END_TAG
            else
              call push_chars(fb, c)
            endif
          else
            tempString => fx%token
            fx%token => vs_str_alloc(str_vs(tempString)//c)
            deallocate(tempString)
          endif
        endif

      end select
    end subroutine tokenizeDTD

  end subroutine sax_tokenize


  recursive function normalize_attribute_text(fx, s_in) result(s_out)
    type(sax_parser_t), intent(inout) :: fx
    character, dimension(:), intent(in) :: s_in
    character, dimension(:), pointer :: s_out

    character, dimension(:), pointer :: s_temp, s_temp2, s_ent, tempString
    character :: dummy
    integer :: i, i2, j
    type(entity_t), pointer :: ent
#ifdef PGF90
    type(URI), pointer :: nullURI
    nullURI => null()
#endif

    ! Condense all whitespace, only if we are validating,
    ! Expand all &
    ! Complain about < and &

    allocate(s_temp(size(s_in))) ! in the first instance
    allocate(s_out(0)) ! in case we return early ...
    s_ent => null()
    tempString => null()

    i2 = 1
    i = 1
    do 
      if (i > size(s_in)) exit
      ! Firstly, all whitespace must become 0x20
      if (verify(s_in(i),XML_WHITESPACE)==0) then
        s_temp(i2) = " "
        ! Then, < is always illegal
        i = i + 1
        i2 = i2 + 1
      elseif (s_in(i)=='<') then
        call add_error(fx%error_stack, "Illegal < found in attribute.")
        goto 100
        ! Then, expand <
      elseif (s_in(i)=='&') then
        j = index(str_vs(s_in(i+1:)), ';')
        if (j==0) then
          call add_error(fx%error_stack, "Illegal & found in attribute")
          goto 100
        elseif (j==1) then
          call add_error(fx%error_stack, "No entity reference found")
          goto 100
        endif
        allocate(tempString(j-1))
        tempString = s_in(i+1:i+j-1)
        if (existing_entity(fx%predefined_e_list, str_vs(tempString))) then
          ! Expand immediately
          s_temp(i2) = expand_entity_text(fx%predefined_e_list, str_vs(tempString))
          i = i + j + 1
          i2 = i2 + 1
        elseif (checkCharacterEntityReference(str_vs(tempString), fx%xds%xml_version)) then
          ! Expand all character entities
          s_temp(i2) = expand_char_entity(str_vs(tempString))
          i = i + j  + 1
          i2 = i2 + 1 ! fixme
        elseif (checkName(str_vs(tempString), fx%xds%xml_version)) then
          ent => getEntityByName(fx%forbidden_ge_list, str_vs(tempString))
          if (associated(ent)) then
            call add_error(fx%error_stack, 'Recursive entity expansion')
            goto 100
          else
            ent => getEntityByName(fx%xds%entityList, str_vs(tempString))
          endif
          if (associated(ent)) then
            if (ent%wfc.and.fx%xds%standalone) then
              call add_error(fx%error_stack, 'Externally declared entity referenced in standalone document')
              goto 100
            endif
            !is it the right sort of entity?
            if (ent%external) then
              call add_error(fx%error_stack, "External entity forbidden in attribute")
              goto 100
            endif
#ifdef PGF90
            call add_internal_entity(fx%forbidden_ge_list, str_vs(tempString), "", nullURI, .false.)
#else
            call add_internal_entity(fx%forbidden_ge_list, str_vs(tempString), "", null(), .false.)
#endif
            ! Recursively expand entity, checking for errors.
            s_ent => normalize_attribute_text(fx, &
              vs_str(expand_entity_text(fx%xds%entityList, str_vs(tempString))))
            dummy = pop_entity_list(fx%forbidden_ge_list)
            if (in_error(fx%error_stack)) then
              goto 100
            endif
            allocate(s_temp2(size(s_temp)+size(s_ent)-j))
            s_temp2(:i2-1) = s_temp(:i2-1)
            s_temp2(i2:i2+size(s_ent)-1) = s_ent
            deallocate(s_temp)
            s_temp => s_temp2
            nullify(s_temp2)
            i = i + j + 1
            i2 = i2 + size(s_ent)
            deallocate(s_ent)
          else
            s_temp(i2:i2+j) = s_in(i:i+j)
            i = i + j + 1
            i2 = i2 + j + 1
            if (.not.fx%skippedExternal.or.fx%xds%standalone) then
              call add_error(fx%error_stack, "Undeclared entity encountered in standalone document.")
              goto 100
            endif
          endif
        else
          call add_error(fx%error_stack, "Illegal entity reference")
          goto 100
        endif
        deallocate(tempString)
      else
        s_temp(i2) = s_in(i)
        i = i + 1
        i2 = i2 + 1
      endif
    enddo

    deallocate(s_out)
    allocate(s_out(i2-1))
    s_out = s_temp(:i2-1)
100 deallocate(s_temp)
    if (associated(s_ent))  deallocate(s_ent)
    if (associated(tempString)) deallocate(tempString)

  end function normalize_attribute_text

  recursive function expand_pe_text(fx, s_in, fb) result(s_out)
    type(sax_parser_t), intent(inout) :: fx
    character, dimension(:), intent(in) :: s_in
    type(file_buffer_t), intent(inout) :: fb
    character, dimension(:), pointer :: s_out

    character, dimension(:), pointer :: s_temp, s_temp2, s_ent, tempString
    character :: dummy
    integer :: i, i2, j, iostat
    type(entity_t), pointer :: ent
#ifdef PGF90
    type(URI), pointer :: nullURI
    nullURI => null()
#endif

    ! Expand all %PE;

    allocate(s_temp(size(s_in))) ! in the first instance
    allocate(s_out(0)) ! in case we return early ...
    s_ent => null()
    tempString => null()
    s_temp2 => null()

    i2 = 1
    i = 1
    do while (i <= size(s_in))
      if (s_in(i)=='%') then
        j = index(str_vs(s_in(i+1:)), ';')
        if (j==0) then
          call add_error(fx%error_stack, "Illegal % found in attribute")
          goto 100
        elseif (j==1) then
          call add_error(fx%error_stack, "No entity reference found")
          goto 100
        endif
        allocate(tempString(j-1))
        tempString = s_in(i+1:i+j-1)
        if (checkName(str_vs(tempString), fx%xds%xml_version)) then
          ent => getEntityByName(fx%forbidden_pe_list, str_vs(tempString))
          if (associated(ent)) then
            call add_error(fx%error_stack, 'Recursive entity expansion')
            goto 100
          endif
          ent => getEntityByName(fx%xds%peList, str_vs(tempString))
          if (associated(ent)) then
            if (ent%wfc.and.fx%xds%standalone) then
              call add_error(fx%error_stack, &
                "Externally declared entity used in standalone document")
              goto 100
            elseif (str_vs(ent%notation)/="") then
              call add_error(fx%error_stack, "Unparsed entity reference forbidden in entity value")
              goto 100
            endif
#ifdef PGF90
            call add_internal_entity(fx%forbidden_pe_list, str_vs(tempString), "", nullURI, .false.)
#else
            call add_internal_entity(fx%forbidden_pe_list, str_vs(tempString), "", null(), .false.)
#endif
            ! Recursively expand entity, checking for errors.
            if (ent%external) then
              call open_new_file(fb, ent%baseURI, iostat)
              if (iostat/=0) then
                call add_error(fx%error_stack, "Unable to access external parameter entity")
                goto 100
              endif
              call parse_text_declaration(fb, fx%error_stack)
              if (in_error(fx%error_stack)) goto 100
              s_temp2 => get_all_characters(fb, fx%error_stack)
              call pop_buffer_stack(fb)
              if (in_error(fx%error_stack)) goto 100
              s_ent => expand_pe_text(fx, s_temp2, fb)
              deallocate(s_temp2)
            else
              s_ent => expand_pe_text(fx, &
                vs_str(expand_entity_text(fx%xds%peList, str_vs(tempString))), fb)
            endif
            dummy = pop_entity_list(fx%forbidden_pe_list)
            if (in_error(fx%error_stack)) then
              goto 100
            endif
            allocate(s_temp2(size(s_temp)+size(s_ent)-j))
            s_temp2(:i2-1) = s_temp(:i2-1)
            s_temp2(i2:i2+size(s_ent)-1) = s_ent
            deallocate(s_temp)
            s_temp => s_temp2
            s_temp2 => null()
            i = i + j + 1
            i2 = i2 + size(s_ent)
            deallocate(s_ent)
          else
            s_temp(i2:i2+j) = s_in(i:i+j)
            i = i + j + 1
            i2 = i2 + j + 1
            if (.not.fx%skippedExternal.or.fx%xds%standalone) then
              call add_error(fx%error_stack, "Reference to undeclared parameter entity encountered in standalone document.")
              goto 100
            endif
          endif
        else
          call add_error(fx%error_stack, "Illegal parameter entity reference")
          goto 100
        endif
        deallocate(tempString)
      else
        s_temp(i2) = s_in(i)
        i = i + 1
        i2 = i2 + 1
      endif
    enddo

    deallocate(s_out)
    allocate(s_out(i2-1))
    s_out = s_temp(:i2-1)
100 deallocate(s_temp)
    if (associated(s_temp2))  deallocate(s_temp2)
    if (associated(s_ent))  deallocate(s_ent)
    if (associated(tempString)) deallocate(tempString)

  end function expand_pe_text
#endif

end module m_sax_tokenizer
