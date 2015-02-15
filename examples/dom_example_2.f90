program dom_example

  use FoX_dom
  implicit none

  type(Node), pointer :: myDoc, p
  type(NodeList), pointer :: parameterList, children
  integer :: i, j

  real :: energy

  ! Load in the document
  myDoc => parseFile("h2o.xml")

  ! Find all the parameters:
  parameterList => getElementsByTagNameNS(myDoc, &
    "http://www.xml-cml.org/schema", "parameter")

  print*, "Found ", getLength(parameterList), " parameters."

  ! Loop over the parameter list. Note that the DOM
  ! counts from zero, not from one.
  do i = 0, getLength(parameterList)-1
    p => item(parameterList, i)
    ! Check for the existence of the attribute we're looking for
    if (hasAttribute(p, "name")) then
      if (getAttribute(p, "name")=="DM.EnergyTolerance") then
        ! The energy is in the text node which is the child of the <scalar> element under this node ...
        ! Check all the children of the node for the <scalar> element.
        children => getChildNodes(p)
        do j = 0, getLength(children)-1
          p => item(children, j)
          if (getLocalName(p) =="scalar") then
            ! This is the scalar node whose child we want:
            call extractDataContent(p, energy)
            print*, "Energy Tolerance is ", energy
          endif
        enddo
      endif
    endif
  enddo

  ! Clear up all allocated memory
  call destroy(myDoc)
end program dom_example
