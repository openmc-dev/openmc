TOHW_m_dom_publics(`
  
  public :: getNotationName

  public :: getIllFormed
  public :: setIllFormed

')`'dnl
dnl
TOHW_m_dom_contents(`

  TOHW_m_dom_get(DOMString, notationName, np%dtdExtras%notationName, (ENTITY_NODE))

!Internally-used getters/setters:

  TOHW_m_dom_get(logical, illFormed, np%dtdExtras%illFormed, (ENTITY_NODE))
  TOHW_m_dom_set(logical, illFormed, np%dtdExtras%illFormed, (ENTITY_NODE))

  TOHW_m_dom_get(DOMString, stringValue, np%nodeValue, (ENTITY_NODE))
  TOHW_m_dom_set(DOMString, stringValue, np%nodeValue, (ENTITY_NODE))

')`'dnl
