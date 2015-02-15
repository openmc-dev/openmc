# DOM

## Overview

The FoX DOM interface exposes an API as specified by the W3C DOM Working group.

FoX implements essentially all of DOM Core Levels 1 and 2, (there are a number of minor exceptions which are listed below) and a substantial portion of DOM Core Level 3.

* [Quick overview of how to map the DOM interface to Fortran](#DomQuickOverview)  
* [More detailed explanation of Fortran interface](#DomDetailedInterface)  
* [Additional (non-DOM) utility functions](#DomUtilityFunctions)  
* [String handling](#DomString)  
* [Exception handling](#DomException)  
* [Live nodelists](#DomLiveNodelists)  
* [DOM Configuration](#DomConfiguration)  
* [Miscellanea](#DomMiscellanea)

## Interface Mapping

<a name="DomQuickOverview"/>

FoX implements all objects and methods mandated in DOM Core Level 1 and 2. (A listing of supported DOM Core Level 3 interfaces is given below.)

In all cases, the mapping from DOM interface to Fortran implementation is as follows:

1. All DOM objects are available as Fortran types, and should be referenced only as pointers (though see 7 and 8 below). Thus, to use a Node, it must be declared first as:    
`type(Node), pointer :: aNode`
2. A flat (non-inheriting) object hierarchy is used. All DOM objects which inherit from Node are represented as Node types.   
3. All object method calls are modelled as functions or subroutines with the same name, whose first argument is the object. Thus:  
`aNodelist = aNode.getElementsByTagName(tagName)`  
should be converted to Fortran as:  
`aNodelist => getElementsByTagName(aNode, tagName)`  
4. All object method calls whose return type is void are modelled as subroutines. Thus:   
`aNode.normalize()`  
becomes
`call normalize(aNode)`   
5. All object attributes are modelled as a pair of get/set calls (or only get where the attribute is readonly), with the naming convention being merely to prepend get or set to the attribute name. Thus:   
`name = node.nodeName`  
`node.nodeValue = string`  
should be converted to Fortran as    
`name = getnodeName(node)`  
`call setnodeValue(string)`  
6. Where an object method or attribute getter returns a DOM object, the relevant Fortran function must always be used as a pointer function. Thus:  
`aNodelist => getElementsByTagName(aNode, tagName)`  
7. No special DOMString object is used - all string operations are done on the standard Fortran character strings, and all functions that return DOMStrings return Fortran character strings.
8. Exceptions are modelled by every DOM subroutine/function allowing an optional additional argument, of type DOMException. For further information see (|DOM Exceptions|) below. 

### String handling


<a name="DomString"/>

The W3C DOM requires that a `DOMString` object exist, capable of holding Unicode strings; and that all DOM functions accept and emit DOMString objects when string data is to be transferred.

FoX does not follow this model. Since (as mentioned elsewhere) it is impossible to perform Unicode I/O in standard Fortran, it would be obtuse to require users to manipulate additional objects merely to transfer strings. Therefore, wherever the DOM mandates use of a `DOMString`, FoX merely uses standard Fortran character strings.

All functions or subroutines which expect DOMString input arguments should be used with normal character strings.  
All functions which should return DOMString objects will return Fortran character strings.

### Using the FoX DOM library.

All functions are exposed through the module `FoX_DOM`. `USE` this in your program:

    program dom_example

      use FoX_DOM
      type(Node) :: myDoc

      myDoc => parseFile("fileIn.xml")
      call serialize(myDoc, "fileOut.xml")
    end program dom_example



## Documenting DOM functions

<a name="DomDetailedInterface"/>

This manual will not exhaustively document the functions available through the `Fox_DOM` interface. Primary documentation may be found in the W3C DOM specifications:`

* [DOM Core Level 1](http://www.w3.org/TR/REC-DOM-Level-1/)
* [DOM Core Level 2](http://www.w3.org/TR/DOM-Level-2-Core/)
* [DOM Core Level 3](http://www.w3.org/TR/DOM-Level-3-Core/)

The systematic rules for translating the DOM interfaces to Fortran are given in the previous section. For completeness, though, there is a list here. The W3C specifications should be consulted for the use of each.


DOMImplementation:  
`type(DOMImplementation), pointer`

* `hasFeature(impl, feature, version)`  
* `createDocumentType(impl, qualifiedName, publicId, systemId)`  
* `createDocument(impl, qualifiedName, publicId, systemId)`  

Document: 
`type(Node), pointer`

* `getDocType(doc)`
* `getImplementation(doc)`    
* `getDocumentElement(doc)`    
* `createElement(doc, tagname)`    
* `createDocumentFragment(doc)`    
* `createTextNode(doc, data)`    
* `createComment(doc, data)`    
* `createCDataSection(doc, data)`    
* `createProcessingInstruction(doc, target, data)`    
* `createAttribute(doc, name)`    
* `createEntityReference(doc, name)`    
* `getElementsByTagName(doc, tagname)`    
* `importNode(doc, importedNode, deep)`    
* `createElementNS(doc, namespaceURI, qualifiedName)`    
* `createAttributeNS(doc, namespaceURI, qualifiedName)`    
* `getElementsByTagNameNS(doc, namespaceURI, qualifiedName)`    
* `getElementById(doc, elementId)`</li>

Node:   
`type(Node), pointer`

* `getNodeName(arg)`    
* `getNodeValue(arg)`    
* `setNodeValue(arg, value)`    
* `getNodeType(arg)`    
* `getParentNode(arg)`    
* `getChildNodes(arg)`    
* `getFirstChild(arg)`    
* `getLastChild(arg)`    
* `getPreviousSibling(arg)`    
* `getNextSibling(arg)`    
* `getAttributes(arg)`    
* `getOwnerDocument(arg)`    
* `insertBefore(arg, newChild, refChild)`    
* `replaceChild(arg, newChild, refChild)`    
* `removeChild(arg, oldChild)`    
* `appendChild(arg, newChild)`    
* `hasChildNodes(arg)`  
* `cloneNode(arg, deep)`  
* `normalize`    
* `isSupported(arg, feature, version)`    
* `getNamespaceURI(arg)`    
* `getPrefix(arg)`    
* `setPrefix(arg, prefix)`  
* `getLocalName(arg)`    
* `hasAttributes(arg)`    
</ul>

NodeList:  
`type(NodeList), pointer`

* `item(arg, index)`    
* `getLength(arg)`  

NamedNodeMap:  
`type(NamedNodeMap), pointer`

* `getNamedItem(map, name)`    
* `setNamedItem(map, arg)`    
* `removeNamedItem(map, name)`    
* `item(map, index)`    
* `getLength(map)`    
* `getNamedItemNS(map, namespaceURI, qualifiedName)`    
* `setNamedItemNS(map, arg)`    
* `removeNamedItemNS(map, namespaceURI, qualifiedName)`    

CharacterData:  
`type(Node), pointer`

* `getData(np)`    
* `setData(np, data)`    
* `getLength(np)`    
* `substringData(np, offset, count)`    
* `appendData(np, arg)`    
* `deleteData(np, offset, count)`    
* `replaceData(np, offset, count, arg)`

Attr:  
`type(Node), pointer`

* `getName(np)`    
* `getSpecified(np)`    
* `getValue(np)`    
* `setValue(np, value)`    
* `getOwnerElement(np)`

Element:  
`type(Node), pointer`

* `getTagName(np)`    
* `getAttribute(np, name)`    
* `setAttribute(np, name, value)`    
* `removeAttribute(np, name)`    
* `getAttributeNode(np, name)`    
* `setAttributeNode(np, newAttr)`    
* `removeAttributeNode(np, oldAttr)`    
* `getElementsByTagName(np, name)`  
* `getAttributeNS(np, namespaceURI, qualifiedName)`    
* `setAttributeNS(np, namespaceURI, qualifiedName, value)`    
* `removeAttributeNS(np, namespaceURI, qualifiedName)`    
* `getAttributeNode(np, namespaceURI, qualifiedName)`    
* `setAttributeNode(np, newAttr)`    
* `removeAttributeNode(np, oldAttr)`    
* `getElementsByTagNameNS(np, namespaceURI, qualifiedName)`    
* `hasAttribute(np, name)`    
* `hasAttributeNS(np, namespaceURI, qualifiedName)`  

Text:  
`type(Node), pointer`

* `splitText(np, offset)`    

DocumentType:  
`type(Node), pointer`

* `getName(np)`    
* `getEntites(np)`    
* `getNotations(np)`    
* `getPublicId(np)`    
* `getSystemId(np)`    
* `getInternalSubset(np)`    

Notation:  
`type(Node), pointer` 

* `getPublicId(np)`    
* `getSystemId(np)`    

Entity:  
`type(Node), pointer`

* `getPublicId(np)`    
* `getSystemId(np)`  
* `getNotationName(np)`    

ProcessingInstruction:  
`type(Node), pointer`

* `getTarget(np)`    
* `getData(np)`    
* `setData(np, data)` 

In addition, the following DOM Core Level 3 functions are available:

Document:

* `getDocumentURI(np)`  
* `setDocumentURI(np, documentURI)`  
* `getDomConfig(np)`  
* `getInputEncoding(np)`  
* `getStrictErrorChecking(np)`  
* `setStrictErrorChecking(np, strictErrorChecking)`  
* `getXmlEncoding(np)`  
* `getXmlStandalone(np)`  
* `setXmlStandalone(np, xmlStandalone)`  
* `getXmlVersion(np)`  
* `setXmlVersion(np, xmlVersion)`  
* `adoptNode(np, source)`   
* `normalizeDocument(np)`  
* `renameNode(np, namespaceURI, qualifiedName)`  


Node:

* `getBaseURI(np)`  
* `getTextContent(np)`  
* `setTextContent(np, textContent)`  
* `isEqualNode(np, other)`
* `isSameNode(np)`  
* `isDefaultNamespace(np, namespaceURI)`  
* `lookupPrefix(np, namespaceURI)`  
* `lookupNamespaceURI(np, prefix)`  

Attr:

* `getIsId(np)`

Entity:  

* `getInputEncoding(np)`  
* `getXmlVersion(np)`  
* `getXmlEncoding(np)`  

Text:

* `getIsElementContentWhitespace(np)`  

DOMConfiguration:  
`type(DOMConfiguration)`

* `canSetParameter(arg, name, value)`  
* `getParameter(arg, name)` 
* `getParameterNames(arg)`  
* `setParameter(arg, name)`  

*NB For details on DOMConfiguration, see [below](#DomConfiguration)*


### Object Model


The DOM is written in terms of an object model involving inheritance, but also permits a flattened model. FoX implements this flattened model - all objects descending from the Node are of the opaque type `Node`. Nodes carry their own type, and attempts to call functions defined on the wrong nodetype (for example, getting the `target` of a node which is not a PI) will result in a `FoX_INVALID_NODE` exception.

The other types available through the FoX DOM are:

* `DOMConfiguration`  
* `DOMException`   
* `DOMImplementation`  
* `NodeList`  
* `NamedNodeMap`

### FoX DOM and pointers 

All DOM objects exposed to the user may only be manipulated through pointers. Attempts to access them directly will result in compile-time or run-time failures according to your environment.

This should have little effect on the structure of your programs, except that you must always remember, when calling a DOM function, to perform pointer assignment, not direct assignment, thus:  
`child => getFirstChild(parent)`  
and _not_  
`child = getFirstChild(parent)`  


### Memory handling

Fortran offers no garbage collection facility, so unfortunately a small degree of memory
handling is necessarily exposed to the user.

However, this has been kept to a minimum. FoX keeps track of all memory allocated and used when calling DOM routines, and keeps references to all DOM objects created.

The only memory handling that the user needs to take care of is destroying any
DOM Documents (whether created manually, or by the `parse()` routine.) All other nodes or node structures created will be destroyed automatically by the relevant `destroy()` call.

As a consequence of this, all DOM objects which are part of a given document will become inaccessible after the document object is destroyed.

## Additional functions.

<a name="DomUtilityFunctions"/>

Several additional utility functions are provided by FoX.

### Input and output of XML data

Firstly, to construct a DOM tree, from either a file or a string containing XML data.

* `parseFile`  
**filename**: *string*  
(**configuration**): *DOMConfiguration*  
(**ex**): *DOMException*  

**filename** should be an XML document. It will be opened and parsed into a DOM tree. The parsing is performed by the FoX SAX parser; if the XML document is not well-formed, a `PARSE_ERR` exception will be raised. **configuration** is an optional argument - see [DOMConfiguration](#DomConfiguration) for its meaning.

* `parseString`  
**XMLstring**: *string*  
(**configuration**): *DOMConfiguration*  
(**ex**): *DOMException*

**XMLstring** should be a string containing XML data. It will be parsed into a DOM tree. The parsing is performed by the FoX SAX parser; if the XML document is not well-formed, a `PARSE_ERR` exception will be raised. **configuration** is an optional argument - see [DOMConfiguration](#DomConfiguration) for its meaning.

Both `parseFile` and `parseString` return a pointer to a `Node` object containing the Document Node.`

Secondly, to output an XML document:

* `serialize`  
**arg**: *Node, pointer*
**fileName**: *string* 

This will open `fileName` and serialize the DOM tree by writing into the file. If `fileName` already exists, it will be overwritten. If an problem arises in serializing the document, then a fatal error will result.

(Control over serialization options is done through the configuration of the **arg**'s ownerDocument, see [below](#DomConfiguration).)


Finally, to clean up all memory associated with the DOM, it is necessary to call:

* `destroy`  
**np**: *Node, pointer*

This will clear up all memory usage associated with the document (or documentType) node passed in.


### Extraction of data from an XML file.

<a name="dataExtraction"/>

The standard DOM functions only deal with string data. When dealing with numerical (or logical) data,
the following functions may be of use.

* `extractDataContent`  
* `extractDataAttribute`  
* `extractDataAttributeNS`

These extract data from, respectively, the text content of an element, from one of its attributes, or from one of its namespaced attributes.
They are used like so:

(where `p` is an element which has been selected by means of the other DOM functions)

    call extractDataContent(p, data)

The subroutine will look at the text contents of the element, and interpret according to the type of `data`. That is, if `data` has been declared as an `integer`, then the contents of `p` will be read as such an placed into `data`.

`data` may be a `string`, `logical`, `integer`, `real`, `double precision`, `complex` or `double complex` variable.

In addition, if `data` is supplied as a rank-1 or rank-2 variable (ie an array or a matrix) then the data will be read in assuming it to be a space- or comma-separated list of such data items.

Thus, the array of integers within the XML document:

    <element> 1 2 3 4 5 6 </element>

could be extracted by the following Fortran program:

    type(Node), pointer :: doc, p
    integer :: i_array(6)

    doc => parseFile(filename)
    p => item(getElementsByTagName(doc, "element"), 0)
    call extractDataContent(p, i_array)

#### Contents and Attributes

For extracting data from text content, the example above suffices. For data in a non-namespaced attribute (in this case, a 2x2 matrix of real numbers)

    <element att="0.1, 2.3 7.56e23, 93"> Some uninteresting text </element>

then use a Fortran program like:

    type(Node), pointer :: doc, p
    real :: r_matrix(2,2)

    doc => parseFile(filename)
    p => item(getElementsByTagName(doc, "element"), 0)
    call extractDataAttribute(p, "att", r_matrix)

or for extracting from a namespaced attribute (in this case, a length-2 array of complex numbers):

    <myml xmlns:ns="http://www.example.org">
      <element ns:att="0.1,2.3  3.4e2,5.34"> Some uninteresting text </element>
    </myml>

then use a Fortran program like:

    type(Node), pointer :: doc, p
    complex :: c_array(2)

    doc => parseFile(filename)
    p => item(getElementsByTagName(doc, "element"), 0)
    call extractDataAttributeNS(p, &
         namespaceURI="http://www.example.org", localName="att", &
         data=c_array)

#### Error handling

The extraction may fail of course, if the data is not of the sort specified, or if there are not enough elements to fill the array or matrix. In such a case, this can be detected by the optional arguments `num` and `iostat`. 

`num` will hold the number of items successfully read. Hopefully this should be equal to the expected number of items; but it may be less if reading failed for some reason, or if there were less items than expected in the element.

`iostat` will hold an integer - this will be `0` if the extraction went ok; `-1` if too few elements were found, `1` if although the read went ok, there were still some elements left over, or `2` if the extraction failed due to either a badly formatted number, or due to the wrong data type being found.

#### String arrays

For all data types apart from strings, arrays and matrices are specified by space- or comma-separated lists. For strings, some additional options are available. By default, arrays will be extracted assuming that separators are spaces (and multiple spaces are ignored). So:

    <element> one two     three </element>

will result in the string array `(/"one", "two", "three"/)`.

However, you may specify an optional argument `separator`, which specifies another single-character separator to use (and does not ignore multiple spaces). So:

    <element>one, two, three </element>

will result in the string array `(/"one", " two", " three "/)`. (note the leading and trailing spaces).

Finally, you can also specify an optional logical argument, `csv`. In this case, the `separator` is ignored, and the extraction proceeds assuming that the data is a list of comma-separated values. (see: [CSV](http://en.wikipedia.org/wiki/Comma-separated_values))

### Other utility functions

* `setFoX_checks`  
**FoX_checks**: *logical*

This affects whether additional FoX-only checks are made (see [DomExceptions](#DomException) below). 

* `getFoX_checks`  
**arg**: *DOMImplementation, pointer*  

Retrieves the current setting of FoX_checks.

Note that FoX_checks can only be turned on and off globally, not on a per-document basis.

* `setLiveNodeLists`  
**arg**: *Node, pointer*  
**liveNodeLists**: *logical*

**arg** must be a Document Node. Calling this function affects whether any nodelists active on the document are treated as live - ie whether updates to the documents are reflected in the contents of nodelists (see [DomLiveNodelists](#DomLiveNodelists) below).

* `getLiveNodeLists`  
**arg**: *Node, pointer*

Retrieves the current setting of liveNodeLists.

Note that the live-ness of nodelists is a per-document setting.


### Exception handling

<a name="DomException"/>

Exception handling is important to the DOM. The W3C DOM standards provide not only interfaces to the DOM, but also specify the error handling that should take place when invalid calls are made.

The DOM specifies these in terms of a `DOMException` object, which carries a numeric code whose value reports the kind of error generated. Depending upon the features available in a particular computer language, this DOMException object should be generated and thrown, to be caught by the end-user application.

Fortran of course has no mechanism for throwing and catching exceptions. However, the behaviour of an exception can be modelled using Fortran features.

FoX defines an opaque `DOMException` object.
Every DOM subroutine and function implemented by FoX will take an optional argument, 'ex', of type `DOMException`. 

If the optional argument is not supplied, any errors within the DOM will cause an immediate abort, with a suitable error message. However, if the optional argument *is* supplied, then the error will be captured within the `DOMException` object, and returned to the caller for inspection. It is then up to the application to decide how to proceed.

Functions for inspecting and manipulating the `DOMException` object are described below:

* `inException`:   
**ex**: *DOMException*


A function returning a logical value, according to whether `ex` is in exception - that is, whether the last DOM function or subroutine, from which `ex` returned, caused an error. Note that this will not change the status of the exception.

* `getExceptionCode`   
**ex**: *DOMException*

A function returning an integer value, describing the nature of the exception reported in `ex`. If the integer is 0, then `ex` does not hold an exception. If the integer is less than 200, then the error encountered was of a type specified by the DOM standard; for a full list, see below, and for explanations, see the various DOM standards. If the integer is 200 or greater, then the code represents a FoX-specific error. See the list below.

Note that calling `getExceptionCode` will clean up all memory associated with the DOMException object, and reset the object such that it is no longer in exception.

#### Exception handling and memory usage.

Note that when an Exception is thrown, memory is allocated within the DOMException object. Calling `getExceptionCode` on a DOMEXception will clean up this memory. If you use the exception-handling interfaces of FoX, then you must check every exception, and ensure you check its code, otherwise your program will leak memory.


#### FoX exceptions.

The W3C DOM interface allows the creation of unserializable XML document in various ways. For example, it permits characters to be added to a text node which would be invalid XML. FoX performs multiple additional checks on all DOM calls to prevent the creation of unserializable trees. These are reported through the DOMException mechanisms noted above, using additional exception codes. However, if for some reason, you want to create such trees, then it is possible to switch off all FoX-only checks. (DOM-mandated checks may not be disabled.) To do this, use the `setFoX_checks` function described in [DomUtilityFunctions](#DomUtilityFunctions).

Note that FoX does not yet currently check for all ways that a tree may be made non-serializable.

#### List of exceptions.

The following is the list of all exception codes (both specified in the W3C DOM and those related to FoX-only checks) that can be generated by FoX:

* INDEX\_SIZE\_ERR = 1
* DOMSTRING\_SIZE\_ERR = 2
* HIERARCHY\_REQUEST\_ERR = 3
* WRONG\_DOCUMENT\_ERR = 4
* INVALID\_CHARACTER\_ERR = 5
* NO\_DATA\_ALLOWED\_ERR = 6
* NO\_MODIFICATION\_ALLOWED\_ERR = 7
* NOT\_FOUND\_ERR = 8
* NOT\_SUPPORTED\_ERR = 9
* INUSE\_ATTRIBUTE\_ERR = 10
* INVALID\_STATE\_ERR = 11
* SYNTAX\_ERR = 12
* INVALID\_MODIFICATION\_ERR = 13
* NAMESPACE\_ERR = 14
* INVALID\_ACCESS\_ERR = 15
* VALIDATION\_ERR = 16
* TYPE\_MISMATCH\_ERR = 17
* INVALID\_EXPRESSION\_ERR = 51
* TYPE\_ERR = 52
* PARSE\_ERR = 81
* SERIALIZE\_ERR = 82
* FoX\_INVALID\_NODE = 201
* FoX\_INVALID\_CHARACTER = 202
* FoX\_NO\_SUCH\_ENTITY = 203
* FoX\_INVALID\_PI\_DATA = 204
* FoX\_INVALID\_CDATA\_SECTION = 205
* FoX\_HIERARCHY\_REQUEST\_ERR = 206
* FoX\_INVALID\_PUBLIC\_ID = 207
* FoX\_INVALID\_SYSTEM\_ID = 208
* FoX\_INVALID\_COMMENT = 209
* FoX\_NODE\_IS\_NULL = 210
* FoX\_INVALID\_ENTITY = 211
* FoX\_INVALID\_URI = 212
* FoX\_IMPL\_IS\_NULL = 213
* FoX\_MAP\_IS\_NULL = 214
* FoX\_LIST\_IS\_NULL = 215
* FoX\_INTERNAL\_ERROR = 999

### Live nodelists

<a name="DomLiveNodelists"/>

The DOM specification requires that all NodeList objects are <em>live</em> - that is, that any change in the document structure is immediately reflected in the contents of any nodelists.

For example, any nodelists returned by getElementsByTagName or getElementsByTagNameNS must be updated whenever nodes are added to or removed from the document; and the order of nodes in the nodelists must be changed if the document structure changes.

Though FoX does keep all nodelists live, this can impose a significant performance penalty when manipulating large documents. Therefore, FoX can be instructed to inly use 'dead' nodelists - that is, nodelists which reflect a snapshot of the document structure at the point they were created. To do this, call `setLiveNodeLists` (see API documentation).

However, note that the nodes within the nodelist remain live - any changes made to the nodes will be reflected in accessing them through the nodelist.

Furthermore, since the nodelists are still associated with the document, they and their contents will be rendered inaccessible when the document is destroyed.

## DOM Configuration

<a name="DomConfiguration"/>

Multiple valid DOM trees may be produced from a single document. When parsing input, some of these choices are made available to the user.

By default, the DOM tree presented to the user will be produced according to the following criteria:

* there will be no adjacent text nodes  
* Cdata nodes will appear as such in the DOM tree  
* EntityReference nodes will appear in the DOM tree.

However, if another tree is desired, the user may change this. For example, very often you would rather be working with the fully canonicalized tree, with all cdata sections replaced by text nodes and merged, and all entity references replaced with their contents.

The mechanism for doing this is the optional `configuration` argument to `parseFile` and `parseString`. `configuration` is a `DOMConfiguration` object, which may be manipulated by `setParameter` calls.

Note that FoX's implementation of `DOMConfiguration` does not follow the specification precisely. One `DOMConfiguration` object controls all of parsing, normalization and serialization. It can be used like so:

    use FoX_dom
    implicit none
    type(Node), pointer :: doc
    ! Declare a new configuration object
    type(DOMConfiguration), pointer :: config
    ! Request full canonicalization
    ! ie convert CDATA sections to text sections, remove all entity references etc.
    config => newDOMConfig()
    call setParameter(config, "canonical-form", .true.)
    ! Turn on validation
    call setParameter(config, "validate", .true.)
    ! parse the document
    doc => parseFile("doc.xml", config)

    ! Do a whole lot of DOM processing ...

    ! change the configuration to allow cdata-sections to be preserved.
    call setParameter(getDomConfig(doc), "cdata-sections", .true.)
    ! normalize the document again 
    call normalizeDocument(doc)
    ! change the configuration to influence the output - make sure there is an XML declaration
    call setParameter(getDomConfig(doc), "xml-declaration", .true.)
    ! and write the document out.
    call serialize(doc)
    ! once everything is done, destroy the doc and config
    call destroy(doc)
    call destroy(config)

The available configuration options are fully explained in:

* [DOM Core 3](http://www.w3.org/TR/DOM-Level-3-Core/core.html#DOMConfiguration)  
* [DOM Core LSParser](http://www.w3.org/TR/2004/REC-DOM-Level-3-LS-20040407/load-save.html#LS-LSParser)  
* [DOM Core LSSerializer](http://www.w3.org/TR/2004/REC-DOM-Level-3-LS-20040407/load-save.html#LS-LSSerializer)  

and are all implemented, with the exceptions of: `error-handler`, `schema-location`, and `schema-type`.   
In total there are 24 implemented configuration options (`schema-location` and `schema-type` are not
implemented). The options known by FoX are as follows:

* `canonical-form` default: false, can be set to true. See note below.
* `cdata-sections` default: true, can be changed.
* `check-character-normalization` default: false, cannot be changed.
* `comments` default: true, can be changed.
* `datatype-normalization` default: false, cannot be changed.
* `element-content-whitespace` default: true, can be changed.
* `entities` default: true, can be changed.
* `error-handler` default: false, cannot be changed. This is a breach of the DOM specification.
* `namespaces` default: true, can be changed.
* `namespace-declarations` default: true, can be changed.
* `normalize-characters` default: false, cannot be changed.
* `split-cdata-sections` default: true, can be changed.
* `validate` default: false, can be changed. See note below.
* `validate-if-schema` default: false, can be changed.
* `well-formed` default true, cannot be changed.
* `charset-overrides-xml-encoding` default false, cannot be changed.
* `disallow-doctype` default false, cannot be changed.
* `ignore-unknown-character-denormalizations` default true, cannot be changed.
* `resource-resolver` default false, cannot be changed.
* `supported-media-types-only` default false, cannot be changed.
* `discard-default-content` default: true, can be changed.
* `format-pretty-print` default: false, cannot be changed.
* `xml-declaration` default: true, can be changed.
* `invalid-pretty-print` default: false, can be changed. This is a FoX specific extension which works like `format-pretty-print` but does not preseve the validity of the document.

Setting `canonical-form` changes the value of `entities`, `cdata-sections`, `discard-default-content`, `invalid-pretty-print`, and `xml-declaration`to false and changes `namespaces`, `namespace-declarations`, and `element-content-whitespace` to true. Unsetting `canonical-form` causes these options to revert to the defalt settings. Changing the values of any of these options has the side effect of unsetting `canonical-form` (but does not cause the other options to be reset). Setting `validate` unsets `validate-if-schema` and vica versa.

## DOM Miscellanea

<a name="DomMiscellanea"/>

Other issues


* As mentioned in the documentation for WXML, it is impossible within Fortran to reliably output lines longer than 1024 characters. While text nodes containing such lines may be created in the DOM, on serialization newlines will be inserted as described in the documentation for WXML.
* All caveats with regard to the FoX SAX processor apply to reading documents through the DOM interface. In particular, note that documents containing characters beyond the US-ASCII set will not be readable.

It was decided to implement W3C DOM interfaces primarily because they are specified in a language-agnostic fashion, and thus made Fortran implementation possible. A number of criticisms have been levelled at the W3C DOM, but many apply only from the perspective of Java developers. However, more importantly, the W3C DOM suffers from a lack of sufficient error checking so it is very easy to create a DOM tree, or manipulate an existing DOM tree into a state, that cannot be serialized into a legal XML document.

(Although the Level 3 DOM specifications finally addressed this issue, they did so in a fashion that was neither very useful, nor easily translatable into a Fortran API.)

Therefore, FoX will by default produce errors about many attempts to manipulate the DOM in such a way as would result in invalid XML. These errors can be switched off if standards-compliant behaviour is wanted. Although extensive, these checks are not complete.
In particular, the way the W3C DOM mandates namespace handling makes it trivial to produce namespace non-well-formed document trees, and very difficult for the processor to automatically detect the non-well-formedness. Thus a fully well-formed tree is only guaranteed after a suitable `normalizeDocument` call.
