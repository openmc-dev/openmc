# WXML

`wxml` is a general Fortran XML output library. It offers a Fortran interface, in the form of a number of subroutines,  to generate well-formed XML documents. Almost all of the XML features described in [XML11](#XML11)  and [Namespaces](#Namespaces) are available, and `wxml` will diagnose almost all attempts to produce an invalid document. [Exceptions](#Exceptions) below describes where `wxml` falls short of these aims.

First, [Conventions](#Conventions) describes the conventions use in this document.

Then, [Functions](#Functions) lists all of `wxml`'s publically exported functions, in three sections:

1. [Firstly](#simple), the very few functions necessary to create the simplest XML document, containing only elements, attributes, and text. 
2. [Secondly](#NSfunctions), those functions concerned with XML Namespaces, and how Namespaces affect the behaviour of the first tranche of functions.  
3. [Thirdly](#obscure), a set of more rarely used functions required to access some of the more esoteric corners of the XML specification.

Please note that where the documentation below is not clear, it may be useful to look at some of the example files. There is a very simple example in the `examples/` subdirectory, but which nevertheless shows the use of most of the features you will use.

A more elaborate example, using almost all of the XML features found here, is available in the top-level directory as `wxml_example.f90`. It will be automatically compiled as part of the build porcess.

<a name="Conventions"/>

##Conventions and notes:

####Conventions used below.

* Function names are in `monospace`
* argument names are in **bold**
* optional argument names are in (**parenthesized bold**)
* argument types are in *italic* and may consist of:
* *string*: string of arbitrary (unless otherwise specified) length
* *integer*: default integer
* *real(sp)*: single precision real number
* *real(dp)*: double precision real number
* *logical*: default logical 
* *real*: either of *real(sp)* or *real(dp)*
* *anytype*: any of *logical*, *integer*, *real(sp)*, *real(dp)*, *string*

Note that where *strings* are passed in, they will be passed through entirely unchanged to the output file - no truncation of whitespace will occur.

It is strongly recommended that the functions be used with keyword arguments rather than replying on implicit ordering.


####Derived type: `xmlf_t`
This is an opaque type representing the XML file handle. Each function requires this as an argument, so it knows which file to operate on. (And it is an output of the xml_OpenFile subroutine) Since all subroutines require it, it is not mentioned below.

<a name="Functions"/>

## Function listing


<a name="simple"/>

### Frequently used functions

* `xml_OpenFile`  
**filename**: *string*: Filename to be opened  
**xf**: *xmlf_t*: XML File handle  
(**channel**): *integer*: What Fortran file handle should the XML file be attached to? 
  *default: picked by the library at runtime*  
(**pretty_print**): *logical*: Should the XML output be formatted to look pretty? (This implies that whitespace is not significant) 
  *default: false*  
(**replace**): *logical*: Should the file be replaced if it already exists? 
  *default: no, stop at runtime if file already exists*  
(**addDecl**): *logical*: Should an XML declaration be added at the start of the file? 
  *default: yes*  
(**namespace**): *logical*: Should wxml prevent the output of namespace-ill-formed documents? 
 *default: yes*  
(**validate**): *logical*: Should wxml carry out any checks on the optional VC constraints specified by XML? 
 *default: no*  
(**warning**): *logical*: Should wxml emit warnings when it is unable to guarantee well-formedness? 
 *default: no*  

Open a file for writing XML

By default, the XML will have no extraneous text nodes. This can have the effect of it
looking slightly ugly, since there will be no newlines inserted between tags.

This behaviour can be changed to produce slightly nicer looking XML, by switching
on pretty_print. This will insert newlines and spaces between some tags where
they are unlikely to carry semantics. Note, though, that this does result in 
the XML produced being not quite what was asked for, since extra characters and
text nodes have been inserted.

NB: The **replace** option should be noted. By default, xml_OpenFile will fail with a runtime error if you try and write to an existing file. If you are sure you want to continue on in such a case, then you can specify `**replace**=.true.` and any existing files will be overwritten. If finer granularity is required over how to proceed in such cases, use the Fortran `inquire` statement in your code. There is no 'append' functionality by design - any XML file created by appending to an existing file would be invalid.

* `xml_Close`   
**xf**: *xmlf_t*: XML File handle
(**empty**): Can the file be empty? *default: .false.*

Close an opened XML file, closing all still-opened tags so that it is well-formed.

In the normal run of event, trying to close an XML file with no root element will cause an error, since this is not well-formed. However, an optional argument, **empty** is provided in case it is desirable to close files which may be empty. In this case, a warning will still be emitted, but no fatal error generated.

* `xml_NewElement`  
**name**: *string*:
 Name of tag (for namespaced output, you need to include the prefix)

Open a new element tag

* `xml_EndElement`  
**name**: *string*: 
 Name of tag to be closed (if it doesn't match currently open tag, you'll get an error)

Close an open tag

* `xml_AddAttribute`  
**name**: *string*: Name of attribute  
**value**: *anytype*: Value of attribute  
(**escape**): *logical*: if the attribute value is a string, should the attribute value be escaped?
  *default: true*  
(**type**): *string*: the type of the attribute. This must be one of `CDATA`, `ID`, `IDREF`, `IDREFS`, `NMTOKEN`, `NMTOKENS`, `ENTITY`, `ENTITIES`, or `NOTATION` (always upper case). If specified, this must match any attribute declarations that have been previously declared in the DTD. If unspecified this (as the XML standard requires) defaults to `CDATA`.

Add an attribute to the currently open tag.

By default, if the attribute value contains markup characters, they will be escaped automatically by
wxml before output.

However, in rare cases you may not wish this to happen - if you wish to output Unicode
characters, or entity references. In this case, you should set `escape=.false.` for the relevant
subroutine call. Note that if you do this, no checking on the validity of the output string iis performed; the onus is on you to ensure well-formedness

The value to be added may be of any type; it will be converted to text according to FoX's [formatting rules](str.html),
and if it is a 1- or 2-dimensional array, the elements will all be output, separated by spaces (except if it is a character array, in which
case the delimiter may be changed to any other single character using an optional argument).

NB The **type** option is only provided so that in the case of an external DTD which FoX is unaware of, the attribute type can be specified (which gives FoX more information to ensure well-formedness and validity). Specifying the type incorrectly may result in spurious error messages)

* `xml_AddCharacters`  
**chars** *anytype*:
 The text to be output  
(**parsed**): *logical*: Should the output characters be parsed (ie should the library replace '&' with '&amp;' etc?) or unparsed (in which case
the characters will be surrounded by CDATA tags.
 *default: yes*  
(**delimiter**): *character(1)*: If **data** is a character array, what should the delimiter between elements be on output?
 *default: a single space*  
(**ws_significant**): *logical*: Is any whitespace in the string significant? *default: unknown*

Add text data. The data to be added may be of any type; they will be converted to text according to FoX's [formatting rules](str.html),
and if they are a 1- or 2-dimensional array, the elements will all be output, separated by spaces (except if it is a character array, in which
case the delimiter may be changed to any other single character using an optional argument).

* `xml_AddNewline`

Within the context of character output, add a (system-dependent) newline character. This function can only
be called wherever `xml_AddCharacters` can be called. (Newlines outside of character context are under
FoX's control, and cannot be manipulated by the user.)

<a name="NSfunctions"/>

### Namespace-aware functions:

* `xml_DeclareNamespace`  
**nsURI** *string*: The URI of the namespace   
(**prefix**) *string*: The namespace prefix to be used in the document. If absent, then the default namespace is affected.

Add an XML Namespace declaration. This function may be called at any time, and its precise effect depends on when it is called; see below

* `xml_UndeclareNamespace`  
(**prefix**) *string*: The namespace prefix to be used in the document. If absent, then the default namespace is affected.

Undeclare an XML namespace. This is equivalent to declaring an namespace with an empty URI, and renders the namespace ineffective for the scope of the declaration. For explanation of its scope, see below.

**NB** Use of `xml_UndeclareNamespace` implies that the resultant document will be compliant with XML Namespaces 1.1, but not 1.0; wxml will issue an error when trying to undeclare namespaces under XML 1.0.

#### Scope of namespace functions

If  `xml_[Un]declareNamespace` is called immediately prior to an `xml_NewElement` call, then the namespace will be declared in that next element, and will therefore take effect in all child elements.

If it is called prior to an `xml_NewElement` call, but that element has namespaced attributes 

To explain by means of example: In order to generate the following XML output:

     <cml:cml xmlns:cml="http://www.xml-cml.org/schema"/>

then the following two calls are necessary, in the prescribed order:

      xml_DeclareNamespace(xf, 'cml', 'http://www.xml-cml.org')
      xml_NewElement(xf, 'cml:cml')

However, to generate XML input like so:
      <cml xhtml:class="symbol" xmlns:xhtml="http://www.w3.org/1999/xhtml"/>
that is, where the namespace refers to an attribute at the same level,
then as long as the `xml_DeclareNamespace` call is made before the element tag is closed (either by `xml_EndElement`, or by a new element tag being opened, or some text being added etc.) the correct XML will be generated.

Two previously mentioned functions are affected when used in a namespace-aware fashion.

* `xml_NewElement`, `xml_AddAttribute`

The element or attribute name is checked, and if it is a QName (ie if it is of the form prefix:tagName) then wxml will check that prefix is a
registered namespace prefix, and generate an error if not.


<a name="obscure"/>

### More rarely used functions:

If you don't know the purpose of any of these, then you don't need to. 


* `xml_AddXMLDeclaration`  
(**version**) *string*: XML version to be used.
 *default: 1.0*  
(**encoding**) *string*: character encoding of the document
  *default: absent*  
(**standalone**) *logical*: is this document standalone?
  *default: absent*  

Add XML declaration to the first line of output. If used, then the file must have been opened with `addDecl = .false.`, and this must be the first wxml call to the document.o


NB The only XML versions available are 1.0 and 1.1. Attempting to specify anything else will result in an error. Specifying version 1.0 results in additional output checks to ensure the resultant document is XML-1.0-conformant.

NB Note that if the encoding is specified, and is specified to not be UTF-8, then if the specified encoding does not match that supported by the Fortran processor, you may end up with output you do not expect.

* `xml_AddDOCTYPE`  
**name** *string*: DOCTYPE name  
(**system**) *string*: DOCTYPE SYSTEM ID  
(**public**) *string*: DOCTYPE PUBLIC ID  

Add an XML document type declaration. If used, this must be used prior to first `xml_NewElement` call, and only one such call must be made.

* `xml_AddInternalEntity`  
**name** *string*: name of internal entity  
**value** *string*: value of internal entity  

Define an internal entity for the document. If used, this call must be made after `xml_AddDOCTYPE` and before the first `xml_NewElement` call.

* `xml_AddExternalEntity`  
**name** *string*: name of external entity  
**system** *string*: SYSTEM ID of external entity  
(**public**) *string*: PUBLIC ID of external entity
  *default: absent*  
(**notation**) *string*: notation for external entity
  *default: absent*  

Define an external entity for the document. If used, this call must be made after `xml_AddDOCTYPE` and before the first `xml_NewElement` call.

* `xml_AddParameterEntity`  
**name** *string*: name of parameter entity  
(**PEdef**) *string*: definition of parameter entity
  *default: absent*  
(**system**) *string*: SYSTEM ID of parameter entity
  *default: absent*  
(**public**) *string*: PUBLIC ID of parameter entity
  *default: absent*  

Define a parameter entity for the document. If used, this call must be made after `xml_AddDOCTYPE` and before the first `xml_NewElement` call.

* `xml_AddNotation`  
**name** *string*: name of notation  
(**system**) *string*: SYSTEM ID of notation
  *default: absent*  
(**public**) *string*: PUBLIC ID of notation
  *default: absent*  

Define a notation for the document. If used, this call must be made after `xml_AddDOCTYPE` and before the first `xml_NewElement` call.

* `xml_AddElementToDTD`  
**name** *string*: name of element  
**declaration** *string*: declaration of element  

Add an ELEMENT declaration to the DTD. The syntax of the declaration is not checked in any way, nor does this affect how elements may be added in the content of the XML document.

If used, this call must be made after `xml_AddDOCTYPE` and before the first `xml_NewElement` call.

* `xml_AddAttlistToDTD`  
**name** *string*: name of element  
**declaration** *string*: declaration of element  

Add an ATTLIST declaration to the DTD. The syntax of the declaration is not checked in any way, nor does this affect how attributes may be added in the content of the XML document.

If used, this call must be made after `xml_AddDOCTYPE` and before the first `xml_NewElement` call.

* `xml_AddPEreferenceToDTD`  
**name** *string*: name of PEreference

Add a reference to a Parameter Entity in the DTD. No check is made according to whether the PE exists, has been declared, or may legally be used.

If used, this call must be made after `xml_AddDOCTYPE` and before the first `xml_NewElement` call.


* `xml_AddXMLStylesheet`  
**href** :*string*: 
 address of stylesheet    
**type**: *string*:
 type of stylesheet (generally "text/xsl")  
(**title**): *string*:
 title of stylesheet
  *default: none*  
(**media**): *string:*
 output media type
*default: none*  
(**charset**): *string*
 charset of media type
 *default:none*  
 (**alternate**): *string*:
 alternate
*default:none*  

Add XML stylesheet processing instruction, as described in [Stylesheets]. If used, this call must be made before the first `xml_NewElement` call.


* `xml_AddXMLPI`   
**name**: *string*:
name of PI    
(**data**): *string*:
data for PI   
(**xml**): *logical*: (see below)
  *default: false*  
(**ws_significant**): *logical*: if this is a PI containing only **data**, then is any whitespace in the data significant? *default: unknown*

Add an XML Processing Instruction.

If data is present, nothing further can be added to the PI. If it is *not* present, then pseudoattributes may be added using the call below.
Normally, the **name** is checked to ensure that it is XML-compliant. This requires that PI targets not start with `[Xx][Mm][Ll]`, because such names are reserved. However, some are defined by later W3 specificataions. If you wish to use such PI targets, then set `xml=.true.` when outputting them.

The output PI will look like:
`<?name data?>`

* `xml_AddPseudoAttribute`  
**name**: *string*:
 Name of pseudoattribute  
**value**: *anytype*:
 Value of pseudoattribute
(**ws_significant**): *logical*: If there is any whitespace in the value of this pseudoattribute, is is significant?

Add a pseudoattribute to the currently open PI.


* `xml_AddComment`  
**comment**: *string*
 Contents of comment  
(**ws_significant**): *logical*: is any whitespace in the comment string significant? *default: unknown*

Add an XML comment.


* `xml_AddEntityReference`  
**entityref**: Entity reference.

This may be used anywhere that `xml_AddCharacters` may be, and will insert an entity reference into the contents of the XML document at that point. Note that if the entity inserted is a character entity, its validity well be checked according to the rules of XML-1.1, not 1.0.

If the entity reference is not a character entity, then no check is made of its validity, and a warning will be issued

### Functions to query XML file objects

These functions may be of use in building wrapper libraries:

* `xmlf_Name` result(*string*)

Return the filename of an open XML file

* `xmlf_OpenTag` result(*string*)

Return the currently open tag of the current XML file (or the empty string if none is open)

* `xmlf_GetPretty_print` result(*logical*)

Return the current value of pretty_print.

* `xmlf_SetPretty_print`
**NewValue**: *logical*

Set the current value of pretty_print to the NewValue. This may be useful in a mixed namespace
document where pretty printing the output may change the meaning under one of the namespaces.

##Exceptions

<a name="Exceptions"/>

Below are explained areas where wxml fails to implement the whole of XML 1.0/1.1. These are divided into two lists; where wxml **does not** permit the generation of a particular well-formed XML document, and where it **does** permit the generation of a particular non-well-formed document.

Ways in which wxml renders it impossible to produce a certain sort of well-formed XML document:

 1. Unicode support is limited. Due to the limitations of Fortran, wxml is unable to manipulate characters outwith 7-bit US-ASCII. wxml will ensure that characters corresponding to those in 7-bit ASCII are output correctly within the constraints of the version of XML in use, for a UTF-8 encoding. Attempts to directly output any other characters will have undefined effects.  Output of other unicode characters is possible through the use of character entities.
 1. Due to the constraints of the Fortran IO specification, it is impossible to output arbitrary long strings without carriage returns. The size of the limit varies between processors, but may be as low as 1024 characters. To avoid overrunning this limit, wxml will by default insert carriage returns before every new element, and if an unbroken string of attribute or text data is requested greater than 1024 characters, then carriage returns will be inserted as appropriate; within whitespace if possible; to ensure it is broken up into smaller sections to fit within the limits.

wxml will try very hard to ensure that output is well-formed. However, it is possible to fool wxml into producing ill-formed XML documents. Avoid doing so if possible; for completeness these ways are listed here. In all cases where ill-formedness is a possibility, a warning can be issued. These warnings can be verbose, so are off by default, but if they are desired, they can be switched on by manipulating the `warning` argument to `xml_OpenFile`.

 1. If you specify a non-default text encoding, and then run FoX on a platform which does not use this encoding, then the result will be nonsense, and more than likely ill-formed. FoX will issue a warning in this case.
 4. When adding any text, if any characters are passed in (regardless of character set) which do not have equivalants within 7-bit ASCII, then the results are processor-dependent, and may result in an invalid document on output. A warning will be issued if this occurs. If you need a guarantee that such characters will be passed correctly, use character entities.
 2. If any parameter entities are referenced, no checks are made that the document after parameter-entity-expansion is well-formed. A warning will be issued. 

###Validity constraints

Finally, note that constraints on XML documents are divided into two sets - well-formedness constraints (WFC) and validity constraints (VC). The above only applies to WFC checks. wxml can make some minimal checks on VCs, but this is by no means complete, nor is it intended to be. These checks are off by default, but may be switched on by manipulating the `validate` argument to `xml_OpenFile`.
