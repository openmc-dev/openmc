# Attributes dictionaries.

When parsing XML using the FoX SAX module, attributes are returned contained within a dictionary object.

This dictionary object implements all the methods described by the SAX interfaces Attributes and Attributes2. Full documentation is available from the SAX Javadoc, but is reproduced here for ease of reference.

All of the attribute dictionary objects and functions are exported through FoX\_sax - you must USE the module to enable them. The dictionary API is described here.

An attribute dictionary consists of a list of entries, one for each attribute. The entries all have the following pieces of data:

* qName - the attribute's full name  
* value - the attribute's value

and for namespaced attributes:

* uri - the namespace URI (if any) of the attribute  
* localName - the local name of the attribute

In addition, the following pieces of data will be picked up from a DTD if present:

* declared - is the attribute declared in the DTD?  
* specified - is this instance of the attribute specified in the XML document, or is it a default from the DTD?  
* type - the type of the attribute (if declared)

------

## Derived types

There is one derived type of interest, `dictionary_t`.

It is opaque - that is, it should only be manipulated through the functions described here.

## Functions

### Inspecting the dictionary

* `getLength  
   type(dictionary_t), intent(in) :: dict`

Returns an integer with the length of the dictionary, _ie_ the number of dictionary entries.

* `hasKey  
    type(dictionary_t), intent(in) :: dict  
    character(len=*), intent(in) :: key`

Returns a logical value according to whether the dictionary contains an attribute named `key` or not.

* `hasKey  
    type(dictionary_t), intent(in) :: dict  
    character(len=*), intent(in) :: uri  
    character(len=*), intent(in) :: localname`

Returns a logical value according to whether the dictionary contains an attribute with the correct `URI` and `localname`.

### Retrieving data from the dictionary

* `getQName  
    type(dictionary_t), intent(in) :: dict   
    integer, intent(in) :: i`

Return the full name of the `i`th dictionary entry.

* `getValue  
    type(dictionary_t), intent(in)  
    integer, intent(in) :: i`

If an integer is passed in - the value of the `i`th attribute. 

* `getValue  
    type(dictionary_t), intent(in)  
    character(len=*), intent(in) :: qName`

If a single string is passed in, the value of the attribute with that name.

* `getValue  
    type(dictionary_t), intent(in)  
    character(len=*), intent(in) :: uri, localname`

If two strings are passed in, the value of the attribute with that uri and localname.

* `getURI  
    type(dictionary_t), intent(in)  
    integer, intent(in) :: i`

Returns a string containing the nsURI of the `i`th attribute.

* `getlocalName  
    type(dictionary_t), intent(in)  
    integer, intent(in) :: i`

Returns a string containing the localName of the `i`th attribute.

### DTD-driven functions

The following functions are only of interest if you are using DTDs.

* `getType 
    type(dictionary_t), intent(in)  
    integer, intent(in), optional :: i`

If an integer is passed in, returns the type of the `i`th attribute.

* `getType  
    type(dictionary_t), intent(in)  
    character(len=*), intent(in) :: qName`

If a single string is passed in, returns the type of the attribute with that QName.

* `getType  
    type(dictionary_t), intent(in)  
    character(len=*), intent(in) :: uri  
    character(len=*), intent(in) :: localName`

If a single string is passed in, returnsthe type of the attribute with that {uri,localName}.

* `isDeclared  
    type(dictionary_t), intent(in)  
    integer, intent(in), optional :: i`

If an integer is passed in, returns false unless the `i`th attribute is declared in the DTD.

* `isDeclared 
    type(dictionary_t), intent(in)  
    character(len=*), intent(in) :: qName`

If a single string is passed in, returns false unless the attribute with that QName is declared in the DTD.

* `isDeclared  
    type(dictionary_t), intent(in)  
    character(len=*), intent(in) :: uri  
    character(len=*), intent(in) :: localName`

If a single string is passed in, returns false unless the attribute with that {uri,localName} is declared in the DTD.

* `isSpecified  
    type(dictionary_t), intent(in)  
    integer, intent(in), optional :: i`

If an integer is passed in, returns true unless the `i`th attribute is a default value from the DTD.

* `isSpecified 
    type(dictionary_t), intent(in)  
    character(len=*), intent(in) :: qName`

If a single string is passed in, returns true unless the attribute with that QName is a default value from the DTD.

* `isSpecified  
    type(dictionary_t), intent(in)  
    character(len=*), intent(in) :: uri  
    character(len=*), intent(in) :: localName`

If a single string is passed in, returns true unless the attribute with that {uri,localName} is a default value from the DTD.
