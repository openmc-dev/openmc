from typing import get_type_hints
from inspect import getmro
from xml.etree import ElementTree as ET


def clean_indentation(element, level=0, spaces_per_level=2, trailing_indent=True):
    """Set indentation of XML element and its sub-elements.
    Copied and pasted from https://effbot.org/zone/element-lib.htm#prettyprint.
    It walks your tree and adds spaces and newlines so the tree is
    printed in a nice way.

    Parameters
    ----------
    level : int
        Indentation level for the element passed in (default 0)
    spaces_per_level : int
        Number of spaces per indentation level (default 2)
    trailing_indent : bool
        Whether or not to add indentation after closing the element

    """
    i = "\n" + level*spaces_per_level*" "

    # ensure there's always some tail for the element passed in
    if not element.tail:
        element.tail = ""

    if len(element):
        if not element.text or not element.text.strip():
            element.text = i + spaces_per_level*" "
        if trailing_indent and (not element.tail or not element.tail.strip()):
            element.tail = i
        for sub_element in element:
            # `trailing_indent` is intentionally not forwarded to the recursive
            # call. Any child element of the topmost element should add
            # indentation at the end to ensure its parent's indentation is
            # correct.
            clean_indentation(sub_element, level+1, spaces_per_level)
        if not sub_element.tail or not sub_element.tail.strip():
            sub_element.tail = i
    else:
        if trailing_indent and level and (not element.tail or not element.tail.strip()):
            element.tail = i


def get_text(elem, name, default=None):
    """Retrieve text of an attribute or subelement.

    Parameters
    ----------
    elem : xml.etree.ElementTree.Element
        Element from which to search
    name : str
        Name of attribute/subelement
    default : object
        A defult value to return if matching attribute/subelement exists

    Returns
    -------
    str
        Text of attribute or subelement

    """
    if name in elem.attrib:
        return elem.get(name, default)
    else:
        child = elem.find(name)
        return child.text if child is not None else default


def reorder_attributes(root):
    """Sort attributes in XML to preserve pre-Python 3.8 behavior

    Parameters
    ----------
    root : xml.etree.ElementTree.Element
        Root element

    """
    for el in root.iter():
        attrib = el.attrib
        if len(attrib) > 1:
            # adjust attribute order, e.g. by sorting
            attribs = sorted(attrib.items())
            attrib.clear()
            attrib.update(attribs)


def xml_repr(x):
    if hasattr(x, 'xml_repr'):
        return x.xml_repr()
    elif isinstance(x, dict):
        return x
    elif isinstance(x, list):
        return ' '.join([str(xi) for xi in x])
    elif isinstance(x, tuple):
        return ' '.join([str(xi) for xi in x])
    elif isinstance(x, bool):
        return str(x).lower()
    else:
        return str(x)


def from_xml_repr(desired_type, value):
    '''
    TODO: this could probably be improved.
    Also need to make it work on dicts.

    Notably, while our mapping from class -> XML layout
    isn't consistent, we could make it so that exports happen
    in a new way, but only maintain read compatibility for
    old-style XML.
    '''
    if hasattr(desired_type, 'from_xml_repr'):
        return desired_type.from_xml_repr(value)
    if desired_type == str:
        return value
    elif desired_type == int:
        return int(value)
    elif desired_type == float:
        return float(value)
    elif desired_type == bool:
        return value == 'true' or value == '1' or value == 'True' or value == 'y'
    elif desired_type == list[float]:
        return [float(x) for x in value.split(' ')]
    elif desired_type == list[int]:
        return [int(x) for x in value.split(' ')]
    elif desired_type == tuple[int, int, int]:
        splitvalue = value.split(' ')
        assert len(splitvalue) == 3
        return (int(splitvalue[0]), int(splitvalue[1]), int(splitvalue[2]))
    elif desired_type == tuple[int, int]:
        splitvalue = value.split(' ')
        assert len(splitvalue) == 2
        return (int(splitvalue[0], int(splitvalue[1])))
    else:
        raise Exception("Cannot convert type {} from XML text to data".format(
            desired_type))


def xmlinator(cls):
    """Class decorator that automatically adds to_xml_element
    and from_xml_element methods to a class

    _xml_name -- element names in XML file
    pre_xml_export -- define for pre-export logic
    to_xml_element_finalize -- additional modifications to perform on
      the newly created element. Can be used to handle logic not otherwise
      covered.


    Parameters
    ----------
    cls : class
        Class to be passed in
    """
    mro = getmro(cls)

    cls._xml_attributes = []
    cls._xml_optional_attributes = []
    cls._xml_elements = []
    cls._xml_optional_elements = []

    # Loop over base classes and include variables defined therein
    for cc in mro:
        if cc == object:
            continue
        if cc.__name__ in _xml_class_attributes.keys():
            cls._xml_attributes.extend(_xml_class_attributes[cc.__name__])
        if cc.__name__ in _xml_class_optional_attributes.keys():
            cls._xml_optional_attributes.extend(
                _xml_class_optional_attributes[cc.__name__])
        if cc.__name__ in _xml_class_elements.keys():
            cls._xml_elements.extend(_xml_class_elements[cc.__name__])
        if cc.__name__ in _xml_class_optional_elements.keys():
            cls._xml_optional_elements.extend(
                _xml_class_optional_elements[cc.__name__])

    # The XML element name is assumed to be the class name in lower
    # case by default. However, the developer can set it using _xml_name.
    if hasattr(cls, '_xml_name'):
        xmlname = cls._xml_name
    else:
        xmlname = cls.__name__.lower()

    # Create the to_xml_element method and stick it to the class.
    def to_xml_element(self):
        """Return XML representation of {}

        Returns
        -------
        element : xml.etree.ElementTree.Element
            XML element containing mesh data
        """.format(cls.__name__)

        # TODO run pre_xml_export from base classes too!
        # or should this be left to the user to call super(...?
        if hasattr(cls, 'pre_xml_export'):
            print('running pre xml export')
            self.pre_xml_export()

        element = ET.Element(xmlname)

        for attrname, attrtype in cls._xml_attributes:
            element.set(attrname, xml_repr(getattr(self, attrname)))

        for attrname, attrtype in cls._xml_optional_attributes:
            if not hasattr(self, '_'+attrname):
                raise Exception("You wrote code which specifies that {0} is a\
                        private member of {1}, but it does not have a member\
                        named _{0}. xmlinator requires this naming convention,\
                        sorry!".format(attrname, cls.__name__))
            value = getattr(self, '_'+attrname)
            if value is not None:
                element.set(attrname, xml_repr(value))

        for elemname, elemtype in cls._xml_elements:
            subelement = ET.SubElement(element, elemname)
            the_xml_repr = xml_repr(getattr(self, elemname))
            if isinstance(the_xml_repr, str):
                subelement.text = the_xml_repr
            elif isinstance(the_xml_repr, dict):
                for key, value in the_xml_repr.items():
                    if value is not None:
                        subelement.set(key, xml_repr(value))
            else:
                raise Exception("No valid XML representation for type \
                        {}".format(type(getattr(self, elemname))))

        for elemname, elemtype in cls._xml_optional_elements:
            if not hasattr(self, '_'+elemname):
                raise Exception("You wrote code which specifies that {0} is a\
                        private member of {1}, but it does not have a member\
                        named _{0}. xmlinator requires this naming convention,\
                        sorry!".format(elemname, cls.__name__))
            value = getattr(self, '_'+elemname)
            if value is not None:
                subelement = ET.SubElement(element, elemname)
                subelement.text = xml_repr(value)

        if hasattr(cls, 'to_xml_element_finalize'):
            self.to_xml_element_finalize(element)

        return element

    cls.to_xml_element = to_xml_element

    def from_xml_element(cls, elem):

        # avoid gettinga warning about deplicate IDs
        if hasattr(cls, 'id'):
            instance = cls(9999999999)
        else:
            instance = cls()

        for attrname, attrtype in cls._xml_attributes:
            value = elem.get(attrname, None)
            if value is None:
                raise Exception("Attribute {} is mandatory;\
                        failed to read XML file.")
            setattr(instance, attrname, from_xml_repr(attrtype, value))

        for attrname, attrtype in cls._xml_optional_attributes:
            value = elem.get(attrname, None)
            if value is not None:
                setattr(instance, attrname, from_xml_repr(attrtype, value))

        for elemname, elemtype in cls._xml_elements:
            child = elem.find(elemname)
            if child is None:
                raise Exception("Element {} is mandatory;\
                        failed to read XML file.")
            if elemtype == dict:
                # TODO TODO
                raise Exception(
                    'This is an edge case that needs to be addressed later...')
            else:
                setattr(instance, elemname, from_xml_repr(elemtype, child.text))

        for elemname, elemtype in cls._xml_optional_elements:
            child = elem.find(elemname)
            if child is None:
                continue
            if elemtype == dict:
                # TODO TODO
                raise Exception('This is an edge case to address later...')
            else:
                setattr(instance, elemname, from_xml_repr(elemtype, child.text))

        if hasattr(cls, 'from_xml_element_finalize'):
            instance.from_xml_element_finalize(elem)

        return instance

    cls.from_xml_element = classmethod(from_xml_element)

    return cls


# These are temporary variables which are used to externally
# mark variables as being xml-compatible.
_xml_class_attributes = {}
_xml_class_optional_attributes = {}
_xml_class_elements = {}
_xml_class_optional_elements = {}


def _split_qualname(qualname):
    splitname = qualname.split('.')
    if len(splitname) > 2:
        raise Exception("Cannot define XML-interfaceable classes inside of\
                 XML-interfaceable classes.")
    if len(splitname) == 1:
        raise Exception("Must use toxml inside of classes only.")
    return splitname[0], splitname[1]


def xml_attribute(method):
    '''
    Marker decorator. This doesn't alter the function that's passed in.
    It just makes note that this is a property inside the current class
    that we will export to xml in the future, specifically putting it as
    an attribute on the element.
    '''
    hints = get_type_hints(method)
    if 'return' not in hints.keys():
        raise Exception(
            "Must give type hint on XML-able property {}".format(method.__qualname__))
    rt = hints['return']
    classname, methodname = _split_qualname(method.__qualname__)
    if classname not in _xml_class_attributes.keys():
        _xml_class_attributes[classname] = []
    _xml_class_attributes[classname].append((methodname, rt))
    return method


def optional_xml_attribute(method):
    '''
    Marker decorator. This doesn't alter the function that's passed in.
    It just makes note that this is a property inside the current class
    that we will export to xml in the future, specifically putting it as
    an attribute on the element.

    This assumes that the attribute is represented by a private variable
    of the same name prefixed by an underscore. If that private value is
    None, nothing is written.
    '''
    hints = get_type_hints(method)
    if 'return' not in hints.keys():
        raise Exception(
            "Must give type hint on XML-able property {}".format(method.__qualname__))
    rt = hints['return']
    classname, methodname = _split_qualname(method.__qualname__)
    if classname not in _xml_class_optional_attributes.keys():
        _xml_class_optional_attributes[classname] = []
    _xml_class_optional_attributes[classname].append((methodname, rt))
    return method


def xml_element(method):
    '''
    Marker decorator. This doesn't alter the function that's passed in.
    It just makes note that this is a property inside the current class
    that we will export to xml in the future, specifically putting it as
    a subelement on the element.
    '''
    hints = get_type_hints(method)
    if 'return' not in hints.keys():
        raise Exception(
            "Must give type hint on XML-able property {}".format(method.__qualname__))
    rt = hints['return']
    classname, methodname = _split_qualname(method.__qualname__)
    if classname not in _xml_class_elements.keys():
        _xml_class_elements[classname] = []
    _xml_class_elements[classname].append((methodname, rt))
    return method


def optional_xml_element(method):
    '''
    Marker decorator. This doesn't alter the function that's passed in.
    It just makes note that this is a property inside the current class.
    HOWEVER, it assumes that this property is represented by a private
    variable prefixed with an underscore.

    If the variable is None, than nothing is written. Otherwise,
    it is written as a usual XML element.
    '''
    hints = get_type_hints(method)
    if 'return' not in hints.keys():
        raise Exception(
            "Must give type hint on XML-able property {}".format(method.__qualname__))
    rt = hints['return']
    classname, methodname = _split_qualname(method.__qualname__)
    if classname not in _xml_class_optional_elements.keys():
        _xml_class_optional_elements[classname] = []
    _xml_class_optional_elements[classname].append((methodname, rt))
    return method
