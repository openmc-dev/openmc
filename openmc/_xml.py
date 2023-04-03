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
    return str(x)


def xmlinator(cls):
    """Class decorator that automatically adds to_xml_element
    and from_xml_element methods to a class

    Parameters
    ----------
    cls : class
        Class to be passed in
    """
    mro = getmro(cls)

    cls._xml_attributes = []
    cls._xml_elements = []
    cls._xml_optional_elements = []

    # Loop over base classes and include variables defined therein
    for cc in mro:
        if cc == object:
            continue
        if cc.__name__ in _xml_class_attributes.keys():
            cls._xml_attributes.extend(_xml_class_attributes[cc.__name__])
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

        element = ET.Element(xmlname)

        for attrname in cls._xml_attributes:
            element.set(attrname, xml_repr(getattr(self, attrname)))

        for elemname in cls._xml_elements:
            subelement = ET.SubElement(element, elemname)
            subelement.text = xml_repr(getattr(self, elemname))

        for elemname in cls._xml_optional_elements:
            value = getattr(self, '_'+elemname)
            if value is not None:
                subelement = ET.SubElement(element, elemname)
                subelement.text = xml_repr(value)
        return element

    cls.to_xml_element = to_xml_element

    return cls


# These are temporary variables which are used to externally
# mark variables as being xml-compatible.
_xml_class_attributes = {}
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
    classname, methodname = _split_qualname(method.__qualname__)
    if classname not in _xml_class_attributes.keys():
        _xml_class_attributes[classname] = []
    _xml_class_attributes[classname].append(methodname)
    return method


def xml_element(method):
    '''
    Marker decorator. This doesn't alter the function that's passed in.
    It just makes note that this is a property inside the current class
    that we will export to xml in the future, specifically putting it as
    a subelement on the element.
    '''
    classname, methodname = _split_qualname(method.__qualname__)
    if classname not in _xml_class_elements.keys():
        _xml_class_elements[classname] = []
    _xml_class_elements[classname].append(methodname)
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
    classname, methodname = _split_qualname(method.__qualname__)
    if classname not in _xml_class_optional_elements.keys():
        _xml_class_optional_elements[classname] = []
    _xml_class_optional_elements[classname].append(methodname)
    return method
