def clean_indentation(element, level=0, spaces_per_level=2):
    """
    copy and paste from https://effbot.org/zone/element-lib.htm#prettyprint
    it basically walks your tree and adds spaces and newlines so the tree is
    printed in a nice way
    """
    i = "\n" + level*spaces_per_level*" "

    if len(element):
        if not element.text or not element.text.strip():
            element.text = i + spaces_per_level*" "
        if not element.tail or not element.tail.strip():
            element.tail = i
        for sub_element in element:
            clean_indentation(sub_element, level+1, spaces_per_level)
        if not sub_element.tail or not sub_element.tail.strip():
            sub_element.tail = i
    else:
        if level and (not element.tail or not element.tail.strip()):
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
