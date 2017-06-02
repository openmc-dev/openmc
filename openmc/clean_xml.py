def sort_xml_elements(tree):

    # Retrieve all children of the root XML node in the tree
    elements = list(tree)

    # Initialize empty lists for the sorted and comment elements
    sorted_elements = []

    # Initialize an empty set of tags (e.g., Surface, Cell, and Lattice)
    tags = set()

    # Find the unique tags in the tree
    for element in elements:
        tags.add(element.tag)

    # Initialize an empty list for the comment elements
    comment_elements = []

    # Find the comment elements and record their ordering within the
    # tree using a precedence with respect to the subsequent nodes
    for index, element in enumerate(elements):
        next_element = None

        if 'Comment' in str(element.tag):

            if index < len(elements)-1:
                next_element = elements[index+1]

            comment_elements.append((element, next_element))

    # Now iterate over all tags and order the elements within each tag
    for tag in sorted(list(tags)):

        # Retrieve all of the elements for this tag
        try:
            tagged_elements = tree.findall(tag)
        except:
            continue

        # Initialize an empty list of tuples to sort (id, element)
        tagged_data = []

        # Retrieve the IDs for each of the elements
        for element in tagged_elements:
            key = element.get('id')

            # If this element has an "ID" tag, append it to the list to sort
            if key is not None:
                tagged_data.append((int(key), element))

        # Sort the elements according to the IDs for this tag
        tagged_data.sort()
        sorted_elements.extend(list(item[-1] for item in tagged_data))

    # Add the comment elements while preserving the original precedence
    for element, next_element in comment_elements:
        index = sorted_elements.index(next_element)
        sorted_elements.insert(index, element)

    # Remove all of the sorted elements from the tree
    for element in sorted_elements:
        tree.remove(element)

    # Add the sorted elements back to the tree in the proper order
    tree.extend(sorted_elements)


def clean_xml_indentation(element, level=0, spaces_per_level=2):
    """
    copy and paste from http://effbot.org/zone/elementent-lib.htm#prettyprint
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
            clean_xml_indentation(sub_element, level+1, spaces_per_level)

        if not sub_element.tail or not sub_element.tail.strip():
            sub_element.tail = i

    else:
        if level and (not element.tail or not element.tail.strip()):
            element.tail = i
