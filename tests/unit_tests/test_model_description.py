import openmc

DESCRIPTION_TEXT = "This is a test model."


def test_model_description():
    """Test the description attribute on the Model class."""
    model = openmc.Model()
    model.description = DESCRIPTION_TEXT

    # Check that the description is set on the underlying settings object
    assert model.settings.description == DESCRIPTION_TEXT


def test_settings_description_xml():
    """Test the XML representation of the description in Settings."""
    settings = openmc.Settings()
    settings.description = DESCRIPTION_TEXT

    # Generate XML element
    elem = settings.to_xml_element()

    # Check for the presence and content of the description tag
    desc_elem = elem.find('description')
    assert desc_elem is not None
    assert desc_elem.text == DESCRIPTION_TEXT

    # Test from_xml_element
    new_settings = openmc.Settings.from_xml_element(elem)
    assert new_settings.description == DESCRIPTION_TEXT
