import openmc
import openmc.stats


def test_source():
    space = openmc.stats.Point()
    energy = openmc.stats.Discrete([1.0e6], [1.0])
    angle = openmc.stats.Isotropic()

    src = openmc.Source(space=space, angle=angle, energy=energy)
    assert src.space == space
    assert src.angle == angle
    assert src.energy == energy
    assert src.strength == 1.0

    elem = src.to_xml_element()
    assert 'strength' in elem.attrib
    assert elem.find('space') is not None
    assert elem.find('angle') is not None
    assert elem.find('energy') is not None


def test_source_file():
    filename = 'source.h5'
    src = openmc.Source(filename=filename)
    assert src.file == filename

    elem = src.to_xml_element()
    assert 'strength' in elem.attrib
    assert 'file' in elem.attrib
