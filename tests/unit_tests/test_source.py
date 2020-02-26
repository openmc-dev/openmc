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

    elem = src.to_xml_element()
    assert 'strength' in elem.attrib
    assert elem.find('space') is not None
    assert elem.find('angle') is not None
    assert elem.find('energy') is not None

    src = openmc.Source.from_xml_element(elem)
    assert isinstance(src.angle, openmc.stats.Isotropic)
    assert src.space.xyz == [0.0, 0.0, 0.0]
    assert src.energy.x == [1.0e6]
    assert src.energy.p == [1.0]
    assert src.strength == 1.0


def test_source_file():
    filename = 'source.h5'
    src = openmc.Source(filename=filename)
    assert src.file == filename

    elem = src.to_xml_element()
    assert 'strength' in elem.attrib
    assert 'file' in elem.attrib

def test_source_dlopen():
    library = './libsource.so'
    src = openmc.Source(library=library)
    assert src.library == library

    elem = src.to_xml_element()
    assert 'library' in elem.attrib
