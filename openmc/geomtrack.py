
class ElementTracker:
    class __ElementTracker:
        def __init__(self):
            self.cells = set()
            self.surfaces = set()
            self.lattices = set()
            self.universes = set()

        def update(self, element_type, element_id):
            if element_type == "cell":
                self.cells.add(element_id)
            elif element_type == "surface":
                self.surfaces.add(element_id)
            elif element_type == "lattice":
                self.lattices.add(element_id)
            elif element_type == "universe":
                self.universes.add(element_id)

        def add_cell(self, cell_id):
            self.update("cell", cell_id)

        def add_surface(self, surface_id):
            self.update("surface", surface_id)

        def add_lattice(self, lattice_id):
            self.update("lattice", lattice_id)

        def add_universe(self, universe_id):
            self.update("universe", universe_id)

        def reset(self):
            self.cells = set()
            self.surfaces = set()
            self.lattices = set()
            self.universes = set()

    instance = None

    def __init__(self):
        if not ElementTracker.instance:
            ElementTracker.instance = ElementTracker.__ElementTracker()
        else:
            pass
    def __getattr__(self, name):
        return getattr(self.instance, name)
