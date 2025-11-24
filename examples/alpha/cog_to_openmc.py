#!/usr/bin/env python3
"""
COG to OpenMC Python Converter
Translates COG Monte Carlo input files to OpenMC Python scripts.

This converter produces a Python script that creates:
- Materials with proper nuclide/element definitions
- Geometry with surfaces, cells, universes, and lattices
- Settings with run parameters and source definition

Author: William Zywiec
"""

import re
from pathlib import Path
from dataclasses import dataclass, field
from typing import Optional


# =============================================================================
# Element Data
# =============================================================================

ATOMIC_MASSES = {
    'H': 1.008, 'He': 4.003, 'Li': 6.941, 'Be': 9.012, 'B': 10.81,
    'C': 12.01, 'N': 14.01, 'O': 16.00, 'F': 19.00, 'Ne': 20.18,
    'Na': 22.99, 'Mg': 24.31, 'Al': 26.98, 'Si': 28.09, 'P': 30.97,
    'S': 32.07, 'Cl': 35.45, 'Ar': 39.95, 'K': 39.10, 'Ca': 40.08,
    'Sc': 44.96, 'Ti': 47.87, 'V': 50.94, 'Cr': 52.00, 'Mn': 54.94,
    'Fe': 55.85, 'Co': 58.93, 'Ni': 58.69, 'Cu': 63.55, 'Zn': 65.38,
    'Ga': 69.72, 'Ge': 72.63, 'As': 74.92, 'Se': 78.97, 'Br': 79.90,
    'Kr': 84.80, 'Rb': 85.47, 'Sr': 87.62, 'Y': 88.91, 'Zr': 91.22,
    'Nb': 92.91, 'Mo': 95.95, 'Tc': 98.00, 'Ru': 101.1, 'Rh': 102.9,
    'Pd': 106.4, 'Ag': 107.9, 'Cd': 112.4, 'In': 114.8, 'Sn': 118.7,
    'Sb': 121.8, 'Te': 127.6, 'I': 126.9, 'Xe': 131.3, 'Cs': 132.9,
    'Ba': 137.3, 'La': 138.9, 'Ce': 140.1, 'Pr': 140.9, 'Nd': 144.2,
    'Pm': 145.0, 'Sm': 150.4, 'Eu': 152.0, 'Gd': 157.3, 'Tb': 158.9,
    'Dy': 162.5, 'Ho': 164.9, 'Er': 167.3, 'Tm': 168.9, 'Yb': 173.0,
    'Lu': 175.0, 'Hf': 178.5, 'Ta': 180.9, 'W': 183.8, 'Re': 186.2,
    'Os': 190.2, 'Ir': 192.2, 'Pt': 195.1, 'Au': 197.0, 'Hg': 200.6,
    'Tl': 204.4, 'Pb': 207.2, 'Bi': 209.0, 'Po': 209.0, 'At': 210.0,
    'Rn': 222.0, 'Fr': 223.0, 'Ra': 226.0, 'Ac': 227.0, 'Th': 232.0,
    'Pa': 231.0, 'U': 238.0, 'Np': 237.0, 'Pu': 244.0, 'Am': 243.0,
    'Cm': 247.0, 'Bk': 247.0, 'Cf': 251.0, 'Es': 252.0, 'Fm': 257.0,
}

ELEMENTS = set(ATOMIC_MASSES.keys())


# =============================================================================
# Data Classes
# =============================================================================

@dataclass
class Nuclide:
    name: str
    density: float
    is_element: bool
    sab: Optional[str] = None


@dataclass
class Material:
    mat_id: int
    name: str
    nuclides: list = field(default_factory=list)
    sab_tables: list = field(default_factory=list)
    is_fissile: bool = False


@dataclass
class Surface:
    surf_id: int
    surf_type: str
    params: list
    boundary: Optional[str] = None
    comment: str = ""


@dataclass
class Cell:
    cell_id: int
    name: str
    material_id: Optional[int]
    region: str
    universe_id: Optional[int] = None
    fill_universe: Optional[int] = None
    translation: Optional[tuple] = None


@dataclass
class Lattice:
    lattice_id: int
    dimension: tuple
    lower_left: tuple
    pitch: tuple
    universes: list


@dataclass
class Universe:
    universe_id: int
    cells: list = field(default_factory=list)
    lattice: Optional[Lattice] = None


# =============================================================================
# COG Parser
# =============================================================================

class COGParser:
    """Parser for COG input files."""

    def __init__(self, filename):
        self.filename = filename
        self.title = ""
        self.materials = {}
        self.surfaces = {}
        self.cells = []
        self.universes = {}
        self.lattices = {}
        self.boundary_surfaces = {}

        self.npart = 10000
        self.nbatch = 100
        self.nfirst = 20
        self.source_points = [(0.0, 0.0, 0.0)]

        self._current_section = None
        self._current_unit = None
        self._cell_counter = 1
        self._lines = []
        self._line_idx = 0
        self._last_material_id = None

    def parse(self):
        """Parse the COG input file."""
        with open(self.filename, 'r') as f:
            self._lines = f.readlines()

        if self._lines:
            self.title = self._lines[0].strip()

        self._line_idx = 0
        while self._line_idx < len(self._lines):
            line = self._lines[self._line_idx].strip()

            if not line or line.startswith('$'):
                self._line_idx += 1
                continue

            lower_line = line.lower()
            if lower_line == 'basic':
                self._current_section = 'basic'
            elif lower_line == 'criticality':
                self._current_section = 'criticality'
            elif lower_line.startswith('mix'):
                self._current_section = 'mix'
            elif lower_line == 'geometry':
                self._current_section = 'geometry'
            elif lower_line == 'surfaces':
                self._current_section = 'surfaces'
            elif lower_line.startswith('assign-mc'):
                self._current_section = 'assign-mc'
            elif lower_line in ['end', 'end.']:
                break
            elif self._current_section:
                if self._current_section == 'geometry' and lower_line.startswith('define unit'):
                    self._parse_define_unit()
                    continue
                else:
                    self._parse_section_line(line)

            self._line_idx += 1

        self._identify_boundary_surfaces()

    def _parse_section_line(self, line):
        if self._current_section == 'criticality':
            self._parse_criticality(line)
        elif self._current_section == 'mix':
            self._parse_mix(line)
        elif self._current_section == 'geometry':
            self._parse_geometry(line)
        elif self._current_section == 'surfaces':
            self._parse_surface(line)

    def _parse_criticality(self, line):
        lower = line.lower()

        match = re.search(r'npart=(\d+)', lower)
        if match:
            self.npart = int(match.group(1))

        match = re.search(r'nbatch=(\d+)', lower)
        if match:
            self.nbatch = int(match.group(1))

        match = re.search(r'nfirst=(\d+)', lower)
        if match:
            self.nfirst = int(match.group(1))

        match = re.search(r'nsource=(\d+)', lower)
        if match:
            nsource = int(match.group(1))
            coords_text = line[match.end():]
            all_coords = []
            nums = re.findall(r'[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?', coords_text)
            all_coords.extend([float(n) for n in nums])

            while len(all_coords) < nsource * 3:
                self._line_idx += 1
                if self._line_idx >= len(self._lines):
                    break
                next_line = self._lines[self._line_idx].strip()
                if not next_line or next_line.startswith('$'):
                    continue
                if any(kw in next_line.lower() for kw in ['npart', 'nbatch', 'mix', 'assign', 'geometry']):
                    self._line_idx -= 1
                    break
                nums = re.findall(r'[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?', next_line)
                all_coords.extend([float(n) for n in nums])

            self.source_points = []
            for i in range(0, min(len(all_coords), nsource * 3), 3):
                if i + 2 < len(all_coords):
                    self.source_points.append((all_coords[i], all_coords[i+1], all_coords[i+2]))

            if not self.source_points:
                self.source_points = [(0.0, 0.0, 0.0)]

    def _parse_mix(self, line):
        lower = line.lower()
        if 'mat=' not in lower:
            if self._last_material_id is not None and self._last_material_id in self.materials:
                material = self.materials[self._last_material_id]
                self._parse_bunches_material(line, material, is_continuation=True)
            return

        mat_match = re.search(r'mat=(\d+)', lower)
        if not mat_match:
            return

        mat_id = int(mat_match.group(1))
        self._last_material_id = mat_id

        comment = ""
        if '$' in line:
            comment = line.split('$')[1].strip()
            line = line.split('$')[0]

        material = Material(mat_id=mat_id, name=comment)

        if 'a-f' in lower:
            self._parse_atom_fraction_material(line, material)
        elif 'w-f' in lower:
            self._parse_weight_fraction_material(line, material)
        else:
            self._parse_bunches_material(line, material)

        self.materials[mat_id] = material

    def _parse_bunches_material(self, line, material, is_continuation=False):
        if is_continuation:
            content = line
        else:
            content = re.sub(r'mat=\d+\s+bunches\s*', '', line, flags=re.IGNORECASE)

        pattern = r'\(?([\w.]+)\)?\s+([-+]?\d*\.?\d+(?:[-+eE]\d+)?)'
        matches = re.findall(pattern, content)

        for name, density_str in matches:
            nuclide = self._convert_nuclide(name, density_str)
            if nuclide:
                material.nuclides.append(nuclide)
                if nuclide.sab and nuclide.sab not in material.sab_tables:
                    material.sab_tables.append(nuclide.sab)
                if any(f in nuclide.name for f in ['U233', 'U235', 'Pu239', 'Pu241']):
                    material.is_fissile = True

    def _parse_atom_fraction_material(self, line, material):
        if '$' in line:
            line = line.split('$')[0]

        match = re.search(r'a-f\s+([\d.]+)\s+(.+)', line, re.IGNORECASE)
        if not match:
            return

        mass_density = float(match.group(1))
        content = match.group(2).strip()

        pattern = r'([a-zA-Z]+\d*)\s+([\d.]+)'
        matches = re.findall(pattern, content)

        if not matches:
            return

        total_at_pct = sum(float(pct) for _, pct in matches)
        if total_at_pct == 0:
            return

        avg_mass = 0.0
        for elem_or_isotope, pct in matches:
            at_frac = float(pct) / total_at_pct
            iso_match = re.match(r'([a-zA-Z]+)(\d*)', elem_or_isotope)
            if iso_match:
                elem = iso_match.group(1).capitalize()
                mass_num = iso_match.group(2)
                if mass_num:
                    avg_mass += at_frac * float(mass_num)
                elif elem in ATOMIC_MASSES:
                    avg_mass += at_frac * ATOMIC_MASSES[elem]

        if avg_mass > 0:
            total_atom_density = mass_density * 0.6022 / avg_mass

            for elem_or_isotope, pct in matches:
                at_frac = float(pct) / total_at_pct
                atom_density = total_atom_density * at_frac

                nuclide = self._convert_nuclide(elem_or_isotope, str(atom_density))
                if nuclide:
                    material.nuclides.append(nuclide)
                    if any(f in nuclide.name for f in ['U233', 'U235', 'Pu239', 'Pu241']):
                        material.is_fissile = True

    def _parse_weight_fraction_material(self, line, material):
        if '$' in line:
            line = line.split('$')[0]

        match = re.search(r'w-f\s+([\d.]+)\s+(.+)', line, re.IGNORECASE)
        if not match:
            return

        mass_density = float(match.group(1))
        content = match.group(2).strip()

        pattern = r'([a-zA-Z]+\d*)\s+([\d.]+)'
        matches = re.findall(pattern, content)

        if not matches:
            return

        total_wt_pct = sum(float(pct) for _, pct in matches)
        if total_wt_pct == 0:
            return

        for elem_or_isotope, pct in matches:
            wt_frac = float(pct) / total_wt_pct
            iso_match = re.match(r'([a-zA-Z]+)(\d*)', elem_or_isotope)
            if iso_match:
                elem = iso_match.group(1).capitalize()
                mass_num = iso_match.group(2)
                if mass_num:
                    atomic_mass = float(mass_num)
                elif elem in ATOMIC_MASSES:
                    atomic_mass = ATOMIC_MASSES[elem]
                else:
                    continue

                atom_density = mass_density * wt_frac * 0.6022 / atomic_mass
                nuclide = self._convert_nuclide(elem_or_isotope, str(atom_density))
                if nuclide:
                    material.nuclides.append(nuclide)
                    if any(f in nuclide.name for f in ['U233', 'U235', 'Pu239', 'Pu241']):
                        material.is_fissile = True

    def _convert_nuclide(self, cog_name, density_str):
        density_str = density_str.strip()
        if density_str and density_str[0].isdigit():
            if re.match(r'[\d.]+[-+]\d+$', density_str):
                density_str = re.sub(r'([\d.]+)([-+])(\d+)$', r'\1e\2\3', density_str)
        try:
            density = float(density_str)
        except ValueError:
            return None

        cog_lower = cog_name.lower().strip('()')
        sab = None
        is_element = False

        if 'h.h2o' in cog_lower or cog_lower == 'h2o':
            return Nuclide(name='H1', density=density, is_element=False, sab='c_H_in_H2O')

        if 'h.ch2' in cog_lower or cog_lower == 'ch2':
            return Nuclide(name='H1', density=density, is_element=False, sab='c_H_in_CH2')

        if cog_lower in ['c.graphite', 'c30p']:
            return Nuclide(name='C', density=density, is_element=True, sab='c_Graphite')

        if cog_lower == 'c':
            return Nuclide(name='C', density=density, is_element=True, sab=None)

        match = re.match(r'([a-z]+)(\d*)', cog_lower)
        if not match:
            return None

        elem = match.group(1).capitalize()
        mass = match.group(2)

        if elem not in ELEMENTS:
            return None

        if mass:
            name = f'{elem}{mass}'
            is_element = False
        else:
            name = elem
            is_element = True

        return Nuclide(name=name, density=density, is_element=is_element, sab=sab)

    def _parse_geometry(self, line):
        lower = line.lower()

        if any(cmd in lower for cmd in ['picture', 'volume']):
            return

        if lower.startswith('use unit'):
            self._parse_use_unit(line)
        elif 'sector' in lower:
            self._parse_sector(line, universe_id=0)
        elif lower.startswith('boundary'):
            self._parse_boundary(line)

    def _parse_sector(self, line, universe_id):
        match = re.search(r'sector\s+(\d+)\s+(\S+)\s+(.*)', line, re.IGNORECASE)
        if not match:
            return

        mat_id = int(match.group(1))
        name = match.group(2)
        surf_expr = match.group(3).strip()

        if '$' in surf_expr:
            surf_expr = surf_expr.split('$')[0].strip()

        cell = Cell(
            cell_id=self._cell_counter,
            name=name,
            material_id=mat_id,
            region=surf_expr,
            universe_id=universe_id
        )

        self._cell_counter += 1

        if universe_id == 0:
            self.cells.append(cell)
        else:
            if universe_id not in self.universes:
                self.universes[universe_id] = Universe(universe_id=universe_id)
            self.universes[universe_id].cells.append(cell)

    def _parse_use_unit(self, line):
        match = re.search(r'use\s+unit\s+(\d+)\s+(\S+)\s+(.*)', line, re.IGNORECASE)
        if not match:
            return

        unit_id = int(match.group(1))
        name = match.group(2)
        rest = match.group(3).strip()

        translation = None
        tr_match = re.search(r'tr\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)', rest, re.IGNORECASE)
        if tr_match:
            translation = (float(tr_match.group(1)), float(tr_match.group(2)), float(tr_match.group(3)))
            rest = rest[:tr_match.start()].strip()

        if '$' in rest:
            rest = rest.split('$')[0].strip()

        cell = Cell(
            cell_id=self._cell_counter,
            name=name,
            material_id=None,
            region=rest,
            universe_id=0,
            fill_universe=unit_id,
            translation=translation
        )

        self._cell_counter += 1
        self.cells.append(cell)

    def _parse_boundary(self, line):
        parts = line.split()
        if len(parts) < 3:
            return

        if parts[1].lower() in ['vacuum', 'reflecting', 'reflective']:
            boundary_type = parts[1].lower()
            if boundary_type == 'reflective':
                boundary_type = 'reflecting'
            try:
                surf_id = int(parts[2])
                self.boundary_surfaces[surf_id] = boundary_type
            except ValueError:
                pass
        else:
            try:
                surf_id = int(parts[1])
                boundary_type = parts[2].lower() if len(parts) > 2 else 'vacuum'

                if boundary_type == 'periodic' and len(parts) >= 4:
                    paired_surf = int(parts[3].lstrip('-'))
                    self.boundary_surfaces[surf_id] = f'periodic:{paired_surf}'
                elif boundary_type == 'reflective':
                    self.boundary_surfaces[surf_id] = 'reflecting'
                else:
                    self.boundary_surfaces[surf_id] = boundary_type
            except ValueError:
                pass

    def _parse_define_unit(self):
        line = self._lines[self._line_idx].strip()

        match = re.search(r'define\s+unit\s+(\d+)', line, re.IGNORECASE)
        if not match:
            self._line_idx += 1
            return

        unit_id = int(match.group(1))
        self._current_unit = unit_id

        if unit_id not in self.universes:
            self.universes[unit_id] = Universe(universe_id=unit_id)

        self._line_idx += 1
        lattice_lines = []
        in_lattice = False

        while self._line_idx < len(self._lines):
            line = self._lines[self._line_idx].strip()

            if not line or line.startswith('$'):
                self._line_idx += 1
                continue

            lower = line.lower()

            if (lower.startswith('define unit') or
                lower.startswith('picture') or
                lower.startswith('volume') or
                lower == 'surfaces' or
                lower in ['end', 'end.']):
                self._line_idx -= 1
                break

            if '{' in line:
                in_lattice = True
                lattice_lines.append(line)
            elif in_lattice:
                lattice_lines.append(line)
                if '}' in line:
                    self._parse_lattice(unit_id, lattice_lines)
                    in_lattice = False
                    lattice_lines = []
            elif 'sector' in lower:
                self._parse_sector(line, unit_id)

            self._line_idx += 1

        self._current_unit = None

    def _parse_lattice(self, unit_id, lines):
        full_text = ' '.join(lines)
        match = re.search(r'\{(.*?)\}', full_text, re.DOTALL)
        if not match:
            return

        content = match.group(1)
        tokens = content.split()

        x_data = None
        y_data = None
        fill_pattern = []

        i = 0
        while i < len(tokens):
            token = tokens[i].lower()

            if token == 'x' and i + 4 <= len(tokens):
                try:
                    x_data = {
                        'n': int(tokens[i+2]),
                        'range': (float(tokens[i+3].strip('[]')), float(tokens[i+4].strip('[]')))
                    }
                except (ValueError, IndexError):
                    pass
                i += 5
            elif token == 'y' and i + 4 <= len(tokens):
                try:
                    y_data = {
                        'n': int(tokens[i+2]),
                        'range': (float(tokens[i+3].strip('[]')), float(tokens[i+4].strip('[]')))
                    }
                except (ValueError, IndexError):
                    pass
                i += 5
            elif token == 'z' and i + 2 <= len(tokens):
                i += 4 if i + 3 < len(tokens) and not tokens[i+3].lower() == 'fill' else 3
            elif token == 'fill':
                i += 1
                while i < len(tokens):
                    fill_token = tokens[i]
                    if '*' in fill_token:
                        try:
                            parts = fill_token.split('*')
                            count = int(parts[0])
                            unit = int(parts[1])
                            fill_pattern.extend([unit] * count)
                        except (ValueError, IndexError):
                            pass
                    else:
                        try:
                            fill_pattern.append(int(fill_token))
                        except ValueError:
                            pass
                    i += 1
                break
            else:
                i += 1

        if x_data and y_data:
            nx = x_data['n'] - 1
            ny = y_data['n'] - 1
            xmin, xmax = x_data['range']
            ymin, ymax = y_data['range']

            pitch_x = (xmax - xmin) / nx if nx > 0 else 1.0
            pitch_y = (ymax - ymin) / ny if ny > 0 else 1.0

            universes = []
            idx = 0
            for j in range(ny):
                row = []
                for k in range(nx):
                    if idx < len(fill_pattern):
                        row.append(fill_pattern[idx])
                    else:
                        row.append(fill_pattern[0] if fill_pattern else 1)
                    idx += 1
                universes.append(row)

            lattice = Lattice(
                lattice_id=unit_id,
                dimension=(nx, ny),
                lower_left=(xmin, ymin),
                pitch=(pitch_x, pitch_y),
                universes=universes
            )

            self.universes[unit_id].lattice = lattice

    def _parse_surface(self, line):
        comment = ""
        if '$' in line:
            comment = line.split('$')[1].strip()
            line = line.split('$')[0].strip()

        parts = line.split()
        if len(parts) < 2:
            return

        try:
            surf_id = int(parts[0])
        except ValueError:
            return

        surf_type = parts[1].lower()
        params = parts[2:]

        # Handle multi-line prism definitions
        if surf_type == 'prism' and len(params) > 0:
            n_vertices = int(params[0])
            needed_coords = n_vertices * 2
            all_params = params[1:]

            # Keep reading lines until we have all coordinates
            while len([p for p in all_params if self._is_number(p)]) < needed_coords:
                self._line_idx += 1
                if self._line_idx >= len(self._lines):
                    break
                next_line = self._lines[self._line_idx].strip()
                if not next_line or next_line.startswith('$'):
                    continue
                if '$' in next_line:
                    next_line = next_line.split('$')[0].strip()
                # Check if this is a new surface definition
                next_parts = next_line.split()
                if next_parts and next_parts[0].isdigit() and len(next_parts) > 1:
                    if next_parts[1].lower() in ['sphere', 'sph', 'so', 'c', 'cyl', 'cylinder',
                                                   'px', 'py', 'pz', 'p', 'plane', 'pla',
                                                   'box', 'rpp', 'c/x', 'c/y', 'c/z', 'prism']:
                        self._line_idx -= 1
                        break
                all_params.extend(next_line.split())

            params = [str(n_vertices)] + all_params

        surface = Surface(
            surf_id=surf_id,
            surf_type=surf_type,
            params=params,
            comment=comment
        )

        self.surfaces[surf_id] = surface

    def _is_number(self, s):
        try:
            float(s)
            return True
        except ValueError:
            return False

    def _identify_boundary_surfaces(self):
        for surf_id, boundary_type in self.boundary_surfaces.items():
            if surf_id in self.surfaces:
                self.surfaces[surf_id].boundary = boundary_type


# =============================================================================
# OpenMC Python Generator
# =============================================================================

class OpenMCPythonGenerator:
    """Generate OpenMC Python script from parsed COG data."""

    def __init__(self, parser):
        self.parser = parser
        self._prism_planes = {}  # Maps prism surf_id to list of plane variable names

    def generate(self):
        """Generate OpenMC Python script."""
        lines = []

        # Header
        lines.append('"""')
        lines.append(f'{self.parser.title}')
        lines.append('Converted from COG to OpenMC')
        lines.append('"""')
        lines.append('')
        lines.append('import openmc')
        lines.append('')

        # Materials
        lines.append('# ' + '=' * 78)
        lines.append('# Materials')
        lines.append('# ' + '=' * 78)
        lines.append('')
        lines.extend(self._generate_materials())

        # Geometry
        lines.append('# ' + '=' * 78)
        lines.append('# Geometry')
        lines.append('# ' + '=' * 78)
        lines.append('')
        lines.extend(self._generate_surfaces())
        lines.append('')
        lines.extend(self._generate_universes())
        lines.extend(self._generate_cells())

        # Settings
        lines.append('# ' + '=' * 78)
        lines.append('# Settings')
        lines.append('# ' + '=' * 78)
        lines.append('')
        lines.extend(self._generate_settings())

        # Export
        lines.append('# ' + '=' * 78)
        lines.append('# Export and Run')
        lines.append('# ' + '=' * 78)
        lines.append('')
        lines.append('materials.export_to_xml()')
        lines.append('geometry.export_to_xml()')
        lines.append('settings.export_to_xml()')
        lines.append('')

        return '\n'.join(lines)

    def _generate_materials(self):
        lines = []

        for mat_id in sorted(self.parser.materials.keys()):
            mat = self.parser.materials[mat_id]
            var_name = f'mat{mat_id}'

            if mat.name:
                lines.append(f'# {mat.name}')

            lines.append(f'{var_name} = openmc.Material(material_id={mat_id})')
            lines.append(f'{var_name}.set_density("sum")')

            for nuclide in mat.nuclides:
                if nuclide.is_element:
                    lines.append(f'{var_name}.add_element("{nuclide.name}", {nuclide.density:.6e})')
                else:
                    lines.append(f'{var_name}.add_nuclide("{nuclide.name}", {nuclide.density:.6e})')

            for sab in mat.sab_tables:
                lines.append(f'{var_name}.add_s_alpha_beta("{sab}")')

            lines.append('')

        mat_list = ', '.join(f'mat{i}' for i in sorted(self.parser.materials.keys()))
        lines.append(f'materials = openmc.Materials([{mat_list}])')
        lines.append('')

        return lines

    def _generate_surfaces(self):
        lines = []

        for surf_id in sorted(self.parser.surfaces.keys()):
            surface = self.parser.surfaces[surf_id]
            surf_code = self._convert_surface(surface)
            if surf_code:
                if surface.comment:
                    lines.append(f'# {surface.comment}')
                lines.append(surf_code)

        return lines

    def _convert_surface(self, surface):
        surf_type = surface.surf_type.lower()
        params = surface.params
        surf_id = surface.surf_id
        var_name = f'surf{surf_id}'

        bc = ''
        if surface.boundary:
            if surface.boundary.startswith('periodic:'):
                bc = ', boundary_type="periodic"'
            elif surface.boundary in ['vacuum', 'reflecting']:
                bc = f', boundary_type="{surface.boundary}"'

        try:
            if surf_type in ['sphere', 'sph', 'so']:
                return self._gen_sphere(var_name, surf_id, params, bc)
            elif surf_type in ['c', 'cyl', 'cylinder']:
                return self._gen_cylinder(var_name, surf_id, params, bc)
            elif surf_type in ['px']:
                return f'{var_name} = openmc.XPlane(surface_id={surf_id}, x0={params[0]}{bc})'
            elif surf_type in ['py']:
                return f'{var_name} = openmc.YPlane(surface_id={surf_id}, y0={params[0]}{bc})'
            elif surf_type in ['pz']:
                return f'{var_name} = openmc.ZPlane(surface_id={surf_id}, z0={params[0]}{bc})'
            elif surf_type in ['p', 'plane', 'pla']:
                return self._gen_plane(var_name, surf_id, params, bc)
            elif surf_type in ['box']:
                return self._gen_box(var_name, surf_id, params, bc)
            elif surf_type in ['rpp']:
                return self._gen_rpp(var_name, surf_id, params, bc)
            elif surf_type in ['c/x']:
                return self._gen_axial_cyl(var_name, surf_id, 'x', params, bc)
            elif surf_type in ['c/y']:
                return self._gen_axial_cyl(var_name, surf_id, 'y', params, bc)
            elif surf_type in ['c/z']:
                return self._gen_axial_cyl(var_name, surf_id, 'z', params, bc)
            elif surf_type == 'prism':
                return self._gen_prism(var_name, surf_id, params, bc)
            else:
                return f'# {var_name}: Unsupported surface type "{surf_type}" with params {params}'
        except (ValueError, IndexError):
            return f'# {var_name}: Error converting surface type "{surf_type}"'

    def _gen_sphere(self, var_name, surf_id, params, bc):
        if len(params) == 1:
            return f'{var_name} = openmc.Sphere(surface_id={surf_id}, r={params[0]}{bc})'
        elif len(params) >= 4:
            return f'{var_name} = openmc.Sphere(surface_id={surf_id}, x0={params[0]}, y0={params[1]}, z0={params[2]}, r={params[3]}{bc})'
        else:
            return f'{var_name} = openmc.Sphere(surface_id={surf_id}, r={params[0]}{bc})'

    def _gen_cylinder(self, var_name, surf_id, params, bc):
        if not params:
            return None

        axis = params[0].lower() if params[0].lower() in ['x', 'y', 'z'] else 'z'
        idx = 1 if params[0].lower() in ['x', 'y', 'z'] else 0

        cyl_class = {'x': 'XCylinder', 'y': 'YCylinder', 'z': 'ZCylinder'}[axis]
        remaining = params[idx:]

        if not remaining:
            return None

        radius = float(remaining[0])

        if len(remaining) >= 5:
            x0, y0 = remaining[1], remaining[2]
            return f'{var_name} = openmc.{cyl_class}(surface_id={surf_id}, x0={x0}, y0={y0}, r={radius}{bc})'
        else:
            return f'{var_name} = openmc.{cyl_class}(surface_id={surf_id}, r={radius}{bc})'

    def _gen_axial_cyl(self, var_name, surf_id, axis, params, bc):
        cyl_class = {'x': 'XCylinder', 'y': 'YCylinder', 'z': 'ZCylinder'}[axis]
        if len(params) >= 3:
            return f'{var_name} = openmc.{cyl_class}(surface_id={surf_id}, x0={params[0]}, y0={params[1]}, r={params[2]}{bc})'
        elif params:
            return f'{var_name} = openmc.{cyl_class}(surface_id={surf_id}, r={params[0]}{bc})'
        return None

    def _gen_plane(self, var_name, surf_id, params, bc):
        if not params:
            return None
        if params[0].lower() in ['x', 'y', 'z']:
            axis = params[0].lower()
            value = params[1] if len(params) > 1 else '0'
            plane_class = {'x': 'XPlane', 'y': 'YPlane', 'z': 'ZPlane'}[axis]
            coord = {'x': 'x0', 'y': 'y0', 'z': 'z0'}[axis]
            return f'{var_name} = openmc.{plane_class}(surface_id={surf_id}, {coord}={value}{bc})'
        elif len(params) >= 4:
            return f'{var_name} = openmc.Plane(surface_id={surf_id}, a={params[0]}, b={params[1]}, c={params[2]}, d={params[3]}{bc})'
        return None

    def _gen_box(self, var_name, surf_id, params, bc):
        if len(params) < 3:
            return None

        dx, dy, dz = float(params[0]), float(params[1]), float(params[2])
        x0, y0, z0 = 0.0, 0.0, 0.0

        for i, p in enumerate(params):
            if str(p).lower() == 'tr' and i + 3 < len(params):
                x0, y0, z0 = float(params[i+1]), float(params[i+2]), float(params[i+3])
                break

        xmin, xmax = x0 - dx/2, x0 + dx/2
        ymin, ymax = y0 - dy/2, y0 + dy/2
        zmin, zmax = z0 - dz/2, z0 + dz/2

        return f'{var_name} = openmc.model.RectangularParallelepiped({xmin}, {xmax}, {ymin}, {ymax}, {zmin}, {zmax}{bc})'

    def _gen_rpp(self, var_name, surf_id, params, bc):
        if len(params) >= 6:
            return f'{var_name} = openmc.model.RectangularParallelepiped({params[0]}, {params[1]}, {params[2]}, {params[3]}, {params[4]}, {params[5]}{bc})'
        return None

    def _gen_prism(self, var_name, surf_id, params, bc):
        """Generate planes for a prism (polygon extruded in Z)."""
        lines = []

        # Parse prism: n_vertices x1 y1 x2 y2 ... xn yn [tr tx ty tz ...]
        n_vertices = int(params[0])
        coords = []
        tr_offset = (0.0, 0.0)

        i = 1
        while i < len(params) and len(coords) < n_vertices * 2:
            if params[i].lower() == 'tr':
                if i + 2 < len(params):
                    tr_offset = (float(params[i+1]), float(params[i+2]))
                break
            try:
                coords.append(float(params[i]))
            except ValueError:
                pass
            i += 1

        # Extract vertices
        vertices = []
        for j in range(0, len(coords) - 1, 2):
            vertices.append((coords[j] + tr_offset[0], coords[j+1] + tr_offset[1]))

        if len(vertices) < 3:
            return f'# {var_name}: Prism needs at least 3 vertices'

        # Calculate centroid for inside/outside determination
        cx = sum(v[0] for v in vertices) / len(vertices)
        cy = sum(v[1] for v in vertices) / len(vertices)

        lines.append(f'# Prism {surf_id}: {len(vertices)}-sided polygon')

        # Generate a plane for each edge
        plane_vars = []
        for j in range(len(vertices)):
            x1, y1 = vertices[j]
            x2, y2 = vertices[(j + 1) % len(vertices)]

            # Plane equation: a*x + b*y = d
            a = y2 - y1
            b = x1 - x2
            d = a * x1 + b * y1

            # Normalize
            length = (a*a + b*b) ** 0.5
            if length > 0:
                a, b, d = a/length, b/length, d/length

            plane_var = f'{var_name}_{j}'
            plane_vars.append(plane_var)

            # Determine sign: centroid should be on negative side
            centroid_val = a * cx + b * cy - d
            if centroid_val > 0:
                a, b, d = -a, -b, -d

            lines.append(f'{plane_var} = openmc.Plane(a={a:.10f}, b={b:.10f}, c=0, d={d:.10f})')

        # Store the component planes for region construction
        self._prism_planes[surf_id] = plane_vars

        return '\n'.join(lines)

    def _generate_universes(self):
        lines = []

        if not self.parser.universes:
            return lines

        lines.append('# ' + '-' * 78)
        lines.append('# Universes')
        lines.append('# ' + '-' * 78)
        lines.append('')

        for unit_id in sorted(self.parser.universes.keys()):
            universe = self.parser.universes[unit_id]

            if universe.lattice:
                lines.extend(self._gen_lattice_universe(unit_id, universe))
            else:
                lines.extend(self._gen_simple_universe(unit_id, universe))

            lines.append('')

        return lines

    def _gen_simple_universe(self, unit_id, universe):
        lines = []
        cell_vars = []

        for i, cell in enumerate(universe.cells):
            cell_var = f'u{unit_id}_cell{i}'
            cell_vars.append(cell_var)

            region = self._convert_region(cell.region)

            if cell.material_id is not None and cell.material_id in self.parser.materials:
                lines.append(f'{cell_var} = openmc.Cell(fill=mat{cell.material_id})')
            else:
                lines.append(f'{cell_var} = openmc.Cell()')

            if region:
                lines.append(f'{cell_var}.region = {region}')

        lines.append(f'universe{unit_id} = openmc.Universe(universe_id={unit_id}, cells=[{", ".join(cell_vars)}])')

        return lines

    def _gen_lattice_universe(self, unit_id, universe):
        lines = []
        lat = universe.lattice

        lines.append(f'# Lattice {unit_id}: {lat.dimension[0]}x{lat.dimension[1]} array')
        lines.append(f'lattice{unit_id} = openmc.RectLattice(lattice_id={unit_id})')
        lines.append(f'lattice{unit_id}.lower_left = [{lat.lower_left[0]}, {lat.lower_left[1]}]')
        lines.append(f'lattice{unit_id}.pitch = [{lat.pitch[0]:.6f}, {lat.pitch[1]:.6f}]')

        rows_str = []
        for row in reversed(lat.universes):
            row_str = '[' + ', '.join(f'universe{u}' for u in row) + ']'
            rows_str.append(row_str)

        lines.append(f'lattice{unit_id}.universes = [')
        for row_str in rows_str:
            lines.append(f'    {row_str},')
        lines.append(']')

        lines.append(f'universe{unit_id} = openmc.Universe(universe_id={unit_id})')
        lines.append(f'universe{unit_id}.add_cell(openmc.Cell(fill=lattice{unit_id}))')

        return lines

    def _generate_cells(self):
        lines = []

        lines.append('# ' + '-' * 78)
        lines.append('# Root Cells')
        lines.append('# ' + '-' * 78)
        lines.append('')

        cell_vars = []
        for cell in self.parser.cells:
            cell_var = f'cell{cell.cell_id}'
            cell_vars.append(cell_var)

            if cell.name:
                lines.append(f'# {cell.name}')

            if cell.fill_universe is not None:
                if cell.translation:
                    lines.append(f'{cell_var} = openmc.Cell(cell_id={cell.cell_id}, fill=universe{cell.fill_universe})')
                    lines.append(f'{cell_var}.translation = ({cell.translation[0]}, {cell.translation[1]}, {cell.translation[2]})')
                else:
                    lines.append(f'{cell_var} = openmc.Cell(cell_id={cell.cell_id}, fill=universe{cell.fill_universe})')
            elif cell.material_id is not None and cell.material_id in self.parser.materials:
                lines.append(f'{cell_var} = openmc.Cell(cell_id={cell.cell_id}, fill=mat{cell.material_id})')
            else:
                lines.append(f'{cell_var} = openmc.Cell(cell_id={cell.cell_id})')

            region = self._convert_region(cell.region)
            if region:
                lines.append(f'{cell_var}.region = {region}')

            lines.append('')

        lines.append('root_universe = openmc.Universe(cells=[' + ', '.join(cell_vars) + '])')
        lines.append('geometry = openmc.Geometry(root_universe)')
        lines.append('')

        return lines

    def _convert_region(self, cog_region):
        if not cog_region:
            return ""

        tokens = cog_region.split()
        parts = []

        for token in tokens:
            token = token.strip()
            if not token:
                continue

            clean = token.lstrip('-+')
            if clean.isdigit():
                surf_num = int(clean)
                is_negative = token.startswith('-')

                # Check if this is a prism surface
                if surf_num in self._prism_planes:
                    plane_vars = self._prism_planes[surf_num]
                    if is_negative:
                        # Inside prism: intersection of all negative half-spaces
                        prism_region = ' & '.join(f'-{pv}' for pv in plane_vars)
                    else:
                        # Outside prism: union of positive half-spaces
                        prism_region = ' | '.join(f'+{pv}' for pv in plane_vars)
                    parts.append(f'({prism_region})')
                else:
                    if is_negative:
                        parts.append(f'-surf{surf_num}')
                    else:
                        parts.append(f'+surf{surf_num}')

        return ' & '.join(parts) if parts else ""

    def _generate_settings(self):
        lines = []

        lines.append('settings = openmc.Settings()')
        lines.append(f'settings.particles = {self.parser.npart}')
        lines.append(f'settings.batches = {self.parser.nbatch}')
        lines.append(f'settings.inactive = {self.parser.nfirst}')
        lines.append('settings.run_mode = "eigenvalue"')
        lines.append('')

        if len(self.parser.source_points) == 1:
            pt = self.parser.source_points[0]
            lines.append(f'source = openmc.IndependentSource()')
            lines.append(f'source.space = openmc.stats.Point(({pt[0]}, {pt[1]}, {pt[2]}))')
        else:
            xs = [p[0] for p in self.parser.source_points]
            ys = [p[1] for p in self.parser.source_points]
            zs = [p[2] for p in self.parser.source_points]
            lines.append(f'source = openmc.IndependentSource()')
            lines.append(f'source.space = openmc.stats.Box(({min(xs)-1}, {min(ys)-1}, {min(zs)-1}), ({max(xs)+1}, {max(ys)+1}, {max(zs)+1}))')

        lines.append('settings.source = source')
        lines.append('')

        return lines


# =============================================================================
# Main Functions
# =============================================================================

def translate_cog_to_openmc(cog_file, output_file=None):
    """Translate a COG input file to OpenMC Python format."""
    cog_path = Path(cog_file)

    if output_file is None:
        output_file = cog_path.parent / (cog_path.stem + '_openmc.py')

    parser = COGParser(cog_file)
    parser.parse()

    generator = OpenMCPythonGenerator(parser)
    python_code = generator.generate()

    with open(output_file, 'w') as f:
        f.write(python_code)

    print(f"Translated {cog_file} -> {output_file}")
    return output_file


def is_cog_file(filepath):
    """Detect if a file is a COG input file."""
    cog_keywords = {
        'basic', 'criticality', 'geometry', 'surfaces',
        'mix', 'define unit', 'sector', 'use unit',
        'npart=', 'nbatch=', 'nfirst=', 'nsource=',
        'assign-mc', 'boundary'
    }

    try:
        with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
            content = f.read(50000).lower()
            keyword_count = sum(1 for keyword in cog_keywords if keyword in content)
            return keyword_count >= 3
    except (IOError, OSError, PermissionError):
        return False


def translate_directory(input_dir, output_dir=None, recursive=False):
    """Translate all COG files in a directory."""
    input_path = Path(input_dir)
    if not input_path.is_dir():
        print(f"Error: {input_dir} is not a directory")
        return

    if output_dir:
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
    else:
        output_path = input_path

    if recursive:
        all_files = [f for f in input_path.rglob('*') if f.is_file()]
    else:
        all_files = [f for f in input_path.iterdir() if f.is_file()]

    print(f"Scanning {len(all_files)} file(s) for COG input format...")

    cog_files = [f for f in all_files if is_cog_file(f)]

    if not cog_files:
        print(f"No COG input files found in {input_dir}")
        return

    print(f"Found {len(cog_files)} COG file(s) to process")
    print("-" * 60)

    successful = 0
    failed = 0

    for cog_file in cog_files:
        try:
            out_file = output_path / (cog_file.stem + '_openmc.py')
            translate_cog_to_openmc(str(cog_file), str(out_file))
            successful += 1
        except Exception as e:
            print(f"Error processing {cog_file}: {e}")
            failed += 1

    print("-" * 60)
    print(f"Translation complete: {successful} successful, {failed} failed")


if __name__ == "__main__":
    import sys

    if len(sys.argv) < 2:
        print("COG to OpenMC Python Converter")
        print()
        print("Usage: python cog_to_openmc.py <cog_file_or_directory> [output] [-r]")
        print()
        print("Options:")
        print("  -r, --recursive    Search subdirectories recursively")
        print()
        print("Examples:")
        print("  python cog_to_openmc.py input.cog output.py")
        print("  python cog_to_openmc.py ./cog_files ./output -r")
        sys.exit(1)

    input_path = sys.argv[1]
    output_path = None
    recursive = False

    for arg in sys.argv[2:]:
        if arg in ['-r', '--recursive']:
            recursive = True
        elif output_path is None:
            output_path = arg

    if Path(input_path).is_dir():
        translate_directory(input_path, output_path, recursive=recursive)
    else:
        translate_cog_to_openmc(input_path, output_path)
