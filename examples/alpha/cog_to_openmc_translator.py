#!/usr/bin/env python3
"""
COG to OpenMC Input Deck Translator
Translates COG Monte Carlo input files to OpenMC Python format
Author: William Zywiec
"""

import re
import numpy as np
from pathlib import Path


class COGParser:
    """Parser for COG input files"""
    
    def __init__(self, filename):
        self.filename = filename
        self.title = ""
        self.materials = {}
        self.surfaces = {}
        self.cells = {}
        self.boundary_type = "vacuum"
        self.source = {"x": 0.0, "y": 0.0, "z": 0.0}
        self.npart = 10000
        self.nbatch = 100
        self.nfirst = 20
        self.units = "cm"
        self.libraries = {}
        self.universes = {}  # COG defined units
        self.lattices = {}  # COG lattice fills
        self.current_section = None
        self.current_unit = None
        
    def parse(self):
        """Parse the COG input file"""
        with open(self.filename, 'r') as f:
            lines = f.readlines()
        
        # Get title (first line)
        if lines:
            self.title = lines[0].strip()
        
        i = 0
        while i < len(lines):
            line = lines[i].strip()
            
            # Skip empty lines and comments
            if not line or line.startswith('$'):
                i += 1
                continue
            
            # Detect sections
            if line.lower() == 'basic':
                self.current_section = 'basic'
            elif line.lower() == 'criticality':
                self.current_section = 'criticality'
            elif line.lower().startswith('mix'):
                self.current_section = 'mix'
                self._parse_mix_header(line)
            elif line.lower() == 'geometry':
                self.current_section = 'geometry'
            elif line.lower() == 'surfaces':
                self.current_section = 'surfaces'
            elif line.lower() == 'assign-mc':
                self.current_section = 'assign-mc'
            elif line.lower() in ['end', 'end.']:
                break
            elif self.current_section:
                # Check for define unit
                if self.current_section == 'geometry' and line.lower().startswith('define unit'):
                    i = self._parse_define_unit(lines, i)
                    continue
                else:
                    self._parse_section_line(line, lines, i)
            
            i += 1
    
    def _parse_define_unit(self, lines, start_idx):
        """Parse a complete define unit block"""
        line = lines[start_idx].strip()
        parts = line.split()
        
        if len(parts) >= 3:
            unit_id = int(parts[2])
            self.current_unit = unit_id
            self.universes[unit_id] = {
                'cells': [],
                'lattice': None,
                'comment': ' '.join(parts[3:]).strip('$').strip()
            }
        
        # Parse lines until next define unit, picture, or end of geometry
        idx = start_idx + 1
        lattice_lines = []
        in_lattice = False
        
        while idx < len(lines):
            line = lines[idx].strip()
            
            if not line or line.startswith('$'):
                idx += 1
                continue
            
            # Check for end of unit definition
            if (line.lower().startswith('define unit') or 
                line.lower().startswith('picture') or
                line.lower().startswith('surfaces') or
                line.lower() in ['end', 'end.']):
                self.current_unit = None
                return idx - 1
            
            # Check for lattice start
            if '{' in line:
                in_lattice = True
                lattice_lines.append(line)
            elif in_lattice:
                lattice_lines.append(line)
                if '}' in line:
                    # Parse complete lattice
                    self._parse_lattice_fill(unit_id, lattice_lines)
                    in_lattice = False
                    lattice_lines = []
            elif 'sector' in line.lower():
                # Parse sector within unit
                match = re.search(r'sector\s+(\d+)\s+(\S+)\s+(.*)', line, re.IGNORECASE)
                if match:
                    mat_id = int(match.group(1))
                    name = match.group(2)
                    surf_expr = match.group(3).strip()
                    self.universes[unit_id]['cells'].append({
                        'mat_id': mat_id,
                        'name': name,
                        'surfaces': surf_expr
                    })
            
            idx += 1
        
        self.current_unit = None
        return idx
    
    def _parse_lattice_fill(self, unit_id, lattice_lines):
        """Parse COG lattice fill specification"""
        # Join all lines and extract content between braces
        full_text = ' '.join(lattice_lines)
        match = re.search(r'\{(.*?)\}', full_text, re.DOTALL)
        if not match:
            return
        
        content = match.group(1)
        
        # Parse lattice specification
        lattice_data = {
            'x': None,
            'y': None,
            'z': None,
            'fill': []
        }
        
        # Split into tokens
        tokens = content.split()
        i = 0
        while i < len(tokens):
            token = tokens[i].lower()
            
            if token == 'x' and i + 4 <= len(tokens):
                # x start_id n_planes [xmin xmax]
                lattice_data['x'] = {
                    'start_id': int(tokens[i+1]),
                    'n': int(tokens[i+2]),
                    'range': [float(tokens[i+3].strip('[]')), float(tokens[i+4].strip('[]'))]
                }
                i += 5
            elif token == 'y' and i + 4 <= len(tokens):
                lattice_data['y'] = {
                    'start_id': int(tokens[i+1]),
                    'n': int(tokens[i+2]),
                    'range': [float(tokens[i+3].strip('[]')), float(tokens[i+4].strip('[]'))]
                }
                i += 5
            elif token == 'z' and i + 2 <= len(tokens):
                # z might have range or just bounds
                if '[' in tokens[i+2] or i + 3 < len(tokens):
                    lattice_data['z'] = {
                        'start_id': int(tokens[i+1]),
                        'range': [float(tokens[i+2].strip('[]')), float(tokens[i+3].strip('[]'))]
                    }
                    i += 4
                else:
                    lattice_data['z'] = {
                        'start_id': int(tokens[i+1]),
                        'range': [float(tokens[i+2]), float(tokens[i+2])]
                    }
                    i += 3
            elif token == 'fill':
                # Rest is fill pattern
                i += 1
                while i < len(tokens):
                    # Parse fill pattern like "361*1" or individual numbers
                    fill_token = tokens[i]
                    if '*' in fill_token:
                        parts = fill_token.split('*')
                        count = int(parts[0])
                        unit = int(parts[1])
                        lattice_data['fill'].extend([unit] * count)
                    else:
                        try:
                            lattice_data['fill'].append(int(fill_token))
                        except:
                            pass
                    i += 1
                break
            else:
                i += 1
        
        self.universes[unit_id]['lattice'] = lattice_data
    
    def _parse_mix_header(self, line):
        """Parse mix block header for library information"""
        if 'nlib=' in line.lower():
            match = re.search(r'nlib=(\S+)', line, re.IGNORECASE)
            if match:
                self.libraries['neutron'] = match.group(1)
        if 'ptlib=' in line.lower():
            match = re.search(r'ptlib=(\S+)', line, re.IGNORECASE)
            if match:
                self.libraries['photon'] = match.group(1)
        if 'sablib=' in line.lower():
            match = re.search(r'sablib=(\S+)', line, re.IGNORECASE)
            if match:
                self.libraries['sab'] = match.group(1)
    
    def _parse_section_line(self, line, all_lines, index):
        """Parse line based on current section"""
        if self.current_section == 'basic':
            self._parse_basic(line)
        elif self.current_section == 'criticality':
            self._parse_criticality(line)
        elif self.current_section == 'mix':
            self._parse_mix(line)
        elif self.current_section == 'geometry':
            self._parse_geometry(line)
        elif self.current_section == 'surfaces':
            self._parse_surfaces(line)
    
    def _parse_basic(self, line):
        """Parse basic block"""
        parts = line.lower().split()
        if 'centimeters' in parts or 'cm' in parts:
            self.units = 'cm'
    
    def _parse_criticality(self, line):
        """Parse criticality parameters"""
        if 'npart=' in line.lower():
            match = re.search(r'npart=(\d+)', line, re.IGNORECASE)
            if match:
                self.npart = int(match.group(1))
        
        if 'nbatch=' in line.lower():
            match = re.search(r'nbatch=(\d+)', line, re.IGNORECASE)
            if match:
                self.nbatch = int(match.group(1))
        
        if 'nfirst=' in line.lower():
            match = re.search(r'nfirst=(\d+)', line, re.IGNORECASE)
            if match:
                self.nfirst = int(match.group(1))
        
        if 'nsource=' in line.lower():
            # Parse source location: nsource=1 x y z (can be on same or next line)
            match = re.search(r'nsource=\d+\s+([\d.-]+)\s+([\d.-]+)\s+([\d.-]+)', line, re.IGNORECASE)
            if match:
                self.source = {
                    'x': float(match.group(1)),
                    'y': float(match.group(2)),
                    'z': float(match.group(3))
                }
        
        # Handle source coordinates on separate line (no keywords, just numbers)
        if not line.lower().startswith(('npart', 'nbatch', 'nfirst', 'nsource', 'alpha', 'sdt', 'norm')):
            # Try to parse as three floats (source coordinates)
            parts = line.split()
            if len(parts) == 3:
                try:
                    self.source = {
                        'x': float(parts[0]),
                        'y': float(parts[1]),
                        'z': float(parts[2])
                    }
                except ValueError:
                    pass  # Not source coordinates
    
    def _parse_mix(self, line):
        """Parse material definitions"""
        # Match material definition: mat=N bunches ...
        if 'mat=' in line.lower():
            mat_match = re.search(r'mat=(\d+)', line, re.IGNORECASE)
            if mat_match:
                mat_id = int(mat_match.group(1))
                
                # Extract nuclides and densities
                # Handle special cases: (h.h2o), (c), (c30p), etc.
                # Pattern: element/isotope (with optional parentheses) followed by scientific notation number
                nuclide_pattern = r'\(?([\w.]+)\)?[\s]+([\d.eE+-]+)'
                matches = re.findall(nuclide_pattern, line, re.IGNORECASE)
                
                # Filter out non-nuclide entries (keywords)
                nuclides = []
                for nuclide, density in matches:
                    # Skip if it's a keyword like 'mat', 'bunches'
                    if nuclide.lower() not in ['mat', 'bunches']:
                        nuclides.append((nuclide, density))
                
                self.materials[mat_id] = {
                    'nuclides': nuclides,
                    'comment': self._extract_comment(line)
                }
    
    def _extract_comment(self, line):
        """Extract comment from line"""
        if '$' in line:
            return line.split('$')[1].strip()
        return ""
    
    def _parse_geometry(self, line):
        """Parse geometry definitions"""
        # Check for use unit
        if line.lower().startswith('use unit'):
            match = re.search(r'use\s+unit\s+(\d+)\s+(\S+)\s+(.*)', line, re.IGNORECASE)
            if match:
                unit_id = int(match.group(1))
                instance_name = match.group(2)
                surf_expr = match.group(3).strip()
                
                # Store use unit info in cells
                if unit_id not in self.cells:
                    self.cells[unit_id] = []
                self.cells[unit_id].append({
                    'name': instance_name,
                    'surfaces': surf_expr,
                    'is_unit_instance': True,
                    'unit_id': unit_id
                })
            return
        
        if 'sector' in line.lower():
            # sector N NAME surface_list
            match = re.search(r'sector\s+(\d+)\s+(\S+)\s+(.*)', line, re.IGNORECASE)
            if match:
                mat_id = int(match.group(1))
                name = match.group(2)
                surf_expr = match.group(3).strip()
                
                if mat_id not in self.cells:
                    self.cells[mat_id] = []
                
                self.cells[mat_id].append({
                    'name': name,
                    'surfaces': surf_expr
                })
        
        elif 'boundary' in line.lower():
            # boundary type surface_id
            parts = line.split()
            if len(parts) >= 2:
                self.boundary_type = parts[1].lower()
    
    def _parse_surfaces(self, line):
        """Parse surface definitions"""
        # Remove comments
        if '$' in line:
            comment = line.split('$')[1].strip()
            line = line.split('$')[0].strip()
        else:
            comment = ""
        
        parts = line.split()
        if len(parts) < 2:
            return
        
        surf_id = int(parts[0])
        surf_type = parts[1].lower()
        
        self.surfaces[surf_id] = {
            'type': surf_type,
            'params': parts[2:],
            'comment': comment
        }


class OpenMCGenerator:
    """Generate OpenMC input deck from parsed COG data"""
    
    def __init__(self, parser):
        self.parser = parser
        
    def generate(self):
        """Generate OpenMC Python script"""
        lines = []
        
        # Header
        lines.append('"""')
        lines.append(f'{self.parser.title}')
        lines.append('Translated from COG to OpenMC')
        lines.append('"""')
        lines.append('')
        lines.append('import openmc')
        lines.append('import numpy as np')
        lines.append('')
        
        # Materials
        lines.append('# ' + '='*78)
        lines.append('# Materials')
        lines.append('# ' + '='*78)
        lines.append('')
        
        for mat_id in sorted(self.parser.materials.keys()):
            mat_data = self.parser.materials[mat_id]
            lines.extend(self._generate_material(mat_id, mat_data))
        
        # Create materials object
        lines.append('materials = openmc.Materials([{}])'.format(
            ', '.join(f'mat{i}' for i in sorted(self.parser.materials.keys()))
        ))
        lines.append('materials.export_to_xml()')
        lines.append('')
        
        # Geometry
        lines.append('# ' + '='*78)
        lines.append('# Geometry')
        lines.append('# ' + '='*78)
        lines.append('')
        
        # Surfaces
        for surf_id in sorted(self.parser.surfaces.keys()):
            surf_data = self.parser.surfaces[surf_id]
            lines.extend(self._generate_surface(surf_id, surf_data))
        
        lines.append('')
        
        # Generate universes (from define unit blocks)
        if self.parser.universes:
            lines.extend(self._generate_universes())
        
        # Cells and geometry
        lines.extend(self._generate_geometry())
        
        # Settings
        lines.append('# ' + '='*78)
        lines.append('# Settings')
        lines.append('# ' + '='*78)
        lines.append('')
        lines.extend(self._generate_settings())
        
        # Tallies
        lines.append('# ' + '='*78)
        lines.append('# Tallies')
        lines.append('# ' + '='*78)
        lines.append('')
        lines.append('tallies = openmc.Tallies()')
        lines.append('tallies.export_to_xml()')
        lines.append('')
        
        # Run info
        lines.append('# ' + '='*78)
        lines.append('# Run OpenMC')
        lines.append('# ' + '='*78)
        lines.append('')
        lines.append('openmc.run()')
        lines.append('')
        
        return '\n'.join(lines)
    
    def _generate_material(self, mat_id, mat_data):
        """Generate OpenMC material definition"""
        lines = []
        
        comment = mat_data.get('comment', '')
        if comment:
            lines.append(f'# {comment}')
        
        lines.append(f'mat{mat_id} = openmc.Material(material_id={mat_id}, name="{comment}")')
        
        # Calculate total atom density from nuclides (excluding thermal scattering markers)
        total_density = 0.0
        for nuclide, density in mat_data['nuclides']:
            openmc_nuclide = self._convert_nuclide_name(nuclide)
            if openmc_nuclide is not None:
                # Convert density string to float
                density_str = density.replace('-', 'e-').replace('+', 'e+')
                if 'e' not in density_str.lower():
                    try:
                        total_density += float(density_str)
                    except:
                        pass
                else:
                    try:
                        total_density += float(density_str)
                    except:
                        pass
        
        # Set density - COG uses atom density in units of 1e24 atoms/cm³
        if total_density > 0:
            lines.append(f'mat{mat_id}.set_density("atom/b-cm", {total_density:.6e})')
        else:
            lines.append(f'mat{mat_id}.set_density("atom/b-cm", 1.0)  # TODO: Verify density')
        
        # Add nuclides
        has_water_sab = False
        has_graphite_sab = False
        
        for nuclide, density in mat_data['nuclides']:
            # Convert COG nuclide notation to OpenMC format
            openmc_nuclide = self._convert_nuclide_name(nuclide)
            
            if openmc_nuclide is None:
                # Check what thermal scattering to add
                if 'h.h2o' in nuclide.lower() or 'h2o' in nuclide.lower():
                    has_water_sab = True
                elif 'c30p' in nuclide.lower() or (nuclide.lower() in ['(c)', 'c.graphite']):
                    has_graphite_sab = True
                continue
            
            # Convert density (COG uses scientific notation like 4.9184-4 for 4.9184e-4)
            density_str = density.replace('-', 'e-').replace('+', 'e+')
            if 'e' not in density_str.lower():
                density_str = density_str
            
            lines.append(f'mat{mat_id}.add_nuclide("{openmc_nuclide}", {density_str})')
        
        # Add thermal scattering if water is present
        if has_water_sab:
            lines.append(f'mat{mat_id}.add_s_alpha_beta("c_H_in_H2O")')
        
        # Add thermal scattering for graphite
        if has_graphite_sab:
            lines.append(f'mat{mat_id}.add_s_alpha_beta("c_Graphite")')
        
        lines.append('')
        return lines
    
    def _convert_nuclide_name(self, cog_name):
        """Convert COG nuclide name to OpenMC format"""
        # Handle special thermal scattering cases
        cog_lower = cog_name.lower()
        
        if 'h.h2o' in cog_lower or 'h2o' in cog_lower:
            return None  # Mark for S(α,β) addition
        
        if cog_lower in ['c', 'c30p', 'c.graphite']:
            # Graphite - return as carbon but mark for potential S(α,β)
            return 'C0'  # Natural carbon
        
        if cog_lower.startswith('(') or cog_lower.endswith(')'):
            # Remove parentheses from notation like (c) or (h.h2o)
            cleaned = cog_lower.strip('()')
            if '.' in cleaned or cleaned in ['c30p']:
                return None  # Thermal scattering marker
        
        # Convert element names to proper format
        # COG: u235 -> OpenMC: U235
        # COG: fe -> OpenMC: Fe (natural)
        
        # Extract element and mass number
        match = re.match(r'([a-z]+)(\d*)', cog_name, re.IGNORECASE)
        if match:
            element = match.group(1).capitalize()
            mass = match.group(2)
            
            if mass:
                return f'{element}{mass}'
            else:
                # Natural element - OpenMC typically needs mass number
                # For now, return the element for natural composition
                # Special case for carbon
                if element.lower() == 'c':
                    return 'C0'
                return element
        
        return cog_name
    
    def _generate_surface(self, surf_id, surf_data):
        """Generate OpenMC surface definition"""
        lines = []
        
        surf_type = surf_data['type'].lower()
        params = surf_data['params']
        comment = surf_data.get('comment', '')
        
        if comment:
            lines.append(f'# {comment}')
        
        # Convert COG surface types to OpenMC
        if surf_type in ['sphere', 'sph', 'so']:
            # COG: SPHERE r [x0 y0 z0] or SO r (sphere at origin)
            if len(params) >= 1:
                if len(params) == 1:
                    # Sphere at origin
                    radius = params[0]
                    lines.append(f'surf{surf_id} = openmc.Sphere(surface_id={surf_id}, r={radius})')
                elif len(params) >= 4:
                    # Sphere with center
                    x0, y0, z0, radius = params[0], params[1], params[2], params[3]
                    lines.append(f'surf{surf_id} = openmc.Sphere(surface_id={surf_id}, x0={x0}, y0={y0}, z0={z0}, r={radius})')
        
        elif surf_type in ['cylinder', 'cyl', 'c']:
            # COG: CYLINDER [axis] radius [bounds]
            if len(params) >= 1:
                # Determine axis
                axis = 'z'  # default
                radius_idx = 0
                if params[0].lower() in ['x', 'y', 'z']:
                    axis = params[0].lower()
                    radius_idx = 1
                
                if radius_idx < len(params):
                    radius = params[radius_idx]
                    
                    # Create appropriate cylinder
                    if axis == 'x':
                        lines.append(f'surf{surf_id} = openmc.XCylinder(surface_id={surf_id}, r={radius})')
                    elif axis == 'y':
                        lines.append(f'surf{surf_id} = openmc.YCylinder(surface_id={surf_id}, r={radius})')
                    else:  # z
                        lines.append(f'surf{surf_id} = openmc.ZCylinder(surface_id={surf_id}, r={radius})')
        
        elif surf_type in ['plane', 'pla', 'p']:
            # COG: PLANE axis value or PLANE A B C D (general plane)
            if len(params) >= 2:
                if params[0].lower() in ['x', 'y', 'z']:
                    # Coordinate plane
                    axis = params[0].lower()
                    value = params[1]
                    if axis == 'x':
                        lines.append(f'surf{surf_id} = openmc.XPlane(surface_id={surf_id}, x0={value})')
                    elif axis == 'y':
                        lines.append(f'surf{surf_id} = openmc.YPlane(surface_id={surf_id}, y0={value})')
                    else:  # z
                        lines.append(f'surf{surf_id} = openmc.ZPlane(surface_id={surf_id}, z0={value})')
                elif len(params) >= 4:
                    # General plane: Ax + By + Cz = D
                    a, b, c, d = params[0], params[1], params[2], params[3]
                    lines.append(f'surf{surf_id} = openmc.Plane(surface_id={surf_id}, a={a}, b={b}, c={c}, d={d})')
        
        elif surf_type in ['box']:
            # COG: BOX dx dy dz [TR x0 y0 z0]
            if len(params) >= 3:
                dx = float(params[0])
                dy = float(params[1])
                dz = float(params[2])
                
                # Check for translation (TR keyword)
                x0, y0, z0 = 0.0, 0.0, 0.0
                for i, p in enumerate(params):
                    if str(p).lower() == 'tr' and i + 3 < len(params):
                        x0 = float(params[i + 1])
                        y0 = float(params[i + 2])
                        z0 = float(params[i + 3])
                        break
                
                # COG box is centered, convert to corners for OpenMC
                xmin = x0 - dx/2
                xmax = x0 + dx/2
                ymin = y0 - dy/2
                ymax = y0 + dy/2
                zmin = z0 - dz/2
                zmax = z0 + dz/2
                
                lines.append(f'surf{surf_id} = openmc.model.RectangularParallelepiped(')
                lines.append(f'    {xmin}, {xmax}, {ymin}, {ymax}, {zmin}, {zmax}, surface_id={surf_id})')
        
        elif surf_type in ['rpp']:
            # COG: RPP x1 x2 y1 y2 z1 z2 (rectangular parallelepiped)
            if len(params) >= 6:
                xmin, xmax = params[0], params[1]
                ymin, ymax = params[2], params[3]
                zmin, zmax = params[4], params[5]
                lines.append(f'surf{surf_id} = openmc.model.RectangularParallelepiped(')
                lines.append(f'    {xmin}, {xmax}, {ymin}, {ymax}, {zmin}, {zmax}, surface_id={surf_id})')
        
        elif surf_type in ['cone', 'c/x', 'c/y', 'c/z']:
            # COG cone surfaces
            # Determine axis from type
            if 'x' in surf_type:
                axis = 'x'
            elif 'y' in surf_type:
                axis = 'y'
            else:
                axis = 'z'
            
            if len(params) >= 3:
                # Typical format: x0 y0/z0 apex t^2
                if axis == 'x' and len(params) >= 4:
                    y0, z0, apex, t2 = params[0], params[1], params[2], params[3]
                    lines.append(f'surf{surf_id} = openmc.XCone(surface_id={surf_id}, y0={y0}, z0={z0}, x0={apex}, r2={t2})')
                elif axis == 'y' and len(params) >= 4:
                    x0, z0, apex, t2 = params[0], params[1], params[2], params[3]
                    lines.append(f'surf{surf_id} = openmc.YCone(surface_id={surf_id}, x0={x0}, z0={z0}, y0={apex}, r2={t2})')
                else:  # z-axis
                    x0, y0, apex, t2 = params[0], params[1], params[2], params[3]
                    lines.append(f'surf{surf_id} = openmc.ZCone(surface_id={surf_id}, x0={x0}, y0={y0}, z0={apex}, r2={t2})')
        
        elif surf_type in ['hex', 'hexprism']:
            # Hexagonal prism
            if len(params) >= 1:
                edge_length = params[0]
                lines.append(f'surf{surf_id} = openmc.model.HexagonalPrism(')
                lines.append(f'    edge_length={edge_length}, surface_id={surf_id})')
        
        elif surf_type in ['gq', 'quadric']:
            # General quadric surface: Ax^2 + By^2 + Cz^2 + Dxy + Eyz + Fxz + Gx + Hy + Iz + J = 0
            if len(params) >= 10:
                a, b, c, d, e, f, g, h, i, j = params[0:10]
                lines.append(f'surf{surf_id} = openmc.Quadric(surface_id={surf_id},')
                lines.append(f'    a={a}, b={b}, c={c}, d={d}, e={e}, f={f}, g={g}, h={h}, i={i}, j={j})')
        
        elif surf_type in ['px']:
            # Plane normal to X axis
            if len(params) >= 1:
                x0 = params[0]
                lines.append(f'surf{surf_id} = openmc.XPlane(surface_id={surf_id}, x0={x0})')
        
        elif surf_type in ['py']:
            # Plane normal to Y axis
            if len(params) >= 1:
                y0 = params[0]
                lines.append(f'surf{surf_id} = openmc.YPlane(surface_id={surf_id}, y0={y0})')
        
        elif surf_type in ['pz']:
            # Plane normal to Z axis
            if len(params) >= 1:
                z0 = params[0]
                lines.append(f'surf{surf_id} = openmc.ZPlane(surface_id={surf_id}, z0={z0})')
        
        else:
            # For any unrecognized surface, provide a clear translation note
            lines.append(f'# COG surface type "{surf_type}" with parameters: {" ".join(str(p) for p in params)}')
            lines.append(f'# This surface type requires manual translation to OpenMC')
            # Create a placeholder that won't break the code
            lines.append(f'surf{surf_id} = openmc.Sphere(surface_id={surf_id}, r=1.0)  # PLACEHOLDER - REPLACE THIS')
        
        lines.append('')
        return lines

    
    def _generate_cells(self):
        """Generate OpenMC cell definitions"""
        lines = []
        
        cell_counter = 0
        for mat_id in sorted(self.parser.cells.keys()):
            for cell_data in self.parser.cells[mat_id]:
                name = cell_data['name']
                surf_expr = cell_data['surfaces']
                
                lines.append(f'# Cell: {name}')
                
                # Parse surface expression
                # COG uses: -1 means inside surface 1, +1 means outside
                region_parts = []
                for surf_token in surf_expr.split():
                    if surf_token.lstrip('-+').isdigit():
                        surf_num = int(surf_token.lstrip('-+'))
                        if surf_token.startswith('-'):
                            region_parts.append(f'-surf{surf_num}')
                        else:
                            region_parts.append(f'+surf{surf_num}')
                
                region_str = ' & '.join(region_parts) if region_parts else ''
                
                lines.append(f'cell{cell_counter} = openmc.Cell(cell_id={cell_counter}, fill=mat{mat_id}, name="{name}")')
                if region_str:
                    # Generate working region assignment
                    lines.append(f'cell{cell_counter}.region = {region_str}')
                lines.append('')
                
                cell_counter += 1
        
        return lines
    
    def _generate_universes(self):
        """Generate OpenMC universes from COG define unit blocks"""
        lines = []
        lines.append('# ' + '='*78)
        lines.append('# Universes (from COG define unit blocks)')
        lines.append('# ' + '='*78)
        lines.append('')
        
        for unit_id in sorted(self.parser.universes.keys()):
            unit_data = self.parser.universes[unit_id]
            
            if unit_data.get('comment'):
                lines.append(f"# Unit {unit_id}: {unit_data['comment']}")
            
            # Check if this unit has a lattice
            if unit_data.get('lattice'):
                lines.extend(self._generate_lattice_universe(unit_id, unit_data))
            else:
                lines.extend(self._generate_simple_universe(unit_id, unit_data))
            
            lines.append('')
        
        return lines
    
    def _generate_simple_universe(self, unit_id, unit_data):
        """Generate a simple universe (pin cell) from define unit"""
        lines = []
        
        # Generate cells for this universe
        cells_in_universe = []
        for i, cell_info in enumerate(unit_data['cells']):
            mat_id = cell_info['mat_id']
            name = cell_info['name']
            surf_expr = cell_info['surfaces']
            
            cell_var = f'u{unit_id}_cell{i}'
            
            # Parse surface expression
            region_parts = []
            for surf_token in surf_expr.split():
                if surf_token.lstrip('-+').isdigit():
                    surf_num = int(surf_token.lstrip('-+'))
                    if surf_token.startswith('-'):
                        region_parts.append(f'-surf{surf_num}')
                    else:
                        region_parts.append(f'+surf{surf_num}')
            
            region_str = ' & '.join(region_parts) if region_parts else ''
            
            lines.append(f'{cell_var} = openmc.Cell(fill=mat{mat_id}, name="{name}")')
            if region_str:
                lines.append(f'{cell_var}.region = {region_str}')
            
            cells_in_universe.append(cell_var)
        
        # Create universe
        cell_list = ', '.join(cells_in_universe)
        lines.append(f'universe{unit_id} = openmc.Universe(universe_id={unit_id}, cells=[{cell_list}])')
        
        return lines
    
    def _generate_lattice_universe(self, unit_id, unit_data):
        """Generate a lattice universe from COG lattice fill"""
        lines = []
        lattice_info = unit_data['lattice']
        
        if not lattice_info:
            return lines
        
        # Calculate lattice parameters
        x_data = lattice_info.get('x')
        y_data = lattice_info.get('y')
        z_data = lattice_info.get('z')
        fill_pattern = lattice_info.get('fill', [])
        
        if not (x_data and y_data):
            lines.append(f'# TODO: Incomplete lattice data for unit {unit_id}')
            return lines
        
        # Calculate dimensions
        # COG: n is number of planes, which creates n-1 divisions
        nx = x_data['n'] - 1
        ny = y_data['n'] - 1
        xmin, xmax = x_data['range']
        ymin, ymax = y_data['range']
        
        pitch_x = (xmax - xmin) / nx
        pitch_y = (ymax - ymin) / ny
        
        lines.append(f'# Lattice for unit {unit_id}: {nx}x{ny} array')
        lines.append(f'# Pitch: ({pitch_x:.6f}, {pitch_y:.6f}) cm')
        lines.append(f'# Lower left: ({xmin}, {ymin})')
        lines.append('')
        
        # Determine which universe fills the lattice
        if fill_pattern and len(set(fill_pattern)) == 1:
            # All same unit
            fill_unit = fill_pattern[0]
            lines.append(f'# Lattice filled with universe{fill_unit}')
            lines.append(f'lattice{unit_id} = openmc.RectLattice(lattice_id={unit_id})')
            lines.append(f'lattice{unit_id}.lower_left = [{xmin}, {ymin}]')
            lines.append(f'lattice{unit_id}.pitch = [{pitch_x:.6f}, {pitch_y:.6f}]')
            lines.append(f'lattice{unit_id}.universes = [[universe{fill_unit}]*{nx}]*{ny}')
        else:
            # Mixed fill - need to create array
            lines.append(f'# Creating {nx}x{ny} lattice with mixed universes')
            lines.append(f'lattice{unit_id} = openmc.RectLattice(lattice_id={unit_id})')
            lines.append(f'lattice{unit_id}.lower_left = [{xmin}, {ymin}]')
            lines.append(f'lattice{unit_id}.pitch = [{pitch_x:.6f}, {pitch_y:.6f}]')
            lines.append(f'# TODO: Set up universe array for mixed lattice')
            lines.append(f'lattice{unit_id}.universes = [[universe{fill_pattern[0] if fill_pattern else 1}]*{nx}]*{ny}')
        
        # Store lattice as universe reference
        lines.append(f'universe{unit_id} = lattice{unit_id}  # Lattice can be used as universe')
        
        return lines
    
    def _generate_geometry(self):
        """Generate final geometry with universes and lattices"""
        lines = []
        
        # Generate root-level cells
        cell_counter = 0
        root_cells = []
        
        for mat_id in sorted(self.parser.cells.keys()):
            for cell_data in self.parser.cells[mat_id]:
                name = cell_data['name']
                surf_expr = cell_data['surfaces']
                
                # Check if this is a use unit instance
                if cell_data.get('is_unit_instance'):
                    unit_id = cell_data['unit_id']
                    lines.append(f'# Cell using unit {unit_id}: {name}')
                    
                    # Parse surface expression
                    region_parts = []
                    for surf_token in surf_expr.split():
                        if surf_token.lstrip('-+').isdigit():
                            surf_num = int(surf_token.lstrip('-+'))
                            if surf_token.startswith('-'):
                                region_parts.append(f'-surf{surf_num}')
                            else:
                                region_parts.append(f'+surf{surf_num}')
                    
                    region_str = ' & '.join(region_parts) if region_parts else ''
                    
                    lines.append(f'cell{cell_counter} = openmc.Cell(cell_id={cell_counter}, fill=universe{unit_id}, name="{name}")')
                    if region_str:
                        lines.append(f'cell{cell_counter}.region = {region_str}')
                else:
                    lines.append(f'# Cell: {name}')
                    
                    # Parse surface expression
                    region_parts = []
                    for surf_token in surf_expr.split():
                        if surf_token.lstrip('-+').isdigit():
                            surf_num = int(surf_token.lstrip('-+'))
                            if surf_token.startswith('-'):
                                region_parts.append(f'-surf{surf_num}')
                            else:
                                region_parts.append(f'+surf{surf_num}')
                    
                    region_str = ' & '.join(region_parts) if region_parts else ''
                    
                    lines.append(f'cell{cell_counter} = openmc.Cell(cell_id={cell_counter}, fill=mat{mat_id}, name="{name}")')
                    if region_str:
                        lines.append(f'cell{cell_counter}.region = {region_str}')
                
                root_cells.append(f'cell{cell_counter}')
                lines.append('')
                cell_counter += 1
        
        # Create root universe and geometry
        lines.append('# Create root universe and geometry')
        cell_list = ', '.join(root_cells)
        lines.append(f'root_universe = openmc.Universe(cells=[{cell_list}])')
        lines.append('geometry = openmc.Geometry(root_universe)')
        lines.append('geometry.export_to_xml()')
        lines.append('')
        
        return lines
    
    def _generate_settings(self):
        """Generate OpenMC settings"""
        lines = []
        
        lines.append('settings = openmc.Settings()')
        lines.append(f'settings.particles = {self.parser.npart}')
        lines.append(f'settings.batches = {self.parser.nbatch}')
        lines.append(f'settings.inactive = {self.parser.nfirst}')
        lines.append('settings.run_mode = "eigenvalue"')
        lines.append('')
        
        # Source
        x, y, z = self.parser.source['x'], self.parser.source['y'], self.parser.source['z']
        lines.append('# Source definition')
        lines.append(f'source = openmc.IndependentSource()')
        lines.append(f'source.space = openmc.stats.Point(({x}, {y}, {z}))')
        lines.append('source.angle = openmc.stats.Isotropic()')
        lines.append('source.energy = openmc.stats.Watt(a=0.988e6, b=2.249e-6)')
        lines.append('settings.source = source')
        lines.append('')
        lines.append('settings.export_to_xml()')
        lines.append('')
        
        return lines


def translate_cog_to_openmc(cog_file, output_file=None):
    """Main translation function"""
    parser = COGParser(cog_file)
    parser.parse()
    
    generator = OpenMCGenerator(parser)
    openmc_code = generator.generate()
    
    if output_file is None:
        # Generate output filename
        input_path = Path(cog_file)
        output_file = input_path.stem + '_openmc.py'
    
    with open(output_file, 'w') as f:
        f.write(openmc_code)
    
    print(f"Translated {cog_file} -> {output_file}")
    return output_file


def is_cog_file(filepath):
    """
    Detect if a file is a COG input file by examining its content.
    Returns True if the file contains COG-specific keywords/patterns.
    """
    # COG-specific keywords that identify the file format
    cog_keywords = {
        'basic', 'criticality', 'geometry', 'surfaces', 
        'mix', 'define unit', 'sector', 'use unit',
        'npart=', 'nbatch=', 'nfirst=', 'nsource=',
        'assign-mc', 'boundary'
    }
    
    try:
        with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
            # Read first reasonable chunk of file (don't read huge files entirely)
            content = f.read(50000).lower()
            
            # Check for multiple COG keywords to reduce false positives
            keyword_count = sum(1 for keyword in cog_keywords if keyword in content)
            
            # If we find at least 3 COG-specific keywords, it's likely a COG file
            return keyword_count >= 3
            
    except (IOError, OSError, PermissionError):
        return False


def translate_directory(input_dir, output_dir=None, recursive=False):
    """
    Translate all COG files in a directory by detecting them via content analysis.
    
    Args:
        input_dir: Directory to search for COG files
        output_dir: Optional output directory for translated files
        recursive: If True, search subdirectories recursively
    """
    input_path = Path(input_dir)
    if not input_path.is_dir():
        print(f"Error: {input_dir} is not a directory")
        return
    
    # Create output directory if specified
    if output_dir:
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
    else:
        output_path = input_path
    
    # Find all files (not directories) in the path
    if recursive:
        all_files = [f for f in input_path.rglob('*') if f.is_file()]
    else:
        all_files = [f for f in input_path.iterdir() if f.is_file()]
    
    print(f"Scanning {len(all_files)} file(s) for COG input format...")
    
    # Detect COG files by content
    cog_files = [f for f in all_files if is_cog_file(f)]
    
    if not cog_files:
        print(f"No COG input files found in {input_dir}")
        print("Searched for files containing COG-specific keywords (basic, criticality, geometry, etc.)")
        return
    
    print(f"Found {len(cog_files)} COG file(s) to process")
    print("-" * 60)
    
    successful = 0
    failed = 0
    
    for cog_file in cog_files:
        try:
            output_file = output_path / (cog_file.stem + '_openmc.py')
            translate_cog_to_openmc(str(cog_file), str(output_file))
            successful += 1
        except Exception as e:
            print(f"Error processing {cog_file}: {e}")
            failed += 1
    
    print("-" * 60)
    print(f"Translation complete: {successful} successful, {failed} failed")


if __name__ == "__main__":
    import sys
    
    if len(sys.argv) < 2:
        print("Usage: python cog_to_openmc.py <cog_input_file_or_directory> [output_file_or_directory] [-r]")
        print()
        print("The script intelligently detects COG files by content (keywords like 'basic',")
        print("'criticality', 'geometry', etc.) regardless of file extension.")
        print()
        print("Examples:")
        print("  Single file:  python cog_to_openmc.py input.cog output.py")
        print("  Directory:    python cog_to_openmc.py ./input_dir ./output_dir")
        print("  Current dir:  python cog_to_openmc.py .")
        print("  Recursive:    python cog_to_openmc.py ./input_dir -r")
        sys.exit(1)
    
    input_path = sys.argv[1]
    output_path = None
    recursive = False
    
    # Parse arguments
    for arg in sys.argv[2:]:
        if arg in ['-r', '--recursive']:
            recursive = True
        elif output_path is None:
            output_path = arg
    
    if Path(input_path).is_dir():
        translate_directory(input_path, output_path, recursive=recursive)
    else:
        translate_cog_to_openmc(input_path, output_path)
