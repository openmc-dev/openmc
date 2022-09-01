from collections.abc import MutableMapping
import os
from pathlib import Path

from openmc.data import DataLibrary


class _Config(MutableMapping):
    def __init__(self, data=()):
        self._mapping = {}
        self.update(data)

    def __getitem__(self, key):
        return self._mapping[key]

    def __delitem__(self, key):
        del self._mapping[key]

    def __setitem__(self, key, value):
        if key == 'cross_sections':
            # Force environment variable to match
            self._mapping[key] = Path(value)
            os.environ['OPENMC_CROSS_SECTIONS'] = str(value)
        elif key == 'mg_cross_sections':
            self._mapping[key] = Path(value)
            os.environ['OPENMC_MG_CROSS_SECTIONS'] = str(value)
        elif key == 'chain_file':
            self._mapping[key] = Path(value)
        else:
            raise KeyError(f'Unrecognized config key: {key}')

    def __iter__(self):
        return iter(self._mapping)

    def __len__(self):
        return len(self._mapping)

    def __repr__(self):
        return repr(self._mapping)


# Create default configuration
config = _Config()

# Set cross sections using environment variable
if "OPENMC_CROSS_SECTIONS" in os.environ:
    config['cross_sections'] = os.environ["OPENMC_CROSS_SECTIONS"]
if "OPENMC_MG_CROSS_SECTIONS" in os.environ:
    config['mg_cross_sections'] = os.environ["OPENMC_MG_CROSS_SECTIONS"]

# Check for depletion chain in cross_sections.xml
# Set depletion chain
chain_file = os.environ.get("OPENMC_DEPLETE_CHAIN")
if chain_file is None and config['cross_sections'] is not None:
    data = DataLibrary.from_xml(config['cross_sections'])
    for lib in reversed(data.libraries):
        if lib['type'] == 'depletion_chain':
            chain_file = lib['path']
            break
if chain_file is not None:
    config['chain_file'] = chain_file
