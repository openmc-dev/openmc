from collections.abc import MutableMapping
from contextlib import contextmanager
import os
from pathlib import Path
import warnings

from openmc.data import DataLibrary
from openmc.data.decay import _DECAY_ENERGY, _DECAY_PHOTON_ENERGY

__all__ = ["config"]


class _Config(MutableMapping):
    def __init__(self, data=()):
        self._mapping = {"resolve_paths": True}
        self.update(data)

    def __getitem__(self, key):
        return self._mapping[key]

    def __delitem__(self, key):
        del self._mapping[key]
        if key == "cross_sections":
            del os.environ["OPENMC_CROSS_SECTIONS"]
        elif key == "mg_cross_sections":
            del os.environ["OPENMC_MG_CROSS_SECTIONS"]
        elif key == "chain_file":
            del os.environ["OPENMC_CHAIN_FILE"]
            # Reset photon source data since it relies on chain file
            _DECAY_PHOTON_ENERGY.clear()

    def __setitem__(self, key, value):
        if key == "cross_sections":
            # Force environment variable to match
            self._set_path(key, value)
            os.environ["OPENMC_CROSS_SECTIONS"] = str(value)
        elif key == "mg_cross_sections":
            self._set_path(key, value)
            os.environ["OPENMC_MG_CROSS_SECTIONS"] = str(value)
        elif key == "chain_file":
            self._set_path(key, value)
            os.environ["OPENMC_CHAIN_FILE"] = str(value)
            # Reset photon source data since it relies on chain file
            _DECAY_PHOTON_ENERGY.clear()
            _DECAY_ENERGY.clear()
        elif key == "resolve_paths":
            self._mapping[key] = value
        else:
            raise KeyError(
                f"Unrecognized config key: {key}. Acceptable keys "
                'are "cross_sections", "mg_cross_sections", '
                '"chain_file", and "resolve_paths".'
            )

    def __iter__(self):
        return iter(self._mapping)

    def __len__(self):
        return len(self._mapping)

    def __repr__(self):
        return repr(self._mapping)

    def _set_path(self, key, value):
        self._mapping[key] = p = Path(value)
        if not p.exists():
            warnings.warn(f"'{value}' does not exist.")

    @contextmanager
    def patch(self, key, value):
        """Temporarily change a value in the configuration.

        Parameters
        ----------
        key : str
            Key to change
        value : object
            New value
        """
        previous_value = self.get(key)
        self[key] = value
        yield
        if previous_value is None:
            del self[key]
        else:
            self[key] = previous_value


def _default_config():
    """Return default configuration"""
    config = _Config()

    # Set cross sections using environment variable
    if "OPENMC_CROSS_SECTIONS" in os.environ:
        config["cross_sections"] = os.environ["OPENMC_CROSS_SECTIONS"]
    if "OPENMC_MG_CROSS_SECTIONS" in os.environ:
        config["mg_cross_sections"] = os.environ["OPENMC_MG_CROSS_SECTIONS"]

    # Set depletion chain
    chain_file = os.environ.get("OPENMC_CHAIN_FILE")
    if (
        chain_file is None
        and config.get("cross_sections") is not None
        and config["cross_sections"].exists()
    ):
        # Check for depletion chain in cross_sections.xml
        data = DataLibrary.from_xml(config["cross_sections"])
        for lib in reversed(data.libraries):
            if lib["type"] == "depletion_chain":
                chain_file = lib["path"]
                break
    if chain_file is not None:
        config["chain_file"] = chain_file

    return config


config = _default_config()
