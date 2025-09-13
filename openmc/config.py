"""Module for handling global configuration in OpenMC.

This module exports a single object, `config`, that can be used to control
various settings, primarily paths to data files. It acts like a dictionary but
with special behaviors.

Examples
--------
>>> import openmc
>>> openmc.config['cross_sections'] = '/path/to/my/cross_sections.xml'
>>> print(openmc.config)
{'resolve_paths': True, 'cross_sections': PosixPath('/path/to/my/cross_sections.xml')}

"""
from collections.abc import MutableMapping
from contextlib import contextmanager
import os
from pathlib import Path
import warnings
from typing import Any, Dict, Iterator

from openmc.data import DataLibrary
from openmc.data.decay import _DECAY_ENERGY, _DECAY_PHOTON_ENERGY

__all__ = ["config"]


class _Config(MutableMapping):
    """A configuration dictionary for OpenMC with special handling for path-like values.

    This class enforces valid configuration keys and synchronizes path-related
    settings with their corresponding environment variables.

    Attributes
    ----------
    cross_sections : pathlib.Path
        Path to a cross_sections.xml file. Also sets/unsets the
        OPENMC_CROSS_SECTIONS environment variable.
    mg_cross_sections : pathlib.Path
        Path to a multi-group cross_sections.h5 file. Also sets/unsets
        the OPENMC_MG_CROSS_SECTIONS environment variable.
    chain_file : pathlib.Path
        Path to a depletion chain XML file. Also sets/unsets the
        OPENMC_CHAIN_FILE environment variable. Setting or deleting this
        clears internal decay data caches.
    resolve_paths : bool
        If True (default), all paths assigned are resolved to absolute
        paths. If False, paths are stored as they are provided.

    """
    _PATH_KEYS: Dict[str, str] = {
        'cross_sections': 'OPENMC_CROSS_SECTIONS',
        'mg_cross_sections': 'OPENMC_MG_CROSS_SECTIONS',
        'chain_file': 'OPENMC_CHAIN_FILE'
    }

    def __init__(self, data: dict = ()):
        self._mapping: Dict[str, Any] = {'resolve_paths': True}
        self.update(data)

    def __getitem__(self, key: str) -> Any:
        return self._mapping[key]

    def __delitem__(self, key: str):
        """Delete a configuration key.

        This also deletes the corresponding environment variable if the key is a
        path-like key, and clears decay data caches if 'chain_file' is deleted.
        'resolve_paths' cannot be deleted.

        """
        if key == 'resolve_paths':
            raise KeyError("'resolve_paths' cannot be deleted.")
        del self._mapping[key]
        if key in self._PATH_KEYS:
            env_var = self._PATH_KEYS[key]
            if env_var in os.environ:
                del os.environ[env_var]
        if key == 'chain_file':
            _DECAY_PHOTON_ENERGY.clear()
            _DECAY_ENERGY.clear()

    def __setitem__(self, key: str, value: Any):
        """Set a configuration key and its corresponding value.

        For path-like keys, this method performs several actions:
        1. Resolves the path to an absolute path if `resolve_paths` is True.
        2. Stores the `pathlib.Path` object.
        3. Sets the corresponding environment variable (e.g., OPENMC_CROSS_SECTIONS).
        4. For 'chain_file', clears internal decay data caches.
        5. Issues a `UserWarning` if the final path does not exist.

        """
        if key in self._PATH_KEYS:
            p = Path(value)
            # Use .get() for robustness, defaulting to True
            if self._mapping.get('resolve_paths', True):
                stored_path = p.resolve(strict=False)
            else:
                stored_path = p

            self._mapping[key] = stored_path
            os.environ[self._PATH_KEYS[key]] = str(stored_path)

            if key == 'chain_file':
                _DECAY_PHOTON_ENERGY.clear()
                _DECAY_ENERGY.clear()

            if not stored_path.exists():
                warnings.warn(f"Path '{stored_path}' does not exist.", UserWarning)

        elif key == 'resolve_paths':
            if not isinstance(value, bool):
                raise TypeError("'resolve_paths' must be a boolean.")
            self._mapping[key] = value
        else:
            valid_keys = list(self._PATH_KEYS.keys()) + ['resolve_paths']
            raise KeyError(
                f"Unrecognized config key: {key}. Acceptable keys are: "
                f"{', '.join(repr(k) for k in valid_keys)}."
            )

    def __iter__(self) -> Iterator[str]:
        return iter(self._mapping)

    def __len__(self) -> int:
        return len(self._mapping)

    def __repr__(self) -> str:
        return repr(self._mapping)

    def clear(self):
        """Clear all configuration keys except for 'resolve_paths'.

        This ensures that the path resolution behavior is not accidentally reset
        when clearing the configuration.

        """
        # Create a copy of keys to iterate over for safe deletion
        keys_to_delete = [k for k in self._mapping if k != 'resolve_paths']
        for key in keys_to_delete:
            del self[key]

    @contextmanager
    def patch(self, key: str, value: Any):
        """Context manager to temporarily change a configuration value.

        After the `with` block, the configuration is restored to its original
        state.

        Parameters
        ----------
        key : str
            The key of the configuration value to change.
        value
            The new temporary value.

        Examples
        --------
        >>> openmc.config['cross_sections'] = 'endf71.xml'
        >>> with openmc.config.patch('cross_sections', 'fendl32.xml'):
        ...     # Code in this block sees the new value
        ...     print(f"Inside with block: {openmc.config['cross_sections']}")
        >>> # Outside the block, the value is restored
        >>> print(f"Outside with block: {openmc.config['cross_sections']}")
        Inside with block: fendl32.xml
        Outside with block: endf71.xml

        """
        previous_value = self.get(key)
        self[key] = value
        try:
            yield
        finally:
            if previous_value is None:
                del self[key]
            else:
                self[key] = previous_value


def _default_config(**kwargs) -> _Config:
    """Create a configuration initialized from environment variables.

    This function checks for OPENMC_CROSS_SECTIONS, OPENMC_MG_CROSS_SECTIONS,
    and OPENMC_CHAIN_FILE environment variables. It also has logic to find
    a chain file within a `cross_sections.xml` file if one is not
    explicitly set.

    Returns
    -------
    _Config
        A new configuration object.

    """
    config = _Config(kwargs)
    for key,var in _Config._PATH_KEYS.items():
        if var in os.environ:
            config[key] = os.environ[var]

    chain_file = config.get("chain_file")
    xs_path = config.get("cross_sections")
    if chain_file is None and xs_path is not None and xs_path.exists():
        try:
            data = DataLibrary.from_xml(xs_path)
        except Exception:
            # Let this pass silently if cross_sections.xml can't be parsed
            # or if a dependency like lxml is not available.
            pass
        else:
            for lib in reversed(data.libraries):
                if lib['type'] == 'depletion_chain':
                    config['chain_file'] = xs_path.parent / lib['path']
                    break
    return config


# Global configuration dictionary for OpenMC settings.
config = _default_config()
