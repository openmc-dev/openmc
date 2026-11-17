"""Adapter from the shared sensitivity runner to the native tokamak model."""
from dataclasses import replace
import numpy as np
import openmc
from sandbox.tokomak_analysis.config import Config as TokamakConfig, PFCoil
from sandbox.tokomak_analysis.geometry.assembly import build_model as build_tokamak
from sandbox.analysis_codes.config import Response
from . import config as c


def build_model(particles=1000, batches=10):
    """Build native CSG with just the verified reaction-production tally."""
    settings_dict = dict(c.TOKAMAK_CONFIG)
    settings_dict["pf_coils"] = tuple(PFCoil(**coil) for coil in settings_dict.get("pf_coils", []))
    settings = replace(TokamakConfig(**settings_dict), particles=particles, batches=batches,
                       threads=c.TRANSPORT_THREADS, enable_mesh=False,
                       enable_heating=False, enable_exposure=False).validate()
    assembly = build_tokamak(settings, tallies=False)
    return assembly.model


def responses(model):
    """Return evaluation-verified, reaction-specific production responses."""
    mapping = {name: mts for name, mts in c.PRODUCTION_CHANNELS.items()
               if any(name in material.get_nuclides()
                      for material in model.geometry.get_all_materials().values())}
    names = list(mapping)
    scores = sorted({mt for mts in mapping.values() for mt in mts})
    tally = openmc.Tally(name='tokamak_verified_tritium_production')
    tally.estimator = 'tracklength'
    tally.nuclides = names
    tally.scores = list(map(str, scores))
    model.tallies.append(tally)
    weights = {}
    for name in names:
        w = np.zeros((len(names), len(scores)))
        for mt in mapping[name]:
            w[names.index(name), scores.index(mt)] = 1.0
        weights[name] = w.tolist()
    result = {'total': sum((np.asarray(w) for w in weights.values()),
                           np.zeros((len(names), len(scores))))}
    result.update({name: np.asarray(w) for name, w in weights.items()})
    return {name: Response(tally.id, units='tritons/source neutron',
                           score_weights=w.tolist(), label=c.RESPONSE_LABEL if name == 'total' else name)
            for name, w in result.items()}


def verify_response(output):
    """Run the shared ENDF-yield and loaded MT205 reconstruction audit."""
    from sandbox.analysis_codes.production import prepare, verify
    prepare(output.parent)
    verify(output)
