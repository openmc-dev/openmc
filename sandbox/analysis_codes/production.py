"""Evaluation-verified tritium production with explicit constituent scope.

Be9 MT105 is a redundant parent: score MT700+701 instead. Never score MT205
under partial-XS perturbation. This adapter fails closed on library changes.
"""
from pathlib import Path
from tempfile import TemporaryDirectory
from io import StringIO
import hashlib
import json
import os
import numpy as np
import pandas as pd
import openmc
import openmc.lib as lib
from openmc.data.endf import Evaluation,get_head_record,get_tab1_record
from sandbox.analysis_codes import Response,ModelSpec
from sandbox.analysis_codes.runner import Runner
from sandbox.analysis_codes.storage import save_json

from .context import settings
VERIFIED_MAPPING={'Li6':[105],'Li7':list(range(52,83)),'Be9':[700,701],'F19':[105]}


def mapping():
    requested=settings().PRODUCTION_CHANNELS
    if any(n not in VERIFIED_MAPPING or mts!=VERIFIED_MAPPING[n] for n,mts in requested.items()):
        raise ValueError('Unverified production channels: provide a reviewed model.responses adapter')
    return requested


def response_mapping(present):
    """Select configured QoI channels, independently of transport nuclides.

    PRODUCTION_CHANNELS is the response whitelist. Nuclides not listed there
    remain in transport and discovery without needing a production-yield map.
    Verification applies only to the requested scored channels.
    """
    selected={name:mts for name,mts in mapping().items() if name in present}
    if not selected:
        raise ValueError('No configured production nuclides are present in the model')
    return selected



def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def data_files():
    index=Path(os.environ.get('OPENMC_CROSS_SECTIONS',openmc.config.get('cross_sections',''))).resolve()
    library=openmc.data.DataLibrary.from_xml(index);result={}
    for name in mapping():
        entry=library.get_by_material(name,data_type='neutron')
        if entry is None:raise ValueError(f'Missing {name}')
        path=Path(entry['path']);path=path if path.is_absolute() else index.parent/path
        result[name]=dict(path=str(path.resolve()),sha256=sha(path))
    return result


def fetch_references(proof):
    """Fetch only the bounded archive prefix, checking ZIP CRC and pinned SHA."""
    import struct,zlib,urllib.request
    pinned={'n-003_Li_007.endf':'fbc9df82be7509fdf3de3c77510ae6671a2f045221e9a457edc9c1a068f7d208','n-003_Li_006.endf':'b2bbcfcb1ba10ec16120813d304aa6d622b61b0602a361df49cad759dad65641','n-004_Be_009.endf':'2878ae42a4f47347396418ed99fa892e2f173155945c0891e02bf477844a6ed6',
            'n-009_F_019.endf':'9b5f0b72a66aa65bf65a4f26fd2d0d345f5c30651ebbe798c0f2cba185576a8a'}
    if all((proof/name).exists() and sha(proof/name)==digest for name,digest in pinned.items()):return
    url='https://www.nndc.bnl.gov/endf-releases/releases/B-VIII.0/neutrons/ENDF-B-VIII.0_neutrons.zip'
    with urllib.request.urlopen(url,timeout=60) as stream:data=stream.read(20000000)
    offset=0
    while data[offset:offset+4]==b'PK\x03\x04':
        header=struct.unpack_from('<4s5H3I2H',data,offset)
        crc,compressed,uncompressed=header[6:9];n,e=header[-2:]
        name=Path(data[offset+30:offset+30+n].decode()).name;start=offset+30+n+e
        if start+compressed>len(data):break
        if name in pinned:
            if header[2]!=0 or header[3]!=8:raise ValueError('Unexpected ZIP representation')
            raw=zlib.decompress(data[start:start+compressed],-15)
            if len(raw)!=uncompressed or zlib.crc32(raw)!=crc or hashlib.sha256(raw).hexdigest()!=pinned[name]:
                raise ValueError('Reference ENDF integrity failure')
            (proof/name).write_bytes(raw)
        offset=start+compressed
        if all((proof/name).exists() for name in pinned):break
    if not all((proof/name).exists() and sha(proof/name)==digest for name,digest in pinned.items()):raise ValueError('Reference archive prefix did not contain verified Be/F evaluations')


def prepare(output):
    """Verify source yield metadata, native Q, and full-grid MT205 reconstruction."""
    from .context import build_model
    MAPPING=mapping()
    openmc.reset_auto_ids();model=build_model()
    present={n for material in model.geometry.get_all_materials().values() for n in material.get_nuclides()}
    MAPPING=response_mapping(present)
    output=Path(output);proof=output/'provenance';proof.mkdir(parents=True,exist_ok=True)
    existing=output/'provenance/production_mapping.json'
    if existing.exists():
        bound=json.loads(existing.read_text())['files'];current=data_files()
        if any(current[n]['sha256']!=bound[n]['sha256'] for n in MAPPING):
            raise ValueError('Installed data changed: the verified manifest is locked; create a separately reviewed adapter/output')
    fetch_references(proof)
    sources={name:proof/f'n-003_Li_{mass:03}.endf' for name,mass in [('Li6',6),('Li7',7)]}
    sources.update(Be9=proof/'n-004_Be_009.endf',F19=proof/'n-009_F_019.endf')
    sources={n:p for n,p in sources.items() if n in MAPPING}
    if not all(p.exists() for p in sources.values()):
        raise ValueError('Reference ENDF files missing; obtain the recorded NNDC VIII.0 references first')
    rows=[];evaluations={name:Evaluation(path) for name,path in sources.items()}
    # Explicit MF6 products independently establish the two Be triton yields
    # and confirm that F MF6 channels contain no additional triton products.
    from openmc.data import Reaction
    for name in ('Be9','F19'):
        if name not in evaluations:continue
        for mf,mt in evaluations[name].section:
            if mf!=6:continue
            reaction=Reaction.from_endf(evaluations[name],mt)
            tritons=[p for p in reaction.products if p.particle=='H3']
            if tritons:
                if mt not in MAPPING[name]:raise ValueError(f'Unmapped explicit triton: {name} MT{mt}')
                for p in tritons:
                    if not np.all(np.asarray(p.yield_.y)==1):raise ValueError('Non-unit triton yield needs an extended adapter')
    t=openmc.Tally(name='mapping_flux');t.scores=['flux'];model.tallies.append(t)
    with TemporaryDirectory() as directory:
        with Runner(ModelSpec('FLiBe_mapping',model,{'flux':Response(t.id)},tuple(MAPPING)),directory):
            checks=[]
            for name,evaluation in evaluations.items():
                index=lib.nuclides[name]._index;grid=lib.nuclide_energy_grid(index)
                def full(mt):
                    xs=lib.get_nuclide_xs(index,mt);values=np.zeros_like(grid);values[-len(xs):]=xs;return values
                for mf,mt in evaluation.section:
                    if mf!=3:continue
                    stream=StringIO(evaluation.section[mf,mt]);get_head_record(stream);p,_=get_tab1_record(stream)
                    if mt in MAPPING[name]:
                        info=lib.get_reaction_info(name,mt)
                        if info['redundant']:raise ValueError('Refusing a redundant production score')
                        if abs(info['q_value']-p[1])>.01:raise ValueError('ENDF/native Q mismatch')
                        if name=='Li7' and p[3]!=33:raise ValueError('Li7 breakup mapping mismatch')
                    if 51<=mt<=91 and p[3] in (33,36) and mt not in MAPPING[name]:raise ValueError('Unmapped triton breakup')
                    rows.append(dict(nuclide=name,mt=mt,QM_eV=p[0],QI_eV=p[1],LR=p[3],included=mt in MAPPING[name]))
                native=full(205);reconstructed=sum(full(mt) for mt in MAPPING[name]);mask=native>0
                np.testing.assert_allclose(native,reconstructed,rtol=1e-6,atol=1e-12)
                checks.append(dict(nuclide=name,mts=MAPPING[name],yield_per_reaction=1,
                    max_absolute_error_barn=float(max(abs(native-reconstructed))),
                    max_relative_error=float(max(abs(native[mask]-reconstructed[mask])/native[mask],default=0))))
    manifest=dict(files=data_files(),sources={n:dict(path=str(p),sha256=sha(p)) for n,p in sources.items()},
        channels=MAPPING,checks=checks,
        scope='Configured production QoI: sum of the listed channels for '+', '.join(MAPPING),
        evidence='VIII.0 MF3 MT105/LR33 and Be MF6 H3 yields; no F MF6 tritons; full-grid native MT205 comparison')
    save_json(proof/'production_mapping.json',manifest)
    pd.DataFrame(rows).to_csv(proof/'endf_channel_metadata.csv',index=False)
    return manifest


def responses(model):
    from .context import output_directory
    manifest=json.loads((output_directory()/'provenance/production_mapping.json').read_text())
    MAPPING=mapping()
    present={n for m in model.geometry.get_all_materials().values() for n in m.get_nuclides()}
    MAPPING=response_mapping(present)
    if manifest['channels']!=MAPPING:
        raise ValueError('Production response scope differs from its verified manifest')
    current=data_files()
    if any(current[n]['sha256']!=manifest['files'][n]['sha256'] for n in MAPPING):
        raise ValueError('Nuclear data differ from the verified FLiBe production mapping')
    names=list(MAPPING);scores=sorted({mt for mts in MAPPING.values() for mt in mts})
    tally=openmc.Tally(name='all_constituent_tritium');tally.estimator='tracklength'
    tally.nuclides=names;tally.scores=list(map(str,scores));model.tallies.append(tally)
    weights={}
    for i,name in enumerate(names):
        w=np.zeros((len(names),len(scores)))
        for mt in MAPPING[name]:w[i,scores.index(mt)]=1.
        weights[name]=w
    total=sum(weights.values());result={}
    for name,w in [('total',total),*weights.items()]:
        result[name]=Response(tally.id,units='tritons/source neutron',score_weights=w.tolist(),
            label=settings().RESPONSE_LABEL if name=='total' else name+' tritium production')
    return result


def verify(output):
    """Controlled native mutations for every constituent, plus exact replay."""
    from .context import build_model,register
    from .context import settings
    VERIFICATION_HISTORY_CAP=settings().VERIFICATION_HISTORY_CAP
    from sandbox.analysis_codes.budget import BudgetedRunner
    from .fingerprint import implementation_sha
    from sandbox.analysis_codes.registry import Parameter
    output=Path(output);output.mkdir(parents=True,exist_ok=True)
    MAPPING=mapping()
    openmc.reset_auto_ids();model=build_model(settings().VERIFICATION_PARTICLES,settings().VERIFICATION_BATCHES);rsp=responses(model)
    MAPPING={n:mts for n,mts in MAPPING.items() if n in rsp}
    cached=openmc.Tally(name='baseline_only_MT205_check');cached.estimator='tracklength'
    cached.nuclides=list(MAPPING);cached.scores=['205'];model.tallies.append(cached)
    for j,name in enumerate(MAPPING):
        w=np.zeros((len(MAPPING),1));w[j]=1
        rsp['cached_'+name]=Response(cached.id,score_weights=w.tolist())
    params=[Parameter(name,mt,'xs',limit=settings().PERTURBATION_LIMIT) for name,mt in [('Li6',105),('Li7',56),('Be9',700),('Be9',701),('F19',105),('Li7',2)] if name in MAPPING]
    rows=[];baseline=None
    with BudgetedRunner(register(model,rsp),output/'checkpoints',VERIFICATION_HISTORY_CAP) as runner:
        runner.bind(params)
        cases=[('baseline',np.zeros(len(params)))]+[(p.nuclide+'_MT'+str(p.mt),np.eye(len(params))[j]*settings().PERTURBATION_LIMIT) for j,p in enumerate(params)]+[('restored',np.zeros(len(params)))]
        for label,values in cases:
            record=runner.evaluate(values,stage='FLiBe_verify_'+label,seed=settings().SEED)
            batch={name:np.asarray(record['batch_responses'][name])[:,0] for name in rsp}
            np.testing.assert_allclose(batch['total'],sum(batch[name] for name in MAPPING),rtol=1e-13,atol=1e-13)
            if label=='baseline':
                baseline=batch
                for name in MAPPING:np.testing.assert_allclose(batch[name],batch['cached_'+name],rtol=1e-6,atol=1e-12)
            if label=='restored':
                for name in rsp:np.testing.assert_array_equal(batch[name],baseline[name])
            runner.state.verify()
            rows.append(dict(case=label,total_mean=float(batch['total'].mean()),
                paired_change=float((batch['total']-baseline['total']).mean()),exact_restore=True,
                **{n+'_mean':float(batch[n].mean()) for n in MAPPING},
                **{n+'_change':float((batch[n]-baseline[n]).mean()) for n in MAPPING}))
        histories=runner.history_budget.used()
    pd.DataFrame(rows).to_csv(output/'controlled_tests.csv',index=False)
    result=dict(status='passed',implementation_sha256=implementation_sha(),adapter_sha256=sha(__file__),
        perturbation_limit=settings().PERTURBATION_LIMIT,exact_restoration=True,baseline_MT205_agreement=True,restored_batch_replay_exact=True,histories=histories)
    save_json(output/'verification.json',result)
    return result
