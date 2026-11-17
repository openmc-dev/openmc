"""Fail closed if the evaluation-specific production mapping is not verified."""
from pathlib import Path
import hashlib,json,os,subprocess
import openmc
import openmc.lib as lib


def sha(path):
    digest=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(1024*1024),b''):digest.update(block)
    return digest.hexdigest()


def verify_library(a,out):
    index=Path(os.environ.get('OPENMC_CROSS_SECTIONS',openmc.config.get('cross_sections',''))).resolve()
    library=openmc.data.DataLibrary.from_xml(index)
    pinned=json.loads((Path(__file__).parents[1]/'data/verified_neutron_files.json').read_text())
    files={}
    names=sorted({n for material in a.model.geometry.get_all_materials().values() for n in material.get_nuclides()})
    for name in names:
        entry=library.get_by_material(name,data_type='neutron')
        if entry is None:raise ValueError(f'No neutron data for {name}')
        path=Path(entry['path']);path=path if path.is_absolute() else index.parent/path
        digest=sha(path);files[name]=dict(path=str(path.resolve()),sha256=digest)
        if name in pinned['sha256'] and digest!=pinned['sha256'][name]:
            raise ValueError(f'{name} differs from the verified production evaluation: review the channel/yield mapping before updating pins')
    code=hashlib.sha256(b''.join(p.read_bytes() for p in sorted(Path(__file__).parents[1].rglob('*.py')))).hexdigest()
    report=dict(geometry_code_sha256=code,index=str(index),index_sha256=sha(index),neutron_files=files,
        compiled_library=str(lib._dll._name),compiled_library_sha256=sha(lib._dll._name),
        openmc_version=openmc.__version__,git_revision=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),
        dirty_worktree=subprocess.check_output(['git','status','--short'],text=True),verified_production_channels=pinned['channels'])
    (Path(out)/'nuclear_data_provenance.json').write_text(json.dumps(report,indent=2))
    return report
