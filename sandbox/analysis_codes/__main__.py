"""Run a study directory using the shared analysis engine."""
import argparse
from pathlib import Path
from .context import load_study
from .workflow import run


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('study',type=Path,help='Directory containing model.py and config.py')
    modes=parser.add_mutually_exclusive_group(required=True)
    for mode in ('run','resume','dry-run','reanalyse'):
        modes.add_argument('--'+mode,action='store_true')
    parser.add_argument('--output',type=Path,help='Alternative output (otherwise config.RESULTS_DIRECTORY)')
    args=parser.parse_args()
    study_path = args.study
    if not study_path.exists():
        # Allow the concise study names used in the documentation while
        # keeping every editable study definition under sandbox/studies/.
        candidate = Path(__file__).resolve().parents[1] / 'studies' / study_path.name
        if candidate.is_dir():
            study_path = candidate
    study=load_study(study_path,args.output)
    if args.output is None:study.output=study.directory/study.config.RESULTS_DIRECTORY
    run(study,dry_run=args.dry_run,reanalyse=args.reanalyse)


if __name__=='__main__':main()
