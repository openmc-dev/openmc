"""Direct-script equivalent of ``python -m sandbox.analysis_codes STUDY --run``."""
if __name__ == '__main__':
    from pathlib import Path
    import sys
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
    from sandbox.analysis_codes.__main__ import main
    main()
