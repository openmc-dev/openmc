import os
import glob

import tempfile

import openmc


def test_statepoint_latest():
    """Test that negative batch numbers keep running statepoints for last N batches."""
    # Create a temporary directory for this test
    with tempfile.TemporaryDirectory() as tmpdir:
        # Save current directory and change to temp
        original_dir = os.getcwd()
        os.chdir(tmpdir)
        
        try:
            # Copy all supporting files from the test directory into temp
            import shutil
            src_dir = os.path.dirname(__file__)
            for fname in os.listdir(src_dir):
                if fname == '__pycache__':
                    continue
                src_path = os.path.join(src_dir, fname)
                # copy files (XML, HDF5, etc.) into the temp working dir
                if os.path.isfile(src_path):
                    shutil.copy(src_path, fname)
            
            # Load the model from the single `model.xml` file
            model = openmc.Model.from_model_xml('model.xml')
            
            # Run the model
            model.run(output=False)
            
            # Check that exactly 2 running statepoint files exist
            running = sorted(glob.glob('statepoint.running.*.h5'))
            assert len(running) == 2, f'Expected 2 running statepoint files, found {len(running)}: {running}'
            
            # Verify the batch numbers are 4 and 5 (last 2 of 5 batches)
            batch_nums = sorted([int(f.split('.')[2]) for f in running])
            assert batch_nums == [4, 5], f'Expected batches [4, 5], found {batch_nums}'
        
        finally:
            os.chdir(original_dir)
