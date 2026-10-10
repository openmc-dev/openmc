from pathlib import Path

from openmc.executor import _process_CLI_arguments


def test_process_cli_arguments_path_input_str(tmp_path):
    # _process_CLI_arguments must return a list of str: openmc.lib.init
    # calls .encode() on each argument (#4185)
    args = list(_process_CLI_arguments(path_input=str(tmp_path)))
    assert args[-1] == str(tmp_path)
    assert all(isinstance(arg, str) for arg in args)


def test_process_cli_arguments_path_input_pathlike(tmp_path):
    args = list(_process_CLI_arguments(path_input=Path(tmp_path)))
    assert args[-1] == str(tmp_path)
    assert all(isinstance(arg, str) for arg in args)
