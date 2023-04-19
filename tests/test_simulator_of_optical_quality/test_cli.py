from click.testing import CliRunner

from simulator_of_optical_quality.cli import cli


def _cli_output_lines(args):
    result = CliRunner().invoke(cli, args)
    assert result.exit_code == 0
    return result.output.splitlines()


def test_help():
    lines = _cli_output_lines(['--help'])
    assert len(lines) > 0
