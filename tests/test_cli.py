from pathlib import Path

from typer.testing import CliRunner

from pharminform.cli import app

FIXTURE_VCF = Path(__file__).parent / "fixtures" / "synthetic_sample.vcf"

runner = CliRunner()


def test_version_command():
    result = runner.invoke(app, ["version"])
    assert result.exit_code == 0
    assert result.stdout.strip()


def test_report_command_writes_html(tmp_path):
    out = tmp_path / "report.html"
    result = runner.invoke(
        app,
        ["report", "--vcf", str(FIXTURE_VCF), "--out", str(out), "--min-depth", "10"],
    )
    assert result.exit_code == 0, result.stdout
    assert out.exists()
    assert "PharmInform" in out.read_text()


def test_report_command_missing_vcf_fails_cleanly(tmp_path):
    out = tmp_path / "report.html"
    result = runner.invoke(
        app,
        ["report", "--vcf", str(tmp_path / "does_not_exist.vcf"), "--out", str(out)],
    )
    assert result.exit_code != 0
