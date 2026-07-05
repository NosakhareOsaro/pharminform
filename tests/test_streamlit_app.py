from pathlib import Path

from streamlit.testing.v1 import AppTest

APP_PATH = Path(__file__).parent.parent / "streamlit_app" / "app.py"
FIXTURE_VCF = Path(__file__).parent / "fixtures" / "synthetic_sample.vcf"


def test_app_loads_without_error():
    at = AppTest.from_file(str(APP_PATH))
    at.run(timeout=30)
    assert not at.exception


def test_app_shows_upload_prompt_with_no_file():
    at = AppTest.from_file(str(APP_PATH))
    at.run(timeout=30)
    assert any("Upload a VCF file" in md.value for md in at.markdown)
