import hashlib

from tests.getrm_oncology_consensus import (
    CYP2C_TSV_PATH,
    TSV_PATH,
    load_cyp2c_consensus,
    load_oncology_consensus,
    translate_cyp2c9,
    translate_cyp2c19,
    translate_dpyd,
    translate_tpmt,
    translate_ugt1a1,
)


def test_extract_is_unmodified():
    digest = hashlib.sha256(TSV_PATH.read_bytes()).hexdigest()
    assert digest == "3a329a683efde66af1d1ea58eeb84b58508879aace4b9a697ecae78f55e0cb5f"


def test_dpyd_translation():
    assert translate_dpyd("*1/*9") == ("*1/*1", "")
    assert translate_dpyd("c.1905+1G>A / no variants detected") == ("*2A/*1", "")
    assert translate_dpyd("c.1601G>A / c.2846A>T") == ("c.2846A>T/*1", "")
    # Two HapB3 markers are one allele.
    assert translate_dpyd("c.85T>C het, c.1129-5923C>G het, c.1236G>A het, c.2194G>A het") == (
        "c.1236G>A(HapB3)/*1", "",
    )
    assert translate_dpyd("exon 4 deletion/ no variants detected")[0] is None
    assert translate_dpyd("c.85T>C, c.557A>G / no variants detected")[0] is None
    assert translate_dpyd("*1/*2")[0] is None


def test_tpmt_translation():
    assert translate_tpmt("*1/*3A") == ("*3A/*1", "")
    assert translate_tpmt("*3C/*3C") == ("*3C/*3C", "")
    assert translate_tpmt("*1/*16")[0] is None


def test_ugt1a1_translation():
    assert translate_ugt1a1("*60/*28") == ("*28/*1", "")
    assert translate_ugt1a1("(*28 + *60)/(*28 + *60)") == ("*28/*28", "")
    assert translate_ugt1a1("(*28 + *60)/(*37 + *60)") == ("*28/*37", "")
    assert translate_ugt1a1("*28/(*60)") == ("*28/*1", "")
    assert translate_ugt1a1("*1/(*7)")[0] is None
    assert translate_ugt1a1("(*37)/*60")[0] is None


def test_every_entry_is_translated_or_explained():
    entries = load_oncology_consensus()
    assert {e.gene for e in entries} == {"DPYD", "TPMT", "UGT1A1"}
    assert all(e.consensus_diplotype or e.scope_note for e in entries)


def test_cyp2c_extract_is_unmodified():
    digest = hashlib.sha256(CYP2C_TSV_PATH.read_bytes()).hexdigest()
    assert digest == "ae0769c0bd5dbfd9e3735a82efbb815b73337e7c3c9292ea1cb045c65836eebf"


def test_cyp2c_translation():
    assert translate_cyp2c9("*2/*3") == ("*2/*3", "")
    assert translate_cyp2c9("*1/*13")[0] is None
    assert translate_cyp2c9("*1/*1 (*1/*9)")[0] is None
    assert translate_cyp2c9("*1/undefined")[0] is None
    assert translate_cyp2c19("*2/*17") == ("*17/*2", "")
    assert translate_cyp2c19("*1/*1 (*1/*17)")[0] is None
    assert translate_cyp2c19("*1/*4 or *4/*17")[0] is None


def test_every_cyp2c_entry_is_translated_or_explained():
    entries = load_cyp2c_consensus()
    assert {e.gene for e in entries} == {"CYP2C9", "CYP2C19"}
    assert all(e.consensus_diplotype or e.scope_note for e in entries)
