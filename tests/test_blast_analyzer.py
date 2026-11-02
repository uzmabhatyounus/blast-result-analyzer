
import unittest
import tempfile
import csv
from pathlib import Path
from unittest.mock import patch

from blast_analyzer import (
    read_blast_results,
    convert_values,
    rank_hits,
    display_results,
)


class TestReadBlastResults(unittest.TestCase):

    def test_read_valid_file(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            filename = Path(temp_dir) / "blast.tsv"

            filename.write_text(
                "sacc\tpident\tqcovs\tevalue\tbitscore\tstitle\n"
                "ABC123\t98.5\t95\t1e-80\t250\tProtein A\n",
                encoding="utf-8"
            )

            hits, columns = read_blast_results(filename)

            self.assertEqual(len(hits), 1)
            self.assertEqual(hits[0]["sacc"], "ABC123")
            self.assertIn("pident", columns)

    def test_empty_file(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            filename = Path(temp_dir) / "empty.tsv"
            filename.write_text("", encoding="utf-8")

            hits, columns = read_blast_results(filename)

            self.assertEqual(hits, [])
            self.assertEqual(columns, set())

    def test_missing_file(self):
        with self.assertRaises(FileNotFoundError):
            read_blast_results("nonexistent.tsv")


class TestConvertValues(unittest.TestCase):

    def setUp(self):
        self.hits = [
            {
                "sacc": "ABC123",
                "pident": "98.5",
                "qcovs": "95",
                "evalue": "1e-80",
                "bitscore": "250",
                "stitle": "Protein A",
            }
        ]

    def test_convert_numeric_values(self):
        hits = convert_values(self.hits)

        self.assertIsInstance(hits[0]["pident"], float)
        self.assertIsInstance(hits[0]["qcovs"], float)
        self.assertIsInstance(hits[0]["evalue"], float)
        self.assertIsInstance(hits[0]["bitscore"], float)

    def test_missing_numeric_value(self):
        self.hits[0]["evalue"] = ""

        with self.assertRaises(ValueError):
            convert_values(self.hits)

    def test_invalid_numeric_value(self):
        self.hits[0]["pident"] = "invalid"

        with self.assertRaises(ValueError):
            convert_values(self.hits)

    def test_missing_accession(self):
        self.hits[0]["sacc"] = ""

        with self.assertRaises(ValueError):
            convert_values(self.hits)

    def test_missing_description(self):
        self.hits[0]["stitle"] = ""

        hits = convert_values(self.hits)

        self.assertEqual(
            hits[0]["stitle"],
            "No description available"
        )


class TestRankHits(unittest.TestCase):

    def test_rank_by_evalue(self):
        hits = [
            {"evalue": 1e-20, "bitscore": 100,
             "pident": 90, "qcovs": 80},

            {"evalue": 1e-50, "bitscore": 150,
             "pident": 95, "qcovs": 90},
        ]

        ranked = rank_hits(hits)

        self.assertEqual(ranked[0]["evalue"], 1e-50)

    def test_rank_by_bitscore_when_evalue_equal(self):
        hits = [
            {"evalue": 1e-20, "bitscore": 100,
             "pident": 90, "qcovs": 80},

            {"evalue": 1e-20, "bitscore": 200,
             "pident": 95, "qcovs": 90},
        ]

        ranked = rank_hits(hits)

        self.assertEqual(ranked[0]["bitscore"], 200)

    def test_empty_hits(self):
        self.assertEqual(rank_hits([]), [])

    def test_original_list_not_modified(self):
        hits = [
            {"evalue": 1e-20, "bitscore": 100,
             "pident": 90, "qcovs": 80},

            {"evalue": 1e-50, "bitscore": 200,
             "pident": 95, "qcovs": 90},
        ]

        original_first = hits[0]

        rank_hits(hits)

        self.assertIs(hits[0], original_first)


class TestDisplayResults(unittest.TestCase):

    def test_display_results(self):
        hits = [
            {
                "sacc": "ABC123",
                "pident": 98.5,
                "qcovs": 95.0,
                "evalue": 1e-80,
                "bitscore": 250.0,
                "stitle": "Protein A",
            }
        ]

        with patch("builtins.print") as mock_print:
            display_results(hits)

        printed_text = " ".join(
            str(call) for call in mock_print.call_args_list
        )

        self.assertIn("ABC123", printed_text)
        self.assertIn("Protein A", printed_text)

    def test_display_empty_results(self):
        with patch("builtins.print") as mock_print:
            display_results([])

        mock_print.assert_any_call(
            "\nBLAST RESULT ANALYZER"
        )


if __name__ == "__main__":
    unittest.main()
