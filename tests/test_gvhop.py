"""Regression tests for GVHoP data handling without model or database downloads."""

import importlib.machinery
import importlib.util
import subprocess
import sys
import tempfile
import unittest
from contextlib import chdir
from pathlib import Path
from types import ModuleType
from unittest.mock import patch

import pandas as pd

REPOSITORY = Path(__file__).resolve().parents[1]


def load_gvhop() -> ModuleType:
    """Import the executable and its sibling module without running its CLI."""
    loader = importlib.machinery.SourceFileLoader(
        "gvhop_under_test", str(REPOSITORY / "gvhop")
    )
    spec = importlib.util.spec_from_loader(loader.name, loader)
    if spec is None:
        raise ImportError("Cannot load the GVHoP executable")
    module = importlib.util.module_from_spec(spec)
    with (
        patch.object(sys, "path", [str(REPOSITORY), *sys.path]),
        patch.dict(sys.modules, {loader.name: module}),
    ):
        loader.exec_module(module)
    return module


GVHOP = load_gvhop()


class ParserTests(unittest.TestCase):
    """Preserve real search evidence through the main-used feature builder."""

    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory(prefix="gvhop parser tests ")
        self.addCleanup(self.temporary.cleanup)
        self.directory = Path(self.temporary.name)

    def test_hmm_counts_survive_feature_and_sample_reindexing(self) -> None:
        (self.directory / "sample.domout.extra.domout").write_text(
            "protein1 - 100 100 - 90 1e-20 80 0 1 1 1e-20 1e-20 80 0 "
            "1 90 1 90 1 90 0.99 description\n"
            "protein1 - 100 200 - 90 1e-10 50 0 1 1 1e-10 1e-10 50 0 "
            "1 90 1 90 1 90 0.99 description\n"
            "protein2 - 100 100 - 90 1e-15 60 0 1 1 1e-15 1e-15 60 0 "
            "1 90 1 90 1 90 0.99 description\n# [ok]\n",
            encoding="utf-8",
        )
        (self.directory / "zero.domout").write_text("# [ok]\n", encoding="utf-8")
        (self.directory / "GVHoP_GVOGs_all.tsv").write_text(
            "missing\t100\t200\n", encoding="utf-8"
        )
        (self.directory / "GVHoP_GVEUKs_all.tsv").write_text(
            "ref2\tref1\n", encoding="utf-8"
        )
        parsed = GVHOP.hmmparser(self.directory)
        parsed_table = pd.read_csv(
            self.directory / "sample.domout.extra.domout.parsed", sep="\t"
        )
        self.assertEqual(
            parsed_table.columns.tolist(),
            [
                "protein_id",
                "accession",
                "best_hit",
                "aln_start",
                "aln_end",
                "aln_length",
                "score",
                "type",
            ],
        )
        self.assertEqual(parsed_table["protein_id"].tolist(), ["protein1", "protein2"])
        self.assertEqual(parsed_table["score"].tolist(), [80.0, 60.0])
        hmm, diamond = GVHOP.prepare_features(
            parsed,
            {"ref1": {"sample.domout.extra": 200.0}},
            ["zero", "sample.domout.extra"],
            self.directory,
        )
        self.assertEqual(hmm.index.tolist(), ["zero", "sample.domout.extra"])
        self.assertEqual(hmm.columns.tolist(), ["missing", "100", "200"])
        self.assertEqual(hmm.to_numpy().tolist(), [[0, 0, 0], [0, 2, 0]])
        self.assertEqual(diamond.index.tolist(), ["zero", "sample.domout.extra"])
        self.assertEqual(diamond.columns.tolist(), ["ref2", "ref1"])
        self.assertEqual(diamond.to_numpy().tolist(), [[0, 0], [0, 200]])

    def test_single_diamond_hit_is_not_a_header(self) -> None:
        (self.directory / "sample.txt.extra.txt").write_text(
            "q1\t00123\t99\t100\t0\t0\t1\t100\t1\t100\t1e-20\t200\n",
            encoding="utf-8",
        )
        self.assertEqual(
            GVHOP.blastparser(self.directory), {"00123": {"sample.txt.extra": 200.0}}
        )

    def test_diamond_retains_maximum_score_for_each_reference(self) -> None:
        (self.directory / "sample.txt").write_text(
            "q1\tref1\t99\t100\t0\t0\t1\t100\t1\t100\t1e-20\t200\n"
            "q2\tref1\t98\t100\t0\t0\t1\t100\t1\t100\t1e-15\t100\n"
            "q3\tref2\t97\t100\t0\t0\t1\t100\t1\t100\t1e-30\t300\n",
            encoding="utf-8",
        )
        self.assertEqual(
            GVHOP.blastparser(self.directory),
            {"ref1": {"sample": 200.0}, "ref2": {"sample": 300.0}},
        )

    def test_empty_diamond_output_is_valid_zero_hit_evidence(self) -> None:
        (self.directory / "sample.txt").write_text("", encoding="utf-8")
        self.assertEqual(GVHOP.blastparser(self.directory), {})

    def test_parsers_accept_relative_search_directories(self) -> None:
        search = self.directory / "search"
        search.mkdir()
        (search / "sample.domout").write_text(
            "protein1 - 100 100 - 90 1e-20 80 0 1 1 1e-20 1e-20 80 0 "
            "1 90 1 90 1 90 0.99 description\n# [ok]\n",
            encoding="utf-8",
        )
        (search / "sample.txt").write_text(
            "q1\tref1\t99\t100\t0\t0\t1\t100\t1\t100\t1e-20\t200\n",
            encoding="utf-8",
        )
        with chdir(self.directory):
            hmm = GVHOP.hmmparser(Path("search"))
            diamond = GVHOP.blastparser(Path("search"))
        self.assertEqual(hmm, {"sample": {"100": 1}})
        self.assertEqual(diamond, {"ref1": {"sample": 200.0}})


class PredictionCountTests(unittest.TestCase):
    """Count threshold-accepted paths without losing named sample rows."""

    def test_hierarchy_exports_preserve_counts_and_abstentions(self) -> None:
        hierarchy = [
            ["broad1", "middle1", "leafA"],
            ["broad1", "middle1", "leafB"],
            ["broad2", "middle2", "leafC"],
        ]
        nodes = {
            "leafA": [0],
            "leafB": [1],
            "leafC": [2],
            "middle1": [0, 1],
            "middle2": [2],
            "broad1": [0, 1],
            "broad2": [2],
        }
        probabilities = pd.DataFrame(
            [
                ["alpha", 0, 0.75, 0.0, 0.25],
                ["beta", 0, 0.4, 0.4, 0.2],
                ["alpha", 0, 0.4, 0.3, 0.3],
                ["beta", 2, 0.1, 0.1, 0.8],
            ]
            * 50,
            columns=["testset", "pred_label", "leafA", "leafB", "leafC"],
        )
        accepted = GVHOP.climbing_inference(probabilities, hierarchy, nodes, 0.75)
        self.assertEqual(
            accepted["pred_hierarchy"].iloc[:4].tolist(),
            [
                "leafA,middle1,broad1",
                "-,middle1,broad1",
                "-,-,-",
                "leafC,middle2,broad2",
            ],
        )
        with tempfile.TemporaryDirectory(prefix="gvhop vote tests ") as directory:
            output = Path(directory)
            GVHOP.write_prediction_counts(
                accepted, ["beta", "alpha"], hierarchy, output
            )
            leaves = pd.read_csv(output / "pred_out_level0.tsv", sep="\t", index_col=0)
            middle = pd.read_csv(output / "pred_out_level1.tsv", sep="\t", index_col=0)
            broad = pd.read_csv(output / "pred_out_level2.tsv", sep="\t", index_col=0)
        self.assertEqual(leaves.index.tolist(), ["beta", "alpha"])
        self.assertEqual(leaves.index.name, "testset")
        self.assertEqual(
            leaves.to_dict(orient="index"),
            {
                "beta": {"leafA": 0, "leafB": 0, "leafC": 50},
                "alpha": {"leafA": 50, "leafB": 0, "leafC": 0},
            },
        )
        self.assertEqual(
            middle.to_dict(orient="index"),
            {
                "beta": {"middle1": 50, "middle2": 50},
                "alpha": {"middle1": 50, "middle2": 0},
            },
        )
        self.assertEqual(
            broad.to_dict(orient="index"),
            {
                "beta": {"broad1": 50, "broad2": 50},
                "alpha": {"broad1": 50, "broad2": 0},
            },
        )


class InputTests(unittest.TestCase):
    """Validate the same filesystem inputs consumed by the CLI."""

    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory(prefix="gvhop input tests ")
        self.addCleanup(self.temporary.cleanup)
        self.directory = Path(self.temporary.name)
        self.inputs = self.directory / "input space"
        self.inputs.mkdir()
        self.database = self.directory / "database space"
        self.database.mkdir()
        self.output = self.directory / "output space"
        self.sample_file = self.directory / "samples.ls"
        self.sample_file.write_text("sample alpha\n", encoding="utf-8")
        (self.inputs / "sample alpha.faa").write_text(
            ">protein\nMKKL\n", encoding="utf-8"
        )
        (self.database / "GVHoP_GVOGs.hmm").write_text("fixture", encoding="utf-8")
        (self.database / "GVHoP_GVEUKs.dmnd").write_text("fixture", encoding="utf-8")

    def validate(self, cpu: int = 2) -> list[str]:
        """Call the production input validator with this test's paths."""
        return GVHOP.validate_inputs(
            self.inputs, "faa", self.sample_file, self.database, self.output, cpu
        )

    def test_search_commands_use_only_requested_files_and_preserve_spaces(self) -> None:
        (self.inputs / "unused.faa").write_text(">unused\nMAAA\n", encoding="utf-8")
        samples = self.validate()
        self.assertEqual(samples, ["sample alpha"])
        hmmer = GVHOP.cmd_hmmer(
            self.inputs, "faa", self.output, self.database, 2, samples
        )
        diamond = GVHOP.cmd_diamond(
            self.inputs, "faa", self.output, self.database, 2, samples
        )
        self.assertEqual(len(hmmer), 1)
        self.assertEqual(len(diamond), 1)
        self.assertEqual(hmmer[0][-1], str(self.inputs / "sample alpha.faa"))
        self.assertEqual(hmmer[0][hmmer[0].index("--cpu") + 1], "1")
        self.assertEqual(
            hmmer[0][hmmer[0].index("--domtblout") + 1],
            str(self.output / "hmmsearch/sample alpha.domout"),
        )
        self.assertEqual(
            diamond[0][diamond[0].index("-q") + 1],
            str(self.inputs / "sample alpha.faa"),
        )
        self.assertEqual(diamond[0][diamond[0].index("-p") + 1], "2")

    def test_one_cpu_uses_hmmer_without_extra_worker_threads(self) -> None:
        command = GVHOP.cmd_hmmer(
            self.inputs, "faa", self.output, self.database, 1, self.validate(cpu=1)
        )[0]
        self.assertEqual(command[command.index("--cpu") + 1], "0")

    def test_empty_sample_list_is_rejected(self) -> None:
        self.sample_file.write_text("", encoding="utf-8")
        with self.assertRaises(ValueError):
            self.validate()

    def test_duplicate_sample_ids_are_rejected(self) -> None:
        self.sample_file.write_text("sample alpha\nsample alpha\n", encoding="utf-8")
        with self.assertRaises(ValueError):
            self.validate()

    def test_sample_path_traversal_is_rejected(self) -> None:
        self.sample_file.write_text("../outside\n", encoding="utf-8")
        with self.assertRaises(ValueError):
            self.validate()

    def test_missing_requested_proteins_are_rejected(self) -> None:
        (self.inputs / "sample alpha.faa").unlink()
        with self.assertRaises(ValueError):
            self.validate()

    def test_empty_requested_proteins_are_rejected(self) -> None:
        (self.inputs / "sample alpha.faa").write_text("", encoding="utf-8")
        with self.assertRaises(ValueError):
            self.validate()

    def test_empty_database_is_rejected(self) -> None:
        (self.database / "GVHoP_GVOGs.hmm").write_text("", encoding="utf-8")
        with self.assertRaises(ValueError):
            self.validate()

    def test_nonpositive_cpu_is_rejected(self) -> None:
        with self.assertRaises(ValueError):
            self.validate(cpu=0)

    def test_nonempty_output_is_rejected_without_changing_contents(self) -> None:
        self.output.mkdir()
        sentinel = self.output / "existing.txt"
        sentinel.write_text("keep", encoding="utf-8")
        with self.assertRaises(ValueError):
            self.validate()
        self.assertEqual(sentinel.read_text(encoding="utf-8"), "keep")

    def test_existing_empty_output_is_accepted(self) -> None:
        self.output.mkdir()
        self.assertEqual(self.validate(), ["sample alpha"])

    def test_cli_rejects_invalid_cpu_before_creating_outputs(self) -> None:
        command = [
            sys.executable,
            str(REPOSITORY / "gvhop"),
            "--in_dir",
            str(self.inputs),
            "--in_file_ext",
            "faa",
            "--sample_ls",
            str(self.sample_file),
            "--db_dir",
            str(self.database),
            "--out_dir",
            str(self.output),
            "--cpu",
            "0",
        ]
        completed = subprocess.run(
            command, cwd=self.directory, capture_output=True, text=True, check=False
        )
        self.assertNotEqual(completed.returncode, 0)
        self.assertIn("cpu", completed.stderr.lower())
        self.assertFalse(self.output.exists())


class ExecutionTests(unittest.TestCase):
    """Use real lightweight child processes to check failure and argument handling."""

    def test_failed_search_stops_before_the_next_command(self) -> None:
        with tempfile.TemporaryDirectory(prefix="gvhop execution tests ") as directory:
            output = Path(directory)
            marker = output / "second command ran"
            commands = [
                [sys.executable, "-c", "print('search failed'); raise SystemExit(7)"],
                [
                    sys.executable,
                    "-c",
                    "import pathlib, sys; pathlib.Path(sys.argv[1]).touch()",
                    str(marker),
                ],
            ]
            with self.assertRaises(RuntimeError) as raised:
                GVHOP.run_jobs(commands, output, ["first", "second"])
            self.assertIsInstance(
                raised.exception.__cause__, subprocess.CalledProcessError
            )
            self.assertIn(str(output / "first.log"), str(raised.exception))
            self.assertFalse(marker.exists())
            self.assertFalse((output / "second.log").exists())
            self.assertIn(
                "search failed", (output / "first.log").read_text(encoding="utf-8")
            )

    def test_execution_preserves_literal_arguments_and_separate_logs(self) -> None:
        with tempfile.TemporaryDirectory(prefix="gvhop execution tests ") as directory:
            output = Path(directory)
            commands = [
                [
                    sys.executable,
                    "-c",
                    "import sys; print(sys.argv[1])",
                    "a b;$(literal)",
                ],
                [sys.executable, "-c", "print('second search')"],
            ]
            GVHOP.run_jobs(commands, output, ["first", "second"])
            self.assertEqual(
                (output / "first.log").read_text(encoding="utf-8"), "a b;$(literal)\n"
            )
            self.assertEqual(
                (output / "second.log").read_text(encoding="utf-8"), "second search\n"
            )


if __name__ == "__main__":
    unittest.main()
