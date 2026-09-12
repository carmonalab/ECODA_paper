#!/usr/bin/env python3
"""Standalone deterministic checks for corrected multi-batch contracts.

The fixtures are deliberately small pandas frames held entirely in memory.  No
H5AD, scheduler, submitter, or production artifact is involved.  The optional
R check uses the repository's Pixi interpreter and is skipped only when that
runtime or its fingerprint dependency is unavailable.
"""
from __future__ import annotations

import shutil
import struct
import subprocess
import sys
from pathlib import Path
from typing import Callable, Iterable

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.utils.py.batch_contract import (  # noqa: E402
    BatchContractError,
    COMPOSITE_SCALARIZATION,
    DIRECT_SCALARIZATION,
    METHOD_IDS,
    RESERVED_OBS_NAME,
    batch_contract_fingerprint,
    build_batch_composite,
    canonicalize_batch_value,
    canonicalize_batch_values,
    composite_token,
    normalize_batch_keys,
    serialize_batch_metadata,
    validate_batch_metadata,
)


CORRECTED_METHOD_IDS = frozenset(
    {
        "ECODA_authors_HR",
        "ECODA_seuratres_2",
        "ECODA_authors_HR_NULL",
        "Pseudobulk",
        "GloScope",
        "PILOT",
        "MrVI",
        "QOT",
    }
)

# This vector intentionally contains separators, a literal backslash, and
# multibyte key/value text.  Lengths are UTF-8 byte lengths, not character
# counts, and no delimiter is ever interpreted as an escape sequence.
GOLDEN_KEYS = ("site|x", "技术;")
GOLDEN_VALUES = ("A|B;C,\\", "é/β")
GOLDEN_TOKEN = (
    "ecoda_batch_composite_v1|2|"
    "6:736974657c78,9:733a417c423b432c5c;"
    "7:e68a80e69caf3b,7:733ac3a92fceb2"
)
THREE_GOLDEN_TOKEN = (
    "ecoda_batch_composite_v1|3|"
    "4:73697465,3:733a41;4:74656368,3:733a78;"
    "4:6c616e65,4:733a4c31"
)

# Fixed cross-language fingerprint vector: keys are ordered as site, tech;
# method and model are policy identifiers, not values derived by the test.
EXPECTED_FINGERPRINT_PAYLOAD = (
    "ecoda_batch_contract_v1\0"
    "8:656e636f64696e67,24:65636f64615f62617463685f636f6d706f736974655f7631;"
    "4:6b657973,24:327c343a37333639373436353b343a37343635363336383b;"
    "13:7363616c6172697a6174696f6e,12:636f6d706f736974655f7631;"
    "6:6d6574686f64,5:50494c4f54;"
    "5:6d6f64656c,16:6876675f636f6d706f736974655f7631;"
)
EXPECTED_FINGERPRINT = "d0d08e80561fa7c013f1495b18b8b713fb5dec74278949d0b90312810c9e7cc3"


def expect_failure(callback: Callable[[], object], *needles: str) -> None:
    """Require a contract error whose message identifies the violated rule."""

    try:
        callback()
    except (BatchContractError, ValueError) as exc:
        message = str(exc)
        assert any(needle in message for needle in needles), (
            f"error {message!r} did not contain any of {needles!r}"
        )
    else:
        raise AssertionError(f"expected failure containing one of {needles!r}")


def _cell_frame(sample_rows: Iterable[tuple[str, dict[str, str]]]) -> pd.DataFrame:
    """Expand sample-level factors to two cells per sample deterministically."""

    cells: list[dict[str, str]] = []
    for sample, factors in sample_rows:
        for _cell in range(2):
            cells.append({"Sample": sample, **factors, "cell_type": "T"})
    return pd.DataFrame(cells)


def scalar_fixture() -> pd.DataFrame:
    frame = _cell_frame(
        [
            ("s1", {"site": "A"}),
            ("s2", {"site": "A"}),
            ("s3", {"site": "A"}),
            ("s4", {"site": "B"}),
            ("s5", {"site": "B"}),
            ("s6", {"site": "B"}),
        ]
    )
    # Exercise factor-label handling while retaining the exact category text.
    frame["site"] = pd.Categorical(frame["site"], categories=["A", "B"])
    return frame


def two_key_fixture() -> pd.DataFrame:
    return _cell_frame(
        [
            ("s1", {"site": "A", "tech": "x"}),
            ("s2", {"site": "A", "tech": "x"}),
            ("s3", {"site": "A", "tech": "y"}),
            ("s4", {"site": "A", "tech": "y"}),
            ("s5", {"site": "B", "tech": "x"}),
            ("s6", {"site": "B", "tech": "x"}),
            ("s7", {"site": "B", "tech": "y"}),
            ("s8", {"site": "B", "tech": "y"}),
        ]
    )


def three_key_fixture() -> pd.DataFrame:
    sample_rows: list[tuple[str, dict[str, str]]] = []
    sample_number = 1
    for site in ("A", "B"):
        for tech in ("x", "y"):
            for lane in ("L1", "L2"):
                for _replicate in range(2):
                    sample_rows.append(
                        (
                            f"s{sample_number}",
                            {"site": site, "tech": tech, "lane": lane},
                        )
                    )
                    sample_number += 1
    return _cell_frame(sample_rows)


def _independent_python_builder(frame: pd.DataFrame, keys: Iterable[str]) -> tuple[str, ...]:
    """Build per-cell tokens through an independent row-oriented path.

    This intentionally does not consume any contract normalizer, canonicalizer,
    or encoder.  It models the separate MRVI/HVG in-memory builder directly
    from the written byte-level contract and is compared with the
    frame-oriented production builder.
    """

    ordered_keys = tuple(keys)
    assert len(ordered_keys) >= 2
    assert all(isinstance(key, str) and key for key in ordered_keys)

    def canonical(value: object) -> str:
        if isinstance(value, (str, np.str_)):
            return f"s:{str(value)}"
        if isinstance(value, (bool, np.bool_)):
            return "b:true" if bool(value) else "b:false"
        if isinstance(value, int) and not isinstance(value, bool):
            return f"i:{int(value):d}"
        if isinstance(value, np.integer) and getattr(value.dtype, "kind", "") == "i":
            return f"i:{int(value):d}"
        if isinstance(value, float) or (
            isinstance(value, np.floating)
            and np.dtype(value.dtype).itemsize == np.dtype(np.float64).itemsize
        ):
            numeric = float(value)
            assert np.isfinite(numeric)
            return f"f64:{struct.pack('>d', numeric).hex()}"
        raise AssertionError(f"independent builder received unsupported value {value!r}")

    tokens: list[str] = []
    for _index, row in frame.iterrows():
        pairs = []
        for key in ordered_keys:
            key_bytes = key.encode("utf-8")
            value_bytes = canonical(row[key]).encode("utf-8")
            pairs.append(
                f"{len(key_bytes)}:{key_bytes.hex()},"
                f"{len(value_bytes)}:{value_bytes.hex()}"
            )
        tokens.append(
            f"ecoda_batch_composite_v1|{len(ordered_keys)}|" + ";".join(pairs)
        )
    return tuple(tokens)


def check_normalization_and_scalar_fixture() -> None:
    # Raw scalar/list/null shapes are not rewritten by this execution-boundary
    # helper: only the validated tuple is returned for corrected consumers.
    raw_scalar = "site"
    raw_list = ["site", "tech", "lane"]
    raw_null = None
    assert raw_scalar == "site"
    assert raw_list == ["site", "tech", "lane"]
    assert raw_null is None
    raw_snapshot = {
        "scalar": raw_scalar,
        "list": list(raw_list),
        "null": raw_null,
    }

    assert normalize_batch_keys(raw_scalar) == ("site",)
    assert normalize_batch_keys(raw_list) == ("site", "tech", "lane")
    assert normalize_batch_keys(["site"]) == ("site",)
    assert normalize_batch_keys(("tech", "site")) == ("tech", "site")
    assert raw_scalar == raw_snapshot["scalar"]
    assert raw_list == raw_snapshot["list"]
    assert raw_null is raw_snapshot["null"]
    expect_failure(lambda: normalize_batch_keys(raw_null), "null")
    expect_failure(lambda: normalize_batch_keys([]), "at least one")
    expect_failure(lambda: normalize_batch_keys(["site", "site"]), "duplicate")
    expect_failure(lambda: normalize_batch_keys(["site", ""]), "nonblank")
    expect_failure(lambda: normalize_batch_keys(["site", 4]), "must be a string")
    expect_failure(lambda: normalize_batch_keys({"site", "tech"}), "ordered list")
    expect_failure(
        lambda: normalize_batch_keys(["site", "cell_type"], biological_column="cell_type"),
        "biological",
    )
    expect_failure(lambda: normalize_batch_keys("Sample"), "sample column")

    assert canonicalize_batch_value(" A ") == "s: A "
    assert canonicalize_batch_value(True) == "b:true"
    assert canonicalize_batch_value(np.bool_(False)) == "b:false"
    assert canonicalize_batch_value(np.int64(-7)) == "i:-7"
    assert canonicalize_batch_value(1.5) == "f64:3ff8000000000000"
    assert canonicalize_batch_value(1, factor=True) == "s:1"
    assert canonicalize_batch_values(["A", "B"]) == ("s:A", "s:B")
    expect_failure(lambda: canonicalize_batch_value(np.float32(1.0)), "unsupported")
    expect_failure(lambda: canonicalize_batch_value(np.uint64(1)), "unsupported")
    expect_failure(lambda: canonicalize_batch_value(pd.Timestamp("2024-01-01")), "date/time")
    expect_failure(lambda: canonicalize_batch_value(["A"]), "unsupported")

    scalar = scalar_fixture()
    validation = validate_batch_metadata(
        scalar,
        "site",
        biological_column="cell_type",
    )
    assert validation.keys == ("site",)
    assert validation.key_count == 1
    assert validation.scalarization == DIRECT_SCALARIZATION
    assert validation.n_obs == 12
    assert validation.sample_ids == ("s1", "s2", "s3", "s4", "s5", "s6")
    assert validation.levels == {"site": ("s:A", "s:B")}
    assert validation.composite_levels == ()
    assert validation.row_tokens == ()
    assert validation.scalarized_values == ("s:A", "s:A", "s:A", "s:A", "s:A", "s:A", "s:B", "s:B", "s:B", "s:B", "s:B", "s:B")
    assert validation.design_rank == validation.design_columns == 2
    assert validation.composite_design_rank == validation.composite_design_columns == 2

    one_key_list = validate_batch_metadata(scalar, ["site"], biological_column="cell_type")
    assert one_key_list.keys == validation.keys
    assert one_key_list.scalarization == validation.scalarization
    assert one_key_list.scalarized_values == validation.scalarized_values
    expect_failure(
        lambda: build_batch_composite(scalar, "site"),
        "at least two",
        "requires",
    )


def check_tokens_and_fixtures() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    assert composite_token(GOLDEN_KEYS, GOLDEN_VALUES) == GOLDEN_TOKEN
    # Mapping iteration order is deliberately reversed; configured key order
    # still controls the token order.
    assert (
        composite_token(
            GOLDEN_KEYS,
            {GOLDEN_KEYS[1]: GOLDEN_VALUES[1], GOLDEN_KEYS[0]: GOLDEN_VALUES[0]},
        )
        == GOLDEN_TOKEN
    )
    reverse_token = composite_token(
        (GOLDEN_KEYS[1], GOLDEN_KEYS[0]),
        (GOLDEN_VALUES[1], GOLDEN_VALUES[0]),
    )
    assert reverse_token != GOLDEN_TOKEN
    assert reverse_token.startswith("ecoda_batch_composite_v1|2|")
    assert composite_token(("site", "tech"), ("A", "x")) == (
        "ecoda_batch_composite_v1|2|"
        "4:73697465,3:733a41;4:74656368,3:733a78"
    )

    scalar = scalar_fixture()
    two = two_key_fixture()
    three = three_key_fixture()
    two_validation = validate_batch_metadata(two, ["site", "tech"], biological_column="cell_type")
    three_validation = validate_batch_metadata(
        three,
        ["site", "tech", "lane"],
        biological_column="cell_type",
    )

    assert two_validation.n_obs == 16
    assert two_validation.n_samples == 8
    assert two_validation.keys == ("site", "tech")
    assert two_validation.scalarization == COMPOSITE_SCALARIZATION
    assert two_validation.levels == {
        "site": ("s:A", "s:B"),
        "tech": ("s:x", "s:y"),
    }
    assert two_validation.composite_level_count == 4
    assert two_validation.design_rank == two_validation.design_columns == 3
    assert two_validation.composite_design_rank == two_validation.composite_design_columns == 4
    assert two_validation.scalarized_values[0] == (
        "ecoda_batch_composite_v1|2|"
        "4:73697465,3:733a41;4:74656368,3:733a78"
    )

    assert three_validation.n_obs == 32
    assert three_validation.n_samples == 16
    assert three_validation.keys == ("site", "tech", "lane")
    assert three_validation.scalarization == COMPOSITE_SCALARIZATION
    assert three_validation.composite_level_count == 8
    assert three_validation.design_rank == three_validation.design_columns == 4
    assert three_validation.composite_design_rank == three_validation.composite_design_columns == 8
    assert three_validation.scalarized_values[0] == (
        "ecoda_batch_composite_v1|3|"
        "4:73697465,3:733a41;4:74656368,3:733a78;"
        "4:6c616e65,4:733a4c31"
    )

    # The separate row-oriented builder and the frame-oriented builder must
    # agree exactly for both two- and three-key fixtures.
    for frame, keys in ((two, ("site", "tech")), (three, ("site", "tech", "lane"))):
        frame_builder = build_batch_composite(frame, keys)
        row_builder_values = _independent_python_builder(frame, keys)
        assert frame_builder.values == row_builder_values
        assert frame_builder.validation.scalarized_values == row_builder_values
        assert frame_builder.column_name == RESERVED_OBS_NAME

    return scalar, two, three


def check_metadata_identity_and_temporary_column(
    scalar: pd.DataFrame,
    two: pd.DataFrame,
    three: pd.DataFrame,
) -> None:
    scalar_validation = validate_batch_metadata(scalar, "site", biological_column="cell_type")
    scalar_metadata = serialize_batch_metadata(
        scalar_validation,
        method_id="ECODA_authors_HR",
        model_id="ecoda_additive_random_intercepts_v1",
    )
    assert scalar_metadata["ordered_keys"] == ["site"]
    assert scalar_metadata["scalarization"] == DIRECT_SCALARIZATION
    assert scalar_metadata["composite_levels"] == []
    assert scalar_metadata["composite_values"] == list(scalar_validation.scalarized_values)
    assert scalar_metadata["reserved_obs_name"] == RESERVED_OBS_NAME

    two_validation = validate_batch_metadata(two, ["site", "tech"], biological_column="cell_type")
    metadata = serialize_batch_metadata(
        two_validation,
        method_id="PILOT",
        model_id="hvg_composite_v1",
        include_tokens=True,
    )
    assert metadata["contract_version"] == "ecoda_batch_contract_v1"
    assert metadata["token_version"] == "ecoda_batch_composite_v1"
    assert metadata["ordered_keys"] == ["site", "tech"]
    assert metadata["keys"] == metadata["ordered_keys"]
    assert metadata["scalarization"] == COMPOSITE_SCALARIZATION
    assert metadata["sample_ids"] == list(two_validation.sample_ids)
    assert metadata["sample_constancy"] is True
    assert metadata["composite_values"] == list(two_validation.scalarized_values)
    assert metadata["tokens"] == metadata["composite_values"]
    assert metadata["fingerprint_payload"] == EXPECTED_FINGERPRINT_PAYLOAD
    assert metadata["fingerprint_payload_hex"] == EXPECTED_FINGERPRINT_PAYLOAD.encode("utf-8").hex()
    assert metadata["fingerprint"] == EXPECTED_FINGERPRINT
    assert (
        batch_contract_fingerprint(
            ["site", "tech"],
            scalarization=COMPOSITE_SCALARIZATION,
            method_id="PILOT",
            model_id="hvg_composite_v1",
        )
        == EXPECTED_FINGERPRINT
    )

    three_validation = validate_batch_metadata(
        three,
        ["site", "tech", "lane"],
        biological_column="cell_type",
    )
    three_metadata = serialize_batch_metadata(
        three_validation,
        method_id="MrVI",
        model_id="mrvi_composite_v1",
    )
    assert three_metadata["ordered_keys"] == ["site", "tech", "lane"]
    assert three_metadata["composite_level_count"] == 8
    assert three_metadata["fingerprint"] == batch_contract_fingerprint(
        ["site", "tech", "lane"],
        scalarization=COMPOSITE_SCALARIZATION,
        method_id="MrVI",
        model_id="mrvi_composite_v1",
    )

    original = two.copy(deep=True)
    before_columns = list(original.columns)
    built = build_batch_composite(original, ["site", "tech"], biological_column="cell_type")
    assert RESERVED_OBS_NAME not in original.columns
    assert list(original.columns) == before_columns
    pd.testing.assert_frame_equal(original, two)
    assert built.frame is not original
    assert RESERVED_OBS_NAME in built.frame.columns
    assert list(built.frame[RESERVED_OBS_NAME]) == list(built.values)
    assert list(built.frame["site"]) == list(original["site"])
    assert list(built.frame["tech"]) == list(original["tech"])

    # The temporary field is deleted before the caller's corrected artifact
    # boundary; the source and the final frame both have no persisted temp key.
    del built.frame[RESERVED_OBS_NAME]
    assert RESERVED_OBS_NAME not in built.frame.columns
    assert RESERVED_OBS_NAME not in original.columns
    assert RESERVED_OBS_NAME not in two.columns

    already_temporary = two.copy(deep=True)
    already_temporary[RESERVED_OBS_NAME] = "preexisting"
    expect_failure(
        lambda: validate_batch_metadata(already_temporary, ["site", "tech"]),
        "reserved",
    )
    expect_failure(
        lambda: build_batch_composite(already_temporary, ["site", "tech"]),
        "reserved",
    )


def check_rejections(two: pd.DataFrame, three: pd.DataFrame) -> None:
    for key in ("site", "tech"):
        missing = two.drop(columns=[key])
        expect_failure(
            lambda missing=missing: validate_batch_metadata(missing, ["site", "tech"]),
            "missing",
        )

    for sentinel in (None, np.nan, "", " ", "NA", " nan ", "None", "<NA>", "N/A", "null", " UNKNOWN "):
        invalid = two.copy(deep=True)
        invalid.loc[invalid.index[0], "tech"] = sentinel
        expect_failure(
            lambda invalid=invalid: validate_batch_metadata(invalid, ["site", "tech"]),
            "missing",
            "blank",
            "sentinel",
        )

    # Every declared key is checked across every cell, not just the first row
    # retained for a sample-level collapse.
    for frame, keys in (
        (two, ("site", "tech")),
        (three, ("site", "tech", "lane")),
    ):
        for key in keys:
            disagreement = frame.copy(deep=True)
            disagreement.loc[disagreement.index[1], key] = "__different__"
            expect_failure(
                lambda disagreement=disagreement, keys=keys: validate_batch_metadata(
                    disagreement, keys
                ),
                "disagrees within Sample",
                "disagrees",
            )

    confounded = _cell_frame(
        [
            ("s1", {"site": "A", "tech": "x"}),
            ("s2", {"site": "A", "tech": "x"}),
            ("s3", {"site": "A", "tech": "x"}),
            ("s4", {"site": "B", "tech": "y"}),
            ("s5", {"site": "B", "tech": "y"}),
            ("s6", {"site": "B", "tech": "y"}),
        ]
    )
    expect_failure(
        lambda: validate_batch_metadata(confounded, ["site", "tech"]),
        "rank",
        "deficient",
        "disconnected",
        "non-estimable",
    )

    composite_near_unique = _cell_frame(
        [
            ("s1", {"site": "A", "tech": "x"}),
            ("s2", {"site": "A", "tech": "y"}),
            ("s3", {"site": "A", "tech": "z"}),
            ("s4", {"site": "B", "tech": "x"}),
            ("s5", {"site": "B", "tech": "y"}),
            ("s6", {"site": "B", "tech": "z"}),
        ]
    )
    expect_failure(
        lambda: validate_batch_metadata(composite_near_unique, ["site", "tech"]),
        "near-unique",
    )

    individually_near_unique = _cell_frame(
        [
            ("s1", {"site": "A", "tech": "x"}),
            ("s2", {"site": "B", "tech": "x"}),
            ("s3", {"site": "C", "tech": "y"}),
            ("s4", {"site": "D", "tech": "y"}),
        ]
    )
    expect_failure(
        lambda: validate_batch_metadata(individually_near_unique, ["site", "tech"]),
        "near-unique",
    )


def check_method_policy(two: pd.DataFrame) -> None:
    assert set(METHOD_IDS) == set(CORRECTED_METHOD_IDS) | {"preprocess"}
    assert "pilot-gm-vae" not in METHOD_IDS
    assert "PILOT-GM-VAE" not in METHOD_IDS
    assert "pilotgm" not in METHOD_IDS

    validation = validate_batch_metadata(two, ["site", "tech"])
    for method_id in sorted(CORRECTED_METHOD_IDS):
        metadata = serialize_batch_metadata(
            validation,
            method_id=method_id,
            model_id="hvg_composite_v1",
        )
        assert metadata["method_id"] == method_id

    for forbidden in ("pilot-gm-vae", "PILOT-GM-VAE", "pilotgm"):
        expect_failure(
            lambda forbidden=forbidden: serialize_batch_metadata(
                validation,
                method_id=forbidden,
                model_id="hvg_composite_v1",
            ),
            "unsupported",
            "recognized",
        )
        expect_failure(
            lambda forbidden=forbidden: batch_contract_fingerprint(
                ["site", "tech"],
                method_id=forbidden,
                model_id="hvg_composite_v1",
            ),
            "unsupported",
            "recognized",
        )


def _r_utf8_expression(value: str) -> str:
    """Encode a Python string as an R UTF-8 expression without quoted escapes."""

    codepoints = ", ".join(f"{ord(character)}L" for character in value)
    return f"intToUtf8(c({codepoints}))"


def _r_optional_unavailable(output: str) -> bool:
    lowered = output.casefold()
    if "digest" in lowered and (
        "not available" in lowered
        or "no package" in lowered
        or "required" in lowered
        or "cannot load" in lowered
    ):
        return True
    if "rscript" in lowered and ("not found" in lowered or "no such file" in lowered):
        return True
    if "environment" in lowered and ("not found" in lowered or "unavailable" in lowered):
        return True
    return any(
        marker in lowered
        for marker in ("failed to solve", "failed to fetch", "no solution found", "could not resolve host")
    )


def check_optional_r_parity() -> None:
    """Compare the same golden token/fingerprint through Pixi's R runtime."""

    r_source = _r_utf8_expression(str(ROOT / "src/utils/batch_contract.R"))
    r_keys = ", ".join(_r_utf8_expression(key) for key in ("site|x", "技术;"))
    r_values = ", ".join(_r_utf8_expression(value) for value in GOLDEN_VALUES)
    r_three_keys = ", ".join(_r_utf8_expression(key) for key in ("site", "tech", "lane"))
    r_three_values = ", ".join(_r_utf8_expression(value) for value in ("A", "x", "L1"))
    r_fingerprint_keys = ", ".join(_r_utf8_expression(key) for key in ("site", "tech"))
    r_code = f"""
source({r_source})
composite <- ecoda_batch_composite_token(
  list({r_keys}),
  list({r_values})
)
cat(composite, "\\n", sep = "")
three_composite <- ecoda_batch_composite_token(
  list({r_three_keys}),
  list({r_three_values})
)
cat(three_composite, "\\n", sep = "")
if (!requireNamespace("digest", quietly = TRUE)) {{
  cat("__NO_DIGEST__\\n", sep = "")
  quit(status = 0)
}}
cat(ecoda_batch_fingerprint_payload_hex(
  list({r_fingerprint_keys}), "PILOT", "hvg_composite_v1"
), "\\n", sep = "")
cat(ecoda_batch_fingerprint(
  list({r_fingerprint_keys}), "PILOT", "hvg_composite_v1"
), "\\n", sep = "")
"""
    pixi = shutil.which("pixi")

    try:
        result = subprocess.run(
            [pixi, "run", "-e", "default", "Rscript", "-e", r_code],
            cwd=ROOT,
            text=True,
            encoding="utf-8",
            capture_output=True,
            check=False,
            timeout=90,
        )
    except (OSError, subprocess.TimeoutExpired) as exc:
        print(f"multibatch contracts: R parity skipped (Pixi R unavailable: {exc})")
        return

    if result.returncode != 0:
        details = f"stdout={result.stdout}\nstderr={result.stderr}"
        if _r_optional_unavailable(details):
            print("multibatch contracts: R parity skipped (optional R runtime/dependency unavailable)")
            return
        raise AssertionError(f"Pixi R parity command failed:\n{details}")

    lines = result.stdout.splitlines()
    if lines == [GOLDEN_TOKEN, THREE_GOLDEN_TOKEN, "__NO_DIGEST__"]:
        print("multibatch contracts: R parity skipped (digest package unavailable)")
        return
    expected_lines = [
        GOLDEN_TOKEN,
        THREE_GOLDEN_TOKEN,
        EXPECTED_FINGERPRINT_PAYLOAD.encode("utf-8").hex(),
        EXPECTED_FINGERPRINT,
    ]
    assert lines == expected_lines, f"unexpected Pixi R contract vectors: {lines!r}"
    print("multibatch contracts: R parity OK")


def main() -> None:
    check_normalization_and_scalar_fixture()
    scalar, two, three = check_tokens_and_fixtures()
    check_metadata_identity_and_temporary_column(scalar, two, three)
    check_rejections(two, three)
    check_method_policy(two)
    check_optional_r_parity()
    print("multibatch contracts: OK")


if __name__ == "__main__":
    main()
