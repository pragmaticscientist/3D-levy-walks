"""Versioned, auditable on-disk artifacts for pipeline stages.

CSV is used deliberately: every intermediate table remains inspectable without
special tooling.  A small JSON sidecar records the schema version, row/column
counts, and a SHA-256 digest so stale or partially replaced intermediates fail
clearly instead of being consumed silently.
"""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
from typing import Any, Dict, List, Optional

import pandas as pd

from .. import SCHEMA_VERSIONS
from ..errors import LtdbLevyError
from ..logging_utils import get_logger

logger = get_logger(__name__)


class ArtifactError(LtdbLevyError):
    """An artifact is absent, corrupt, or has an unsupported schema version."""


class IncrementalFrameWriter:
    """Atomically build a large CSV artifact from bounded DataFrame chunks."""

    def __init__(
        self,
        store: "ArtifactStore",
        stage: str,
        name: str,
        schema_name: Optional[str],
        extra_metadata: Optional[Dict[str, Any]],
    ) -> None:
        self.store = store
        self.stage = stage
        self.name = name
        self.schema_name = schema_name
        self.extra_metadata = dict(extra_metadata or {})
        directory = store.stage_dir(stage)
        self.path = directory / "{}.csv".format(name)
        self.tmp_path = directory / ".{}.csv.incremental.tmp".format(name)
        self.columns: Optional[List[str]] = None
        self.rows = 0
        self.closed = False
        if self.tmp_path.exists():
            self.tmp_path.unlink()

    def append(self, frame: pd.DataFrame) -> None:
        """Append one non-empty chunk, enforcing a stable column schema."""
        if self.closed:
            raise ArtifactError("Cannot append to a closed artifact writer")
        if frame.empty:
            return
        columns = list(frame.columns)
        if self.columns is None:
            self.columns = columns
        elif columns != self.columns:
            raise ArtifactError(
                "Incremental artifact columns changed from {} to {}".format(
                    self.columns, columns
                )
            )
        frame.to_csv(
            self.tmp_path,
            mode="a",
            header=self.rows == 0,
            index=False,
        )
        self.rows += int(len(frame))

    def close(self) -> Path:
        """Commit the completed CSV and write its integrity sidecar."""
        if self.closed:
            return self.path
        if self.rows == 0 or self.columns is None:
            raise ArtifactError(
                "Incremental artifact {} received no rows".format(self.name)
            )
        os.replace(str(self.tmp_path), str(self.path))
        metadata: Dict[str, Any] = {
            "artifact": self.name,
            "format": "csv",
            "rows": self.rows,
            "columns": self.columns,
            "sha256": sha256_file(self.path),
        }
        if self.schema_name is not None:
            if self.schema_name not in SCHEMA_VERSIONS:
                raise ArtifactError(
                    "Unknown artifact schema {!r}".format(self.schema_name)
                )
            metadata["schema_name"] = self.schema_name
            metadata["schema_version"] = SCHEMA_VERSIONS[self.schema_name]
        metadata.update(self.extra_metadata)
        self.store._write_json_path(self.path.with_suffix(".meta.json"), metadata)
        self.closed = True
        logger.info("Wrote %s (%d rows, incremental)", self.path, self.rows)
        return self.path


def sha256_file(path: Path, chunk_size: int = 1024 * 1024) -> str:
    """Return the SHA-256 digest of ``path`` without loading it all at once."""
    digest = hashlib.sha256()
    with Path(path).open("rb") as fh:
        while True:
            chunk = fh.read(chunk_size)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()


class ArtifactStore:
    """Read and write stage artifacts beneath one output root."""

    def __init__(self, root: Path) -> None:
        self.root = Path(root)

    def stage_dir(self, stage: str) -> Path:
        path = self.root / stage
        path.mkdir(parents=True, exist_ok=True)
        return path

    def frame_path(self, stage: str, name: str) -> Path:
        return self.root / stage / "{}.csv".format(name)

    def json_path(self, stage: str, name: str) -> Path:
        return self.root / stage / "{}.json".format(name)

    def write_frame(
        self,
        stage: str,
        name: str,
        frame: pd.DataFrame,
        schema_name: Optional[str] = None,
        extra_metadata: Optional[Dict[str, Any]] = None,
    ) -> Path:
        """Atomically write a CSV table and its integrity sidecar."""
        directory = self.stage_dir(stage)
        path = directory / "{}.csv".format(name)
        tmp_path = directory / ".{}.csv.tmp".format(name)
        frame.to_csv(tmp_path, index=False)
        os.replace(str(tmp_path), str(path))

        metadata: Dict[str, Any] = {
            "artifact": name,
            "format": "csv",
            "rows": int(len(frame)),
            "columns": list(frame.columns),
            "sha256": sha256_file(path),
        }
        if schema_name is not None:
            if schema_name not in SCHEMA_VERSIONS:
                raise ArtifactError("Unknown artifact schema {!r}".format(schema_name))
            metadata["schema_name"] = schema_name
            metadata["schema_version"] = SCHEMA_VERSIONS[schema_name]
        if extra_metadata:
            metadata.update(extra_metadata)
        self._write_json_path(path.with_suffix(".meta.json"), metadata)
        logger.info("Wrote %s (%d rows)", path, len(frame))
        return path

    def incremental_frame_writer(
        self,
        stage: str,
        name: str,
        schema_name: Optional[str] = None,
        extra_metadata: Optional[Dict[str, Any]] = None,
    ) -> IncrementalFrameWriter:
        """Create an atomic chunk writer for an artifact too large for memory."""
        return IncrementalFrameWriter(
            self, stage, name, schema_name, extra_metadata
        )

    def read_frame(
        self,
        stage: str,
        name: str,
        schema_name: Optional[str] = None,
        verify_hash: bool = True,
        expected_config_hash: Optional[str] = None,
    ) -> pd.DataFrame:
        """Read a table after checking its sidecar and schema version."""
        path = self.frame_path(stage, name)
        meta_path = path.with_suffix(".meta.json")
        if not path.is_file():
            raise ArtifactError(
                "Required artifact is missing: {}. Run the preceding pipeline stage.".format(path)
            )
        if not meta_path.is_file():
            raise ArtifactError("Artifact sidecar is missing: {}".format(meta_path))
        with meta_path.open("r", encoding="utf-8") as fh:
            metadata = json.load(fh)
        if schema_name is not None:
            expected = SCHEMA_VERSIONS.get(schema_name)
            found = metadata.get("schema_version")
            found_name = metadata.get("schema_name")
            if found_name != schema_name or found != expected:
                raise ArtifactError(
                    "{} has schema {} version {}; expected {} version {}".format(
                        path, found_name, found, schema_name, expected
                    )
                )
        if (
            expected_config_hash is not None
            and metadata.get("config_hash") != expected_config_hash
        ):
            raise ArtifactError(
                "{} was produced with config hash {}; current hash is {}. "
                "Rerun the preceding pipeline stages.".format(
                    path, metadata.get("config_hash"), expected_config_hash
                )
            )
        if verify_hash:
            actual = sha256_file(path)
            if actual != metadata.get("sha256"):
                raise ArtifactError(
                    "Artifact hash mismatch for {}; the CSV or sidecar was modified".format(path)
                )
        frame = pd.read_csv(path)
        if int(metadata.get("rows", -1)) != len(frame):
            raise ArtifactError(
                "Artifact row-count mismatch for {} (sidecar {}, CSV {})".format(
                    path, metadata.get("rows"), len(frame)
                )
            )
        return frame

    def write_json(self, stage: str, name: str, value: Any) -> Path:
        """Atomically write a JSON artifact."""
        path = self.stage_dir(stage) / "{}.json".format(name)
        self._write_json_path(path, value)
        logger.info("Wrote %s", path)
        return path

    def read_json(self, stage: str, name: str) -> Any:
        """Read a JSON artifact written by :meth:`write_json`."""
        path = self.json_path(stage, name)
        if not path.is_file():
            raise ArtifactError(
                "Required JSON artifact is missing: {}".format(path)
            )
        with path.open("r", encoding="utf-8") as fh:
            return json.load(fh)

    def require_run_context(
        self, stage: str, expected_config_hash: str
    ) -> Dict[str, Any]:
        """Validate that a completed stage belongs to the current analysis."""
        context = self.read_json(stage, "run_context")
        found = context.get("config_hash") if isinstance(context, dict) else None
        if found != expected_config_hash:
            raise ArtifactError(
                "{} stage was produced with config hash {}; current hash is {}. "
                "Rerun that stage.".format(stage, found, expected_config_hash)
            )
        return context

    @staticmethod
    def _write_json_path(path: Path, value: Any) -> None:
        tmp_path = path.with_name(".{}.tmp".format(path.name))
        with tmp_path.open("w", encoding="utf-8") as fh:
            json.dump(value, fh, indent=2, sort_keys=True, allow_nan=False)
            fh.write("\n")
        os.replace(str(tmp_path), str(path))
