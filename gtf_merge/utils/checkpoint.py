"""Checkpoint manager for resume support."""
import json
import logging
import os
from pathlib import Path
from typing import Set, Optional
from datetime import datetime

logger = logging.getLogger(__name__)

CHECKPOINT_FILE = "checkpoint.json"


class CheckpointManager:
    def __init__(self, checkpoint_dir: str):
        self.checkpoint_dir = Path(checkpoint_dir)
        self.checkpoint_dir.mkdir(parents=True, exist_ok=True)
        self.checkpoint_path = self.checkpoint_dir / CHECKPOINT_FILE
        self._data = self._load()

    def _load(self) -> dict:
        if self.checkpoint_path.exists():
            try:
                with open(self.checkpoint_path) as f:
                    return json.load(f)
            except Exception as e:
                logger.warning(f"Could not load checkpoint: {e}")
        return {"completed_chroms": [], "started_at": str(datetime.now()), "params": {}}

    def _save(self):
        with open(self.checkpoint_path, "w") as f:
            json.dump(self._data, f, indent=2)

    def is_completed(self, chrom: str) -> bool:
        return chrom in self._data.get("completed_chroms", [])

    def mark_completed(self, chrom: str, output_gtf: str, output_tracking: str):
        if chrom not in self._data["completed_chroms"]:
            self._data["completed_chroms"].append(chrom)
        if "chrom_outputs" not in self._data:
            self._data["chrom_outputs"] = {}
        self._data["chrom_outputs"][chrom] = {
            "gtf":      output_gtf,
            "tracking": output_tracking,
            "time":     str(datetime.now()),
        }
        self._save()

    def get_completed_chroms(self) -> Set[str]:
        return set(self._data.get("completed_chroms", []))

    def get_chrom_output(self, chrom: str) -> Optional[dict]:
        return self._data.get("chrom_outputs", {}).get(chrom)

    def save_params(self, params: dict):
        self._data["params"] = params
        self._save()

    def reset(self):
        self._data = {"completed_chroms": [], "started_at": str(datetime.now()), "params": {}}
        self._save()
        logger.info("Checkpoint reset")
