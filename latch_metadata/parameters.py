from dataclasses import dataclass
import typing

from latch.types.metadata import SnakemakeParameter, SnakemakeFileParameter
from latch.types.file import LatchFile
from latch.types.directory import LatchDir
from pathlib import Path
from latch.types.metadata import SnakemakeMetadata, SnakemakeFileParameter
from latch.types import LatchFile, LatchDir
from latch.types.metadata import LatchAuthor

SnakemakeMetadata(
    display_name='tamirna',
    author=LatchAuthor(
        name="LatchBio",
    ),
    parameters={
        "data": SnakemakeFileParameter(
            display_name="FastQ Directory",
            type=LatchDir,
            path=Path("data")
        ),
        "sampleSheet": SnakemakeFileParameter(
            display_name="Sample Sheet Excel",
            type=LatchFile,
            path=Path("data/E-MTAB-6885_SampleContrastSheet.xlsx"),
            download=True
        ),
        "repoPath": SnakemakeFileParameter(
            display_name="Genome",
            type=LatchDir,
            path=Path("repository/data")
        )
    },
)
