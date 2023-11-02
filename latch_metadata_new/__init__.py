from latch.types.metadata import SnakemakeMetadata, LatchAuthor
from latch.types.directory import LatchDir

from .parameters import generated_parameters

SnakemakeMetadata(
    output_dir=LatchDir("latch:///your_output_directory"),
    display_name="Your Workflow Name",
    author=LatchAuthor(
        name="Your Name",
    ),
    # Add more parameters
    parameters=generated_parameters,
)
