from dataclasses import dataclass
import typing

from latch.types.metadata import SnakemakeParameter, SnakemakeFileParameter
from latch.types.file import LatchFile
from latch.types.directory import LatchDir

@dataclass
class threads:
    high: int
    medium: int
    low: int


@dataclass
class Homo_sapiens:
    genomeID: str
    code: str
    genomeVersion: str
    txid: str


@dataclass
class Mus_musculus:
    genomeID: str
    code: str
    genomeVersion: str
    txid: str


@dataclass
class Sus_scrofa:
    genomeID: str
    code: str
    genomeVersion: str
    txid: str


@dataclass
class Equus_caballus:
    genomeID: str
    code: str
    genomeVersion: str
    txid: str


@dataclass
class Caenorhabditis_elegans:
    genomeID: str
    code: str
    genomeVersion: str
    txid: str


@dataclass
class Felis_catus:
    genomeID: str
    code: str
    genomeVersion: str
    txid: str


@dataclass
class Canis_lupus_familiaris:
    genomeID: str
    code: str
    genomeVersion: str
    txid: str


@dataclass
class Bos_taurus:
    genomeID: str
    code: str
    genomeVersion: str
    txid: str


@dataclass
class Rattus_norvegicus:
    genomeID: str
    code: str
    genomeVersion: str
    txid: str


@dataclass
class Oryctolagus_cuniculus:
    genomeID: str
    code: str
    genomeVersion: str
    txid: str


@dataclass
class Danio_rerio:
    genomeID: str
    code: str
    genomeVersion: str
    txid: str


@dataclass
class Drosophila_melanogaster:
    genomeID: str
    code: str
    genomeVersion: str
    txid: str


@dataclass
class Cricetulus_griseus:
    genomeID: str
    code: str
    genomeVersion: str
    txid: str


@dataclass
class Cavia_porcellus:
    genomeID: str
    code: str
    genomeVersion: str
    txid: str


@dataclass
class speciesData:
    Homo_sapiens: Homo_sapiens
    Mus_musculus: Mus_musculus
    Sus_scrofa: Sus_scrofa
    Equus_caballus: Equus_caballus
    Caenorhabditis_elegans: Caenorhabditis_elegans
    Felis_catus: Felis_catus
    Canis_lupus_familiaris: Canis_lupus_familiaris
    Bos_taurus: Bos_taurus
    Rattus_norvegicus: Rattus_norvegicus
    Oryctolagus_cuniculus: Oryctolagus_cuniculus
    Danio_rerio: Danio_rerio
    Drosophila_melanogaster: Drosophila_melanogaster
    Cricetulus_griseus: Cricetulus_griseus
    Cavia_porcellus: Cavia_porcellus




# Import these into your `__init__.py` file:
#
# from .parameters import generated_parameters
#
generated_parameters = {
    'sampleSheet': SnakemakeFileParameter(
        display_name='sampleSheet',
        type=LatchFile,
        config=True,
    ),
    'projectID': SnakemakeParameter(
        display_name='projectID',
        type=str,
        default='E-MTAB-6885',
    ),
    'projectComment': SnakemakeParameter(
        display_name='projectComment',
        type=str,
        default='Samples from ENA project PRJEB27261 ',
    ),
    'speciesName': SnakemakeParameter(
        display_name='speciesName',
        type=str,
        default='Homo sapiens',
    ),
    'deAnalysis': SnakemakeParameter(
        display_name='deAnalysis',
        type=int,
        default=1,
    ),
    'alpha': SnakemakeParameter(
        display_name='alpha',
        type=float,
        default=0.05,
    ),
    'includeSpikeIns': SnakemakeParameter(
        display_name='includeSpikeIns',
        type=int,
        default=0,
    ),
    'spikeInVersion': SnakemakeParameter(
        display_name='spikeInVersion',
        type=str,
        default="No"
    ),
    'includeSequence': SnakemakeParameter(
        display_name='includeSequence',
        type=int,
        default=0,
    ),
    'adapter': SnakemakeParameter(
        display_name='adapter',
        type=str,
        default='-a IlluminaSmallRNA3p=TGGAATTCTC',
    ),
    'adapterName': SnakemakeParameter(
        display_name='adapterName',
        type=str,
        default="Illumina Small RNA 3'",
    ),
    'readMinLength': SnakemakeParameter(
        display_name='readMinLength',
        type=int,
        default=17,
    ),
    'qualityCutoff': SnakemakeParameter(
        display_name='qualityCutoff',
        type=int,
        default=30,
    )
}
