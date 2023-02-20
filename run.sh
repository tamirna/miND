#!/bin/bash
source ~/.bashrc

MINDVERSION="1.2 RC1"
# read command line arguments

POSITIONAL=()
count=0
while [[ $# -gt 0 ]]
do
  key="$1"
  case $key in
      -i|--input)
      (( count++ ))
      INPUT="$2"
      shift
      shift
      ;;
      -)
      (( count++ ))
      STDIN=YES
      shift
      ;;
      -o|--output)
      (( count++ ))
      OUTPUT="$2"
      shift
      shift
      ;;
      -k|--keeptmp)
      (( count++ ))
      KEEPTMP=YES
      shift
      ;;
      --help)
      (( count++ ))
      HELP=YES
      shift
      ;;
      --version)
      (( count++ ))
      VERSION=YES
      shift
      ;;
      *)    # unknown option
      (( count++ ))
      POSITIONAL+=("$1") # save it in an array for later
      shift
      ;;
  esac
done
set -- "${POSITIONAL[@]}" # restore positional parameters

# set to read from stdin if no arguments are passed
if (( ${count} == 0 )); then
  STDIN=YES
fi

if [[ "${STDIN}" == "YES" ]] || [[ "${INPUT}" == "-" ]]; then
  INPUT=$(</dev/stdin)
fi

if [[ -n $1 ]]; then
  INPUT=$1
fi
if [[ -n $2 ]]; then
  OUTPUT=$2
fi

if [[ "${HELP}" == "YES" ]]; then
  echo "Usage: ./run.sh [OPTION]... [FILE]..."
  echo "Run MIND pipeline based on configuration in the -i FILE and generate output"
  echo ""
  echo "With no FILE, or when FILE is -, read configurtaion file path from standard input."
  echo ""
  echo "  -i, --input       Path to configuration SampleContrastSheet.xlsx file"
  echo "  -o, --output      Output subfolder name"
  echo "  -k, --keeptmp     Keep temporary outputfiles (like trimmed fastq files)"
  echo "      --help        display this help and exit"
  echo ""
  echo "MIND pipeline online help: <https://www.tamirna.com/mind>"
  exit 0
fi

if [[ "${VERSION}" == "YES" ]]; then
  echo "MIND ${MINDVERSION}"
  echo "Copyright (C) 2021 TAmiRNA GmbH"
  echo "License GNU AGPLv3: GNU Affero General Public License v3.0 or later <https://www.gnu.org/licenses/agpl-3.0.html>."
  echo "This is free software: you are free to change and redistribute it."
  echo "There is NO WARRANTY, to the extent permitted by law."
  echo ""
  echo "Written and developed by Andreas B. Diendorfer"
  exit 0
fi

if [[ -f "~/miniconda3/etc/profile.d/conda.sh" ]]; then
    source ~/miniconda3/etc/profile.d/conda.sh
fi

# run snakemake and archive the workflow
if [[ $(type -P "conda") ]]; then
  printf "Conda path is set\n"
else
  printf "Conda path no set, adding to PATH variable\n"
  wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh
  bash ~/miniconda.sh -u -b
  export PATH="$HOME/miniconda3/bin:$PATH"
  ~/miniconda3/bin/conda init
fi

if ! [[ $(type -P "conda") ]]; then
  printf "Conda binary not found in path\nexiting\n"
  exit 1
fi

ENVNAME="mind-workflow-env"

if conda env list | awk '{print $1}' | grep -xqFe "$ENVNAME"; then
   source activate $ENVNAME
else
   echo "Creating conda environment $ENVNAME"
   conda create -y --name $ENVNAME
   conda activate $ENVNAME
   conda install -y mamba -c conda-forge
fi;
echo "Current conda env: $CONDA_DEFAULT_ENV"

if hash wc 2>/dev/null; then
  echo "wc command found"
else
  printf "wc binary not found in path\nexiting\n"
  exit 1
fi

if hash snakemake 2>/dev/null; then
  echo "snakemake command found"
else
  echo "snakemake command NOT found"
  conda install -y mamba -c conda-forge
  mamba install -y -c bioconda -c conda-forge snakemake=6.*
  mamba install -y -c conda-forge pip xorg-libxrender xorg-libxpm python=3.9.* pandas=1.1.* xlrd=1.*
fi

if hash dot 2>/dev/null; then
  echo "dot command found"
else
  echo "dot command NOT found"
  mamba install -y -c conda-forge graphviz=2.42.*
fi

if hash aria2c 2>/dev/null; then
  echo "aria2c command found"
else
  echo "aria2c command NOT found"
  mamba install -y -c conda-forge aria2
fi

# read condaPath from config file
CONDA_PREFIX=`grep 'condaPath:' config.yaml | tail -n1 | awk '{gsub(/"/, "", $2); print $2}'`

# check command line arguments
cmdlineargs="--config sampleSheet='${INPUT}'"
if [[ "${OUTPUT}" != "" ]]; then
 cmdlineargs="${cmdlineargs} outputSubfolder='${OUTPUT}'"
fi
if [[ "${KEEPTMP}" == "YES" ]]; then
 cmdlineargs="${cmdlineargs} --notemp"
fi

echo "Unlocking snakemake"
snakemake --unlock $cmdlineargs
echo "Running pipeline"
snakemake --use-conda --conda-prefix=$CONDA_PREFIX --conda-frontend=mamba --restart-times=1 --cores $cmdlineargs

exit 0
