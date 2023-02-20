#!/usr/bin/env bash

source ~/.bashrc

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

source ~/miniconda3/etc/profile.d/conda.sh

ENVNAME="mind-repository-env"

if conda env list | awk '{print $1}' | grep -xqFe "$ENVNAME"; then
   source activate $ENVNAME
else
   echo "Creating conda environment $ENVNAME"
   conda create -y --name $ENVNAME
   conda activate $ENVNAME
   conda install -y mamba -c conda-forge
fi;
echo "Current conda env: $CONDA_DEFAULT_ENV"

if hash snakemake 2>/dev/null; then
  echo "snakemake command found"
else
  echo "snakemake command NOT found"
  conda install -y mamba -c conda-forge
  mamba install -y -c bioconda -c conda-forge snakemake=6.*
  mamba install -y -c conda-forge pip xorg-libxrender xorg-libxpm python=3.9.* pandas=1.1.* xlrd=1.*
fi

if hash aria2c 2>/dev/null; then
  echo "aria2c command found"
else
  echo "aria2c command NOT found"
  mamba install -y -c conda-forge aria2
fi

# `find . -type d`
for source in `find . -type d`
do
  if [ -f "$source/Snakefile" ]; then
    cd $source
    echo "Snakemake $source"
    snakemake --unlock
    snakemake --use-conda --conda-frontend=mamba --conda-prefix=../envs/tmp  --restart-times=1 --cores --resources network=100
    cd ../
  fi
done
