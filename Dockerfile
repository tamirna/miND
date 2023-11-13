# DO NOT CHANGE
from 812206152185.dkr.ecr.us-west-2.amazonaws.com/latch-base:fe0b-main

workdir /tmp/docker-build/work/

shell [ \
    "/usr/bin/env", "bash", \
    "-o", "errexit", \
    "-o", "pipefail", \
    "-o", "nounset", \
    "-o", "verbose", \
    "-o", "errtrace", \
    "-O", "inherit_errexit", \
    "-O", "shift_verbose", \
    "-c" \
]
env TZ='Etc/UTC'
env LANG='en_US.UTF-8'

arg DEBIAN_FRONTEND=noninteractive

# Latch SDK
# DO NOT REMOVE
run mkdir /opt/latch

# Install Mambaforge
run apt-get update --yes && \
    apt-get install --yes curl && \
    curl \
        --location \
        --fail \
        --remote-name \
        https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh && \
    `# Docs for -b and -p flags: https://docs.anaconda.com/anaconda/install/silent-mode/#linux-macos` \
    bash Mambaforge-Linux-x86_64.sh -b -p /opt/conda -u && \
    rm Mambaforge-Linux-x86_64.sh

# Set conda PATH
env PATH=/opt/conda/bin:$PATH

# Build conda environment
copy environment.yaml /opt/latch/environment.yaml
run mamba env create \
    --file /opt/latch/environment.yaml \
    --name workflow
env PATH=/opt/conda/envs/workflow/bin:$PATH

# Latch SDK
# DO NOT REMOVE

copy envs/* /opt/latch/
run mamba env create --file /opt/latch/bbmap.yml
run mamba env create --file /opt/latch/bowtie1.yml
run mamba env create --file /opt/latch/cutadapt.yml
run mamba env create --file /opt/latch/fastqc.yml
run mamba env create --file /opt/latch/mirdeep2.yml
run mamba env create --file /opt/latch/multiqc.yml
run mamba env create --file /opt/latch/ranalysis.yml
run mamba env create --file /opt/latch/samtools.yml
run mamba env create --file /opt/latch/seqcluster.yml
run mamba env create --file /opt/latch/seqtk.yml

# Copy workflow data (use .dockerignore to skip files)
copy . .latch/* /root/

run cd /root/envs/latch && pip install -e .

# Latch workflow registration metadata
# DO NOT CHANGE
arg tag
# DO NOT CHANGE
env FLYTE_INTERNAL_IMAGE $tag

workdir /root
