import pandas as pd
import pathlib
import xlrd
xlrd.xlsx.ensure_elementtree_imported(False, None)
xlrd.xlsx.Element_has_iter = True

configfile: "config.yaml"

sampleSheet = config['sampleSheet']
sampleSheetPathObj = pathlib.Path(sampleSheet)

sampleSheetPath = str(sampleSheetPathObj.parent) # folder containing the sample sheet
sampleSheetName = str(sampleSheetPathObj.name) # filename of the sample sheet without path
sampleSheetStem = str(sampleSheetPathObj.stem) # filename without extension

# read configuration XLS sheet
xlsConfig = {}
xls = pd.read_excel(sampleSheet, header = None, index_col = 0, sheet_name = 'Project Details', converters={0:str,1:str}, usecols = "A,B")

for index, row  in xls.itertuples():
    xlsConfig[index] = row

if('Project ID' in xlsConfig and isinstance(xlsConfig['Project ID'], str)): config['projectID'] = xlsConfig['Project ID']
else: config['projectID'] = 'none'

if('Project comment' in xlsConfig and isinstance(xlsConfig['Project comment'], str)): config['projectComment'] = xlsConfig['Project comment']
else: config['projectComment'] = ""

# set default to hsa
config['speciesName'] = "Homo sapiens"
config['speciesTxid'] = config['speciesData'][config['speciesName']]['txid']
config['speciesCode'] = config['speciesData'][config['speciesName']]['code']
config['genomeID'] = config['speciesData'][config['speciesName']]['genomeID']
config['genomeVersion'] = config['speciesData'][config['speciesName']]['genomeVersion']

if xlsConfig['Species']:
    if xlsConfig['Species'] in config['speciesData']:
        config['speciesName'] = xlsConfig['Species']
        config['speciesCode'] = config['speciesData'][xlsConfig['Species']]['code']
        config['speciesTxid'] = config['speciesData'][xlsConfig['Species']]['txid']
        config['genomeID'] = config['speciesData'][xlsConfig['Species']]['genomeID']
        config['genomeVersion'] = config['speciesData'][xlsConfig['Species']]['genomeVersion']
    else:
        sys.exit("Set species is not valid")

if(xlsConfig['Differential expression analysis'] == "Yes"):
    config['deAnalysis'] = 1
else:
    config['deAnalysis'] = 0

if('Significance level' in xlsConfig and isinstance(xlsConfig['Significance level'], str)): config['alpha'] = xlsConfig['Significance level']
else: config['alpha'] = "0.05"

if(xlsConfig['Spikein analysis'] != "No"):
    config['includeSpikeIns'] = 1
    config['spikeInVersion'] = xlsConfig['Spikein analysis']
else:
    config['includeSpikeIns'] = 0

if(xlsConfig['Report spikein sequences'] == "Yes"):
    config['includeSequence'] = 1
else:
    config['includeSequence'] = 0

if(xlsConfig['Adapter'] == "Illumina Universal Adapter"):
    config['adapter'] = "-a IlluminaUniversal=AGATCGGAAGAG"
    config['adapterName'] = xlsConfig['Adapter']
elif (xlsConfig['Adapter'] == "Illumina Small RNA 3'"):
    config['adapter'] = "-a IlluminaSmallRNA3p=TGGAATTCTC"
    config['adapterName'] = xlsConfig['Adapter']
elif (xlsConfig['Adapter'] == "RealSeq"):
    config['adapter'] = "-a Realseq3P=TGGAATTCTC -u1"
    config['adapterName'] = xlsConfig['Adapter']
else:
    config['adapter'] = xlsConfig['Custom adapter (cutadapt cfg)']
    config['adapterName'] = ""

if('Minimum read length' in xlsConfig and isinstance(xlsConfig['Minimum read length'], str)): config['readMinLength'] = xlsConfig['Minimum read length']
else: config['readMinLength'] = "17"

if('Reads quality cutoff' in xlsConfig and isinstance(xlsConfig['Reads quality cutoff'], str)): config['qualityCutoff'] = xlsConfig['Reads quality cutoff']
else: config['qualityCutoff'] = "30"

#read samples from XLS sheet

def get_ids(sampleSheet):
    xlsSamples = []
    xls = pd.read_excel(sampleSheet, index_col = 0, sheet_name = 'Sample Group Matrix', converters={0:str,1:str,2:str}, usecols = "A:G")
    for index, row in xls.iterrows():
        filename = row["Filename"]
        if isinstance(filename, str):
            # remove file extensions should there be any
            filename = filename.rstrip(".fastq.gz")
            filename = filename.rstrip(".fq.gz")
            filename = filename.rstrip(".fastq")
            filename = filename.rstrip(".fq")
            filename = filename.rstrip(".bam")
            xlsSamples.append(filename)

    if(len(xlsSamples) > 0):
        IDS = xlsSamples
    else:
        print("No input files found!")
        sys.exit()

    return IDS

if 'outputSubfolder' in config and config['outputSubfolder'] == 'date':
    from datetime import datetime
    outputSubfolder = datetime.now().strftime("%d%m%y_%H%M%S")
elif 'outputSubfolder' in config and str(config['outputSubfolder']) != '':
    outputSubfolder = str(config['outputSubfolder']).strip("/").replace("/", "_")
else:
    if 'sampleSheetStem' == 'SampleContrastSheet':
        outputSubfolder = 'default'
    else:
        outputSubfolder = sampleSheetStem.replace("-SampleContrastSheet", "").replace("_SampleContrastSheet", "")

outPath = "output/" + outputSubfolder

def getInput_all(wildcards):
    inputList = []
    inputList.append(expand("%s/analysis/report.html" % outPath))
    return inputList

rule all:
    input:
        unpack(getInput_all)
    shell:
        '''
            rm -fr {outPath}/rundata
            mkdir {outPath}/rundata
            echo '{sampleSheetPath}' > {outPath}/rundata/inputfiles.txt
            ls -l "{sampleSheetPath}" >> {outPath}/rundata/inputfiles.txt
            cp {sampleSheet} {outPath}/rundata/
            tar cfz {outPath}/rundata/snakemake.tar.gz Snakefile run.sh config.yaml scripts/
            rm -fr dir_mapper_seq* mapper_logs Rplots.pdf

            # generate report zip file
            cp {outPath}/analysis/report.html "{outPath}/$(date -I)_{config[projectID]}_MINDreport.html"
            zip -q -j "{outPath}/$(date '+%Y-%m-%d')_{config[projectID]}_{config[speciesCode]}_{config[alpha]}_MINDreport.zip" "{outPath}/$(date -I)_{config[projectID]}_MINDreport.html" {outPath}/multiqc/multiqc_report.html "{sampleSheet}"
            rm "{outPath}/$(date -I)_{config[projectID]}_MINDreport.html"

            echo ""
            echo "####################################################"
            echo "#          miND - miRNA NGS data pipeline          #"
            echo "#         Copyright (C) 2021 TAmiRNA GmbH          #"
            echo "#  Written and developed by Andreas B. Diendorfer  #"
            echo "#                                                  #"
            echo "#                  RUN FINISHED!                   #"
            echo "####################################################"
            echo ""
            echo "Your results file (including the HTML report and multiQC details) is located at"
            echo "{outPath}/$(date '+%Y-%m-%d')_{config[projectID]}_{config[speciesCode]}_{config[alpha]}_MINDreport.zip"
            echo ""
            echo "The HTML report can also be found at"
            echo "{outPath}/analysis/report.html"
            echo ""
        '''

def getInput_makeHTMLReport(wildcards):
    inputList = {}
    inputList['sampleSheetFile'] = sampleSheet
    inputList['mirnaMappingStats'] = "%s/analysis_sampleMappingStats/all_samples.csv" % outPath
    inputList['mirnaMappingStatsRPM'] = "%s/analysis_sampleMappingStats/all_samples_rpm_miRNA.csv" % outPath
    inputList['sampleMappingQuantifications'] = "%s/analysis_quantifySampleMappings/all_samples.csv" % outPath
    inputList['mirnaMappingStatsRPMLib'] = "%s/analysis_sampleMappingStats/all_samples_rpm_lib.csv" % outPath
    inputList['mirnaMappingStatsRPMDistDat'] = "%s/analysis_sampleMappingStats/all_samples_rpm_dist.dat" % outPath
    inputList['sampleMappingQuantificationsDat'] = "%s/analysis_quantifySampleMappings/all_samples.dat" % outPath
    inputList['sampleMappingQuantificationsDatPerc'] = "%s/analysis_quantifySampleMappings/all_samples_perc.dat" % outPath
    inputList['sampleMappingQuantificationsCsvPerc'] = "%s/analysis_quantifySampleMappings/all_samples_perc.csv" % outPath
    inputList['csvGenomeMapping'] = "%s/analysis_quantifySampleMappings/all_samples_genome_mapping.csv" % outPath
    inputList['multiQC'] = "%s/multiqc/multiqc_report.html" % outPath

    if config["deAnalysis"]:
        inputList['deDir'] = "%s/analysis_differentialExpression/" % outPath
    return inputList

rule makeHTMLReport:
    input:
        unpack(getInput_makeHTMLReport)
    output:
        "%s/analysis/report.html" % outPath
    log:    "%s/logs/analysis/report.log" % outPath
    conda:  "ranalysis"
    threads: 1
    script:
        "scripts/report.Rmd"

rule runDEAnalysis:
    input:
        readsFile = "%s/analysis_sampleMappingStats/all_samples.csv" % outPath,
        sampleSheetFile = sampleSheet
    output:
        directory("%s/analysis_differentialExpression/" % outPath)
    params:
        alpha=config["alpha"]
    log:    "%s/logs/analysis/differentialExpression.log" % outPath
    conda:  "ranalysis"
    threads: 1
    shell:
        "sleep 3 && Rscript scripts/differentialExpression.R --outDir='{output}' --readsFile='{input.readsFile}' --sampleSheetFile='{input.sampleSheetFile}' --alpha='{params.alpha}' > {log} 2>&1"


rule quantifySampleMappings:
    input:
        sampleSheetFile = sampleSheet,
        fastqLibSize = expand("%s/fastq_trimmed/{filename}.fastq.libsize" % outPath, filename=get_ids(sampleSheet)),
        spikeins = expand("%s/mapping_bowtie_spikeins/{filename}.map" % outPath, filename=get_ids(sampleSheet)),
        genome = expand("%s/mapping_bowtie_genome/{filename}.map" % outPath, filename=get_ids(sampleSheet)),
        mirna = expand("%s/mapping_bowtie_mirna/{filename}.map" % outPath, filename=get_ids(sampleSheet)),
        rnacentral = expand("%s/mapping_bowtie_rnacentral/{filename}.map" % outPath, filename=get_ids(sampleSheet)),
        cdna = expand("%s/mapping_bowtie_cdna/{filename}.map" % outPath, filename=get_ids(sampleSheet)),
        mappedSeqSizeDist = expand("%s/mapping_bowtie_genome/{filename}.mapped.ssd" % outPath, filename=get_ids(sampleSheet)),
        unmappedSeqSizeDist = expand("%s/mapping_bowtie_genome/{filename}.unmapped.ssd" % outPath, filename=get_ids(sampleSheet))
    output: csv = "%s/analysis_quantifySampleMappings/all_samples.csv" % outPath,
            dat = "%s/analysis_quantifySampleMappings/all_samples.dat" % outPath,
            csvPerc = "%s/analysis_quantifySampleMappings/all_samples_perc.csv" % outPath,
            datPerc = "%s/analysis_quantifySampleMappings/all_samples_perc.dat" % outPath,
            csvGenomeMapping = "%s/analysis_quantifySampleMappings/all_samples_genome_mapping.csv" % outPath
    log:    "%s/logs/analysis/quantifySampleMappings.log" % outPath
    conda:  "ranalysis"
    threads: config['threads']['high']
    shell:
        "sleep 3 && Rscript scripts/quantifySampleMappings.R --outFile='{output.csv}' --outFileGenomeMapping='{output.csvGenomeMapping}' --sampleSheetFile='{input.sampleSheetFile}' --cores={threads} {input.genome} 2> {log}"

rule getFastQLibSize:
    input:  "%s/fastq_trimmed/{filename}.fastq" % outPath
    output: "%s/fastq_trimmed/{filename}.fastq.libsize" % outPath
    log:    "%s/logs/getFastQLibSize/{filename}.log" % outPath
    threads: 1
    shell:
        "wc -l {input} | awk '{{print $1/4}}' > {output} 2> {log}"

rule combineSampleMappingStats:
    input:
        miRNA = "%s/mapping_mirdeep2_miRNA/{filename}/miRNAs_expressed_all_samples_default.csv" % outPath,
        spikeins = "%s/mapping_bbduk_spikeins_core/{filename}.txt" % outPath
    output: "%s/analysis_sampleMappingStats/{filename}.csv" % outPath
    log:    "%s/logs/analysis/sampleMappingStats/{filename}.log" % outPath
    conda:  "ranalysis"
    threads: 1
    shell:
        "sleep 3 && Rscript scripts/combineSampleMappingStats.R --outFile='{output}' --includeSpikeIns={config[includeSpikeIns]} --input.miRNA='{input.miRNA}' --input.spikeins='{input.spikeins}' 2> {log}"


rule combineAllMappingStats:
    input:  csv = expand("%s/analysis_sampleMappingStats/{filename}.csv" % outPath, filename=get_ids(sampleSheet)),
            sampleSheetFile = sampleSheet
    output: csv = "%s/analysis_sampleMappingStats/all_samples.csv" % outPath,
            csvRPMmiRNA = "%s/analysis_sampleMappingStats/all_samples_rpm_miRNA.csv" % outPath,
            csvRPMLib = "%s/analysis_sampleMappingStats/all_samples_rpm_lib.csv" % outPath,
            plot = "%s/analysis_sampleMappingStats/all_samples_rpm_dist.dat" % outPath
    log:    "%s/logs/analysis/combineAllMappingStats.log" % outPath
    conda:  "ranalysis"
    threads: 1
    shell:
        "sleep 3 && Rscript scripts/combineAllMappingStats.R --outFile='{output.csv}' --includeSequence={config[includeSequence]}  --sampleSheetFile='{input.sampleSheetFile}' {input.csv} 2> {log}"

rule sortBowtieMapping:
    input: "%s/{mapping}/{filename}.map" % outPath
    output: temp("%s/{mapping}/{filename}.sorted.map" % outPath)
    log:    "%s/logs/{mapping}/{filename}.log" % outPath
    threads: 1
    shell:
        "cat '{input}' | sort -nr -k2,2 -t \"x\" > '{output}' 2> {log}"

rule mappingBowtieSpikeIns:
    input:  "%s/mapping_mirdeep2_miRNA/{filename}.fasta" % outPath
    output: map = "%s/mapping_bowtie_spikeins/{filename}.map" % outPath,
    #        mapped = temp("%s/mapping_bowtie_spikeins/{filename}.mapped.fasta" % outPath),
            unmapped = temp("%s/mapping_bowtie_spikeins/{filename}.unmapped.fasta" % outPath)
    threads: config['threads']['medium']
    log:    "%s/logs/mapping_bowtie_spikeins/{filename}.log" % outPath
    conda:
        "bowtie1"
    shell:
        "bowtie --threads {threads} -f -k1 --fullref --best -v0 --norc --un '{output.unmapped}' libs/spikeins/spikeins_full '{input}' > '{output.map}' 2> {log}"


rule mappingBowtieGenome:
    input:  "%s/mapping_bowtie_spikeins/{filename}.unmapped.fasta" % outPath
    output: map = "%s/mapping_bowtie_genome/{filename}.map" % outPath,
            mapped = temp("%s/mapping_bowtie_genome/{filename}.mapped.fasta" % outPath),
            unmapped = temp("%s/mapping_bowtie_genome/{filename}.unmapped.fasta" % outPath),
            mappedSeqSizeDist = "%s/mapping_bowtie_genome/{filename}.mapped.ssd" % outPath,
            unmappedSeqSizeDist = "%s/mapping_bowtie_genome/{filename}.unmapped.ssd" % outPath
    log:    "%s/logs/mapping_bowtie_genome/{filename}.log" % outPath
    conda:  "bowtie1"
    threads: config['threads']['medium']
    shell:
        """
            latch cp latch://23349.account/snakemake/mind/{config[repoPath]}/{config[genomeID]}/{config[genomeVersion]} {config[genomeVersion]} &&\
            bowtie --threads {threads} -f -k1 -v2 --fullref --un '{output.unmapped}' --al '{output.mapped}' {config[genomeVersion]}/bowtiedb/genome '{input}' > '{output.map}' 2> {log}
            cat '{output.mapped}' | awk '{{if(NR%2==1) {{printf "%s\t",$0}} else {{printf "%i\\n",length($1)}} }}' | cut -d "x" -f2 | awk '{{i[$2]+=$1}} END{{for(x in i){{print i[x]" "x}}}}' | sort -k2 -n > {output.mappedSeqSizeDist} 2> {log}
            cat '{output.unmapped}' | awk '{{if(NR%2==1) {{printf "%s\t",$0}} else {{printf "%i\\n",length($1)}} }}' | cut -d "x" -f2 | awk '{{i[$2]+=$1}} END{{for(x in i){{print i[x]" "x}}}}' | sort -k2 -n > {output.unmappedSeqSizeDist} 2> {log}
        """

rule mappingBowtieMirna:
    input:  "%s/mapping_bowtie_genome/{filename}.mapped.fasta" % outPath
    output: map = "%s/mapping_bowtie_mirna/{filename}.map" % outPath,
            # mapped = temp("%s/mapping_bowtie_mirna/{filename}.mapped.fasta" % outPath),
            unmapped = temp("%s/mapping_bowtie_mirna/{filename}.unmapped.fasta" % outPath)
    threads: config['threads']['medium']
    log:    "%s/logs/mapping_bowtie_mirna/{filename}.log" % outPath
    conda:
        "bowtie1"
    shell:
        """
        latch cp latch://23349.account/snakemake/mind/repository/data/mirbase mirbase &&\
        bowtie --threads {threads} -f -k1 --fullref --best -v1 --un '{output.unmapped}' mirbase/{config[miRBaseVersion]}/hairpin/bowtiedb/hairpin-{config[speciesCode]} '{input}' > '{output.map}' 2> {log}
        """

rule mappingBowtieRNAcentral:
    input:  "%s/mapping_bowtie_mirna/{filename}.unmapped.fasta" % outPath
    output: map = "%s/mapping_bowtie_rnacentral/{filename}.map" % outPath,
            # mapped = temp("%s/mapping_bowtie_rnacentral/{filename}.mapped.fasta" % outPath),
            unmapped = temp("%s/mapping_bowtie_rnacentral/{filename}.unmapped.fasta" % outPath)
    threads: config['threads']['medium']
    log:    "%s/logs/mapping_bowtie_rnacentral/{filename}.log" % outPath
    conda:
        "bowtie1"
    shell:
        """
        latch cp latch://23349.account/snakemake/mind/{config[repoPath]}/rnacentral/{config[rnacentralVersion]}/bowtiedb bowtiedb &&\
        bowtie --threads {threads} -f -k1 --fullref --best -v1 --norc --un '{output.unmapped}' bowtiedb/rnacentral_species_specific_ids-{config[speciesTxid]} '{input}' > '{output.map}' 2> {log}
        """

rule mappingBowtieCDNA:
    input:  "%s/mapping_bowtie_rnacentral/{filename}.unmapped.fasta" % outPath
    output: map = "%s/mapping_bowtie_cdna/{filename}.map" % outPath,
            # mapped = temp("%s/mapping_bowtie_cdna/{filename}.mapped.fasta" % outPath),
            unmapped = temp("%s/mapping_bowtie_cdna/{filename}.unmapped.fasta" % outPath)
    threads: config['threads']['medium']
    log:    "%s/logs/mapping_bowtie_cdna/{filename}.log" % outPath
    conda:
        "bowtie1"
    shell:
        """
        latch cp latch://23349.account/snakemake/mind/{config[repoPath]}/{config[genomeID]}/{config[genomeVersion]} {config[genomeVersion]} &&\
        bowtie --threads {threads} -f -k1 --fullref --best -v1 --norc --un '{output.unmapped}' {config[genomeVersion]}/bowtiedb/cdna '{input}' > '{output.map}' 2> {log}
        """

rule miRDeepPrep:
    input: "%s/fastq_trimmed/{filename}.fastq" % outPath
    output:
        collapsedReads = temp("%s/mapping_mirdeep2_miRNA/{filename}.fasta" % outPath)
#        readsVsGenome = temp("%s/mapping_mirdeep2_miRNA/{filename}.arf" % outPath)
    threads: config['threads']['medium']
    log:    "%s/logs/mapping_mirdeep2_miRNA/{filename}.log" % outPath
    conda:
        "mirdeep2"
    shell:
        """
            mapper.pl '{input}' -e -h -j -l %s -o {threads} -m -s '{output.collapsedReads}' > '{log}' 2>&1
            rm -fr mapper_logs
            rm -fr bowtie.log
        """ % config['readMinLength']

rule mappingMiRDeep2MiRNA:
    input:
        collapsedReads = "%s/mapping_bowtie_genome/{filename}.mapped.fasta" % outPath
        # readsVsGenome = "%s/miRDeep2/genomemapped/{filename}.arf"
    output:
        map = "%s/mapping_mirdeep2_miRNA/{filename}/miRNAs_expressed_all_samples_default.csv" % outPath
    params:
        outputDir = "%s/mapping_mirdeep2_miRNA/{filename}" % outPath
    threads: config['threads']['medium']
    log:    "%s/logs/mapping_mirdeep2_miRNA/{filename}.log" % outPath
    conda:
        "mirdeep2"
    shell:
        """
        latch cp latch://23349.account/snakemake/mind/repository/data/mirbase mirbase &&\
        subDirs=$(awk -F"/" '{{print NF-1}}' <<< "{input.collapsedReads}")
        subDirPath=$(seq ${{subDirs}} | awk '{{printf "../"}}')
        mkdir -p '{params.outputDir}'
        cd '{params.outputDir}'
        quantifier.pl -p '{workflow.basedir}/mirbase/{config[miRBaseVersion]}/hairpin/uncompressed/hairpin-{config[speciesCode]}.dna.fa' \
            -y 'default' \
            -T {threads} \
            -d \
            -m '{workflow.basedir}/mirbase/{config[miRBaseVersion]}/mature/uncompressed/mature-{config[speciesCode]}.dna.fa' \
            -r "../${{subDirPath}}{input.collapsedReads}" > "../${{subDirPath}}{log}" 2>&1
        rm -fr expression_analyses/expression_analyses_default/*.ebwt
        """

rule mappingBBDukSpikeInsCore:
    input:  "%s/fastq_trimmed/{filename}.fastq" % outPath
    output: stats = temp("%s/mapping_bbduk_spikeins_core/{filename}.txt" % outPath)
    #        mapped = temp("%s/mapping_bbduk_spikeins_core/{filename}.mapped.fastq" % outPath)
    threads: config['threads']['high']
    log:    "%s/logs/mapping_bbduk_spikeins_core/{filename}.log" % outPath
    conda:
        "bbmap"
    shell:
        "bbduk.sh threads={threads} -Xmx1024m in='{input}' outm=stdout.fq ref='libs/spikeins/spikeins_core.fa' stats='{output.stats}' statscolumns=5 k=13 maskmiddle=f rcomp=f hdist=0 edist=0 rename=t > /dev/null 2> {log}"

rule mappingBBDukSpikeIns:
    input:  "%s/fastq_trimmed/{filename}.fastq" % outPath
    output: stats = temp("%s/mapping_bbduk_spikeins/{filename}.txt" % outPath)
    #        mapped = temp("%s/mapping_bbduk_spikeins/{filename}.mapped.fastq" % outPath)
    threads: config['threads']['high']
    log:    "%s/logs/mapping_bbduk_spikeins/{filename}.log" % outPath
    conda:
        "bbmap"
    shell:
        "bbduk.sh threads={threads} -Xmx1024m in='{input}' outm=stdout.fq ref='libs/spikeins/spikeins.fa' stats='{output.stats}' statscolumns=5 k=21 maskmiddle=f rcomp=f hdist=0 edist=0 rename=t > /dev/null 2> {log}"

rule filterCollapsedReads:
    input: "%s/fastq_collapsed/{filename}.fastq" % outPath
    output:temp("%s/fastq_collapsed_filtered/{filename}.fastq" % outPath)
    threads: 1
    run:
        shell("awk 'BEGIN {{FS = \"_x\" ; OFS = \"\\n\"}} {{header = $0 ; getline seq ; getline qheader ; getline qseq ; if ($2 >= {config[minReadCount]}) {{print header, seq, qheader, qseq}}}}' < '{input}' > '{output}'")

rule collapseReadsFastq:
    input:  "%s/fastq_trimmed/{filename}.fastq" % outPath
    output: temp("%s/fastq_collapsed/{filename}.fastq" % outPath)
    log:    "%s/logs/fastq_collapsed/{filename}.log" % outPath
    params:
        outPath = outPath
    conda:  "seqcluster"
    threads: 1
    shell: """
        seqcluster collapse -f '{input}' -o {params.outPath}/fastq_collapsed/{wildcards.filename}/ > {log} 2>&1
        rename 's/_trimmed//' {params.outPath}/fastq_collapsed/{wildcards.filename}/*.fastq
        mv {params.outPath}/fastq_collapsed/{wildcards.filename}/*.fastq {params.outPath}/fastq_collapsed/ &&
        rm -fr {params.outPath}/fastq_collapsed/{wildcards.filename}/
    """

rule uncompressing:
    input: sampleSheetPath + "/{filename}.{ext}.gz"
    output: temp("%s/fastq_raw/{filename}.{ext}" % outPath)
    threads: config['threads']['low']
    log: "%s/logs/fastq_raw/{filename}.{ext}.log" % outPath
    shell:
        "gzip -dc '{input}' > '{output}' 2> {log}"

rule moveFqToFastq:
    input: "%s/fastq_raw/{filename}.fq" % outPath
    output: temp("%s/fastq_raw/{filename}.fastq" % outPath)
    threads: config['threads']['low']
    log: "%s/logs/fastq_raw/mv_{filename}.log" % outPath
    shell:
        "mv '{input}' '{output}' 2> {log}"

rule copyFq2ToFastq:
    input: sampleSheetPath + "/{filename}.fq"
    output: temp("%s/fastq_raw/{filename}.fastq" % outPath)
    threads: config['threads']['low']
    log: "%s/logs/fastq_raw/cp_{filename}.log" % outPath
    shell:
        "cp '{input}' '{output}' 2> {log}"

rule copyFastqToFastq:
    input: sampleSheetPath + "/{filename}.fastq"
    output: temp("%s/fastq_raw/{filename}.fastq" % outPath)
    threads: config['threads']['low']
    log: "%s/logs/fastq_raw/cp_{filename}.log" % outPath
    shell:
        "cp '{input}' '{output}' 2> {log}"

rule bam2fastq:
    input: sampleSheetPath + "/{filename}.bam"
    output: temp("%s/fastq_raw/{filename}.fastq" % outPath)
    threads: config['threads']['low']
    log: "%s/logs/fastq_raw/{filename}.log" % outPath
    conda:  "samtools"
    shell:
        "samtools bam2fq {input} > {output} 2> {log}"

rule multiQC:
    input:
        fastqc = expand("%s/{type}/{filename}_fastqc.zip" % outPath, type={'fastqc_raw', 'fastqc_trimmed'}, filename=get_ids(sampleSheet)),
        cutadapt = expand("%s/logs/fastq_trimmed/{filename}.log" % outPath, filename=get_ids(sampleSheet))
    output: "%s/multiqc/multiqc_report.html" % outPath
    log:    "%s/multiqc/multiqc.log" % outPath
    params:
        outPath = outPath
    conda:  "multiqc"
    threads: 1
    shell:
        """
        rm -fr '{params.outPath}/multiqc/' && mkdir '{params.outPath}/multiqc/'
        multiqc '{params.outPath}/fastqc_raw/' '{params.outPath}/fastqc_trimmed/' '{params.outPath}/logs/fastq_trimmed/' --config multiqc_config.yaml -o '{params.outPath}/multiqc/' > '{log}' 2>&1
        """

rule fastQC:
    input:  "%s/fastq_{mod}/{filename}.fastq" % outPath
    output:
        html = "%s/fastqc_{mod}/{filename}_fastqc.html" % outPath,
        zip = temp("%s/fastqc_{mod}/{filename}_fastqc.zip" % outPath)
    log:    "%s/logs/fastqc_{mod}/{filename}.log" % outPath
    params:
        outPath = outPath
    conda:  "fastqc"
    threads: config['threads']['low']
    shell:
        """
        fastqc '{input}' --outdir='{params.outPath}/fastqc_{wildcards.mod}/' > {log} 2>&1
        """

rule adapterTrimming:
    input:  "%s/fastq_raw/{filename}.fastq" % outPath
    output: temp("%s/fastq_trimmed/{filename}.fastq" % outPath)
    log:    "%s/logs/fastq_trimmed/{filename}.log" % outPath
    conda:  "cutadapt"
    threads: config['threads']['high']
    shell:
        "cutadapt {config[adapter]} --minimum-length {config[readMinLength]} --quality-cutoff {config[qualityCutoff]} --discard-untrimmed --cores={threads} -o '{output}' '{input}' > {log} 2>&1"
