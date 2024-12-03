# mamba activate cut_tag

import cmdlogtime
import os
import json
from subprocess import run

COMMAND_LINE_DEF_FILE = "./getpeaksCommandLine.txt"

def main():
    (start_time_secs, pretty_start_time, my_args, addl_logfile) = cmdlogtime.begin(COMMAND_LINE_DEF_FILE)
    addl_logfile.write("starting visualization and peak finding run\n")
    out_dir = my_args["out_dir"]
    bowtie2_dir = my_args["bowtie2_dir"]

    # Get the directory where the script is located
    script_dir = os.path.dirname(os.path.abspath(__file__))

    # Open and read the config.json
    config_path = os.path.join(bowtie2_dir, "config.json")
    with open(config_path, "r") as f:
        config = json.load(f)
    print(config)
    # copy config to out_dir
    os.system("cp " + bowtie2_dir + "/config.json " + out_dir)

    # get bowtie2 directory
    #bowtie2_dir = config["bowtie2dir"]
    # get list of samples
    samples = []
    for key in config.keys():
        if "SAMPLE" in key:
            samples.append(config[key])
    histList = [sample["histName"] for sample in samples]
    sampleList=[]
    for sample in samples:
        rep = sample["histName"]+"Rep"+str(sample["rep"])
        sampleList.append(rep)
    print(sampleList)
    # run alignment_summary.r for alignment visualization
    print("Running alignment_summary.r")
    ralignsumcmd = ("Rscript --vanilla " + script_dir+"/alignment_summary.r " + bowtie2_dir + " " + out_dir)
    addl_logfile.write("\n\nR align cmd: " + ralignsumcmd + "\n")
    result = run(ralignsumcmd, check=True, capture_output=True, text=True, shell=True)
    if (result.returncode):
        print("RC:", result.returncode, "\nOUT:", result.stdout, "\nERR:", result.stderr)
    addl_logfile.write("RC:" + str(result.returncode) + "\nOUT:" + result.stdout + "\nERR:" + result.stderr)

    # get sequencing depth
    print("Getting sequencing depth")
    for sample in sampleList:
        seqdepthcmd = ("samtools view -F 0x04 " + bowtie2_dir + "/alignment/sam/" + sample + "_bowtie2.sam | wc -l")
        addl_logfile.write("\n\nseq depth cmd: " + seqdepthcmd + "\n")
        seqDepthDouble = os.popen(seqdepthcmd).read()
        seqDepth = int(seqDepthDouble)/2
        with open(bowtie2_dir + "/alignment/sam/bowtie2_summary/" + sample + "_bowtie2.seqDepth", "w") as f:
            f.write(str(seqDepth))
        print("Sequencing depth for " + sample + " is " + str(seqDepth))

    
    spikeIn = config["spikeIn"]
    if spikeIn == "true":
        print("Getting spike-in sequencing depth")
        for sample in sampleList:
            seqdepthcmd = ("samtools view -F 0x04 " + bowtie2_dir + "/alignment/sam/" + sample + "_bowtie2_spikeIn.sam | wc -l")
            addl_logfile.write("\n\nseq depth cmd: " + seqdepthcmd + "\n")
            seqDepthDouble = os.popen(seqdepthcmd).read()
            seqDepth2 = int(seqDepthDouble)/2
            with open(bowtie2_dir + "/alignment/sam/bowtie2_summary/" + sample + "_bowtie2_spikeIn.seqDepth", "w") as f:
                f.write(str(seqDepth2))
            print("Spike-in sequencing depth for " + sample + " is " + str(seqDepth2))

    # assess mapped fragment size distribution
    os.makedirs(bowtie2_dir + "/alignment/sam/fragmentLen", exist_ok=True)

    print("Assessing mapped fragment size distribution")
    for sample in sampleList:
        frsizecmd = ("samtools view -F 0x04 " + bowtie2_dir + "/alignment/sam/" + sample + "_bowtie2.sam | awk -F'\t'" +
        " 'function abs(x){return ((x < 0.0) ? -x : x)} {print abs($9)}' | sort | uniq -c | awk -v OFS='\t' '{print $2, $1/2}'")
        addl_logfile.write("\n\nfragment size cmd: " + frsizecmd + "\n")
        fragSizeDist = os.popen(frsizecmd).read()
        with open(bowtie2_dir + "/alignment/sam/fragmentLen/" + sample + "_fragmentLen.txt", "w") as f:
            f.write(fragSizeDist)
        print("Fragment size distribution for " + sample + " is " + fragSizeDist)

    # # run alignment_summary.r for alignment visualization
    print("Running fragment_length.r")
    rfrsizecmd = ("Rscript --vanilla " + script_dir+"/fragment_length.r " + bowtie2_dir + " " + out_dir)
    addl_logfile.write("\n\nR align cmd: " + rfrsizecmd + "\n")
    result = run(rfrsizecmd, check=True, capture_output=True, text=True, shell=True)
    if (result.returncode):
        print("RC:", result.returncode, "\nOUT:", result.stdout, "\nERR:", result.stderr)
    addl_logfile.write("RC:" + str(result.returncode) + "\nOUT:" + result.stdout + "\nERR:" + result.stderr)

    # filter mapped reads by minimal quality score
    minQualityScore = config["minQualityScore"]
    print("Filtering mapped reads by minimal quality score")
    for sample in sampleList:
        filtercmd = ("samtools view -h -q " + str(minQualityScore) + " " + bowtie2_dir + "/alignment/sam/" + sample + "_bowtie2.sam " +
        "> " + bowtie2_dir + "/alignment/sam/" + sample + "_bowtie2.qualityScore" + str(minQualityScore) + ".sam")
        addl_logfile.write("\n\nfilter cmd: " + filtercmd + "\n")
        result = run(filtercmd, check=True, capture_output=True, text=True, shell=True)
        if (result.returncode):
            print("RC:", result.returncode, "\nOUT:", result.stdout, "\nERR:", result.stderr)
        addl_logfile.write("RC:" + str(result.returncode) + "\nOUT:" + result.stdout + "\nERR:" + result.stderr)

    ## Filter and keep the mapped read pairs
    print("Filtering and keeping the mapped read pairs")
    for sample in sampleList:
        filtercmd = ("samtools view -bS -F 0x04 " + bowtie2_dir + "/alignment/sam/" + sample + "_bowtie2.qualityScore" + str(minQualityScore) + ".sam " +
        "> " + bowtie2_dir + "/alignment/bam/" + sample + "_bowtie2.mapped.bam")
        addl_logfile.write("\n\nfilter cmd: " + filtercmd + "\n")
        result = run(filtercmd, check=True, capture_output=True, text=True, shell=True)
        if (result.returncode):
            print("RC:", result.returncode, "\nOUT:", result.stdout, "\nERR:", result.stderr)
        addl_logfile.write("RC:" + str(result.returncode) + "\nOUT:" + result.stdout + "\nERR:" + result.stderr)

    ## Convert into bed file format
    print("Converting into bed file format")
    os.makedirs(bowtie2_dir + "/alignment/bed/", exist_ok=True)
    for sample in sampleList:
        bedcmd = ("bedtools bamtobed -i " + bowtie2_dir + "/alignment/bam/" + sample + "_bowtie2.mapped.bam -bedpe " +
        "> " + bowtie2_dir + "/alignment/bed/" + sample + "_bowtie2.bed")
        addl_logfile.write("\n\nbed cmd: " + bedcmd + "\n")
        result = run(bedcmd, check=True, capture_output=True, text=True, shell=True)
        if (result.returncode):
            print("RC:", result.returncode, "\nOUT:", result.stdout, "\nERR:", result.stderr)
        addl_logfile.write("RC:" + str(result.returncode) + "\nOUT:" + result.stdout + "\nERR:" + result.stderr)

    ## Keep the read pairs that are on the same chromosome and fragment length less than 1000bp.
    print("Keeping the read pairs that are on the same chromosome and fragment length less than 1000bp")
    for sample in sampleList:
        cleanbedcmd = ("awk '$1==$4 && $6-$2 < 1000 {print $0}' " + bowtie2_dir + "/alignment/bed/" + sample + "_bowtie2.bed " +
                       "> " + bowtie2_dir + "/alignment/bed/" + sample + "_bowtie2.clean.bed")
        addl_logfile.write("\n\nclean bed cmd: " + cleanbedcmd + "\n")
        result = run(cleanbedcmd, check=True, capture_output=True, text=True, shell=True)
        if (result.returncode):
            print("RC:", result.returncode, "\nOUT:", result.stdout, "\nERR:", result.stderr)
        addl_logfile.write("RC:" + str(result.returncode) + "\nOUT:" + result.stdout + "\nERR:" + result.stderr)

    ## Only extract the fragment related columns
    print("Extracting the fragment related columns")
    for sample in sampleList:
        fragcmd = ("cut -f 1,2,6 " + bowtie2_dir + "/alignment/bed/" + sample + "_bowtie2.clean.bed | sort -k1,1 -k2,2n -k3,3n " +
        "> " + bowtie2_dir + "/alignment/bed/" + sample + "_bowtie2.fragments.bed")
        addl_logfile.write("\n\nfrag cmd: " + fragcmd + "\n")
        result = run(fragcmd, check=True, capture_output=True, text=True, shell=True)
        if (result.returncode):
            print("RC:", result.returncode, "\nOUT:", result.stdout, "\nERR:", result.stderr)
        addl_logfile.write("RC:" + str(result.returncode) + "\nOUT:" + result.stdout + "\nERR:" + result.stderr)

    ## Assess fragment counts for each bin- normally 500bp 
    print("Assessing fragment counts for each bin")
    binLen = config["binLen"]
    for sample in sampleList:
        fragcountcmd = ("awk -v w=" + str(binLen) + " '{print $1, int(($2 + $3)/(2*w))*w + w/2}' " + bowtie2_dir + "/alignment/bed/" + sample + "_bowtie2.fragments.bed " +
        "| sort -k1,1V -k2,2n | uniq -c | awk -v OFS='\t' '{print $2, $3, $1}' | sort -k1,1V -k2,2n " +
        "> " + bowtie2_dir + "/alignment/bed/" + sample + "_bowtie2.fragmentsCount.bin" + str(binLen) + ".bed")
        addl_logfile.write("\n\nfrag count cmd: " + fragcountcmd + "\n")
        result = run(fragcountcmd, check=True, capture_output=True, text=True, shell=True)
        if (result.returncode):
            print("RC:", result.returncode, "\nOUT:", result.stdout, "\nERR:", result.stderr) 
        addl_logfile.write("RC:" + str(result.returncode) + "\nOUT:" + result.stdout + "\nERR:" + result.stderr)

    # assess spike-in calibration
    print("Assessing spike-in calibration")
    os.makedirs(bowtie2_dir + "/alignment/bedgraph/", exist_ok=True)
    spikeIn = config["spikeIn"]
    if spikeIn == "true":
        for sample in sampleList:
            seqDepth = float(open(bowtie2_dir + "/alignment/sam/bowtie2_summary/" + sample + "_bowtie2_spikeIn.seqDepth").read())
            seqDepth = int(seqDepth)
            if seqDepth > 1:
                scale_factor = 10000 / seqDepth
                print("Scaling factor for " + sample + " is " + str(scale_factor))
                scale_factor = round(scale_factor, 2)
                spikeInNormcmd = ("bedtools genomecov -bg -scale " + str(scale_factor) + " -i " + bowtie2_dir + "/alignment/bed/" 
                                  + sample + "_bowtie2.fragments.bed -g " + config["chromSize"] +
                                " > " + bowtie2_dir + "/alignment/bedgraph/" + sample + "_bowtie2.fragments.normalized.bedgraph")
                addl_logfile.write("\n\nspike-in normalization cmd: " + spikeInNormcmd + "\n")
                result = run(spikeInNormcmd, check=True, capture_output=True, text=True, shell=True)
                if (result.returncode):
                    print("RC:", result.returncode, "\nOUT:", result.stdout, "\nERR:", result.stderr)
                addl_logfile.write("RC:" + str(result.returncode) + "\nOUT:" + result.stdout + "\nERR:" + result.stderr)
            else:
                print("No spike-in reads found for " + sample, ", normalize to genome")
                spikeInNormcmd = ("bedtools genomecov -bg -i " + bowtie2_dir + "/alignment/bed/" + sample + "_bowtie2.fragments.bed -g " + 
                                  config["chromSize"] + " > " + bowtie2_dir + "/alignment/bedgraph/" + sample + 
                                  "_bowtie2.fragments.normalized.bedgraph")
                addl_logfile.write("\n\nspike-in normalization cmd: " + spikeInNormcmd + "\n")
                result = run(spikeInNormcmd, check=True, capture_output=True, text=True, shell=True)
                if (result.returncode):
                    print("RC:", result.returncode, "\nOUT:", result.stdout, "\nERR:", result.stderr)
                addl_logfile.write("RC:" + str(result.returncode) + "\nOUT:" + result.stdout + "\nERR:" + result.stderr)
    else:
        print("Spike-in calibration skipped, normalize to genome")
        for sample in sampleList:
            Normcmd = ("bedtools genomecov -bg -i " + bowtie2_dir + "/alignment/bed/" + sample + "_bowtie2.fragments.bed -g " + 
                   config["chromSize"] + " > " + bowtie2_dir + "/alignment/bedgraph/" + sample + 
                   "_bowtie2.fragments.normalized.bedgraph")
            addl_logfile.write("\n\nspike-in normalization cmd: " + Normcmd + "\n")
            result = run(Normcmd, check=True, capture_output=True, text=True, shell=True)
            if (result.returncode):
                print("RC:", result.returncode, "\nOUT:", result.stdout, "\nERR:", result.stderr)
            addl_logfile.write("RC:" + str(result.returncode) + "\nOUT:" + result.stdout + "\nERR:" + result.stderr)

    print("running scale.r")
    rscalecmd = ("Rscript --vanilla " + script_dir + "/scale.r " + bowtie2_dir + " " + out_dir)
    addl_logfile.write("\n\nR scale cmd: " + rscalecmd + "\n")
    result = run(rscalecmd, check=True, capture_output=True, text=True, shell=True)
    if (result.returncode):
        print("RC:", result.returncode, "\nOUT:", result.stdout, "\nERR:", result.stderr)
    addl_logfile.write("RC:" + str(result.returncode) + "\nOUT:" + result.stdout + "\nERR:" + result.stderr)

    os.makedirs(bowtie2_dir +  "/peakCalling/SEACR", exist_ok=True)
    print("Running SEACR peak calling")
    seacr = config["seacr"]
    for sample in sampleList:
        # run SEACR peak calling
        # because we don't have a control file, just select top 1% of regions by AUC
        seacr_cmd = (seacr + " " + bowtie2_dir + "/alignment/bedgraph/" + sample + "_bowtie2.fragments.normalized.bedgraph 0.01 non stringent " +
        bowtie2_dir + "/peakCalling/SEACR/" + sample + "_seacr_top0.01.peaks")
        addl_logfile.write("\n\nSEACR cmd: " + seacr_cmd + "\n")
        result = run(seacr_cmd, check=True, capture_output=True, text=True, shell=True)
        if (result.returncode):
            print("RC:", result.returncode, "\nOUT:", result.stdout, "\nERR:", result.stderr)
        addl_logfile.write("RC:" + str(result.returncode) + "\nOUT:" + result.stdout + "\nERR:" + result.stderr)

    # run peak_sum.r
    print("Running peak_sum.r")
    rpeaksumcmd = ("Rscript --vanilla " + script_dir + "/peak_sum.r " + bowtie2_dir + " " + out_dir)
    addl_logfile.write("\n\nR peak sum cmd: " + rpeaksumcmd + "\n")
    result = run(rpeaksumcmd, check=True, capture_output=True, text=True, shell=True)
    if (result.returncode):
        print("RC:", result.returncode, "\nOUT:", result.stdout, "\nERR:", result.stderr)
    addl_logfile.write("RC:" + str(result.returncode) + "\nOUT:" + result.stdout + "\nERR:" + result.stderr)
    
    # make bigwig files
    print("Making bigwig files")
    cores = config["cores"]
    os.makedirs(bowtie2_dir + "/alignment/bigwig", exist_ok=True)
    for sample in sampleList:
        # sort
        sortcmd = ("samtools sort -o " + bowtie2_dir + "/alignment/bam/" + sample + "_bowtie2.sorted.bam " +
        bowtie2_dir + "/alignment/bam/" + sample + "_bowtie2.mapped.bam")
        addl_logfile.write("\n\nsort cmd: " + sortcmd + "\n")
        result = run(sortcmd, check=True, capture_output=True, text=True, shell=True)
        if (result.returncode):
            print("RC:", result.returncode, "\nOUT:", result.stdout, "\nERR:", result.stderr)
        addl_logfile.write("RC:" + str(result.returncode) + "\nOUT:" + result.stdout + "\nERR:" + result.stderr)

        # index, -c for csi for genomes with chromosomes longer than 512Mb
        indexcmd = ("samtools index -c " + bowtie2_dir + "/alignment/bam/" + sample + "_bowtie2.sorted.bam")
        addl_logfile.write("\n\nindex cmd: " + indexcmd + "\n")
        result = run(indexcmd, check=True, capture_output=True, text=True, shell=True)
        if (result.returncode):
            print("RC:", result.returncode, "\nOUT:", result.stdout, "\nERR:", result.stderr)
        addl_logfile.write("RC:" + str(result.returncode) + "\nOUT:" + result.stdout + "\nERR:" + result.stderr)

        # convert to bw and normalize
        bwcmd = ("bamCoverage -b " + bowtie2_dir + "/alignment/bam/" + sample + "_bowtie2.sorted.bam " +
        "-o " + bowtie2_dir + "/alignment/bigwig/" + sample + "_norm.smooth.bw --normalizeUsing CPM " +
        "--binSize 500 --smoothLength 2000 --minMappingQuality 1 --numberOfProcessors " + str(cores))
        addl_logfile.write("\n\nbw cmd: " + bwcmd + "\n")
        result = run(bwcmd, check=True, capture_output=True, text=True, shell=True)
        if (result.returncode):
            print("RC:", result.returncode, "\nOUT:", result.stdout, "\nERR:", result.stderr)
        addl_logfile.write("RC:" + str(result.returncode)+ "\nOUT:" + result.stdout + "\nERR:" + result.stderr)

        # no smoothing
        bwcmd2 = ("bamCoverage -b " + bowtie2_dir + "/alignment/bam/" + sample + "_bowtie2.sorted.bam " +
        "-o " + bowtie2_dir + "/alignment/bigwig/" + sample + "_norm.bw --normalizeUsing CPM " +
        "--minMappingQuality 1 --numberOfProcessors " + str(cores))
        addl_logfile.write("\n\nnosmooth bw cmd: " + bwcmd2 + "\n")
        result = run(bwcmd2, check=True, capture_output=True, text=True, shell=True)
        if (result.returncode):
            print("RC:", result.returncode, "\nOUT:", result.stdout, "\nERR:", result.stderr)
        addl_logfile.write("RC:" + str(result.returncode) + "\nOUT:" + result.stdout + "\nERR:" + result.stderr)

    # use deeptools to compute matrix
    print("Computing matrix")

    projPath = config["projPath"]
    cores = config["cores"]
    gtfFile = config["gtfFile"]
    # on both marks
    # computeMatrixcmd = ("computeMatrix scale-regions -S " + bowtie2_dir + "/alignment/bigwig/" + sampleList[0] + "_norm.smooth.bw " +
    #                     bowtie2_dir + "/alignment/bigwig/" + sampleList[1] + "_norm.smooth.bw -R " +  
    #                     gtfFile + " --beforeRegionStartLength 5000 " +
    #                     "--regionBodyLength 5000 --afterRegionStartLength 10000 --skipZeros -o " + bowtie2_dir +
    #                     "/peakCalling/SEACR/" + histList[0] + "-" + histList[1] + "_gene_cpm_smooth10k.mat.gz " +
    #                     "-p " + str(cores))
    # addl_logfile.write("\n\ncomputeMatrix cmd: " + computeMatrixcmd + "\n")
    # result = run(computeMatrixcmd, check=True, capture_output=True, text=True, shell=True)
    # if (result.returncode):
    #     print("RC:", result.returncode, "\nOUT:", result.stdout, "\nERR:", result.stderr)
    # addl_logfile.write("RC:" + str(result.returncode) + "\nOUT:" + result.stdout + "\nERR:" + result.stderr)

    # on individual marks
    for sample in sampleList:
        computeMatrixcmd = ("computeMatrix scale-regions -S " + bowtie2_dir + "/alignment/bigwig/" + sample + "_norm.smooth.bw " +
                            "-R " + gtfFile  +
                            " --beforeRegionStartLength 5000 --regionBodyLength 5000 --afterRegionStartLength 5000 --skipZeros " +
                            "-o " + bowtie2_dir + "/peakCalling/SEACR/" + sample + "_gene_cpm_smooth10k.mat.gz " +
                            "-p " + str(cores))
        addl_logfile.write("\n\ncomputeMatrix cmd: " + computeMatrixcmd + "\n")
        result = run(computeMatrixcmd, check=True, capture_output=True, text=True, shell=True)
        if (result.returncode):
            print("RC:", result.returncode, "\nOUT:", result.stdout, "\nERR:", result.stderr)
        addl_logfile.write("RC:" + str(result.returncode) + "\nOUT:" + result.stdout + "\nERR:" + result.stderr)

    # plot heatmap
    print("Plotting heatmap")
    stat = "mean"
    # on both marks
    # plotHeatmapcmd = ("plotHeatmap -m " + bowtie2_dir + "/peakCalling/SEACR/" + histList[0] + "-" + histList[1] + "_gene_cpm_smooth10k.mat.gz " +
    #                   "-o " + bowtie2_dir + "/peakCalling/SEACR/AmexT_v47-AmexG_v6.0-DD_" + histList[0] + "-" + histList[1] + "_gene_cpm-smooth10k.pdf " +
    #                   "--averageTypeSummaryPlot " + stat + " --sortUsing sum --heatmapHeight 16 --heatmapWidth 8 --outFileSortedRegions " +
    #                   bowtie2_dir + "/peakCalling/SEACR/AmexT_v47-AmexG_v6.0-DD_" + histList[0] + "-" + histList[1] + "_gene.histone.cpm.smooth10k.bed")
    # addl_logfile.write("\n\nplotHeatmap cmd: " + plotHeatmapcmd + "\n")
    # result = run(plotHeatmapcmd, check=True, capture_output=True, text=True, shell=True)
    # if (result.returncode):
    #     print("RC:", result.returncode, "\nOUT:", result.stdout, "\nERR:", result.stderr)
    # addl_logfile.write("RC:" + str(result.returncode) + "\nOUT:" + result.stdout + "\nERR:" + result.stderr)

    # plot heatmap for individual marks
    for sample in sampleList:
        plotHeatmapcmd = ("plotHeatmap -m " + bowtie2_dir + "/peakCalling/SEACR/" + sample + "_gene_cpm_smooth10k.mat.gz " +
                          "-o " + bowtie2_dir + "/peakCalling/SEACR/" + sample + "_gene_cpm-smooth10k.pdf " +
                          "--averageTypeSummaryPlot " + stat + " --sortUsing sum --heatmapHeight 16 --heatmapWidth 8 --outFileSortedRegions " +
                          bowtie2_dir + "/peakCalling/SEACR/" + sample + "_gene.histone.cpm.smooth10k.bed")
        addl_logfile.write("\n\nplotHeatmap cmd: " + plotHeatmapcmd + "\n")
        result = run(plotHeatmapcmd, check=True, capture_output=True, text=True, shell=True)
        if (result.returncode):
            print("RC:", result.returncode, "\nOUT:", result.stdout, "\nERR:", result.stderr)
        addl_logfile.write("RC:" + str(result.returncode) + "\nOUT:" + result.stdout + "\nERR:" + result.stderr)

    # Heatmap on peaks
    # get summit region
    print("Getting summit region")
    for sample in sampleList:
        summitcmd = ("awk '{split($6, summit, \":\"); split(summit[2], region, \"-\"); print summit[1]\"\\t\"region[1]\"\\t\"region[2]}' " +
                     bowtie2_dir + "/peakCalling/SEACR/" + sample + "_seacr_top0.01.peaks.stringent.bed " +
                     "> " + bowtie2_dir + "/peakCalling/SEACR/" + sample + "_seacr_top0.01.peaks.summitRegion.bed")
        addl_logfile.write("\n\nsummit cmd: " + summitcmd + "\n")
        result = run(summitcmd, check=True, capture_output=True, text=True, shell=True)
        if (result.returncode):
            print("RC:", result.returncode, "\nOUT:", result.stdout, "\nERR:", result.stderr)
        addl_logfile.write("RC:" + str(result.returncode) + "\nOUT:" + result.stdout + "\nERR:" + result.stderr)

    # compute matrix
    print("Computing matrix an plotting heatmap for each sample")
    for sample in sampleList:
        computeMatrixcmd = ("computeMatrix reference-point -S " + bowtie2_dir + "/alignment/bigwig/" + sample + "_norm.smooth.bw " +
                            "-R " + bowtie2_dir + "/peakCalling/SEACR/" + sample + "_seacr_top0.01.peaks.summitRegion.bed " +
                            "--skipZeros -o " + bowtie2_dir + "/peakCalling/SEACR/" + sample + "_SEACR.mat.gz -p " + str(cores) +
                            " -a 3000 -b 3000 --referencePoint center")
        addl_logfile.write("\n\ncomputeMatrix cmd: " + computeMatrixcmd + "\n")
        result = run(computeMatrixcmd, check=True, capture_output=True, text=True, shell=True)
        if (result.returncode):
            print("RC:", result.returncode, "\nOUT:", result.stdout, "\nERR:", result.stderr)
        addl_logfile.write("RC:" + str(result.returncode) + "\nOUT:" + result.stdout + "\nERR:" + result.stderr)
        plotHeatmapcmd = ("plotHeatmap -m " + bowtie2_dir + "/peakCalling/SEACR/" + sample + "_SEACR.mat.gz -out " +
                          bowtie2_dir + "/peakCalling/SEACR/" + sample + "_SEACR_heatmap.pdf --sortUsing sum --startLabel " +
                          "'Peak Start' --endLabel 'Peak End' --xAxisLabel '' --regionsLabel 'Peaks' --samplesLabel " + sample +
                          " --averageTypeSummaryPlot " + stat)
        addl_logfile.write("\n\nplotHeatmap cmd: " + plotHeatmapcmd + "\n")
        result = run(plotHeatmapcmd, check=True, capture_output=True, text=True, shell=True)
        if (result.returncode):
            print("RC:", result.returncode, "\nOUT:", result.stdout, "\nERR:", result.stderr)
        addl_logfile.write("RC:" + str(result.returncode) + "\nOUT:" + result.stdout + "\nERR:" + result.stderr)


    cmdlogtime.end(addl_logfile, start_time_secs)

if __name__ == "__main__":
    main()