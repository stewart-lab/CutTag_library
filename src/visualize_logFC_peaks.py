
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

    # get chromSize file
    chromSize = config["chromSize"]
    # make bigwig files
    # print("Making bigwig files")
    # for file in os.listdir(bowtie2_dir + "/DESeq2_out/"):
    #     if file.endswith(".bedgraph"):
    #         samplename = file.strip().split(".bedgraph")[0]
    #         print(samplename)
    #         bwcmd = ("bedGraphToBigWig " + bowtie2_dir + "/DESeq2_out/" + file + " " + chromSize +
    #              " " + bowtie2_dir + "/DESeq2_out/" + samplename + ".bw")
    #         addl_logfile.write("\n\nbw cmd: " + bwcmd + "\n")
    #         result = run(bwcmd, check=True, capture_output=True, text=True, shell=True)
    #         if (result.returncode):
    #             print("RC:", result.returncode, "\nOUT:", result.stdout, "\nERR:", result.stderr)
    #         addl_logfile.write("RC:" + str(result.returncode)+ "\nOUT:" + result.stdout + "\nERR:" + result.stderr)

    cores = config["cores"]
    gtfFile = config["gtfFile"]
    
    # compute matrix on logFC comparisons
    print("Computing matrices")
    for file in os.listdir(bowtie2_dir + "/DESeq2_out/"):
        if file.endswith(".bw"):
            samplename = file.strip().split(".bw")[0]
            print(samplename)
            computeMatrixcmd = ("computeMatrix scale-regions -S " + bowtie2_dir + "/DESeq2_out/" + file +
                            " -R " + gtfFile  + " --missingDataAsZero "
                            " --beforeRegionStartLength 5000 --regionBodyLength 5000 --afterRegionStartLength 5000 --skipZeros " +
                            "-o " + bowtie2_dir + "/DESeq2_out/" + samplename + "_gene_logFC_5k.mat.gz " +
                            "-p " + str(cores))
            addl_logfile.write("\n\ncomputeMatrix cmd: " + computeMatrixcmd + "\n")
            result = run(computeMatrixcmd, check=True, capture_output=True, text=True, shell=True)
            if (result.returncode):
                print("RC:", result.returncode, "\nOUT:", result.stdout, "\nERR:", result.stderr)
            addl_logfile.write("RC:" + str(result.returncode) + "\nOUT:" + result.stdout + "\nERR:" + result.stderr)
            
    # plot heatmap
    print("Plotting heatmap")
    stat = "mean"
    for file in os.listdir(bowtie2_dir + "/DESeq2_out/"):
        if file.endswith("_gene_logFC_5k.mat.gz"):
            samplename = file.strip().split(".gz")[0]
            print(samplename)
            plotHeatmapcmd = ("plotHeatmap -m " + bowtie2_dir + "/DESeq2_out/" + file +
                          " -o " + bowtie2_dir + "/DESeq2_out/" + samplename + "_hm.pdf " +
                          "--averageTypeSummaryPlot " + stat + " --sortUsing sum --heatmapHeight 16 --heatmapWidth 8 --outFileSortedRegions " +
                          bowtie2_dir + "/DESeq2_out/" + samplename + ".bed")
            addl_logfile.write("\n\nplotHeatmap cmd: " + plotHeatmapcmd + "\n")
            result = run(plotHeatmapcmd, check=True, capture_output=True, text=True, shell=True)
            if (result.returncode):
                print("RC:", result.returncode, "\nOUT:", result.stdout, "\nERR:", result.stderr)
            addl_logfile.write("RC:" + str(result.returncode) + "\nOUT:" + result.stdout + "\nERR:" + result.stderr)

    print("finished!")
    cmdlogtime.end(addl_logfile, start_time_secs)

if __name__ == "__main__":
    main()