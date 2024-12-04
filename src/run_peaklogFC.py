import cmdlogtime
import os
import json
from subprocess import run

COMMAND_LINE_DEF_FILE = "./getpeaksCommandLine.txt"

def main():
    (start_time_secs, pretty_start_time, my_args, addl_logfile) = cmdlogtime.begin(COMMAND_LINE_DEF_FILE)
    addl_logfile.write("starting log fold change and visualization run\n")
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
    projPath = config["projPath"]
    cores = config["cores"]
    gtfFile = config["gtfFile"]
    DEdir = projPath + "/" + "DESeq2_out/"
    # # get list of samples
    # # samples = []
    # # for key in config.keys():
    # #     if "SAMPLE" in key:
    # #         samples.append(config[key])
    # # histList = [sample["histName"] for sample in samples]
    # # sampleList=[]
    # # for sample in samples:
    # #     rep = sample["histName"]+"Rep"+str(sample["rep"])
    # #     sampleList.append(rep)
    # # print(sampleList)
    # # Convert bedgraph to bigwig file
    # print("converting bedgraph to bigwig")
    # bedgrtobwcmd = ("bedGraphToBigWig CtrlIgGvs.VipIgG.bedgraph chromSize CtrlIgGvs.VipIgG.bw")
    # addl_logfile.write("\n\nbedGraphToBigWig cmd: " + bedgrtobwcmd + "\n")
    # result = run(bedgrtobwcmd, check=True, capture_output=True, text=True, shell=True)
    # if (result.returncode):
    #     print("RC:", result.returncode, "\nOUT:", result.stdout, "\nERR:", result.stderr)
    # addl_logfile.write("RC:" + str(result.returncode) + "\nOUT:" + result.stdout + "\nERR:" + result.stderr)
    # "computeMatrix scale-regions -S CtrlIgGvs.VipIgG.bw -R /w5home/bmoore/genomes/human_genome_38/Homo_sapiens.GRCh38.108.gtf  --beforeRegionStartLength 5000 --regionBodyLength 5000 --afterRegionStartLength 5000 --skipZeros --missingDataAsZero -o CtrlIgGvs.VipIgG_gene_logFC_5k.mat.gz -p 16"
    # # Compute gene matrix from bigwig file
    # print("computing gene matrix")
    # computeMatrixcmd = ("computeMatrix scale-regions -S " + bowtie2_dir + "/alignment/bigwig/" + sample + "_norm.smooth.bw " +
    #                         "-R " + gtfFile  +
    #                         " --beforeRegionStartLength 5000 --regionBodyLength 5000 --afterRegionStartLength 5000 --skipZeros " +
    #                         "-o " + bowtie2_dir + "/peakCalling/SEACR/" + sample + "_gene_logFC_5k.mat.gz " +
    #                         "-p " + str(cores))
    # addl_logfile.write("\n\ncomputeMatrix cmd: " + computeMatrixcmd + "\n")
    # result = run(computeMatrixcmd, check=True, capture_output=True, text=True, shell=True)
    # if (result.returncode):
    #     print("RC:", result.returncode, "\nOUT:", result.stdout, "\nERR:", result.stderr)
    # addl_logfile.write("RC:" + str(result.returncode) + "\nOUT:" + result.stdout + "\nERR:" + result.stderr)
    # plot heatmap on logFC
    print("plotting heatmap")
    stat = "mean"
    for file in os.listdir(DEdir):
        if file.endswith("_gene_logFC_5k.mat.gz"):
            samplename = file.strip().split("_")[0]
            print(samplename)
            plotHeatmapcmd = ("plotHeatmap -m " + DEdir + file + " -o " + DEdir + samplename + "_heatmap.pdf " +
                        "--averageTypeSummaryPlot " + stat + " --sortUsing sum --heatmapHeight 16 --heatmapWidth 8 " +
                        "--outFileSortedRegions " + DEdir + samplename + "_gene_logFC_5k.bed")
            addl_logfile.write("\n\nplotHeatmap cmd: " + plotHeatmapcmd + "\n")
            result = run(plotHeatmapcmd, check=True, capture_output=True, text=True, shell=True)
            if (result.returncode):
                print("RC:", result.returncode, "\nOUT:", result.stdout, "\nERR:", result.stderr)
            addl_logfile.write("RC:" + str(result.returncode) + "\nOUT:" + result.stdout + "\nERR:" + result.stderr)

if __name__ == "__main__":
    main()