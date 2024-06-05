# mamba activate cut_tag

import cmdlogtime
import os
import json
from subprocess import run

COMMAND_LINE_DEF_FILE = "./runBowtie2CommandLine.txt"


def run_fastqc(out_dir, histName, projPath, fastq, addl_logfile):
    os.makedirs(out_dir + "/fastqFileQC/" + histName, exist_ok=True)
    fastqccmd = "fastqc -o " + out_dir + "/fastqFileQC/" + histName + " -f fastq " + projPath + "/" + fastq
    addl_logfile.write("\n\nfastqc: " + fastqccmd + "\n")
    result = run(fastqccmd, capture_output=True, text=True, shell=True)
    if (result.returncode):
        print("RC:", result.returncode, "\nOUT:", result.stdout, "\nERR:", result.stderr)
    addl_logfile.write("RC:" + str(result.returncode) + "\nOUT:" + result.stdout + "\nERR:" + result.stderr)


def run_fastp(out_dir, histName, rep, projPath, fastq1, fastq2, filename1, filename2, addl_logfile):
    os.makedirs(out_dir + "/fastqFileQC/" + histName, exist_ok=True)
    fastpcmd = ("fastp -i " + projPath + "/" + fastq1 + " -I " + projPath + "/" + fastq2 + " -o " + 
                projPath + "/" + filename1 + " -O " + projPath  + "/" + filename2 + 
                " -h " + out_dir + "/fastqFileQC/" + histName + "/" + histName + "_rep" + rep + "_fastp.html"
                + " -j " + out_dir + "/fastqFileQC/" + histName + "/" + histName + "_rep" + rep + "_fastp.json")
    addl_logfile.write("\n\nfastp: " + fastpcmd + "\n")
    result = run(fastpcmd, capture_output=True, text=True, shell=True)
    if (result.returncode):
        print("RC:", result.returncode, "\nOUT:", result.stdout, "\nERR:", result.stderr)
    addl_logfile.write("RC:" + str(result.returncode) + "\nOUT:" + result.stdout + "\nERR:" + result.stderr)


def run_bowtie2(out_dir, histName, rep, projPath, ref, filename1, filename2, cores, addl_logfile):
    bowtie2cmd = ("nohup bowtie2 --local --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 700 -p " 
                  + cores + " -x " + projPath + ref + " -1 " + projPath + "/" + filename1 + " -2 " + projPath + "/" +
                  filename2 + " -S " + out_dir + "/alignment/sam/" + histName  + "_rep" + rep + "_bowtie2.sam &> " 
                  + out_dir + "/alignment/sam/bowtie2_summary/" + histName  + "_rep" + rep +  "_bowtie2.txt & > " + 
                  out_dir + "/alignment/" + histName  + "_rep" + rep +  "_bowtie2.out")
    addl_logfile.write("\n\nbowtie2: " + bowtie2cmd + "\n")
    result = run(bowtie2cmd, capture_output=True, text=True, shell=True)
    if (result.returncode):
        print("RC:", result.returncode, "\nOUT:", result.stdout, "\nERR:", result.stderr)
    addl_logfile.write("RC:" + str(result.returncode) + "\nOUT:" + result.stdout + "\nERR:" + result.stderr)


def run_bowtie2sp(out_dir, histName, rep, projPath, ref, filename1, filename2, cores, addl_logfile):
    bowtie2spcmd = ("nohup bowtie2 --local --very-sensitive  --no-overlap --no-dovetail --no-mixed --no-discordant " + 
                    "--phred33 -I 10 -X 700 -p " + cores + " -x " + projPath + ref + " -1 " + projPath + "/" + 
                    filename1 + " -2 " + projPath + "/" + filename2 + " -S " + out_dir + "/alignment/sam/" + histName + 
                    "_rep" + rep + "_bowtie2_spikeIn.sam &> " + out_dir + "/alignment/sam/bowtie2_summary/" + histName + 
                    "_rep" + rep +"_bowtie2_spikeIn.txt & > " + out_dir + "/alignment/" + histName + "_rep" + rep +
                    "_spikein_bowtie2.out")
    addl_logfile.write("\n\nbowtie2 spikein: " + bowtie2spcmd + "\n")
    result = run(bowtie2spcmd, capture_output=True, text=True, shell=True)
    if (result.returncode):
        print("RC:", result.returncode, "\nOUT:", result.stdout, "\nERR:", result.stderr)
    addl_logfile.write("RC:" + str(result.returncode) + "\nOUT:" + result.stdout + "\nERR:" + result.stderr)


def main():
    (start_time_secs, pretty_start_time, my_args, addl_logfile) = cmdlogtime.begin(COMMAND_LINE_DEF_FILE)
    addl_logfile.write("Starting Bowtie2 run\n")
    out_dir = my_args["out_dir"]
    
    # Get the directory where the script is located
    script_dir = os.path.dirname(os.path.abspath(__file__))

    # Open and read the config.json
    config_path = os.path.join(script_dir, "config.json")
    with open(config_path, "r") as f:
        config = json.load(f)
    print(config)
    # copy config to out_dir
    os.system("cp " + script_dir + "/config.json " + out_dir)
    
    # get list of samples
    samples = []
    for key in config.keys():
        if "SAMPLE" in key:
            samples.append(config[key])
    print(samples)

    # run fastqc
    print("Running fastqc")
    for sample in samples:
        run_fastqc(out_dir, sample["histName"], config["projPath"], sample["fastq_PE1"], addl_logfile)
        run_fastqc(out_dir, sample["histName"], config["projPath"], sample["fastq_PE2"], addl_logfile)
    
    # run fastp if needed
    if config["runfastp"]=="true":
        print("Running fastp")
        for sample in samples:
            filename1 = sample["fastq_PE1"].split(".")[0]+ "_trimmed.fastq.gz"
            filename2 = sample["fastq_PE2"].split(".")[0]+ "_trimmed.fastq.gz"
            run_fastp(out_dir, sample["histName"], sample["rep"], config["projPath"], sample["fastq_PE1"], 
                      sample["fastq_PE2"], filename1, filename2, addl_logfile)

    else:
        print("fastp skipped")
    
    # make directories for bowtie2 output
    os.makedirs(out_dir + "/alignment/sam/bowtie2_summary", exist_ok=True)
    os.makedirs(out_dir + "/alignment/bam", exist_ok=True)
    os.makedirs(out_dir + "/alignment/bed", exist_ok=True)
    os.makedirs(out_dir + "/alignment/bedgraph", exist_ok=True)

    # run bowtie2 with parameters from config
    print("Running bowtie2 on reference genome")
    if config["runfastp"]=="true":
        for sample in samples:
            filename1 = sample["fastq_PE1"].split(".")[0]+ "_trimmed.fastq.gz"
            filename2 = sample["fastq_PE2"].split(".")[0]+ "_trimmed.fastq.gz"
            run_bowtie2(out_dir, sample["histName"], sample["rep"], config["projPath"], config["ref"], filename1, filename2, config["cores"], addl_logfile)
    else:
        for sample in samples:
            run_bowtie2(out_dir, sample["histName"], sample["rep"], config["projPath"], config["ref"], sample["fastq_PE1"], sample["fastq_PE2"], config["cores"], addl_logfile)

    # run bowtie2 on spikein
    if config["spikeIn"]=="true":
        print("Running bowtie2 on spikein")
        if config["runfastp"]=="true":
            for sample in samples:
                filename1 = sample["fastq_PE1"].split(".")[0]+ "_trimmed.fastq.gz"
                filename2 = sample["fastq_PE2"].split(".")[0]+ "_trimmed.fastq.gz"
                run_bowtie2sp(out_dir, sample["histName"], sample["rep"], config["projPath"], config["spikeInRef"], filename1, filename2, config["cores"], addl_logfile)
        else:
            for sample in samples:
                run_bowtie2sp(out_dir, sample["histName"], sample["rep"], config["projPath"], config["spikeInRef"], sample["fastq_PE1"], sample["fastq_PE2"], config["cores"], addl_logfile)
    print("Bowtie2 running, check bowtie summary files for results\n")
    # edit config file to add out_dir as bowtie2dir
    config["bowtie2dir"] = out_dir
    # write config file
    # new config path
    config_path = os.path.join(out_dir, "config.json")
    with open(config_path, "w") as f:
        json.dump(config, f, indent=4)

    cmdlogtime.end(addl_logfile, start_time_secs)

if __name__ == "__main__":
    main()