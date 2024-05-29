import cmdlogtime
import os
import json
from subprocess import run

COMMAND_LINE_DEF_FILE = "./runBowtie2CommandLine.txt"

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
    # make directories
    os.makedirs(out_dir + "/fastqFileQC/" + config["histName1"], exist_ok=True)
    os.makedirs(out_dir + "/fastqFileQC/" + config["histName2"], exist_ok=True)
    
    # run fastqc
    print("Running fastqc")
    fastqccmd = "fastqc -o " + out_dir + "/fastqFileQC/" + config["histName1"] + " -f fastq " + config["projPath"] + "/" + config["fastq_PE1"]
    fastqccmd2 = "fastqc -o " + out_dir + "/fastqFileQC/" + config["histName1"] + " -f fastq " + config["projPath"] + "/" + config["fastq_PE2"]
    addl_logfile.write("\n\nfastqc: " + fastqccmd + "\n")
    result = run(fastqccmd, capture_output=True, text=True, shell=True)
    if (result.returncode):
        print("RC:", result.returncode, "\nOUT:", result.stdout, "\nERR:", result.stderr)
    addl_logfile.write("RC:" + str(result.returncode) + "\nOUT:" + result.stdout + "\nERR:" + result.stderr)
    addl_logfile.write("\n\nfastqc2: " + fastqccmd2 + "\n")
    result = run(fastqccmd2, capture_output=True, text=True, shell=True)
    if (result.returncode):
        print("RC:", result.returncode, "\nOUT:", result.stdout, "\nERR:", result.stderr)
        addl_logfile.write("RC:" + str(result.returncode) + "\nOUT:" + result.stdout + "\nERR:" + result.stderr)
    # get new file name
    filename1 = config["fastq_PE1"].split(".")[0]
    filename2 = config["fastq_PE2"].split(".")[0]
    
    # run fastp if needed
    if config["runfastp"]=="true":
        print("Running fastp")
        fastpcmd = ("fastp -i " + config["projPath"] + "/" + config["fastq_PE1"] + " -I " + config["projPath"] + "/" + config["fastq_PE2"] + 
        " -o " + config["projPath"] + "/" + filename1 + "_trimmed.fastq.gz" + " -O " + config["projPath"] + "/" + filename2 + "_trimmed.fastq.gz" +
        " -h " + out_dir + "/fastqFileQC/" + config["histName1"] + "/" + config["histName1"] + "_fastp.html")
        addl_logfile.write("\n\nfastp: " + fastpcmd + "\n")
        result = run(fastpcmd, capture_output=True, text=True, shell=True)
        if (result.returncode):
            print("RC:", result.returncode, "\nOUT:", result.stdout, "\nERR:", result.stderr)
            addl_logfile.write("RC:" + str(result.returncode) + "\nOUT:" + result.stdout + "\nERR:" + result.stderr)
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
        bowtie2cmd = ("nohup bowtie2 --local --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 700 -p " + config["cores"] + " -x " + 
                      config["projPath"] + config["ref"] + " -1 " + config["projPath"] + "/" + filename1 + "_trimmed.fastq.gz" + " -2 " + config["projPath"] + "/" + 
                      filename2 + "_trimmed.fastq.gz" + " -S " + out_dir + "/alignment/sam/" + config["histName1"] + "_bowtie2.sam &> " + out_dir + 
                      "/alignment/sam/bowtie2_summary/" + config["histName1"] + "_bowtie2.txt & > " + out_dir + "/alignment/" + config["histName1"] + "_bowtie2.out")
    else:
        bowtie2cmd = ("nohup bowtie2 --local --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 700 -p " + config["cores"] + " -x " + 
                      config["projPath"] + config["ref"] + " -1 " + config["projPath"] + "/" + config["fastq_PE1"] + " -2 " + config["projPath"] + "/" + 
                      config["fastq_PE2"] + " -S " + out_dir + "/alignment/sam/" + config["histName1"] + "_bowtie2.sam &> " + out_dir + 
                      "/alignment/sam/bowtie2_summary/" + config["histName1"] + "_bowtie2.txt & > " + out_dir + "/alignment/" + config["histName1"] + "_bowtie2.out")
    addl_logfile.write("\n\nbowtie2: " + bowtie2cmd + "\n")
    result = run(bowtie2cmd, capture_output=True, text=True, shell=True)
    if (result.returncode):
        print("RC:", result.returncode, "\nOUT:", result.stdout, "\nERR:", result.stderr)
    addl_logfile.write("RC:" + str(result.returncode) + "\nOUT:" + result.stdout + "\nERR:" + result.stderr)

    # run bowtie2 on spikein
    if config["spikein"]=="true":
        print("Running bowtie2 on spikein")
        if config["runfastp"]=="true":
            bowtie2spcmd = ("nohup bowtie2 --local --very-sensitive  --no-overlap --no-dovetail --no-mixed --no-discordant --phred33 -I 10 -X 700 -p " + config["cores"] + " -x " + 
                          config["projPath"] + config["spikeInRef"] + " -1 " + config["projPath"] + "/" + filename1 + "_trimmed.fastq.gz" + " -2 " + config["projPath"] + "/" + 
                          filename2 + "_trimmed.fastq.gz" + " -S " + out_dir + "/alignment/sam/" + config["histName1"] + "_spikein_bowtie2.sam &> " + out_dir + 
                          "/alignment/sam/bowtie2_summary/" + config["histName1"] + "_spikein_bowtie2.txt & > " + out_dir + "/alignment/" + config["histName1"] + "_spikein_bowtie2.out")
        else:
            bowtie2spcmd = ("nohup bowtie2 --local --very-sensitive  --no-overlap --no-dovetail --no-mixed --no-discordant --phred33 -I 10 -X 700 -p " + config["cores"] + " -x " + 
                          config["projPath"] + config["spikeInRef"] + " -1 " + config["projPath"] + "/" + config["fastq_PE1"] + " -2 " + config["projPath"] + "/" + 
                          config["fastq_PE2"] + " -S " + out_dir + "/alignment/sam/" + config["histName1"] + "_spikein_bowtie2.sam &> " + out_dir + 
                          "/alignment/sam/bowtie2_summary/" + config["histName1"] + "_spikein_bowtie2.txt & > " + out_dir + "/alignment/" + config["histName1"] + "_spikein_bowtie2.out")
        addl_logfile.write("\n\nbowtie2 spikein: " + bowtie2spcmd + "\n")
        result = run(bowtie2spcmd, capture_output=True, text=True, shell=True)
        if (result.returncode):
            print("RC:", result.returncode, "\nOUT:", result.stdout, "\nERR:", result.stderr)
        addl_logfile.write("RC:" + str(result.returncode) + "\nOUT:" + result.stdout + "\nERR:" + result.stderr)

    cmdlogtime.end(addl_logfile, start_time_secs)

if __name__ == "__main__":
    main()