#!/usr/bin/python
## Example subprocess push to cluster job submission in Python

import sys, os, re, subprocess, os.path, datetime
sys.path.append('/home/bst/bin/anaconda2/lib/python2.7')
import argparse

# ncores = "8"
# email = "bst@stowers.org"
# pe = "by_node"
# queue = "all.q"
# molng_path = "/home/bst/Random/qsubing/test_data/"
# bowtie_index = "/n/data1/genomes/bowtie-index/mm10/mm10"

# aligner = "tophat"
# extra_commands = ""
# get_cwd = "/home/bst/Random/rletest/test2"
# sge_out = get_cwd + "/SGE_out"
# bam_dir = get_cwd + "/bam"
# tophat_dir = get_cwd + "/tophat_aln"

# things to add
# add paired end compatibility - yup
# add exclusion of samples
# add subselection of samples
# add help
# add final summary script that runs after indexing!
# auto find user for e-mailing

aligner = "tophat"
space_string = ">>>                                                    <<<\n"


def test_existance(theobject):
    isreal = os.path.isfile(theobject)
    if (isreal is False):
        sys.stderr.write(space_string)
        sys.exit("ERROR: " + theobject + " doesn't exist!   \n")

def analyze_molng(molng_path,args):
    indexed_sample_reports = {}
    ##sample_report_file = molng_path + "/Sample_Report.csv"
    for flowcell in molng_path:
        ## print(flowcell)
        sample_report_file = flowcell + "/Sample_Report.csv"
        test_existance(sample_report_file)
        sample_from_exp  = os.listdir(os.path.dirname(sample_report_file))
        sample_from_exp = filter(lambda x:re.search(r'[ACGT].fastq.gz$', x), sample_from_exp)
        sample_report_file_handle = open(sample_report_file, "r")
        sample_reports = [line.split(",") for line in sample_report_file_handle.read().splitlines()]
        sample_report_file_handle.close()
        for line in sample_reports:                
            if any(sample in sample_from_exp for sample in line):
                sample_fastq = line[0]
                sample_name = line[4]
                sample_fastq = flowcell + sample_fastq              
                if(len(line) < 17):
                ## catch wrong Sample_Report.csv from older LIMS
                    sys.stderr.write(space_string)
                    sys.stderr.write(">>>   WELL... WELL.. WELL... THE NUMBER OF COLUMNS IS NOT CORRECT!\nUSING the OLD LIMS ARE WE?   \n")
                    sys.stderr.write(space_string)
                    sys.exit("The number of columns in Sample_Report.csv is only " + str(len(line)) + " when it should be 17   \n")
                if(line[8] is "2" and args.paired_end is False):
                ## catch paired end
                    sys.stderr.write(space_string)
                    sys.stderr.write(">>>   WOAH BRO GET THOSE PAIRED END READS OUTTA HERE!   \n")
                    sys.stderr.write(space_string)
                    sys.exit("ERROR: " + sample_name + " is a twin!   \n")
                if sample_name not in indexed_sample_reports:
                    indexed_sample_reports[sample_name] = [sample_fastq]
                else:
                    indexed_sample_reports[sample_name].append(sample_fastq)
        if len(indexed_sample_reports) is 0:
            sys.stderr.write(space_string)
            sys.exit("ERROR: Looks like none of the .fastq.gz's in your directory were matched in the sample report!   \n")
        else:
            sys.stderr.write(space_string)
            sys.stderr.write(">>>   Sample reports generated successfully/   \n")
    indexed_sample_reports = indexed_sample_reports
    return(indexed_sample_reports)


def filter_paired_read(files,num):
    files_r = filter(lambda x:re.search(r'_' + num + '_[AGCT]*.[ACGT]*.fastq.gz$', x), files)
    if len(files_r) > 1:
        files_r = ','.join(files_r)
    else:
        files_r = ''.join(files_r)
    return(files_r)


def transcript_test_index(args):
    if args.bowtie_one is True:
        use_bowtie_one = "--bowtie1 "
    else:
        use_bowtie_one = ""
    if args.gtf_index is "":
        use_gtf = ""
    else:
        use_gtf = " -G " + args.gtf_index
        if os.path.isfile(args.gtf_index) is False:
            sys.stderr.write(">>>   GFF3/GTF FILE NOT FOUND in " + args.gtf_index + "    \n")
    if args.transcriptome_index is '':
        index_command = "/n/local/bin/" + aligner + " " + use_gtf + " " + use_bowtie_one + args.bowtie_index
    else:
        index_command = "/n/local/bin/" + aligner + " --transcriptome-index " + args.transcriptome_index + use_gtf + " " + use_bowtie_one + args.bowtie_index
        if use_gtf is "":
            if not os.path.exists(args.transcriptome_index + ".ver"):
                sys.stderr.write("No GFF3/GTF file provided   \n")
                sys.stderr.write("No transcriptome detected in the provided path: " + os.path.dirname(args.transcriptome_index) + "   \n")
                sys.stderr.write(space_string)
                sys.exit("ERROR: Unable to create new transcriptome without annotation file     \n")
        if not os.path.exists(os.path.dirname(args.transcriptome_index)):
            os.makedirs(os.path.dirname(args.transcriptome_index))            
    full_command = "'" + index_command + "'"
    sys.stderr.write(">>>   Checking for your transcriptome-index: \n")
    sys.stderr.write(space_string)
    if not os.path.exists(args.transcriptome_index + ".ver"):
        sys.stderr.write("No transcriptome detected in " + os.path.dirname(args.transcriptome_index) + "   \n")
        sys.stderr.write("Remember! In order to create a transcriptome index you need to have the genome fasta file in teh SAME directory as the gtf/gff!   \n")
        sys.stderr.write(space_string)
    else:
        sys.stderr.write("Previous transcriptome detected!   \n")
        sys.stderr.write(space_string)
    processors = "1"
    sge_out = args.output_dir + "/SGE_out"
    ##if args.gtf_index is not "":
    job_name = "transcriptome_index_" + os.path.basename(args.transcriptome_index)
    if args.email is None:
        job_string = """-m ae -N %s -pe %s %s -q %s -l mem_free=1G -j y -o %s -V -cwd -terse -b y %s""" % (job_name, args.pe, processors, args.queue, sge_out, full_command)
    else:
        job_string = """-M %s -m ae -N %s -pe %s %s -q %s -l mem_free=1G -j y -o %s -V -cwd -terse -b y %s""" % (args.email, job_name, args.pe, processors, args.queue, sge_out, full_command)
    server_command = """qsub %s""" % (job_string)
    if args.gtf_index is not "":
        sys.stderr.write(space_string)
        sys.stderr.write(">>>   Transcriptome Indexing Command:   \n")
        sys.stderr.write(server_command + "\n")
        sys.stderr.write(space_string)
        align_job = subprocess.Popen(server_command,shell=True,stdout=subprocess.PIPE)
        align_job_id = align_job.stdout.read().rstrip()
        sys.stderr.write(space_string)
        sys.stderr.write(">>>   Checking the tx index and building if necessary!   \n")
        sys.stderr.write(space_string)
    else:
        align_job_id = None
    return align_job_id


def tophat_commander(sample_report, sample, args):
    sample_name = sample
    fastq_files = sample_report.get(sample_name)
    bam_dir = args.output_dir + "/bam"
    tophat_dir = args.output_dir + "/tophat_aln"
    if args.paired_end is False:
        if len(fastq_files) > 1:
            fastq_files = ','.join(fastq_files)
        else:
            fastq_files = ''.join(fastq_files)
    else:
        fastq_files_r1 = filter_paired_read(fastq_files,'1')
        fastq_files_r2 = filter_paired_read(fastq_files,'2')
        if len(fastq_files_r2) is 0:
            sys.stderr.write(space_string)
            sys.stderr.write(">>>   WOAH BRO THIS DOESN'T LOOK LIKE PAIRED END DATA!   \n")
            sys.stderr.write(space_string)
            sys.exit("ERROR: " + sample_name + " hase no pair!   \n")
        fastq_files = fastq_files_r1 + " " + fastq_files_r2
    if args.bowtie_one is True:
        use_bowtie_one = "--bowtie1 "
        sys.stderr.write(">>>   --use-bowtie1 option included. Testing for bowtie1 index!   \n")
        test_existance((args.bowtie_index + ".1.ebwt"))
    else:
        use_bowtie_one = ""
        sys.stderr.write(">>>   Testing for bowtie2 index! To use Tophat with a bowtie1 index you MUST specify --use-bowtie1.   \n")
        small_bt = os.path.isfile((args.bowtie_index + ".1.bt2"))
        large_bt = os.path.isfile((args.bowtie_index + ".1.bt2l"))
        if small_bt or large_bt:
            if small_bt:
                bt_type = 'SMALL'
            elif large_bt:
                bt_type = 'LARGE'
            sys.stderr.write(">>>   Found bowtie2 " + bt_type + " index. Great job !   \n")
        else:
            sys.stderr.write(">>>   Can't find .1.bt2 or .1.bt2l   \n")
            test_existance((args.bowtie_index + ".1.bt2"))
    if args.gtf_index is "":
        use_gtf = ""
    else:
        use_gtf = " -G " + args.gtf_index
    if args.transcriptome_index is "":
        use_transcriptome = ""
    else:
        use_transcriptome = " --transcriptome-index " + args.transcriptome_index
    if args.extra_commands is '':
        extra_args = ''.join(args.extra_commands)
    else:
        extra_args = ''.join(args.extra_commands)
        extra_args = " " + extra_args.strip()
    tophat_command = "/n/local/bin/" + aligner + " -p " + args.ncores + use_transcriptome + use_gtf + " -o " +  tophat_dir + "/" + sample + " --no-coverage-search" + extra_args + " " + use_bowtie_one + args.bowtie_index
    fastq_command = " " + fastq_files
    if not os.path.exists(bam_dir):
        os.makedirs(bam_dir)
    if not os.path.exists(tophat_dir):
        os.makedirs(tophat_dir)
    full_command = "'" + tophat_command + fastq_command + "'"
    return(full_command)


def samtools_commander_tophat(sample, args):
    sample_name = sample
    bam_dir = args.output_dir + "/bam"
    tophat_dir = args.output_dir + "/tophat_aln"
    ## add in a symlinking command
    symlink_command = "ln -s " + tophat_dir + "/" + sample + "/accepted_hits.bam" + " " + bam_dir + "/" + sample + ".bam; "
    full_command = "'" + symlink_command + "samtools index "+ bam_dir + "/" + sample + ".bam" + "'"
    return(full_command)


def align_job_command(sample_fastq, job_id, command, args):
    sys.stderr.write(">>>   Preparing some SGE cluster kung fu to run the alignments:   \n")
    sys.stderr.write(space_string)
    processors =  args.ncores
    sge_out = args.output_dir + "/SGE_out"
    job_name = "aln_" + sample_fastq
    ## skip indexing pause if necessary
    if job_id is None:
        if args.email is None:
            job_string = """-m ae -N %s -pe %s %s -q %s -l mem_free=1G,h_vmem=10G -j y -o %s -V -cwd -terse -b y %s""" % (job_name, args.pe, processors, args.queue, sge_out, command)
        else:
            job_string = """-M %s -m ae -N %s -pe %s %s -q %s -l mem_free=1G,h_vmem=10G -j y -o %s -V -cwd -terse -b y %s""" % (args.email, job_name, args.pe, processors, args.queue, sge_out, command)
    else:
        if args.email is None:
            job_string = """-m ae -N %s -pe %s %s -q %s -l mem_free=1G,h_vmem=10G -j y -o %s -hold_jid %s -V -cwd -terse -b y %s""" % (job_name, args.pe, processors, args.queue, sge_out, job_id, command)
        else:
            job_string = """-M %s -m ae -N %s -pe %s %s -q %s -l mem_free=1G,h_vmem=10G -j y -o %s -hold_jid %s -V -cwd -terse -b y %s""" % (args.email, job_name, args.pe, processors, args.queue, sge_out, job_id, command)
    server_command = """qsub %s""" % (job_string)
    sys.stderr.write(space_string)
    sys.stderr.write(">>>   Command:   \n")
    sys.stderr.write(server_command + "\n")
    sys.stderr.write(space_string)
    align_job = subprocess.Popen(server_command,shell=True,stdout=subprocess.PIPE)
    align_job_id = align_job.stdout.read().rstrip()
    sys.stderr.write(space_string)
    sys.stderr.write(">>>   READY! Alignment job parameters confirmed and pushed to server.   \n")
    sys.stderr.write(space_string)
    return align_job_id


def index_job_command(sample_name, job_id, command, args):
    processors = "1"
    sge_out = args.output_dir + "/SGE_out"
    sys.stderr.write(">>>   BEEP BEEP - Setting up post-alignment indexing:   \n")
    sys.stderr.write(space_string)
    job_name = "idx_" + sample_name
    if args.email is None:
        job_string = """-m ae -N %s -pe %s %s -q %s -l mem_free=1G -j y -o %s -hold_jid %s -V -cwd -terse -b y %s""" % (job_name, args.pe, processors, args.queue, sge_out, job_id, command)
    else:
        job_string = """-M %s -m ae -N %s -pe %s %s -q %s -l mem_free=1G -j y -o %s -hold_jid %s -V -cwd -terse -b y %s""" % (args.email, job_name, args.pe, processors, args.queue, sge_out, job_id, command)
    server_command = """qsub %s""" % (job_string)
    sys.stderr.write(space_string)
    sys.stderr.write(">>>   Command:   \n")
    sys.stderr.write(server_command + "\n")
    sys.stderr.write(space_string)
    index_job = subprocess.Popen(server_command,shell=True,stdout=subprocess.PIPE)
    index_job_id = index_job.stdout.read().rstrip()
    final_out = ("Sample " + sample_name + " prepared for alignment as tophat job id: " + job_id + ". Queued prior to index job id: " + index_job_id + ". Pushed to server at " + str(datetime.datetime.now()))
    sys.stderr.write(space_string)
    sys.stderr.write(">>> " + final_out  + " \n")
    sys.stderr.write(space_string)
    sys.stderr.write(">>>   Many congratulations. You did it! Give yourself a pat on the back.   \n")    
    sys.stderr.write(space_string)
    sys.stderr.write(">>>   Results in:" + args.output_dir + " bam folder!   \n")
    sys.stderr.write(space_string)    
    sys.stderr.write(">>>   Check the SGE_out file post-run for logs and potential internal errors.   \n")
    sys.stderr.write(space_string)
    return index_job_id


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument('--directories', action='store', nargs='*', dest='molng_path', required=True,
                        help="Path to the location of a flowcells containting the fastqs and the Sample_Report.csv file describing the samples. Multiple directories may be passed (separated by space)")
    
    parser.add_argument('--transcriptome-index', action='store', dest='transcriptome_index', default='',
                        help='Path to pre-built transcriptome index or location for building an index')

    parser.add_argument('--gtf', action='store', dest='gtf_index', default='',
                        help='Path to GTF file for building a transcriptome index')
    
    parser.add_argument('--bowtie', action='store', dest='bowtie_index', required=True,
                        help='Path to a bowtie2 index *.bt2 -- if specifying a bowtie1 index please include the --use-bowtie1 command')
    
    parser.add_argument('--destdir', action='store', dest='output_dir', required=True,
                        help='Output destination. The directory destination for the alignment output')
    
    parser.add_argument('--use-bowtie1', action='store_true', dest='bowtie_one', default=False,
                        help='Include this if you are using a bowtie one index ')
    
    parser.add_argument('--paired', action='store_true', dest='paired_end', default=False,
                        help='Include this if the run is paired end otherwise paired end files will be ignored')

    parser.add_argument('--email', action='store', dest='email',
                        help='E-mail address for getting status updates during the run')

    parser.add_argument('--cpu', action='store', dest='ncores', default='4',
                        help='Number of threads to use per alignment job (Default: 4)')
    
    parser.add_argument('--queue', action='store', dest='queue', default='all.q',
                        help='SGE queue to use (Default: all.q)')
    
    parser.add_argument('--pe', action='store', dest='pe', default='by_node',
                        help='SGE parallel environment setting (Default: by_node)')

    parser.add_argument('--extra', action='store', nargs='*', dest='extra_commands', default='',
                        help='Provide additional parameters for the aligner. PUT QUOTES AROUND THEM!')

    args = parser.parse_args()
    molng_path = args.molng_path
    if type(molng_path) is str:
        molng_path = [molng_path]
    ##molng_path = molng_path + '/'
    molng_path = [s + '/' for s in molng_path]
    sample_reports  = analyze_molng(molng_path,args)
    ##index_check = check_index(args)
    job0 = transcript_test_index(args)    
    for i in sample_reports:
        com1 = tophat_commander(sample_reports,i,args)
        com2 = samtools_commander_tophat(i,args)
        job1 = align_job_command(i,job0,com1,args)
        job2 = index_job_command(i,job1,com2,args)
