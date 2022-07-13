import subprocess
import pysam

def bowtie_align_read_sample(fastq_path, index_path, number_of_reads=10000):
    '''
    align a sample of reads (number of reads size) to a provided reference and return the alignment output
    ''' 

    alignment_output = subprocess.run(['bowtie', '-v', '3', '-p', '8', f'{index_path}', '-q', f'{fastq_path}', '-u', f'{number_of_reads}'], stderr=subprocess.PIPE, stdout=subprocess.PIPE)
    alignment_output_stderr = alignment_output.stderr
    alignment_output_stdout = alignment_output.stdout
    return alignment_output_stderr.decode('utf-8'), alignment_output_stdout.decode('utf-8')


def bowtie_alignment_output_stderr_to_report_dict(alignment_output):
    '''
    take the string outputted by bowtie with information on read alignments and return a dictionary with the parsed information  
    the function starts by removing the hashes at the start of each line and creating a list of each line in the output 
    then some ugly clean up is done to produce a dictionary with keys ['reads_processed', 'reads_with_at_least_one_reported_alignment', 'reads_that_failed_to_align']
    '''
    cleaned_output_list = [i.strip('# ') for  i in alignment_output.split('\n')]

    report_dict = {}
    for line in cleaned_output_list:
        line_elements_list = line.split(': ')
        try:
            report_dict[line_elements_list[0].replace(' ', '_')] = int(line_elements_list[1].split(' ')[0])
        except: 
            continue

    return report_dict

    
def bowtie_alignment_output_stdout_parsing(alignment_output):
    '''
    Take the default bowtie alignment output and return a dictionary with alignment details
    '''
    cleaned_alignment_output = [i.split('\t') for i in alignment_output.split('\n')]
    print(cleaned_alignment_output)


if __name__ == "__main__":
    not_clipped = "/home/jack/projects/riboseq_database/data/Park2017_Homo_sapiens_GSE97384_SRP103009/fastq/SRR5413155.fastq"
    clipped = "/home/jack/projects/riboseq_database/data/Park2017_Homo_sapiens_GSE97384_SRP103009/fastq/park17_clipped_fq/SRR5413155.fastq.clipped_fastq"
    index_path = "/mnt/data/indices/bowtie/genome/hg38/hg38"
    # print(bowtie_alignment_output_stderr_to_report_dict(bowtie_align_read_sample(not_clipped, index_path)[1]))
    print(bowtie_alignment_output_stdout_parsing(bowtie_align_read_sample(not_clipped, index_path)[1]))
    # print()
    # print(bowtie_align_read_sample(clipped, index_path))


'''
To DO:

- Design brute force method where reads are modified one base at a time for a certain numbe of bases until read mapping meets certain % threshold 
- Function to check if threshold is met by stderr 
- Explore the use of STAR as the aligned where read modificaitons are informed by soft clipping. 
- When the RPF can be identified we can identify UMI/Barcodes from patterns in reads
- Explore non-templated additions 


'''