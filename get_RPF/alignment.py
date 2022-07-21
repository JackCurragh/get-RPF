import subprocess
import os


def is_gz_file(filepath):
    with open(filepath, 'rb') as test_f:
        return test_f.read(2) == b'\x1f\x8b'

def check_and_unzip_fq(fastq_path):
    if is_gz_file(fastq_path):
        print(fastq_path)

        subprocess.run(['gzip', '-k', '-f', '-d', f'{fastq_path}'])
        print(fastq_path.split('.'))
        fastq_path = '.'.join(fastq_path.split('.')[:-1])
    return fastq_path


def make_tmp_output_dir(fastq_path):
    '''
    create the tmp output dir if required
    '''
    tmp_output_dir = '/'.join(fastq_path.split('/')[:-1]) 
    if 'tmp' not in fastq_path:
        if not os.path.exists(tmp_output_dir + '/tmp'):
            os.makedirs(tmp_output_dir + '/tmp')
        tmp_output_dir += '/tmp/'
    return tmp_output_dir


def trim_reads(fastq_path, five_prime=0, three_prime=0, number_of_reads=10000):
    '''
    remove fixed number of bases from 5' and/or 3' of each read 
    '''

    tmp_output_dir = make_tmp_output_dir(fastq_path)

    filename = fastq_path.split('/')[-1]
    if 'less_' in filename and '_fiveprime' in filename:
        filename_list = filename.split('_')
        pervious_fiveprime_count = int(filename_list[filename_list.index('fiveprime') - 1])
        pervious_threeprime_count = int(filename_list[filename_list.index('threeprime') - 1])


        # outpath = f"{tmp_output_dir}less_{pervious_fiveprime_count + five_prime}_fiveprime_less_{pervious_threeprime_count + three_prime}_threeprime_{filename_list[-1]}"

        outpath = f"{tmp_output_dir}/less_{five_prime}_fiveprime_less_{three_prime}_threeprime_{filename_list[-1]}"
    else:
        filename_list = filename.split('_')
        outpath = f"{tmp_output_dir}less_{five_prime}_fiveprime_less_{three_prime}_threeprime_{filename_list[-1]}"


    if five_prime > 0 and three_prime > 0:
        head_output = subprocess.Popen(['head', '-n', f'{number_of_reads * 4}', f'{fastq_path}'], stdout=subprocess.PIPE)
        trimming_output = subprocess.run(['cutadapt', '-u', f'{five_prime}', '-o', f'{outpath}', '-'], stdin=head_output.stdout, stderr=subprocess.PIPE, stdout=subprocess.PIPE)

        head_output.wait()
        head_output2 = subprocess.Popen(['head', '-n', f'{number_of_reads * 4}', f'{outpath}'], stdout=subprocess.PIPE)
        trimming_output = subprocess.run(['cutadapt', '-u', f'-{three_prime}', '-o', f'{outpath}', '-'], stdin=head_output2.stdout, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
        head_output2.wait()

    elif five_prime > 0:
        head_output = subprocess.Popen(['head', '-n', f'{number_of_reads * 4}', f'{fastq_path}'], stdout=subprocess.PIPE)
        trimming_output = subprocess.run(['cutadapt', '-u', f'{five_prime}', '-o', f'{outpath}', '-'], stdin=head_output.stdout, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
        head_output.wait()


    elif three_prime > 0:
        head_output = subprocess.Popen(['head', '-n', f'{number_of_reads * 4}', f'{fastq_path}'], stdout=subprocess.PIPE)
        trimming_output = subprocess.run(['cutadapt', '-u', f'-{three_prime}', '-o', f'{outpath}', '-'], stdin=head_output.stdout, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
        head_output.wait()
    
    # print(trimming_output.stdout.decode())
    # print(trimming_output.stderr.decode())

    return outpath


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
    bowtie_column_names comes from the specificaiton http://bowtie-bio.sourceforge.net/manual.shtml#default-bowtie-output

    '''
    bowtie_output_column_names = [ 
        'read_name',
        'strand',
        'chr',
        'alignment_start',
        'sequence',
        'phred',
        'number_of_other_alignment_sites',
        'mismatch_descriptors'
    ]
    cleaned_alignment_output = [i.split('\t') for i in alignment_output.split('\n')]

    alignments_dict = {}
    for read in cleaned_alignment_output:
        alignments_dict[read[0]] = {i[0]:i[1] for i in zip(bowtie_output_column_names, read)}

    return alignments_dict


def define_searches(max_5=20, max_3=20, increment=1):
    '''
    produce a list of tuples where tup[0] is the num bases to remove 5' and tup[1] is the num 3'

    increment determines how many tupes are produced. 
    '''
    output_list = []
    for i in range(0, max_5 + 1):
        for j in range(max_3 + 1):
            output_list.append((i,j))

    sorted_output_list = sorted(output_list, key = lambda x:x[0] + x[1])

    return sorted_output_list


def find_rpf(fastq_path, index_path, number_of_reads=10000):
    '''
    Find the position for the rpfs in the given fastq path
    '''
    fastq_path = check_and_unzip_fq(fastq_path)
    results_dict = {}
    for five_prime, three_prime in define_searches():

        if (five_prime, three_prime) == (0, 0):
            alignment_stderr, alignment_stdout = bowtie_align_read_sample(fastq_path, index_path, number_of_reads=number_of_reads)
            report_dict = bowtie_alignment_output_stderr_to_report_dict(alignment_stderr)
            iteration_proportion_aligned = report_dict['reads_with_at_least_one_reported_alignment']/report_dict['reads_processed']
            if iteration_proportion_aligned > 0.9:
                return (five_prime, three_prime)
            continue

        new_fastq_path = trim_reads(fastq_path, five_prime=five_prime, three_prime=three_prime, number_of_reads=number_of_reads)
        alignment_stderr, alignment_stdout = bowtie_align_read_sample(new_fastq_path, index_path, number_of_reads=10000)
        report_dict = bowtie_alignment_output_stderr_to_report_dict(alignment_stderr)
        iteration_proportion_aligned = report_dict['reads_with_at_least_one_reported_alignment']/report_dict['reads_processed']

        if iteration_proportion_aligned > 0.95:
            results_dict[(five_prime, three_prime)] = iteration_proportion_aligned
            return results_dict

        results_dict[(five_prime, three_prime)] = iteration_proportion_aligned
    return results_dict


if __name__ == "__main__":
    not_clipped = "/home/jack/projects/riboseq_database/data/Park2017_Homo_sapiens_GSE97384_SRP103009/fastq/SRR5413155.fastq"
    clipped = "/home/jack/projects/riboseq_database/data/Park2017_Homo_sapiens_GSE97384_SRP103009/fastq/park17_clipped_fq/SRR5413155.fastq.clipped_fastq"
    w_umis = "/home/jack/projects/riboseq_data_processing/data/Chatterji2018_Homo_sapiens_GSE112305_SRP136411/fastq/SRR6893921.fastq.gz_clipped.fastq"
    w_umis_zipped = "/home/jack/projects/riboseq_data_processing/data/Chatterji2018_Homo_sapiens_GSE112305_SRP136411/fastq/SRR6893922.fastq.gz"

    index_path = "/mnt/data/indices/bowtie/genome/hg38/hg38"
    results_dict = find_rpf(w_umis, index_path)

    results_dict_top =  sorted(results_dict.keys(), key = lambda x:x[0] + x[1])
    for i in results_dict_top:
        if results_dict[i] > 0.9:
            fq_path = f"/home/jack/projects/riboseq_data_processing/data/Chatterji2018_Homo_sapiens_GSE112305_SRP136411/fastq/tmp/less_{i[0]}_fiveprime_less_{i[1]}_threeprime_clipped.fastq"
            avg_read_length = subprocess.run(['awk', '{if(NR%4==2) {count++; bases += length} } END{print bases/count}', f'{fq_path}'], stderr=subprocess.PIPE, stdout=subprocess.PIPE)
            print(i, results_dict[i], avg_read_length.stdout)
    print(results_dict_top)

    # print()
    # print(bowtie_align_read_sample(clipped, index_path))

'''
To DO:

- Design brute force method where reads are modified one base at a time for a certain number of bases until read mapping meets certain % threshold 
- Function to check if threshold is met by stderr 
- Explore the use of STAR as the aligned where read modificaitons are informed by soft clipping. 
- When the RPF can be identified we can identify UMI/Barcodes from patterns in reads
- Explore non-templated additions 


'''