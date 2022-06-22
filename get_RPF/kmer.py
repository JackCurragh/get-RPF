from heapq import merge
from Bio import SeqIO
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import subprocess

def read_to_kmers(sequence, k):
    '''
    split read sequence into kmers of length k. Returns list 
    '''
    return [sequence[i:i+k] for i in range(0, len(sequence) - (k - 1))]


def read_to_kmer_range(sequence, kmer_range):
    '''
    produce a dict of k:[kmers] where k are all values within the kmer_range 
    '''
    kmer_dict = {}
    for k in range(kmer_range[0], kmer_range[1]):
        kmer_dict[k] = read_to_kmers(sequence, k)
    return kmer_dict

def count_kmer(fastq_path, kmer):
    '''
    count number of occurrences of kmer in provided file
    '''
    kmer_count_raw = subprocess.check_output(
        f"head -2000000 {fastq_path} | sed -n '2~4p' > test.fastq ; agrep -c1 \"{kmer}\" test.fastq ; rm test.fastq",
        shell=True,
    )
    kmer_count = float(kmer_count_raw.decode('utf-8').strip('\n'))
    return kmer_count


def merge_kmer_dicts(kmer_dict1, kmer_dict2):
    '''
    kmer dicts have the following structure k:[list of kmers] for a range of k values
    this funciton combines the lists avoiding redundancy 
    '''
    if kmer_dict1 == {}:
        return kmer_dict2
    else:
        merged_kmer_dict = {}
        for k in kmer_dict2.keys():
            kmer_list = list(set(kmer_dict1[k] + kmer_dict2[k]))
            merged_kmer_dict[k] = kmer_list
        return merged_kmer_dict


def get_kmer_sample(fastq_path, k_range=(15, 19), reads_to_use=10):
    '''
    run the kmer approach for finding rpf in raw reads 
    '''
    records_used = 0
    kmer_dict = {}
    while records_used < reads_to_use:
        record = next(SeqIO.parse(fastq_path, "fastq"))
        new_kmer_dict = read_to_kmer_range(str(record.seq), k_range)
        kmer_dict = merge_kmer_dicts(kmer_dict, new_kmer_dict)

        records_used += 1
    return kmer_dict


def get_number_of_reads(fastq_path):
    '''
    return the number of reads in the fastq file (num lines / 4)
    '''
    fastq_lines_output = subprocess.check_output(f"wc -l {fastq_path}", shell=True)
    fastq_lines = float(fastq_lines_output.split()[0])
    number_of_reads = fastq_lines/4
    return number_of_reads

def get_kmer_counts(fastq_path, kmer_dict, number_of_reads):
    '''
    Count the number of times each kmer occurs in a read in the fastq file
    '''
    kmer_count_dict = {} 
    for k in kmer_dict:
        print(k)
        kmer_count_dict[k] = {}
        for kmer in kmer_dict[k]:
            kmer_count_dict[k][kmer] = count_kmer(fastq_path, kmer)
            print(kmer, kmer_count_dict[k][kmer], round(kmer_count_dict[k][kmer]/number_of_reads * 100, 2))
        print()


if __name__ == "__main__":
    not_clipped = "/home/jack/projects/riboseq_database/data/Park2017_Homo_sapiens_GSE97384_SRP103009/fastq/SRR5413155.fastq"
    clipped = "/home/jack/projects/riboseq_database/data/Park2017_Homo_sapiens_GSE97384_SRP103009/fastq/park17_clipped_fq/SRR5413155.fastq.clipped_fastq"
    kmer_dict = get_kmer_sample(not_clipped)
    number_of_reads = get_number_of_reads(not_clipped)
    get_kmer_counts(not_clipped, kmer_dict, number_of_reads)
    print(number_of_reads)

