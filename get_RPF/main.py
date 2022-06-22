import Bio
import argparse

def get_parser():
    '''
    return the argparse parser for this script 
    '''
    parser = argparse.ArgumentParser()

    parser.add_argument("-f", "--fastq", help="path to the FASTQ file you want to investigate")
    parser.add_argument("-m", "--mode", choices=["kmer", "align"], help="Method to use to identify rpf")
    parser.add_argument("-b", "--bowtie", help="Bowtie index to use for alignment method")

    return parser


if __name__ == "__main__":
    parser = get_parser()
    args = parser.parse_args()

    if args.fastq:
        print(args.fastq)
    print()