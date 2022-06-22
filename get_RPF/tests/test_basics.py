from get_RPF.kmer import read_to_kmers

def test_read_to_dimer():
    assert read_to_kmers("ACGTA", 2) == ['AC', 'CG', 'GT', 'TA']

def test_read_to_trimer():
    assert read_to_kmers("ACGTA", 3) == ['ACG', 'CGT', 'GTA']

def test_read_to_quadmer():
    assert read_to_kmers('ACGTA', 4) == ['ACGT', 'CGTA']