import pytest
import os
import tempfile
from mut_stats import calculate_mutation_stats

def create_temp_tsv(content):
    """Helper to create temporary TSV files."""
    fd, path = tempfile.mkstemp(suffix='.tsv')
    with os.fdopen(fd, 'w') as tmp:
        tmp.write(content)
    return path

def test_basic_counts():
    """Test basic mutation counting."""
    tsv_content = """REGION\tPOS\tREF\tALT\tPASS
chr1\t100\tA\tG\tTRUE
chr1\t200\tC\tT\tTRUE
chr1\t300\tA\tT\tTRUE
chr1\t400\tG\tC\tTRUE
chr1\t500\tA\tAT\tTRUE
chr1\t600\tAT\tA\tTRUE
chr1\t700\tA\tN\tFALSE
"""
    tsv_file = create_temp_tsv(tsv_content)
    stats = calculate_mutation_stats(tsv_file)

    assert(stats['total_mutations'] == 6)
    assert(stats['insertions'] == 1)
    assert(stats['deletions'] == 1)
    assert(stats['snps'] == 4)
    assert(stats['transitions'] == 2)
    assert(stats['transversions'] == 2)
    assert(stats['ti_tv_ratio'] == 1.00)
    os.unlink(tsv_file)



def test_no_transversions():
    """Test case with only transitions."""
    tsv_content = """REGION\tPOS\tREF\tALT\tPASS
chr1\t100\tA\tG\tTRUE
chr1\t200\tC\tT\tTRUE
"""
    tsv_file = create_temp_tsv(tsv_content)
    stats = calculate_mutation_stats(tsv_file)
    #output = capture_stdout.getvalue()

    assert(stats["transversions"] == 0)
    assert(stats["ti_tv_ratio"] == float("inf"))
    
    os.unlink(tsv_file)

def test_no_passing_mutations():
    """Test case where no mutations PASS."""
    tsv_content = """REGION\tPOS\tREF\tALT\tPASS
chr1\t100\tA\tG\tFALSE
chr1\t200\tC\tT\tFALSE
"""
    tsv_file = create_temp_tsv(tsv_content)
    stats = calculate_mutation_stats(tsv_file)
    #output = capture_stdout.getvalue()
    assert(stats["total_mutations"] == 0)

    os.unlink(tsv_file)

def test_complex_indels():
    """Test multi-base insertions/deletions."""
    tsv_content = """REGION\tPOS\tREF\tALT\tPASS
chr1\t100\tA\tAC\tTRUE
chr1\t200\tAC\tA\tTRUE
chr1\t300\tACG\tA\tTRUE
chr1\t400\tA\tACGT\tTRUE
"""
    tsv_file = create_temp_tsv(tsv_content)
    stats = calculate_mutation_stats(tsv_file)
    #output = capture_stdout.getvalue()
    
    assert(stats["insertions"] == 2)
    assert(stats["deletions"] == 2)
    assert(stats["snps"] == 0)
    
    os.unlink(tsv_file)

def test_case_insensitivity():
    """Test that case doesn't affect counting."""
    tsv_content = """REGION\tPOS\tREF\tALT\tPASS
chr1\t100\ta\tg\tTRUE
chr1\t200\tc\tt\tTRUE
chr1\t300\tA\tt\tTRUE
"""
    tsv_file = create_temp_tsv(tsv_content)
    stats = calculate_mutation_stats(tsv_file)

    assert(stats["transitions"] == 2)
    assert(stats["transversions"] == 1)
    
    os.unlink(tsv_file)