#!/usr/bin/env python
"""
liftover_intropolis.py

Lifts intropolis over to hg19 and writes output including both hg38 coordinates
and hg19 coordinates.

Requires

http://hgdownload.cse.ucsc.edu/goldenpath/hg38/
    liftOver/hg38ToHg19.over.chain.gz

liftOver executable available from
    https://genome-store.ucsc.edu/products/

intropolis.v2.hg38.tsv.gz .

Writes to stdout. We ran

pypy liftover_intropolis.py
    --liftover /path/to/liftOver
    --chain /path/to/hg38ToHg19.over.chain
    --intropolis /path/to/intropolis.v2.hg38.tsv.gz
    | gzip >intropolis.v2.hg38_with_liftover_to_hg19.tsv.gz

Tab-separated output fields
1. hg38 chrom
2. hg38 start (1-based, inclusive)
3. hg38 end (1-based, inclusive)
4. hg38 strand
5. left motif (e.g., GT)
6. right motif (e.g., AG)
7. comma-separated list of indexes of samples from
    intropolis in which junction was found
8. comma-separated list of numbers of reads in corresponding samples from
    field 7 overlapping junction
9. hg19 chrom or NA if liftover unsuccessful
10. hg19 start or NA
11. hg19 end or NA
12. hg19 strand or NA
"""
import tempfile
import gzip
import shutil
import atexit
import subprocess
import os

if __name__ == '__main__':
    import argparse
    # Print file's docstring if -h is invoked
    parser = argparse.ArgumentParser(description=__doc__,
                formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--liftover', type=str, required=True,
            help=('path to liftOver executable available from '
                  'https://genome-store.ucsc.edu/products/')
        )
    parser.add_argument('--chain', type=str, required=True,
            help=('path to unzipped liftover chain; this should be '
                  'hg38ToHg19.over.chain')
        )
    parser.add_argument('--intropolis', type=str, required=True,
            help='path to intropolis.v2.hg38.tsv.gz'
        )
    parser.add_argument('--temp-dir', type=str, required=False,
            default=None,
            help='where to store temporary files; defaults to TMPDIR'
        )
    args = parser.parse_args()
    temp_dir = tempfile.mkdtemp(dir=args.temp_dir)
    to_liftover = os.path.join(temp_dir, 'to_liftover.bed')
    temp_hg19 = os.path.join(temp_dir, 'hg19.bed')
    temp_hg38 = os.path.join(temp_dir, 'hg38.bed')
    with open(temp_hg38, 'w') as hg38_stream, gzip.open(
            args.intropolis
        ) as input_stream:
        for i, line in enumerate(input_stream):
            tokens = line.strip().split('\t')
            chrom, start, end, strand = (
                    tokens[0], str(int(tokens[1]) - 1),
                    tokens[2], tokens[3],
                ) # zero-based, half-open coordinates for BED
            # Tack original junction onto junction name
            junction_name = ';'.join([str(i), chrom, start, end, strand])
            print >>hg38_stream, '{}\t{}\t{}\tinfo_{}\t1\t{}'.format(
                    chrom, start, end, junction_name, strand
                )
    # Convert junctions from hg38 to hg19
    liftover_process = subprocess.call(' '.join([
                                            args.liftover,
                                            '-ends=2',
                                            '-minMatch=1.0',
                                            temp_hg38,
                                            args.chain,
                                            temp_hg19,
                                            os.path.join(
                                                    temp_dir,
                                                    'unmapped.bed'
                                                )
                                        ]),
                                        shell=True,
                                        executable='/bin/bash'
                                    )
    to_sort = os.path.join(temp_dir, 'intropolis_and_liftover.tsv.gz')
    with gzip.open(to_sort, 'w') as both_stream:
        with open(temp_hg19) as hg19_stream:
            for line in hg19_stream:
                chrom, start, end, name, score, strand = line.strip().split(
                                                                        '\t'
                                                                    )[:6]
                (_, hg38_chrom, hg38_start,
                        hg38_end, hg38_strand) = name.split(';')
                hg38_start, start = int(hg38_start), int(start)
                print >>both_stream, '\t'.join(
                                [hg38_chrom, str(hg38_start + 1), hg38_end,
                                    hg38_strand, chrom, str(start + 1),
                                    end, strand, 'FAKE']
                            )
        with gzip.open(args.intropolis) as intropolis_stream:
            for line in intropolis_stream:
                print >>both_stream, line,
    sorted_together = os.path.join(temp_dir, 'sorted_together.tsv.gz')
    subprocess.check_call(
            'gzip -cd {} | sort -k1,1 -k2,2n -k3,3n | gzip >{}'.format(
                    to_sort, sorted_together
                ), shell=True, bufsize=-1
        )
    with gzip.open(sorted_together) as sorted_stream:
        junction_1_tokens = sorted_stream.readline().strip().split('\t')
        junction_2_tokens = sorted_stream.readline().strip().split('\t')
        while True:
            if junction_1_tokens[:4] == junction_2_tokens[:4]:
                # Liftover available
                if len(junction_1_tokens) > len(junction_2_tokens):
                    hg19_tokens = junction_1_tokens
                    hg38_tokens = junction_2_tokens
                else:
                    hg19_tokens = junction_2_tokens
                    hg38_tokens = junction_1_tokens
                print '\t'.join(hg38_tokens + hg19_tokens[4:8])
                junction_1_tokens = sorted_stream.readline().strip()
                if not junction_1_tokens:
                    # End of file; nothing to print
                    break
                junction_1_tokens = junction_1_tokens.split('\t')
                junction_2_tokens = sorted_stream.readline().strip()
                if not junction_2_tokens:
                    # End of file; print junction 1 tokens and sign out
                    print '\t'.join(junction_1_tokens + ['NA'] * 4)
                    break
                junction_2_tokens = junction_2_tokens.split('\t')
            else:
                '''Liftover not available for junction 1, but have to check
                junction 2 against next junction.'''
                print '\t'.join(junction_1_tokens + ['NA'] * 4)
                junction_1_tokens = junction_2_tokens
                junction_2_tokens = sorted_stream.readline().strip()
                if not junction_2_tokens:
                    # End of file; print new junction 1 tokens and sign out
                    print '\t'.join(junction_1_tokens + ['NA'] * 4)
                    break
                junction_2_tokens = junction_2_tokens.split('\t')
