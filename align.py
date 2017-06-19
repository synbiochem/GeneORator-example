'''
synbiochem (c) University of Manchester 2016

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
import subprocess
import sys
import tempfile

from Bio import SeqIO, Seq, SeqRecord
import pysam


def align(templ_filename, seqs_filename):
    '''Aligns sequences in barcoded bins.'''
    # Index template:
    subprocess.call(['bwa', 'index', templ_filename])

    # Read sequences:
    with open(seqs_filename, 'rU') as fle:
        seqs = {record.id: str(record.seq)
                for record in SeqIO.parse(fle, 'fasta')}

    # Align and strip indels from file:
    strip_id_seqs = _strip_indels(_mem(seqs, templ_filename), templ_filename)

    # Align and sort strip indels:
    _sort(_mem(strip_id_seqs, templ_filename), 'align.sam')


def _mem(seqs, templ_filename, readtype='ont2d'):
    '''Runs BWA MEM.'''
    out_file = tempfile.NamedTemporaryFile(delete=False)
    seq_file = tempfile.NamedTemporaryFile(delete=False)

    records = [SeqRecord.SeqRecord(Seq.Seq(seq), seq_id, '', '')
               for seq_id, seq in seqs.iteritems()]

    SeqIO.write(records, seq_file.name, 'fasta')

    with open(out_file.name, 'w') as out:
        subprocess.call(['bwa', 'mem', '-x', readtype,
                         templ_filename, seq_file.name], stdout=out)

    return out_file.name


def _sort(in_filename, out_filename):
    '''Custom sorts SAM file.'''
    sam_file = pysam.Samfile(in_filename, 'r')
    out_file = pysam.AlignmentFile(out_filename, 'wh',
                                   template=sam_file,
                                   header=sam_file.header)

    for read in sorted([read for read in sam_file],
                       key=lambda x: (-x.query_length,
                                      x.reference_start)):
        out_file.write(read)

    out_file.close()

    return out_filename


def _strip_indels(in_filename, templ_filename):
    '''Generator to convert sam files into Biopython SeqRecords.'''
    seqs = {}
    sam_file = pysam.Samfile(in_filename, 'r')

    with open(templ_filename, 'rU') as fle:
        templ_seq = [str(record.seq)
                     for record in SeqIO.parse(fle, 'fasta')][0]

    for read in sam_file:
        # Perform mapping of nucl indices to remove spurious indels:
        seq = ''.join([read.seq[pair[0]]
                       if pair[0]
                       else templ_seq[pair[1]]
                       for pair in read.aligned_pairs
                       if pair[1] is not None])

        if len(seq):
            seqs[read.qname] = seq

    return seqs


def main(args):
    '''main method.'''
    align(args[0], args[1])


if __name__ == '__main__':
    main(sys.argv[1:])
