'''
This is a quick script that runs as follows:
python3 simulate_low_coverage_genome.py KY075937.fasta

It takes in a reference genome file and simulates reference genomes with Ns
inserted to reach the coverage breadths specified within the script
(currently .9, .75, .5, .25)
'''

from Bio import SeqIO
from collections import Counter
from random import randint
import sys

input_filename = sys.argv[1]
# simulated_coverage = float(sys.argv[2])

records = []
for r in SeqIO.parse(input_filename, "fasta"):
    records.append(r)

this_sequence = records[0].seq
len_of_seq = len(this_sequence)
nt_counts = Counter(this_sequence)
existing_coverage = 1 - nt_counts['N']/len_of_seq


for simulated_coverage in [.90, .75, .5, .25]:
    print(simulated_coverage)
    delta_to_simulated_coverage = existing_coverage - simulated_coverage
    ns_to_create_sim_cov = round(delta_to_simulated_coverage*len_of_seq)
    print(ns_to_create_sim_cov)
    start_pos = randint(0, (len_of_seq-ns_to_create_sim_cov-nt_counts['N']))

    print('desired simulated cov:', simulated_coverage)
    print('existing cov:', existing_coverage)

    # generate new sequence with appropriate number of Ns
    new_sequence = list(str(this_sequence))
    i = 0
    position = 0
    while i < ns_to_create_sim_cov:
        position += 1
        if new_sequence[start_pos + position] == 'N':
            continue
        else:
            new_sequence[start_pos + position] = 'N'
            i += 1

    final_coverage = 1 - Counter(new_sequence)['N']/len_of_seq
    print("final coverage:", final_coverage)

    output_filename = '.'.join(input_filename.split('.')[:-1]) + '_' + str(int(simulated_coverage*100)) + '.fasta'
    with open(output_filename, 'w') as the_file:
        the_file.write('>'+str(records[0].id))
        the_file.write(''.join(new_sequence))
