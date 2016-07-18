import sys
import subprocess as sp

# with open(sys.argv[1], 'r') as starOutput:
#     for line in starOutput:
#         print line
#         jxn_data = line.strip().split()
#         chr_list = ["chr1", "chr2", "chr3", "chr4", "chr5"]
#         myvals = {}
#         for chr in chr_list:
#             if jxn_data[0] == chr_list[0] or jxn_data[3] == chr_list[0]:
#                 if key in myvals:
#                     my_dict[key] += 1
#                 else:
#                     my_dict[key] = 1
#         bedpe = '\t'.join([chrL, str(int(posL) - 1), posL, chrR, str(int(posR) - 1), posR, '-'.join([geneL, geneR]), '0', strandL, strandR])
#         print bedpe


def sortbedpe(input_path, output_path):
    outbedpe_fh = open(output_path, "w")
    sp.call(["sort", "-k1,1", "-k2,2n", "-k4,4", "-k5,5n", input_path], stdout=outbedpe_fh)
    outbedpe_fh.close()

sortbedpe(sys.argv[1], "test.bedpe")
