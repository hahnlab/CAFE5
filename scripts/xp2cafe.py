import sys

# sys.argv[1]: cols to keep
# sys.argv[2]: input data path
# sys.argv[3]: output data path

tips = sys.argv[1].split(',')
spp_names = list()
print_header = False
with open(sys.argv[3]) as output_file:
    with open(sys.argv[2]) as input_data:
        for line in input_data:
            if line.startswith('#'):
                sp = line.lstrip('#').rstrip()
                spp_names.append(sp)
            else:
                if not print_header:
                    output_file.write('\t'.join(spp_names) + '\n')
                    print_header = True
                tokens = line.rstrip().split('\t')
                output_file.write('\t'.join([t for t in tokens if str(t+1) in tips]) + '\n')

