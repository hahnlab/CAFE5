import sys

# sys.argv[1]: cols to keep
# sys.argv[2]: input data path
# sys.argv[3]: output data path

tips = sys.argv[1].split(',')
spp_names = list()
print_header = False
with open(sys.argv[3], 'w') as output_file:
    with open(sys.argv[2], 'r') as input_data:
        for line in input_data:
            if line.startswith('#'):
                sp = line.lstrip('#').rstrip()
                spp_names.append(sp)
            else:
                if not print_header:
                    output_file.write('\t'.join([n for i,n in enumerate(spp_names) if str(i+1) in tips]) + '\n')
                    print_header = True
                tokens = line.rstrip().split('\t')
                output_file.write('\t'.join([t for j,t in enumerate(tokens) if str(j+1) in tips]) + '\n')

