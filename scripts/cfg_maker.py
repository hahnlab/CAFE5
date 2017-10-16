import argparse

# def cfg_writer(cafexp_folder_opt, input_files_opt, options_opt, parameters_opt, output_opt, output_path):
#     """Write .cfg files from list of CAFExp option-value tuples"""

class ArgStash:

    def __init__(self, my_args):
        # dealing with multiple simulation/inference runs (either through -n, or through -i)
        print my_args
        self.n_runs = int()
        self.instructions_list = list()

        try:
            if my_args.instr_path and my_args.number_runs != "1":
                exit('You can only run many simulations with the same parameters, or specify a file with instructions. Exiting...\n')
            elif my_args.instr_path:
                self.instr_path = my_args.instr_path
                self.read_instruction_file()
            else:
                self.n_runs = int(my_args.number_runs)
        except:
            exit('Either you did no specify the path to the instruction file (-i), or you did not specify the number of simulations (-n). Exiting...')

        # where to write .cfg files 
        self.output_path = './'
        if my_args.output_path:
            self.output_path = my_args.output_path

        # CAFExp options
        self.cafexp_path = my_args.cafexp_path
        self.input_args = self.tokenize_string(my_args.input_args)
        self.input_values = self.tokenize_string(my_args.input_values)
        self.options_args = self.tokenize_string(my_args.options_args)
        self.options_values = self.tokenize_string(my_args.options_values)
        self.parameters_args = self.tokenize_string(my_args.parameters_args)
        self.parameters_values = self.tokenize_string(my_args.parameters_values)
        self.output_args, self.output_values = list(), list()
        
        if my_args.output_args:
            self.output_args = self.tokenize_string(my_args.output_args)
            self.output_values = self.tokenize_string(my_args.output_values)

        # finally writing stuff
        self.cfg_writer()

    def read_instruction_file(self):
        with open(self.instr_path, 'r') as instr_file:
            for idx, line in enumerate(instr_file):
                line = line.rstrip()

                if not line.startswith('#'):
                    tokens = line.split()
                    self.instructions_list.append(tokens)

        self.n_runs = len(self.instructions_list)
    
    def cfg_writer(self):
        for r in xrange(self.n_runs):
            with open(self.output_path + 'config_file_' + str(r) + '.cfg', 'w') as config_file:
                config_file.write('[CAFExp folder]\ncafexp = ' + self.cafexp_path)
                
                config_file.write('\n\n[input files]\ngene families = ')
                input_dict = self.dictionize_lists(self.input_args, self.input_values)            
                if 'i' in input_dict: config_file.write(input_dict['i'])
                config_file.write('\ntree = ')
                if 't' in input_dict: config_file.write(input_dict['t'])
                config_file.write('\nlambda tree = ')
                if 'y' in input_dict: config_file.write(input_dict['y'])
                config_file.write('\nroot distribution = ')
                if 'f' in input_dict: config_file.write(input_dict['f'])

                config_file.write('\n\n[options]\nsimulate = ')
                options_dict = self.dictionize_lists(self.options_args, self.options_values)
                if 's' in options_dict: config_file.write('True')
                else: config_file.write('False')
                config_file.write('\nsimulationsN = ')
                if 'n' in options_dict: config_file.write(options_dict['n'])
                config_file.write('\ngammacatN = ')
                if 'k' in options_dict: config_file.write(options_dict['k'])
                config_file.write('\ngammacatMultiplier = ')
                config_file.write('\ngammacat = ')

                par_dict = self.dictionize_lists(self.parameters_args, self.parameters_values)
                config_file.write('\n\n[parameters]\nlambda = ')
                if 'l' in par_dict and not self.instructions_list: config_file.write(par_dict['l'])
                else: config_file.write(self.instructions_list[r][1])
                config_file.write('\nalpha = ')
                if 'a' in par_dict and not self.instructions_list: config_file.write(par_dict['a'])
                elif self.instructions_list and self.instructions_list[r][2] != 'N/A': config_write(self.instructions_list[r][2])
                config_file.write('\npoisson = ')
                if 'p' in par_dict: config_file.write(par_dict['p'])

                output_dict = self.dictionize_lists(self.output_args, self.output_values)
                config_file.write('\n\n[output]\noutput folder = ')

                config_file.write('\noutput suffix = ')
                # for multiple simulations with same parameter values
                if self.n_runs > 1 and not self.instructions_list:
                    config_file.write(str(r+1))
                # for multiple simulations with instruction file
                elif self.instructions_list:
                    config_file.write(self.instructions_list[r][3])
                # for inference
                elif self.n_runs == 1 and not self.instructions_list and 's' in output_dict:
                    config_file.write(output_dict['s'])    
            
    @staticmethod
    def tokenize_string(a_string):
        if a_string == ',':
            return []

        return a_string.split(',')

    @staticmethod
    def dictionize_lists(list1, list2):
        if len(list1) != len(list2): exit("The list of arguments was longer or shorter than the list of values")
        
        return dict((i, list2[idx]) for idx, i in enumerate(list1))
    
if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--cafexp-bin-path", action="store", dest="cafexp_path", default='./', type=str, help="path to CAFExp bin", required=True)
    parser.add_argument("-ia", "--input-args", action="store", dest="input_args", default='', type=str, help="arguments for input options, separated by comma", required=True)
    parser.add_argument("-iv", "--input-values", action="store", dest="input_values", default='', type=str, help="values for input section, separated by comma", required=True)
    parser.add_argument("-oa", "--options-args", action="store", dest="options_args", default='', type=str, help="arguments for options section, separated by comma", required=False)
    parser.add_argument("-ov", "--options-values", action="store", dest="options_values", default='', type=str, help="values for options section, separated by comma", required=False)
    parser.add_argument("-pa", "--parameters-args", action="store", dest="parameters_args", default='', type=str, help="arguments for parameters section, separated by comma", required=False)
    parser.add_argument("-pv", "--parameters-values", action="store", dest="parameters_values", default='', type=str, help="values for parameters section, separated by comma", required=False)
    parser.add_argument("-oua", "--output-args", action="store", dest="output_args", default='', type=str, help="arguments for output section, separated by comma", required=False)
    parser.add_argument("-ouv", "--output-values", action="store", dest="output_values", default='', type=str, help="values for output section, separated by comma", required=False)
    parser.add_argument("-o", "--output-path", action="store", dest="output_path", default='./', type=str, help="path to write cfg files")
    parser.add_argument("-n", "--n-runs", action="store", dest="number_runs", default="1", type=str, help="number of runs (will go into output suffix, one per cfg file, 1:number_runs)")
    parser.add_argument("-i", "--intructions-path", action="store", dest="instr_path", default='', type=str, help="path to file with simulation instructions", required=False)
    
    args = parser.parse_args()

    my_arg_stash = ArgStash(args)
