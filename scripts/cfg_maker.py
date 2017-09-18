import argparse

# def cfg_writer(cafexp_folder_opt, input_files_opt, options_opt, parameters_opt, output_opt, output_path):
#     """Write .cfg files from list of CAFExp option-value tuples"""

class ArgStash:

    def __init__(self, my_args):
        self.output_path = './'
        if my_args.output_path:
            self.output_path = my_args.output_path
            
        self.n_runs = int(my_args.number_runs)
        self.cafexp_path = my_args.cafexp_path
        self.input_args = self.tokenize_string(args.input_args)
        self.input_values = self.tokenize_string(args.input_values)
        self.options_args = self.tokenize_string(args.options_args)
        self.options_values = self.tokenize_string(args.options_values)
        self.parameters_args = self.tokenize_string(args.parameters_args)
        self.parameters_values = self.tokenize_string(args.parameters_values)
        self.output_args, self.output_values = list(), list()
        
        if args.output_args:
            self.output_args = self.tokenize_string(args.output_args)
            self.output_values = self.tokenize_string(args.output_values)

        self.cfg_writer()
        
    def cfg_writer(self):
        for r in xrange(self.n_runs):
            with open(self.output_path + 'config_file_' + str(r+1) + '.cfg', 'w') as config_file:
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
                if 'l' in par_dict: config_file.write(par_dict['l'])
                config_file.write('\nalpha = ')
                # if 'a' in parameters_dict: config_file.write(parameters_dict['a'])

                output_dict = self.dictionize_lists(self.output_args, self.output_values)
                config_file.write('\n\n[output]\noutput folder = ')

                config_file.write('\noutput suffix = ')
                if self.n_runs > 1:
                    config_file.write(str(r+1))
                    # config_file.write('\noutput prefix = '+str(r+1))
            
            
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
    parser.add_argument("-c", "--cafexp-bin-path", action="store", dest="cafexp_path", default="./", type=str, help="path to CAFExp bin", required=True)
    parser.add_argument("-ia", "--input-args", action="store", dest="input_args", default=None, type=str, help="arguments for input options, separated by comma", required=True)
    parser.add_argument("-iv", "--input-values", action="store", dest="input_values", default=None, type=str, help="values for input section, separated by comma", required=True)
    parser.add_argument("-oa", "--options-args", action="store", dest="options_args", default=None, type=str, help="arguments for options section, separated by comma", required=True)
    parser.add_argument("-ov", "--options-values", action="store", dest="options_values", default=None, type=str, help="values for options section, separated by comma", required=True)
    parser.add_argument("-pa", "--parameters-args", action="store", dest="parameters_args", default=None, type=str, help="arguments for parameters section, separated by comma", required=False)
    parser.add_argument("-pv", "--parameters-values", action="store", dest="parameters_values", default=None, type=str, help="values for parameters section, separated by comma", required=False)
    parser.add_argument("-oua", "--output-args", action="store", dest="output_args", default=None, type=str, help="arguments for output section, separated by comma", required=False)
    parser.add_argument("-ouv", "--output-values", action="store", dest="output_values", default=None, type=str, help="values for output section, separated by comma", required=False)
    parser.add_argument("-o", "--output-path", action="store", dest="output_path", default="./", type=str, help="path to write cfg files")
    parser.add_argument("-n", "--n-runs", action="store", dest="number_runs", default="1", type=str, help="number of runs (will go into output suffix, one per cfg file, 1:number_runs)")
    
    args = parser.parse_args()

    my_arg_stash = ArgStash(args)
