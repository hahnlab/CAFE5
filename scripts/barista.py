import sys
import os
import ConfigParser
import subprocess

class Barista:
    """Read config file and execute CAFExp"""

    def __init__(self, config_file):
        self.command_line = list()
        self.config = ConfigParser.SafeConfigParser()
        self.config.read(config_file)

        self.binary_path = self.config.get('CAFExp folder', 'cafexp')

        # input files
        self.tree_path = self.config.get('input files', 'tree') # -t
        self.gene_families_path = self.config.get('input files', 'gene families', None) # -i
        self.lambda_tree_path = self.config.get('input files', 'lambda tree', None) # -y
        self.root_dist_path = self.config.get('input files', 'root distribution', None) # -f

        # options
        self.simulate = self.config.getboolean('options', 'simulate') # -s
        self.simulationsN = self.config.get('options', 'simulationsN', 0) # -n
        self.gammacatN = self.config.get('options', 'gammacatN', 1) # -k

        # parameters
        self.a_lambda = self.config.get('parameters', 'lambda') # -l
        # self.alpha = self.config.get('parameters', 'alpha') # -a
        
        # building command_line
        self.build_command_line()
        
    def check_options(self):
        if simulate and gene_families_path:
            exit('You can only infer/compute a likelihood or simulate. Exiting\n')

    def add_option(self, option, value):
        """option (e.g., -n), value (e.g., 10)"""

        if value:
            self.command_line += [option, value]

    def build_command_line(self):
        self.command_line.append(os.path.join(self.binary_path, 'cafexp'))

        self.add_option('-l', self.a_lambda)
        self.add_option('-t', self.tree_path)
        self.add_option('-n', self.simulationsN)
        self.add_option('-f', self.root_dist_path)
        self.add_option('-k', self.gammacatN)
        self.add_option('-i', self.gene_families_path)
        
        if self.simulate:
            self.command_line += ['-s']
        
    def pour_me_one(self):
        # print " ".join(self.command_line)
        return subprocess.check_output(" ".join(self.command_line), shell=True)
        
if __name__ == "__main__":
    # sys.argv[1] = config file

    barista = Barista(sys.argv[1])
    print barista.pour_me_one()
