import os
import sys
from cfg_maker import *

class MyArgs:
  def __init__(self, bin_path='./', input_file_path=None, tree_file_path=None, lambdas_list=[], alphas_list=[]):
      self.cafexp_path = bin_path
      self.input_args = 'i,t'
      self.input_values = input_file_path+','+tree_file_path
      self.options_args, self.options_values = '', ''
      self.parameters_args = 'l'
      self.parameters_values = ''
      self.output_args, self.output_values = 's', ''

      # default values
      self.number_runs = 1
      self.instr_path = ''
      self.output_path = './'

def batch_cfg_maker(bin_path, input_file_path, tree_file_path, lambdas_list, alphas_list):
    myargs = MyArgs(bin_path, input_file_path, tree_file_path, lambdas_list, alphas_list)
    # args = {'cafexp_path':bin_path, 'input_args':'i,t', 'input_values':input_file_path+','+tree_path, 'parameters_args':'l', 'instr_path':'', 'number_runs':'1'}

    for idx, l in enumerate(lambdas_list):

        with_alpha = False
        if alphas_list:
            myargs.parameters_args = 'l,a'
            with_alpha = True
            
            for idx2, a in enumerate(alphas_list):
                myargs.parameters_values = l+','+a
                myargs.output_values = 'l'+l+'a'+a
                myarg_stash = ArgStash(myargs)
                new_config_file_name = 'config_file_' + str(idx+1)+'-'+str(idx2+1)+'.txt'
                os.system('mv config_file_1.cfg ' + new_config_file_name)
                # os.system('python scripts/barista.py ' + new_config_file_name) # running CAFExp

        if not with_alpha:
            myargs.parameters_values = l
            myargs.output_values = 'l'+l
            new_config_file_name = 'config_file_' + str(idx+1)+'.txt'
            myarg_stash = ArgStash(myargs)
            os.system('mv config_file_1.cfg ' + new_config_file_name)
            # os.system('python scripts/barista.py ' + new_config_file_name) # running CAFExp
    
if __name__ == "__main__":
    # sys.argv[1]: bin path
    # sys.argv[2]: input data
    # sys.argv[3]: tree
    # sys.argv[4]: lambda values separated by ','
    # sys.argv[5]: alpha values separated by ','

    lambdas_list = sys.argv[4].split(',')
    alphas_list = list()
    try:
        alphas_list = sys.argv[5].split(',')
    except:
        print 'No alphas were provided.'
        
    batch_cfg_maker(sys.argv[1], sys.argv[2], sys.argv[3], lambdas_list, alphas_list)
