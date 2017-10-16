import os
import sys
from cfg_maker import *

class MyArgs:
  def __init__(self, bin_path='./', input_file_path=None, tree_file_path=None, poisson_lambda=None, lambdas_list=[], alphas_list=[]):
      self.cafexp_path = bin_path
      self.input_args = 'i,t'
      self.input_values = input_file_path+','+tree_file_path
      self.options_args, self.options_values = '', ''
      self.output_args, self.output_values = 's', ''
      self.parameters_values, self.parameters_args = '', ''
      
      if poisson_lambda:
          self.parameters_args += 'p,'
          self.parameters_values = poisson_lambda+','
          self.output_values = 'p'+poisson_lambda

      self.parameters_args += 'l'

      # default values
      self.number_runs = 1
      self.instr_path = ''
      self.output_path = './'

def batch_cfg_maker(bin_path, input_file_path, tree_file_path, poisson_lambda, lambdas_list, alphas_list):
    myargs = MyArgs(bin_path, input_file_path, tree_file_path, poisson_lambda, lambdas_list, alphas_list)
    # args = {'cafexp_path':bin_path, 'input_args':'i,t', 'input_values':input_file_path+','+tree_path, 'parameters_args':'l', 'instr_path':'', 'number_runs':'1'}

    for idx, l in enumerate(lambdas_list):

        with_alpha = False
        if alphas_list:
            myargs.parameters_args += 'l,a'
            with_alpha = True
            
            for idx2, a in enumerate(alphas_list):
                par_value = l+','+a
                output_value = 'l'+l+'a'+a
                myargs.parameters_values += par_value
                myargs.output_values += output_value
                myarg_stash = ArgStash(myargs)
                myargs.parameters_values = myargs.parameters_values[:-len(par_value)]
                myargs.output_values = myargs.output_values[:-len(output_value)]
                new_config_file_name = 'config_file_' + str(idx)+'-'+str(idx+1)+'.cfg'
                os.system('mv config_file_0.cfg ' + new_config_file_name)
                os.system('python scripts/barista.py ' + new_config_file_name) # running CAFExp

        if not with_alpha:
            par_value = l
            output_value = 'l'+l
            myargs.parameters_values += par_value
            myargs.output_values += output_value
            myarg_stash = ArgStash(myargs)
            myargs.parameters_values = myargs.parameters_values[:-len(par_value)]
            myargs.output_values = myargs.output_values[:-len(output_value)]
            new_config_file_name = 'config_file_' + str(idx+1)+'.cfg'
            os.system('mv config_file_0.cfg ' + new_config_file_name)
            os.system('python scripts/barista.py ' + new_config_file_name) # running CAFExp
    
if __name__ == "__main__":
    # sys.argv[1]: bin path
    # sys.argv[2]: input data
    # sys.argv[3]: tree
    # sys.argv[4]: poisson lambda (or 'None')
    # sys.argv[4]: lambda values separated by ','
    # sys.argv[5]: alpha values separated by ','

    poisson_lambda = None
    if not sys.argv[4] == 'None':
        poisson_lambda = sys.argv[4]
    lambdas_list = sys.argv[5].split(',')
    alphas_list = list()
    try:
        alphas_list = sys.argv[6].split(',')
    except:
        print 'No alphas were provided.'
        
    batch_cfg_maker(sys.argv[1], sys.argv[2], sys.argv[3], poisson_lambda, lambdas_list, alphas_list)
