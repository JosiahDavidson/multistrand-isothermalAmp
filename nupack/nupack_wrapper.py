
# Code contributors:  Erik Winfree, Chris Thachuk, Justin Bois, Joseph Berleant.
#
# The following functions are currently wrapped:
#  pfunc
#  pairs
#  mfe
#  energy
#  prob
#  defect
#  sample
#
# The following functions may be wrapped in a future release:
#  complexes
#  count
#  concentrations
#  design
#  distributions
#  subopt



import math
import subprocess as sub
import os

def dGadjust(T,N):
    """Adjust NUPACK's native free energy (with reference to mole fraction units) to be appropriate for molar units, assuming N strands in the complex."""
    R=0.0019872041 # Boltzmann's constant in kcal/mol/K 
    water=55.14    # molar concentration of water at 37 C, ignore temperature dependence, which is about 5%
    K=T+273.15     # Kelvin
    adjust = R*K*math.log(water) # converts from NUPACK mole fraction units to molar units, per association
    return adjust*(N-1)

def get_nupack_exec_path(exec_name):
  """ If the NUPACKHOME environment variable is set, use that as the directory
  of the Nupack executables. Otherwise, have Python search the PATH directly. """
  if 'NUPACKHOME' in os.environ:
    return os.environ['NUPACKHOME'] + '/bin/' + exec_name;
  else:
    return exec_name;
    
def call_with_file(args, cmd_input, outsuffix):
  """ Performs a NUPACK call, returning the lines of the output in a temporary
  output file. The output file is assumed to have the suffix 'outsuffix'.
  outsuffix includes the period (.) delimiter.
    Ex:
      call_with_file(args, input, '.sample')
  """
  
  import tempfile

  ## Preliminaries
  # Set up temporary output file
  outfile = tempfile.NamedTemporaryFile(delete=False, suffix=outsuffix)
  outprefix = outfile.name[:-len(outsuffix)]
  
  # Close the output file so sample can open/write to it.
  # Will reopen it later to get the output.
  outfile.close()
  
  ## Perform executable call, ignoring pipe output
  args = [str(s) for s in args] # all argument elements must be strings
  cmd_input = outprefix + '\n' + cmd_input # prepend the output file prefix to the input for NUPACK
  p = sub.Popen(args, stdin=sub.PIPE, stdout=sub.PIPE, stderr=sub.STDOUT)
  p.communicate(cmd_input)
  
  ## Process and return output
  # Read output file and clean it up
  # Note that it was created by us, so it won't be cleaned up automatically
  out = open(outfile.name, "rt")
  output_lines = out.readlines()
  out.close()
  os.remove(outfile.name)
  
  return output_lines
  

def call_with_pipe(args, cmd_input):
  """ Performs a NUPACK call, returning the lines of the output from the pipe.
  """
  args = [str(s) for s in args] # all argument elements must be strings
  
  p=sub.Popen(args,stdin=sub.PIPE,stdout=sub.PIPE,stderr=sub.PIPE)
  output,error = p.communicate(cmd_input)
  output_lines = output.split('\n')
  return (output_lines, error)
    
def pfunc(sequences, material, dangles = 'some', T = 37, multi = True, pseudo = False, sodium = 1.0, magnesium = 0.0):
  """Calls NUPACK's pfunc on a complex consisting of the unique strands in sequences, returns dG.
       sequences is a list of the strand sequences
       See NUPACK User Manual for information on other arguments. """
  
  ## Set up command-line arguments and input
  args = [get_nupack_exec_path('pfunc'),
          '-material', material,   '-sodium', sodium,
          '-magnesium', magnesium, '-dangles', dangles, '-T', T]
  if multi: args += ['-multi']
  if pseudo: args += ['-pseudo']
  
  if not multi:
    cmd_input = '+'.join(sequences)
  else:
    n_seqs = len(sequences)
    seq_order = ' '.join([str(i) for i in range(1, n_seqs+1)])
    cmd_input = str(n_seqs) + '\n' + ('\n'.join(sequences)) + '\n' + seq_order
  
  ## Perform call until it works (should we have a max # of tries before quitting?)
  output, error = call_with_pipe(args, cmd_input)
  while len(output) < 4 : # can't figure out why, but occasionally NUPACK returns empty-handed.  Subsequent tries seem to work...
      print 'Retrying pfunc: NUPACK failed with output ' + `output` + ' and error ' + `error` +" ."
      output, error = call_with_pipe(args, cmd_input)

  ## Parse and return output
  if output[-4] != "% Free energy (kcal/mol) and partition function:" :
      raise NameError('NUPACK output parsing problem')

  if float(output[-3])==float('inf') : return 0               # if these strands can't base-pair
  else: return float(output[-3]) + dGadjust(T,len(sequences))
  
  
def pairs(sequences, material, dangles = 'some',
          T = 37, multi = True, pseudo = False,
          sodium = 1.0, magnesium = 0.0, cutoff = 0.001):
  """Calls NUPACK's pairs executable on a complex consisting of the unique strands in sequences.
     Returns the probabilities of pairs of bases being bound, only including those pairs
     with probability greater than cutoff.
       sequences is a list of the strand sequences
       See NUPACK User Manual for information on other arguments.
  """
  
  ## Set up command-line arguments and input
  args = [get_nupack_exec_path('pairs'),
          '-material', material,   '-cutoff', cutoff,   '-sodium', sodium,
          '-magnesium', magnesium, '-dangles', dangles, '-T', T]
  if multi: args += ['-multi']
  if pseudo: args += ['-pseudo']
  
  if not multi:
    cmd_input = '+'.join(sequences)
    suffix = '.ppairs'
  else:
    n_seqs = len(sequences)
    seq_order = ' '.join([str(i) for i in range(1, n_seqs+1)])
    cmd_input = str(n_seqs) + '\n' + ('\n'.join(sequences)) + '\n' + seq_order
    suffix = '.epairs'
  
  ## Perform call
  output = call_with_file(args, cmd_input, suffix)

  ## Parse and return output
  pair_probs = []
  for l in filter(lambda x: x[0].isdigit(), output):
    if len(l.split()) > 1:
      pair_probs.append(tuple(l.split()))
  
  return pair_probs

  
def mfe(sequences, material, dangles = 'some',
        T = 37, multi = True, pseudo = False,
        sodium = 1.0, magnesium = 0.0, degenerate = False):
  """Calls NUPACK's mfe executable on a complex consisting of the unique strands in sequences.
     Returns the minimum free energy structure, or multiple mfe structures if the degenerate
     option is specified
       sequences is a list of the strand sequences
       degenerate is a boolean specifying whether to include degenerate mfe structures
       See NUPACK User Manual for information on other arguments.
  """
  
  ## Set up command-line arguments and input
  args = [get_nupack_exec_path('mfe'),
          '-material', material,   '-sodium', sodium,
          '-magnesium', magnesium, '-dangles', dangles, '-T', T]
  if multi: args += ['-multi']
  if pseudo: args += ['-pseudo']
  if degenerate: args += ['-degenerate']
  
  if not multi:
    cmd_input = '+'.join(sequences)
  else:
    n_seqs = len(sequences)
    seq_order = ' '.join([str(i) for i in range(1, n_seqs+1)])
    cmd_input = str(n_seqs) + '\n' + ('\n'.join(sequences)) + '\n' + seq_order
  
  ## Perform call
  output = call_with_file(args, cmd_input, '.mfe')

  ## Parse and return output
  structs = []
  for i, l in enumerate(output):
    if l[0] == '.' or l[0] == '(':
      s = l.strip()
      e = output[i-1].strip()
      structs.append((s,e))
  
  return structs
  
  
def energy(sequences, structure, material, dangles = 'some', T = 37, multi = True, pseudo = False, sodium = 1.0, magnesium = 0.0):
  """Calls NUPACK's energy executable. Returns the microstate dG.
       sequences is a list of the strand sequences
       structure is a string with the dot-paren structure notation
         (pair-list notation for structures is not currently supported)
       See NUPACK User Manual for information on the other arguments.
  """
  
  ## Set up command-line arguments and input
  args = [get_nupack_exec_path('energy'),
          '-material', material,   '-sodium', sodium,
          '-magnesium', magnesium, '-dangles', dangles, '-T', T]
  if multi: args += ['-multi']
  if pseudo: args += ['-pseudo']
  
  if not multi:
    cmd_input = '+'.join(sequences) + '\n' + structure
  else:
    n_seqs = len(sequences)
    seq_order = ' '.join([str(i) for i in range(1, n_seqs+1)])
    cmd_input = str(n_seqs) + '\n' + ('\n'.join(sequences)) + '\n' + seq_order + '\n' + structure
                
  ## Perform call
  output, error = call_with_pipe(args, cmd_input)

  ## Parse and return output
  if output[-3] != "% Energy (kcal/mol):" :
     raise ValueError('NUPACK output parsing problem')

  return float(output[-2])


def prob(sequences, structure, material, dangles = 'some', T = 37, multi = True, pseudo = False, sodium = 1.0, magnesium = 0.0):
  """Calls NUPACK's prob executable. Returns the probability of the given structure.
       sequences is a list of the strand sequences
       structure is a string with the dot-paren structure notation
         (pair-list notation for structures is not currently supported)
       See NUPACK User Manual for information on the other arguments.
  """
  
  ## Set up command-line arguments and input
  args = [get_nupack_exec_path('prob'),
          '-material', material,   '-sodium', sodium,
          '-magnesium', magnesium, '-dangles', dangles, '-T', T]
  if multi: args += ['-multi']
  if pseudo: args += ['-pseudo']
  
  if not multi:
    cmd_input = '+'.join(sequences) + '\n' + structure
  else:
    n_seqs = len(sequences)
    seq_order = ' '.join([str(i) for i in range(1, n_seqs+1)])
    cmd_input = str(n_seqs) + '\n' + ('\n'.join(sequences)) + '\n' + seq_order + '\n' + structure
                
  ## Perform call
  output, error = call_with_pipe(args, cmd_input)

  ## Parse and return output
  if output[-3] != "% Probability:" :
     raise ValueError('NUPACK output parsing problem')

  return float(output[-2])


def defect(sequences, structure, material,
           dangles = 'some', T = 37, multi = True, pseudo = False,
           sodium = 1.0, magnesium = 0.0, mfe = False):
  """Calls NUPACK's defect executable. Returns the ensemble defect (default) or the mfe defect.
       sequences is a list of the strand sequences
       structure is a string with the dot-paren structure notation
         (pair-list notation for structures is not currently supported)
       See NUPACK User Manual for information on the other arguments.
  """
  
  ## Set up command-line arguments and input
  args = [get_nupack_exec_path('defect'),
          '-material', material,   '-sodium', sodium,
          '-magnesium', magnesium, '-dangles', dangles, '-T', T]
  if multi: args += ['-multi']
  if pseudo: args += ['-pseudo']
  if mfe: args += ['-mfe']
  
  if not multi:
    cmd_input = '+'.join(sequences) + '\n' + structure
  else:
    n_seqs = len(sequences)
    seq_order = ' '.join([str(i) for i in range(1, n_seqs+1)])
    cmd_input = str(n_seqs) + '\n' + ('\n'.join(sequences)) + '\n' + seq_order + '\n' + structure
                
  ## Perform call
  output, error = call_with_pipe(args, cmd_input)

  ## Parse and return output
  if "% Ensemble defect" not in output[-4] and \
        "% Fraction of correct nucleotides vs. MFE" not in output[-4]:
    raise ValueError('NUPACK output parsing problem')

  # We don't return the normalized ensemble defect, because that is easily calculable on your own
  return float(output[-3])


def sample(samples, sequences, material,
           dangles = 'some', T = 37, multi = True,
           pseudo = False, sodium = 1.0, magnesium = 0.0):
  """ Calls the NUPACK sample executable.
        samples is the number of Boltzmann samples to produce.
        See NUPACK User Manual for information on the other arguments.
      This only works with NUPACK 3.0.2+
      
      Note that if using OS X and sample is not in your $PATH, this will try
      to run the standard BSD tool 'sample'. """
      
  ## Set up command-line arguments and input
  args = [get_nupack_exec_path('sample'),
          '-material', material, '-sodium', sodium, '-magnesium', magnesium,
          '-dangles', dangles, '-T', T, '-samples', samples]
  if multi: args += ['-multi']
  if pseudo: args += ['-pseudo']
  
  if not multi:
    cmd_input = '+'.join(sequences);
  else:
    n_seqs = len(sequences)
    seq_order = ' '.join([str(i) for i in range(1, n_seqs+1)])
    cmd_input = str(n_seqs)          + '\n' + \
                '\n'.join(sequences) + '\n' + \
                seq_order

  # Call executable
  output = call_with_file(args, cmd_input, '.sample')

  # Check NUPACK version
  if not "NUPACK 3.0" in output[0]:
    raise IOError("Boltzmann sample function is not up to date. NUPACK 3.0.2 or greater needed.")

  # Parse and return output
  sampled = [l.strip() for l in output[14:]]
  return sampled
