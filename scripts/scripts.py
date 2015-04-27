import sys, traceback, logging, os, sys
from pathlib import Path

"""A bunch of utility functions which make it convenient for creating scripts using thompson."""
def setup_script(__file__):
	"""A naughty hack which makes thompson importable by examining __file__ and adding to sys.path."""
	path = Path(__file__)
	try:
		parent = str(path.parents[1])
	except IndexError:
		parent = str(Path('..').resolve())
	sys.path.insert(0, parent)

def choose_from_enum(enum, desc = None):
	"""An input loop which has a console user select an element of an enumeration."""
	if desc is None:
		desc = {}
	num_choices = len(enum)
	print('Choices are:')
	for i, x in enumerate(enum):
		print('\t[{}] {}. {}'.format(i + 1, x.name, desc.get(x, '')))

		ans = None
	while ans is None:
		try:
			ans = int(input('Please make a choice by entering an integer in [1--{}]: '.format(
			  num_choices
			)))
		except Exception as e:
			print(e)
		else:
			if not (1 <= ans <= num_choices):
				ans = None
	return enum(ans)

def prepare_logging(log_filepath):
	"""This is how I prefer to setup the logging module."""
	logging.basicConfig(level   =logging.DEBUG,
	                    format  ='%(asctime)s %(levelname)-8s %(message)s',
	                    filename=log_filepath,
	                    filemode='wt')
	console = logging.StreamHandler()
	console.setLevel('INFO') #show info and above
	formatter = logging.Formatter('%(levelname)-8s %(message)s')
	console.setFormatter(formatter)
	logging.getLogger().addHandler(console)

def find_examples_passing(test_functions,
                          automorphism_generator,
                          test_name,
                          save_examples=True,
                          description='',
                          max_examples=float('inf')):
	"""Generates random examples using the specified *automorphism_generator* and applies the user-supplied *test_functions* to see if they match certain conditions. If a *test_function* decides that the example matches the condition, it should return a **non-empty** string giving the details. Otherwise the condition is not met, and the *test_function* should return an empty string.
	
	When an example is found, its details are saved to disk if *save_examples* is True. At most *max_examples* are found in this way.
	"""
	os.makedirs(test_name, exist_ok=True)
	log_filepath = os.path.join(test_name, test_name + '.log')
	output_path  = os.path.join(test_name, test_name + '_{}.aut')
	time_series_path = os.path.join(test_name, test_name + '.csv')
	
	prepare_logging(log_filepath)
	
	logging.info('Starting a search for examples: {}'.format(test_name))
	if description:
		logging.info(description)
	
	num_attempts = 0
	num_passes   = [0]  * len(test_functions)
	results      = [''] * len(test_functions)
	num_found    = 0
	
	with open(time_series_path, 'wt') as series:
		headers = ["Attempt"]
		for i in range(len(test_functions)):
			headers.append( 'Test{}Passes'.format(i+1) )
		print(', '.join(headers), file=series)
		
		while num_found < max_examples:
			num_attempts += 1
			if num_attempts % 4000 == 0:
				logging.info('Attempt number {}'.format(num_attempts))
				logging.debug('Test results so far: {}'.format(num_passes))
			aut = automorphism_generator()
			
			#call all the test functions
			try:
				passed_all = True
				for index, test in enumerate(test_functions):
					results[index] = test(aut)
					assert isinstance(results[index], str), "test_functions should return the empty string for a fail \
					  and a non-empty string for a pass."
					if results[index]:
						num_passes[index] += 1
					else:
						passed_all = False
						break
			except Exception as e:
				logging.error("An exception occured when calling test function {} on attempt {}.\n\t{}".format(
				  index, num_attempts, e))
				logging.debug(aut)
				logging.debug(traceback.format_exc())
			
			if passed_all:
				num_found += 1
				logging.info('Found example number {} on attempt {}. Success proportion: {:.2%}'.format(
				  num_found, num_attempts, num_found / num_attempts))
			if passed_all or index > 0:
				entry = str(num_attempts) + ', '
				entry += ', '.join(  str(passes) for passes in num_passes  )
				print(entry, file = series)
				
			if not (passed_all and save_examples):
				continue
			
			path = output_path.format(num_found)
			logging.debug('Saving the details of example number {} to {}.'.format(
			  num_found, path))
			
			result = '\n\n'.join(r for r in results)
			aut.save_to_file(path, result)
	logging.info('Maxmium number of examples found; ending search.')
"""Here is a snippet of R code to plot the change in the success ratio.

data = read.csv('FILENAME_GOES_HERE.csv', header=FALSE)
attempts = data[,1]
success = data[,2]
plot(attempts, success, type='l')
abline(h=mean(success), col='red', lty=2)
"""