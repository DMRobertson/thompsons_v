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

def find_examples_passing(test_func,
                          automorphism_generator,
                          test_name,
                          save_examples=True,
                          description=''):
	"""Generates random examples using the specified *automorphism_generator* and applies the user-supplied *test_func* to see if they match a certain condition. If *test_func* determines the example matches the condition, it should return a **non-empty** string giving the details. Otherwise if the condition is not met, *test_func* should return an empty string.
	
	When an example is found, its details are saved to disk if *save_examples* is True.
	"""
	os.makedirs(test_name, exist_ok=True)
	log_filepath = os.path.join(test_name, test_name + '.log')
	output_path  = os.path.join(test_name, test_name + '_{}.aut')
	
	prepare_logging(log_filepath)
	
	logging.info('Starting a search for examples: {}'.format(test_name))
	if description:
		logging.info(description)
	
	num_attempts = 0
	num_found = 0
	while True:
		num_attempts += 1
		if num_attempts % 1000 == 0:
			logging.debug('Attempt number {}'.format(num_attempts))
		
		aut = automorphism_generator()
		try:
			result = test_func(aut)
		except Exception as e:
			logging.error("An exception occured when calling the test function on attempt {}.\n\t{}").format(
			  num_attempts, e)
			logging.debug(traceback.format_exc())
		
		assert isinstance(result, str), "test_func should return the empty string for a fail \
		  and a non-empty string for a pass."
		if result:
			num_found += 1
			logging.info('Found example number {} on attempt {}. Success proportion: {:.2%}'.format(
			  num_found, num_attempts, num_found / num_attempts))
		if not (result and save_examples):
			continue
		
		path = output_path.format(num_found)
		logging.debug('Saving the details of example number {} to {}.'.format(
		  num_found, path))
		aut.save_to_file(path, result)