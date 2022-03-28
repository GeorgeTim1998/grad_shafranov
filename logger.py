import logging
from termcolor import colored

logging.basicConfig(filename='test.log', level = logging.INFO,
                    format='%(asctime)s:%(message)s')

def log_n_output(message, color):
    logging.info(message)
    print_colored(message, color)
    
def log_n_output_return_carriage(message, color):
    logging.info(message)
    print_colored(message + "\n", color)

def info(message):
    logging.info(message)
    
def print_colored(text, color):
    print(colored(text, color))