import numpy
from termcolor import colored

class PointSource:
    def __init__(self):
        # watch out for missing ',' in arrays
        disp = 0.01
        self.r = [
            [100, 0],
        ] #array of vector point sources position 
        self.i_disp = [
            [1, disp],
        ] #array of [current of the point source, characteristic decay distance]
        
        poor_input = (numpy.shape(self.r) != numpy.shape(self.i_disp)) # check if you input all
        #data concerning each point source
        
        if poor_input:
            print(colored("WARNING! Size of self.r and self.i_disp arrays are not the same!\n", 'red'))
        