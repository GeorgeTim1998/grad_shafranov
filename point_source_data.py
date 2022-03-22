import numpy
from termcolor import colored

class PointSource:
    def __init__(self, alpha):
        # watch out for missing ',' in arrays
        disp = 0.01
        self.r = [
            [0.1, 0.34],
            [0.1, -0.34],
            [0.17, 0.37],
            [0.17, -0.37],
            [0.4, 0.37],
            [0.4, -0.37],
            [0.53, 0.185],
            [0.53, -0.185]
        ] #array of vector point sources position 
        self.i_disp = [
            [-alpha * 0, disp],
            [-alpha * 0, disp],
            [-alpha * 39.1, disp],
            [-alpha * 39.1, disp],
            [-alpha * 62.6, disp],
            [-alpha * 62.6, disp],
            [-alpha * 5.3, disp],
            [-alpha * 5.3, disp]
        ] #array of [current of the point source, characteristic decay distance]
        
        poor_input = (numpy.shape(self.r) != numpy.shape(self.i_disp)) # check if you input all
        #data concerning each point source
        
        if poor_input:
            print(colored("WARNING! Size of self.r and self.i_disp arrays are not the same!\n", 'red'))
        