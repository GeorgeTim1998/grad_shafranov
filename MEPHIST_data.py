class MEPhIST:
    # in SI
    def __init__(self):
        self.R0 = 0.25
        self.a = 0.13
        self.k = 2
        self.e = self.R0/self.a
        
        self.B0 = 1
        self.I = 160e3
        
        self.p_axis = 12.5e3
        
        self.psi_axis = pow(self.R0, 2)*self.B0
        
        