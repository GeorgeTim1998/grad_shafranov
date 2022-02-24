class MEPhIST:
    # in SI
    def __init__(self):
        self.R0 = 0.25 # Большой радиус плазмы
        self.a = 0.13 # Малый радиус плазмы
        self.k = 2 # Вытянутость плазмы
        self.delta = 0.4 # вытянутость
        self.e = self.R0/self.a # Аспектное соотношение
        
        self.B0 = 1 # Тороидальное поле на оси
        self.I = 160e3 # Ток плазмы
        
        self.p_axis = 12.5e3 # Давление плазмы на оси
        
        self.psi_axis = pow(self.R0, 2)*self.B0 # Пси на оси
