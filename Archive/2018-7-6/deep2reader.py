import numpy as np




class d2r:

    def __init__(self):

        self.OIIEW = np.loadtxt('./DEEP2/O2EW_histogram.txt')
        self.OIIIHb = np.loadtxt('./DEEP2/O3Hb_histogram.txt')
        self.OIIIOII = np.loadtxt('./DEEP2/O3O2_histogram.txt')

