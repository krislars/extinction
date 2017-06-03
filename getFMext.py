def getFMext(x,R):
    
    """ Fitzpatrick 1999 
        Input: wavelength in microns, nominal R = 3
        Output: Al/EBV, so divide by R to get Al/AV
    """
    
    from scipy.interpolate import interp1d
    
    import numpy as np

    f99_anchor = 1.0E4 / np.array([np.inf, 26500., 12200., 6000., 5470., 4670., 4110.]) #microns
    
    a6000= -0.426 +1.0044*R
    a5470= -0.050 +1.0016*R
    a4670=  0.701 +1.0016*R
    a4110=  1.208 +1.0032*R -0.00033*R**2. # typo in the paper!

    af99_anchor = np.array([0.0, 0.265, 0.829, a6000, a5470, a4670, a4110])
        
    f=interp1d(f99_anchor,af99_anchor, kind='cubic') 
    
    return f(1.0/x)

