# Various functions
import os

# ===============================================================================
# determine expnos automatically
# ===============================================================================
def autoexpnos(path):
    expnos = sorted(os.listdir(path))
    for expno in expnos[:]:  # [:] is to create a copy of expnos 
        # which will not be modified when the original list is changed
        if expno.startswith('.'):
            expnos.remove(expno)
    return expnos

