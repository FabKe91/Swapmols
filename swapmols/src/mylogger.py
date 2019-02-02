''' Handles all the logging '''

import logging
logging.getLogger().addHandler(logging.NullHandler())
LOGGER = logging.getLogger("mylogger")
LOGGER.setLevel(logging.INFO)
SH = logging.StreamHandler()
SH.setLevel(logging.WARNING)
SH.setFormatter(logging.Formatter('%(asctime)s %(name)s - %(levelname)s - %(message)s'))
FH = logging.FileHandler("SwapMols.log")
FH.setLevel(logging.WARNING)
FH.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s\n'))
LOGGER.addHandler(SH)
LOGGER.addHandler(FH)

print("In logger:", LOGGER.handlers)

def set_verbosity(slevel: int):
    ''' Sets level of streamhandler
        slevels:
            0 Warning
            1 Info
            2 Debug
    '''
    if slevel is None:
        slevel = 0
    elif slevel > 2:
        slevel = 2
    slevel = 30 -slevel*10
    SH.setLevel(slevel)#logging.__getattribute__(slevel))


def set_loglevel(slevel: int):
    ''' Sets level of FileHandler
        see set_verbosity
    '''
    slevel = 30 -slevel*10
    FH.setLevel(slevel)#logging.__getattribute__(slevel))
