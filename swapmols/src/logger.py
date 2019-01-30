''' Handles all the logging '''

import logging

LOGGER = logging.getLogger()
LOGGER.setLevel(logging.DEBUG)
SH = logging.StreamHandler()
SH.setLevel(logging.INFO)
SH.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
LOGGER.addHandler(SH)
FH = logging.FileHandler("SwapMols.log")
FH.setLevel(logging.INFO)
FH.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s\n'))
LOGGER.addHandler(FH)

def set_verbosity(slevel: str):
    ''' Sets level of streamhandler
        slevel can be DEBUG, INFO, WARNING, ...
    '''
    SH.setLevel(logging.__getattribute__(slevel))


def set_loglevel(slevel: str):
    ''' Sets level of FileHandler
        slevel can be DEBUG, INFO, WARNING, ...
    '''
    FH.setLevel(logging.__getattribute__(slevel))
