''' Handles all the logging '''

import logging
logging.getLogger().addHandler(logging.NullHandler())
LOGGER = logging.getLogger("Swapmols")
LOGGER.setLevel(logging.DEBUG)

def add_handlers(log):
    sh = logging.StreamHandler()
    sh.setLevel(logging.INFO)
    sh.setFormatter(logging.Formatter('%(asctime)s %(name)s - %(levelname)s - %(message)s'))
    log.addHandler(sh)

    fh = logging.FileHandler("SwapMols.log")
    fh.setLevel(logging.INFO)
    fh.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s\n'))
    log.addHandler(fh)
    return sh, fh


def set_verbosity(slevel: int, handler):
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

    handler.setLevel(slevel)#logging.__getattribute__(slevel))


#def set_loglevel(slevel: int):
#    ''' Sets level of FileHandler
#        see set_verbosity
#    '''
#    slevel = 30 -slevel*10
#    FH.setLevel(slevel)#logging.__getattribute__(slevel))
