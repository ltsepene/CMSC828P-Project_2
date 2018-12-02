#!/usr/bin/env python

# Load required modules
import sys, os

################################################################################
# LOGGING
################################################################################
import logging

# Logging format
FORMAT = '%(asctime)s ReproducingLiu2017 %(levelname)-10s: %(message)s'
logging.basicConfig(format=FORMAT)

def get_logger(verbosity=logging.INFO):
    '''
    Returns logger object
    '''
    logger = logging.getLogger(__name__)
    logger.setLevel(verbosity)
    return logger