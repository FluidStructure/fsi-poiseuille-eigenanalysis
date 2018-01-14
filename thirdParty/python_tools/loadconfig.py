#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Loads in the relative C-fmm and pythontools paths.

AUTHOR - Jarrad Kapor (jarrad.kapor@postgrad.curtin.edu.au)
MODIFIED - 2009.08

"""

import os
import sys
BASE_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
CFMM_PATH = os.path.join(BASE_PATH, 'cfmm-2d.git')
PYTOOLS_PATH = os.path.join(BASE_PATH, 'pythontools.git')
sys.path += [CFMM_PATH, PYTOOLS_PATH]
