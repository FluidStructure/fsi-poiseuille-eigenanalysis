#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Miscellaneous useful functions

AUTHOR - Jarrad Kapor (jarrad.kapor@postgrad.curtin.edu.au)
MODIFIED - 2010.01

"""
import inspect

def where_am_i(debug_level=1):
    frame_info = inspect.getframeinfo(inspect.currentframe(1))

    prepend = ' ' * debug_level
    message = prepend + frame_info.filename + "->" + frame_info.function + "():"+ str(frame_info.lineno) + ':'
    
    return message 

def main():
    pass

if __name__ == '__main__':
    main()
