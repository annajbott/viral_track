""" =====================
Pipeline viral track
=========================

Overview 
=========

This pipeline was developed to ...


Requires:

* fastq files etc...
* Reference genome 
* ... 

Pipeline output
================
<What comes out of pipeline...>


Code
=====

"""

from ruffus import *

import sys
import os

import cgatcore.pipeline as P
import cgatcore.experiment as E


# load options from the config file
PARAMS = P.get_parameters(
    ["%s/pipeline.yml" % os.path.splitext(__file__)[0],
     "../pipeline.yml",
     "pipeline.yml"])






