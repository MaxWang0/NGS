#!/bin/sh
#Last-modified: 29 Jan 2015 03:07:03 PM

####################### Module/Scripts Description ######################
#  
#  Copyright (c) 2008 Yu Wang <yxw124430@utdallas.com>
#  
#  This code is free software; you can redistribute it and/or modify it
#  under the terms of the BSD License (see the file COPYING included with
#  the distribution).
#  
#  @status:  experimental
#  @version: $Revision$
#  @author:  Yu Wang
#  @contact: yxw124430@utdallas.com
#
#########################################################################


g++ -I/usr/include/gsl -lgsl -lgslcblas eSNVlm.cpp -o eSNVlm
