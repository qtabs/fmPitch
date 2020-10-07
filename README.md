# fwPitch
FM-Feedback model for pitch and FM-direction encoding 

This folder contains everything you need to reproduce the results of the manuscript
"Neural modelling of the encoding of fast frequency modulation"

The code is written in a mix-up of python and MATLAB (sorry!).

-- Contents --

jibia.m : Wrapper for the plots of all figures in the manuscript

octopus.py : Wrapper for all the modelling results reported in the manuscript

deltaPitchDB.mat : Experimental data for the single-sweeps
deltaPitchSSDB.mat : Experimental data for the sweep-trains

expResAvgSw.mat : Average sweep pitch shift for the single-sweeps
expResAvgTr.mat : Average sweep pitch shift for the sweep-trains

hierField : Model libraries (this folder will be published in github)
	startHere.py : installation instructions, licence, and running examples
	hierField.py : model library
	soch.py : main python librariy for stimulus generation

sacf : Libraries necessary for running the temporal model
	sweepSacf.m : main wrapper, run this to get the results
	sacf.py : python wrapper for the SACF
	moch.py : main python library running the SACF
	soch.py : same as in hierField
