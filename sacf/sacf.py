import numpy as np 
import moch
import soch
import os
import sys
import scipy.io
import thorns

def main(parseID):
    
    parseIn  = parseID + 'In.mat'
    parseOut = parseID + 'Out.mat' 

    parse = scipy.io.loadmat(parseIn)

    os.remove(parseIn)

    lagSpace = 1. * parse['lagSpace'] / 1000
    
    parsStruct = parse['pars'][0, 0]

    # Parametres 
    est = {'duration' : 1. * parsStruct['est'][0,0]['dur'][0][0] / 1000,
           'loudness' : 70.,
           'onset'    : 1. * parsStruct['est'][0,0]['onset' ][0][0] / 1000,  
           'tail'     : 1. * parsStruct['est'][0,0]['tail'][0][0] / 1000,
           'maskN'    : 0,
           'filename' : -1,
           'bandpass' : False,
           'save'     : False
          }

    if est['filename'] == -1:
        est['type']     = parsStruct['est'][0,0]['type'][0]
        est['freq']     = parsStruct['est'][0,0]['f'][0][0]
        est['harms']    = 0
        est['harmFact'] = 0
        est['shift']    = parsStruct['est'][0,0]['shift'][0][0]
        est['nOfIts']   = parsStruct['est'][0,0]['nOfIts'][0][0]
        est['notes']    = 0
        est['tuning']   = 0
        est['noiseOff'] = 1. * parsStruct['est'][0,0]['noiseOff'][0][0] / 1000
    else:
        est['type']     = 'external'


    par = {'periphFs'   : 100000,
           'cochChanns' : (125, 10000, 30),
           'SACFTau'    : 1. * parsStruct['tauSACF'][0,0] / 1000,
           'subCortTau' : 1. * parsStruct['tauSubthal'][0,0] / 1000,
           'solvOnset'  : 1. * parsStruct['solvOnset'][0] / 1000,
           'subCortFs'  : 100000,
           'subCortAff' : parsStruct['subCortAff'][0,0],
           'regularise' : parsStruct['regularise'][0,0],
           'mu0'        : parsStruct['mu0'][0,0],
           'SACFGround' : parsStruct['SACFGround'][0,0],
           'cortFs'     : parsStruct['cortFs'][0,0],
          }

    [A, n, b] = thalamicInput(lagSpace, par, est)

    duration = 1.* len(A) / par['cortFs']

    dti = 1./par['cortFs'] 
    timeSpace = np.arange(start = dti, stop = duration + dti, step = dti)
    
    if 'off' in est.keys():
        timeSpace = timeSpace - est['off']

    scipy.io.savemat(parseOut, {'A':A, 'n':n, 'b':b, 'timeSpace': timeSpace})



def thalamicInput(lagSpace, par, est, raster = False):

    fs = par['periphFs']

    # Subcortical processing
    sound = soch.createStimulus(est, par['periphFs'])
    prob = moch.peripheral(sound, par)

    [A, n, b] = moch.subcortical(prob, lagSpace, par)

    for i in range(1, par['subCortAff']):
        sound = soch.createStimulus(est, par['periphFs'])
        prob = moch.peripheral(sound, par)
        [A0, n0, b0] = moch.subcortical(prob, lagSpace, par)
        A = A + A0
        n = n + n0
        b = b + b0

    A = (1. / par['subCortAff']) * A
    n = (1. / par['subCortAff']) * n
    b = (1. / par['subCortAff']) * b

    if raster:    
        anfTrains = moch.peripheralSpikes(sound, par, fs = -1)
        thorns.plot_raster(anfTrains)
        thorns.show()

    return [A, n, b]



main(sys.argv[1])


