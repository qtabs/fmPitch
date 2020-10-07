import numpy as np
import scipy.integrate
import scipy.signal
import soch
import cochlea
import time
import scipy.io
import matplotlib.pyplot as plt

plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['path.simplify_threshold'] = 0.5
plt.rcParams['savefig.transparent'] = False


def swField(par = 0, est = 0, plotSignal = False, saveMat = ''):

    if par == 0:
        par = loadParameters()

    if est == 0:
        est = stimulusPars()

    W = connectivityMatrices(par)

    p = peripheral(par, est)

    h, s, timeSpace = initialiseStates(p, par)
    
    for tIx in range(1, len(timeSpace)):
        h, s = updateStates(h, s, tIx, p[tIx], par, W)

    if plotSignal:
        print("E[chann] = %.1f %s %.1f".format(tonePred, u"\u00B1", toneVar))
        plotPops(timeSpace, p, h)

    if saveMat:
        saveResultsInMatFile(timeSpace, p, h, est, saveMat)

    return(h)



def swFieldFast(par, p):

    W = connectivityMatrices(par)

    h, s, timeSpace = initialiseStates(p, par)
    
    for tIx in range(1, len(timeSpace)):
        h, s = updateStates(h, s, tIx, p[tIx], par, W)

    return(h)



def updateStates(h, s, t, p, par, W):
    
    x  = synInput(s, t - 1, par, W)
    dh = dhdt(h, x, t, par)
    ds = dsdt(s, h, p, t, par)

    for key in s:
        s[key][t] = s[key][t - 1] + par['dt'] * ds[key]
        s[key][t][s[key][t] < 0] = 0

    for key in h:
        h[key][t] = h[key][t - 1] + par['dt'] * dh[key]
        h[key][t][h[key][t] < 0] = 0

    return([h, s])



def synInput(s, t, par, W):

    x = dict()

    xUp      = par['JsN'] * np.dot(W['upf'],    s['upN'][t])
    xDown    = par['JsN'] * np.dot(W['downf'],  s['downN'][t])
    xPeriph  = par['Jin'] * np.dot(W['periph'], s['inA'][t])
    x['fre'] = xPeriph + xUp + xDown # I0 = 0

    # Delay lines
    wDsUp, wDsDown = delayInputs(s, t, par, W)

    xFreq      = par['JfA'] * wDsUp
    xInh       = par['JsG'] * np.dot(W['downie'], s['downG'][t])
    x['upe']   = xFreq - xInh + par['I0e']

    xFreq      = par['JfA'] * wDsDown
    xInh       = par['JsG'] * np.dot(W['upie'], s['upG'][t])
    x['downe'] = xFreq - xInh + par['I0e']    

    x['upi']   = par['JsA'] * np.dot(W['upei'],   s['upA'][t])   + par['I0i']  
    x['downi'] = par['JsA'] * np.dot(W['downei'], s['downA'][t]) + par['I0i']  

    return(x)



def dhdt(h, x, t, par):

    dh = dict()

    for key in h:
        phix, tauPop = phi(x[key], h[key][t - 1], par, key[-1])
        dh[key] = (phix - h[key][t - 1]) / tauPop
        
    return(dh)



def dsdt(s, h, p, t, par):

    ds = dict()

    ds['inA']   = par['inGain'] * p
    ds['freqA'] = h['fre'][t - 1]
    ds['upA']   = h['upe'][t - 1]
    ds['upN']   = par['gamma'] * (1 - s['upN'][t - 1]  ) * h['upe'][t - 1]
    ds['upG']   = h['upi'][t - 1]
    ds['downA'] = h['downe'][t - 1]
    ds['downN'] = par['gamma'] * (1 - s['downN'][t - 1]) * h['downe'][t - 1]
    ds['downG'] = h['downi'][t - 1]

    tau = {'A': 'tauAMPA', 'N': 'tauNMDA', 'G': 'tauGABA'}

    for key in s:
        leak    = - s[key][t - 1] / par[tau[key[-1]]]
        noise   = par['sigma'] * np.random.randn(len(leak))
        ds[key] = leak + .001 * ds[key] + noise

    return(ds)



def phi(x, h, par, nt):

    y   = par['a'+nt] * x - par['b'+nt]
    phi = y / (1 - np.exp(-par['d'+nt] * y))
    
    if par['tauEff']:
        dphidx  = phi * (1 / y + par['d'+nt] / (1 - np.exp(par['d'+nt] * y)))
        tauFact = np.minimum(1, par['Delta'] * dphidx / (h + 1E-10))
    else:
        tauFact = 1.

    return([phi, tauFact * par['tauPop'+nt]])



def delayInputs(s, t, par, W):

    D0 = int(par['tShift'] / par['dt'])
    wDsUpMat   = np.zeros((par['N'], par['N']))
    wDsDownMat = np.zeros((par['N'], par['N']))
    for i in range(par['N']):
        for j in range(par['N']):
            DT = int((t - abs(i - j) * D0))
            if DT >= 0:
                if i >= j:
                    wDsUpMat[i, j]   = s['freqA'][DT, j] * W['fup'][j, i]
                if j >= i:
                    wDsDownMat[i, j] = s['freqA'][DT, j] * W['fdown'][j, i]

    return([wDsUpMat.sum(1), wDsDownMat.sum(1)])



def initialiseStates(p, par):

    timeSpace = np.linspace(par['dt'], p.shape[0] * par['dt'], p.shape[0])

    h = {'fre'   : np.zeros(p.shape),
         'upe'   : np.zeros(p.shape),
         'upi'   : np.zeros(p.shape),
         'downe' : np.zeros(p.shape),
         'downi' : np.zeros(p.shape),
    }

    s = { 'inA'   : np.zeros(p.shape), # AMPA synaptic gating in the inputs
          'freqA' : np.zeros(p.shape), # AMPA synaptic gating in freq
          'upA'   : np.zeros(p.shape), # AMPA synaptic gating in up
          'upN'   : np.zeros(p.shape), # NMDA synaptic gating in up
          'upG'   : np.zeros(p.shape), # GABA synaptic gating in up
          'downA' : np.zeros(p.shape), # AMPA synaptic gating in down
          'downN' : np.zeros(p.shape), # NMDA synaptic gating in down
          'downG' : np.zeros(p.shape), # GABA synaptic gating in down
        }

    return([h, s, timeSpace])



def peripheral(par = 0, est = 0):

    if est == 0:
        est = stimulusPars()

    if par == 0:
        par = loadParameters()

    outFs = int(1000. / par['dt']) # Factor 1000 to transform seconds in ms

    sound = soch.createStimulus(est, par['cochFs'])
    ANRts = cochlea.run_zilany2014_rate(sound,  par['cochFs'], 
                                        anf_types = ('lsr', 'msr', 'hsr'), 
                                        cf = (125, 10000, par['N']), 
                                        species = 'human', 
                                        cohc = 1,
                                        cihc = 1)

    ANRts = .6 * ANRts['hsr'] + .25 * ANRts['msr'] + .15 * ANRts['lsr']

    if outFs == par['cochFs']:
        p = np.asarray(ANRts)
    else:
        resampleN = int(len(sound) * outFs / par['cochFs'])
        p = scipy.signal.resample(ANRts, num = resampleN)
        p[p < 0] = 0

    return(p)



def plotPops(timeSpace, p, h):

 
    minAct = min(h['fre'].min(), h['upe'].min(), h['upi'].min(), \
                                           h['downe'].min(), h['downi'].min())
    maxAct = max(h['fre'].max(), h['upe'].max(), h['upi'].max(), \
                                           h['downe'].max(), h['downi'].max())


    f, ax = plt.subplots(2, 2, figsize=(15,8))
    

    im = ax[0,0].imshow(h['fre'].T, aspect='auto', origin = 'lower', \
                                                vmin = minAct, vmax = maxAct)
    ax[0,0].set(ylabel = 'freq')
    

    ax[1,0].imshow(h['downe'].T, aspect='auto', origin = 'lower', \
                                                vmin = minAct, vmax = maxAct)
    ax[1,0].set(xlabel = 'time (ms)', ylabel = 'down exc')
    

    ax[0,1].imshow(h['fre'].T, aspect='auto', origin = 'lower', \
                                                vmin = minAct, vmax = maxAct)
    ax[0,1].set(ylabel = 'freq')
    

    ax[1,1].imshow(h['downe'].T, aspect='auto', origin = 'lower', \
                                                vmin = minAct, vmax = maxAct)
    ax[1,1].set(xlabel = 'time (ms)', ylabel = 'down exc')

    f.colorbar(im, ax=ax.ravel().tolist())

    f.canvas.draw()
    labels = [item.get_text() for item in ax[0,0].get_xticklabels()]
    for ix in range(len(labels)):
        if labels[ix] != '':
            labels[ix] = '{}'.format(int(labels[ix]) / 10)
    for axis in [axrow[i] for axrow in ax for i in range(len(axrow))]:
        axis.set_xticklabels(labels)

    plt.show(block=False)



def saveResultsInMatFile(timeSpace, p, h, est, matfile):

    fStim = np.zeros(timeSpace.shape)
    fs = int(1000 / (timeSpace[1] - timeSpace[0]))

    if est['type'] == 'FMsweep':
        f0 = est['freq'] - 0.5 * est['shift']
        f1 = est['freq'] + 0.5 * est['shift']
        tOnset  = int(fs * est['onset'])
        tSweep0 = tOnset + int(fs * est['noiseOff']) 
        tSweep1 = tSweep0 + int(fs * est['duration']) 
        tOffset = tSweep1 + int(fs * est['noiseOff']) 
        fSweep  = 1. / np.linspace(1/f0, 1/f1, int(fs * est['duration']))

        fStim[0:tOnset]        = -1
        fStim[tOnset:tSweep0]  = f0
        fStim[tSweep0:tSweep1] = fSweep
        fStim[tSweep1:tOffset] = f1
        fStim[tOffset:-1]      = -1


    matDict = {'peripheral' : p,
               'hFrequency' : h['fre'],
               'hUpExcit'   : h['upe'],
               'hUpInhib'   : h['upi'],
               'hDownExcit' : h['downe'],
               'hDownInhib' : h['downi'],
               'fStim'      : fStim,
               'timeSpace'  : timeSpace}

    scipy.io.savemat(matfile, matDict)



def loadParameters():

    par = { # Peripheral
            'dt'      : 0.1,    # ms
            'tShift'  : 1.,     # ms per channel
            'cochFs'  : 100000, # Hz; cochlea and soch are in seconds
            'N'       : 100,    # number of cochlear channels
            'tauPope' : 20.,    # ms; temporal constant of excit neural pops
            'tauPopi' : 10.,    # ms; temporal constant of inhib neural pops
            'tauNMDA' : 100.,   # ms; temporal constant of NMDA gating
            'tauAMPA' : 2.,     # ms; temporal constant of AMPA gating
            'tauGABA' : 5.,     # ms; temporal constant of GABA gating
            'sigma'   : 0.0007, # multiplicative factor for synaptic noise
            'gamma'   : 0.641,  # NMDA coupling factor
            'ae'      : 310,    # transfer function parameter  (exc) 
            'be'      : 125,    # transfer function parameter  (exc) 
            'de'      : 0.16,   # transfer function parameter  (exc)
            'ai'      : 615,    # transfer function parameter  (inh) 
            'bi'      : 177,    # transfer function parameter  (inh) 
            'di'      : 0.087,  # transfer function parameter  (inh)
            'I0e'     : 0.23,   # nA, contant background input
            'I0i'     : 0.10,   # nA, contant background input
            'inGain'  : 1.,     # input gain
            'Jin'     : 0.38,   # AMPA Conductivity peripheral -> freq
            'JfA'     : 0.55,   # AMPA Conductivity freq -> sweep
            'JsN'     : 0.05,   # NMDA Conductivity sweep -> freq
            'JsA'     : 0.67,   # AMPA Conductivity sweep exc -> sweep inh
            'JsG'     : 0.30,   # GABA Conductivity sweep inh -> sweep exc
            'Delta'   : 1.,     # action potential initiation (mV)
            'tauEff'  : True,   # False for a constant tauPop
            'wFeed'   : 0.05,   # as a fraction of the tonotopic range
            'dFeed'   : 0.05,   # as a fraction of the tonotopic range
            'wDelta'  : 0.12,   # as a fraction of the tonotopic range
            'inhSig'  : 0.5,    # as a fraction of the tonotopic range
            'excSig'  : 0.03,   # as a fraction of the tonotopic range
            'inpSig'  : 0.1,    # as a fraction of the tonotopic range
    }

    return(par)



def connectivityMatrices(par = 0, plotMatrices = False):

    if par == 0:
        par = loadParameters()

    W = { # Connectivity
          'periph': np.eye(par['N']),
          'fup'   : np.zeros((par['N'], par['N'])),
          'fdown' : np.zeros((par['N'], par['N'])),
          'upf'   : np.zeros((par['N'], par['N'])),
          'downf' : np.zeros((par['N'], par['N'])),
          'upei'  : np.eye(par['N']),
          'downei': np.eye(par['N']),
          'upie'  : np.ones((par['N'], par['N'])),
          'downie': np.ones((par['N'], par['N']))
    }

    wDelta = par['wDelta'] * par['N'] 
    wFeed  = par['wFeed']  * par['N'] 
    dFeed  = par['dFeed']  * par['N']
    inhSig = par['inhSig'] * par['N']
    inpSig = par['inpSig'] * par['N']
    excSig = par['excSig'] * par['N'] 
    j = 1. * np.arange(par['N'])

    for i in range(par['N']):
        W['periph'][i, :] = np.exp(-(i-j)**2 / (2.*inpSig)) / np.sqrt(inpSig)
        W['fup'  ][i, abs(i-j + wDelta/2) - wDelta/2 <= 0] = 1
        W['fdown'][i, abs(i-j - wDelta/2) - wDelta/2 <= 0] = 1
        
        W['downf'][i, abs(i-j + dFeed + wFeed/2)  - wFeed/2 < 0] = 1
        W['upf'  ][i, abs(i-j - dFeed - wFeed/2)  - wFeed/2 < 0] = 1
        
        W['upie'  ][i, :] = np.exp(-(i-j)**2 / (2.*inhSig))
        W['downie'][i, :] = np.exp(-(i-j)**2 / (2.*inhSig))
        W['upei'  ][i, :] = np.exp(-(i-j)**2 / (2.*excSig))
        W['downei'][i, :] = np.exp(-(i-j)**2 / (2.*excSig))


    if plotMatrices: # Show connectivity matrices
        f, ax = plt.subplots(2, 4)
        ax[0,0].imshow(W['fup'], origin = 'lower')
        ax[0,0].set(title = 'f -> up')
        ax[1,0].imshow(W['fdown'], origin = 'lower')
        ax[1,0].set(title = 'f -> down')
        ax[0,1].imshow(W['upf'], origin = 'lower')
        ax[0,1].set(title = 'up -> f')
        ax[1,1].imshow(W['downf'], origin = 'lower')
        ax[1,1].set(title = 'down -> f')
        ax[0,2].imshow(W['upie'], origin = 'lower')
        ax[0,2].set(title = 'up i - > down e')
        ax[1,2].imshow(W['downie'], origin = 'lower')
        ax[1,2].set(title = 'down i - > up e')
        ax[0,3].imshow(W['upei'], origin = 'lower')
        ax[0,3].set(title = 'up e - > up i')
        ax[1,3].imshow(W['downei'], origin = 'lower')
        ax[1,3].set(title = 'down e - > down i')
        f.show()


    return(W)



def stimulusPars():

    # Stimulus parameters (soch measures time in seconds, not ms!)
    est = {'duration' : 0.05,      # seconds
           'loudness' : 70,        # decibels
           'onset'    : 0.00,      # seconds, silece before onset  
           'tail'     : 0.00,      # seconds, silence after offset
           'maskN'    : 0,         # [0--1], weight of masking noise 
           'filename' : -1,        # string or -1, to use a wavefile as input
           'intv'     : -1,        # seconds or -1, fragment of the file 
           'bandpass' : False,     # [,]Hz or False, bandpass filter
           'save'     : 0,         # 1 to save the stimulus, 0 otherwise 
           'type'     : 'FMsweep', # stimulus type to generate ('PT')
           'freq'     : 1200,       # Hz, average/fundamental frequency
           'shift'    : -300,       # Hz, frequency gap (delta) of the sweep
           'nOfIts'   : 1,         # number of repetitions
           'noiseOff' : 0.005,     # seconds, constant freq onset/offset
           }

    return est



if __name__ == "__main__":
    swField(plotSignal = True)



