import numpy as np
import sys
import scipy.io
import matplotlib.pyplot as plt
sys.path.append('./hierField')
import hierField
import pickle
import time


def calculateSpectralModel(): # Predictions of the spectral model (Fig 3)

    experiment = scipy.io.loadmat('expResAvgSw.mat')
    percFr = experiment['fAvg']

    f0 = [900., 1200., 1500.]
    deltas = np.linspace(-600, 600, 10)

    est = hierField.stimulusPars()
    est['type'] = 'FMsweep'
    est['duration'] = 0.05   # This is the total duration, including noiseOff
    est['noiseOff'] = 0.005  # This is the duration of the f=cte segments
    est['nOfIts']   = 1      # This is the number of sweeps

    par = hierField.loadParameters()

    p = dict()
    for fIx, freq in enumerate(f0):
        print('F = {}'.format(freq), end = ' ')
        for dIx, delta in enumerate(deltas):
            print('.', end = ' ')
            est['freq'] = freq
            est['shift'] = delta
            p['f{:.0f}d{:.0f}'.format(freq, delta)] = hierField.peripheral(par, est)
        print('')

    expFr = dict()
    for fIx, freq in enumerate(f0):
        print('F = {}'.format(freq), end = ' ')
        for dIx, delta in enumerate(deltas):
            print('.', end = ' ')
            est['freq'] = percFr[fIx, dIx]
            est['shift'] = 0
            expFr['f{:.0f}d{:.0f}'.format(freq, delta)] = hierField.peripheral(par, est)
        print('')

    scipy.io.savemat('spectSims.mat', {'p': p, 'expFr': expFr})



def exampleResponses(): # Example of the population responses (Fig 6)

    est = hierField.stimulusPars()
    est['onset'] = 0.02 # <- we want to see the onset/offset in the plots!
    est['tail']  = 0.04 # 
    est['freq']  = 1200

    est['shift'] = +300
    drop = hierField.swField(est = est, saveMat = 'popresUp.mat')

    est['shift'] = -300
    drop = hierField.swField(est = est, saveMat = 'popresDown.mat')



def directionSelectivity(): # DSI of the model (fig 7)

    deltas = np.linspace(-600, 600, 10)
    f0     = np.array([900, 1200, 1500])

    est = hierField.stimulusPars()
    est['type'] = 'FMsweep'
    est['duration'] = 0.05
    est['noiseOff'] = 0.005
    est['nOfIts']   = 1

    downExc, upExc = np.zeros([3, 10]), np.zeros([3, 10])
    dsiDown, dsiUp = np.zeros([3, 5]), np.zeros([3, 5])
    for fIx, freq in enumerate(f0):
        for dIx, delta in enumerate(deltas):
            est['shift'] = delta
            est['freq']  = freq
            h = hierField.swField(est = est)
            upExc[fIx, dIx]   = h['upe'].mean()
            downExc[fIx, dIx] = h['downe'].mean()
        for dIx in range(5):
            d = downExc[fIx,9-dIx] - downExc[fIx,dIx]
            dsiDown[fIx,dIx] = d / (downExc[fIx,9-dIx] + downExc[fIx,dIx])
            d = upExc[fIx,9-dIx] - upExc[fIx,dIx]
            dsiUp[fIx,dIx] 

    # Data storing
    matDict = dict()
    matDict['topDownOn'] = {'downExc' : downExc,
                            'upExc'   : upExc,
                            'dsiDown' : dsiDown,
                            'dsiUp'   : dsiUp}


    deltas = np.linspace(-600, 600, 10)
    f0     = np.array([900, 1200, 1500])

    est = hierField.stimulusPars()
    est['type'] = 'FMsweep'
    est['duration'] = 0.05
    est['noiseOff'] = 0.005
    est['nOfIts']   = 1

    par = hierField.loadParameters()
    par['JsN'] = 0

    downExc, upExc = np.zeros([3, 10]), np.zeros([3, 10])
    dsiDown, dsiUp = np.zeros([3, 5]), np.zeros([3, 5])
    for fIx, freq in enumerate(f0):
        for dIx, delta in enumerate(deltas):
            est['shift'] = delta
            est['freq']  = freq
            h = hierField.swField(par, est)
            upExc[fIx, dIx]   = h['upe'].mean()
            downExc[fIx, dIx] = h['downe'].mean()
        for dIx in range(5):
            d = downExc[fIx,9-dIx] - downExc[fIx,dIx]
            dsiDown[fIx,dIx] = d / (downExc[fIx,9-dIx] + downExc[fIx,dIx])
            d = upExc[fIx,9-dIx] - upExc[fIx,dIx]
            dsiUp[fIx,dIx] 

    # Data storing
    matDict['topDownOff'] = {'downExc' : downExc,
                             'upExc'   : upExc,
                             'dsiDown' : dsiDown,
                             'dsiUp'   : dsiUp}

    scipy.io.savemat('dsi.mat', matDict)



def calculateFMPitchModel(): # Predictions of the FM-feedback model for sweeps

    experiment = scipy.io.loadmat('expResAvgSw.mat')
    percFr = experiment['fAvg']

    f0 = [900., 1200., 1500.]
    deltas = np.linspace(-600, 600, 10)
    sSpace = np.arange(100)

    est = hierField.stimulusPars()
    est['type'] = 'FMsweep'
    est['duration'] = 0.05   # This is the total duration, including noiseOff
    est['noiseOff'] = 0.005  # This is the duration of the f=cte segments
    est['nOfIts']   = 1      # This is the number of sweeps

    par = hierField.loadParameters()

    sp = dict()
    sPeaks = np.zeros([len(f0), len(deltas)])
    for fIx, freq in enumerate(f0):
        for dIx, delta in enumerate(deltas):
            est['freq'] = freq
            est['shift'] = delta
            h = hierField.swField(par, est)
            sp['f{:.0f}d{:.0f}'.format(freq, delta)] = h['fre']
            c = h['fre'].mean(0)
            fDist = np.exp(c) / np.exp(c).sum()
            sPeaks[fIx, dIx] = sum(sSpace * fDist);

    ep = dict()
    ePeaks = np.zeros([len(f0), len(deltas)])
    for fIx, freq in enumerate(f0):
        for dIx, delta in enumerate(deltas):
            est['freq'] = percFr[fIx, dIx]
            est['shift'] = 0
            h = hierField.swField(par, est)
            ep['f{:.0f}d{:.0f}'.format(freq, delta)] = h['fre']
            c = h['fre'].mean(0)
            fDist = np.exp(c) / np.exp(c).sum()
            ePeaks[fIx, dIx] = sum(sSpace * fDist);

    asymm = np.array([[sPeaks[fIx][9-i] - sPeaks[fIx][i] 
                                    for fIx in range(3)] for i in range(5)])

    rr = {'sp':sp, 'ep':ep, 'sPeaks':sPeaks, 'ePeaks':ePeaks, 'asymm':asymm}
    scipy.io.savemat('pitchPredictionsSw.mat', rr)



def exploreParameterSpace(): # Predictions of the spectral model (Fig 3)

    def explore(par):

        f0 = [900., 1200., 1500.]
        deltas = np.linspace(-600, 600, 10)
        sSpace = np.arange(100)

        est = hierField.stimulusPars()
        est['type'] = 'FMsweep'
        est['duration'] = 0.05   # This is the total duration, including noiseOff
        est['noiseOff'] = 0.005  # This is the duration of the f=cte segments
        est['nOfIts']   = 1      # This is the number of sweeps

        sweepPeak = np.zeros([len(f0), len(deltas)])
        for fIx, freq in enumerate(f0):
            for dIx, delta in enumerate(deltas):
                est['freq'] = freq
                est['shift'] = delta
                h = hierField.swField(par, est)
                c = h['fre'].mean(0)
                fDist = np.exp(c) / np.exp(c).sum()
                sweepPeak[fIx, dIx] = sum(sSpace * fDist);

        experiment = scipy.io.loadmat('expResAvgSw.mat')
        percFr = experiment['fAvg']
        expPeak = np.zeros([len(f0), len(deltas)])
        for fIx, freq in enumerate(f0):
            for dIx, delta in enumerate(deltas):
                est['freq'] = percFr[fIx, dIx]
                est['shift'] = 0
                h = hierField.swField(par, est)
                c = h['fre'].mean(0)
                fDist = np.exp(c) / np.exp(c).sum()
                expPeak[fIx, dIx] = sum(sSpace * fDist);

        mo = np.reshape(sweepPeak, len(deltas)*len(f0))
        ex = np.reshape(expPeak, len(deltas)*len(f0))
        
        return(1 - sum((mo-ex)**2) / sum((mo-mo.mean())**2) )


    results = []

    # 1. Jn VS tau
    TAU = np.arange(2., 22., 2.)
    JN  = np.arange(0., 0.105, 0.005)

    R2 = np.nan * np.ones((len(TAU), len(JN))) 
    for tIx, tau in enumerate(TAU):
        tt = time.time()
        print("it {}/{}".format(tIx, len(TAU)), end = ' ')
        for jIx, jn in enumerate(JN):
            par = hierField.loadParameters()
            par['tauPope'], par['JsN'] = tau, jn
            R2[tIx, jIx] = explore(par)
        print("done! Time: {:.0f}".format(time.time()-tt))

    results.append(dict([('contrast', ['tauPope VS JsN']),
                         ('par1', TAU),
                         ('par2', JN), 
                         ('R2'  , R2)]))
    
    # 2. Jn VS tau with constant Tau
    TAU = np.arange(2., 22., 2.)
    JN  = np.arange(0., 0.105, 0.005)

    R2 = np.nan * np.ones((len(TAU), len(JN))) 
    for tIx, tau in enumerate(TAU):
        tt = time.time()
        print("it {}/{}".format(tIx, len(TAU)), end = ' ')
        for jIx, jn in enumerate(JN):
            par = hierField.loadParameters()
            par['tauEff'] = False
            par['tauPope'], par['JsN'] = tau, jn
            R2[tIx, jIx] = explore(par)
        print("done! Time: {:.0f}".format(time.time()-tt))

    results.append(dict([('contrast', ['tauPope VS JsN (Tau = cte)']),
                         ('par1', TAU),
                         ('par2', JN), 
                         ('R2'  , R2)]))
    with open('cache.pickle', 'wb') as handle:
        pickle.dump(results, handle)

    # 3. wFeed VS dFeed with constant Tau
    WF  = np.arange(0.01, 0.11, 0.01)
    DF  = np.arange(0.01, 0.11, 0.01)

    R2 = np.nan * np.ones((len(WF), len(DF))) 
    for wIx, wf in enumerate(WF):
        tt = time.time()
        print("it {}/{}".format(wIx, len(WF)), end = ' ')
        for dIx, df in enumerate(DF):
            par = hierField.loadParameters()
            par['wFeed'], par['dFeed'] = wf, df
            R2[wIx, dIx] = explore(par)
        print("done! Time: {:.0f}".format(time.time()-tt))

    results.append(dict([('contrast', ['wFeed VS dFeed']),
                         ('par1', WF),
                         ('par2', DF), 
                         ('R2'  , R2)]))

    scipy.io.savemat('parameterSpace.mat', {'r2': results})



def calculateBradys():

    par = hierField.loadParameters()

    # Brady expII -- different onsets
    expFreqII = np.array([[1507, 1517, 1511, 1499, 1230, 1029],
                        [ 969,  981,  973, 1006, 1455, 1481]])

    onsetsII  = np.arange(0.010, 0.070, 0.010)

    est = hierField.stimulusPars()
    est['type']     = 'FMsweep'
    est['nOfIts']   = 1
    est['duration'] = 0.1

    bradyIImod = [[0] * len(onsetsII) for ix in range(2)]
    bradyIIExp = [[0] * len(onsetsII) for ix in range(2)]
    for jx in range(2):
        for ix in range(len(onsetsII)):
            est['freq']     = 1250
            est['noiseOff'] = onsetsII[ix]
            est['transOff'] = 0.1 - onsetsII[ix] - 0.02
            est['shift']    = (-1)**jx * 250
            bradyIImod[jx][ix] = hierField.swField(par, est)['fre']
            est['freq']     = expFreqII[jx, ix]
            est['shift']    = 0
            bradyIIExp[jx][ix] = hierField.swField(par, est)['fre']


    # Brady expIII -- different durations with the same offset
    expFreqIII    = np.array([[ 1458, 1460, 1427, 1194, 1036],
                            [ 1149, 1183, 1240, 1466, 1477]])

    onsetsIII  = np.arange(0.010, 0.060, 0.010)

    est = hierField.stimulusPars()
    est['type']     = 'FMsweep'
    est['nOfIts']   = 1
    est['duration'] = 0.08
    est['transOff'] = 0.01

    bradyIIImod = [[0] * len(onsetsIII) for ix in range(2)]
    bradyIIIExp = [[0] * len(onsetsIII) for ix in range(2)]
    for jx in range(2):
        for ix in range(len(onsetsIII)):
            est['freq']     = 1250
            est['noiseOff'] = onsetsIII[ix]
            est['shift']    = (-1)**jx * 250
            bradyIIImod[jx][ix] = hierField.swField(par, est)['fre']
            est['freq']     = expFreqIII[jx, ix]
            est['shift']    = 0
            bradyIIIExp[jx][ix] = hierField.swField(par, est)['fre']

    # Data storing
    matDict = {'onsetsII'    : onsetsII,
               'onsetsIII'   : onsetsIII,
               'bradyIImod'  : bradyIImod,
               'bradyIIImod' : bradyIIImod,
               'bradyIIExp'  : bradyIIExp,
               'bradyIIIExp' : bradyIIIExp}

    scipy.io.savemat('bradys.mat', matDict)



def calculateFMPitchModelTrains():  # Predictions for sweep trains

    par = hierField.loadParameters()
    
    experiment = scipy.io.loadmat('expResAvgTr.mat')
    percFr = experiment['fAvg']

    f0 = [900., 1200., 1500.]
    deltas = np.linspace(-333, 333, 6)
    sSpace = np.arange(100)

    est = hierField.stimulusPars()
    est['type'] = 'FMsweep'
    est['duration'] = 0.05   # This is the total duration, including noiseOff
    est['noiseOff'] = 0.005  # This is the duration of the f=cte segments
    est['nOfIts']   = 5      # This is the number of sweeps

    sp = dict()
    sPeaks = np.zeros([len(f0), len(deltas)])
    for fIx, freq in enumerate(f0):
        for dIx, delta in enumerate(deltas):
            est['freq'] = freq
            est['shift'] = delta
            h = hierField.swField(par, est)
            sp['f{:.0f}d{:.0f}'.format(freq, delta)] = h['fre']
            c = h['fre'].mean(0)
            fDist = np.exp(c) / np.exp(c).sum()
            sPeaks[fIx, dIx] = sum(sSpace * fDist);

    ep = dict()
    ePeaks = np.zeros([len(f0), len(deltas)])
    for fIx, freq in enumerate(f0):
        for dIx, delta in enumerate(deltas):
            est['freq'] = percFr[fIx, dIx]
            est['shift'] = 0
            h = hierField.swField(par, est)
            ep['f{:.0f}d{:.0f}'.format(freq, delta)] = h['fre']
            c = h['fre'].mean(0)
            fDist = np.exp(c) / np.exp(c).sum()
            ePeaks[fIx, dIx] = sum(sSpace * fDist);

    asymm = np.array([[sPeaks[fIx][5-i] - sPeaks[fIx][i] 
                                    for fIx in range(3)] for i in range(3)])

    rr = {'sp':sp, 'ep':ep, 'sPeaks':sPeaks, 'ePeaks':ePeaks, 'asymm':asymm}
    scipy.io.savemat('pitchPredictionsTr.mat', rr)



def connectivityMatrices(): # Connectivity matrices
    
    w = hierField.connectivityMatrices()

    scipy.io.savemat('connectivity.mat', w)



tt = time.time()

t0 = time.time()
print('example responses', end = ' ')
exampleResponses()
print("done! Time {:.0f}seconds".format(time.time() - t0))

t0 = time.time()
print('DSI', end = ' ')
directionSelectivity()
print("done! Time {:.0f}seconds".format(time.time() - t0))

t0 = time.time()
print('Sweeps', end = ' ')
calculateFMPitchModel()
print("done! Time {:.0f}seconds".format(time.time() - t0))

t0 = time.time()
print('R2 across parametrizations', end = ' ')
exploreParameterSpace()
print("done! Time {:.0f}seconds".format(time.time() - t0))

t0 = time.time()
print('Bradys', end = ' ')
calculateBradys()
print("done! Time {:.0f}seconds".format(time.time() - t0))

t0 = time.time()
print('Trains', end = ' ')
calculateFMPitchModelTrains()
print("done! Time {:.0f}seconds".format(time.time() - t0))

print('Connectivity matrices')
connectivityMatrices()

print("Total time: {:.0f}seconds".format(time.time() - tt))
