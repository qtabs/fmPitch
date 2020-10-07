import numpy as np 
import cochlea
import scipy.signal
#import matplotlib.pyplot as plt


def peripheralSpikes(sound, par, fs = -1):

    if fs == -1:
        fs = par['periphFs']

    anfTrains = cochlea.run_zilany2014(sound, fs, 
                                   anf_num = [60, 25, 15], 
                                   cf = par['cochChanns'], 
                                   species = 'human', seed = 0);

    return(anfTrains)



def peripheral(sound, par, fs = -1):

    if fs == -1:
        fs = par['periphFs']

    ANRts = cochlea.run_zilany2014_rate(sound, fs, 
                                        anf_types = ('lsr', 'msr', 'hsr'), 
                                        cf = par['cochChanns'], 
                                        species = 'human', 
                                        cohc = 1,
                                        cihc = 1)

    ANRts = .6 * ANRts['hsr'] + .25 * ANRts['msr'] + .15 * ANRts['lsr']


    if par['subCortFs'] == fs:
        p = np.asarray(ANRts) / par['subCortFs']
    else:
        resampleN = len(sound) * par['subCortFs'] / fs 
        p = scipy.signal.resample(ANRts, num = resampleN) / par['subCortFs']  
        p[p < 0] = 0
        p[p > 1] = 1

    return(p)



def subcortical(prob, lagSpace, par):

    # Processing constants
    dt = 1./par['subCortFs']
    timeSpace = np.arange(start = dt, stop = (len(prob) + 1) * dt, step = dt)
    
    if par['SACFTau'] <= 0:
        tauFactor = 2 # [Wiegriebe2000]
        taus = tauFactor * lagSpace
        taus = np.maximum(taus, 0.0025) # [Wiegriebe2000]
    else:
        taus = par['SACFTau'] * np.ones(lagSpace.shape)


    # Initalising variables
    a  = np.zeros((len(lagSpace), par['cochChanns'][2]))
    B0 = np.zeros((len(lagSpace), par['cochChanns'][2]))
    B  = np.zeros(len(lagSpace))
    z0 = np.zeros(par['cochChanns'][2])
    z  = np.zeros(len(prob))
    k  = 0.5 * np.ones(len(prob))
    C  = np.zeros((len(prob), len(lagSpace)))

    for ti in range(1, len(prob)):
        
        # SACF 
        for li in range(len(lagSpace)):
            if (timeSpace[ti - 1] - lagSpace[li] - par['solvOnset']) >  dt:
                tiL = int(max((ti - par['subCortFs'] * lagSpace[li]).round(), 1))
                a[li] = prob[ti] * prob[tiL] * par['subCortFs'] 
                
        B0 = B0 * np.exp(-dt / np.tile(taus, (1, par['cochChanns'][2])))
        B0 = B0 * dt / np.tile(taus, (1, par['cochChanns'][2])) + a 
        B0[B0 < 0] = 0
        B = B + (dt / par['subCortTau']) * (B0.sum(1) - B)

        if par['regularise'] == 1:

            # Normalisation factor (multiplicative)
            a0 = (prob[ti]**2) * par['subCortFs']
            z0 = z0 * np.exp(-dt / taus.min()) * (dt / taus.min()) + a0
            z0[z0 < 0] = 0    
            z[ti] = z[ti-1] + (dt/par['subCortTau']) * (z0.sum(0) - z[ti-1])
            
            # Normalisation factor (additive) 
            sInd = np.argmin((timeSpace[ti - 1] - 1.25 * lagSpace)**2)
            if sInd > (len(lagSpace)):
                k[ti] = 0.5
            else:
                k[ti] = B[sInd:].mean() / (z[ti] + 0.01)

            if z[ti] > 5:
                C[ti] = B / (z[ti] + 0.01)
            else:
                C[ti] = (B - 1 / (z[ti] + 0.1)) / (z[ti] + 0.01)

        else:
            C[ti] = B
            z[ti] = 0
            k[ti] = 0


    # Recomputing additive normalisation factor k and multiplicative gain A0
    if (par['regularise'] == 1):
        
        if (par['SACFGround'] < 0):
            if (len(prob) * dt < 0.075):
                print('Warning: dur(stim) < 75ms; Using default baseline 0.5')
                k0 = 0.5
            else:
                k0 = np.mean(k[int(0.05/dt) : int(min(0.10/dt, len(prob)))])            
            k[0:int(np.ceil(0.075 / dt))] = k0
            for ti in range(1, len(prob)):
                k[ti] = k[ti-1] + (dt/par['subCortTau']) * (k[ti] - k[ti-1])
        else:
            k[:] = par['SACFGround']

        kMat = np.transpose(np.tile(k, [len(lagSpace), 1]))
        A0 = par['mu0'] / np.maximum(0.1, 1 - kMat)
        C = np.maximum(0, A0 * (C - kMat))

    resampleN = int(len(prob) * par['cortFs'] / par['subCortFs'])
    A = scipy.signal.resample(C, num = resampleN, axis = 0) 
    n = scipy.signal.resample(z, num = resampleN, axis = 0) 
    b = scipy.signal.resample(k, num = resampleN, axis = 0)

    # Resampling introduces non-existent negative values
    A[A < 0] = 0
    # Resampling introduces early activity not observed in the original series
    A[0:int((par['solvOnset'] * par['cortFs']) + 1), :] = 0

    return [A, n, b]



