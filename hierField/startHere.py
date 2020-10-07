# ---------------------------------------------------------
# | FM-Feedback model for pitch and FM-direction encoding |
# ---------------------------------------------------------

# This code is shared under a GNU General Public License v3.0
# License, and it's absolutely not guaranteed! If you use this for your 
# research please cite: 
# [PLACE HOLDER FOR CITATION]
# Author: Alejandro Tabas <tabas--AT--cbs--DOT--mpg--DOT--de>
# Coded and tested on python3


#  --- INSTALLATION ---
# The model requires cochlea and thorns to run:
# >>> pip3 install cochlea
# >>> pip3 install thorns 
#
# At the time of writting there is an unresolved bug in the python3 
# package for cochlea that will manifest itself when running swField
# bellow with the following error message:  
#      AttributeError: module 'numpy.fft' has no attribute 'fftpack'
# You can patch this easily. First locate where is zilany2014_rate.py 
# installed in your system. In linux you can look for this using pip:
# >>> pip3 show cochlea | grep Location
# The script you need to modify is located in:
# <PATH-TO-COCHLEA>/zilany2014/zilany2014_rate.py
# In my system this is:
# ~/.local/lib/python3.8/site-packages/cochlea/zilany2014/zilany2014_rate.py
# Open zilany2014_rate.py and comment out the lines (in my system 87 and 88):
#  if isinstance(np.fft.fftpack._fft_cache, dict):
#      np.fft.fftpack._fft_cache = {} 
# Save the file and you are done. 


#  --- GUIDED VISIT ---

# hierFiled wraps all the commands you need to run the model
import hierField

# Just for plotting...
import matplotlib.pyplot as plt

# To run the model you just need swField()
h = hierField.swField()

# The output is a dictionary with the timecourses of all dynamic variables:
print(h.keys())

# h['fre'] is the activity at the f populations
# h['fre'][:, 0] is the time course of the first population
# h['fre'].mean(0) is the activity across all cochlear channels
# h['downe'] is the activity at the down-excitatory populations
# similarly: downi, upe, upi -> down-inhibitory, up-excitatory, up-inhibitory

# The parameters and stimulus properties of the model can be accessed with:
pars = hierField.loadParameters()
est  = hierField.stimulusPars()

# pars and est are dictionaries containing all relevant paramters
print(pars.keys())
print(est.keys())

# Lets say e.g. you want to run the model with a finer temporal resolution:
pars['dt'] = 0.01
h = hierField.swField(pars)

# You can modify the stimulus parameters similary. Some important properties:
# est['nOfIts']:   number of consecutive sweeps
# est['duration']: total duration of the sweep, including the f=cte flanks
# est['freq']:     average frequency of the sweep (Hz)
# est['shift']:    frequency shift (Hz), negative for down sweeps
# est['noiseOff']: duration of the constant frequency flanks

# We use soch.py (github.com/qtabs/moch) to generate the sounds, so in theory
# any stimulus that can be created with soch can be used to run the model. 
# You can modify the stimulus type with est['type'], check out soch.py for an
# overview of the available sounds. You can also use your own soundfiles with
# the field est['filename']. We haven't tested the libraries with stimuli 
# other than pure tones and FM-sweeps, but they should work just fine. Let us
# know if the don't!


# --- EXAMPLES ---

# 1. Population responses to a down sweep
pars = hierField.loadParameters()
est = hierField.stimulusPars()
est['freq']   = 1200 
est['shift']  = -300 
est['nOfIts'] = 1    

h = hierField.swField(pars, est)

minAct = 0
maxAct = 100

f, ax = plt.subplots(3, 1, figsize=(15,8))
im = ax[0].imshow(h['fre'].T, aspect='auto', vmin = minAct, vmax = maxAct)
ax[0].set(ylabel = 'freq')
ax[1].imshow(h['upe'].T, aspect='auto', vmin = minAct, vmax = maxAct)
ax[1].set(ylabel = 'up exc')
ax[2].imshow(h['downe'].T, aspect='auto', vmin = minAct, vmax = maxAct)
ax[2].set(xlabel = 'time (ms)', ylabel = 'down exc')
f.colorbar(im, ax=ax.ravel().tolist())
f.canvas.draw()
plt.show(block=False)


# 2. Direction selectiviy in the model
par = hierField.loadParameters()
est = hierField.stimulusPars()   
est['freq'] = 1200
    
deltas = [200.0, 400.0, 600.0]
dsiUp, dsiDown = [], []

for dIx, delta in enumerate(deltas):
    est['shift'] = delta
    hUp   = hierField.swField(par, est)
    est['shift'] = -delta
    hDown = hierField.swField(par, est)     
    DSI = hUp['upe'].mean() - hDown['upe'].mean()
    DSI = DSI / (hUp['upe'].mean() + hDown['upe'].mean())
    dsiUp.append(DSI)
    DSI = hUp['downe'].mean() - hDown['downe'].mean() 
    DSI = DSI / (hUp['downe'].mean() + hDown['downe'].mean())
    dsiDown.append(DSI)

f, ax = plt.subplots(1, 2, sharex = True)
ax[0].plot(deltas, dsiUp, marker='o')
ax[0].set(title='up-sweep network', xlabel='|Delta f| (Hz)', ylabel='DSI')
ax[1].plot(deltas, dsiDown, marker='o')
ax[1].set(title='down-sweep network', xlabel='|Delta f| (Hz)', ylabel='DSI')
plt.show()


# Code used to produce all modelling figures in the paper are wrapped in 
# octopus.py. Check it out for much more many examples!