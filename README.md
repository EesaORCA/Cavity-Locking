# Cavity-Locking

Scripts written in pyrpl 0.9.8.0 to control a STEMlab 125-14 Gen 2 Red Pitaya for:
(1) Scanning the transmission of the cavity, sweep different modulation depths for lock-in detection. 
(2) Script to analyse the data, finding the optimal modulation depth and phase.
(3) Lockbox script which scans the cavity, identifies the error signal and attempts to hold a lock. If lock is lost then will attempt to relock around last known locking point. Also has toggle-able network analyser functionality such that the resonances in the system can be identified.

Physical layout is as follows:

Laser ---> Cavity ---> Photodiode --
  |          |                       |
  |          |                        --> in1
  |          ----------Amplifier------------ out2
   -------------------------------------- out1 

Where in1 monitors the laser transmission through the cavity, out1 is used to AC modulate the laser frequency (which is then demodulated by the IQ module in the RP) and out2 is used to both scan the piezo and feedback on the cavity in locking mode. 

All scripts written in VScode and meant to be used with interactive window. 

Workflow starting from scratch with an arbitrary cavity is:
(1) Tune laser onto resonance with the cavity, sweep across cavity resonances using bench top equipment.
(2) Park the laser and do the same scans with cavities.
(3) Unplug bench-top function generator and scope, switch to Redpitaya and do initial scans with Lockinscript.py to identify the differential lock signal and check that input voltage ranges/ modulation depths and phase are broadly acceptable.
(3a) Optimal is perform various scans and feed into the ANALYSIS script to find optimal modulation depth and phase.
(4) Run LOCKBOX.py and toggle network analyser on to identify resonances and recommend a particular integral gain setting. 
(5) C

##Lockinscript.py and LockinscriptPIDscan.py 

Both are used for scanning the cavity and performing lock-in detection across a range of modulation depths. Difference is how the piezo is scanned and subsequently measure. For Lockinscript.py the piezo is modulated by an asg ramp whilst a scope trace reads for the duration of the ramp. This is a fast measurement a la scop (as fast as you set the ramp to be) 

