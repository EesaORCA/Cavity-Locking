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

##Lockinscript.py and LockinscriptPIDscan.py 

Both are used for scanning the cavity and performing lock-in detection across a range of modulation depths. Difference is how the piezo is scanned and subsequently measure. For Lockinscript.py

