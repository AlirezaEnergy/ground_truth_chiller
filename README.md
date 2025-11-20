# Ground Truth Chiller Model

This repository contains a detailed physics-based model of a vapor compression chiller.  
It includes sub-models for:

- supply fan  
- cooling coil  
- compressor  
- condenser and condenser fan  
- refrigerant cycle  
- chiller ON/OFF logic  

The code uses **CoolProp** for thermodynamic properties and follows assumptions from  
EnergyPlus and ASHRAE fundamentals.

---

## 📦 Installation

```bash
pip install CoolProp
```

## 📘 Documentation
Main functions

fan(...)
  Models the supply fan and calculates:
  fan power
  outlet temperature
  humidity and enthalpy

cooling_coil(...)
  Models:
  sensible and latent cooling
  bypass factor
  condensate mass flow rate
  coil outlet state

chiller_on_off_status(...)
  Checks:
  if the chiller is allowed to run
  coil load relative to minimum and maximum capacity

condenser_fan(...)
  Models:
  condenser airflow
  condenser fan power
  ground_truth_chiller(...)

Main wrapper model:
  combines coil, fan, compressor, condenser, and refrigerant cycle
  returns all outputs in a Python dictionary

## How to Cite

This software is archived on Zenodo.
If you use this model in research, please cite:
```
Ghadertootoonchi, A., & Lee, S. (2025). Ground Truth Chiller (Version 1.0.0) [Computer software]. https://doi.org/10.5281/zenodo.17664545
