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
