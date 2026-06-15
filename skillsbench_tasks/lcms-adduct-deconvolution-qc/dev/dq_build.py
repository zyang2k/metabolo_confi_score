#!/usr/bin/env python3
"""Rebuilt dataset: crux = data-quality dubious-detection (charge-state / multi-
alkali / notation), NOT raw-formula notation (that was the artifact).

A no-skill model defaults to "chemically parseable => ok", so it labels
[M-2H]2-, [M-3H]3-, [M+3Na-6H]3-, [M+Na+K-4H]2-, [M-2H]-, M-CO+H as ok and keeps
the flag-bins built from them. The skill teaches the QC red-flags (Rule D).
"""
import pandas as pd
PROTON=1.007276; E=0.000549
H2O=18.010565; CO2=43.989830; CO=27.994915
CL=34.968853+E; FORMATE=44.998203; ACETATE=59.013851
NA=22.989770-E; K=38.963707-E

def mz(M,ion):
    return {
        "MH":M-PROTON, "Cl":M+CL, "acetate":M+ACETATE, "formate":M+FORMATE,
        "MH_H2O":M-PROTON-H2O, "MH_CO2":M-PROTON-CO2,
        "M2H1":M-2*PROTON,              # [M-2H]- : 2 deprotonations, 1- charge -> inconsistent
        "M2H2":(M-2*PROTON)/2,          # [M-2H]2-
        "M3H3":(M-3*PROTON)/3,          # [M-3H]3-
        "NaK2":(M+NA+K-4*PROTON)/2,     # [M+Na+K-4H]2-
        "Na3_3":(M+3*NA-6*PROTON)/3,    # [M+3Na-6H]3-
        "MCOH":M-CO+PROTON,             # M-CO+H : unbracketed loss-then-add (notation error)
    }[ion]

# (compound, M, rt, gate, [(annotation, ion, category, is_isf, is_crux)])
FAM=[
 ("citric acid",192.027003,9.21,"keep",[("[M-H]-","MH","ok",0,0),("[M+CH3COO]-","acetate","ok",0,0),("[M-H-CO2]-","MH_CO2","isf",1,0)]),
 ("phenylalanine",165.078979,5.22,"keep",[("[M-H]-","MH","ok",0,0),("[M+Cl]-","Cl","ok",0,0),("[M-H-H2O]-","MH_H2O","isf",1,0)]),
 ("malonic acid",104.010959,7.50,"keep",[("[M-H]-","MH","ok",0,0),("[M-H-CO2]-","MH_CO2","isf",1,0)]),
 ("malic acid",134.021524,8.47,"keep",[("[M-H]-","MH","ok",0,0),("[M-2H]2-","M2H2","dubious",0,1),("[M-3H]3-","M3H3","dubious",0,1),("[M-H-H2O]-","MH_H2O","isf",1,0)]),
 ("glucose",180.063388,10.12,"keep",[("[M-H]-","MH","ok",0,0),("[M+HCOO]-","formate","ok",0,0),("[M-2H]2-","M2H2","dubious",0,1)]),
 ("glutaric acid",132.042259,8.01,"keep",[("[M-H]-","MH","ok",0,0),("[M+Na+K-4H]2-","NaK2","dubious",0,1),("[M-H-CO2]-","MH_CO2","isf",1,0)]),
 ("tartaric acid",150.016438,9.03,"keep",[("[M-H]-","MH","ok",0,0),("[M-2H]-","M2H1","dubious",0,1),("[M-H-H2O]-","MH_H2O","isf",1,0)]),
 ("succinic acid",118.026609,7.88,"flag",[("[M-2H]2-","M2H2","dubious",0,1),("[M-2H]-","M2H1","dubious",0,1),("[M-H-H2O]-","MH_H2O","isf",1,0)]),
 ("fumaric acid",116.010959,7.05,"flag",[("[M-3H]3-","M3H3","dubious",0,1),("[M-H-H2O]-","MH_H2O","isf",1,0)]),
 ("2-oxoglutarate",146.021524,8.90,"flag",[("[M+3Na-6H]3-","Na3_3","dubious",0,1),("[M-H-CO2]-","MH_CO2","isf",1,0)]),
 ("aconitic acid",174.016438,9.55,"flag",[("[M+Na+K-4H]2-","NaK2","dubious",0,1),("[M-H-H2O]-","MH_H2O","isf",1,0)]),
 ("gluconic acid",196.058303,11.20,"flag",[("[M-2H]-","M2H1","dubious",0,1),("M-CO+H","MCOH","dubious",0,1)]),
]
INTF=[("noise1",250.1812,4.10),("noise2",311.0827,13.40),("noise3",288.9934,3.05)]

rows,gold=[],[]; fid=1
for comp,M,rt,gate,feats in FAM:
    for j,(ann,ion,cat,isf,crux) in enumerate(feats):
        rows.append({"feature_id":f"F{fid:02d}","mz":round(mz(M,ion),4),"rt":round(rt+(j-1)*0.012,3),
                     "intensity":120000//(j+1),"adduct_annotation":ann})
        gold.append({"feature_id":f"F{fid:02d}","compound":comp,"neutral_mass":round(M,4),
                     "adduct_category":cat,"is_isf":bool(isf),"gate":gate,"is_crux":bool(crux)})
        fid+=1
for n,m,rt in INTF:
    rows.append({"feature_id":f"F{fid:02d}","mz":m,"rt":rt,"intensity":60000,"adduct_annotation":"unknown"})
    gold.append({"feature_id":f"F{fid:02d}","compound":"(interferent)","neutral_mass":None,
                 "adduct_category":"unassigned","is_isf":False,"gate":"n/a","is_crux":False})
    fid+=1

pd.DataFrame(rows).sort_values("rt").to_csv("dq_input.csv",index=False)
pd.DataFrame(gold).to_csv("dq_gold.csv",index=False)
g=pd.DataFrame(gold)
print(g[g.adduct_category!='unassigned'].groupby('adduct_category').size().to_string())
print(f"\ncrux: {int(g.is_crux.sum())} | gates keep={sum(1 for f in FAM if f[3]=='keep')} flag={sum(1 for f in FAM if f[3]=='flag')} | features {len(rows)}")
