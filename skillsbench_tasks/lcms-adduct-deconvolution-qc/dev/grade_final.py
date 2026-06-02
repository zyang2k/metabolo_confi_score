import pandas as pd, io
gold = pd.read_csv('final_gold.csv')
pred = pd.read_csv('pred_final.csv')
m = gold.merge(pred, on='feature_id', suffixes=('_g','_p'))
# deconvolution (informational)
def nm(a,b):
    if pd.isna(a) and pd.isna(b): return True
    if pd.isna(a) or pd.isna(b): return False
    return abs(a-b)<0.02
m['nm_ok']=[nm(a,b) for a,b in zip(m.neutral_mass_g,m.neutral_mass_p)]
m['cat_ok']=m.adduct_category_g.str.strip()==m.adduct_category_p.str.strip()
# SCORED: crux category + gate
crux = m[m.is_crux]
gates_g = gold[gold.gate!='n/a'][['compound','gate']].drop_duplicates()
gates_p = m[['compound','gate_p']].drop_duplicates(subset='compound')
gm = gates_g.merge(gates_p, on='compound'); gm['ok']=gm.gate==gm.gate_p
crux_pass, gate_pass = int(crux.cat_ok.sum()), int(gm.ok.sum())
total = len(crux)+len(gm)
print('PLUMBING (not scored): neutral_mass %d/%d = %.0f%%'%(m.nm_ok.sum(),len(m),100*m.nm_ok.mean()))
print('SCORED  crux category: %d/%d'%(crux_pass,len(crux)))
print('SCORED  bin gate:      %d/%d'%(gate_pass,len(gm)))
print('==> NO-SKILL REWARD (crux+gate): %d/%d = %.1f%%'%(crux_pass+gate_pass,total,100*(crux_pass+gate_pass)/total))
print('\ncrux misses:')
print(crux[~crux.cat_ok].merge(pd.read_csv('final_input.csv')[['feature_id','adduct_annotation']],on='feature_id')[['feature_id','adduct_annotation','adduct_category_g','adduct_category_p']].to_string(index=False))
print('\ngate misses:'); print(gm[~gm.ok].to_string(index=False))
