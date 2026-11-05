"""Search previously checked trivial intermediate backgrounds for OBC stability."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
from pathlib import Path
import json
import numpy as np
from scipy.linalg import eigvals
from bulk_stability import matrix

def main():
    data=np.load('data/generated/intermediate_topology_general/scan.npz')['records']
    rows=data[data[:,9]==0]; rates=[]
    for i,r in enumerate(rows):
        rates.append(float(eigvals(matrix(50,1,r[2],r[5],r[6],r[1]/2),check_finite=False).imag.max()))
        if i%500==0:print(i,len(rows),flush=True)
    rates=np.array(rates);out=Path('data/generated/stable_trivial');out.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(out/'search.npz',records=rows,growth=rates)
    best=np.argsort(rates)[:10]
    print(json.dumps(dict(stable=int((rates < -1e-7).sum()),best=[dict(parameters=rows[j].tolist(),growth=float(rates[j])) for j in best]),indent=2))

if __name__=='__main__':main()
