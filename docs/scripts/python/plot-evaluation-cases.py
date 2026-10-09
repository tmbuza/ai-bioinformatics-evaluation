"""Generate supplied figure previews. Requires matplotlib and numpy.
Run from repository root; base-R alternative: scripts/R/plot-evaluation-cases.R.
All values are computed from the synthetic input files, not entered as results.
"""
from pathlib import Path
import csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
out=Path('results/figures');out.mkdir(parents=True,exist_ok=True)
teal,gold,ink='#036281','#f7c546','#0f172a'
plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False,'text.color':ink,'axes.labelcolor':ink,'axes.titleweight':'bold','figure.facecolor':'white'})
def read(p):
 with open(p) as f:return list(csv.DictReader(f,delimiter='\t'))
def save(fig,name,note):
 fig.text(.5,.025,note,ha='center',fontsize=9,color='#475569');fig.tight_layout(rect=(0,.065,1,1));fig.savefig(out/name,dpi=180);plt.close(fig)
m=read('cases/microbiome/data/sample-metadata.tsv');f=read('cases/microbiome/data/feature-counts.tsv');feature=next(x for x in f if x['feature_id']=='F1')
fig,ax=plt.subplots(figsize=(10,6))
for id,color in zip(sorted({x['participant_id'] for x in m}),[teal,'#577d34','#98622b','#7463a0']):
 rows=[next(x for x in m if x['participant_id']==id and x['visit']==v) for v in ['before','after']]
 vals=[100*float(feature[x['sample_id']])/sum(float(y[x['sample_id']]) for y in f) for x in rows]
 ax.plot([0,1],vals,'o-',color=color,lw=2);ax.text(1.04,vals[1],id,va='center',color=color)
ax.set(xticks=[0,1],xticklabels=['Before','After'],ylim=(0,50),xlim=(-.15,1.35),ylabel='F1 relative abundance (%)',xlabel='Visit',title='Microbiome: follow each participant');ax.grid(axis='y',alpha=.15)
save(fig,'04-microbiome-paired.png','Synthetic counts; relative abundance does not measure absolute load.')
m=read('cases/rnaseq/data/sample-metadata.tsv');a=np.array([[sum(x['condition']==c and x['batch']==b for x in m) for b in ['B1','B2']] for c in ['control','treated']])
fig,ax=plt.subplots(figsize=(10,6));from matplotlib.colors import ListedColormap
ax.imshow(a>0,cmap=ListedColormap(['#eef2f6',teal]),vmin=0,vmax=1,aspect='auto')
for i in range(2):
 for j in range(2):ax.text(j,i,f'{a[i,j]} samples',ha='center',va='center',fontsize=19,color='white' if a[i,j] else ink)
ax.set(xticks=[0,1],xticklabels=['B1','B2'],yticks=[0,1],yticklabels=['Control','Treated'],xlabel='Batch',title='RNA-seq: condition and batch are confounded')
save(fig,'05-rnaseq-design.png','No condition comparison is available within either batch. Synthetic study design.')
d=read('cases/singlecell/data/donors.tsv');cells=read('cases/singlecell/data/cells.tsv')
for x in d:
 z=[c for c in cells if c['donor_id']==x['donor_id']];x['n']=len(z);x['p']=100*sum(c['cluster']=='C1' for c in z)/len(z)
fig,(ax,bx)=plt.subplots(1,2,figsize=(11,6));bars=ax.bar([x['donor_id'] for x in d],[x['p'] for x in d],color=[teal if x['condition']=='control' else gold for x in d]);ax.bar_label(bars,labels=[f"n={x['n']}" for x in d],padding=4)
from matplotlib.patches import Patch
ax.legend(handles=[Patch(color=teal,label='Control'),Patch(color=gold,label='Treated')],loc='upper center',frameon=False,ncol=2);ax.set(ylim=(0,110),ylabel='C1 cells (%)',title='Each donor')
pooled=[];means=[]
for group in ['control','treated']:
 z=[x for x in d if x['condition']==group];pooled.append(sum(x['p']*x['n'] for x in z)/sum(x['n'] for x in z));means.append(np.mean([x['p'] for x in z]))
for offset,vals,col,lab in [(-.18,pooled,teal,'Pooled cells'),(.18,means,gold,'Donor mean')]:
 bars=bx.bar(np.arange(2)+offset,vals,.36,color=col,label=lab);bx.bar_label(bars,padding=4,fmt='%.0f')
bx.set(xticks=[0,1],xticklabels=['Control','Treated'],ylim=(0,110),ylabel='C1 cells (%)',title='Choice of summary');bx.legend(frameon=False,loc='upper center')
save(fig,'06-singlecell-weighting.png','Synthetic retained cells. Descriptive summaries do not establish a treatment effect.')
m=read('cases/proteomics/data/sample-metadata.tsv');p=read('cases/proteomics/data/log2-intensities.tsv');a=np.array([[float(x[y['sample_id']]) if x[y['sample_id']]!='NA' else np.nan for y in m] for x in p])
fig,(ax,bx)=plt.subplots(1,2,figsize=(11,6));ax.imshow(np.isnan(a),cmap=ListedColormap([teal,gold]),vmin=0,vmax=1,aspect='auto')
for i in range(a.shape[0]):
 for j in range(a.shape[1]):ax.text(j,i,'NA' if np.isnan(a[i,j]) else f'{a[i,j]:.0f}',ha='center',va='center',color=ink if np.isnan(a[i,j]) else 'white')
ax.set(xticks=range(6),xticklabels=[x['sample_id'] for x in m],yticks=range(3),yticklabels=[x['protein_id'] for x in p],xlabel='Sample',title='Observed or missing?\nNumbers are log2 intensities')
v=a[next(i for i,x in enumerate(p) if x['protein_id']=='P1')];delta=[]
for fill in [np.nan,18,22]:
 z=np.where(np.isnan(v),fill,v);delta.append(np.nanmean(z[[x['condition']=='treated' for x in m]])-np.nanmean(z[[x['condition']=='control' for x in m]]))
bx.scatter(range(3),delta,color=teal,s=75);bx.axhline(0,color='gray',ls='--',lw=1)
for i,v in enumerate(delta):bx.text(i,v+.3,f'{v:.3f}',ha='center')
bx.set(xticks=range(3),xticklabels=['Observed\nonly','Fill 18','Fill 22'],xlim=(-.4,2.4),ylim=(-4.5,3),ylabel='Treated − control (mean log2)',xlabel='Missing-value assumption',title='P1: sensitivity of direction')
save(fig,'07-proteomics-missingness.png','Synthetic protein groups. Fixed replacements illustrate sensitivity; they are not recommended imputation methods.')
g=read('cases/gwas/data/carrier-counts.tsv');ors=[]
for group in ['A','B','pooled']:
 z=[x for x in g if group=='pooled' or x['stratum']==group]
 a,b=[[sum(int(x[k]) for x in z if x['status']==s) for k in ['carriers','noncarriers']] for s in ['case','control']]
 ors.append(a[0]*b[1]/(a[1]*b[0]))
fig,ax=plt.subplots(figsize=(10,6));ax.scatter(ors,[2,1,0],s=100,color=[teal,teal,'#aa7100']);ax.axvline(1,ls='--',color='gray',lw=1)
for v,y in zip(ors,[2,1,0]):ax.text(v,y+.18,f'{v:.3f}',ha='center')
ax.set(xscale='log',xlim=(.5,8),ylim=(-.5,2.6),yticks=[2,1,0],yticklabels=['Stratum A','Stratum B','Pooled'],xlabel='Carrier odds ratio (log scale)',title='GWAS: pooling changes the association');ax.set_xticks([.5,1,2,4,8],labels=['0.5','1','2','4','8']);ax.minorticks_off()
save(fig,'08-gwas-stratification.png','Synthetic counts. Descriptive point estimates; no confidence intervals or significance tests.')
print('Generated five figures from case inputs.')
