"""Build version-consistent manuscript tables from frozen inputs and model outputs."""
from __future__ import annotations
import argparse, ast, hashlib, json, re
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

PATHWAYS = {
 'ISC': ['iscR','iscS','iscU','iscA','hscB','hscA','fdx'],
 'SUF': ['sufA','sufB','sufC','sufD','sufS','sufE'],
 'Cysteine': ['cysK','cysM','cysE','cysA','cysW','cysU','cysJ','cysI','cysH','cysD','cysN','cysC','cysP','sbp','cysB'],
 'Glutathione': ['gshA','gshB','gor','ggt','gloA','gloB','gloC','yajG','grxA','grxB','grxC','grxD','trxA','trxB','trxC'],
 'Electrophile': ['nemA','nemR','yqhD','yqhC','dkgA','dkgB','ahr','ydgI','yafC','yqhE'],
 'FeS/iron': ['fur','ryhB','sodB','fumA','acnA','acnB','sdhC','sdhD','sdhA','sdhB','bfr','ftnA','ftnB','nfuA','erpA','yggX'],
 'Other sulfur': ['tauA','tauB','tauC','tauD','ssuA','ssuB','ssuC','ssuD','ssuE','metA','metB','metC','metE','metH','metF','metK','metL'],
}

def network_tables(data: Path):
    raw = pd.read_excel(data/'clash/elife-54655-supp3-v2.xlsx', header=None).iloc[9:].copy()
    raw.columns=['idx','sRNA','mRNA','sRNA_frag','mRNA_frag','total_hybrids','min_MFE','in_RILseq','Exp','Transition','Stationary']
    raw['source_excel_row'] = raw.index + 1
    raw=raw.dropna(subset=['sRNA','mRNA'])
    for col in ['sRNA','mRNA']: raw[col]=raw[col].astype(str).str.strip()
    raw['total_hybrids']=pd.to_numeric(raw.total_hybrids,errors='raise')
    raw['min_MFE']=pd.to_numeric(raw.min_MFE,errors='coerce')
    edges=[]
    for (s,t),g in raw.groupby(['sRNA','mRNA'],sort=True):
        count=int(g.total_hybrids.sum()); ril=g.in_RILseq.eq('Yes').any()
        phases=sum(any(g[col].notna() & ~g[col].astype(str).isin(['No','0','0.0','nan','-'])) for col in ['Exp','Transition','Stationary'])
        edges.append(dict(sRNA=s,mRNA=t,total_hybrids=count,best_MFE=g.min_MFE.min(),in_RILseq='Yes' if ril else 'No',
            n_phases=phases,source_excel_rows=';'.join(g.source_excel_row.astype(str)),
            evidence_tier='high' if ril or count>=5 else 'standard' if count>=2 else 'exploratory'))
    edges=pd.DataFrame(edges)
    sulfur=set(sum(PATHWAYS.values(),[])); fes=set(PATHWAYS['ISC']+PATHWAYS['SUF']+PATHWAYS['FeS/iron'])
    assert len(sulfur)==86 and len(fes)==29
    t1=edges[edges.mRNA.isin(sulfur)].copy()
    t1['target_pathway']=t1.mRNA.map({g:k for k,gs in PATHWAYS.items() for g in gs})
    assert len(t1)==52
    # Expression proxy is explicitly a separate study, not same-study matched expression.
    x=[]
    for p in [data/'GSM5024438_8h-WT-minus-1.txt.gz',data/'GSM5024439_8h-WT-minus-2.txt.gz']:
        z=pd.read_csv(p,sep='\t',header=None,names=['tx','count']);z.tx=z.tx.str.replace('transcript:','',regex=False)
        x.append(z.set_index('tx')['count'])
    counts=pd.concat(x,axis=1).sum(axis=1)
    ann=pd.read_csv(data/'GSE126176_counts_genes_formatted_metadata.txt.gz',sep='\t')
    proxy=pd.DataFrame({'tx':counts.index,'count':counts.values}).merge(ann[['Transcript','GeneName']],left_on='tx',right_on='Transcript')
    logcpm=np.log2((proxy['count']+.5)/proxy['count'].sum()*1e6)
    ge=pd.Series(logcpm.values,index=proxy.GeneName).groupby(level=0).max()
    detectable=set(ge[ge>1].index)
    targets=set(edges.mRNA)
    records=[]
    def add(ts,univ,bg,label,subsystem):
        hit=len(ts&univ); n_t=len(ts); n_u=len(univ)
        tab=[[hit,n_t-hit],[n_u-hit,bg-n_t-n_u+hit]]
        assert min(sum(tab,[]))>=0
        odds,p=fisher_exact(tab)
        se=np.sqrt(sum(1/c for c in sum(tab,[]))) if min(sum(tab,[]))>0 else np.nan
        records.append(dict(analysis=label,subsystem=subsystem,universe_size=n_u,hit=hit,n_targets=n_t,background=bg,
            OR=odds,CI95_low=np.exp(np.log(odds)-1.96*se),CI95_high=np.exp(np.log(odds)+1.96*se),P=p,contingency=json.dumps(tab)))
    for label,ts,bg in [('Primary',targets,4400),('Multi-chimera pairs',set(edges.loc[edges.total_hybrids>=2,'mRNA']),4400),
                        ('RIL-supported pairs',set(edges.loc[edges.in_RILseq.eq('Yes'),'mRNA']),4400),('Expression-proxy background',targets&detectable,len(detectable))]:
        for name,u in [('Sulfur network',sulfur),('FeS/iron',fes)]:
            add(ts,u&detectable if label=='Expression-proxy background' else u,bg,label,name)
    t5=pd.DataFrame(records);t5['BH_q_eight_Fisher_tests']=multipletests(t5.P,method='fdr_bh')[1]
    det=sorted(detectable); bins=pd.qcut(ge[det],10,labels=False)
    pools={i:np.array([g for g in det if bins[g]==i]) for i in range(10)}
    need=pd.Series([bins[g] for g in sulfur&detectable]).value_counts().sort_index()
    rng=np.random.default_rng(42); observed=len(sulfur&detectable&targets)
    null=np.array([len(set(np.concatenate([rng.choice(pools[int(b)],int(n),replace=False) for b,n in need.items()]))&targets) for _ in range(10000)])
    extra=dict(analysis='Expression-decile permutation',subsystem='Sulfur network',universe_size=len(sulfur&detectable),hit=observed,
               n_targets=len(targets&detectable),background=len(detectable),P=(1+(null>=observed).sum())/10001,
               contingency=f'mean_null_overlap={null.mean():.4f};n_permutations=10000;seed=42')
    t5=pd.concat([t5,pd.DataFrame([extra])],ignore_index=True)
    return t1,t5


def conservation_tables(compgen: Path):
    align=json.loads((compgen/'alignments_v3.json').read_text())
    legacy=pd.read_csv(compgen/'Table_S3_conservation_v3.csv').set_index('genome')
    rows=[]; positions=[]
    for genome,a in sorted(align.items()):
        g=genome.replace('UMN026_APEC','UMN026_UPEC')
        species='E. coli' if g.startswith('E.coli') else 'Salmonella'
        metrics={}
        for locus,qlen in [('iscS',210),('ryhB',95)]:
            qa,sa=a[locus+'_q'],a[locus+'_s'];assert len(qa)==len(sa)
            qidx=-1; first=next(i for i,c in enumerate(qa) if c!='-');last=max(i for i,c in enumerate(qa) if c!='-')
            for qc,sc in zip(qa,sa):
                if qc=='-': continue
                qidx+=1
                positions.append(dict(genome=g,species=species,locus=locus,query_index=qidx,relative_ATG=qidx-150 if locus=='iscS' else None,
                    reference_base=qc,query_base=sc,match=int(qc==sc and sc!='-'),gap=int(sc=='-')))
            assert qidx+1==qlen
            qcore,score=qa[first:last+1],sa[first:last+1]
            metrics[locus+'_aligned_span_columns']=len(qcore)
            metrics[locus+'_insertion_columns']=sum(qc=='-' for qc in qcore)
        records=pd.DataFrame([r for r in positions if r['genome']==g])
        row=dict(genome=g,species=species,**metrics)
        for locus,reg,sub in [('ryhB','query95',records[records.locus.eq('ryhB')]),
                              ('iscS','query210',records[records.locus.eq('iscS')]),
                              ('iscS','pairing_window',records[records.locus.eq('iscS') & records.query_index.between(130,148)]),
                              ('iscS','coding60',records[records.locus.eq('iscS') & records.query_index.ge(150)])]:
            tag=locus+'_'+reg
            row.update({tag+'_n_positions':len(sub),tag+'_n_matches':int(sub.match.sum()),tag+'_n_deletions':int(sub.gap.sum()),
                tag+'_identity_gaps_included':100*sub.match.sum()/len(sub),tag+'_identity_gap_excluded':100*sub.match.sum()/(len(sub)-sub.gap.sum())})
        lr=legacy.loc[g if g in legacy.index else genome]
        row['ryhB_hit_cluster_count']=int(lr.ryhB_n_loci)
        row['ryhB_contig']=lr.ryhB_contig;row['ryhB_HSP_start']=int(lr.ryhB_smin);row['ryhB_HSP_end']=int(lr.ryhB_smax)
        row['interpretation']='query-anchored similarity; hit clusters are not proof of functional orthologs'
        rows.append(row)
    pos=pd.DataFrame(positions)
    window=pos[pos.locus.eq('iscS') & pos.query_index.between(130,148)]
    summary=[]
    for species,z in [('All',window),('E. coli',window[window.species.eq('E. coli')]),('Salmonella',window[window.species.eq('Salmonella')])]:
        summary.append(dict(group=species,n_genomes=z.genome.nunique(),n_positions=len(z),matches=int(z.match.sum()),deletions=int(z.gap.sum()),
            identity_gap_inclusive=100*z.match.sum()/len(z),identity_nongap=100*z.match.sum()/(len(z)-z.gap.sum())))
    return pd.DataFrame(rows),pos,pd.DataFrame(summary)


def build(data, stats, compgen, output):
    data,stats,compgen,output=map(Path,[data,stats,compgen,output]);supp=output/'supplementary';supp.mkdir(parents=True,exist_ok=True)
    tab_dir=output/'audit';tab_dir.mkdir(exist_ok=True)
    t1,t5=network_tables(data)
    t1.to_csv(supp/'Table_S1_CLASH_interactions_v4.csv',index=False)
    t5.to_csv(supp/'Table_S5_network_sensitivity_v4.csv',index=False)
    t3,pos,cs=conservation_tables(compgen)
    t3.to_csv(supp/'Table_S3_conservation_v4.csv',index=False)
    pos.to_csv(tab_dir/'conservation_positions_v4.csv',index=False)
    cs.to_csv(tab_dir/'pairing_window_summary_v4.csv',index=False)
    de={c:pd.read_csv(stats/f'DE_{c}.csv') for c in ['Reuterin','HOCl','Ferrate','Hpx']}
    de['Reuterin'].to_csv(supp/'Table_S4_Reuterin_DE_v4.csv',index=False)
    p=pd.read_csv(stats/'genesets_primary.csv')
    second=pd.read_csv(stats/'genesets_seed_check.csv')
    p=p.merge(second[['contrast','gene_set','roast_q_up']].rename(columns={'roast_q_up':'second_seed_q_up'}),on=['contrast','gene_set'])
    p['seed_threshold_stable']=(p.roast_q_up<.05)==(p.second_seed_q_up<.05)
    all_tables=[]
    eligibility=pd.read_csv(stats/'array_rotation_eligibility.csv')
    allowed_reu=set(eligibility.loc[eligibility.eligible_rotation_primary,'symbol'].dropna())
    fc={c:d.dropna(subset=['symbol','logFC']).groupby('symbol').logFC.mean() for c,d in de.items()}
    rng=np.random.default_rng(42)
    permp=[]
    for row in p.itertuples():
        pool=fc[row.contrast]
        if row.contrast=='Reuterin': pool=pool[pool.index.isin(allowed_reu)]
        if row.contrast=='Hpx': pool=pool.drop(['katE','katG','ahpC'],errors='ignore')
        vals=pool.to_numpy();obs=row.raw_mean_logFC
        null=np.array([rng.choice(vals,row.n_genes,replace=False).mean() for _ in range(10000)])
        permp.append((1+(np.abs(null)>=abs(obs)).sum())/10001)
        expected=fc[row.contrast].loc[row.used_genes.split(';')].mean()
        assert np.isclose(obs,expected,atol=1e-10)
    p['permutation_P']=permp;p['permutation_BH_q']=multipletests(permp,method='fdr_bh')[1];p['analysis']='primary'
    all_tables.append(p)
    second['analysis']='second_seed';all_tables.append(second)
    for label in ['array_unweighted','array_masked_probe','rnaseq_trend']:
        z=pd.read_csv(stats/f'genesets_sensitivity_{label}.csv');z['analysis']=label;all_tables.append(z)
    t2=pd.concat(all_tables,ignore_index=True)
    t2.to_csv(supp/'Table_S2_gene_set_statistics_v4.csv',index=False)
    # Resample genes only, preserving paired cross-study membership. Descriptive,
    # not biological confidence intervals and not evidence for chemical causality.
    bootstrap=[]
    for sn in p.gene_set.unique():
        ra=p[(p.contrast=='Reuterin') & p.gene_set.eq(sn)].iloc[0]
        for c in ['HOCl','Ferrate','Hpx']:
            rb=p[p.contrast.eq(c) & p.gene_set.eq(sn)].iloc[0]
            common=sorted(set(ra.used_genes.split(';'))&set(rb.used_genes.split(';')))
            if len(common)<3: continue
            delta=fc['Reuterin'].loc[common].values-fc[c].loc[common].values
            boots=rng.choice(delta,size=(10000,len(delta)),replace=True).mean(axis=1)
            lo,hi=np.quantile(boots,[.025,.975])
            bootstrap.append(dict(gene_set=sn,comparison='Reuterin vs '+c,n_common_genes=len(common),common_genes=';'.join(common),
                mean_logFC_difference=delta.mean(),resampling_2_5_percentile=lo,resampling_97_5_percentile=hi,
                resampling_unit='gene (not biological sample)',interpretation='descriptive membership sensitivity, not a cross-study causal test',n_resamples=10000,seed=42))
    pd.DataFrame(bootstrap).to_csv(supp/'Table_S6_gene_resampling_v4.csv',index=False)
    pd.read_csv(stats/'batch_sensitivity.csv').to_csv(supp/'Table_S7_block_sensitivity_v4.csv',index=False)
    old_map={'ISC operon':'iscRSUA operon','SUF operon':'suf operon','Cysteine regulon':'Cys regulon','Glutathione/glyoxalase':'GSH/glyoxalase','Electrophile detoxification':'Electrophile detox'}
    details={'pairing_window':cs.to_dict(orient='records'),'network':{'n_edges':len(t1),'n_sRNAs':t1.sRNA.nunique(),'n_targets':t1.mRNA.nunique(),'tiers':t1.evidence_tier.value_counts().to_dict()},
        'seed_flips':int((~p.seed_threshold_stable).sum()),'n_primary_gene_set_tests':len(p),
        'array_rotation_eligible_genes':int(eligibility.eligible_rotation_primary.sum()),
        'principal_model':'native limma weighted array DE; complete-positive-weight array roast; TMM/voom RNA-seq',
        'genotype_exclusion':'Hpx: katE, katG, ahpC, from gene-set scoring only; raw counts and DE preserved'}
    (tab_dir/'result_summary_v4.json').write_text(json.dumps(details,indent=2,ensure_ascii=False)+'\n')
    return details

if __name__=='__main__':
    q=argparse.ArgumentParser();q.add_argument('--data',required=True);q.add_argument('--stats',required=True);q.add_argument('--compgen',required=True);q.add_argument('--out',required=True)
    a=q.parse_args();print(json.dumps(build(a.data,a.stats,a.compgen,a.out),indent=2,ensure_ascii=False))
