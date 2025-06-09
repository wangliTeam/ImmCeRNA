import os 
import pandas as pd 
import numpy as np 
import pyreadr 
from config import RNR_DIR, CERNAS_PATH 
 
def load_rdata_file(rdata_path):
    """Load single RData file"""
    try:
        result = pyreadr.read_r(rdata_path) 
        data_dict = {}
        
        # Extract expression data 
        expr_key = [k for k in result.keys()  if 'exp_data' in k.lower()] 
        if expr_key:
            expr_df = result[expr_key[0]]
            expr_df.index  = expr_df.index.astype(str) 
            expr_df.columns  = expr_df.columns.astype(str) 
            data_dict['expression'] = expr_df 
        
        # Extract sample annotations 
        anno_key = [k for k in result.keys()  if 'sample_anno' in k.lower()] 
        if anno_key:
            sample_anno = result[anno_key[0]]
            if 'RNR' in sample_anno.columns: 
                sample_anno['label'] = sample_anno['RNR'].map({'NR': 0, 'R': 1})
                data_dict['sample_anno'] = sample_anno 
                
        return data_dict 
    except Exception as e:
        print(f"Error loading {os.path.basename(rdata_path)}:  {str(e)}")
        return {}
 
def load_all_data():
    """Load all RData files and ceRNA pairs"""
    ceRNA_df = pd.read_csv(CERNAS_PATH) 
    all_features = []
    biomarker_data = []
    
    for filename in os.listdir(RNR_DIR): 
        if not filename.endswith(".RData"): 
            continue 
            
        # Extract cancer type from filename 
        cancer_type = filename.split('_')[1] 
        rdata_path = os.path.join(RNR_DIR,  filename)
        rdata_dict = load_rdata_file(rdata_path)
        
        if not rdata_dict or 'expression' not in rdata_dict or 'sample_anno' not in rdata_dict:
            continue 
            
        expr_df = rdata_dict['expression']
        sample_anno = rdata_dict['sample_anno']
        
        # Find common samples 
        common_samples = list(set(expr_df.columns)  & set(sample_anno['sample_id']))
        if not common_samples:
            continue 
            
        # Process PD-L1 biomarker if available 
        if 'CD274' in expr_df.index: 
            biomarker_data.append(pd.DataFrame({ 
                'sample_id': common_samples,
                'cancer_type': cancer_type,
                'PD_L1': expr_df.loc['CD274',  common_samples].astype(float),
                'label': sample_anno.set_index('sample_id').loc[common_samples,  'label']
            }))
            
        # Process ceRNA pairs 
        cancer_ceRNAs = ceRNA_df[ceRNA_df['Cancer'] == cancer_type]
        for _, row in cancer_ceRNAs.iterrows(): 
            lnc, gene = row['lncRNA'], row['gene']
            if lnc not in expr_df.index  or gene not in expr_df.index: 
                continue 
                
            # Calculate features 
            lnc_expr = expr_df.loc[lnc,  common_samples].astype(float)
            gene_expr = expr_df.loc[gene,  common_samples].astype(float)
            
            features = pd.DataFrame({
                'sample_id': common_samples,
                'cancer_type': cancer_type,
                'ceRNA_pair': f"{lnc}_{gene}",
                'label': sample_anno.set_index('sample_id').loc[common_samples,  'label'],
                f"{lnc}_{gene}_lnc_expr": lnc_expr,
                f"{lnc}_{gene}_gene_expr": gene_expr,
                f"{lnc}_{gene}_comp_index": lnc_expr / (lnc_expr + gene_expr + 1e-5),
                f"{lnc}_{gene}_corr": np.corrcoef(lnc_expr,  gene_expr)[0, 1],
                f"{lnc}_{gene}_expr_diff": lnc_expr - gene_expr,
                f"{lnc}_{gene}_interaction": lnc_expr * gene_expr * abs(np.corrcoef(lnc_expr,  gene_expr)[0, 1])
            })
            all_features.append(features) 
    
    return pd.concat(all_features),  pd.concat(biomarker_data)  if biomarker_data else None 