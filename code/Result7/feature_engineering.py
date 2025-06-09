import pandas as pd 
import numpy as np 
from sklearn.preprocessing  import StandardScaler 
from sklearn.impute  import SimpleImputer 
from sklearn.feature_selection  import VarianceThreshold, SelectKBest, f_classif 
from config import VARIANCE_THRESHOLD, SELECT_K_FEATURES 
 
def preprocess_features(full_data):
    """Clean and preprocess features"""
    # Select feature columns 
    feature_cols = [col for col in full_data.columns  if col not in 
                   ['sample_id', 'cancer_type', 'label', 'ceRNA_pair']]
    X = full_data[feature_cols]
    y = full_data['label']
    
    # 1. Handle missing values 
    X.replace([np.inf,  -np.inf],  np.nan,  inplace=True)
    imputer = SimpleImputer(strategy='median')
    X_imputed = pd.DataFrame(imputer.fit_transform(X),  columns=X.columns,  index=X.index) 
    
    # 2. Variance threshold 
    selector = VarianceThreshold(threshold=VARIANCE_THRESHOLD)
    X_var = selector.fit_transform(X_imputed) 
    retained_cols = X_imputed.columns[selector.get_support()] 
    X = pd.DataFrame(X_var, columns=retained_cols, index=X_imputed.index) 
    
    # 3. Feature scaling 
    scaler = StandardScaler()
    X_scaled = pd.DataFrame(scaler.fit_transform(X),  columns=X.columns,  index=X.index) 
    
    # 4. Feature selection 
    selector = SelectKBest(f_classif, k=min(SELECT_K_FEATURES, X_scaled.shape[1])) 
    X_selected = selector.fit_transform(X_scaled,  y)
    selected_cols = X_scaled.columns[selector.get_support()] 
    
    return pd.DataFrame(X_selected, columns=selected_cols, index=X_scaled.index),  y 
 
def prepare_pdl1_data(pdl1_data):
    """Prepare PD-L1 biomarker data"""
    if pdl1_data is None:
        return None 
    
    # Create wide format for PD-L1 
    pdl1_wide = pdl1_data.pivot(index='sample_id',  columns='cancer_type', values='PD_L1')
    pdl1_wide.columns  = [f"PD_L1_{col}" for col in pdl1_wide.columns] 
    
    # Fill missing values with median 
    pdl1_wide = pdl1_wide.fillna(pdl1_wide.median()) 
    
    return pdl1_wide 