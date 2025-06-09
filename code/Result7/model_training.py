import xgboost as xgb 
import numpy as np 
from sklearn.model_selection  import StratifiedKFold 
from sklearn.metrics  import roc_auc_score, roc_curve 
from imblearn.pipeline  import Pipeline 
from imblearn.over_sampling  import SMOTE 
from scipy.stats  import ttest_rel 
from config import XGB_PARAMS, CV_FOLDS, OUTPUT_DIR 
import joblib 
import os 
 
def train_xgb_model(X_train, y_train):
    """Train XGBoost model with SMOTE"""
    pipeline = Pipeline([
        ('smote', SMOTE(random_state=42)),
        ('xgb', xgb.XGBClassifier(**XGB_PARAMS))
    ])
    pipeline.fit(X_train,  y_train)
    return pipeline 
 
def cross_validate(X, y, pdl1_data=None):
    """Perform stratified 5-fold cross-validation"""
    model_names = ['ceRNA_only', 'PD_L1_only', 'ceRNA_PDL1']
    results = {name: {'auc': [], 'fpr': [], 'tpr': []} for name in model_names}
    
    # Prepare PD-L1 data 
    X_pdl1 = pdl1_data.values  if pdl1_data is not None else np.random.normal(size=(len(y),  1))
    
    # Initialize KFold 
    skf = StratifiedKFold(n_splits=CV_FOLDS, shuffle=True, random_state=42)
    
    for fold, (train_idx, test_idx) in enumerate(skf.split(X,  y)):
        print(f"\nFold {fold + 1}/{CV_FOLDS}")
        
        # Split data 
        X_train, X_test = X.iloc[train_idx],  X.iloc[test_idx] 
        y_train, y_test = y.iloc[train_idx],  y.iloc[test_idx] 
        
        # 1. ceRNA-only model 
        model_ceRNA = train_xgb_model(X_train, y_train)
        proba_ceRNA = model_ceRNA.predict_proba(X_test)[:,  1]
        
        # 2. PD-L1-only model 
        model_pdl1 = train_xgb_model(X_pdl1[train_idx], y_train)
        proba_pdl1 = model_pdl1.predict_proba(X_pdl1[test_idx])[:,  1]
        
        # 3. Combined model 
        X_train_comb = np.hstack([X_train,  X_pdl1[train_idx]])
        X_test_comb = np.hstack([X_test,  X_pdl1[test_idx]])
        model_comb = train_xgb_model(X_train_comb, y_train)
        proba_comb = model_comb.predict_proba(X_test_comb)[:,  1]
        
        # Store results 
        for name, proba in zip(model_names, [proba_ceRNA, proba_pdl1, proba_comb]):
            auc = roc_auc_score(y_test, proba)
            fpr, tpr, _ = roc_curve(y_test, proba)
            results[name]['auc'].append(auc)
            results[name]['fpr'].append(fpr)
            results[name]['tpr'].append(tpr)
            print(f"{name}: AUC = {auc:.4f}")
    
    # Calculate statistics 
    metrics = {
        'ceRNA_mean_auc': np.mean(results['ceRNA_only']['auc']), 
        'ceRNA_std_auc': np.std(results['ceRNA_only']['auc']), 
        'PDL1_mean_auc': np.mean(results['PD_L1_only']['auc']), 
        'PDL1_std_auc': np.std(results['PD_L1_only']['auc']), 
        'Combined_mean_auc': np.mean(results['ceRNA_PDL1']['auc']), 
        'Combined_std_auc': np.std(results['ceRNA_PDL1']['auc']) 
    }
    
    # Statistical tests 
    pval_ceRNA = ttest_rel(results['ceRNA_PDL1']['auc'], results['ceRNA_only']['auc']).pvalue 
    pval_pdl1 = ttest_rel(results['ceRNA_PDL1']['auc'], results['PD_L1_only']['auc']).pvalue 
    
    # Save final models 
    final_models = {
        'ceRNA_model': train_xgb_model(X, y),
        'pdl1_model': train_xgb_model(X_pdl1, y) if pdl1_data is not None else None,
        'combined_model': train_xgb_model(np.hstack([X,  X_pdl1]), y)
    }
    joblib.dump(final_models,  os.path.join(OUTPUT_DIR,  'final_models.pkl')) 
    
    return results, metrics, pval_ceRNA, pval_pdl1 