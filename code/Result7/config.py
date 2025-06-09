import os 
 
# Path configurations 
RNR_DIR = "data/RNR_DIR/bulk/RNR"
CERNAS_PATH = "data/all-ceRNA list.csv" 
OUTPUT_DIR = "./OUTPUT_DIR"
 
# XGBoost parameters 
XGB_PARAMS = {
    'objective': 'binary:logistic',
    'eval_metric': 'auc',
    'random_state': 42,
    'learning_rate': 0.1,
    'max_depth': 7,
    'subsample': 0.8,
    'colsample_bytree': 0.8,
    'min_child_weight': 1,
    'gamma': 0,
    'n_estimators': 100 
}
 
# Cross-validation settings 
CV_FOLDS = 5 
TEST_SIZE = 0.2 
STRATIFIED = True 
 
# Feature processing 
FEATURE_TYPES = ['lnc_expr', 'gene_expr', 'comp_index', 'corr', 'expr_diff', 'interaction']
VARIANCE_THRESHOLD = 0.01 
SELECT_K_FEATURES = 150 
 
# Visualization settings 
PLOT_PARAMS = {
    'font.family':  'sans-serif',
    'font.sans-serif':  ['Arial', 'Helvetica'],
    'font.size':  10,
    'axes.labelsize':  11,
    'axes.titlesize':  12,
    'figure.dpi':  300,
    'savefig.dpi':  300,
    'savefig.format':  'pdf'
}
 
# Create output directory 
os.makedirs(OUTPUT_DIR,  exist_ok=True)