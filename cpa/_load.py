"""Functionalities to load data from splits to CPA AnnData with the correct fields."""

import scanpy as sc
from pathlib import Path
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def _validate_cols(adata: sc.AnnData) -> bool:
    """Validate that CPA-required columns are present."""
    required_cols = ['split', 'condition', 'dose_val']
    for col in required_cols:
        if col not in adata.obs.columns:
            logging.error(f'Column {col} not found in adata.obs')
            return False
    return True

def load_split(train_path: str,
               valid_path: str,
               test_path: str,
               control_condition: str = 'GFP') -> sc.AnnData:
    """Load a train/valid/test split and format for CPA."""
    adata_train = sc.read_h5ad(train_path)
    adata_valid = sc.read_h5ad(valid_path)
    adata_test = sc.read_h5ad(test_path)
    adata = adata_train.concatenate([adata_valid, adata_test],
                                      batch_key='split',
                                        batch_categories=['train', 'test', 'ood'])
    assert _validate_cols(adata), "Required columns not found in adata.obs"
    if 'control' not in adata.obs.columns:
        adata.obs['control'] = (adata.obs['condition'].values == control_condition).astype(int)
    else:
        logger.info('control column already present in adata.obs')
    return adata
    