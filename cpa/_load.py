"""Functionalities to load data from splits to CPA AnnData with the correct fields."""

import scanpy as sc
from pathlib import Path
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def _validate_cols(adata: sc.AnnData) -> bool:
    """Validate that CPA-required columns are present."""
    required_cols = ['condition', 'dose_val']
    for col in required_cols:
        if col not in adata.obs.columns:
            logging.error(f'Column {col} not found in adata.obs')
            return False
    return True

def _add_condition_dosage(adata: sc.AnnData,
                        condition_col: str = None,
                        dosage_col: str = None) -> sc.AnnData:
    """Add condition and dosage columns to adata.obs."""
    if condition_col is not None:
        adata.obs['condition'] = adata.obs[condition_col]
    else:
        adata.obs['condition'] = 'condition'
    if dosage_col is not None:
        adata.obs['dose_val'] = adata.obs[dosage_col]
    else:
        adata.obs['dose_val'] = 1.
    return adata

def load_split(train_path: str,
               valid_path: str,
               test_path: str,
               control_condition: str = 'GFP',
               condition_col: str = None,
               dosage_col: str = None) -> sc.AnnData:
    """Load a train/valid/test split and format for CPA."""
    adata_train = sc.read_h5ad(train_path)
    adata_train = _add_condition_dosage(adata_train, condition_col, dosage_col)
    assert _validate_cols(adata_train), "Required columns not found in adata.obs"
    adata_valid = sc.read_h5ad(valid_path)
    adata_valid = _add_condition_dosage(adata_valid, condition_col, dosage_col)
    assert _validate_cols(adata_valid), "Required columns not found in adata.obs"
    adata_test = sc.read_h5ad(test_path)
    adata_test = _add_condition_dosage(adata_test, condition_col, dosage_col)
    assert _validate_cols(adata_test), "Required columns not found in adata.obs"
    adata = adata_train.concatenate([adata_valid, adata_test],
                                    batch_key='split',
                                    batch_categories=['train', 'test', 'ood'])
    if 'control' not in adata.obs.columns:
        adata.obs['control'] = (adata.obs['condition'].values == control_condition).astype(int)
    else:
        logger.info('control column already present in adata.obs')
    return adata
    