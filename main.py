"""
@author: yohanna
@email: yohanna.wang0924@gmail.com
"""

from utils import args
from data import load_data
from evaluate import count_accuracy
from dcov.dcov_known_cc import DCOV_known
from dcov.dcov_unknown_cc import DCOV_unknown


if __name__ == '__main__':
    
    simu, G_true, _, cc, _ = load_data()

    if args.algorithm == "known_dcov":
        G_est = DCOV_known(simu, G_true, cc)    
        
    if args.algorithm == "unknown_dcov":
        G_est = DCOV_unknown(simu, G_true, cc)
        

    fdr_npcov, tpr_npcov, fpr_npcov, shd_npcov, pred_size_npcov = count_accuracy(G_true, G_est)
    print('DCOV Accuracy: fdr', fdr_npcov, ' tpr ', tpr_npcov, ' fpr ', fpr_npcov, 'shd', shd_npcov, 'pred_size', pred_size_npcov)


    