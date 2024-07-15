import argparse
from utils import prepare_save_dir
from STELLAR import STELLAR
import numpy as np
import os
import torch
from datasets import CodexGraphDataset, load_hubmap_data

def main():
    parser = argparse.ArgumentParser(description='STELLAR')
    parser.add_argument('--dataset', default='Hubmap', help='dataset setting')
    parser.add_argument('--seed', type=int, default=1, metavar='S', help='random seed (default: 1)')
    parser.add_argument('--name', type=str, default='STELLAR')
    parser.add_argument('--epochs', type=int, default=20)
    parser.add_argument('--lr', type=float, default=1e-3)
    parser.add_argument('--wd', type=float, default=5e-2)
    parser.add_argument('--input-dim', type=int, default=48)
    parser.add_argument('--num-heads', type=int, default=22)
    parser.add_argument('--num-seed-class', type=int, default=0)
    #parser.add_argument('--sample-rate', type=float, default=0.5)
    parser.add_argument('-b', '--batch-size', default=1, type=int,
                    metavar='N',
                    help='mini-batch size')
    parser.add_argument('--distance_thres', default=50, type=int)
    parser.add_argument('--savedir', type=str, default='./')
    args = parser.parse_args()
    args.cuda = torch.cuda.is_available()
    args.device = torch.device("cuda" if args.cuda else "cpu")

    # Seed the run and create saving directory
    args.name = '_'.join([args.dataset, args.name])
    args = prepare_save_dir(args, __file__)
    
    labeled_X, labeled_y, unlabeled_X, labeled_edges, unlabeled_edges = load_hubmap_data('./data/'+args.dataset+'_training.csv', './data/'+args.dataset+'_all.csv', args.distance_thres)
    dataset = CodexGraphDataset(labeled_X, labeled_y, unlabeled_X, labeled_edges, unlabeled_edges)
    
    stellar = STELLAR(args, dataset)
    stellar.train()
    _, results = stellar.pred()
    np.savetxt(os.path.join(args.savedir, args.dataset + '_results.csv'), results, delimiter=",")

if __name__ == '__main__':
    main()
