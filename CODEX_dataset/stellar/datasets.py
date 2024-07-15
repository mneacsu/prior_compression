import torch
import numpy as np
from builtins import range
from torch_geometric.data import InMemoryDataset, Data
from sklearn.metrics import pairwise_distances
import pandas as pd
import csv

def get_hubmap_edge_index(pos, regions, distance_thres):
    edge_list = []
    regions_unique = np.unique(regions)
    for reg in regions_unique:
        locs = np.where(regions == reg)[0]
        pos_region = pos[locs, :]
        dists = pairwise_distances(pos_region)
        dists_mask = dists < distance_thres
        np.fill_diagonal(dists_mask, 0)
        region_edge_list = np.transpose(np.nonzero(dists_mask)).tolist()
        for (i, j) in region_edge_list:
            edge_list.append([locs[i], locs[j]])
    return edge_list


def load_hubmap_data(labeled_file, unlabeled_file, distance_thres):
    train_df = pd.read_csv(labeled_file)
    test_df = pd.read_csv(unlabeled_file)
    train_X = train_df.iloc[:, :-4].values
    test_X = test_df.iloc[:, :-3].values
    labeled_pos = train_df[['x','y']].values
    unlabeled_pos = test_df[['x','y']].values
    labeled_regions = train_df['unique_region'].values
    unlabeled_regions = test_df['unique_region'].values
    train_y = train_df['cell_type']
    cell_types = np.sort(list(set(train_df['cell_type'].values))).tolist()
    cell_type_dict = {}
    for i, cell_type in enumerate(cell_types):
        cell_type_dict[cell_type] = i
    with open(labeled_file.split('_')[0] + '.csv', 'w') as csvfile:
        for key in cell_type_dict.keys():
            csvfile.write("%s,%s\n"%(key,cell_type_dict[key]))
    train_y = np.array([cell_type_dict[x] for x in train_y])
    labeled_edges = get_hubmap_edge_index(labeled_pos, labeled_regions, distance_thres)
    unlabeled_edges = get_hubmap_edge_index(unlabeled_pos, unlabeled_regions, distance_thres)
    return train_X, train_y, test_X, labeled_edges, unlabeled_edges


class CodexGraphDataset(InMemoryDataset):

    def __init__(self, labeled_X, labeled_y, unlabeled_X, labeled_edges, unlabeled_edges, transform=None,):
        self.root = '.'
        super(CodexGraphDataset, self).__init__(self.root, transform)
        self.labeled_data = Data(x=torch.FloatTensor(labeled_X), edge_index=torch.LongTensor(labeled_edges).T, y=torch.LongTensor(labeled_y))
        self.unlabeled_data = Data(x=torch.FloatTensor(unlabeled_X), edge_index=torch.LongTensor(unlabeled_edges).T)

    def __len__(self):
        return 2

    def __getitem__(self, idx):
        return self.labeled_data, self.unlabeled_data