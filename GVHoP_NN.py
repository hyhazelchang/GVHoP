#!/usr/bin/env python

# GVHoP_NN.py

import torch
from torch.utils.data import Dataset
import torch.nn as nn
from torch.utils.data import DataLoader
from sklearn.compose import ColumnTransformer
from sklearn.preprocessing import StandardScaler
from torch.utils.data import Sampler, RandomSampler, SequentialSampler
import joblib


class MetaDataset(Dataset):
    def __init__(self, X, y):
        self.X = torch.tensor(X, dtype=torch.float32)
        self.y = torch.tensor(y, dtype=torch.long)
    def __len__(self):
        return len(self.X)
    def __getitem__(self, idx):
        return self.X[idx], self.y[idx]

class MetaNN(nn.Module):
    def __init__(self, n_classes):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(965, 64),
            nn.BatchNorm1d(64),
            nn.ReLU(),
            nn.Dropout(p=0.3),
            nn.Linear(64, 32),
            nn.BatchNorm1d(32),
            nn.ReLU(),
            nn.Dropout(p=0.2),
            nn.Linear(32, 16),
            nn.ReLU(),
            nn.Linear(16, n_classes)
        )
    def forward(self, x):
        return self.net(x)

class SafeBatchSampler(Sampler):
    def __init__(self, dataset, batch_size, shuffle=True):
        self.batch_size = batch_size
        # Use RandomSampler for training, SequentialSampler for validation
        self.sampler = RandomSampler(dataset) if shuffle else SequentialSampler(dataset)
    def __iter__(self):
        batch = []
        prev_batch = None
        for idx in self.sampler:
            batch.append(idx)
            if len(batch) == self.batch_size:
                # If we already have a buffered batch, yield it now
                if prev_batch is not None:
                    yield prev_batch
                prev_batch = batch
                batch = []
        # Time to handle the leftover elements at the end of the epoch
        if batch:
            # CRITICAL CHECK: If the last batch size is exactly 1, 
            # merge it into the previous batch.
            if len(batch) == 1 and prev_batch is not None:
                prev_batch.extend(batch)
                yield prev_batch
            else:
                if prev_batch is not None:
                    yield prev_batch
                yield batch
        else:
            if prev_batch is not None:
                yield prev_batch
    def __len__(self):
        # Approximate length
        return (len(self.sampler) + self.batch_size - 1) // self.batch_size
