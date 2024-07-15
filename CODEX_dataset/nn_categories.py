import pandas as pd
import numpy as np
import torch
import torch.nn as nn
from torch.utils.data import DataLoader, TensorDataset
from sklearn.preprocessing import OneHotEncoder
import random
from imblearn.over_sampling import RandomOverSampler, SMOTE
import sys
from sklearn.utils.class_weight import compute_class_weight

class NetCellCategory(nn.Module):
   
    def __init__(self, no_features):
        
        super().__init__()
        self.fc1 = nn.Linear(no_features, 5)
        self.fc2 = nn.Linear(5, 3)
        self.softmax = nn.Softmax(dim=1)
        
    def forward(self, x):
        x = self.fc1(x)
        x= torch.asinh(x)
        x=self.fc2(x)
        x= self.softmax(x)
        
        return x


df = pd.read_csv('categories_training.csv')
x= df.iloc[:,:-1].values
y= df.iloc[:,-1].values

os = RandomOverSampler(random_state=0, shrinkage=1)
x, y = os.fit_resample(x, y)


torch.manual_seed(100)
random.seed(100)
np.random.seed(100)

x= np.float32(x)
y = y.reshape(-1,1)

ohe = OneHotEncoder(sparse_output=False).fit(y)
y = ohe.transform(y)

model = NetCellCategory(int(sys.argv[1]))
model.train()

train_loader = DataLoader(dataset=TensorDataset(torch.tensor(x), torch.tensor(y)), batch_size=1000, shuffle=True)

n_epochs=5
loss_fn = nn.CrossEntropyLoss(weight = torch.tensor([25,1,1]))
optimizer = torch.optim.Adam(model.parameters(), lr=0.1) 

for epoch in range(n_epochs):   
    for [X,Y] in train_loader:
        y_pred = model(X)
        optimizer.zero_grad()
        loss = loss_fn(y_pred, Y)
        loss.backward()
        optimizer.step()

df = pd.read_csv('categories_all.csv')

x= torch.tensor(np.float32(df.values))

model.eval()

pd.DataFrame(model(x).detach().numpy(), columns=ohe.categories_).to_csv('probabilities.csv', index=False)
