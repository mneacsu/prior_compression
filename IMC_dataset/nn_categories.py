import pandas as pd
import numpy as np
import torch
import torch.nn as nn
from torch.utils.data import DataLoader, TensorDataset
from sklearn.preprocessing import OneHotEncoder
import random
from imblearn.over_sampling import RandomOverSampler
import sys

class NetCellCategory(nn.Module):
   
    def __init__(self, no_features):
        
        super().__init__()
        self.fc1 = nn.Linear(no_features, 4)
        self.fc2 = nn.Linear(4, 4)
        self.softmax = nn.Softmax(dim=1)
        
    def forward(self, x):
        x = self.fc1(x)
        x= torch.relu(x)
        x=self.fc2(x)
        x= self.softmax(x)
        
        return x


df = pd.read_csv('df_train.csv')
x= df.iloc[:,:-1].values
y= df.iloc[:,-1].values

torch.manual_seed(100)
random.seed(100)
np.random.seed(100)

os = RandomOverSampler()
x, y = os.fit_resample(x, y)

x= np.float32(x)
y = y.reshape(-1,1)

ohe = OneHotEncoder(sparse_output=False).fit(y)
y = ohe.transform(y)

model = NetCellCategory(int(sys.argv[1]))
model.train()

train_loader = DataLoader(dataset=TensorDataset(torch.tensor(x), torch.tensor(y)), batch_size=1000, shuffle=True)

n_epochs=5
loss_fn = nn.CrossEntropyLoss()
optimizer = torch.optim.Adam(model.parameters(), lr=0.1) 

for epoch in range(n_epochs):   
    for [X,Y] in train_loader:
        y_pred = model(X)
        optimizer.zero_grad()
        loss = loss_fn(y_pred, Y)
        loss.backward()
        optimizer.step()

df = pd.read_csv('df_all.csv')

x= torch.tensor(np.float32(df.values))

model.eval()

pd.DataFrame(model(x).detach().numpy(), columns=ohe.categories_).to_csv('probabilities.csv', index=False)
