import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.neighbors import KNeighborsClassifier
from imblearn.over_sampling import RandomOverSampler

df = pd.read_csv('types_train.csv')

df2 = pd.read_csv('types_all.csv')
df2['cell_type']=""

for t in ['immune', 'exocrine', 'islet', 'stromal']:

    match t:
        case "immune":
            markers=["CD8a", "CD20",  "CD68",  "MPO", "CD4", "CD3e"]
        case "exocrine":
            markers=["PDX1", "AMY2A","KRT19"]
        case "islet":
            markers=["GCG", "PIN", "PPY", "SST", "PCSK2", "NKX6-1", "IAPP", "PDX1", "INS"]
        case "stromal":
            markers=["CD31","SMA"]

    x= df.loc[df['cell_category']==t,markers].values
    y= df.loc[df['cell_category']==t,'cell_type'].values

    os = RandomOverSampler(random_state=0)
    x, y = os.fit_resample(x, y)

    classifier = KNeighborsClassifier()
    classifier.fit(x, y)


    df2.loc[df2['cell_category']==t,'cell_type']=classifier.predict(df2.loc[df2['cell_category']==t,markers].values)

df2.to_csv('cell_types.csv')


