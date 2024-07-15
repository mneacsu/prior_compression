import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.neighbors import KNeighborsClassifier
from sklearn.ensemble import RandomForestClassifier
from imblearn.over_sampling import SMOTE, RandomOverSampler

df = pd.read_csv('types_training.csv')

df2 = pd.read_csv('types_all.csv')
df2['cell_type']=""

for t in ['Immune', 'Epithelial', 'Stromal']:

    match t:
        case "Immune":
            markers=["CD8",   "CD3",   "CD4",   "CD19",  "CD21" , "CD38", "CD68",  "HLADR", "CD163", "CD206", "CD15",  "CD11c" ,"NKG2D", "CD161",  "CD66" ]
        case "Epithelial":
            markers=["CHGA" , "MUC2" , "aDef5", "SOX9" ]
        case "Stromal":
            markers=["aSMA" ,      "CD31"  ,     "Synapto"   , "Podoplanin" ,"CD117" ]

    x= df.loc[df['cell_category']==t,markers].values
    y= df.loc[df['cell_category']==t,'cell_type'].values

    classifier = KNeighborsClassifier()
    classifier.fit(x, y)

    df2.loc[df2['cell_category']==t,'cell_type']=classifier.predict(df2.loc[df2['cell_category']==t,markers].values)

df2.to_csv('cell_types.csv')


