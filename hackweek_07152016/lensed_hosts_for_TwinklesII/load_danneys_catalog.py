import numpy as np
import pickle
df = pickle.load(open('./data/glsn_protodc2.pkl','rb'))

# print type(df)
# print df['zl'][7]
# print df['zs']
# print df['LSST_filters/magnitude:LSST_u:observed']
# print df['e']
# print df['imno']
print df['xhost']
print df['yhost']
