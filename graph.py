import matplotlib.pyplot as py
import pandas as pd

data = pd.import_csv("output_g.csv")
analytic = data.iloc[0]

h1_f=abs(data.iloc[1] -analytic)
h1_c=abs(data.iloc[2]-analytic)
h01_f=abs(analytic-data.iloc[3])
h01_c=abs(data.iloc[4]- analytic)
h001_f=abs(data.iloc[5]-analytic)
h001_c=abs(data.iloc[6]-analytic)
h0001_f=abs(data.iloc[7]-analytic)
h0001_c=abs(data.iloc[8]-analytic)
analytic = abs(data.iloc[0]-analytic)

h = [ 0.1, 0.01, 0.001, 0.001]

