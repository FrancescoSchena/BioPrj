import os
import pandas as pd
import statsmodels.api as sm
from scipy import stats
import seaborn as sn
import matplotlib.pyplot as plt
from sklearn.metrics import r2_score


# Aminoacids
path = "results/aminoacids"
if os.path.isdir(path):
	print(f"The directory {path} exists already")
else:
	try:
		os.mkdir(path)
	except OSError:
		print("Creation of the directory %s failed" % path)
		quit()
	else:
		print("Successfully created the directory %s" % path)
df = pd.read_csv("results/data.csv")C

print("Plotting Aminoacids (x23)")
aminoacids = [
			"PhePer",
			"LeuPer",
			"MetPer",
			"ValPer",
			"SerPer",
			"IsoPer",
			"ProPer",
			"ThrPer",
			"AlaPer",
			"TyrPer",
			"HisPer",
			"GlnPer",
			"AsnPer",
			"LysPer",
			"AspPer",
			"GluPer",
			"CysPer",
			"TrpPer",
			"ArgPer",
			"GlyPer"
			]
x = df["IntergenomicGC"]
for aminoacid in aminoacids:
	y = df[aminoacid]    
	slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
	R_square = r2_score(x, y)
	sn.regplot(x=x, y=y, line_kws={'label':"y = {0:.2f}+{1:.2f}x | R {2:.2f} | R squared {3:.2f}".format(intercept, slope, r_value,R_square)})
	plt.legend()
	plt.savefig(path + "/" + aminoacid + ".jpg")
	plt.close("all")