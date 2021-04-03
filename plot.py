import os
import pandas as pd
import statsmodels.api as sm
from scipy import stats
import seaborn as sn
import matplotlib.pyplot as plt
from sklearn.metrics import r2_score

df = pd.read_csv("results/data.csv")


# Aminoacids
path = "results/aminoacids"
# 	Creating the target directory
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
print("Plotting Aminoacids (x23)")
# 	List of aminoacids (columns of the df)
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
# 	Setting the common x
x = df["IntergenomicGC"]
for aminoacid in aminoacids:
	y = df[aminoacid]
	# Calculating the values to show in the label
	slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
	R_square = r2_score(x, y)
	# Plotting the graph and the label
	sn.regplot(x=x, y=y, line_kws={'label':"y = {0:.2f}+{1:.2f}x | R {2:.2f} | R squared {3:.2f}".format(intercept, slope, r_value, po
	R_square)})
	plt.legend()
	# Saving in the predetermined path
	plt.savefig(path + "/" + aminoacid + ".jpg")
	# Getting ready for the next plot
	plt.close("all")