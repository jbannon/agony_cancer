import seaborn as sns
import matplotlib.pyplot as plt 
import pandas as pd 
import numpy as np 
from scipy.stats import wilcoxon

def plot_agony_dist(df):
	ax= sns.kdeplot(data=df, x="agony", hue="sample", multiple="stack")
	ax.set_title("Agony Distributions")
	plt.savefig("../figs/global_agony.png")
	plt.close()
	ax= sns.kdeplot(data=df, x="log_agony", hue="sample", multiple="stack")
	ax.set_title("Log Agony Distributions")
	plt.savefig("../figs/global_log_agony.png")


def reshape_dataset(df):
	tumor_df = pd.DataFrame({"agony":df["tumor_agonies"].tolist(),"log_agony":df["log_tumor_agonies"].tolist(), "sample":["tumor" for x in range(df.shape[0])] })
	normal_df = pd.DataFrame({"agony":df["normal_agony"].tolist(),"log_agony":df["log_normal_agonies"].tolist(), "sample":["normal" for x in range(df.shape[0])] })
	full_df = pd.concat((tumor_df,normal_df), ignore_index=True)
	return full_df, tumor_df, normal_df



df = pd.read_csv("../data/results_data/tcga/ALL/global_agonies.csv",index_col=0)
df,tdf, ndf = reshape_dataset(df)
print("two sided")
w, p = wilcoxon(tdf["agony"].tolist(), ndf["agony"].tolist())
print(w)
print(p)
print("less")
w, p = wilcoxon(tdf["agony"].tolist(), ndf["agony"].tolist(),alternative="less")
print("greater")
print(w)
print(p)
w, p = wilcoxon(tdf["agony"].tolist(), ndf["agony"].tolist(),alternative="greater")
print(w)
print(p)
plot_agony_dist(df)