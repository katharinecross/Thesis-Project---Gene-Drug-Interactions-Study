# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import pickle
import matplotlib.pyplot as plt
import numpy as np
with open("/Users/tingtingzhao/Documents/Research/data/rubicinAndOther/readcombinedDatavorinostat.pkl", 'rb') as f:
    GeneMtx1, combined_df1 = pickle.load(f)
#with open(r'C:\Users\student\OneDrive - Bryant University\Desktop\Honors Thesis\readcombinedDatavorinostat.pkl', 'rb') as f:
#    GeneMtx1, combined_df1 = pickle.load(f)
GeneMtx1.shape
GeneMtx1Sub = GeneMtx1[0:3107,:]

GeneMtx1Sub.shape

from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score

# Kmeans clustering
kmeans_per_k = [KMeans(n_clusters=k, random_state=42).fit(GeneMtx1Sub)
                for k in range(1, 11)]
inertias = [model.inertia_ for model in kmeans_per_k]
silhouette_scores = [silhouette_score(GeneMtx1Sub, model.labels_)
                     for model in kmeans_per_k[1:]]      

plt.figure(figsize=(8, 4))
plt.plot(range(2, 11), silhouette_scores, "bo-")
plt.xlabel("$k$", fontsize=14)
plt.ylabel("Silhouette score", fontsize=14)
plt.show()

#The elbow curve

wcss=[]

from sklearn.metrics import silhouette_score
silhouette_scores = []
for i in range(1,11):
    km=KMeans(n_clusters=i)
    km.fit(GeneMtx1Sub)
    wcss.append(km.inertia_)
    

plt.figure(figsize=(12,6))

plt.plot(range(1,11),wcss)

plt.plot(range(1,11),wcss, linewidth=2, color="red", marker ="8")

plt.xlabel("K Value")
plt.xticks(np.arange(1,11,1))
plt.ylabel("WCSS")

plt.show()


from sklearn.metrics import silhouette_samples
from matplotlib.ticker import FixedLocator, FixedFormatter
import matplotlib as mpl



km1=KMeans(n_clusters=2)

#Fitting the input data

km1.fit(GeneMtx1Sub)

#predicting the labels of the input data

y=km1.predict(GeneMtx1Sub)

#adding the labels to a column named label
import pandas as pd

df = pd.DataFrame(data=GeneMtx1Sub)
df['y'] = y

## Exploratory Data Analysis

## explore the difference between the two groups using y labels
# Group data by label and calculate mean for each class
mean_by_class = df.groupby(y).mean()
mean_by_class
# Compare mean values between the two classes
print(mean_by_class)

import matplotlib.pyplot as plt
import seaborn as sns

# Separate data for each class
class_0_data = mean_by_class.iloc[0]
class_1_data = mean_by_class.iloc[1]

# Number of bins for the histogram
num_bins = 20

# Create subplots with one row and two columns
fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(12, 6))

# Create a histogram for class 0 in the first subplot
axes[0].hist(class_0_data.values.flatten(), bins=num_bins, color='blue')
axes[0].set_title('Custer 0')
axes[0].set_xlabel('Feature Value')
axes[0].set_ylabel('Frequency')
axes[0].set_xlim(-20, 20)  # Set x-axis range

# Create a histogram for class 1 in the second subplot
axes[1].hist(class_1_data.values.flatten(), bins=num_bins, color='orange')
axes[1].set_title('Cluster 1')
axes[1].set_xlabel('Feature Value')
axes[1].set_ylabel('Frequency')
axes[1].set_xlim(-20, 20)  # Set x-axis range

# Adjust layout for better spacing between subplots
plt.tight_layout()

# Display the plot
plt.show()


import pandas as pd

# Combine class_0_data and class_1_data vertically
combined_data = pd.concat([class_0_data, class_1_data], axis=0)

# Create cluster labels (0 for class 0 and 1 for class 1)
cluster_labels = [0] * len(class_0_data) + [1] * len(class_1_data)

# Create a DataFrame
combined_df = pd.DataFrame({'Feature': combined_data.values.flatten(), 'Cluster': cluster_labels})

print(combined_df.head())

# Kernel Density Plot by Seaborn
#  https://seaborn.pydata.org/generated/seaborn.kdeplot.html

plt.figure(figsize=(10, 8))

sns.kdeplot(
   data=combined_df, x="Feature", hue="Cluster",
   fill=True, common_norm=False,
   alpha=.5, linewidth=0,
)

# Add labels and title
plt.xlabel('Cluster', fontsize=16, fontweight='bold')
plt.ylabel('', fontsize=16, fontweight='bold')
plt.title('Density Comparison of Landmark Genes Mean Z-score Values by Cluster', fontsize=16, fontweight='bold')

# Customize tick labels
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)

# Save the plot as a PDF file
plt.savefig('/Users/tingtingzhao/Documents/Teaching/HornorsThesis/HonorThesisKatie/density_comparison.pdf', format='pdf')

# Display the plot
plt.show()

# Create a boxplot using Seaborn
plt.figure(figsize=(10, 8))
sns.boxplot(x='Cluster', y='Feature', data=combined_df)

# Add labels and title
plt.xlabel('Cluster', fontsize=16, fontweight='bold')
plt.ylabel('Feature Means', fontsize=16, fontweight='bold')
plt.title('Boxplot Comparison of Landmark Genes Mean Z-score Values by Cluster', fontsize=16, fontweight='bold')

# Customize tick labels
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)

# Save the plot as a PDF file
plt.savefig('/Users/tingtingzhao/Documents/Teaching/HornorsThesis/HonorThesisKatie/boxplot_comparison.pdf', format='pdf')
# Display the plot
plt.show()

## Randomly select 40 genes from class 0 and 
# randomly select 40 genes from class 1 and do a heatmap
import numpy as np
import pandas as pd

# Separate data for class 0 and class 1
class_0_data = df[df['y'] == 0].sample(n=20, random_state=0)
class_1_data = df[df['y'] == 1].sample(n=20, random_state=0)

combined_data = pd.concat([class_0_data, class_1_data], axis=0)

sns.heatmap(class_0_data, cmap='coolwarm')

# Extract only gene columns (first 978 columns)
gene_data = combined_data.iloc[:, :978]

# Create a heatmap
plt.figure(figsize=(10, 8))
sns.heatmap(gene_data, cmap='coolwarm')

# Customize axis labels and title
plt.xlabel('Genes', fontsize=14)
plt.ylabel('Observations', fontsize=14)
plt.title('Heatmap of Gene Expression by Observations', fontsize=16)
# Display the plot
plt.show()

######################################################
## write df to a csv file and save it for later 
## so that we can generate the Gaussian mixture knockoffs
## get the 978 gene names and put it as column names
## Find the gene names after feature selection
gene_info_path="/Users/tingtingzhao/Documents/Research/data/GSE92742_Broad_LINCS_gene_info.txt"
gene_info = pd.read_csv(gene_info_path, sep="\t")
numOfGenes = GeneMtx1Sub.shape[1]
geneNames = gene_info['pr_gene_symbol'][0:numOfGenes]
columnNames = list(geneNames)
columnNames.append("cluster")
df.columns = columnNames
df = df.rename(columns=columnNames)
df.to_csv("/Users/tingtingzhao/Documents/Teaching/HornorsThesis/HonorThesisKatie/TwoClusterGenes.csv", index=False)








# Statistical Test to Check Significant Difference
#from scipy.stats import ttest_ind

#group_0 = combined_df[combined_df['Cluster'] == 0]['Feature']
#group_1 = combined_df[combined_df['Cluster'] == 1]['Feature']

#t_statistic, p_value = ttest_ind(group_0, group_1)
#print("t-statistic:", t_statistic)
#print("p-value:", p_value)

# Pairplot to visualize relationships between features by class
#sns.pairplot(data=data, hue='y', diag_kind='kde')
#plt.title('Pairplot of Features by Class')
#plt.show()
