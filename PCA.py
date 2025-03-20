import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

import Aim1
import DifferenceVector as dv

cosmic_percent_rank = Aim1.cosmic_percent_rank


def plot_pca(provided_df, x_axis='PC1', y_axis='PC2', use_kde=1, dot_alpha=0.3, dot_size=10, use_name_to_save='plot',
             draw_h_line=0, draw_v_line=0, plot_points=None, use_marginals=0,
             remove_legend=0, bw_adjust=1):

    plt.figure(figsize=(20, 20))

    g = sns.JointGrid(data=provided_df, x=x_axis, y=y_axis, height=8, ratio=5)

    if use_kde:
        # plot kde
        g = g.plot_joint(sns.kdeplot, alpha=0.9, bw_adjust=bw_adjust)

    # plot scatterplot on top
    g.plot_joint(sns.scatterplot, s=dot_size, edgecolor='black', alpha=dot_alpha)

    if use_marginals:
        g.plot_marginals(sns.kdeplot)

    if remove_legend:
        # remove the lefend from the joint plot
        g.ax_joint.legend_.remove()

    g.ax_joint.set_xlim(-40, 40)
    g.ax_joint.set_ylim(-40, 40)

    # adjust font size for ax_joint
    g.ax_joint.set_xlabel(x_axis, fontsize=18)
    g.ax_joint.set_ylabel(y_axis, fontsize=18)
    g.ax_joint.tick_params(axis='both', which='major', labelsize=14)

    if draw_h_line:
        # draw a horizontal line
        g.ax_joint.axhline(y=draw_h_line, color='darkgrey', linestyle='--')

    if draw_v_line:
        # draw a vertical line
        g.ax_joint.axhline(x=draw_v_line, color='darkgrey', linestyle='--')

    if plot_points is not None:
        additional_x = plot_points['PC-1'].values
        additional_y = plot_points['PC-2'].values
        # plot the additional points on top of the existing JointGrid plot
        g.ax_joint.plot(additional_x, additional_y, 'ko', markersize=3)
        #'ko' means black color and circle markers

    plt.tight_layout()
    plt.savefig(use_name_to_save, transparent=True, dpi=300)


def get_pca(percent_rank_df):

    #threshold = input('Set threshold for difference vector: ')
    #print(f'Threshold set to {threshold}')

    #dv_dict = dv.percent_change(percent_rank_df, threshold)
    #dv_df = dv.change_df(dv_dict)
    dv_df = Aim1.readCSV('dv_cosmic_flank_st_threshold_90.csv')
    #dv_df = dv_df_1.reset_index(drop=True)

    # Example df
    columns = dv_df.columns.values.tolist()
    index = dv_df.index.values.tolist()
    num_samples = len(index)
    X = pd.DataFrame(dv_df, columns=columns, index=index)
    # Remove columns with any missing values
    # X_cleaned = X.dropna(axis=1)

    # Remove any columns where every column is equal to 0
    X_cleaner = X.drop(columns=[col for col in X if (X[col] == 0).all()])
    X_cleaned = X_cleaner.set_index('Unnamed: 0')
    print(X_cleaned)

    clean_col = X_cleaned.columns.values.tolist()
    num_cols = len(clean_col)

    # Feature Scaling
    scaler = StandardScaler()
    df_scaled = scaler.fit_transform(X_cleaned)

    # PCA Execution
    n_components = min(num_samples, num_cols)
    #n_components = 10
    pca = PCA(n_components=n_components, svd_solver='full')

    principalComponents = pca.fit_transform(df_scaled)
    explained_variance_ratio = pca.explained_variance_ratio_
    cumulative_explained_variance = np.cumsum(explained_variance_ratio)

    # Cumulative Explained Variance to decide on number of components
    pc_index = 1
    for i in cumulative_explained_variance:
        print(f'PC{pc_index} = {i}')
        pc_index += 1

    # Plotting Cumulative Explained Variance to decide on number of components
    plt.figure(figsize=(8, 5))
    plt.plot(range(1, len(cumulative_explained_variance) + 1), cumulative_explained_variance, marker='o',
             linestyle='--')
    plt.title('Cumulative Explained Variance by PCA Components', fontsize=16)
    plt.xlabel('Number of Components', fontsize=14)
    plt.ylabel('Cumulative Explained Variance', fontsize=14)
    plt.grid(True)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.savefig("pca_cumvar_st_flank_mut_negative.png")

    # Extracting the first 2 principal components
    pc1 = principalComponents[:, 0]
    pc2 = principalComponents[:, 1]

    M = pd.DataFrame({'PC1': pc1, 'PC2': pc2})
    pca_data = pd.DataFrame(pca.fit_transform(X),columns=['PC1','PC2']) 
    Aim1.to_csv(M, "pca_data_st_flank")
    plot_pca(M, use_marginals=1, dot_size=30, dot_alpha=1, use_name_to_save="pca_flank_st_mut_negative.png")


get_pca()

