import time,re,random
import pandas as pd
import numpy as np
import scipy.stats as stats
from scipy.stats import spearmanr
from scipy.cluster import hierarchy
from scipy.spatial.distance import squareform
from scipy.interpolate import splev, splrep 
import matplotlib.pyplot as plt
from matplotlib import rcParams
import seaborn as sns

import sklearn
from sklearn.ensemble import RandomForestRegressor
from sklearn.utils import shuffle
from sklearn.datasets import fetch_openml
from sklearn.impute import SimpleImputer
from sklearn.compose import ColumnTransformer
from sklearn.model_selection import train_test_split
from sklearn.model_selection import cross_val_score
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder
from sklearn.preprocessing import OrdinalEncoder
from sklearn import metrics
from sklearn.metrics import mean_squared_error
from sklearn.metrics import r2_score
from sklearn.metrics import accuracy_score
from sklearn.metrics import mean_absolute_error, mean_squared_error
from sklearn.inspection import permutation_importance
from sklearn.inspection import PartialDependenceDisplay
from sklearn.inspection import partial_dependence 
from sklearn.compose import ColumnTransformer
from collections import defaultdict

config = {
    "font.family":'serif',  
    "font.size": 14,
    'font.weight' : 'bold',
    "mathtext.fontset":'stix',
    "font.serif": ['SimSun'],
    "axes.unicode_minus": False
}
rcParams.update(config)
###########################
#['scale_slope','scale_slope_deficit','scale_slope_surplus','Sensitivity_slope','Sensitivity_slope_deficit','Sensitivity_slope_surplus']
th=0
writer = pd.ExcelWriter(r'D:\Droughts_scales\Output_data\Regressor_30drivers_sig_permute_decolinear_regions.xlsx',engine='openpyxl')
for name_sheet in ['scale_slope_deficit','Sensitivity_slope_deficit','scale_slope_surplus','Sensitivity_slope_surplus'] :  
    #name_sheet='scale_slope'
    print(name_sheet)
    th=th+1;
    if name_sheet.find('scale')>=0:
      df = pd.read_excel(r'D:\Droughts_scales\Output_data\30drivers_noclass_sig.xls',sheet_name='scale_slope')
    else: 
      df = pd.read_excel(r'D:\Droughts_scales\Output_data\30drivers_noclass_sig.xls',sheet_name='Sensitivity_slope')
    columns=df.columns.values
    cls=pd.Index([1,1,1,1,1,2,2,3,2,3,3,3,3,3,3,2,3,3,3,3,2,2,1,1,1,1,1,3,2,1]) #.astype('object').T
    x_columns = df.columns.tolist()
    x_columns.remove(columns[0])
    X = df[x_columns]
    y = df[columns[0]]
    X.columns=x_columns
    x_columns=[name.replace('_','') for name in x_columns]
    x_columns=pd.Index(x_columns)
    start =time.time()
    X_train,X_test,y_train,y_test = train_test_split(X,y,test_size=0.1,random_state=0)
    corr = spearmanr(X).correlation
    # Ensure the correlation matrix is symmetric
    corr = (corr + corr.T) / 2
    np.fill_diagonal(corr, 1)
    # We convert the correlation matrix to a distance matrix before performing
    # hierarchical clustering using Ward's linkage.
    distance_matrix = 1 - np.abs(corr)
    dist_linkage = hierarchy.ward(squareform(distance_matrix))

    cluster_ids = hierarchy.fcluster(dist_linkage, 0.75, criterion="distance") #chnage the thresholds 0.5 075 1
    cluster_id_to_feature_ids = defaultdict(list)
    for idx, cluster_id in enumerate(cluster_ids):
        cluster_id_to_feature_ids[cluster_id].append(idx)
    selected_features = [v[0] for v in cluster_id_to_feature_ids.values()]
    typ=cls[selected_features]
    selected_features_names = X.columns[selected_features]
    X_columns_sel=selected_features_names
    print(X_columns_sel)
    print(X_columns_sel.shape[0])
    ########################
    ########################
    df = pd.read_excel(r'D:\Droughts_scales\Output_data\30drivers_noclass_sig.xls',sheet_name=name_sheet)
    columns=df.columns.values
    cls=pd.Index([1,1,1,1,1,2,2,3,2,3,3,3,3,3,3,2,3,3,2,3,2,2,1,1,1,1,1,2,2,1]) #.astype('object').T
    x_columns = df.columns.tolist()
    x_columns.remove(columns[0])
    X = df[x_columns]
    y = df[columns[0]]
    X.columns=x_columns
    x_columns=[name.replace('_','') for name in x_columns]
    x_columns=pd.Index(x_columns)
    ########################
    X_train,X_test,y_train,y_test = train_test_split(X,y,test_size=0.1, random_state=0)
    X_train_sel = X_train[selected_features_names]
    X_test_sel = X_test[selected_features_names]
    ########################
    rf_sel = RandomForestRegressor(bootstrap=True,n_estimators=200,n_jobs=-1,oob_score=True,random_state=0,min_samples_leaf=1)
    rf_sel.fit(X_train_sel, y_train)
    print('Running time: %s Seconds'%( time.time()-start))

    # Error printing
    print("Baseline accuracy on train data with features removed:"f" {rf_sel.score(X_train_sel, y_train):.4}") # calculate the Gini scores of train data 
    print("Baseline accuracy on test data with features removed:"f" {rf_sel.score(X_test_sel, y_test):.4}")  # calculate the gini scores of test data
    ############
    ########################Importances  plot
    #####Permutation Importances   
    print('Permutation Importances (train set)...')  
    start =time.time()
    result =permutation_importance(rf_sel, X_train_sel, y_train, n_repeats=20, random_state=0, n_jobs=-1)
    print('Running time: %s Seconds'%( time.time()-start))
    importances=result.importances_mean 
    sorted_idx = importances.argsort() [::-1]
    names=['$%s$'%re.sub(r'(\s+)', r'\0$ $', i) for i in X_columns_sel[sorted_idx]]
    permute_importances =pd.Series(importances[sorted_idx], index=names)

    #text (chr(97+i))
    if name_sheet.find('scale')>=0:
      titl='$TS_d$ $Slope$'
    else:
      titl='$R_d$ $Slope$'
    ########################Importances  plot    
    ########################    
    # Save key parameters
    df_imp = pd.DataFrame()
    df_imp['RI_'+ name_sheet] = list(result.importances_mean)
    df_imp.index =X_columns_sel
    df_imp.to_excel(excel_writer=writer, sheet_name=name_sheet)
    ########################
    print("Computing partial dependence plots...")
    start=time.time()
    sorted_idx = result.importances_mean.argsort() [::-1]
    lb=X_columns_sel[sorted_idx]
    im=result.importances_mean[sorted_idx]
    typc=typ[sorted_idx]
    features = lb#list([lb[0],lb[1],lb[2],lb[3]])
    #########
    if name_sheet.find('scale')>=0:
      name='Time-Scale'
    else:
      name='R'

    if name_sheet.find('surplus')>=0:
      loc='Surplus'
    elif name_sheet.find('deficit')>=0:
      loc='Deficit'
    else:
      loc='Globe'
    #########
    count=0
    sns.set_theme(style="ticks", palette="deep", font_scale = 1.1)

    fig = plt.figure(figsize=(18, 10)) #34 12  
    ax  = plt.subplot(3,4,(1,2))

    x1=[x+1 for x in range(10) if typc[x]==1]
    y1=[permute_importances[x] for x in range(10) if typc[x]==1]
    x2=[x+1 for x in range(10) if typc[x]==2]
    y2=[permute_importances[x] for x in range(10) if typc[x]==2]
    x3=[x+1 for x in range(10) if typc[x]==3]
    y3=[permute_importances[x] for x in range(10) if typc[x]==3]
    ax.bar(x1,y1, label="Climate",color="blue")
    ax.bar(x2,y2, label="Environment",color="red")
    ax.bar(x3,y3, label="Plant",color="green")
    ax.set_title("%s $(%s$ $Zones)$"%(titl,loc), fontsize=14, fontweight='bold',fontname='Arial')
    ax.set_ylabel("Permutation Importances", fontsize=14, fontweight='bold',fontname='Arial')
    ax.set_xticks(range(1,11))
    ax.set_xticklabels(names,rotation=30, fontsize=12, fontweight='bold',fontname='Arial') 
    ax.legend(loc=0,title='Drivers types')  
    fig.text(0.22, 0.93,'Gini scores of 90%% train data: %.3f'%rf_sel.score(X_train_sel, y_train),fontsize=12, fontweight='bold',horizontalalignment='left')
    fig.text(0.22, 0.91,'Gini scores of 10%% test data: %.3f'%rf_sel.score(X_test_sel, y_test),fontsize=12, fontweight='bold',horizontalalignment='left')
    ######################################
    for i in features:
        count=count+1
        pdp = partial_dependence(rf_sel, X_train_sel, [i], kind="both",response_method='auto', method='brute', grid_resolution=50) 

        ax  = plt.subplot(3,4,count+2)
        plot_x = pd.Series(pdp['values'][0]).rename('x')
        plot_i = pdp['individual']
        plot_y = pdp['average'][0]
        tck = splrep(plot_x, plot_y, s=30)
        xnew = np.linspace(plot_x.min(),plot_x.max(),300)
        ynew = splev(xnew, tck, der=0)

        plot_df = pd.DataFrame(columns=['x','y']) 
        for a in plot_i[0]: 
            a2 = pd.Series(a)
            df_i = pd.concat([plot_x, a2.rename('y')], axis=1) 
            plot_df = pd.concat([plot_df,df_i],  ignore_index=True)

        sns.lineplot(data=plot_df, x="x", y="y", color='k', linewidth = 1.5, linestyle='--', alpha=0.6) 
        plt.plot(xnew, ynew, linewidth=2) 
        sns.rugplot(df[i].sample(100), x = i, height=.05, color='k', alpha = 0.3) 

        x_min = plot_x.min()-(plot_x.max()-plot_x.min())*0.1
        x_max = plot_x.max()+(plot_x.max()-plot_x.min())*0.1

        i=re.sub(r'(\s+)', r'\0$ $', i) 
        plt.title('Rank %d: %.3f'%(count,im[count-1]), fontsize=14, fontweight='bold',fontname='Arial')
        plt.ylabel('%s'%titl, fontsize=14, fontweight='bold',fontname='Arial')
        plt.xlabel('$%s$'%i, fontsize=14, fontweight='bold',fontname='Arial')
        plt.xlim(x_min,x_max)

    if th%2==1:
         fig.text(0.01, 0.98, '(a)', fontsize=16, fontweight='bold',fontname='Arial')
    else:
          fig.text(0.01, 0.98, '(b)', fontsize=16, fontweight='bold',fontname='Arial')

    fig.tight_layout() 


    fig.savefig('D:\Droughts_scales\figures\Regressor_30drivers_sig_Partial_dependence_decolinear_global_%s_%s.jpg'%(name,loc),dpi=350,bbox_inches='tight')
    print('Finished...'+name_sheet)

writer.save()
writer.close() 
print('Running time: %s Seconds'%(time.time()-start))
print('All Finished!')
    ########################
