## This code is adapted from Brian Rogers's benchmarks.ipynb notebook
## as well as from the PyOD documentation

## M. Ellie White ~ 11 Feb 2025

from pyod.models.knn import KNN
from pyod.models.iforest import IForest
from pyod.models.ecod import ECOD
from pyod.models.lof import LOF
from pyod.models.abod import ABOD
from pyod.models.cblof import CBLOF
from pyod.models.cof import COF
from pyod.models.cd import CD
from pyod.models.copod import COPOD
from pyod.models.feature_bagging import FeatureBagging
from pyod.models.hbos import HBOS
from pyod.models.inne import INNE
from pyod.models.kde import KDE
from pyod.models.loci import LOCI
from pyod.models.pca import PCA

def unsupervised_runsort(clf, clf_name, x_test, x_train, df_initial):

    ## this module runs an unsupervised ML algorithm on a 
    ## dataset, then sorts the output, ranking in order from
    ## most to least outlying data points. The function then
    ## returns a dataframe containing all of the information 
    ## about the given objects of interest from the original 
    ## simulation dataset, with the addition of a column for 
    ## anomaly scores. 
    
    subspace = ["a*", "i-z", "a", "sini", "e", "v-vk", "r"]
    X_train = x_train[subspace]
    X_test = x_test[subspace]
    
    clf.fit(X_train)

    # get the prediction labels and outlier scores of the training data
    y_train_pred = clf.labels_ # binary labels (0: inliers, 1: outliers)
    y_train_scores = clf.decision_scores_ # raw outlier scores
    
    # get the prediction on the test data
    y_test_pred = clf.predict(X_test) # outlier labels (0 or 1)
    y_test_scores = clf.decision_function(X_test) # outlier scores
    
    # it is possible to get the prediction confidence as well
    # outlier labels (0 or 1) and confidence in the range of [0,1]
    y_test_pred, y_test_pred_confidence = clf.predict(X_test, return_confidence=True) 
    
    X_test_modified = X_test
    
    X_test_modified['{}_score'.format(clf_name)] = y_test_scores

    ## sort the test samples from highest to lowest outlier scores
    xtest_sorted = X_test_modified.sort_values(by='{}_score'.format(clf_name), ascending=False) 

    ## retrieve the indices of each sample, in a list where the indices are sorted from those
    ## corresponding to the sample with the highest outlier scores to the lowest 
    indices = xtest_sorted.index.to_list() #retrieve indices of each sample

    # get the sorted version of the dataframe with *all* of the fields, not just the subset fields
    df_full_table_sorted = df_initial.loc[indices] 

    #drop duplicate "observations" of same object
    df_full_table_sorted = df_full_table_sorted.drop_duplicates(subset='ssObjectId', keep='first')

    return df_full_table_sorted
