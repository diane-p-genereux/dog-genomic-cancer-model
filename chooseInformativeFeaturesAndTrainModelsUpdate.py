'''
Editing to run on master file from code directory. Previously, the input file already had folds assigned. So, need to add that here,
as well as a way of which samples are in each fold so subsequent runs and downstream analysis maintain this. This script selects features and trains
a model for only one fold.
KJM update 03/2025: edited to allow selection of more subsets, updated master file path.
'''


#source activate /seq/vgb/swofford/conda/my3env <- must be active before running.
#qsub -V -cwd -b y -pe smp 10 -binding linear:10 -l h_vmem=5G -l h_rt=240:00:00 python chooseInformativeFeatures.py -o /seq/vgb/swofford/R01/convertingToSKLInput/testVCFs/sigOnlyDog101422/ -s 100 -rfco 20 -rfci 3 -xgbo 20 -xgbi 3 -tr dog -f 4
import pandas as pd
import numpy as np
import xgboost as xgb
from sklearn.utils import class_weight
import pickle
from sklearn.feature_selection import RFECV
import argparse
import os
from sklearn.metrics import accuracy_score
from sklearn.metrics import balanced_accuracy_score
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import GridSearchCV
from bayes_opt import BayesianOptimization
from sklearn.model_selection import cross_val_score
import time


parser = argparse.ArgumentParser()
parser.add_argument('-re', '--reRun', action = 'store_true', default = False) #rerun with only the feature types available to original classifier 
parser.add_argument('-i', '--inFile', default = '/seq/vgb/cancer_r01/ML_annotations/remove_duplicates/allxall_master_features_samples_duplicates_removed.txt.gz') #allxall master file
parser.add_argument('-f', '--fold') #which fold to use as hold out
parser.add_argument('-rhm', '--removeHighMut', action='store_true', default=False) #remove human cancer types defined as high mutation rate
parser.add_argument('-gpsm', '--genesPathsSigsMutCount', action='store_true', default=False) #Select one subset only.
parser.add_argument('-nt', '--noTrinuc', action='store_true', default=False)
parser.add_argument('-so', '--sigsOnly', action='store_true', default=False)
parser.add_argument('-do', '--driversOnly', action='store_true', default=False)
parser.add_argument('-o', '--outDir') #output directory
parser.add_argument('-s', '--stepSize', default = 100) #number of features to eliminate per round of random forest feature selection
parser.add_argument('-rfco', '--rfcOptRounds', default = 20) #number of optimization rounds for random forest classifier
parser.add_argument('-rfci', '--rfcInitPoints', default = 3) #number of initiation points for random forest optimizer
parser.add_argument('-xgbo', '--xgbOptRounds', default = 20) #number of optimization rounds for xgboost classifier
parser.add_argument('-xgbi', '--xgbInitPoints', default = 3) #number of initiation points for xgboost optimizer
parser.add_argument('-tr', '--trainOn', nargs = '+', default = ['human']) #adding option to choose which (or both) species to train on
args = parser.parse_args()

'''
After reading about feature selection pitfalls, it's clear that
feature selection should be performed on a per-fold level. xgboost is
prohibitively slow with the number of features we're using for this.
The below is a proposal to use a standard random forest classifier as the
feature selection method. The entire human data set is split into 5 folds.
For each fold, recursive feature selection is performed using randomforest as
the classifier, and one of the folds as a holdout set. The recursive feature
selection itself has 5 recursive holdouts of the subset train data. Currently
using default values to see how quickly it will run, and how many features per fold
will be selected. Subsequently could run xgboost on each of the folds with the selected
feature set for that fold vs the test set for that fold for an average fold accuracy.
k-folds is to estimate real-world accuracy, not to make the best model possible,
after the accuracy is estimated, the full process of feature selection and training
is performed on the full human data, and that is what is used to predict dog values.
Alternatively, each fold could predict dog values and we could look at an average.
'''

def rfc_cv(data, label, min_samples_split, max_features, criterion, max_depth, n_estimators):  
    estimator = RandomForestClassifier(n_jobs=-1, 
                                       class_weight = "balanced",
                                       max_depth = max_depth,
                                       min_samples_split = min_samples_split,
                                       criterion = criterion, 
                                       max_features = max_features,
                                       n_estimators = n_estimators
                                       )
    results = cross_val_score(estimator, data, label, scoring = 'neg_log_loss', cv = 3)
    return results.mean()
    
                          
'''
optimization function for random forest. models are iteratively trained with varying parameter
values using bayesian optimization to identify optimal parameters from limited testing.
'''                                  
def rfc_optimize(data, label):
    data = data
    label = label
    def rfc_crossval(min_samples_split, max_features, criterion, max_depth, n_estimators):
        if int(criterion) == 0:
            criterion = 'gini'
        elif int(criterion) == 1:
            criterion = 'entropy'
        if int(max_features) == 0:
            max_features = 'sqrt'
        elif int(max_features) == 1:
            max_features = 'log2'
        return rfc_cv(data,
                      label,
                      min_samples_split = int(min_samples_split),
                      max_features = max_features,
                      criterion = criterion, 
                      max_depth=int(max_depth),
                      n_estimators = int(n_estimators) 
                      )
    
    optimizer = BayesianOptimization(
        f = rfc_crossval,      
        pbounds = {
            'min_samples_split': (2, 25),
            'max_features':(0,1.99),
            'criterion':(0,1.99),
            'max_depth': (1, 15),
            'n_estimators': (200, 1000)
        },
        random_state=42,
        verbose = 10
    )
    #optimizer.maximize(n_iter=20, init_points=3)
    #optimizer.maximize(n_iter=60, init_points=10)
    optimizer.maximize(n_iter=int(args.rfcOptRounds), init_points=int(args.rfcInitPoints))
    min_samples_split = int(optimizer.max['params']['min_samples_split'])
    if int(optimizer.max['params']['max_features']) == 0:
            max_features = 'sqrt'
    elif int(optimizer.max['params']['max_features']) == 1:
            max_features = 'log2'
    if int(optimizer.max['params']['criterion']) == 0:
            criterion = 'gini'
    elif int(optimizer.max['params']['criterion']) == 1:
            criterion = 'entropy'
    max_depth = int(optimizer.max['params']['max_depth'])
    n_estimators = int(optimizer.max['params']['n_estimators'])
    
    params = {
        'min_samples_split': min_samples_split,
        'max_features': max_features,
        'criterion': criterion,
        'max_depth': max_depth,
        'n_estimators': n_estimators
        }
    return(params)


'''
uses a random forest classifier with parameters determined via parameter search
to rank and choose features for inclusion in the xgboost model.
'''
def chooseFeatures(data, label, params):
    clf = RandomForestClassifier(n_jobs = -1, class_weight = 'balanced', 
                                 max_features = params['max_features'], 
                                 n_estimators = params['n_estimators'], 
                                 max_depth = params['max_depth'], 
                                 criterion = params['criterion']
                                 )
    #selector = RFECV(clf, cv=5, step=200, verbose=10, min_features_to_select=400)
    selector = RFECV(clf, cv=5, step=int(args.stepSize), verbose=10, min_features_to_select=400)
    selector.fit(data, label)
    importantFeatures = []
    for x, val in enumerate(selector.support_):
        if val:
            importantFeatures.append(features[x])
    return(importantFeatures)

'''
for parameter tuning. Creates and trains a 3-fold xgboost.cv instance returning the mean mlogloss
of the final iteration allowing for variable input values to be tested.
''' 
def xgb_cv(data, label, max_depth, learning_rate, gamma, reg_alpha, 
           reg_lambda, colsample_bytree, min_child_weight):
    dtrain = xgb.DMatrix(data=data, label=label)
    params = {
        'objective':'multi:softprob', 
        'random_state':42,
        'subsample':.8, 
        'nthread':-1, 
        'max_depth':max_depth,
        'num_class':len(pd.unique(dtrain.get_label())), 
        'learning_rate':learning_rate,
        'gamma':gamma, 
        'reg_alpha':reg_alpha, 
        'reg_lambda':reg_lambda, 
        'colsample_bytree':colsample_bytree, 
        'min_child_weight':min_child_weight
        }
    cv_results = xgb.cv(dtrain = dtrain, params=params, nfold = 3, 
                        num_boost_round = 1000, early_stopping_rounds = 10, 
                        metrics = 'mlogloss', as_pandas=True, seed=42)
    results = float(cv_results['test-mlogloss-mean'][-1:])
    #lower mlogloss is better, so we return the negative value so the optimizer
    #optimizes in the correct direction.
    return (-1 * results)
                          
'''
optimization function for xgboost. models are iteratively trained with varying parameter
values using bayesian optimization to identify optimal parameters from limited testing.
'''                                  
def xgb_optimize(data, label):
    data = data
    label = label
    #def xgb_crossval(max_depth=6, learning_rate=.3, gamma=0, reg_alpha=0, 
    #                 reg_lambda=1, colsample_bytree=1, min_child_weight=1):
    def xgb_crossval(max_depth='', learning_rate='', gamma='', reg_alpha=0, 
                     reg_lambda='', colsample_bytree='', min_child_weight=''):
        return xgb_cv(data,
                      label, 
                      max_depth=int(max_depth), 
                      learning_rate = learning_rate,
                      gamma = gamma, 
                      reg_alpha = reg_alpha, 
                      reg_lambda = reg_lambda, 
                      colsample_bytree = colsample_bytree, 
                      min_child_weight = min_child_weight
                      )
    
    optimizer = BayesianOptimization(
        f = xgb_crossval,
        pbounds={
            'max_depth': (3,18),#18
            'learning_rate':(.01,1),
            'gamma': (1,9),
            #'reg_alpha' : (0, 180),
            'reg_lambda': (0, 100),
            'colsample_bytree': (0.5, 1),
            'min_child_weight': (0, 10),
        },
        random_state=42,
        verbose = 10
    )
    #optimizer.maximize(n_iter=20, init_points=3)
    #optimizer.maximize(n_iter=60, init_points=10)
    optimizer.maximize(n_iter=int(args.xgbOptRounds), init_points=int(args.xgbInitPoints))
    
    max_depth = int(optimizer.max['params']['max_depth'])
    learning_rate = optimizer.max['params']['learning_rate']
    gamma = optimizer.max['params']['gamma']
    reg_lambda = optimizer.max['params']['reg_lambda']
    colsample_bytree = optimizer.max['params']['colsample_bytree']
    min_child_weight = optimizer.max['params']['min_child_weight']
    
    params = {
        'objective':'multi:softprob', 
        'subsample':.8, 
        'nthread':-1, 
        'max_depth':max_depth,
        'learning_rate':learning_rate,
        'gamma':gamma, 
        'reg_lambda':reg_lambda, 
        'colsample_bytree':colsample_bytree, 
        'min_child_weight':min_child_weight
        }
    return(params)

def trainModel(data, label, test, testLabel, params):
    data = data
    label = label
    xgb_model = xgb.XGBClassifier(objective='multi:softprob', random_state=42, 
                              subsample = .8, nthread = -1, max_depth = params['max_depth'], 
                              n_estimators = 3000, learning_rate = params['learning_rate'],
                              gamma = params['gamma'], reg_lambda = params['reg_lambda'],
                              colsample_bytree = params['colsample_bytree'], 
                              min_child_weight = params['min_child_weight'])
    eval_set = [(data, label), (test, testLabel)]
    xgb_model.fit(data, label, early_stopping_rounds=10,
                  eval_metric=["merror", "mlogloss"], eval_set=eval_set, verbose=True, 
                  sample_weight=class_weight.compute_sample_weight("balanced", label))
    return xgb_model

if not os.path.exists(args.outDir):
    os.mkdir(args.outDir)


df = pd.read_csv(args.inFile, sep = '\t')
#drop any columns with na values
#df.dropna(axis = 0, inplace = True)
df.fillna(0, inplace = True)

#Drop high mutation rate cancers if option is selected
if args.removeHighMut:
    highMutList=["BLCA","CESC","CORE","DLBCL","ESCA","HNSC","KIRC","KIRP","LICA","LUNG","MELA","ORCA","OV","STAD","UCEC"]
    filtered_df = df.loc[~df['label'].isin(highMutList)]
    df = filtered_df
df.sort_values(['label','species'], inplace = True)
df.set_index('ID', inplace = True)
labels = pd.unique(df['label'])
#df['factorLabel'] = pd.factorize(df['label'])[0]

'''
check to see if folds have been assigned previously, and if not assign them and save fold/ID mapping to file
'''
try:
    folds = pd.read_csv(args.outDir + '/fold_assignments.txt', sep = '\t')
    folds.set_index('ID', inplace = True)
    df['is_train'] = folds['is_train']
    df['factorLabel'] = folds['factorLabel']
except FileNotFoundError:
    foldBalance = [0]
    #repeat until all folds contain all cancer classes
    while np.min(foldBalance) != len(labels):    
        foldBalance = []
        df['is_train'] = np.random.randint(0,5, len(df))
        for i in range(5):
            foldBalance.append(len(pd.unique(df[df['is_train'] == i]['label'])))
    df['factorLabel'] = pd.factorize(df['label'])[0]
    df[['is_train','factorLabel']].to_csv(args.outDir + '/fold_assignments.txt', sep = '\t', index = True)
            
            
    
    
    
features = list(df.columns[2:])
#drop any columns that sum to less than 5, arbitrary filter to minimize number of likely uniformative features.
#features = df[features][df[features].columns[df[features.sum()>5]]].columns
try:
    features.remove('is_train')
except ValueError:
    pass
try:
    features.remove('factorLabel')
except ValueError:
    pass
try:
    features.remove('label.1')
except ValueError:
    pass
try:
    features.remove('species.1')
except ValueError:
    pass

if args.reRun:
    keep = []
    for i in features:
        if 'chr' and '@' in i:
            continue
        elif 'amag' in i:
            continue
        elif '>' in i:
            continue
        elif 'lower' in i:
            continue
        else:
            keep.append(i)
    features = keep

# Subset features. Note that options cannot be combined - select one subset only.
if args.noTrinuc: # Skip the trinucleotide contexts, which have ">" in the name
    features = features[96:] #dropping trinucleotide contexts and id, label, species
elif args.genesPathsSigsMutCount: #Select features included in the previous "full" models
    features = features[96:40290]
elif args.sigsOnly:
    features = features[40247:40290]
elif args.driversOnly:
    keepFeatures=[]
    binary = features[96:15411]
    with open('/seq/vgb/cancer_r01/ML_annotations/driver_genes.txt', 'r') as d:
        for line in d:
            driverList.add(line.strip())
    for i in binary:
        if i in driverList:
            keepFeatures.append(i)
            features = keepFeatures # Set features to the subset of selected features

labelAdjust = 0
#define training and test sets based on inputs
if args.trainOn == ['human']:
    dogs = df[df['species'] == 'dog']
    train = df[df['is_train'] != int(args.fold)]
    test = df[df['is_train'] == int(args.fold)]
    train = train[train['species'] == 'Human']
    test = test[test['species'] == 'Human']
elif args.trainOn == ['dog']:
    train = df[df['is_train'] != int(args.fold)]
    train = train[train['species'] == 'dog']
    test = df[df['is_train'] == int(args.fold)]
    test = test[test['species'] == 'dog']
    labelAdjust = min(train['factorLabel'])
    train['factorLabel'] = train['factorLabel'] - labelAdjust
    test['factorLabel'] = test['factorLabel'] - labelAdjust
    print('here')
elif 'dog' in args.trainOn and 'human' in args.trainOn:
    train = df[df['is_train'] != int(args.fold)]
    test = df[df['is_train'] == int(args.fold)]

#run steps, if step hasn't previously run        
if not os.path.exists(args.outDir + '/rfc_params_fold_' + str(args.fold) + '.pickle'):
    rfcParams = rfc_optimize(train[features].copy(deep = True), train['factorLabel'])
    with open(args.outDir + '/rfc_params_fold_' + str(args.fold) + '.pickle', 'wb') as outFile:
        pickle.dump(rfcParams, outFile)
else:
    with open(args.outDir + '/rfc_params_fold_' + str(args.fold) + '.pickle', 'rb') as rfcInFile:
        rfcParams = pickle.load(rfcInFile)
if not os.path.exists(args.outDir + '/importantFeatures_fold_' + str(args.fold) + '.pickle'):
    importantFeatures = chooseFeatures(train[features],train['factorLabel'], rfcParams)
    with open(args.outDir + '/importantFeatures_fold_' + str(args.fold) + '.pickle', 'wb') as outFile:
        pickle.dump(importantFeatures, outFile)
else:
    with open(args.outDir + '/importantFeatures_fold_' + str(args.fold) + '.pickle', 'rb') as featureFile:
             importantFeatures = pickle.load(featureFile)
if not os.path.exists(args.outDir + '/xgb_params_fold_' + str(args.fold) + '.pickle'):             
    xgbParams = xgb_optimize(train[importantFeatures], train['factorLabel'])
    with open(args.outDir + '/xgb_params_fold_' + str(args.fold) + '.pickle', 'wb') as outFile:
        pickle.dump(xgbParams, outFile)
else:
    with open(args.outDir + '/xgb_params_fold_' + str(args.fold) + '.pickle', 'rb') as xgbParamFile:
        xgbParams = pickle.load(xgbParamFile)
xgb_model = trainModel(train[importantFeatures], train['factorLabel'], test[importantFeatures], test['factorLabel'], xgbParams)
with open(args.outDir + '/xgb_model_fold_' + str(args.fold) + '.pickle', 'wb') as outFile:
    pickle.dump(xgb_model, outFile)
xgb_model.save_model(args.outDir + '/xgb_model_fold_' + str(args.fold) + '.json')
    
'''        
#from xgboost import XGBClassifier
keep = {}
preds = {}
probs = {}
df['is_train'] = np.random.randint(5, size=len(df))
for i in range(5):
    train, test = df[df['is_train']!= i], df[df['is_train']==i]
    train = train[train['species'] != 'dog']
    test = test[test['species'] != 'dog']
    y_train = train['factorLabel']
    y_test = test['factorLabel']
    clf = RandomForestClassifier(n_jobs=-1, class_weight = "balanced")
    param_grid = { 
    'n_estimators': [200, 500],
    'max_features': ['sqrt', 'log2'],
    'max_depth' : [1,15],
    'criterion' :['gini', 'entropy']
    }
    griddy = GridSearchCV(estimator=clf, param_grid=param_grid, cv= 5, verbose = 10)
    griddy.fit(train[features], y_train)
    crit = griddy.best_params_['criterion'] #entropy
    nEst = griddy.best_params_['n_estimators'] #500
    maxFeat = griddy.best_params_['max_features'] #auto
    maxDepth = griddy.best_params_['max_depth'] #8
    clf = RandomForestClassifier(n_jobs = -1, class_weight = 'balanced', 
                                 max_features = maxFeat, 
                                 n_estimators = nEst, 
                                 max_depth = maxDepth, criterion = crit)
                
    selector = RFECV(clf, cv=5, step=200, verbose=10, min_features_to_select=400)
    selector.fit(train[features], y_train)
    importantFeatures = []
    for x, val in enumerate(selector.support_):
        if val:
            importantFeatures.append(features[x])
    keep[i] = importantFeatures
    
    params = {
        'min_child_weight': [1, 5, 10],
        'gamma': [0.5, 1, 1.5, 2, 5],
        'subsample': [0.6, 0.8, 1.0],
        'colsample_bytree': [0.6, 0.8, 1.0],
        'max_depth': [3, 4, 5]
        }
    
    
    xgb_model = xgb.XGBClassifier(objective='multi:softprob', random_state=42, 
                              subsample = .8, nthread = -1, max_depth = 10, 
                              n_estimators = 3000, learning_rate = 0.05)
    eval_set = [(train[importantFeatures], y_train), (test[importantFeatures], y_test)]
    xgb_model.fit(train[importantFeatures], y_train, early_stopping_rounds=10,
                  eval_metric=["merror", "mlogloss"], eval_set=eval_set, verbose=True, 
                  sample_weight=class_weight.compute_sample_weight("balanced", y_train))
    preds[i] = xgb_model.predict(test[importantFeatures])
    probs[i] = xgb_model.predict_prob(test[importantFeatures])

plot_confusion_matrix(xgb_model, test[importantFeatures], y_test, normalize = 'true', 
                      xticks_rotation = 'vertical', display_labels = labels, include_values = False)


start = time.time()
importantFeatures = chooseFeatures(train[features],train['factorLabel'], rfcParams)
end = time.time()
print(end - start)

class CustomModel(keras.Model):
    def train_step(self, data):
        if len(data) == 3:
            x, y, sample_weight = data
        else:
            sample_weight = None
            x, y = data
        with tf.GradientTape() as tape:
            y_pred = self(x, training=True)  # Forward pass
            loss = self.compiled_loss(
                y,
                y_pred,
                sample_weight=sample_weight,
                regularization_losses=self.losses
            )
        trainable_vars = self.trainable_variables
        gradients = tape.gradient(loss, trainable_vars)
        self.optimizer.apply_gradients(zip(gradients, trainable_vars))
        self.compiled_metrics.update_state(y, y_pred, sample_weight=sample_weight)
        return {m.name: m.result() for m in self.metrics}   

shapley export
for i in range(5):
    for x in range(labels):
        name = labels[x]
        tempDF = pd.DataFrame(shap_vals[i][x], columns=features)
#        tempDF.columns = features[i]
        cols = np.array(tempDF.columns)
        tempDF['ID'] = df['ID']
        tempDF['label'] = df['label']
        tempDF['is_train'] = df['is_train']
        cols = np.append(['ID','label','is_train'], cols)
        tempDF[cols].to_csv(outDir + 'fold_' + str(i) + '/' + name + '.txt', sep = '\t', index = False)
 
'''
