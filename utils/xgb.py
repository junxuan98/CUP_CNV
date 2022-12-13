from sklearn.metrics import accuracy_score, roc_auc_score, confusion_matrix, plot_confusion_matrix
from sklearn.preprocessing import LabelEncoder
from xgboost.sklearn import XGBClassifier
from sklearn.feature_selection import SelectFromModel
from matplotlib import pyplot as plt
from sklearn.model_selection import train_test_split, StratifiedShuffleSplit, RepeatedKFold
import pandas as pd
import xgboost as xgb
import numpy as np
import itertools
from collections import Counter

from sklearn import metrics
from sklearn.model_selection import GridSearchCV

umap = pd.read_csv(r"../data/UMAP.csv")
#'UCS': 46, 'DLBC': 34, 'CHOL': 33, 'ESAD': 27, 'CEAD': 25
# cancerType = ['GBM','KICH','KIRC','KIRP','LUSC','PRAD','READ','TGCT','THCA','UVM']
cancerType = ["ACC","BLCA","BRCA","CESC","COAD","ESSC",
              "GBM","HNSC","KICH","KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC",
              "MESO","OV","PAAD","PCPG","PRAD","READ","SARC","SKCM","STAD","TGCT",
              "THCA","THYM","UCEC","UVM"]
#"UCS", "UVM" ,"MESO", "SKCM"
nec = ["GBM", "LAML", "LGG", "PCPG", "SARC", "TGCT", "THYM" ]
oec = ["ACC", "BLCA", "KICH", "KIRC", "KIRP", "LIHC"]
#, "CHOL","CEAD", "ESAD"
adc = ["BRCA", "COAD", "LUAD", "OV", "PAAD", "PRAD", "STAD", "THCA", "UCEC"]
#"ESSC"
scc = ["CESC","HNSC", "LUSC"]
classes = scc
type = 'SCC'
num = 6

umap = umap[umap['Cancer'].isin(cancerType)]
umap_data = umap.iloc[:,0:umap.shape[1]-1].values
umap_target = umap.iloc[:,umap.shape[1]-1].values
print(Counter(umap_target))

Ly = LabelEncoder()
umap_label = Ly.fit_transform(umap_target)
# Xtrain, Xtest, Ytrain, Ytest = train_test_split(umap_data,umap_label,test_size=0.1,random_state=66)
# xgb_train = xgb.DMatrix(Xtrain,label = Ytrain)
# xgb_test = xgb.DMatrix(Xtest,label = Ytest)
ss = StratifiedShuffleSplit(n_splits=1,test_size=0.2,random_state=42)
for train_index, test_index in ss.split(umap_data,umap_target):
    Xtrain = umap_data[train_index]
    Xtest = umap_data[test_index]
    Ytrain = umap_label[train_index]
    Ytest = umap_label[test_index]
xgb_train = xgb.DMatrix(Xtrain,label = Ytrain)
xgb_test = xgb.DMatrix(Xtest,label = Ytest)


xgb_model = XGBClassifier(
                     learning_rate=0.1,
                     objective='multi:softmax',
                     num_class = num,
                     nthread = 4,
                     n_estimators=454,
                     max_depth=5,
                     gamma=0,
                     min_child_weight=2,
                     subsample=0.5,
                     colsample_bytree=0.7,
                     reg_alpha=0.3,
                     seed=666)

def get_param(xgb_model,xgb_train):
    rkf=RepeatedKFold(n_splits=10,n_repeats=5,random_state=88)#设置分割策略
    cv_result = xgb.cv(xgb_model.get_xgb_params(),
                       xgb_train,
                       num_boost_round=xgb_model.get_params()['n_estimators'],
                       folds=rkf,
                       metrics='acc',
                       early_stopping_rounds=50,
                       callbacks=[xgb.callback.early_stop(50),
                                  xgb.callback.print_evaluation(period=1,show_stdv=True)])
    print(cv_result)

# get_param(xgb_model,xgb_train)

def get_cv(xgb_model,xgb_train):
    cv_results = xgb.cv(xgb_model.get_params(),
                        xgb_train,
                        nfold=5,
                        num_boost_round=100,
                        metrics='auc',
                        early_stopping_rounds=100,
                        seed=666)
    print(cv_results)

# param_grid={'reg_alpha':[i/10 for i in range(0,10)]}

def get_grid_search(xgb_model,param_grid):
    grid_search=GridSearchCV(xgb_model,param_grid,scoring='accuracy',cv=5,iid=False)
    grid_search.fit(Xtrain,Ytrain)
    print('best_params:',grid_search.best_params_)
    print('best_score:',grid_search.best_score_)

# get_grid_search(xgb_model,param_grid)

def multiclass_metrics(xgb_model):
    xgb_model.fit(Xtrain,Ytrain)
    y_prediction = xgb_model.predict(Xtest)
    print(y_prediction)
    y_predprob = xgb_model.predict_proba(Xtest)
    print(Ytest)
    print(y_predprob)
    print("Accuracy (Test):",round(accuracy_score(y_true=Ytest, y_pred=y_prediction), 4))
    print("Recall Score (Test):",round(metrics.recall_score(Ytest,y_prediction,average='macro'),4))
    print("Prediction Score (Test):",round(metrics.precision_score(Ytest,y_prediction,average='macro'),4))
    print("F1 Score (Test):",round(metrics.f1_score(Ytest,y_prediction,average='macro'),4))
    print("ROC_AUC Score (Test,ovr):",round(roc_auc_score(Ytest,y_predprob,multi_class='ovr'),4))
    print("ROC_AUC (Test,ovo):",round(roc_auc_score(Ytest,y_predprob,multi_class='ovo'),4))
multiclass_metrics(xgb_model)

def fea_selection(model):
    model.fit(umap_data,umap_label)
    xgb.plot_importance(model)
    plt.show()
    thresholds = np.sort(model.feature_importances_)
    print(model.feature_importances_)
    y_prediction = xgb_model.predict(Xtest)
    y_predprob = xgb_model.predict_proba(Xtest)
    for threshold in thresholds:
        # select features using threshold
        selection = SelectFromModel(model, threshold=threshold, prefit=True)
        select_X_train = selection.transform(Xtrain)
        # train model
        selection_model = XGBClassifier(learning_rate=0.1,
                     objective='multi:softmax',
                     num_class = num,
                     nthread = 4,
                     n_estimators=454,
                     max_depth=5,
                     gamma=0,
                     min_child_weight=2,
                     subsample=0.5,
                     colsample_bytree=0.7,
                     reg_alpha=0.3,
                     seed=666)
        selection_model.fit(select_X_train, Ytrain)
        # eval model
        select_X_test = selection.transform(Xtest)
        y_pred = selection_model.predict(select_X_test)
        y_predprob = selection_model.predict_proba(select_X_test)
        predictions = [round(value) for value in y_pred]
        accuracy = accuracy_score(Ytest, predictions)
        roc_auc_ovo = roc_auc_score(Ytest,y_predprob,multi_class='ovo')
        roc_auc_ovr = roc_auc_score(Ytest, y_predprob, multi_class='ovr')
        print("Thresh=%.3f, n=%d, Accuracy: %.2f%%" % (threshold, select_X_train.shape[1], accuracy * 100.0))
        print("Thresh=%.3f, n=%d, roc_auc_ovo: %.2f%%" % (threshold, select_X_train.shape[1], roc_auc_ovo * 100.0))
        print("Thresh=%.3f, n=%d, roc_auc_ovr: %.2f%%" % (threshold, select_X_train.shape[1], roc_auc_ovr * 100.0))


def plot_confuse(xgb_model,type):
    xgb_model.fit(Xtrain, Ytrain)
    y_pred = xgb_model.predict(Xtest)
    cm = confusion_matrix(y_true=Ytest, y_pred=y_pred)
    cm = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis]
    cmap = plt.cm.Blues
    plt.imshow(cm, interpolation='nearest', cmap=cmap)
    plt.title(type,fontsize = 14, fontweight = 'heavy')
    plt.colorbar()
    tick_marks = np.arange(len(classes))
    plt.xticks(tick_marks, classes, rotation=45, fontsize = 10)
    plt.yticks(tick_marks, classes, fontsize = 10)
    # thresh = cm.max() / 2.
    for i, j in itertools.product(range(cm.shape[0]), range(cm.shape[1])):
        plt.text(j, i, '{:.2f}'.format(cm[i, j]) if i == j else '',
                 horizontalalignment="center", verticalalignment = 'center',
                 color="white" if i == j else "white",fontsize = 10,weight = "bold")
    plt.tight_layout()
    plt.ylabel('True label',fontsize = 14)
    plt.xlabel('Predicted label', fontsize = 14)
    plt.savefig("newconfusion"+type+".png",dpi=300,figsize=(12,6),bbox_inches='tight')

plot_confuse(xgb_model,type)