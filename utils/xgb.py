from sklearn.metrics import accuracy_score, roc_auc_score, confusion_matrix, plot_confusion_matrix
from sklearn.preprocessing import LabelEncoder
from xgboost.sklearn import XGBClassifier
from sklearn.feature_selection import SelectFromModel
from matplotlib import pyplot as plt
from sklearn.model_selection import train_test_split
import pandas as pd
import xgboost as xgb
import numpy as np
import itertools

from sklearn import metrics
from sklearn.model_selection import GridSearchCV

umap4 = pd.read_csv(r"./data/UMAP.csv")
cancerType = ["ACC","BLCA","BRCA","CEAD","CESC","CHOL","COAD","DLBC","ESAD","ESSC",
              "GBM","HNSC","KICH","KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC",
              "MESO","OV","PAAD","PCPG","PRAD","READ","SARC","SKCM","STAD","TGCT",
              "THCA","THYM","UCEC","UCS","UVM"]
umap4 = umap4[umap4['Cancer'].isin(cancerType)]
umap4_data = umap4.iloc[:,0:umap4.shape[1]-1].values
umap4_target = umap4.iloc[:,umap4.shape[1]-1].values
Ly = LabelEncoder()
umap4_label = Ly.fit_transform(umap4_target)

Xtrain, Xtest, Ytrain, Ytest = train_test_split(umap4_data,umap4_label,test_size=0.1,random_state=666)
xgb_train = xgb.DMatrix(Xtrain,label = Ytrain)
xgb_test = xgb.DMatrix(Xtest,label = Ytest)

xgb_model = XGBClassifier(
                     learning_rate=0.1,
                     objective='multi:softmax',
                     num_class = 35,
                     nthread = 4,
                     n_estimators=400,
                     max_depth=5,
                     gamma=0,
                     min_child_weight=3,
                     subsample=0.7,
                     colsample_bytree=0.3,
                     reg_alpha=0.5,
                     early_stopping_rounds = 50,
                     seed=666)

def get_cv(xgb_model,xgb_train):
    cv_results = xgb.cv(xgb_model.get_xgb_params(),
                        xgb_train,
                        nfold=5,
                        num_boost_round=1000,
                        metrics='auc',
                        early_stopping_rounds=100,
                        seed=666)
    print(cv_results)

def get_grid_search(xgb_model,param_grid):
    grid_search=GridSearchCV(xgb_model,param_grid,scoring='accuracy',cv=5,iid=False)
    grid_search.fit(Xtrain,Ytrain)
    print('best_params:',grid_search.best_params_)
    print('best_score:',grid_search.best_score_)

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

def fea_selection(model):
    model.fit(umap4_data,umap4_label)
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
                         num_class = 4,
                         nthread = 4,
                         n_estimators= 399,
                         max_depth=8,
                         gamma=0,
                         min_child_weight=2,
                         subsample=0.8,
                         colsample_bytree=0.8,
                         reg_alpha=0.48,
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

def plot_confuse(xgb_model):
    xgb_model.fit(Xtrain, Ytrain)
    y_pred = xgb_model.predict(Xtest)
    cm = confusion_matrix(y_true=Ytest, y_pred=y_pred)
    classes = cancerType
    title='Confusion matrix'
    cmap=plt.cm.Blues
    cm = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis]
    plt.imshow(cm, interpolation='nearest', cmap=cmap)
    plt.title(title)
    plt.colorbar()
    tick_marks = np.arange(len(classes))
    plt.xticks(tick_marks, classes, rotation=45, fontsize = 6)
    plt.yticks(tick_marks, classes, fontsize = 6)
    thresh = cm.max() / 2.
    for i, j in itertools.product(range(cm.shape[0]), range(cm.shape[1])):
        plt.text(j, i, '{:.2f}'.format(cm[i, j]) if i == j else '',
                 horizontalalignment="center", verticalalignment = 'center',
                 color="red" if i == j else "white",fontsize = 5,weight = "bold")
    plt.tight_layout()
    plt.ylabel('True label')
    plt.xlabel('Predicted label')
    plt.savefig("confusion.png",dpi=300,figsize=(12,8))