import warnings

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sklearn

from sklearn.metrics import explained_variance_score, mean_absolute_error, mean_squared_error, r2_score
from sklearn.model_selection import cross_val_score, GridSearchCV, KFold, learning_curve, LearningCurveDisplay, RandomizedSearchCV, RepeatedKFold, ShuffleSplit, train_test_split
from sklearn.preprocessing import MinMaxScaler, PolynomialFeatures, PowerTransformer, StandardScaler


DPI = None
NUM_JOBS = 4
RANDOM_SEED = 42
TEST_SET_FRACTION = 0.2
PLOT_COLOURS = ['#212121', '#3F51B5', '#303F9F', '#FF5252', '#D32F2F']  # black, blue, deep blue, red, deep red from https://www.materialpalette.com
mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color=['#303F9F', '#FF5252', '#D32F2F'])


def rmLowVarFeats(featsDF, varThresh=0.0, verbose=False):
    if verbose:
        print(f"Removing the features with null values...")
        print(f"  Original number of features: {len(featsDF.columns)}")
    featsToDrop = []
    for feat in featsDF.columns:
        if featsDF[feat].isnull().any():
            featsToDrop.append(feat)
    featsDF.drop(featsToDrop, axis=1, inplace=True)
    if verbose:
        print(f"  Total number of features left: {len(featsDF.columns)}\n")
    
    if verbose:
        print(f"Removing the features with variance below {varThresh:.2f}...")
        print(f"  Original number of features: {len(featsDF.columns)}")
    featsToDrop = []
    for feat in featsDF.columns:
        if featsDF[feat].var() <= varThresh:
            if verbose:
                print(f"    {feat}:    {featsDF[feat].var():.3f}")
            featsToDrop.append(feat)
    featsNoLowVarDF = featsDF.drop(featsToDrop, axis=1, inplace=False)
    if verbose:
        print(f"  Total number of features left: {len(featsNoLowVarDF.columns)}\n")
    return featsNoLowVarDF


def rmHighCorrFeats(featsDF, corrThresh=0.9, verbose=False):
    if verbose:
        print(f"Removing the second feature from every pair of features with correlation above {corrThresh:.2f}...")
        print(f"  Original number of features: {len(featsDF.columns)}")
    corrMat = featsDF.corr()
    featsToDrop = []
    for (i, feat1) in enumerate(featsDF.columns):
        for (j, feat2) in enumerate(featsDF.columns):
            if feat1 == feat2: continue
            if j <= i: continue
            if abs(corrMat.loc[feat1][feat2]) >= corrThresh:
                if feat1 in featsToDrop or feat2 in featsToDrop:
                    continue
                print(f"    {feat1} {feat2}:    {abs(corrMat.loc[feat1][feat2]):.3f}")
                featsToDrop.append(feat2)
    featsNoHighCorrDF = featsDF.drop(featsToDrop, axis=1, inplace=False)
    if verbose:
        print(f"  Total number of features left: {len(featsNoHighCorrDF.columns)}\n")
    return featsNoHighCorrDF


def splitScaleData(df, polyFeats=False, scaling='minmax', verbose=False):
    X, y = df.iloc[:, :-1], df.iloc[:, -1]
    if polyFeats:
        X = PolynomialFeatures(degree=2, include_bias=False).fit_transform(X)  # Polynomial features (optional)
        XNoLowVar = rmLowVarFeats(X, varThresh=0.0, verbose=verbose)
        X = rmHighCorrFeats(XNoLowVar, corrThresh=0.9, verbose=verbose)
    Xtrain, Xtest, yTrain, yTest = train_test_split(X, y, test_size=TEST_SET_FRACTION, random_state=RANDOM_SEED, shuffle=True, stratify=None)  # Train-test split

    if scaling:
        if verbose:
            print(f"Applying {scaling} scaling to feature values...")
        if scaling == 'minmax':
            scaler = MinMaxScaler() 
        elif scaling == 'standard':
            scaler = StandardScaler()
        elif scaling == 'yeo-johnson':
            scaler = PowerTransformer(method='yeo-johnson', standardize=True)
        else:
            raise Exception(f"Unknown {scaling} scaling specified!")
        Xtrain = scaler.fit_transform(Xtrain)
        Xtrain = pd.DataFrame(Xtrain, columns=X.columns)
        Xtest = scaler.transform(Xtest)
        Xtest = pd.DataFrame(Xtest, columns=X.columns)
    return Xtrain, Xtest, yTrain, yTest


def evalModel(regModel, Xtrain, yTrain, Xtest, yTest, metric='coefficient of determination'):
    regModel.fit(Xtrain, yTrain)
    yPredTrain = regModel.predict(Xtrain)
    yPredTest = regModel.predict(Xtest)
    print(f"  Trained and tested using {metric}.")
    print(f"    Training score: {regModel.score(Xtrain, yTrain):.3f}")
    print(f"    Testing score:")  # {regModel.score(Xtest, yTest):.3f}
    print(f"      Mean absolute error: {mean_absolute_error(yTest, yPredTest):.6f}")
    print(f"      Mean squared error: {mean_squared_error(yTest, yPredTest, squared=True):.6f}")
    print(f"      Root mean squared error: {mean_squared_error(yTest, yPredTest, squared=False):.6f}")
    print(f"      Coefficient of determination: {r2_score(yTest, yPredTest):.3f}")
    return yPredTrain, yPredTest, regModel


def plotR2(yTrain, yTest, yPredTrain, yPredTest, modelName, 
           dataColour=PLOT_COLOURS[1], bestFitColour=PLOT_COLOURS[3], bestFitLW=1, perfectFitLW=1.5, 
           markerSize=20, markerAlpha=1.0, markerLW=0.7, 
           figSize=(7, 3), figName='figures/DBoxR2.png', showFig=False):
    fig, axes = plt.subplots(1, 2, figsize=figSize, dpi=DPI)

    perfectFitTrain = np.linspace(yTrain.min(), yTrain.max(), num=100)
    mTrain, cTrain = np.polyfit(yTrain, yPredTrain, deg=1)
    axes[0].scatter(yTrain, yPredTrain, color=dataColour, s=markerSize, alpha=markerAlpha, lw=markerLW, edgecolors=PLOT_COLOURS[0], zorder=2)
    axes[0].plot(perfectFitTrain, cTrain + mTrain*perfectFitTrain, color=bestFitColour, linestyle='--', lw=bestFitLW, zorder=3)
    axes[0].plot(perfectFitTrain, perfectFitTrain, color=PLOT_COLOURS[0], lw=perfectFitLW, zorder=1)
    axes[0].grid(linestyle='dotted')
    axes[0].legend(['Training Data', 'Best Fit', 'Perfect Fit'])

    perfectFitTest = np.linspace(yTest.min(), yTest.max(), num=100)
    mTest, cTest = np.polyfit(yTest, yPredTest, deg=1)
    axes[1].scatter(yTest, yPredTest, color=dataColour, s=markerSize, alpha=markerAlpha, lw=markerLW, edgecolors=PLOT_COLOURS[0],  zorder=2)
    axes[1].plot(perfectFitTest, cTest + mTest*perfectFitTest, color=bestFitColour, linestyle='--', lw=bestFitLW, zorder=3)
    axes[1].plot(perfectFitTest, perfectFitTest, color=PLOT_COLOURS[0], lw=perfectFitLW, zorder=1)
    axes[1].grid(linestyle='dotted')
    axes[1].legend(['Testing Data', 'Best Fit', 'Perfect Fit'])

    fig.text(0.53, 0.0, r'Actual $D_{B}$', ha='center')
    fig.text(0.0, 0.4, r'Predicted $D_{B}$', ha='center', rotation=90)
    plt.tight_layout()
    plt.savefig(figName, bbox_inches='tight')
    if showFig:
        plt.show()


def plotR2All(yTrain, yTest, yPredTrain, yPredTest, modelName, 
              trainDataColour=PLOT_COLOURS[1], trainBestFitColour=PLOT_COLOURS[2], testDataColour=PLOT_COLOURS[3], testBestFitColour=PLOT_COLOURS[4], 
              bestFitLW=1.8, perfectFitLW=1.5, 
              markerSize=25, markerAlpha=1.0, markerLW=0.9, 
              figSize=(4.5, 3.5), figName='figures/DBoxR2All.png', showFig=False):
    fig, ax = plt.subplots(figsize=figSize, dpi=DPI)
    
    mTrain, cTrain = np.polyfit(yTrain, yPredTrain, deg=1)
    perfectFitTrain = np.linspace(yTrain.min(), yTrain.max(), num=100)
    plt.scatter(yTrain, yPredTrain, color=trainDataColour, s=markerSize, alpha=markerAlpha, lw=markerLW, edgecolors=PLOT_COLOURS[0], zorder=2)
    plt.plot(perfectFitTrain, cTrain + mTrain*perfectFitTrain, color=trainBestFitColour, linestyle='--', lw=bestFitLW, zorder=3)

    mTest, cTest = np.polyfit(yTest, yPredTest, deg=1)
    perfectFitTest = np.linspace(yTest.min(), yTest.max(), num=100)
    plt.scatter(yTest, yPredTest, color=testDataColour, s=markerSize, alpha=markerAlpha, lw=markerLW, edgecolors=PLOT_COLOURS[0], zorder=2)
    plt.plot(perfectFitTest, cTest + mTest*perfectFitTest, color=testBestFitColour, linestyle='--', lw=bestFitLW, zorder=3)

    perfectFit = np.linspace(min(yTrain.min(), yTest.min()), max(yTrain.max(), yTest.max()), num=100)
    plt.plot(perfectFit, perfectFit, color=PLOT_COLOURS[0], lw=perfectFitLW, zorder=1)
    
    plt.xlabel(r'Actual $D_{B}$')
    plt.ylabel(r'Predicted $D_{B}$')
    plt.grid(linestyle='dotted')
    plt.legend(['Testing Data', 'Testing Best Fit', 'Training Data', 'Training Best Fit', 'Perfect Fit'])
    plt.savefig(figName, bbox_inches='tight')
    if showFig:
        plt.show()


def hyperparamTune(estimator, paramGrid, Xtrain, yTrain, Xtest, yTest, scoringMetric='r2', cv=KFold(n_splits=5, shuffle=True, random_state=RANDOM_SEED), hideWarnings=True):
    gridSearch = GridSearchCV(estimator=estimator, param_grid=paramGrid, cv=cv, scoring=scoringMetric, 
                               n_jobs=1, error_score=np.nan, verbose=0)
    # randomSearch = RandomizedSearchCV(estimator=estimator, n_jobs=-1, cv=cvFold, param_distributions=paramGrid, scoring=scoringMetric)
    if hideWarnings:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            searchResults = gridSearch.fit(Xtrain, yTrain)
    else:
        searchResults = gridSearch.fit(Xtrain, yTrain)
    bestModel = searchResults.best_estimator_
    print(f"Best score: {bestModel.score(Xtest, yTest):.3f}")
    print(f"Best parameter combination: {searchResults.best_params_}")
    return searchResults


def plotLearningCurve(estimator, Xtrain, yTrain, modelName, figSize=(4.5, 3.5), metric='Coefficient of determination', figName='figures/DBoxLearnCurve.png', showFig=False):
    fig, ax = plt.subplots(figsize=figSize, dpi=DPI)
    learningCurveParams = {'X': Xtrain, 'y': yTrain, 'train_sizes': np.linspace(0.1, 1.0, 5), 
                           'cv': ShuffleSplit(n_splits=50, test_size=TEST_SET_FRACTION, random_state=RANDOM_SEED), 
                           'score_type': 'both', 'n_jobs': NUM_JOBS, 'std_display_style': 'fill_between', 'score_name': metric, 
                           'line_kw': {'marker': 'o'}}
    LearningCurveDisplay.from_estimator(estimator, **learningCurveParams, ax=ax)
    handles, label = ax.get_legend_handles_labels()
    plt.xlabel('Number of training samples')
    plt.legend(handles[:2], ["Training Score", "Test Score"])
    plt.grid(linestyle='dotted')
    plt.savefig(figName, bbox_inches='tight')
    if showFig:
        plt.show()


def plotTopImpFeats(model, Xtrain, modelName, isForest=False, numFeat=40, figSize=(13, 3), figName='figures/DBoxFeatImp.png', showFig=False):
    featureNames = model.feature_names_in_ if modelName != 'LGBMR' else model.feature_name_
    if isForest:
        impStds = np.std([tree.feature_importances_ for tree in model.estimators_])
        featImps = pd.DataFrame({'feature': featureNames, 'importance': model.feature_importances_, 'std': impStds})
    elif 'ETR' in modelName:
        featImps = pd.DataFrame({'feature': featureNames, 'importance': model.estimator.feature_importances_})
    else:
        featImps = pd.DataFrame({'feature': featureNames, 'importance': model.feature_importances_})
    featImps.sort_values(by='importance', ascending=False, inplace=True)
    topFeats = featImps.iloc[:numFeat, :]

    fig, ax = plt.subplots(figsize=figSize, dpi=DPI)
    if isForest:
        plt.bar(x=topFeats['feature'], height=topFeats['importance'], yerr=impStds, width=0.9, color=PLOT_COLOURS[-1], edgecolor=PLOT_COLOURS[0], lw=1)
    else:
        plt.bar(x=topFeats['feature'], height=topFeats['importance'], width=0.9, color=PLOT_COLOURS[-1], edgecolor=PLOT_COLOURS[0], lw=1)
    plt.xticks(rotation=45, ha='right')
    plt.ylabel('Feature importance')
    plt.grid(linestyle='dotted')
    plt.savefig(figName, bbox_inches='tight')
    if showFig:
        plt.show()


def runModel(Xtrain, yTrain, Xtest, yTest, model, modelName, datasetName, figNameAppendix, showFeatCoef=True, plotFeatImp=False, verbose=False, showFig=False):
    if verbose: 
        print(f"  Running {modelName}")
    yPredTrain, yPredTest, model = evalModel(model, Xtrain, yTrain, Xtest, yTest)
    plotR2(yTrain, yTest, yPredTrain, yPredTest, modelName, figName=f"figures/{datasetName}DBoxR2{modelName}_{figNameAppendix}.png", showFig=showFig)
    plotR2All(yTrain, yTest, yPredTrain, yPredTest, modelName, figName=f"figures/{datasetName}DBoxR2All{modelName}_{figNameAppendix}.png", showFig=showFig)
    plotLearningCurve(model, Xtrain, yTrain, modelName, figName=f"figures/{datasetName}DBoxLearnCurve{modelName}_{figNameAppendix}.png", showFig=showFig)
    if showFeatCoef:
        print("  Feature coefficients:")
        for (i, coefficient) in enumerate(model.coef_):
            if round(coefficient, 3) >= 0.001:
                print(f"    {model.feature_names_in_[i]}: {coefficient:.3f}")
    if plotFeatImp:
        plotTopImpFeats(model, Xtrain, modelName, figName=f"figures/{datasetName}DBoxFeatImp{modelName}_{figNameAppendix}.png", showFig=showFig)
    return yPredTrain, yPredTest, model