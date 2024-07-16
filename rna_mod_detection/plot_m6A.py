import os, csv, math
import numpy as np
import pandas as pd
from collections import Counter
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc, precision_recall_curve, accuracy_score


def intersection(lst1, lst2):
    return list(set(lst1) & set(lst2))


def compute_roc_auc(y_true, y_pred):
    fpr, tpr, _ = roc_curve(y_true, y_pred)
    roc_auc = auc(fpr, tpr)
    return roc_auc


def compute_pr_auc(y_true, y_pred):
    precision, recall, _ = precision_recall_curve(y_true, y_pred, pos_label=1)
    pr_auc = auc(recall, precision)
    return pr_auc


def plot_roc_auc(model_df, gt):

    plt.figure(dpi=600)
    plt.rc('font', family='Times New Roman')
    to_compare = ["Tombo", "MINES", "Nanom6A", "m6Anet", "Epinano", "SegPore"]
    model_df = model_df[to_compare + [gt]].dropna()
    for col in to_compare:
        y_true, y_pred = model_df[gt].values, model_df[col].values
        fpr, tpr, _ = roc_curve(y_true, y_pred)
        roc_auc = round(auc(fpr, tpr), 3)
        if col == 'Tombo':
            col = 'Tombo_de_novo'
        plt.plot(fpr, tpr, label=col + "(AUC=" + str(roc_auc) + ")", linewidth=1.0)

    plt.title('Ground Truth (' + gt + ')', fontsize=17)
    plt.xlabel("False positive rate", fontsize=15)
    plt.ylabel("True positive rate", fontsize=15)
    plt.legend(fontsize=12)
    x = np.linspace(0, 1)
    plt.plot(x, x, ls="-.", linewidth=0.5, color="gray")
    plt.subplots_adjust(left=0.09, right=0.99, top=0.94, bottom=0.1)
    # plt.show()
    plt.savefig(gt+"_roc.pdf")


def plot_pr_auc(model_df, gt):
    plt.figure(dpi=600)
    plt.rc('font', family='Times New Roman')
    to_compare = ["Tombo", "MINES", "Nanom6A", "m6Anet", "Epinano", "SegPore"]
    model_df = model_df[to_compare + [gt]].dropna()
    for col in to_compare:
        y_true, y_pred = model_df[gt].values, model_df[col].values

        precision, recall, _ = precision_recall_curve(y_true, y_pred, pos_label=1)
        pr_auc = round(auc(recall, precision), 3)

        if col == 'Tombo':
            col = 'Tombo_de_novo'
        plt.plot(recall, precision, label=col + "(AUC=" + str(pr_auc) + ")", linewidth=1.0)

    plt.title('Ground Truth (' + gt + ')', fontsize=17)
    plt.xlabel("Recall", fontsize=15)
    plt.ylabel("Precision", fontsize=15)
    plt.legend(fontsize=12)
    plt.subplots_adjust(left=0.09, right=0.99, top=0.94, bottom=0.1)
    x = np.linspace(0, 1)
    # plt.show()
    plt.savefig(gt+"_pr.pdf")


if __name__ == '__main__':

    model_df = pd.read_csv("m6A.benchmark.csv")
    plot_roc_auc(model_df, "MeRIP-seq")
    plot_roc_auc(model_df, "miCLIP")
    plot_roc_auc(model_df, "miCLIP2")
    plot_pr_auc(model_df, "MeRIP-seq")
    plot_pr_auc(model_df, "miCLIP")
    plot_pr_auc(model_df, "miCLIP2")

