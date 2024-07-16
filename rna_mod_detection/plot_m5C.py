import os, csv, math
import numpy as np
import pandas as pd
from collections import Counter
import matplotlib.pyplot as plt

from sklearn.metrics import roc_curve, auc, precision_recall_curve, accuracy_score



def compute_roc_auc(y_true, y_pred):
    fpr, tpr, _ = roc_curve(y_true, y_pred)
    roc_auc = auc(fpr, tpr)
    return roc_auc


def compute_pr_auc(y_true, y_pred):
    precision, recall, _ = precision_recall_curve(y_true, y_pred, pos_label=1)
    pr_auc = auc(recall, precision)
    return pr_auc


def plot_roc_auc(model_df):
    plt.figure(dpi=600)
    plt.rc('font', family='Times New Roman')
    model_df = model_df[to_compare + [gt]]
    for col in to_compare:
        if col == 'CHEUI-solo':
            _df = model_df[to_compare + [gt]].dropna()
            y_true, y_pred = _df[gt].values, _df[col].values
        else:
            y_true, y_pred = model_df[gt].values, model_df[col].values

        fpr, tpr, _ = roc_curve(y_true, y_pred)
        roc_auc = round(auc(fpr, tpr), round_n)
        plt.plot(fpr, tpr, label=col + "(AUC=" + str(roc_auc) + ")", linewidth=1.0)
    plt.xlabel("False positive rate")
    plt.ylabel("True positive rate")
    plt.legend()
    x = np.linspace(0, 1)
    plt.plot(x, x, ls="-.", linewidth=0.5, color="gray")
    plt.title("ROC AUC")
    # plt.show()
    plt.savefig("m5C_roc.pdf")


def plot_pr_auc(model_df):
    plt.figure(dpi=600)
    plt.rc('font', family='Times New Roman')
    model_df = model_df[to_compare + [gt]]
    for col in to_compare:
        if col == 'CHEUI-solo':
            _df = model_df[to_compare + [gt]].dropna()
            y_true, y_pred = _df[gt].values, _df[col].values
        else:
            y_true, y_pred = model_df[gt].values, model_df[col].values

        precision, recall, _ = precision_recall_curve(y_true, y_pred, pos_label=1)
        pr_auc = round(auc(recall, precision), round_n)
        if col == 'SegPore':
            plt.plot(recall, precision, label=col + "(AUC=" + str(pr_auc) + ")", linewidth=2.0, color="black")
        else:
            plt.plot(recall, precision, label=col + "(AUC=" + str(pr_auc) + ")", linewidth=1.0)

    plt.xlabel("Recall")
    plt.ylabel("Precision")
    plt.legend()
    x = np.linspace(0, 1)
    plt.title("PR AUC")
    # plt.show()
    plt.savefig("m5C_pr.pdf")

# gt = 'gt'
# to_compare = ['tombo']
#
# pos_model_df = pd.read_csv("IVT_m5C_test.fraction_modified_reads_only_data.plus.wig")
# pos_model_df['gt'] = 1
# pos_model_df['tombo'] = pos_model_df['data'].apply(lambda x: float(x.split(" ")[1]))
#
# neg_model_df = pd.read_csv("IVT_normalC_test.fraction_modified_reads._only_data.plus.wig")
# neg_model_df['gt'] = 0
# neg_model_df['tombo'] = neg_model_df['data'].apply(lambda x: float(x.split(" ")[1]))
#
# model_df = pd.concat([neg_model_df, pos_model_df])
# print(model_df)
# plot_roc_auc(model_df)
# plot_pr_auc(model_df)


def generate_dataset(file_name, save_file):

    alter_data = []
    chr = ""
    with open(file_name, 'r') as f:
        for line in f:
            line = line.strip().replace("\n", "").split()
            if line[0] == "track": continue
            if line[0] == "variableStep":
                chr = line[1].split("=")[1]
            else:
                line[0] = chr + ":" + line[0]
                alter_data.append(line)

    df = pd.DataFrame(alter_data, columns=['pos', 'pred'])
    df.to_csv(save_file, index=False)


def summary_data():
    model_df = pd.read_csv("tombo.csv")
    # gt = 'gt'
    # round_n = 3
    # to_compare = ['tombo_alternative', 'tombo_de_novo', 'CHEUI-solo']

    chuei = pd.read_csv("site_level_m5C_predictions.txt", sep='\t')
    chuei['position'] = chuei['position'] + 6
    chuei['position'] = chuei['position'].apply(lambda x: str(x))
    chuei['key'] = chuei['contig'] + ":" + chuei['position']
    chuei_dict = chuei.set_index("key")['probability'].to_dict()

    def get_cheui(x):
        if x not in chuei_dict:
            return None
        else:
            return chuei_dict[x]

    model_df['CHEUI-solo'] = model_df['pos'].apply(lambda x: get_cheui(x))
    model_df.to_csv("m5C.benchmark.csv", index=False)


if __name__ == '__main__':

    model_df = pd.read_csv("m5C.benchmark.csv")
    gt = 'gt'
    round_n = 3
    to_compare = ['Tombo_alternative', 'Tombo_de_novo', 'CHEUI-solo']

    plot_roc_auc(model_df)
    plot_pr_auc(model_df)
