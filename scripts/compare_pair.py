import argparse

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

from scipy.stats import spearmanr
from sklearn.metrics import precision_recall_curve, average_precision_score, auc, precision_score, recall_score, matthews_corrcoef
from sklearn.metrics import PrecisionRecallDisplay

Q_VAL_THR = 0.01
LOR_THR = 0.5

parser = argparse.ArgumentParser(description='Plot a comparison between two Nanocompore results.')
parser.add_argument('--tsv1', help='Path to the first TSV', required=True)
parser.add_argument('--tsv2', help='Path to the second TSV', required=True)
parser.add_argument('--col1', help='Which column from the first TSV to use', required=True)
parser.add_argument('--col2', help='Which column from the second TSV to use', required=True)
parser.add_argument('--name1', help='Name for the first result.', required=True)
parser.add_argument('--name2', help='Name for the second result.', required=True)
# parser.add_argument('--mods', help='BED file with true modified positions.', required=True)
parser.add_argument('--ref', help='Sites with sufficient coverage in both conditions', required=True)
parser.add_argument('--shift', help='Shift the positions in the Nanocompore results before comparing them to the reference.', nargs='?', type=int, const=0, default=0)
args = parser.parse_args()


# ref_mods = pd.read_csv(args.mods, sep='\t', index_col=['chr', 'strand', 'genomicPos'])
# ref_mods['modified'] = True
# ref_mods = pd.read_csv(args.mods,
#                        sep='\t',
#                        usecols=['chr', 'strand', 'genomicPos'])
# ref_mods['modified'] = True
# ref_mods['bin'] = ref_mods['genomicPos'] // 9
# ref_mods_binned = ref_mods.groupby(['chr', 'strand', 'bin']).agg({'modified': 'any'})
# del ref_mods
# print(ref_mods_binned)
ref_raw = pd.read_csv(args.ref,
                      sep='\t',
                      usecols=['chr', 'strand', 'genomicPos', 'modified'])
ref_raw['bin'] = ref_raw['genomicPos'] // 9
ref = ref_raw.groupby(['chr', 'strand', 'bin']).agg({'modified': 'any'})
del ref_raw
print(ref)

# sites = set()
# with open(args.sites) as file:
#     for line in file.readlines():
#         _, _, gpos, chrom, strand = line.split('\t')
#         sites.add((chrom, strand, gpos))
# sites = pd.read_csv(args.sites, sep='\t', usecols=['chr', 'strand', 'genomicPos'])
# sites['bin'] = sites.genomicPos // 9
# del sites['genomicPos']
# sites = sites.set_index(['chr', 'strand', 'bin'])
# sites['modified'] = False
# sites.loc[ref_mods_binned.index, 'modified'] = True

# print(sites)


tsvs = [args.tsv1, args.tsv2]
cols = [args.col1, args.col2]
names = [args.name1, args.name2]


def get_binned(run):
    df = pd.read_csv(tsvs[run], sep='\t') #, nrows=100000)
    if args.shift != 0:
        print(f'Shifting results {tsvs[run]} by {args.shift}')
        df['genomicPos'] += np.where(df['strand'] == '+', args.shift, -args.shift)
    col = cols[run]
    # The q-value for positions we miss is set to 1.0
    df[col] = df[col].fillna(1.0)
    # If q-value is 0, -log10 will be infinity so we
    # add an epsilon.
    df[col] = np.where(df[col] == 0, np.finfo(float).tiny, df[col])

    if LOR_THR:
        df[col] = np.where(np.abs(df['GMM_LOR']) < LOR_THR, 1.0, df[col])

    df['predicted'] = -np.log10(df[col])
    df['bin'] = df['genomicPos'] // 9
    # binned = df.groupby(['chr', 'strand', 'bin']).agg({'predicted': 'max', 'modified': 'any', 'ref_id': 'count', 'gx_pos': list})
    # binned = df.groupby(['chr', 'strand', 'bin']).agg({'predicted': 'max'})
    def best_site(group):
        idx = group.predicted.argmax()
        row = group.iloc[idx]
        return row['predicted'], abs(row['GMM_LOR'])
    binned = pd.DataFrame(
            df.groupby(['chr', 'strand', 'bin']).apply(best_site, include_groups=False)
        ).apply(lambda row: row[0], axis=1, result_type='expand')
    binned.columns = ['predicted', 'LOR']

    binned = ref.join(binned,
                      on=['chr', 'strand', 'bin'],
                      how='left')
    # If the reference is missing a position, we assume
    # it's a true negative.
    # binned['modified'] = binned['modified'].fillna(False)
    # The -log10 q-value for positions we miss is set to 0
    binned['predicted'] = binned['predicted'].fillna(0)
    print(binned)

    return binned

# figure, axes = plt.subplots(3, 2, figsize=(12, 15))
fig = plt.figure(figsize=(14, 12))

gs = fig.add_gridspec(3, 3, height_ratios=[1, 5, 5])
ax0 = fig.add_subplot(gs[0, :])
ax10 = fig.add_subplot(gs[1, 0])
ax11 = fig.add_subplot(gs[1, 1])
ax12 = fig.add_subplot(gs[1, 2])
ax20 = fig.add_subplot(gs[2, 0])
ax21 = fig.add_subplot(gs[2, 1])
ax22 = fig.add_subplot(gs[2, 2])


table_data = []

# fig.patch.set_visible(False)
ax0.axis('off')
ax0.axis('tight')


if __name__ == '__main__':
    df1 = get_binned(0)
    df2 = get_binned(1)

    # PRC
    for i, df in enumerate([df1, df2]):
        ytrue = df['modified']
        ypred = df['predicted']

        precision, recall, thresholds = precision_recall_curve(ytrue, ypred)
        f1_scores = 2 * recall * precision / (recall + precision)
        best_thresh = thresholds[np.argmax(f1_scores)]
        # ap = average_precision_score(ytrue, ypred)
        area = auc(recall, precision)
        print(area)

        display = PrecisionRecallDisplay(
                precision=precision,
                recall=recall)
        best_qval = 10 ** -best_thresh
        display.plot(name=names[i] + f" (AUC={area:.2f}, best_q={best_qval:.2f})", ax=ax10)
        display.ax_.set_title("Precision-recall curve")

        display.ax_.legend()

        ytrue_ = df.modified.astype(int)
        ypred_ = np.where(df.predicted >= -np.log10(Q_VAL_THR), 1, 0)
        detected = f"{np.sum(ytrue_ & ypred_)}/{ytrue_.sum()}"
        fp = np.sum(~ytrue_ & ypred_)
        precision_ = precision_score(ytrue_, ypred_)
        recall_ = recall_score(ytrue_, ypred_)
        f1_ = 2 * recall_ * precision_ / (recall_ + precision_)
        mcc_ = matthews_corrcoef(ytrue_, ypred_)

        table_data.append((detected, fp, precision_, recall_, f1_, mcc_))

    table = pd.DataFrame(table_data, columns=['Detected', 'FP', 'Precision', 'Recall', 'F1', 'MCC']).round(3)
    ax0.table(cellText=table.values, colLabels=table.columns, rowLabels=names, loc='center')

    xmax = max(df1.LOR.max(), df2.LOR.max()) + 0.15
    ymax = max(df1.predicted.max(), df2.predicted.max()) + 15

    # sharkfin plots

    sns.scatterplot(data=df1,
                    x='LOR',
                    y='predicted',
                    hue='modified',
                    s=10,
                    alpha=0.5,
                    linewidth=0,
                    ax=ax11)
    ax11.set_title(f'{names[0]}')
    ax11.set_xlabel('|LOR|')
    ax11.set_ylabel('-log10(q-value)')
    ax11.set_xlim(0, xmax)
    ax11.set_ylim(0, ymax)
    ax11.set_aspect(xmax/ymax)
    sns.scatterplot(data=df2,
                    x='LOR',
                    y='predicted',
                    hue='modified',
                    s=10,
                    alpha=0.5,
                    linewidth=0,
                    ax=ax12)
    ax12.set_title(f'{names[1]}')
    ax12.set_xlabel('|LOR|')
    ax12.set_ylabel('-log10(q-value)')
    ax12.set_xlim(0, xmax)
    ax12.set_ylim(0, ymax)
    ax12.set_aspect(xmax/ymax)

    # correlation plots
    # all
    joined = df1.join(df2,
                      on=['chr', 'strand', 'bin'],
                      how='inner',
                      rsuffix='_2')
    spearman_corr = spearmanr(joined['predicted'],
                              joined['predicted_2']).statistic
    sns.scatterplot(data=joined,
                    x='predicted',
                    y='predicted_2',
                    hue='modified',
                    s=10,
                    alpha=0.5,
                    linewidth=0,
                    ax=ax20)
    ax20.set_xlabel(f'-log10(q-value) {names[0]}')
    ax20.set_ylabel(f'-log10(q-value) {names[1]}')
    ax20.set_title(f'q-value correlation (r={spearman_corr:.2f})')
    ax20.set_xlim(-5, ymax)
    ax20.set_ylim(-5, ymax)
    ax20.set_aspect('equal', adjustable='box')

    # modified
    spearman_corr_mod = spearmanr(joined[joined.modified]['predicted'],
                                  joined[joined.modified]['predicted_2']).statistic
    sns.scatterplot(data=joined[joined.modified],
                    x='predicted',
                    y='predicted_2',
                    s=10,
                    alpha=0.5,
                    linewidth=0,
                    color='tab:orange',
                    ax=ax21)
    ax21.set_xlabel(f'-log10(q-value) {names[0]}')
    ax21.set_ylabel(f'-log10(q-value) {names[1]}')
    ax21.set_title(f'Modified: q-value correlation (r={spearman_corr_mod:.2f})')
    ax21.set_xlim(-5, ymax)
    ax21.set_ylim(-5, ymax)
    ax21.set_aspect('equal', adjustable='box')

    # non-modified
    spearman_corr_nomod = spearmanr(joined[~joined.modified]['predicted'],
                                  joined[~joined.modified]['predicted_2']).statistic
    sns.scatterplot(data=joined[~joined.modified],
                    x='predicted',
                    y='predicted_2',
                    s=10,
                    alpha=0.5,
                    linewidth=0,
                    ax=ax22)
    ax22.set_xlabel(f'-log10(q-value) {names[0]}')
    ax22.set_ylabel(f'-log10(q-value) {names[1]}')
    ax22.set_title(f'Non-modified: q-value correlation (r={spearman_corr_nomod:.2f})')
    ax22.set_xlim(-5, ymax)
    ax22.set_ylim(-5, ymax)
    ax22.set_aspect('equal', adjustable='box')



plt.tight_layout()
fig.savefig('fig.png')

