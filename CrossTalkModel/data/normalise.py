
import pandas, os, glob
import xlrd
import numpy
import matplotlib.pyplot as plt
import seaborn


class GetData:
    new_names = {
        'D': 'DMSO',
        'T': 'TGFB',
        'E': 'TGFB +Everolimus',
        'EA72': 'TGFB + Everolimus+ AZD 3',
        'EA48': 'TGFB + Everolimus+ AZD 2',
        'EA24': 'TGFB + Everolimus+ AZD 1',
        'EA1.25': 'TGFB + Everolimus+ AZD 30 mins',
        'A72': 'TGFB + AZD 3',
        'A48': 'TGFB + AZD  2',
        'A24': 'TGFB + AZD 1',
        'A1.25': 'TGFB + AZD 30mins',
        'EM72': 'TGFB + Everolimus+ MK 3',
        'EM48': 'TGFB + Everolimus+ MK2',
        'EM24': 'TGFB + Everolimus+ MK 1',
        'EM1.25': 'TGFB + Everolimus+ MK30 mins',
        'M72': 'TGFB + MK 3',
        'M48': 'TGFB + MK 2',
        'M24': 'TGFB + MK 1',
        'M1.25': 'TGFB + MK 30mins',
    }
    new_names = {v: k for k, v in new_names.items()}

    col_indices = [1, 3, 5, 7]

    def __init__(self, workbook):
        self.workbook = xlrd.open_workbook(workbook)

    def get_sheet_names(self):
        return self.workbook.sheet_names()

    def get_data_from_sheet(self, sheet):
        sheet = self.workbook.sheet_by_name(sheet)
        names = sheet.col_slice(0, 44, 55)
        print('names', names)

        data = [sheet.col_slice(i, 44, 55) for i in self.col_indices]

        df = pandas.DataFrame(data).transpose()
        df['condition'] = names
        for i in range(df.shape[0]):
            for j in range(df.shape[1]):
                if df.iloc[i, j].ctype == 5:
                    df.iloc[i, j] = numpy.nan
                else:
                    df.iloc[i, j] = df.iloc[i, j].value

        condition_code = []
        for i in df['condition']:
            if i not in self.new_names.keys():
                raise ValueError
            condition_code.append(self.new_names[i])

        df['condition_code'] = condition_code
        df = df.set_index('condition_code')
        df = df.drop('condition', axis=1)
        return df

    def get_data(self):
        df_dct = {i: self.get_data_from_sheet(i) for i in self.get_sheet_names()}
        return pandas.concat(df_dct).stack().unstack(level=1)


def do_stats(data):
    mean_dct = {}
    sd_dct = {}
    for label, df in data.groupby(level=0):
        mean_dct[label] = df.mean()
        sd_dct[label] = df.std()

    mean = pandas.concat(mean_dct, axis=1, sort=False).transpose()
    sd = pandas.concat(sd_dct, axis=1, sort=False).transpose()
    return {
        'mean': mean,
        'sd': sd
    }

def plot_barplot(data, title):
    # seaborn.set_context(context='paper')
    conds = ['D', 'T', 'E', 'A72', 'M72', 'EA72', 'EM72']
    stacked = pandas.DataFrame(data[conds].stack())
    stacked.columns = ['AU']
    fig = plt.figure(figsize=(10, 5))
    seaborn.barplot(data=stacked.reset_index(), hue='condition_code', y='AU', x='protein',
                    n_boot=10000, errcolor='black',
                    errwidth=1, edgecolor='k', capsize=0.04,
                    )
    seaborn.despine(fig=fig, top=True, right=True)
    plt.legend(loc=(1, 0.5))
    plt.xlabel('')
    fname = os.path.join(wd, f'{title}.png')
    plt.title(title)

    plt.savefig(fname, bbox_inches='tight', dpi=300)

    plt.show()


def combine_two_datasets(azd_file, mk_file):
    # azd_v2 = GetData(azd_data_file_v2).get_data()
    # mk_v2 = GetData(mk_data_file_v2).get_data()
    azd_v3 = GetData(azd_file).get_data()
    mk_v3 = GetData(mk_file).get_data()

    azd_v3 = azd_v3.replace('', numpy.nan)
    mk_v3 = mk_v3.replace('', numpy.nan)
    # azd_v3 = azd_v3 * 10
    # mk_v3 = mk_v3 * 10

    mk_dte = mk_v3[['D', 'T', 'E']]
    mk_dte.reset_index(level=1, inplace=True)
    mk_dte.loc[:, 'level_1'] = mk_dte.loc[:, 'level_1'] + 4
    mk_dte = mk_dte.set_index('level_1', append=True)

    # print(mk_dte)
    mk_v3 = mk_v3.drop(['D', 'T', 'E'], axis=1)
    mk_v3 = pandas.concat([mk_dte, mk_v3], axis=1)
    mk_v3.index.names = ['protein', 'repeat']
    azd_v3.index.names = ['protein', 'repeat']

    # mk_v3.iloc[[0, 1, 2, 3], [0, 1, 2]] = azd_v3.iloc[[0, 1, 2, 3], [0, 1, 2]]

    return mk_v3.combine_first(azd_v3)


if __name__ == '__main__':
    wd = os.path.dirname(__file__)
    data_file = os.path.join(wd, 'HardCopy')

    azd_data_file_av = os.path.join(data_file, 'AZD_calculations - v5_norm_to_average.xlsx')
    mk_data_file_av = os.path.join(data_file, 'MK2206_calculations - v5_norm_to_average.xlsx')

    print(azd_data_file_av)


    # data = combine_two_datasets(azd_data_file_av, mk_data_file_av)
    # print(data)
    azd_v3 = GetData(azd_data_file_av).get_data()
    print(azd_v3)






















