import pandas, numpy
from functools import reduce
import os, glob
import xlrd
import matplotlib.pyplot as plt
import seaborn

WORKING_DIRECTORY = r'D:\MesiSTRAT\BreastCancer'
DATA_DIRECTORY = os.path.join(WORKING_DIRECTORY, 'data')
DATA_FILE = os.path.join(DATA_DIRECTORY, 'experimental_data.xlsx')
PLOTS_DIR = os.path.join(DATA_DIRECTORY, 'plots')


class GetData:
    time = [0, 15, 30, 60, 90, 120]
    condition_names = ['MCF-7 0 minutes ins + aa', 'MCF-7 15 minutes ins + aa', 'MCF-7 30 minutes ins + aa',
                       'MCF-7 60 minutes ins + aa', 'MCF-7 90 minutes ins + aa', 'MCF-7 120 minutes ins + aa',
                       'T47D 0 minutes ins + aa', 'T47D 15 minutes ins + aa', 'T47D 30 minutes ins + aa',
                       'T47D 60 minutes ins + aa', 'T47D 90 minutes ins + aa', 'T47D 120 minutes ins + aa']
    condition_codes = [f'MCF{i}' for i in time] + [f'T47D{i}' for i in time]
    assert len(condition_names) == len(condition_codes)

    antibodies = ['Akt', 'AktpT308', 'AktpS473', 'PRAS40', 'PRAS40pT246', 'PRAS40pS183', 'S6K', 'S6KpT389', 'S6KpT229',
                  'TSC2', 'TSC2pT1462', 'IRS1', 'IRS1pS636/639', '4E-BP1', '4E-BP1pT37/46', 'GAPDH', 'ERK',
                  'Coomassie staining']
    antibodies = [i.replace('/', '_') for i in antibodies]
    antibodies = [i.replace('-', '_') for i in antibodies]

    def __init__(self, workbook):
        self.workbook = xlrd.open_workbook(workbook)

    def get_sheet_names(self):
        return self.workbook.sheet_names()

    def get_antibody_names(self):
        sheet = 'Sheet2'
        sheet = self.workbook.sheet_by_name(sheet)
        columns = sheet.row_slice(0)
        columns = [i.value for i in columns]
        columns = [i for i in columns if i != '']
        return columns

    def get_raw_data(self):
        sheet = 'Sheet2'
        sheet = self.workbook.sheet_by_name(sheet)

        data = [sheet.col_slice(i, 2, 14) for i in range(1, 73)]
        new_data = []
        for i in data:
            new_data.append([j.value for j in i])

        df = pandas.DataFrame(new_data).transpose()

        # sanity checks
        top_left = 2599775
        bottom_left = 1503056
        top_right = 5775.57
        bottom_right = 5833.17

        assert top_left == df.loc[0, 0]
        assert bottom_left == df.loc[df.shape[0] - 1, 0], f'{bottom_left} is not {df.loc[df.shape[0] - 1, 0]}'
        assert top_right == df.loc[0, df.shape[1] - 1], f'{top_right} != {df.loc[0, df.shape[1] - 1]}'
        assert bottom_right == df.loc[df.shape[0] - 1,
                                      df.shape[1] - 1], f'{top_right} != {df.loc[df.shape[0] - 1, df.shape[1] - 1]}'

        assert df.shape[1] / 4 == len(self.antibodies)

        df['cell_line'] = ['MCF7'] * 6 + ['T47D'] * 6
        df['time'] = [0, 15, 30, 60, 90, 120] * 2
        df = df.set_index(['cell_line', 'time'])

        df = df.transpose()
        antibody = [[i] * 4 for i in self.antibodies]
        antibody = reduce(lambda x, y: x + y, antibody)
        df['antibody'] = antibody

        repeats = [0, 1, 2, 3] * len(self.antibodies)
        df['repeats'] = repeats
        df = df.set_index(['antibody', 'repeats'])

        df = df.transpose()

        return df

    def get_data_normed_to_average(self):
        data = self.get_raw_data()
        df = data / data.mean()
        return df


    def get_data_normalised_to_coomassie_blue(self):
        # first normalise to average
        data = self.get_data_normed_to_average()
        df_dct = {}
        for ab in self.antibodies:
            df = data[ab]
            df_dct[ab] = df / data['Coomassie staining']
        return pandas.concat(df_dct, axis=1)


def plot(data, prefix, savefig=False):
    data = data.stack().stack()
    data = pandas.DataFrame(data)
    for label, df in data.groupby(level=[3]):
        fig = plt.figure()

        seaborn.lineplot(x='time', y=0, data=df.reset_index(),
                         hue='cell_line', style='cell_line',
                         palette='bright', markers=True, ci=95)
        seaborn.despine(fig, top=True, right=True)
        plt.ylabel('AU')
        plt.title(label)
        fname = os.path.join(PLOTS_DIR, '{}_{}'.format(prefix, label))
        if savefig:
            plt.savefig(fname, dpi=200, bbox_inches='tight')


def principle_component_analysis(data, colourby='cell_line', savefig=False):

    data = data.stack()
    if colourby not in data.index.names:
        raise ValueError(f'colourby variable "{colourby}" not in "{data.index.names}"')
    from sklearn.decomposition import PCA
    from itertools import cycle
    pca = PCA(2)
    x = pca.fit_transform(data)
    explained_variance = pca.explained_variance_ratio_

    df = pandas.DataFrame(x, index=data.index)
    colourby_values = sorted(list(set(df.index.get_level_values(level=colourby))))
    colours = seaborn.color_palette("hls", len(colourby_values))

    fig = plt.figure()
    for i in colourby_values:
        df2 = df.xs(i, level=colourby)
        plt.scatter(df2[0], df2[1], label=i)
    seaborn.despine(fig=fig, top=True, right=True)
    plt.legend(loc='best')
    plt.title('PCA coloured by {}'.format(colourby))
    plt.xlabel('PC1 {}%'.format(round(explained_variance[0], 2)*100))
    plt.ylabel('PC2 {}%'.format(round(explained_variance[1], 2)*100))
    plt_dir = os.path.join(PLOTS_DIR, 'PCAPlots')
    if not os.path.isdir(plt_dir):
        os.makedirs(plt_dir)
    fname = os.path.join(plt_dir, 'PCAColouredBy{}.png'.format(colourby))
    if savefig:
        plt.savefig(fname, dpi=200, bbox_inches='tight')




if __name__ == '__main__':
    gd = GetData(DATA_FILE)
