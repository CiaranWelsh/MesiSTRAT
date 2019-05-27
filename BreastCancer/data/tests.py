import unittest
import os, glob
import xlrd
import matplotlib.pyplot as plt
import seaborn

from .data_analysis import *


class ParseDataTests(unittest.TestCase):

    def setUp(self):
        self.data_dir = os.path.join(WORKING_DIRECTORY, 'data')
        self.data_file = os.path.join(self.data_dir, 'experimental_data.xlsx')
        assert os.path.isfile(self.data_file)
        self.gd = GetData(self.data_file)

    def test_open_workbook(self):
        self.assertIsInstance(self.gd.workbook, xlrd.book.Book)

    def test_get_data_from_sheet(self):
        data = self.gd.get_raw_data()
        self.assertIsInstance(data, pandas.DataFrame)

    def test_raw_data_is_accurate(self):
        data = self.gd.get_raw_data()
        self.assertEqual(2599775, data.loc[('MCF7', 0, 0), 'Akt'])
        self.assertEqual(1156087, data.loc[('MCF7', 0, 0), 'AktpT308'])
        self.assertEqual(8351448, data.loc[('MCF7', 90, 0), 'PRAS40pT246'])
        self.assertEqual(5833.17, data.loc[('T47D', 120, 3), 'Coomassie staining'])

    def test_get_data_normed_to_average(self):
        top_left = 1.422331386
        bottom_left = 0.822318749
        top_right = 0.980851924
        bottom_right = 0.990634001
        df = self.gd.get_data_normed_to_average()
        self.assertAlmostEqual(df.loc[('MCF7', 0), ('Akt', 0)], top_left)
        self.assertAlmostEqual(df.loc[('T47D', 120), ('Akt', 0)], bottom_left)
        self.assertAlmostEqual(df.loc[('MCF7', 0), ('Coomassie staining', 3)], top_right)
        self.assertAlmostEqual(df.loc[('T47D', 120), ('Coomassie staining', 3)], bottom_right)

    def test_get_data_normed_for_coomassie_blue(self):
        top_left = 1.471104582
        bottom_left = 0.770400628
        top_right = 1.650671097
        bottom_right = 0.956924173
        df = self.gd.get_data_normalised_to_coomassie_blue()
        self.assertAlmostEqual(df.loc[('MCF7', 0), ('Akt', 0)], top_left)
        self.assertAlmostEqual(df.loc[('T47D', 120), ('Akt', 0)], bottom_left)
        self.assertAlmostEqual(df.loc[('MCF7', 0), ('GAPDH', 3)], top_right)
        self.assertAlmostEqual(df.loc[('T47D', 120), ('GAPDH', 3)], bottom_right)


class DataTests(unittest.TestCase):

    def setUp(self):
        self.data_dir = os.path.join(WORKING_DIRECTORY, 'data')
        self.data_file = os.path.join(self.data_dir, 'experimental_data.xlsx')
        assert os.path.isfile(self.data_file)
        self.gd = GetData(self.data_file)
        self.raw_data = self.gd.get_raw_data()
        self.normed_to_average_data = self.gd.get_data_normed_to_average()
        self.normed_to_coomassie = self.gd.get_data_normalised_to_coomassie_blue()
        self.plots_dir = os.path.join(self.data_dir, 'plots')
        if not os.path.isdir(self.plots_dir):
            os.makedirs(self.plots_dir)

    def test_plots(self):
        raw = self.gd.get_raw_data()
        av = self.gd.get_data_normed_to_average()
        coomassie = self.gd.get_data_normalised_to_coomassie_blue()
        plot(raw, 'raw', True)
        plot(av, 'average', True)
        plot(coomassie, 'coomassie', True)
        # plt.show()

    def test_pca(self):
        data = self.gd.get_data_normalised_to_coomassie_blue()
        principle_component_analysis(data, colourby='cell_line', savefig=True)
        principle_component_analysis(data, colourby='time', savefig=True)
        principle_component_analysis(data, colourby='repeats', savefig=True)

    def test_get_ic_data(self):
        data = ss_data_to_copasi_format()
        y = data.loc[('MCF7', 0), 'IRS1pS636_639_obs']
        self.assertAlmostEqual(0.8613332330041832, y)

    def test_to_copasi_format_isfile(self):
        fname = os.path.join(DATA_DIRECTORY, 'copasi_formatted_data.csv')
        self.gd.to_copasi_format(fname, '\t')
        self.assertTrue(os.path.isfile(fname))

    def test_to_copasi_format_num_columns(self):
        fname = os.path.join(DATA_DIRECTORY, 'copasi_formatted_data.csv')
        self.gd.to_copasi_format(fname, '\t')
        with open(fname, 'r') as f:
            data = f.read()
        data = data.split('\n')
        columns = data[0].split('\t')
        expected = 20
        actual = len(columns)
        self.assertEqual(expected, actual)



if __name__ == '__main__':
    unittest.main()
