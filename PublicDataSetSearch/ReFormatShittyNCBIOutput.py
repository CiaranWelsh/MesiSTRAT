import glob, pandas, os
import GEOparse as gp
from functools import reduce
from textwrap import wrap
import re
from Bio import Entrez
import webbrowser
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

WORKING_DIR = r'D:\MesiSTRAT\PublicDataSetSearch'

geo_id_list = [
    "GDS1250",
    "GDS1326",
    "GDS1329",
    "GDS1477",
    "GDS1510",
    "GDS1549",
    "GDS1627",
    "GDS1664",
    "GDS1873",
    "GDS1925",
    "GDS2045",
    "GDS2046",
    "GDS2095",
    "GDS2096",
    "GDS2097",
    "GDS2189",
    "GDS2250",
    "GDS2323",
    "GDS2324",
    "GDS2367",
    "GDS2415",
    "GDS2617",
    "GDS2618",
    "GDS2626",
    "GDS2635",
    "GDS2724",
    "GDS2739",
    "GDS2744",
    "GDS2745",
    "GDS2758",
    "GDS2759",
    "GDS2760",
    "GDS2761",
    "GDS2770",
    "GDS2787",
    "GDS2788",
    "GDS2789",
    "GDS2810",
    "GDS2827",
    "GDS2861",
    "GDS2864",
    "GDS3096",
    "GDS3097",
    "GDS3105",
    "GDS3116",
    "GDS3138",
    "GDS3139",
    "GDS3179",
    "GDS3217",
    "GDS3268",
    "GDS3283",
    "GDS3285",
    "GDS3315",
    "GDS3324",
    "GDS3367",
    "GDS3388",
    "GDS3484",
    "GDS360",
    "GDS3604",
    "GDS3609",
    "GDS3638",
    "GDS3716",
    "GDS3721",
    "GDS3722",
    "GDS3853",
    "GDS3910",
    "GDS3952",
    "GDS4051",
    "GDS4052",
    "GDS4053",
    "GDS4054",
    "GDS4055",
    "GDS4056",
    "GDS4057",
    "GDS4059",
    "GDS4061",
    "GDS4063",
    "GDS4065",
    "GDS4066",
    "GDS4067",
    "GDS4069",
    "GDS4070",
    "GDS4071",
    "GDS4074",
    "GDS4078",
    "GDS4079",
    "GDS4080",
    "GDS4082",
    "GDS4083",
    "GDS4084",
    "GDS4085",
    "GDS4088",
    "GDS4089",
    "GDS4090",
    "GDS4091",
    "GDS4092",
    "GDS4093",
    "GDS4095",
    "GDS4096",
    "GDS4110",
    "GDS4114",
    "GDS4296",
    "GDS4359",
    "GDS4388",
    "GDS4760",
    "GDS4761",
    "GDS4762",
    "GDS4763",
    "GDS4764",
    "GDS4765",
    "GDS4766",
    "GDS4797",
    "GDS4800",
    "GDS4814",
    "GDS4818",
    "GDS483",
    "GDS4885",
    "GDS4936",
    "GDS4972",
    "GDS4979",
    "GDS5026",
    "GDS5027",
    "GDS5076",
    "GDS5221",
    "GDS5222",
    "GDS5287",
    "GDS5306",
    "GDS5310",
    "GDS5348",
    "GDS5355",
    "GDS5371",
    "GDS5619",
    "GDS5620",
    "GDS5621",
    "GDS5662",
    "GDS5800",
    "GDS5801",
    "GDS5819",
    "GDS6100",
    "GDS806",
    "GDS807",
    "GDS817",
    "GDS820",
    "GDS823",
    "GDS825",
    "GDS845",
    "GDS846",
    "GDS847",
    "GDS848",
    "GDS88",
    "GDS881",
    "GDS892",
    "GDS893",
    "GDS901",
    "GDS992",
]

all_folders = [os.path.join(WORKING_DIR, i) for i in geo_id_list]


def download_GEOdata(geolist):
    dct = {}
    for i in geolist:
        destdir = os.path.join(WORKING_DIR, i)
        os.makedirs(destdir) if os.path.isdir(destdir) is not True else None
        try:
            fname = os.path.join(destdir, "{}_family.soft.gz".format(i))
            if os.path.isfile(fname):
                gse = gp.get_GEO(filepath=fname, annotate_gpl=True)
            else:
                gse = gp.get_GEO(geo=i, destdir=destdir, annotate_gpl=True)
        except FileNotFoundError:
            continue
        dct[i] = gse
        with open(os.path.join(destdir, 'meta_data.txt'), 'w') as f:
            for j in gse.metadata:
                res = gse.metadata[j][0] if len(gse.metadata[j]) == 1 else reduce(lambda x, y: "{};\t {}".format(x, y),
                                                                                  gse.metadata[j])
                res = '\n'.join(wrap(res, width=100))
                f.write('[{}]\n{}\n\n'.format(j, res))

        gse.table.set_index('ID_REF', inplace=True)
        data_fname = os.path.join(destdir, 'data_file.csv')
        column_fname = os.path.join(destdir, 'column_descriptions.csv')
        gse.table.to_csv(data_fname)
        gse.columns.to_csv(column_fname)

    return dct


def get_unique_platforms(folders):
    platform_list = []
    for i in folders:
        fname = os.path.join(i, 'meta_data.txt')
        l = ''
        with open(fname, 'r') as f:
            data = f.read()

        platform = re.findall('\[platform\]\n(.*)', data)
        platform_list.append(platform[0]) if platform[0] != [] else None
    platform_list = list(set(platform_list))
    print("Number of unique platforms is: {}".format(len(platform_list)))

    return platform_list


#
# for i in sorted(geo_id_list):
#     print("\"{}\",".format(i))

platforms = get_unique_platforms(all_folders)


def google_platforms():
    for i in platforms:
        url = r'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={}'.format(i)
        webbrowser.open(url, new=False)


# geo_dct = download_GEOdata(geo_id_list)
def plot_boxplots(geo_list):
    # current_dataset = geo_id_list[4]
    for i in geo_list:
        print('Plotting "{}"'.format(i))
        current_dir = os.path.join(WORKING_DIR, i)
        current_datafile = os.path.join(current_dir, 'data_file.csv')
        current_column_desc = os.path.join(current_dir, 'column_descriptions.csv')

        data = pandas.read_csv(current_datafile, index_col=[0, 1])
        desc = pandas.read_csv(current_column_desc)

        # print(rma(data))
        sns.set_style('white')
        sns.set_context('talk', font_scale=1)

        data = np.log2(data)
        data = pandas.DataFrame(data.stack())
        data.columns = ['Value']
        data = data.reset_index()

        # print(data)


        fig = plt.figure()
        sns.boxplot(data=data, x='level_2', y='Value')
        plt.xlabel('')
        plt.xticks(rotation=90)
        plt.title(i)
        plt.ylabel('Normalised log2 intensity')
        sns.despine(fig, top=True, right=True)
        fname = os.path.join(current_dir, 'boxplot.png')
        fig.savefig(fname, bbox_inches='tight', dpi=150)


plot_boxplots(geo_id_list)




































