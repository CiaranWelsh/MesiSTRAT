import glob, pandas, os
import GEOparse as gp
from functools import reduce



WORKING_DIR = r'C:\Users\Ciaran\Box Sync\MesiSTRAT\PublicDataSetSearch'
WORKING_DIR = r'C:\Users\cwels\Box Sync\MesiSTRAT\PublicDataSetSearch'
FILENAME = r'BreastCancerCrossRefWithMetabolism.txt'

# filename = os.path.join(WORKING_DIR, FILENAME)
#
#
# with open(filename, 'r') as f:
#     text = f.read()
#
# text = text.split('\n\n')
# # print(text, len(text))
#
# def reformat1(text):
#     text = text.split('\n')
#     text = [i for i in text if i != '']
#     titles = ['Name', 'Description', 'Organism', 'Type', 'Platform', 'FTP download', 'DataSet']
#     try:
#         data = {titles[i]: text[i] for i in range(len(text))}
#     except IndexError:
#         titles += 'fake'
#         data = {titles[i]: text[i] for i in range(len(text))}
#
#         # print('Failed with IndexError. {}'.format(text))
#         # return 0
#     return pandas.DataFrame(data, index=[0])
#
# def reformat(text):
#     df_list = []
#     for i in text:
#         if i != 0:
#             df_list.append(reformat1(i))
#     df = pandas.concat(df_list )
#     df['DataSet'] = df['DataSet'].str.replace('DataSet\t\tAccession: ', '')
#     df['DataSet'] = df['DataSet'].str.replace('\tID: ', '')
#
#     print(df.head())
#     return df

# reformat1(text[0])
# data = reformat(text)
# print(data)





geo_id_list = [
    "GDS1925",
    "GDS2626",
    "GDS2739",
    "GDS2861",
    "GDS3105",
    "GDS3139",
    "GDS3367",
    "GDS3484",
    "GDS3716",
    "GDS3910",
    "GDS4052",
    "GDS4053",
    "GDS4059",
    "GDS4063",
    "GDS4066",
    "GDS4069",
    "GDS4074",
    "GDS4079",
    "GDS4084",
    "GDS4085",
    "GDS4092",
    "GDS4096",
    "GDS4760",
    "GDS4814",
    "GDS5026",
    "GDS5306",
    "GDS5371",
    "GDS5800",
]


def download_GEOdata(geolist):
    dct = {}
    for i in geolist:
        destdir = os.path.join(WORKING_DIR, i)
        os.makedirs(destdir) if os.path.isdir(destdir) is not True else None
        try:
            fname = os.path.join(destdir, "{}_family.soft.gz".format(i))
            if os.path.isfile(fname):
                gse = gp.get_GEO(filepath=fname)
            else:
                gse = gp.get_GEO(geo=i, destdir=destdir)
        except FileNotFoundError:
            continue
        dct[i] = gse.metadata
        with open(os.path.join(destdir, 'meta_data.txt'), 'w') as f:
            for j in gse.metadata:
                res = gse.metadata[j][0] if len(gse.metadata[j]) == 1 else reduce(lambda x, y: "{};\t {}".format(x, y), gse.metadata[j])
                f.write('[{}]\n{}\n\n'.format(j, res))

    return dct

l = ['title', 'geo_accession', 'status', 'submission_date', 'last_update_date', 'pubmed_id',
     'summary', 'overall_design', 'type', 'contributor', 'sample_id', 'contact_name', 'contact_email', 'contact_phone', 'contact_laboratory', 'contact_department', 'contact_institute', 'contact_address', 'contact_city', 'contact_state', 'contact_zip/postal_code', 'contact_country', 'platform_id', 'platform_taxid', 'sample_taxid', 'relation']

# meta  = download_GEOdata(geo_id_list)
#
# print(meta['GSE5800'].keys())
for i in sorted(geo_id_list):
    print(i)









