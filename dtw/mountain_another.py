import json
from threading import Lock

import numpy
import pandas
import pymysql
from dtw.dtw_similarity import measure
from gtfparse import read_gtf
from matplotlib import pyplot as plt

# FIXME
DB_HOST = "localhost"
DB_PORT = 3306
DB_USER = "root"
DB_PASS = ""
DB_NAME = "clickmutation"


def distance_to_similarity(x):
    return 1.0 - (x / (1 + x))


global gtf_data
gtf_data = None

global temp_x
temp_x = 0


class Mountain:
    _commit_lock = Lock()

    def __init__(self, base_path):
        self.base_path = base_path
        self.cancer_type = self.base_path.strip("/").split("/")[-1]
        self.chromosome = "1"
        self.method = "median"
        with open(base_path + "metadata.json", "r") as f:
            self.metadata = json.load(f)

        self.xs = []
        self.ys = []
        self.instance = {
            "tumor": {"median": [], "mean": []},
            "normal": {"median": [], "mean": []},
        }

        self.cyto_band_initial = pandas.read_csv(
            self.base_path + "../cytoBand.txt", sep="\t", low_memory=False
        )
        self.cyto_band = None
        self.center = 0

        global gtf_data
        if gtf_data is None:
            self.df = read_gtf(self.base_path + "../Homo_sapiens.GRCh38.76.gtf")
            self.df = self.df[self.df["strand"] == "+"]
            gtf_data = self.df
        else:
            self.df = gtf_data
        self.df_chr = None
        self.db = pymysql.connect(
            host=DB_HOST,
            port=DB_PORT,
            user=DB_USER,
            password=DB_PASS,
            db=DB_NAME,
            charset="utf8",
            autocommit=True,
        )

    def init(self, chromosome: str, method="mean"):
        self.df_chr = self.df[self.df["seqname"] == chromosome]
        self.chromosome = chromosome
        self.method = method
        self.xs = []
        self.ys = []
        self.instance = {
            "tumor": {"median": [], "mean": []},
            "normal": {"median": [], "mean": []},
        }
        self.cyto_band = self.cyto_band_initial[
            self.cyto_band_initial["chr"] == "chr" + chromosome
        ]
        self.center = self.cyto_band.loc[self.cyto_band["cyto"].str.contains("q")].iloc[
            0
        ]["start"]

    def insert_intergenic_region(self):
        inergenic_regions = self.exec(
            """
        SELECT a.gene_name          AS a_g,
               b.gene_name          AS b_g,
               a.stop + 1           AS start,
               b.start - 1          AS stop,
               b.start - a.stop - 1 AS length,
               a.cyto,
               a.arm
        FROM (SELECT *, @num := @num + 1 AS rn
              FROM (SELECT @num := 0) r, c_gene_grch38_with_intergenic_region
              WHERE chr = %s
              ORDER BY id) a,
             (SELECT *, @num2 := @num2 + 1 AS rn
              FROM (SELECT @num2 := 0) r2, c_gene_grch38_with_intergenic_region
              WHERE chr = %s
              ORDER BY id) b
        WHERE a.rn + 1 = b.rn
          AND b.start - a.stop > 1""",
            (self.chromosome, self.chromosome),
        )

        gene_mean = 34825.7223
        exec_list = []
        for info in inergenic_regions:
            if info["length"] > gene_mean:
                slice_num = round(info["length"] / gene_mean)
                slice_len = round(info["length"] / slice_num)
                start = info["start"]
                i = 0
                for i in range(slice_num - 1):
                    exec_list.append(
                        (
                            "{}_#{}_{}".format(info["a_g"], i, info["b_g"]),
                            self.chromosome,
                            start,
                            start + slice_len,
                            info["cyto"],
                            info["arm"],
                        )
                    )
                    start = start + slice_len + 1
                exec_list.append(
                    (
                        "{}_#{}_{}".format(info["a_g"], i + 1, info["b_g"]),
                        self.chromosome,
                        start,
                        info["stop"],
                        info["cyto"],
                        info["arm"],
                    )
                )
            else:
                exec_list.append(
                    (
                        "{}_#0_{}".format(info["a_g"], info["b_g"]),
                        self.chromosome,
                        info["start"],
                        info["stop"],
                        info["cyto"],
                        info["arm"],
                    )
                )
        self.exec(
            """INSERT INTO c_gene_grch38_with_intergenic_region (gene_name, chr, start, stop, cyto, arm) 
                        VALUES (%s, %s, %s, %s, %s, %s)""",
            exec_list,
        )

    def insert_data(self):
        temp = {"tumor": [], "normal": []}
        sequences = self.exec(
            "SELECT id, start, stop FROM c_gene_grch38_with_intergenic_region WHERE chr = %s",
            self.chromosome,
        )
        for seq in sequences:
            for sample_type in ["tumor", "normal"]:
                if seq["start"] < temp_x and seq["stop"] < temp_x:
                    median_list = []
                    mean_list = []
                elif seq["start"] < temp_x:
                    median_list = self.instance[sample_type]["median"][
                        0 : seq["stop"] - temp_x
                    ]
                    mean_list = self.instance[sample_type]["mean"][
                        0 : seq["stop"] - temp_x
                    ]
                else:
                    median_list = self.instance[sample_type]["median"][
                        seq["start"] - temp_x : seq["stop"] - temp_x
                    ]
                    mean_list = self.instance[sample_type]["mean"][
                        seq["start"] - temp_x : seq["stop"] - temp_x
                    ]
                median = numpy.mean(median_list)
                mean = numpy.mean(mean_list)
                if numpy.isnan(median):
                    median = -1
                if numpy.isnan(mean):
                    mean = -1
                temp[sample_type].append(
                    (
                        seq["id"],
                        self.cancer_type.lower(),
                        round(float(median), 4),
                        round(float(mean), 4),
                    )
                )
            if len(temp["tumor"]) >= 500 or len(temp["normal"]) >= 500:
                for sample_type in ["tumor", "normal"]:
                    self.exec(
                        """INSERT INTO c_{}_y_grch38_with_intergenic_region_moutain
                                    (gene_id, cancer_name, mid, mean) VALUES (%s, %s, %s, %s)""".format(
                            sample_type[0]
                        ),
                        temp[sample_type],
                    )
                temp = {"tumor": [], "normal": []}
        for sample_type in ["tumor", "normal"]:
            self.exec(
                """INSERT INTO c_{}_y_grch38_with_intergenic_region_moutain
                            (gene_id, cancer_name, mid, mean) VALUES (%s, %s, %s, %s)""".format(
                    sample_type[0]
                ),
                temp[sample_type],
            )

    def data_to_mysql(self):
        start = end = 0
        end_temp = 0
        for index, sequence in self.df_chr.iterrows():
            if start <= sequence["start"] <= end and start <= sequence["end"] <= end:
                continue
            start = int(sequence["start"])
            end = int(sequence["end"])
            gene_name = sequence["gene_name"]

            c = c1 = ""
            arm = arm1 = ""
            for cyto_index, cyto in self.cyto_band.iterrows():
                if int(cyto["start"]) <= end_temp <= int(cyto["end"]):
                    c1 = cyto["cyto"]
                    arm1 = "p" if "p" in c1 else "q"
                if int(cyto["start"]) <= start <= int(cyto["end"]):
                    c = cyto["cyto"]
                    arm = "p" if "p" in c else "q"
                if c1 and arm1 and c and arm:
                    break

            self.exec(
                "INSERT INTO c_gene_grch38_with_intergenic_region (gene_name, chr, start, stop, cyto, arm) VALUES (%s, %s, %s, %s, %s, %s)",
                (
                    gene_name,
                    self.chromosome,
                    start,
                    end,
                    c,
                    arm,
                ),
            )

            end_temp = end
            gene_name_temp = gene_name

    def exec(self, sql: str, args=None):
        with self._commit_lock:
            cursor = self.db.cursor(pymysql.cursors.DictCursor)
            if type(args) == list:
                cursor.executemany(sql, args)
            else:
                cursor.execute(sql, args)
            data = cursor.fetchall()
        return data

    def draw(self, ax):
        data_dict = {"tumor": {}, "normal": {}}
        x_axis = {"tumor": [], "normal": []}

        for file in self.metadata:
            if self.cancer_type != "LAML":
                if file["associated_entities"][0]["entity_submitter_id"][
                    15
                ] != "A" and (
                    int(file["associated_entities"][0]["entity_submitter_id"][13:15])
                    == 10
                    or len(file["associated_entities"]) == 1
                ):
                    continue
                if (
                    int(file["associated_entities"][0]["entity_submitter_id"][13:15])
                    == 10
                ):
                    sample_type = "normal"
                elif int(
                    file["associated_entities"][0]["entity_submitter_id"][13:15]
                ) in range(1, 10):
                    sample_type = "tumor"
                else:
                    continue
            else:
                if file["associated_entities"][0]["entity_submitter_id"][
                    15
                ] != "A" and (
                    int(file["associated_entities"][0]["entity_submitter_id"][13:15])
                    == 11
                    or len(file["associated_entities"]) == 1
                ):
                    continue
                if (
                    int(file["associated_entities"][0]["entity_submitter_id"][13:15])
                    == 11
                ):
                    sample_type = "normal"
                elif int(
                    file["associated_entities"][0]["entity_submitter_id"][13:15]
                ) in range(1, 10):
                    sample_type = "tumor"
                else:
                    continue

            file_path = self.base_path + file["file_id"] + "/" + file["file_name"]
            try:
                data = pandas.read_csv(file_path, sep="\t", low_memory=False)
            except FileNotFoundError:
                continue
            except pandas.errors.EmptyDataError:
                print(file_path)
                exit()

            for i, row in data.iterrows():
                if row["Chromosome"] == self.chromosome:
                    x = "{},{},{}".format(row["GDC_Aliquot"], row["Start"], row["End"])
                    x_axis[sample_type].append(row["Start"])
                    x_axis[sample_type].append(row["End"])
                    y = row["Segment_Mean"]
                    data_dict[sample_type][x] = 2 ** y * 2

        for sample_type in ["tumor", "normal"]:
            color = "#2A8BC8" if sample_type == "tumor" else "#2AC376"
            # 先排序，然后从头开始取一个一个最小线段，取包含这些线段的序列的值，取平均/中位
            xx = list(numpy.sort(numpy.unique(x_axis[sample_type])).tolist())
            self.xs.extend(xx)

            x0 = xx.pop(0)
            global temp_x
            temp_x = x0
            while True:
                x1 = xx.pop(0)
                x = [x0, x1]
                values = []
                for key in data_dict[sample_type].keys():
                    k = key.split(",")  # GDC_Aliquot,Start,End
                    start = int(k[1])
                    end = int(k[2])
                    if start <= x0 <= end and start <= x1 <= end:
                        values.append(data_dict[sample_type][key])
                for show_type in ["mean", "median"]:

                    if not values:
                        print(self.cancer_type, self.chromosome, sample_type, x)

                    if show_type == "mean":
                        y = numpy.mean(values)
                    else:
                        y = numpy.median(values)
                    self.ys.append(y)

                    for i in range(x0, x1):
                        self.instance[sample_type][show_type].append(y)

                if sample_type == "tumor":
                    ax.plot(x, [y, y], color, label=sample_type, linewidth=4, alpha=0.5)
                if len(xx):
                    x0 = x1
                else:
                    break

    def main(self):
        fig, ax = plt.subplots()

        self.draw(ax)
        self.insert_data()
        self.data_abstract()

        plt.close()

    def data_abstract(self):
        data_normal = self.exec(
            """
        SELECT a.mid
        FROM c_n_y_grch38_with_intergenic_region_moutain as a,
             c_gene_grch38_with_intergenic_region as b
        WHERE a.cancer_name = %s
          AND b.chr = %s
          AND a.gene_id = b.id
          AND a.mid != -1
        ORDER BY b.start""",
            (
                self.cancer_type.lower(),
                self.chromosome,
            ),
        )
        data_normal = [a["mid"] for a in data_normal]
        data_tumor = self.exec(
            """
        SELECT a.mid
        FROM c_t_y_grch38_with_intergenic_region_moutain as a,
             c_gene_grch38_with_intergenic_region as b
        WHERE a.cancer_name = %s
          AND b.chr = %s
          AND a.gene_id = b.id
          AND a.mid != -1
        ORDER BY b.start""",
            (
                self.cancer_type.lower(),
                self.chromosome,
            ),
        )
        data_tumor = [a["mid"] for a in data_tumor]
        p_arm_length = self.exec(
            """
        SELECT COUNT(*) AS c
        FROM c_n_y_grch38_with_intergenic_region_moutain as a,
             c_gene_grch38_with_intergenic_region as b
        WHERE a.cancer_name = %s
          AND b.chr = %s
          AND a.gene_id = b.id
          AND a.mid != -1
          AND b.arm = 'p'
        ORDER BY b.start""",
            (
                self.cancer_type.lower(),
                self.chromosome,
            ),
        )[0]["c"]
        dtw_similarity_p = measure(
            data_tumor[:p_arm_length], data_normal[:p_arm_length]
        )
        dtw_similarity_q = measure(
            data_tumor[p_arm_length:], data_normal[p_arm_length:]
        )
        dtw_similarity = measure(data_tumor, data_normal)
        print(
            self.cancer_type,
            self.chromosome,
            dtw_similarity_p,
            dtw_similarity_q,
            dtw_similarity,
        )


if __name__ == "__main__":
    chromosomes = list(range(1, 23))
    chromosomes.extend(["X", "Y"])

    # for ct in ['KIRP', 'HNSC', 'THCA', 'DLBC', 'ACC', 'UCS', 'LUSC', 'PRAD', 'GBM', 'STAD', 'READ', 'PAAD', 'KICH',
    #            'KIRC', 'UVM', 'THYM', 'CECC']:
    # for ct in ['LIHC', 'PCPG', 'OV', 'BRCA', 'TGCT', 'CHOL', 'LUAD', 'MESO', 'UCEC', 'SARC', 'ESCC', 'CEAD', 'LGG', 'ESAD', 'SKCM', 'COAD', 'BLCA']:
    # for ct in ['CECC', 'CEAD', 'ESCC', 'ESAD']:
    #     mountain = Mountain('./mountain/{}/'.format(ct))
    #     for i in chromosomes:
    #         mountain.init(str(i), 'median')
    #         mountain.main()
    # for cancer_type in ['ACC', 'BRCA', 'CECC', 'CHOL', 'ESAD', 'ESCC', 'HNSC', 'KICH', 'KIRP', 'LGG', 'LUAD', 'MESO',
    #                     'PAAD', 'PRAD', 'SARC', 'STAD', 'THCA', 'UCEC', 'UVM', 'BLCA', 'CEAD', 'COAD', 'DLBC', 'GBM',
    #                     'KIRC', 'LIHC', 'LUSC', 'OV', 'PCPG', 'READ', 'SKCM', 'TGCT', 'THYM', 'UCS']:
    # ['PCPG', 'GBM', 'KIRP', 'LIHC', 'BLCA', 'SKCM', 'THCA', 'PAAD', 'BRCA']
    # ['STAD', 'LUSC', 'MESO', 'READ', 'DLBC', 'SARC', 'ACC', 'CHOL', 'TGCT']
    # ['PRAD', 'ESAD', 'THYM', 'UCS', 'LGG', 'UCEC', 'LUAD', 'COAD', 'KICH']
    # ['UVM', 'CECC', 'KIRC', 'ESCC', 'CEAD', 'OV', 'HNSC']
    for cancer_type in [
        "PCPG",
        "GBM",
        "KIRP",
        "LIHC",
        "BLCA",
        "SKCM",
        "THCA",
        "PAAD",
        "BRCA",
    ]:
        mountain = Mountain("./mountain/{}/".format(cancer_type))
        print(cancer_type)
        for i in chromosomes:
            mountain.init(str(i), "median")
            # mountain.data_to_mysql()
            # mountain.insert_intergenic_region()
            mountain.main()
