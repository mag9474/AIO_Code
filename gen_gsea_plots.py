import fnmatch
import matplotlib.pyplot as plt
import gseapy
import seaborn as sns
import pandas as pd
import os


def return_title(working_dir):
    for name in os.listdir(working_dir):
        if fnmatch.fnmatch(name, "*.rpt"):
            file_dir = os.path.join(working_dir, name)
    df = pd.read_csv(file_dir)
    title = df[
        df.apply(lambda row: row.astype(str).str.contains("gene_sets").any(), axis=1)
    ]
    title = title.to_string()
    ind = title.find("gene_sets")
    title = title[ind + 10 :]

    return title


def gsea_plots(
    working_dir=r"C:\Users\sorth\OneDrive\NYU\Semester_3\AMI\GSEA_GUI\output",
    output_dir=r"C:\Users\sorth\OneDrive\NYU\Semester_3\AMI\GSEA_GUI\py_out",
):
    for dirs in os.listdir(working_dir):

        sub_dir = os.path.join(working_dir, dirs)
        save_name = dirs + ".png"
        save_dir = os.path.join(output_dir, save_name)
        for name in os.listdir(sub_dir):
            if fnmatch.fnmatch(name, "gsea_report_for_ClassA_*.tsv"):
                file_dir_a = os.path.join(sub_dir, name)
            elif fnmatch.fnmatch(name, "gsea_report_for_Classb_*.tsv"):
                file_dir_b = os.path.join(sub_dir, name)
        title = return_title(sub_dir)
        print("working_on:", title)
        if not "all" in title:
            print('not "all" geneset, skipping')
            continue

        working_dir_a = file_dir_a
        working_dir_b = file_dir_b

        df_gsea_a = pd.read_csv(working_dir_a, sep="\t")

        df_gsea_b = pd.read_csv(working_dir_b, sep="\t")

        df_gsea = df_gsea_a
        df_gsea = df_gsea.append(df_gsea_b)

        df_gsea = df_gsea.rename(columns={"FDR q-val": "FDR_q_val"})

        df_gsea_sort = df_gsea.sort_values(by=["NES"], ascending=False)
        df_gsea_sub = df_gsea_sort.loc[df_gsea_sort["FDR_q_val"] <= 0.25]

        if len(df_gsea_sub) > 20:
            df_gsea_sub = df_gsea_sub.head(10)
            df_gsea_sub = df_gsea_sub.append(df_gsea_sort.tail(10), ignore_index=True)
        elif len(df_gsea_sub) == 0:
            print(title, "had no significant pathways")
            continue

        if "df_gsea_sub_master" in locals():
            print("Adding to GSEA Master")
            df_to_add = df_gsea_sub.copy()
            df_to_add["NAME"] = title[0:2] + "_" + df_to_add["NAME"].astype(str)
            df_to_add["DIR"] = sub_dir
            df_gsea_sub_master = df_gsea_sub_master.append(
                df_to_add[["NAME", "FDR_q_val", "NES", "DIR"]], ignore_index=True
            )
        else:
            df_to_add = df_gsea_sub.copy()
            df_to_add["NAME"] = title[0:2] + df_to_add["NAME"].astype(str)
            df_to_add["DIR"] = sub_dir
            df_gsea_sub_master = df_to_add[["NAME", "FDR_q_val", "NES", "DIR"]]

        plot = plt.scatter(
            df_gsea_sub.NES, df_gsea_sub.NAME, c=df_gsea_sub.FDR_q_val, cmap="winter"
        )
        plt.clf()
        plt.figure(figsize=(10, 6))
        plt.colorbar(plot)

        ax = sns.barplot(
            x=df_gsea_sub.NES,
            y=df_gsea_sub.NAME,
            hue=df_gsea_sub.FDR_q_val,
            palette="winter",
            dodge=False,
        )
        ax.set_title(title)
        ax.set_ylabel("NAME")
        ax.set_xlabel("NES")
        ax.legend_.remove()
        ax2 = ax.twinx()
        ax2.set_ylabel("FDR_q_val")
        ax2.yaxis.set_label_coords(1.2, 0.5)
        ax2.set_yticklabels([])

        plt.savefig(save_dir, bbox_inches="tight")
        print(title, "Generated Plot")

    df_gsea_sub_master = df_gsea_sub_master.sort_values(by=["NES"], ascending=False)
    df_gsea_sub_master = df_gsea_sub_master.loc[df_gsea_sub_master["FDR_q_val"] <= 0.05]

    if len(df_gsea_sub_master) > 20:
        df_gsea_sub_master_sub = df_gsea_sub_master.head(10)
        df_gsea_sub_master_sub = df_gsea_sub_master_sub.append(
            df_gsea_sub_master.tail(10), ignore_index=True
        )
    df_gsea_sub_master = df_gsea_sub_master_sub
    plot = plt.scatter(
        df_gsea_sub_master.NES,
        df_gsea_sub_master.NAME,
        c=df_gsea_sub_master.FDR_q_val,
        cmap="winter",
    )
    plt.clf()
    plt.figure(figsize=(10, 6))
    plt.colorbar(plot)

    ax = sns.barplot(
        x=df_gsea_sub_master.NES,
        y=df_gsea_sub_master.NAME,
        hue=df_gsea_sub_master.FDR_q_val,
        palette="winter",
        dodge=False,
    )
    ax.set_title("GSEA Analysis")
    ax.set_ylabel("NAME")
    ax.set_xlabel("NES")
    ax.legend_.remove()
    ax2 = ax.twinx()
    ax2.set_ylabel("FDR_q_val")
    ax2.yaxis.set_label_coords(1.2, 0.5)
    ax2.set_yticklabels([])

    save_name = "GSEA_all.png"
    save_dir = os.path.join(output_dir, save_name)

    plt.savefig(save_dir, bbox_inches="tight")
    return df_gsea_sub_master


def get_top_genes(
    df, out_dir=r"C:\Users\sorth\OneDrive\NYU\Semester_3\AMI\GSEA_GUI\py_out"
):
    df["NAME"] = df["NAME"].str.slice(3)
    top_genes = []
    for i in range(len(df)):
        gene_name = str(df.iloc[i, 0]) + ".tsv"
        directory = df.iloc[i, -1]
        file_name = os.path.join(directory, gene_name)
        try:
            gene_df = pd.read_csv(file_name, sep="\t")
        except:
            gene_name = "C" + str(df.iloc[i, 0]) + ".tsv"
            file_name = os.path.join(directory, gene_name)
            gene_df = pd.read_csv(file_name, sep="\t")

        top_genes.extend(gene_df["SYMBOL"].values.tolist())

    top_genes_df = pd.DataFrame(data=top_genes, columns=["SYMBOL"])
    save_dir = os.path.join(out_dir, "top_genes_df.csv")
    top_genes_df.to_csv(save_dir, index=False)
    return top_genes_df


df_gsea_sub_master = gsea_plots()

df_gsea_sub_masterbk = df_gsea_sub_master.copy()
get_top_genes(df_gsea_sub_masterbk)
