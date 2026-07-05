"""
Analyse de l'expression du gene MRPS16 dans deux cancers (GBM et LIHC)
compares a leurs tissus sains correspondants (GTEx Cerveau / Foie).

Source des donnees attendue : export CSV depuis UCSC Xena Browser
(xenabrowser.net -> cohorte "TCGA TARGET GTEx"), ou tout autre export
contenant au minimum les colonnes suivantes :

    sample_id, expression, group, cancer

- expression : valeur log2(TPM+1) du gene MRPS16 pour l'echantillon
- group      : "Tumeur" ou "Sain"
- cancer     : "GBM" ou "LIHC"

Adapter la fonction load_data() si le format exporte differe.
"""

import pandas as pd
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns


def load_data(gbm_path: str, lihc_path: str) -> pd.DataFrame:
    """
    Charge et nettoie les fichiers TSV exportes directement de UCSC Xena
    (cohorte TCGA TARGET GTEx). Colonnes source attendues :
    sample, samples, MRPS16, "primary disease or tissue", _primary_site

    - GBM : tumeur = "Glioblastoma Multiforme"
            sain   = toutes les sous-regions GTEx "Brain - ..."
            (on exclut "Brain Lower Grade Glioma", un autre cancer)
    - LIHC: tumeur = "Liver Hepatocellular Carcinoma"
            sain   = "Liver" (GTEx)
    """
    def read_xena_tsv(path):
        df = pd.read_csv(path, sep="\t")
        df = df.rename(columns={
            "MRPS16": "expression",
            "primary disease or tissue": "disease",
        })
        # retire l'eventuelle ligne d'en-tete dupliquee et les valeurs manquantes
        df = df[df["disease"] != "primary disease or tissue"]
        df = df.dropna(subset=["expression"])
        df["expression"] = pd.to_numeric(df["expression"], errors="coerce")
        df = df.dropna(subset=["expression"])
        return df

    gbm = read_xena_tsv(gbm_path)
    gbm = gbm[
        (gbm["disease"] == "Glioblastoma Multiforme")
        | (gbm["disease"].str.startswith("Brain -"))
    ].copy()
    gbm["group"] = np.where(gbm["disease"] == "Glioblastoma Multiforme", "Tumeur", "Sain")
    gbm["cancer"] = "GBM"

    lihc = read_xena_tsv(lihc_path)
    lihc = lihc[lihc["disease"].isin(["Liver Hepatocellular Carcinoma", "Liver"])].copy()
    lihc["group"] = np.where(lihc["disease"] == "Liver Hepatocellular Carcinoma", "Tumeur", "Sain")
    lihc["cancer"] = "LIHC"

    df = pd.concat([gbm, lihc], ignore_index=True)
    return df


def run_statistics(df: pd.DataFrame) -> pd.DataFrame:
    """
    Test de Mann-Whitney (Wilcoxon rank-sum) entre Tumeur et Sain,
    pour chaque cancer separement. Calcule aussi le log2 fold-change
    median comme mesure d'effet.
    """
    results = []
    for cancer in df["cancer"].unique():
        sub = df[df["cancer"] == cancer]
        tumeur = sub.loc[sub["group"] == "Tumeur", "expression"]
        sain = sub.loc[sub["group"] == "Sain", "expression"]

        stat, pval = stats.mannwhitneyu(tumeur, sain, alternative="two-sided")

        median_tumeur = tumeur.median()
        median_sain = sain.median()
        log2fc = median_tumeur - median_sain  # deja en log2(TPM+1)

        results.append({
            "cancer": cancer,
            "n_tumeur": len(tumeur),
            "n_sain": len(sain),
            "median_tumeur": round(median_tumeur, 3),
            "median_sain": round(median_sain, 3),
            "log2FC": round(log2fc, 3),
            "p_value": pval,
            "significatif": "oui" if pval < 0.05 else "non",
        })

    return pd.DataFrame(results)


def plot_boxplots(df: pd.DataFrame, stats_df: pd.DataFrame, output_path: str):
    """
    Genere un boxplot compare Tumeur vs Sain pour chaque cancer,
    avec p-value annotee.
    """
    sns.set_theme(style="whitegrid")
    fig, axes = plt.subplots(1, 2, figsize=(10, 5), sharey=False)

    for ax, cancer in zip(axes, df["cancer"].unique()):
        sub = df[df["cancer"] == cancer]
        sns.boxplot(
            data=sub, x="group", y="expression", hue="group", legend=False,
            order=["Sain", "Tumeur"], palette=["#4C72B0", "#DD8452"], ax=ax
        )
        sns.stripplot(
            data=sub, x="group", y="expression",
            order=["Sain", "Tumeur"], color="black", alpha=0.3, size=2, ax=ax
        )

        pval = stats_df.loc[stats_df["cancer"] == cancer, "p_value"].values[0]
        label = f"p = {pval:.2e}" if pval >= 1e-4 else "p < 1e-4"
        ax.set_title(f"MRPS16 — {cancer}\n{label}")
        ax.set_xlabel("")
        ax.set_ylabel("Expression log2(TPM+1)")

    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    print(f"Figure sauvegardee : {output_path}")


if __name__ == "__main__":
    # Adapter ces chemins vers vos fichiers exportes
    df = load_data("data/gbm_mrps16.tsv", "data/lihc_mrps16.tsv")

    stats_df = run_statistics(df)
    print(stats_df.to_string(index=False))
    stats_df.to_csv("results/mrps16_statistiques.csv", index=False)

    plot_boxplots(df, stats_df, "results/mrps16_boxplots.png")
