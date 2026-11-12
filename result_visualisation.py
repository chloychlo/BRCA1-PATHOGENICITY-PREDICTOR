import matplotlib.pyplot as plt
import pandas as pd

input_file = "brca1_vus_predictions_results_random_forest.csv"
output_image = "brca1_vus_predictions_visualisation.png"

df = pd.read_csv(input_file)

plt.style.use("seaborn-v0_8-whitegrid")
fig, axes = plt.subplots(1, 2, figsize=(16, 6))
fig.suptitle(
    "Distribution of predictions and spatial profile of variant pathogenicity of BRCA1 VUS (n = {})".format(len(df)),
    fontsize=16,
    fontweight="bold",
    y=1.02,
)

# ---------------------------------------------------------
# GRAPHIQUE 1 : Distribution des scores prédits
# ---------------------------------------------------------
ax1 = axes[0]
probs = df["predicted_pathogenic_probability"]

ax1.hist(probs, bins=30, color="#4A777A", edgecolor="white", alpha=0.85)
ax1.axvline(
    0.5, color="#C0392B", linestyle="--", linewidth=2, label="Seuil de décision (0.5)"
)

n_patho = (probs >= 0.5).sum()
n_benign = (probs < 0.5).sum()
ax1.set_title(
    f"Distribution of Predicted probabilities\n(Number of estimated pathogens : {n_patho} | Number of estimated benigns : {n_benign})",
    fontsize=12,
    fontweight="bold",
)
ax1.set_xlabel("Predicted pathogenicity probability", fontsize=11)
ax1.set_ylabel("Number of variants (VUS)", fontsize=11)
ax1.legend(loc="upper center")

# ---------------------------------------------------------
# GRAPHIQUE 2 : Profil le long de la protéine BRCA1 & Domaines
# ---------------------------------------------------------
ax2 = axes[1]

domains = [
    ("RING", 24, 64, "#F39C12"),
    ("Coiled-coil", 1391, 1424, "#8E44AD"),
    ("BRCT1", 1646, 1736, "#27AE60"),
    ("BRCT2", 1756, 1855, "#16A085"),
]

for name, start, end, color in domains:
    ax2.axvspan(start, end, color=color, alpha=0.2, label=f"Domaine {name}")

df_sorted_pos = df.sort_values("aa_pos")
ax2.scatter(
    df_sorted_pos["aa_pos"],
    df_sorted_pos["predicted_pathogenic_probability"],
    alpha=0.25,
    s=12,
    color="#34495E",
    label="Variants",
)

rolling_mean = (
    df_sorted_pos.set_index("aa_pos")["predicted_pathogenic_probability"]
    .rolling(window=50, min_periods=10)
    .mean()
)
ax2.plot(
    rolling_mean.index,
    rolling_mean.values,
    color="#C0392B",
    linewidth=2.5,
    label="Tendancy (sliding mean)",
)

ax2.axhline(0.5, color="grey", linestyle=":", linewidth=1)
ax2.set_xlim(1, 1863)
ax2.set_ylim(-0.02, 1.05)
ax2.set_title(
    "Profile of pathogenicity on protein BRCA1",
    fontsize=12,
    fontweight="bold",
)
ax2.set_xlabel("Position of amino acid (aa_pos)", fontsize=11)
ax2.set_ylabel("Probability of pathogenicity", fontsize=11)
ax2.legend(loc="upper left", fontsize=8, ncol=2)

plt.tight_layout()
plt.savefig(output_image, dpi=300, bbox_inches="tight")
plt.close()

print(f"Visualisation saved in : {output_image}")