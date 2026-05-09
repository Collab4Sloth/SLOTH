import pandas as pd
import matplotlib.pyplot as plt

# =========================
# User parameters
# =========================

file = "./Saves/GrainsProblem/time_specialized.csv"

col_x  = 2
col_y1 = 4
col_yf = 63

n_grains = 30

# =========================
# Read CSV file
# =========================

df = pd.read_csv(file)

# =========================
# Create figure and subplots
# =========================

fig, axes = plt.subplots(1, 2, figsize=(15, 6))

# =========================
# Styles (colors + markers)
# =========================

colors = plt.cm.tab20.colors

markers = [
    'o', 's', '^', 'D', 'v',
    '<', '>', 'p', '*', 'h',
    'H', 'X', 'd', '|', '_'
]

# =========================
# Plot 1 : Normalized fractions
# =========================

for i, col in enumerate(range(col_y1, col_y1 + 2*n_grains, 2)):

    color  = colors[i % len(colors)]
    marker = markers[i % len(markers)]

    axes[0].plot(
        df.iloc[:, col_x],
        df.iloc[:, col] / df.iloc[0, col],
        label=f'grain {i+1}',
        color=color,
        marker=marker,
        markevery=max(len(df)//20, 1),  # espace les symboles
        linewidth=1.5,
        markersize=5
    )

axes[0].set_xlabel("t")
axes[0].set_ylabel(r'Area fractions(t)/Area fractions(0)')
axes[0].set_title("Normalized Area fractions")
axes[0].grid(True)

# =========================
# Better legend placement
# =========================

axes[0].legend(
    loc='center left',
    bbox_to_anchor=(1.02, 0.5),
    fontsize=8,
    ncol=2,
    frameon=True
)

# =========================
# Plot 2 : Normalized energy density
# =========================

axes[1].plot(
    df.iloc[:, col_x],
    df.iloc[:, col_yf] / df.iloc[0, col_yf],
    color='black',
    linewidth=2
)
axes[1].set_xscale('log')
axes[1].set_xlabel("t")
axes[1].set_ylabel(r'$F(t)/F(0)$')
axes[1].set_title("Normalized energy density")
axes[1].grid(True)

# =========================
# Layout and save
# =========================

plt.tight_layout()

plt.savefig(
    "normalized_quantities.png",
    dpi=300,
    bbox_inches='tight'
)

# plt.show()