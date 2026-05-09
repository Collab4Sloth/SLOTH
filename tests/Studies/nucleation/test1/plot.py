import pandas as pd
import matplotlib.pyplot as plt

# =========================
# User parameters
# =========================

file_099 = "./Saves_099/Problem1/time_specialized.csv"
file_100 = "./Saves_1/Problem1/time_specialized.csv"
file_101 = "./Saves_101/Problem1/time_specialized.csv"

col_x = 2
col_y1 = 4
col_y2 = 5

# =========================
# Read CSV files
# =========================

df_099 = pd.read_csv(file_099)
df_100 = pd.read_csv(file_100)
df_101 = pd.read_csv(file_101)

# Column names
x_name  = df_099.columns[col_x]
y1_name = df_099.columns[col_y1]
y2_name = df_099.columns[col_y2]

# =========================
# Create figure and subplots
# =========================

fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# =========================
# Plot 1 : Normalized fraction
# =========================

axes[0].plot(
    df_099.iloc[:, col_x],
    df_099.iloc[:, col_y1] / df_099.iloc[0, col_y1],
    label=r'$r_o = 0.99\, r_c$'
)

axes[0].plot(
    df_100.iloc[:, col_x],
    df_100.iloc[:, col_y1] / df_100.iloc[0, col_y1],
    label=r'$r_o = r_c$'
)

axes[0].plot(
    df_101.iloc[:, col_x],
    df_101.iloc[:, col_y1] / df_101.iloc[0, col_y1],
    label=r'$r_o = 1.01\, r_c$'
)

axes[0].set_xlabel("t")
axes[0].set_ylabel(r'$Y(t)/Y(0)$')
axes[0].set_title("Normalized fraction")
axes[0].grid(True)
axes[0].legend()

# =========================
# Plot 2 : Normalized energy density
# =========================

axes[1].plot(
    df_099.iloc[:, col_x],
    df_099.iloc[:, col_y2] / df_099.iloc[0, col_y2],
    label=r'$r_o = 0.99\, r_c$'
)

axes[1].plot(
    df_100.iloc[:, col_x],
    df_100.iloc[:, col_y2] / df_100.iloc[0, col_y2],
    label=r'$r_o = r_c$'
)

axes[1].plot(
    df_101.iloc[:, col_x],
    df_101.iloc[:, col_y2] / df_101.iloc[0, col_y2],
    label=r'$r_o = 1.01\, r_c$'
)

axes[1].set_xlabel("t")
axes[1].set_ylabel(r'$F(t)/F(0)$')
axes[1].set_title("Normalized energy density")
axes[1].grid(True)
axes[1].legend()

# =========================
# Adjust layout and save figure
# =========================

plt.tight_layout()

# Save figure as PNG
plt.savefig("normalized_quantities.png", dpi=300, bbox_inches='tight')

# # Display figure
# plt.show()