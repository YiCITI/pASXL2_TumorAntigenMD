import pandas as pd
import matplotlib.pyplot as plt

file_a1 = "b1.csv"
file_a2 = "b2.csv"

AXIS_LABEL_FONTSIZE = 28   
TICK_LABEL_FONTSIZE = 22   

TICK_LENGTH = 6
TICK_WIDTH  = 1.2


FIG_SIZE_INCH = 8  


PAD_RATIO = 0.05   #

try:
  
    df1 = pd.read_csv(file_a1, header=None, names=["PC1", "PC2"])
    df2 = pd.read_csv(file_a2, header=None, names=["PC1", "PC2"])

    
    max_abs_x = max(df1["PC1"].abs().max(), df2["PC1"].abs().max())
    max_abs_y = max(df1["PC2"].abs().max(), df2["PC2"].abs().max())
    limit = max(max_abs_x, max_abs_y) * (1.0 + PAD_RATIO)

  
    fig, ax = plt.subplots(figsize=(FIG_SIZE_INCH, FIG_SIZE_INCH))

    ax.scatter(df1["PC1"], df1["PC2"], s=20, alpha=1.0, color="purple")
    ax.scatter(df2["PC1"], df2["PC2"], s=20, alpha=0.5, color="grey")

    ax.set_xlim(-limit, limit)
    ax.set_ylim(-limit, limit)

    ax.set_aspect("equal", adjustable="box")

    ax.set_xlabel("PC1", fontsize=AXIS_LABEL_FONTSIZE)
    ax.set_ylabel("PC2", fontsize=AXIS_LABEL_FONTSIZE)

    ax.tick_params(
        axis="both",
        which="major",
        labelsize=TICK_LABEL_FONTSIZE,
        length=TICK_LENGTH,
        width=TICK_WIDTH,
    )

    ax.grid(True, linestyle="--", alpha=0.5)

    fig.tight_layout()


    output_filename = "beta.svg"
    fig.savefig(output_filename, dpi=300)
    print(f"✅ figure saved: {output_filename}")

except FileNotFoundError as e:
    print(f"❌ error: {e}")
except Exception as e:
    print(f"error: {e}")

