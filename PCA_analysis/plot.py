import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Ellipse

# --- 1. 定义文件和字体大小 ---
file_a1 = 'a1.csv'
file_a2 = 'a2.csv'
FONT_SIZE = 16  # 定义标签和刻度的字体大小

# 置信椭圆设置
CONF_LEVEL = 0.95  # 置信水平：0.68/0.90/0.95/0.99 常用
ELLIPSE_LW = 2.0   # 椭圆线宽

def _chi2_val_2d(conf_level: float) -> float:
    """
    返回二维高斯分布下，给定置信水平对应的卡方分布临界值（df=2）。
    优先用 scipy；如果没有 scipy，则用常见水平的硬编码回退。
    """
    try:
        from scipy.stats import chi2
        return float(chi2.ppf(conf_level, df=2))
    except Exception:
        # 常见置信水平回退（df=2）
        fallback = {
            0.68: 2.279,   # ~1-sigma in 2D
            0.90: 4.605,
            0.95: 5.991,
            0.99: 9.210,
        }
        # 若不是常见值，默认用 0.95
        return fallback.get(round(conf_level, 2), 5.991)

def add_confidence_ellipse(x, y, ax, conf_level=0.95,
                           edgecolor='black', linestyle='-',
                           linewidth=2.0, zorder=3):
    """
    给散点 (x,y) 添加二维置信椭圆（假设近似高斯）。
    椭圆中心=均值，形状由协方差决定；尺度由 chi2.ppf(conf, df=2) 决定。
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    # 去掉 NaN/Inf
    mask = np.isfinite(x) & np.isfinite(y)
    x = x[mask]
    y = y[mask]
    if x.size < 3:
        return  # 点太少不画椭圆

    cov = np.cov(x, y)
    if not np.all(np.isfinite(cov)):
        return

    # 特征分解得到主轴方向与方差
    vals, vecs = np.linalg.eigh(cov)
    order = vals.argsort()[::-1]
    vals = vals[order]
    vecs = vecs[:, order]

    # 椭圆旋转角度（弧度转角度）
    angle = np.degrees(np.arctan2(vecs[1, 0], vecs[0, 0]))

    # 置信尺度：sqrt(chi2_val)
    chi2_val = _chi2_val_2d(conf_level)
    scale = np.sqrt(chi2_val)

    # 半轴长度：scale * sqrt(eigenvalue)
    width = 2.0 * scale * np.sqrt(max(vals[0], 0.0))
    height = 2.0 * scale * np.sqrt(max(vals[1], 0.0))

    mean_x = np.mean(x)
    mean_y = np.mean(y)

    ell = Ellipse(
        xy=(mean_x, mean_y),
        width=width,
        height=height,
        angle=angle,
        facecolor='none',
        edgecolor=edgecolor,
        linestyle=linestyle,
        linewidth=linewidth,
        zorder=zorder
    )
    ax.add_patch(ell)

try:
    # --- 2. 加载数据 ---
    df1 = pd.read_csv(file_a1, header=None, names=['PC1', 'PC2'])
    df2 = pd.read_csv(file_a2, header=None, names=['PC1', 'PC2'])

    # --- 关键修正：计算对称的坐标轴限制 ---
    max_abs_x = max(df1['PC1'].abs().max(), df2['PC1'].abs().max())
    max_abs_y = max(df1['PC2'].abs().max(), df2['PC2'].abs().max())
    limit = max(max_abs_x, max_abs_y) * 1.05

    # --- 3. 创建散点图（正方形画布） ---
    fig, ax = plt.subplots(figsize=(8, 8))

    # 紫色点（更实）
    ax.scatter(df1['PC1'], df1['PC2'], s=20, alpha=1.0, color='purple', zorder=2)

    # 灰色点（更淡）
    ax.scatter(df2['PC1'], df2['PC2'], s=20, alpha=0.5, color='grey', zorder=1)

    # --- 3.1 添加置信椭圆 ---
    # 紫色组：黑色实线椭圆
    add_confidence_ellipse(
        df1['PC1'], df1['PC2'], ax,
        conf_level=CONF_LEVEL,
        edgecolor='black',
        linestyle='-',
        linewidth=ELLIPSE_LW,
        zorder=4
    )

    # 灰色组：黑色虚线椭圆
    add_confidence_ellipse(
        df2['PC1'], df2['PC2'], ax,
        conf_level=CONF_LEVEL,
        edgecolor='black',
        linestyle='--',
        linewidth=ELLIPSE_LW,
        zorder=4
    )

    # --- 4. 设置图表标签和样式 ---
    ax.set_xlim(-limit, limit)
    ax.set_ylim(-limit, limit)

    ax.set_xlabel('Principal Component 1 (PC1)', fontsize=FONT_SIZE)
    ax.set_ylabel('Principal Component 2 (PC2)', fontsize=FONT_SIZE)
    ax.tick_params(axis='both', which='major', labelsize=FONT_SIZE - 2)

    # 网格线：灰色虚线（保留你原来的设置）
    ax.grid(True, linestyle='--', alpha=0.5)

    plt.tight_layout()

    # --- 5. 保存图表 ---
    output_filename = 'alpha.svg'
    plt.savefig(output_filename, dpi=300)
    print(f"✅ 图表已成功保存为: {output_filename}")

except FileNotFoundError as e:
    print(f"❌ 错误: 找不到文件 {e.filename}。请确保文件与脚本在同一目录下。")
except Exception as e:
    print(f"❌ 绘图时发生错误: {e}")
