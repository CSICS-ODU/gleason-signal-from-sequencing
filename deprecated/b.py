# bcrz_comparison.py
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrowPatch, Circle

def panel_left(ax):
    ax.set_title("A) Today: Credential-Dominant ZTA + Overlay Analytics", fontsize=10, pad=8)
    # User
    ax.add_patch(Circle((0.08,0.55),0.03, edgecolor='black', facecolor='none', linewidth=1.5))
    ax.text(0.08,0.46,"User", ha='center', va='center', fontsize=8)
    # ZTA
    ax.add_patch(Rectangle((0.18,0.45),0.18,0.20, fill=False, edgecolor='black', linewidth=2))
    ax.text(0.27,0.55,"ZTA\n(credential)", ha='center', va='center', fontsize=9)
    # Behavioral overlay (dashed)
    ax.add_patch(Rectangle((0.20,0.72),0.16,0.12, fill=False, edgecolor='gray', linewidth=1.5, linestyle='--'))
    ax.text(0.28,0.78,"Behavioral\nanalytics", ha='center', va='center', fontsize=8, color='gray')
    # Arrows
    ax.add_patch(FancyArrowPatch((0.11,0.55),(0.18,0.55), arrowstyle='-|>', mutation_scale=12))
    ax.add_patch(FancyArrowPatch((0.28,0.72),(0.28,0.65), arrowstyle='-|>', mutation_scale=12, color='gray', linestyle='--'))
    # Devices
    dev_x = [0.52,0.62,0.72]
    for i,x in enumerate(dev_x):
        ax.add_patch(Rectangle((x,0.50-0.05),0.08,0.10, fill=False, edgecolor='black'))
        ax.text(x+0.04,0.50,f"D{i+1}", ha='center', va='center', fontsize=9)
        ax.add_patch(FancyArrowPatch((0.36,0.55),(x,0.55), arrowstyle='-|>', mutation_scale=12, color='black'))
    # “revoke user” red X
    ax.text(0.16,0.40,"Default under doubt:\nRevoke user", fontsize=8, color='crimson', ha='left', va='center')
    ax.plot([0.14,0.22],[0.44,0.36], color='crimson', linewidth=2)
    ax.plot([0.14,0.22],[0.36,0.44], color='crimson', linewidth=2)
    ax.set_xlim(0,1); ax.set_ylim(0,1); ax.axis('off')

def panel_right(ax):
    ax.set_title("B) BCRZ: Dual-Identity ZTA (Cut the Edge, Not the User)", fontsize=10, pad=8)
    # Credential pool
    ax.add_patch(Rectangle((0.08,0.18),0.26,0.58, fill=False, edgecolor='purple', linewidth=2))
    ax.text(0.21,0.80,"Credential Identity Pool (C-ID)", ha='center', va='center', fontsize=9, color='purple')
    cid_y = [0.62,0.50,0.38]
    for i,y in enumerate(cid_y):
        ax.add_patch(Rectangle((0.12,y-0.04),0.08,0.08, facecolor='lavender', edgecolor='purple'))
        ax.text(0.16,y,f"C-ID {chr(65+i)}", ha='center', va='center', fontsize=8)
    # New C-ID X
    ax.add_patch(Rectangle((0.12,0.24),0.08,0.08, facecolor='white', edgecolor='purple', linestyle='--'))
    ax.text(0.16,0.28,"New\nC-ID X", ha='center', va='center', fontsize=7)

    # Behavioral pool
    ax.add_patch(Rectangle((0.66,0.18),0.26,0.58, fill=False, edgecolor='green', linewidth=2))
    ax.text(0.79,0.80,"Behavioral Identity Pool (B-ID)", ha='center', va='center', fontsize=9, color='green')
    bid_y = [0.62,0.50,0.38]; labels = ['α','β','γ']
    for i,y in enumerate(bid_y):
        ax.add_patch(Rectangle((0.70,y-0.04),0.08,0.08, facecolor='honeydew', edgecolor='green'))
        ax.text(0.74,y,f"B-ID {labels[i]}", ha='center', va='center', fontsize=8)

    # Policy engine
    ax.add_patch(Rectangle((0.40,0.82),0.20,0.12, fill=False, edgecolor='black', linewidth=2))
    ax.text(0.50,0.88,"Policy Engine (π*)", ha='center', va='center', fontsize=9, fontweight='bold')
    ax.text(0.50,0.84,"ADMIT • CONSTRAIN • MFA • REVOKE", ha='center', va='center', fontsize=7)
    # Inputs to policy
    ax.add_patch(FancyArrowPatch((0.34,0.77),(0.40,0.88), arrowstyle='-|>', mutation_scale=12, color='purple'))
    ax.add_patch(FancyArrowPatch((0.66,0.77),(0.60,0.88), arrowstyle='-|>', mutation_scale=12, color='green'))

    # Mappings (aligned)
    for y in cid_y:
        ax.add_patch(FancyArrowPatch((0.20,y),(0.70,y), arrowstyle='-|>', mutation_scale=12, linewidth=1.5, color='black'))
    # Decayed link (new cred to old behavior)
    ax.add_patch(FancyArrowPatch((0.20,0.28),(0.70,0.38), arrowstyle='-|>', mutation_scale=12, linewidth=1.2, color='gray', linestyle='--'))
    ax.text(0.45,0.41,"soft link p(t)", fontsize=7, ha='center', va='center')

    # Continuity path (behavioral track)
    ax.add_patch(Rectangle((0.42,0.05),0.16,0.08, fill=False, edgecolor='black'))
    ax.text(0.50,0.09,"Behavioral Track\n(constrained continuity)", ha='center', va='center', fontsize=7)
    ax.add_patch(FancyArrowPatch((0.50,0.82),(0.50,0.13), arrowstyle='-|>', mutation_scale=12, color='black'))

    # Red X to show “cut mapping, not user”
    ax.plot([0.46,0.54],[0.50,0.42], color='crimson', linewidth=2)
    ax.plot([0.46,0.54],[0.42,0.50], color='crimson', linewidth=2)
    ax.text(0.50,0.46,"cut mapping,\nnot user", fontsize=7, color='crimson', ha='center')

    ax.set_xlim(0,1); ax.set_ylim(0,1); ax.axis('off')

def draw_comparison(outfile="BCRZ_comparison.png", dpi=300):
    fig, axs = plt.subplots(1,2, figsize=(13,5))
    panel_left(axs[0]); panel_right(axs[1])
    plt.tight_layout()
    plt.savefig(outfile, dpi=dpi, bbox_inches='tight')
    print(f"Saved {outfile}")

if __name__ == "__main__":
    draw_comparison()
