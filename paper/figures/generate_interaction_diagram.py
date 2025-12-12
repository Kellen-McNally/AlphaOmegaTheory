import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches

def generate_color_diagram():
    # Setup plot
    fig, ax = plt.subplots(figsize=(12, 12))
    ax.set_aspect('equal')
    ax.axis('off')
    
    # --- GEOMETRY ---
    short_len = 1.0
    long_len = np.sqrt(3)
    
    # Angles
    angles_short = [30, 90, 150, 210, 270, 330]
    angles_long = [0, 60, 120, 180, 240, 300]
    
    # Coordinates
    def get_coord(r, theta):
        rad = np.radians(theta)
        return np.array([r * np.cos(rad), r * np.sin(rad)])

    short_roots = [get_coord(short_len, a) for a in angles_short]
    long_roots = [get_coord(long_len, a) for a in angles_long]
    
    # --- PHYSICAL MAPPING & COLORS ---
    
    # Matter (Short Roots)
    # We map specific angles to Colors (Red, Green, Blue)
    # q_r (Red) = 90 deg
    # q_g (Green) = 210 deg
    # q_b (Blue) = 330 deg
    
    # Anti-Matter (Short Roots) - Complementary Colors?
    # qbar_r (Cyan) = 270 deg (Opposite Red)
    # qbar_g (Magenta) = 30 deg (Opposite Green)
    # qbar_b (Yellow) = 150 deg (Opposite Blue)
    
    # Node Properties (Index matches angles_short list)
    # Index: 0(30), 1(90), 2(150), 3(210), 4(270), 5(330)
    
    node_props_short = {
        0: {'label': r"$\bar{q}_g$", 'color': 'magenta', 'desc': 'Anti-Green'}, # 30
        1: {'label': r"$q_r$",       'color': 'red',     'desc': 'Red'},        # 90
        2: {'label': r"$\bar{q}_b$", 'color': 'gold',    'desc': 'Anti-Blue'},  # 150
        3: {'label': r"$q_g$",       'color': 'green',   'desc': 'Green'},      # 210
        4: {'label': r"$\bar{q}_r$", 'color': 'cyan',    'desc': 'Anti-Red'},   # 270
        5: {'label': r"$q_b$",       'color': 'blue',    'desc': 'Blue'}        # 330
    }
    
    # Forces (Long Roots)
    # Colors are vector sums. 
    # g(0):   Blue + Anti-Green = Blue + Magenta = Violet/Indigo
    # g(60):  Red + Anti-Green = Red + Magenta = Rose
    # g(120): Red + Anti-Blue = Red + Yellow = Orange
    # g(180): Green + Anti-Blue = Green + Yellow = Lime
    # g(240): Green + Anti-Red = Green + Cyan = Spring Green / Teal
    # g(300): Blue + Anti-Red = Blue + Cyan = Azure
    
    node_props_long = {
        0: {'label': r"$g_{b\bar{g}}$", 'color': 'blueviolet'}, # 0
        1: {'label': r"$g_{r\bar{g}}$", 'color': 'deeppink'},   # 60
        2: {'label': r"$g_{r\bar{b}}$", 'color': 'orange'},     # 120
        3: {'label': r"$g_{g\bar{b}}$", 'color': 'lawngreen'},  # 180
        4: {'label': r"$g_{g\bar{r}}$", 'color': 'mediumspringgreen'}, # 240
        5: {'label': r"$g_{b\bar{r}}$", 'color': 'dodgerblue'}  # 300
    }

    # --- DRAWING INTERACTIONS ---
    # Draw standard triangle interactions
    interactions = [
        (5, 0, 0), # q_b + qbar_g -> g(0)
        (1, 0, 1), # q_r + qbar_g -> g(60)
        (1, 2, 2), # q_r + qbar_b -> g(120)
        (3, 2, 3), # q_g + qbar_b -> g(180)
        (3, 4, 4), # q_g + qbar_r -> g(240)
        (5, 4, 5), # q_b + qbar_r -> g(300)
    ]
    
    for (s1, s2, l) in interactions:
        p1 = short_roots[s1]
        p2 = short_roots[s2]
        pL = long_roots[l]
        
        # Color of the interaction lines follows the particle
        ax.plot([p1[0], pL[0]], [p1[1], pL[1]], color=node_props_short[s1]['color'], lw=1.5, alpha=0.4, linestyle='-')
        ax.plot([p2[0], pL[0]], [p2[1], pL[1]], color=node_props_short[s2]['color'], lw=1.5, alpha=0.4, linestyle='-')
        ax.plot([p1[0], p2[0]], [p1[1], p2[1]], 'gray', lw=0.5, alpha=0.2, linestyle=':')

    # --- DRAW NODES ---
    
    # Short
    for i, props in node_props_short.items():
        x, y = short_roots[i]
        # Halo for visibility
        circle = patches.Circle((x, y), 0.18, facecolor=props['color'], edgecolor='black', zorder=10, lw=1.5, alpha=0.8)
        ax.add_patch(circle)
        # White text for dark colors, Black for light
        text_col = 'white' if props['color'] in ['blue', 'red', 'green', 'blueviolet'] else 'black'
        ax.text(x, y, props['label'], ha='center', va='center', fontsize=14, fontweight='bold', color=text_col, zorder=11)

    # Long
    for i, props in node_props_long.items():
        x, y = long_roots[i]
        circle = patches.Circle((x, y), 0.18, facecolor=props['color'], edgecolor='black', zorder=10, lw=1.5, alpha=0.8)
        ax.add_patch(circle)
        text_col = 'white' if props['color'] in ['blueviolet', 'deeppink', 'dodgerblue'] else 'black'
        ax.text(x, y, props['label'], ha='center', va='center', fontsize=12, fontweight='bold', color=text_col, zorder=11)

    # Center
    ax.plot(0, 0, 'm*', markersize=15)
    
    # --- DECORATIONS ---
    # Outer Red Hexagon
    poly_long = patches.Polygon(long_roots, closed=True, fill=False, edgecolor='gray', linestyle='--', alpha=0.2)
    ax.add_patch(poly_long)

    # Labels for Regions
    ax.text(0, 2.3, r"Color Charge Geometry", fontsize=18, ha='center', fontweight='bold')
    ax.text(0, -2.3, r"Vector Addition = Color Conservation", fontsize=14, ha='center', style='italic')
    
    # Example Annotation
    # Arrow from q_r to g_rg
    # ax.annotate("", xy=long_roots[1]*0.9, xytext=short_roots[1]*1.1, arrowprops=dict(arrowstyle="->", color='gray'))
    
    # Save
    plt.tight_layout()
    plt.savefig('paper/figures/g2_color_interaction_final.pdf')
    print("Generated paper/figures/g2_color_interaction_final.pdf")

def generate_roots_diagram():
    print("Generating CLEAN G2 Roots diagram (Figure 1)...")
    fig, ax = plt.subplots(figsize=(10, 10))
    ax.set_aspect('equal')
    ax.axis('off')
    
    # Geometry
    short_len = 1.0
    long_len = np.sqrt(3)
    
    angles_short = [30, 90, 150, 210, 270, 330]
    angles_long = [0, 60, 120, 180, 240, 300]
    
    # Draw Arrows
    # Long Roots (Red)
    for ang in angles_long:
        rad = np.radians(ang)
        dx = long_len * np.cos(rad)
        dy = long_len * np.sin(rad)
        ax.arrow(0, 0, dx, dy, color='red', width=0.015, head_width=0.08, head_length=0.1, length_includes_head=True, alpha=0.8)
        
    # Short Roots (Blue)
    for ang in angles_short:
        rad = np.radians(ang)
        dx = short_len * np.cos(rad)
        dy = short_len * np.sin(rad)
        ax.arrow(0, 0, dx, dy, color='blue', width=0.015, head_width=0.08, head_length=0.1, length_includes_head=True, alpha=0.8)

    # Origin
    ax.plot(0, 0, 'k.', markersize=8)
    
    # Legend/Labels
    # Place a legend manually or just text
    ax.text(-2, 2, "G2 Root System", fontsize=16, fontweight='bold')
    ax.text(-2, 1.8, "Red: Long Roots (Gluons)", color='red', fontsize=12)
    ax.text(-2, 1.6, "Blue: Short Roots (Matter)", color='blue', fontsize=12)

    # Limits
    ax.set_xlim(-2.2, 2.2)
    ax.set_ylim(-2.2, 2.2)

    plt.tight_layout()
    plt.savefig('paper/figures/g2_roots.pdf')
    print("Generated paper/figures/g2_roots.pdf")

if __name__ == "__main__":
    generate_color_diagram()
    generate_roots_diagram()
