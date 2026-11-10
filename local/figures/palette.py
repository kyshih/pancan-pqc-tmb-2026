COLORS = {
    # Main conditions
    "high":        "#D55E00",  # vermillion - for high TMB, upregulated, etc
    "low":         "#0072B2",  # blue - for low TMB, downregulated, etc
    
    # Significance
    "sig_strong":  "#CC0000",  # FDR < 0.05 dont use this anywhere else - strong red, use for strongest significance
    "ns":          "#999999",  # not significant - mid gray
    
    # Additional categories
    "accent1":     "#009E73",  # teal green
    "accent2":     "#CC79A7",  # mauve/pink
    "accent3":     "#56B4E9",  # sky blue
    "accent4":     "#E69F00",  # amber 
    "accent5":     "#984ea3",  # purple
    "accent6":     "#F0E442",  # yellow (use sparingly, low contrast on white)
    
    # Neutral
    "light_gray":  "#CCCCCC",  # backgrounds, error bars
    "dark_gray":   "#333333",  # text, axes if not black
    "black":       "#000000",  # text
}

# fig1b_kegg_ora_all_tmb_genes.Rmd — functional category colors
FIG1B_KEGG_GROUP_COLORS = {
    "Folding, sorting and degradation": "#D55E00",
    "Translation": "#0072B2",
    "Transcription": "#009E73",
    "DNA replication and repair": "#984ea3",
    "Cell growth and death": "#CC79A7",
    "Metabolism": "#E69F00",
    "Immune": "#56B4E9",
    "Other": "#999999",
}

FIG1B_KEGG_GROUP_ORDER = [
    "Folding, sorting and degradation",
    "Translation",
    "Transcription",
    "DNA replication and repair",
    "Cell growth and death",
    "Metabolism",
    "Immune",
    "Other",
]


def rasterize_figure(fig):
    """Rasterize line/marker artists in vector PDF exports; text stays vector."""
    from matplotlib.collections import LineCollection, PathCollection, PolyCollection
    from matplotlib.lines import Line2D

    for ax in fig.axes:
        for artist in ax.get_children():
            if isinstance(artist, (Line2D, PathCollection, LineCollection, PolyCollection)):
                artist.set_rasterized(True)