import sys
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import os

# Get output filename from command-line argument
output_name = sys.argv[1] if len(sys.argv) > 1 else "sankey_triple"

# Load the three clustering comparisons: Pre-Artificial, Stim-Artificial, Post-Artificial
df_pre = pd.read_csv("sankey_pre.csv")
df_stim = pd.read_csv("sankey_stim.csv")
df_post = pd.read_csv("sankey_post.csv")

def build_sankey(df, label_left='A', label_right='B', stage='1'):
    """
    Given a dataframe with two clustering columns (Ci1 and Ci2),
    build a Sankey diagram comparing the two.
    """
    # Count transitions between clusters
    counts = df.groupby(['Ci1', 'Ci2']).size().reset_index(name='count')

    # Get unique cluster labels from both sides
    source_labels = sorted(df['Ci1'].unique())
    target_labels = sorted(df['Ci2'].unique())

    # Map cluster labels to node indices
    source_map = {val: i for i, val in enumerate(source_labels)}
    target_map = {val: i + len(source_labels) for i, val in enumerate(target_labels)}

    # Build source, target, and values for the links
    source = counts['Ci1'].map(source_map)
    target = counts['Ci2'].map(target_map)
    values = counts['count']

    # Build full label list for all nodes
    labels = [f"{label_left}{val}" for val in source_labels] + [f"{label_right}{val}" for val in target_labels]

    # Return the Sankey trace
    return go.Sankey(
        node=dict(label=labels, pad=15, thickness=20),
        link=dict(source=source, target=target, value=values),
        domain=dict(x=[0,1], y=[0,1]),
        name=stage
    )

# Create a 1-row, 3-column layout for the three Sankeys
fig = make_subplots(
    rows=1,
    cols=3,
    subplot_titles=["Pre → Artificial", "Stim → Artificial", "Post → Artificial"],
    specs=[[{"type": "domain"}, {"type": "domain"}, {"type": "domain"}]]
)

# Add each comparison as a subplot
fig.add_trace(build_sankey(df_pre, 'Pre', 'Art', 'Pre'), row=1, col=1)
fig.add_trace(build_sankey(df_stim, 'Stim', 'Art', 'Stim'), row=1, col=2)
fig.add_trace(build_sankey(df_post, 'Post', 'Art', 'Post'), row=1, col=3)

# Update the figure layout
fig.update_layout(
    title_text="Natural vs Artificial Cluster Assignments",
    font_size=10
)

# Show the figure
fig.show()

# Ensure the output directory exists
output_dir = "SummaryFigures"
os.makedirs(output_dir, exist_ok=True)

# Paths to output files
output_png = os.path.join(output_dir, f"{output_name}.png")
output_html = os.path.join(output_dir, f"{output_name}.html")

# Save the figure with the given output name
fig.write_image(f"{output_name}.png", scale=2)
fig.write_html(f"{output_name}.html")

# Show and save figure
fig.show()
fig.write_image(output_png, scale=2)
fig.write_html(output_html)