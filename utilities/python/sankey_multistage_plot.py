import pandas as pd
import plotly.graph_objects as go
import sys
import os

# Handle command-line argument for filename
output_name = sys.argv[1] if len(sys.argv) > 1 else "sankey_multistage"
output_dir = "SummaryFigures"
os.makedirs(output_dir, exist_ok=True)

output_png = os.path.join(output_dir, f"{output_name}.png")
output_html = os.path.join(output_dir, f"{output_name}.html")

# Load 3-column cluster table: columns = Pre, Stim, Post
df = pd.read_csv("multi_stage_clusters.csv")

# Initialize node and link info
all_labels = []
source = []
target = []
value = []

# We'll process transitions: Pre → Stim, then Stim → Post
n_stages = df.shape[1]
label_offset = 0  # To assign unique node indices across stages
label_maps = []   # Store label-to-index maps per stage

# Create label maps for each stage
for col in df.columns:
    labels = sorted(df[col].unique())
    label_map = {val: i + label_offset for i, val in enumerate(labels)}
    all_labels.extend([f"{col}_{val}" for val in labels])
    label_maps.append(label_map)
    label_offset += len(labels)

# Build links between consecutive stages
for i in range(n_stages - 1):
    col1 = df.columns[i]
    col2 = df.columns[i+1]
    map1 = label_maps[i]
    map2 = label_maps[i+1]
    
    # Count transitions
    grouped = df.groupby([col1, col2]).size().reset_index(name='count')
    
    for _, row in grouped.iterrows():
        src = map1[row[col1]]
        tgt = map2[row[col2]]
        source.append(src)
        target.append(tgt)
        value.append(row['count'])

# Build Sankey
fig = go.Figure(data=[go.Sankey(
    node=dict(label=all_labels, pad=15, thickness=20),
    link=dict(source=source, target=target, value=value)
)])

fig.update_layout(
    # title_text="Pre → Stim → Post Cluster Transitions",
    font_size=10
)

# Show and save
fig.show()
fig.write_image(output_png, scale=2)
fig.write_html(output_html)