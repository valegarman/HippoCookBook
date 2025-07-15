import os
import sys
import pandas as pd
import plotly.graph_objects as go

# Get output filename from command line arguments
# Default name if nothing is provided
if len(sys.argv) > 1:
    output_name = sys.argv[1]
else:
    output_name = "sankey_output.png"

# Load cluster assignments
df = pd.read_csv("cluster_assignments.csv")
counts = df.groupby(['Ci1', 'Ci2']).size().reset_index(name='count')

# Unique labels and maps
source_labels = sorted(df['Ci1'].unique())
target_labels = sorted(df['Ci2'].unique())

source_map = {val: i for i, val in enumerate(source_labels)}
target_map = {val: i + len(source_labels) for i, val in enumerate(target_labels)}

source = counts['Ci1'].map(source_map)
target = counts['Ci2'].map(target_map)
values = counts['count']

labels = [f"A{val}" for val in source_labels] + [f"B{val}" for val in target_labels]

# Build Sankey diagram
fig = go.Figure(data=[go.Sankey(
    node=dict(label=labels, pad=15, thickness=20),
    link=dict(source=source, target=target, value=values)
)])

# Create output directory
output_dir = "SummaryFigures"
os.makedirs(output_dir, exist_ok=True)

# Final path
output_path = os.path.join(output_dir, output_name)
fig.write_image(output_path, width=1000, height=600, scale=3)