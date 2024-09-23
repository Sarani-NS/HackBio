
# Import necessary libraries
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx

# Load the dataset
data = pd.read_csv('preprocessing.csv')

# Ensure data is clean and handle missing values if necessary
data.replace(['-', '?', 'â€“', 'N/A', 'n/a'], np.nan, inplace=True)

# Step 1: Visualization 1 - Distribution of Products Across R&D Phases
rd_phase_counts = data['R&D phase'].value_counts()
rd_phase_order = ['Phase I', 'Phase II', 'Phase III', 'Preregistration']
rd_phase_counts_sorted = rd_phase_counts.reindex(rd_phase_order, fill_value=0)

# Create the bar chart for R&D Phases
plt.figure(figsize=(8, 6))
rd_phase_counts_sorted.plot(kind='bar', color='lightblue', alpha=0.8)
plt.title('Distribution of Products by R&D Phase')
plt.ylabel('Number of Products')
plt.xlabel('R&D Phase')
plt.grid(True, axis='y', linestyle='--', alpha=0.6)
plt.tight_layout()
plt.show()

# Step 2: Visualization 2 - Distribution of Product Types (Antibiotics vs Non-traditional)
product_type_counts = data['Product Type'].value_counts()

# Create the pie chart for Product Types
plt.figure(figsize=(6, 6))
product_type_counts.plot(kind='pie', autopct='%1.1f%%', startangle=90, colors=['lightcoral', 'lightblue'])
plt.title('Distribution of Product Types (Antibiotics vs Non-traditional)')
plt.ylabel('')
plt.tight_layout()
plt.show()

# Step 3: Visualization 3 - Pathogen Targeting Trends by Category
pathogen_counts = data['Pathogen category'].value_counts()

# Create the bar chart for Pathogen Targeting by Category
plt.figure(figsize=(8, 6))
pathogen_counts.plot(kind='bar', color='lightgreen', alpha=0.8)
plt.title('Pathogen Targeting by Category')
plt.ylabel('Number of Products')
plt.xlabel('Pathogen Category')
plt.grid(True, axis='y', linestyle='--', alpha=0.6)
plt.tight_layout()
plt.show()

# Step 4: Visualization 4 - Innovation vs. Product Stage (Stacked Bar Chart)
innovation_stage_data = data.groupby(['R&D phase', 'Innovative?']).size().unstack()

# Stacked Bar Chart
plt.figure(figsize=(10, 6))
innovation_stage_data.plot(kind='bar', stacked=True, color=['lightblue', 'lightpink', 'plum'])
plt.title('Innovation vs. Product Stage (Stacked Bar Chart)')
plt.ylabel('Number of Products')
plt.xlabel('R&D Phase')
plt.legend(title='Innovative Status')
plt.grid(False)
plt.tight_layout()
plt.show()

# Step 5: Visualization 5 - Transition Success Rates Between R&D Phases
# Data for transition rates
transition_rates = pd.Series({
    'Phase I → Phase II': 0.60,
    'Phase II → Phase III': 0.785,
    'Phase III → Preregistration': 0.185
})

# Bar Chart for Transition Success Rates
plt.figure(figsize=(8, 6))
transition_rates.plot(kind='bar', color=['lightblue', 'plum', 'lightcoral'], alpha=0.8)
plt.title('Transition Success Rates Between R&D Phases')
plt.ylabel('Transition Success Rate')
plt.xlabel('R&D Phase')
plt.grid(True, axis='y', linestyle='--', alpha=0.6)
plt.tight_layout()
plt.show()

# Step 6: Visualization 6 - Pathogen-Developer Relationships (Network Diagram)
# Prepare data for Network Visualization
developer_pathogen_data = data.groupby(['Developer', 'Pathogen category']).size().reset_index(name='Product Count')

# Create a graph object
G = nx.Graph()

# Add edges between developers and pathogens, weighted by the number of products
for _, row in developer_pathogen_data.iterrows():
    developer = row['Developer']
    pathogen = row['Pathogen category']
    count = row['Product Count']
    G.add_edge(developer, pathogen, weight=count)

# Set positions for the nodes
pos = nx.spring_layout(G, k=0.5, seed=42)

# Draw the graph
plt.figure(figsize=(12, 10))
nx.draw_networkx_nodes(G, pos, node_size=500, node_color='lightblue', alpha=0.8)
nx.draw_networkx_edges(G, pos, width=[G[u][v]['weight']*0.5 for u, v in G.edges()], alpha=0.6)
nx.draw_networkx_labels(G, pos, font_size=8, font_color='black')
plt.title("Top Pathogen-Developer Relationships in AMR Products")
plt.tight_layout()
plt.show()

# Step 7: Visualization 7 - Bubble Chart for Pathogen Targeting by Product Type
# Prepare data for the bubble chart
bubble_data = data.groupby(['Pathogen category', 'Product Type']).size().unstack()

# Bubble chart
plt.figure(figsize=(8, 6))
for i in range(len(bubble_data.index)):
    for j in range(len(bubble_data.columns)):
        plt.scatter(j, i, s=bubble_data.iloc[i, j] * 100, alpha=0.5)
plt.xticks(range(len(bubble_data.columns)), bubble_data.columns)
plt.yticks(range(len(bubble_data.index)), bubble_data.index)
plt.title('Pathogen Targeting by Product Type (Bubble Chart)')
plt.xlabel('Product Type')
plt.ylabel('Pathogen Category')
plt.tight_layout()
plt.show()

# Step 8: Visualization 8 - Heatmap for Correlation between Product Type and Pathogen Category
# Prepare correlation data
heatmap_data = data.pivot_table(index='Product Type', columns='Pathogen category', aggfunc='size', fill_value=0)

# Create the heatmap
plt.figure(figsize=(8, 6))
sns.heatmap(heatmap_data, annot=True, cmap='coolwarm', linewidths=.5)
plt.title('Correlation between Product Type and Pathogen Category')
plt.tight_layout()
plt.show()
