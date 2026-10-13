
import streamlit as st
import streamlit.components.v1 as components
import networkx as nx
import pandas as pd
import requests
from pyvis.network import Network
from networkx.algorithms.community import greedy_modularity_communities

st.set_page_config(page_title="IPF Protein Network Explorer", layout="wide")

st.title("IPF Protein Interaction Network Explorer")
st.write("Explore functional protein associations related to idiopathic pulmonary fibrosis (IPF).")

# Proteins associated with fibrosis
proteins = [
    "TGFB1", "COL1A1", "COL1A2", "ACTA2", "FN1",
    "CTGF", "MMP2", "MMP9", "TIMP1", "SERPINE1",
    "SMAD2", "SMAD3", "SMAD4", "IL6", "TNF",
    "PDGFA", "PDGFB", "VIM", "SPP1", "CXCL8"
]

# Retrieve real protein associations from STRING
@st.cache_data
def load_network():
    url = "https://string-db.org/api/json/network"
    params = {
        "identifiers": "\r".join(proteins),
        "species": 9606,
        "required_score": 400,
        "network_type": "functional"
    }

    response = requests.get(url, params=params, timeout=30)
    response.raise_for_status()
    data = response.json()

    G = nx.Graph()

    for protein in proteins:
        G.add_node(protein)

    for interaction in data:
        source = interaction["preferredName_A"]
        target = interaction["preferredName_B"]

        if source in proteins and target in proteins:
            G.add_edge(source, target, weight=interaction["score"])

    return G

try:
    G = load_network()
except (requests.RequestException, ValueError, KeyError) as error:
    st.error("Could not load data from STRING. Please try again later.")
    st.stop()

if G.number_of_edges() == 0:
    st.warning("No interactions were returned by STRING.")
    st.stop()

# Network measurements
degrees = dict(G.degree())
communities = list(greedy_modularity_communities(G))

community_map = {}
for i, community in enumerate(communities):
    for node in community:
        community_map[node] = i

palette = [
    "#2E86AB", "#E07A5F", "#81B29A",
    "#F2CC8F", "#9B5DE5", "#F15BB5",
    "#00BBF9", "#6D597A"
]

# Interactive controls
st.sidebar.header("Network Controls")

layout = st.sidebar.selectbox(
    "Choose Layout",
    ["Force-directed", "Circular"]
)

min_degree = st.sidebar.slider(
    "Minimum Degree",
    min_value=0,
    max_value=max(degrees.values()),
    value=0
)

# Filter nodes using original network degree
selected_nodes = [
    node for node in G.nodes()
    if degrees[node] >= min_degree
]

H = G.subgraph(selected_nodes).copy()

# Network summary
col1, col2, col3 = st.columns(3)
col1.metric("Proteins", H.number_of_nodes())
col2.metric("Interactions", H.number_of_edges())
col3.metric("Communities", len(set(community_map[n] for n in H.nodes())))

st.subheader("Interactive Protein Network")

if H.number_of_nodes() == 0:
    st.info("No proteins match the selected degree threshold.")
else:
    net = Network(
        height="600px",
        width="100%",
        bgcolor="#FFFFFF",
        font_color="#222222",
        cdn_resources="in_line"
    )

    for node in H.nodes():
        degree = degrees[node]
        community = community_map[node]

        net.add_node(
            node,
            label=node,
            size=10 + degree * 1.5,
            color=palette[community % len(palette)],
            title=f"Protein: {node}<br>Degree: {degree}<br>Community: {community + 1}"
        )

    for source, target, attributes in H.edges(data=True):
        net.add_edge(
            source,
            target,
            color="#CCCCCC",
            title=f"STRING score: {attributes['weight']:.3f}"
        )

    # Compute actual positions for each layout
    if layout == "Circular":
        positions = nx.circular_layout(H, scale=350)
    else:
        positions = nx.spring_layout(H, seed=42, scale=350)

    for node in net.nodes:
        node_id = node["id"]
        node["x"] = float(positions[node_id][0])
        node["y"] = float(positions[node_id][1])
        node["physics"] = False

    net.toggle_physics(False)

    components.html(net.generate_html(notebook=False), height=620)

st.subheader("About the Network")

st.markdown("""
**Network Description**

This network represents functional associations among human proteins
related to idiopathic pulmonary fibrosis (IPF). The interaction data
are retrieved from the STRING database.

**Visual Encodings**

- Node size represents degree (number of connections).
  Larger nodes indicate proteins with more associations.
- Node color represents communities detected using
  greedy modularity. Different colors highlight groups
  of closely connected proteins.
- The layout selector allows comparison between
  force-directed and circular arrangements.
- The minimum degree slider filters proteins based on
  their number of connections in the original network.

**Layout Comparison**

The force-directed layout places connected proteins closer
together, making clusters easier to identify. The circular
layout provides a more organized arrangement but makes
community structure harder to recognize.

**Limitation**

Dense networks can become difficult to interpret because
edges overlap. Also, STRING associations may represent
functional relationships rather than direct physical
interactions, and they do not necessarily indicate
IPF-specific activity.
""")

st.caption("Data source: STRING database | Organism: Homo sapiens | Minimum STRING score: 0.400")
