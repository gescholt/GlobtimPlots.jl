module GlobtimPlotsGraphMakieExt

# GraphMakie extension for GlobtimPlots: draws a subdivision tree as a graph
# (`plot_subdivision_tree`). Node extraction, layout, styling and the legend are
# backend-free code in src/subdivision_tree_viz.jl; this module adds the graph
# (Graphs.jl), the Buchheim layout (NetworkLayout) and `graphplot!` (GraphMakie).
# GraphMakie depends on Graphs and NetworkLayout, so `using GraphMakie` alone
# loads all three triggers. The docstring is on the stub in src/extension_stubs.jl.

using GlobtimPlots
using GlobtimPlots:
    SubdomainNodeInfo,
    TreeVizStyle,
    BalancedTreeLayout,
    extract_tree_data,
    compute_auto_scale,
    compute_figure_size,
    compute_error_range,
    get_node_color,
    get_node_marker,
    get_dim_color,
    format_node_label,
    draw_depth_bands!,
    _create_tree_legend!
using Graphs: SimpleDiGraph, add_edge!, edges, src
using GraphMakie: graphplot!
using NetworkLayout: Buchheim
using Makie

#                           GRAPH CONVERSION                                    #

"""
    tree_to_graph(nodes::Vector{SubdomainNodeInfo}) -> SimpleDiGraph

Convert extracted node info to Graphs.jl SimpleDiGraph.
"""
function tree_to_graph(nodes::Vector{SubdomainNodeInfo})
    n = length(nodes)
    graph = SimpleDiGraph(n)

    for node in nodes
        if node.children !== nothing
            left_id, right_id = node.children
            add_edge!(graph, node.id, left_id)
            add_edge!(graph, node.id, right_id)
        end
    end

    return graph
end

#                           MAIN PLOTTING FUNCTION                              #

function GlobtimPlots.plot_subdivision_tree(
    tree;
    style::TreeVizStyle = TreeVizStyle(),
    title::String = "",
    show_legend::Bool = true,
    balanced_layout::Bool = true,
)
    # Extract tree data
    nodes, _, _, root_id = extract_tree_data(tree)

    # Compute auto-scaling based on tree size
    scale = compute_auto_scale(nodes, style)

    # Compute figure size
    fig_size = compute_figure_size(nodes, style)

    # Convert to graph
    graph = tree_to_graph(nodes)

    # Compute error range for gradient coloring
    error_range = compute_error_range(nodes)

    # Prepare node attributes with scaled sizes (now with error gradient)
    node_colors = [get_node_color(n, style, error_range) for n in nodes]
    node_sizes = [n.is_leaf ? scale.leaf_node_size : scale.split_node_size for n in nodes]
    node_markers = [get_node_marker(n) for n in nodes]
    node_labels = [
        format_node_label(n, nodes, style; balanced_layout = balanced_layout) for n in nodes
    ]

    # Prepare edge colors (colored by parent's split dimension)
    edge_colors = []
    for e in edges(graph)
        parent_node = nodes[src(e)]
        if parent_node.split_dim !== nothing
            push!(edge_colors, get_dim_color(parent_node.split_dim, style))
        else
            push!(edge_colors, :gray)
        end
    end

    # Create figure with computed size
    fig = Figure(size = fig_size)

    # Main axis for tree
    ax = Axis(
        fig[1, 1],
        title = isempty(title) ? "Subdivision Tree" : title,
        titlesize = style.title_fontsize,
        aspect = DataAspect(),
    )
    hidedecorations!(ax)
    hidespines!(ax)

    # Choose layout: balanced (split-aware) or standard Buchheim
    layout_obj =
        balanced_layout ?
        BalancedTreeLayout(nodes, root_id, scale.h_scale, style.vertical_spacing) :
        Buchheim()

    # Precompute positions for depth bands
    positions = layout_obj(graph)

    # Draw depth bands BEFORE the graph (so they render behind)
    draw_depth_bands!(ax, nodes, positions, style)

    # Plot the graph with scaled attributes and node strokes
    graphplot!(
        ax,
        graph,
        layout = _ -> positions,  # Use precomputed positions
        node_color = node_colors,
        node_size = node_sizes,
        node_marker = node_markers,
        node_strokewidth = style.node_strokewidth,
        node_strokecolor = style.node_strokecolor,
        edge_color = edge_colors,
        edge_width = scale.edge_width,
        nlabels = node_labels,
        nlabels_distance = max(8, round(Int, 12 * scale.font_scale)),
        nlabels_fontsize = scale.label_fontsize,
        arrow_show = false,
    )

    # Add legend if requested
    if show_legend
        _create_tree_legend!(fig, nodes, style, error_range)
    end

    return fig
end

end # module GlobtimPlotsGraphMakieExt
