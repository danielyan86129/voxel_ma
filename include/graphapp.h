#pragma once
#include <boost/graph/undirected_graph.hpp>
#include <string>
#include <vector>
#include "commondefs.h"
#include "geomalgo.h"

using std::string;
using std::vector;
namespace bg = boost::graph;

namespace graphapp
{
	struct GraphNode
	{
		int idx; // node index
		point p; // position
	};
	struct GraphEdge
	{
		int idx; // edge index
		double w; // edge weight
	};
	// typedef boost::property<boost::edge_weight_t, double> EdgeWeightProp;
	typedef boost::undirected_graph<GraphNode, GraphEdge> WeightedGraph;
	typedef boost::graph_traits<WeightedGraph>::vertex_descriptor NodeHandle;
	typedef boost::graph_traits<WeightedGraph>::edge_descriptor EdgeHandle;
	typedef boost::graph_traits<WeightedGraph>::vertex_iterator NodeIter;
	typedef boost::graph_traits<WeightedGraph>::edge_iterator EdgeIter;
	// different weights (to compute/use)
	enum Weight 
	{
		BT3 = 0, BT2 = 1, BT1 = 2
	};
	// method to use to tree-ify a graph
	enum TreeMethod 
	{
		ShortestPath
	};
	// read nodes and edges info from a file
	bool readGraph( const string & _filename,
		vector<point>& _nodes, vector<ivec2>& _edges,
		vector<vector<float>>& _msures );
	// construct a boost graph from file
	bool readGraph( const string& _filename,
		WeightedGraph& _bg,
		vector<vector<float>>& _msures );
	//
	// compute edge weights for the given graph connectivity
	void computeWeights(
		const vector<point>& _nodes, const vector<ivec2>& _edges,
		const vector<vector<float>>& _msures,
		Weight wt_type,
		vector<float>& _weights );
	//
	// construct a weighted graph from given weights and connectivity
	bool makeBoostGraph(
		const vector<point>& _nodes, const vector<ivec2>& _edges,
		const vector<float>& _weights,
		WeightedGraph& _g );
	// 
	// return a tree from the graph, using specified method
	void makeTreeFromGraph(
		const WeightedGraph & _g,
		const vector<vector<float>>& _edge_msures,
		TreeMethod _method,
		vector<NodeHandle> & _t );
	//
	// convert a parent list to a list of edges
	void convertToEdges( const WeightedGraph& _g, const vector<NodeHandle>& _t, vector<ivec2>& _edge_list );
	//
	// write tree to file
	void exportTree(
		const WeightedGraph& _g,
		vector<NodeHandle>& _t,
		const std::string& _filename );
}