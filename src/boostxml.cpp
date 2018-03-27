#include<iostream>
#include <utility>                   // for std::pair
#include <algorithm>                 // for std::for_each
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>

#include<boost/property_tree/ptree.hpp>
#include<boost/property_tree/xml_parser.hpp>

using namespace boost;
using namespace std;
struct VertexData
{
    int node_id;
};
struct EdgeData
{
    double param[7];
};

void usage()
{
    cout<<"argv[1]:"<<"*.xml file path"<<endl;
}

void build_transform_tree(typename boost::adjacency_list<boost::vecS,
                          boost::vecS,
                          boost::bidirectionalS,
                          VertexData,
                          boost::property<boost::edge_weight_t, double, EdgeData>
                          > &g, std::list<int> &transform_tree)
{
 return;
}



int main(int argc, char** argv)
{
    if(argc != 2)
        usage();
std:;string map_filepath = argv[1];
    boost::property_tree::ptree pt;
    boost::property_tree::read_xml(map_filepath, pt);
    if (pt.empty())
    {
        printf("Load %s map.xml failed!", map_filepath.c_str());
        return false;
    }

    // create a typedef for the Graph type
    typedef boost::adjacency_list<boost::vecS,
            boost::vecS,
            boost::bidirectionalS,
            VertexData,
            boost::property<boost::edge_weight_t, double, EdgeData>
            > Graph;
    typedef std::pair<int, int> Edge;
    // declare a graph object
    Graph g;

    std::set<int> node_set;
    std::map<int, boost::graph_traits<Graph>::vertex_descriptor> map_vertex_descriptor;
    boost::property_tree::ptree p_node_list =
            pt.get_child("map.node_list");
    for(auto p_node_it = p_node_list.begin(); p_node_it != p_node_list.end(); p_node_it++)
    {
        // 获取node_id
        int node_id = p_node_it->second.get<int>("id");
        // add vertex property
        if(node_set.find(node_id) == node_set.cend())  // have not add in graph
        {
            VertexData vertexData = {0};
            vertexData.node_id = node_id;
            map_vertex_descriptor[node_id] = boost::add_vertex(vertexData,g);

            node_set.insert(node_id);
        }
        //add_vertex(node_id, g);
        //cout<<"node_id:"<<node_id<<" has neighbours: ";

        boost::property_tree::ptree p_neighbour_node_list =
                p_node_it->second.get_child("neighbour_list");

        // 迭代邻居
        for(auto p_neighbour_node_it = p_neighbour_node_list.begin(); p_neighbour_node_it != p_neighbour_node_list.end(); p_neighbour_node_it++)
        {
            int neighbour_id = p_neighbour_node_it->second.get<int>("id");
            std::string param = p_neighbour_node_it->second.get<std::string>("transform");
            std::string value;
            EdgeData transform_param = {0};

            // split param
            string::size_type pos1, pos2;
            char split_type = ' ';
            int index = 0;
            pos2 = param.find(split_type);
            pos1 = 0;
            while(string::npos != pos2)
            {
                value = param.substr(pos1,pos2-pos1);
                pos1 = pos2 + 1;
                pos2 = param.find(split_type, pos1);
                transform_param.param[index++] = std::stod(value);
                //cout<<transform_param.param[index-1]<<" ";
            }
            if(pos1 != param.length())
            {
                value = param.substr(pos1);
                transform_param.param[index] = std::stod(value);
                //cout<<transform_param.param[index]<<endl;
            }

            // add vertex property
            if(node_set.find(neighbour_id) == node_set.cend())  // have not add in graph
            {
                VertexData vertexData = {0};
                vertexData.node_id = neighbour_id;
                map_vertex_descriptor[neighbour_id] = boost::add_vertex(vertexData,g);
                node_set.insert(neighbour_id);
            }
            // add edge
            // first is the edge. second is a bool telling you whether this is a new edge
            // or an existing one.
            auto e = boost::add_edge(map_vertex_descriptor[neighbour_id], map_vertex_descriptor[node_id], g).first;

            g[e] = transform_param;
        }
    }
    // add weight for edges
    boost::property_map<Graph, boost::edge_weight_t>::type weightmap =
            get(boost::edge_weight, g);
    typedef boost::graph_traits<Graph>::edge_iterator  edge_iter ;
    std::pair<edge_iter, edge_iter> eip;
    eip=boost::edges(g);
    for(auto it = eip.first; it != eip.second; it++)
        weightmap[*it] = 1.0;

    // find the transform param

    typedef boost::graph_traits <Graph>::vertex_descriptor vertex_descriptor;
    std::vector<vertex_descriptor> p(boost::num_vertices(g));
    std::vector<double> d(boost::num_vertices(g));
    boost::dijkstra_shortest_paths(g, map_vertex_descriptor[0],
            boost::predecessor_map(boost::make_iterator_property_map(p.begin(), get(boost::vertex_index, g))).
            distance_map(boost::make_iterator_property_map(d.begin(), get(boost::vertex_index, g))));
    //   auto vertex_property = boost::get(&VertexData::node_id,g);
    //    boost::graph_traits < Graph >::vertex_iterator vi, vend;
    //    for (boost::tie(vi, vend) = boost::vertices(g); vi != vend; ++vi) {
    //      std::cout << "node 0 to " << vertex_property[*vi] << " distance = " << d[*vi] << ", "
    //      << "parent(" << vertex_property[*vi] << ") = " << vertex_property[p[*vi]] << std::endl;
    //   }

    // build the transform tree
    auto vertex_property = boost::get(&VertexData::node_id,g);
    boost::graph_traits < Graph >::vertex_iterator vi, vend;
    std::list<std::list<int>> transform_trees;
    int debug_idx = 0;
    for (boost::tie(vi, vend) = boost::vertices(g); vi != vend; ++vi)
    {
       std::list<int> transform_tree;
        auto cur_vertex = *vi;
        if(vertex_property[cur_vertex] != vertex_property[map_vertex_descriptor[0]])
            transform_tree.push_back(vertex_property[cur_vertex]);

        while(1)
        {
            if(vertex_property[cur_vertex] != vertex_property[map_vertex_descriptor[0]])
            {
                transform_tree.push_front(vertex_property[p[cur_vertex]]);
                cur_vertex = p[cur_vertex];
            }
            else
            {
                debug_idx = 0;
                break;
            }
        }
        transform_trees.push_back(transform_tree);
    }
    //
//    for(auto it = transform_trees.cbegin(); it != transform_trees.cend(); it++)
//    {
//        auto transform_tree = *it;
//         for(auto it_t = transform_tree.cbegin(); it_t != transform_tree.cend(); it_t++)
//         {
//             cout<<*it_t<<"-";
//         }
//         cout<<endl;
//    }




    //遍历边的属性
    //typedef boost::graph_traits<Graph>::edge_iterator  edge_iter ;
    //pair<edge_iter, edge_iter> eip;
    //    eip=boost::edges(g);
    //    auto p_propery_edge = boost::get(&EdgeData::param, g);
    //    for(auto it = eip.first; it != eip.second; it++)
    //    {
    //        cout<<p_propery_edge[*it][0]
    //                <<" "<<p_propery_edge[*it][1]
    //                <<" "<<p_propery_edge[*it][2]
    //                <<" "<<p_propery_edge[*it][3]
    //                <<" "<<p_propery_edge[*it][4]
    //                <<" "<<p_propery_edge[*it][5]
    //                <<" "<<p_propery_edge[*it][6]
    //                <<"  weight:"<<weightmap[*it]<<endl;
    //    }

    //3.遍历顶点和属性
    //    typedef boost::graph_traits<Graph>::vertex_iterator vertex_iter;
    //    pair<vertex_iter, vertex_iter> vip;
    //    cout << "Vertices in g  = [ ";
    //    vip = boost::vertices(g);
    //    auto vertex_property = boost::get(&VertexData::node_id,g);
    //    for(vertex_iter vi = vip.first; vi != vip.second; ++vi) {
    //        cout << vertex_property[*vi] << " ";
    //    }
    //    cout<<"]"<<endl;
    // 遍历边
    //    typedef graph_traits<Graph>::edge_iterator  edge_iter ;
    //       pair<edge_iter, edge_iter> eip;
    //       eip=edges(g);
    //       cout << "Edge in g  = [ ";
    //       for(edge_iter ei = eip.first; ei != eip.second; ++ei) {
    //           //cout << *ei << " ";
    //           cout<<"( source edge="<< source(*ei, g) ;
    //           cout<< " taget edge="<<target(*ei, g) <<")"<<endl;
    //       }
    //       cout<<"]"<<endl;

    // test source and target in edge
    //    auto v_property = boost::get(&VertexData::node_id,g);
    //    VertexData d0,d1;
    //    d0.node_id = 200;
    //    d1.node_id = 201;
    //    auto a0 = boost::add_vertex(d0,g);
    //    auto a1 = boost::add_vertex(d1,g);

    //    auto ex = boost::add_edge(a0,a1,g).first;
    //    auto v0 = boost::source(ex,g);
    //    auto v1 = boost::target(ex,g);
    //    cout<<"source:"<<v_property[v0]
    //          <<"target"<<v_property[v1]<<endl;


}

