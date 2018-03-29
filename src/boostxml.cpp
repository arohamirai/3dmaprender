#include<iostream>
#include <utility>                   // for std::pair
#include <algorithm>                 // for std::for_each
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>

#include<boost/property_tree/ptree.hpp>
#include<boost/property_tree/xml_parser.hpp>

#include <vtkGenericDataObjectReader.h>
#include <vtkStructuredGrid.h>
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <string>
#include <vtkPolyDataMapper.h>
#include <vtkStructuredGridReader.h>
#include <vtkStructuredGridGeometryFilter.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkSmoothPolyDataFilter.h>
#include <vtkPlane.h>
#include <math.h>
#include <limits>
#include <vtkPerspectiveTransform.h>
#include <vtkTransform.h>
#include <vtkMatrix4x4.h>
#include <vtkAppendPolyData.h>
#include <vtkPLYWriter.h>
#include <vtkQuaternion.h>
#include <opencv2/opencv.hpp>
#include <opencv2/core.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/highgui.hpp>
#define RESOLUTION 0.5

using namespace boost;
using namespace std;
using namespace cv;
struct VertexData
{
    int node_id;
};
struct EdgeData
{
    double param[7];
};

int main(int argc, char** argv)
{
std:;string map_filepath = "../frames/map.xml";
    boost::property_tree::ptree pt;
    boost::property_tree::read_xml(map_filepath, pt);
    if (pt.empty())
    {
        printf("Load %s map.xml failed!", map_filepath.c_str());
        return false;
    }

    // create a typedef for the Graph type
    typedef boost::adjacency_list<boost::vecS, boost::vecS,
            boost::bidirectionalS, VertexData,
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
            // add edge and edge property
            // first is the edge. second is a bool telling you whether this is a new edge
            // or an existing one.
            auto e = boost::add_edge(map_vertex_descriptor[neighbour_id], map_vertex_descriptor[node_id], g).first;
            //cout<<"neighbour_id:"<<neighbour_id<<"  node_id"<<node_id<<"  source:"<<g[source(e,g)].node_id<<"  target:"<<g[target(e,g)].node_id<<endl;
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

    // find the shortest path
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
    std::map<int, std::list<int>> transform_trees;
    for (boost::tie(vi, vend) = boost::vertices(g); vi != vend; ++vi)
    {
        std::list<int> transform_tree;
        auto cur_vertex = *vi;
        if(vertex_property[cur_vertex] != vertex_property[map_vertex_descriptor[0]])
            // add self
            transform_tree.push_back(vertex_property[cur_vertex]);
        while(1)
        {
            if(vertex_property[cur_vertex] != vertex_property[map_vertex_descriptor[0]])
            {
                // add parent
                transform_tree.push_back(vertex_property[p[cur_vertex]]);
                cur_vertex = p[cur_vertex];
            }
            else
            {
                break;
            }
        }
        if(transform_tree.size())
            transform_trees[transform_tree.front()] = transform_tree;
    }

//    int debug_index = 0;
//    for(auto it = transform_trees.cbegin(); it != transform_trees.cend(); it++)
//    {
//        auto transform_tree = (*it).second;
//        for(auto it_t = transform_tree.cbegin(); it_t != transform_tree.cend(); it_t++)
//        {
//            cout<<*it_t<<"-";
//        }
//        cout<<"     "<<debug_index++<<endl;
//    }

    // start to merge pointcloud
    // read the first frame
    vtkSmartPointer<vtkGenericDataObjectReader> reader =
            vtkSmartPointer<vtkGenericDataObjectReader>::New();
    string inputFilename = "../frames/" +std::to_string(0)+"/DataPoints.vtk";
    reader->SetFileName(inputFilename.c_str());
    reader->Update();
    // add the first frame
    int insert_idx = 0;
    vtkSmartPointer<vtkPoints> merge_pts =
            vtkSmartPointer<vtkPoints>::New();
    merge_pts->InsertPoints(0,
                            reader->GetPolyDataOutput()->GetNumberOfPoints(),
                            insert_idx,
                            reader->GetPolyDataOutput()->GetPoints());

    auto p_propery_edge = boost::get(&EdgeData::param, g);
    for(auto it = transform_trees.cbegin(); it != transform_trees.cend(); it++)
    {
        vtkSmartPointer<vtkMatrix4x4> matrixA =
                vtkSmartPointer<vtkMatrix4x4>::New();
        vtkSmartPointer<vtkMatrix4x4> matrixB =
                vtkSmartPointer<vtkMatrix4x4>::New();
        vtkSmartPointer<vtkMatrix4x4> matrixC =
                vtkSmartPointer<vtkMatrix4x4>::New();
        matrixA->Identity();


        int node_id =  (*it).first;
        auto transform_tree = (*it).second;

        //if(node_id != 187)  continue;
        //cout<<"node_id:"<<node_id<<endl;
        // parse transform tree
        int pre_frame_id;
        int cur_frame_id = node_id;
        for(auto it_t = transform_tree.cbegin(); it_t != transform_tree.cend(); it_t++)
        {
            if(it_t == transform_tree.cbegin())
                continue;
            pre_frame_id = cur_frame_id;
            cur_frame_id = *it_t;
            auto e = edge(map_vertex_descriptor[pre_frame_id], map_vertex_descriptor[cur_frame_id],g).first;
            vtkQuaternion<double> quaternion_(p_propery_edge[e][3],
                    p_propery_edge[e][4],
                    p_propery_edge[e][5],
                    p_propery_edge[e][6]);
            double q[3][3];
            quaternion_.ToMatrix3x3(q);

            matrixB->SetElement(0,0,q[0][0]);
            matrixB->SetElement(0,1,q[0][1]);
            matrixB->SetElement(0,2,q[0][2]);
            matrixB->SetElement(0,3,p_propery_edge[e][0]);
            matrixB->SetElement(1,0,q[1][0]);
            matrixB->SetElement(1,1,q[1][1]);
            matrixB->SetElement(1,2,q[1][2]);
            matrixB->SetElement(1,3,p_propery_edge[e][1]);
            matrixB->SetElement(2,0,q[2][0]);
            matrixB->SetElement(2,1,q[2][1]);
            matrixB->SetElement(2,2,q[2][2]);
            matrixB->SetElement(2,3,p_propery_edge[e][2]);
            matrixB->SetElement(3,0,0);
            matrixB->SetElement(3,1,0);
            matrixB->SetElement(3,2,0);
            matrixB->SetElement(3,3,1);

            vtkMatrix4x4::Multiply4x4( matrixB,matrixA, matrixC);
            matrixA->DeepCopy(matrixC);
            //cout<<"A:\n";
           // matrixA->Print(std::cout);
        }/* for */
        inputFilename = "../frames/" +std::to_string(node_id)+"/DataPoints.vtk";
        reader->SetFileName(inputFilename.c_str());
        reader->Update();
        vtkSmartPointer<vtkPolyData> output = reader->GetPolyDataOutput();
        vtkSmartPointer<vtkPoints> inPts = output->GetPoints();

        vtkSmartPointer<vtkPerspectiveTransform> perspectiveTransform =
                vtkSmartPointer<vtkPerspectiveTransform>::New();
        perspectiveTransform->SetMatrix(matrixA);

        cout<<"node_id:"<<node_id<<"  A:\n";
        matrixA->Print(std::cout);

        perspectiveTransform->TransformPoints(inPts, merge_pts);
    }

    // 渲染合并后的点云
    vtkSmartPointer<vtkPoints> point_concat = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkPerspectiveTransform> perspectiveTransform_t =
            vtkSmartPointer<vtkPerspectiveTransform>::New();

    double bounds_pre[6];
    merge_pts->GetBounds(bounds_pre);
    perspectiveTransform_t->Translate(-(*bounds_pre),
                                      -(*(bounds_pre + 2)),
                                      -(*(bounds_pre+4)));
    perspectiveTransform_t->TransformPoints(merge_pts,//polydata_contat->GetPoints(),
                                            point_concat);
    double *bounds_dst = point_concat->GetBounds();

    //cout<<*(bounds_dst+1)<<"    "<<*(bounds_dst+3)<<endl;
    int rows, cols, intensity_max;
    rows = *(bounds_dst+1)/RESOLUTION + 0.5;
    cols = *(bounds_dst+3)/RESOLUTION + 0.5;
    intensity_max = *(bounds_dst+5)/RESOLUTION;
    cv::Mat image(rows,cols,CV_8UC1,cv::Scalar(0));

    //cout<<"rows:"<<rows<<"  cols:"<<cols<<endl;
    for(vtkIdType i = 0; i < point_concat->GetNumberOfPoints(); i++)
    {
        double *pts = point_concat->GetPoint(i);
        int x = *(pts)/RESOLUTION;
        int y = *(pts+1)/RESOLUTION;
        int intensity = 255*(*(pts+2)/RESOLUTION)/intensity_max;

        image.at<uchar>(x,y) = intensity>image.at<uchar>(x,y)?
                    intensity:image.at<uchar>(x,y);
        //cout<<"x:"<<x<<"  y:"<<y<<"  intensity:"<<intensity
        //<<"  image.at<uchar>(x,y):"<< int(image.at<uchar>(x,y))<<endl;
    }
    cv::Mat colormap;
    cv::applyColorMap(image,colormap,cv::COLORMAP_JET);

    //无点的地方设置为白色
    Mat copy_image;
    image.copyTo(copy_image);
    copy_image.setTo(255,copy_image);
    Mat full(copy_image.size(),copy_image.type(),cv::Scalar(255));
    Mat mask = full - copy_image;
    colormap.setTo(Vec3b(255,255,255),mask);

    cv::imwrite("intensity.bmp",colormap);
    cv::imwrite("image.bmp",image);
    imshow("intensity",colormap);


    // 渲染和保存点云文件
    vtkSmartPointer<vtkAppendPolyData> appendPolyData = vtkSmartPointer<vtkAppendPolyData>::New();
    appendPolyData->SetOutputPointsPrecision(vtkAlgorithm::SINGLE_PRECISION);
    vtkSmartPointer<vtkPolyData> polyData = vtkSmartPointer<vtkPolyData>::New();
    polyData->SetPoints(point_concat);

    appendPolyData->AddInputData(polyData);
    vtkSmartPointer<vtkPLYWriter> plyWriter =
            vtkSmartPointer<vtkPLYWriter>::New();
    plyWriter->SetFileName("merge.ply");
    plyWriter->SetFileTypeToASCII();
    plyWriter->SetInputConnection(appendPolyData->GetOutputPort());
    plyWriter->Write();


    // Visualize
//    vtkSmartPointer<vtkPolyDataMapper> mapper =
//            vtkSmartPointer<vtkPolyDataMapper>::New();
//    mapper->SetInputConnection(appendPolyData->GetOutputPort());

//    vtkSmartPointer<vtkActor> actor =
//            vtkSmartPointer<vtkActor>::New();
//    actor->SetMapper(mapper);

//    vtkSmartPointer<vtkRenderer> renderer =
//            vtkSmartPointer<vtkRenderer>::New();
//    vtkSmartPointer<vtkRenderWindow> renderWindow =
//            vtkSmartPointer<vtkRenderWindow>::New();
//    renderWindow->AddRenderer(renderer);
//    vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
//            vtkSmartPointer<vtkRenderWindowInteractor>::New();
//    renderWindowInteractor->SetRenderWindow(renderWindow);

//    renderer->AddActor(actor);
//    renderer->SetBackground(1,1,1); // Background color green

//    renderWindow->Render();

//    renderWindowInteractor->Start();

//    cv::waitKey(0);










    //        for(auto it = transform_trees.cbegin(); it != transform_trees.cend(); it++)
    //        {
    //            auto transform_tree = (*it).second;
    //             for(auto it_t = transform_tree.cbegin(); it_t != transform_tree.cend(); it_t++)
    //             {
    //                 cout<<*it_t<<"-";
    //             }
    //             cout<<endl;
    //        }

    // get matrix





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


//    std::ofstream dot_file("dijkstra-eg.dot");

//     dot_file << "digraph D {\n"
//       << "  rankdir=LR\n"
//       << "  size=\"4,3\"\n"
//       << "  ratio=\"fill\"\n"
//       << "  edge[style=\"bold\"]\n" << "  node[shape=\"circle\"]\n";

//     graph_traits < Graph >::edge_iterator ei, ei_end;
//     for (boost::tie(ei, ei_end) = edges(g); ei != ei_end; ++ei) {
//       graph_traits < Graph >::edge_descriptor e = *ei;
//       graph_traits < Graph >::vertex_descriptor
//         u = source(e, g), v = target(e, g);
//       dot_file << g[u].node_id << " -> " << g[v].node_id
//         << "[label=\"" << get(weightmap, e) << "\"";
//       if (p[v] == u)
//         dot_file << ", color=\"black\"";
//       else
//         dot_file << ", color=\"grey\"";
//       dot_file << "]";
//     }
//     dot_file << "}";

}

