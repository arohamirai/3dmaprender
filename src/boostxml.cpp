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
    string map_filepath = "../frames/map.xml";
    boost::property_tree::ptree pt;
    boost::property_tree::read_xml(map_filepath, pt);
    if (pt.empty()){
        printf("Load %s map.xml failed!", map_filepath.c_str());
        return false;
    }

    vtkSmartPointer<vtkGenericDataObjectReader> reader = vtkSmartPointer<vtkGenericDataObjectReader>::New();

    std::set<int> node_set;
    boost::property_tree::ptree p_node_list =
            pt.get_child("map.node_list");
    for(auto p_node_it = p_node_list.begin(); p_node_it != p_node_list.end(); p_node_it++){
        // 获取node_id
        int node_id = p_node_it->second.get<int>("id");
                cout<<node_id<<endl;
        reader->SetFileName(string("../frames/"+ to_string(node_id) + "/DataPoints.vtk").c_str());
        reader->Update();
        vtkSmartPointer<vtkPoints> inPts = reader->GetPolyDataOutput()->GetPoints();
        vtkSmartPointer<vtkPoints> outPts = vtkSmartPointer<vtkPoints>::New();
        vtkSmartPointer<vtkPerspectiveTransform> perspectiveTransform_t =
                vtkSmartPointer<vtkPerspectiveTransform>::New();


        double bounds_pre[6],bounds_aft[6];
        inPts->GetBounds(bounds_pre);
        perspectiveTransform_t->Translate(-bounds_pre[0],
                -bounds_pre[2],
                -bounds_pre[4]);
        perspectiveTransform_t->TransformPoints(inPts, outPts);
        outPts->GetBounds(bounds_aft);

        int rows, cols, intensity_max;
        rows = bounds_aft[1]/RESOLUTION + 1;
        cols = bounds_aft[3]/RESOLUTION + 1;
        intensity_max = bounds_aft[5]/RESOLUTION;
        cv::Mat image(rows,cols,CV_8UC1,cv::Scalar(0));
        for(vtkIdType i = 0; i < outPts->GetNumberOfPoints(); i++)
        {
            double pts[3];
            outPts->GetPoint(i,pts);
            int x = pts[0]/RESOLUTION;
            int y = pts[1]/RESOLUTION;
            int intensity = 255*(pts[2]/RESOLUTION)/intensity_max;

            if(intensity > 128) //丢弃
                continue;

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

        cv::imwrite("./result/intensity"+to_string(node_id)+".bmp",colormap);
        //cv::imwrite("image.bmp",image);
        imshow("intensity",colormap);
        cv::waitKey(33);
        //inPts->Delete();
        //outPts->Delete();

    }
    return EXIT_SUCCESS;
}

