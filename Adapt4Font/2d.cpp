////////////////////////////////////////////////////////////////////////////////////////////
// main1.cpp :
//

#include <stdio.h>
#include <dlfcn.h>

#include <iostream>
#include <fstream>

#include <base/exceptions.h>
#include <lac/vector.h>
#include <lac/sparsity_pattern.h>
#include <lac/sparse_matrix.h>

#include <AFEPack/Geometry.h>
#include <AFEPack/HGeometry.h>
#include <cmath>
#include <algorithm>
#define DIM 2
/*
double point_to_line_segment_distance(double point[2], double segment[2][2]) {
  // Extract the coordinates of the point and line segment
  double Ax = point[0];
  double Ay = point[1];
  double Bx = segment[0][0];
  double By = segment[0][1];
  double Cx = segment[1][0];
  double Cy = segment[1][1];
  
  // Calculate the length of the line segment
  double segment_length = sqrt(pow(Bx - Cx, 2) + pow(By - Cy, 2));
  
  // Calculate the dot products
  double dot_product_1 = (Ax - Bx)*(Cx - Bx) + (Ay - By)*(Cy - By);
  double dot_product_2 = (Ax - Cx)*(Bx - Cx) + (Ay - Cy)*(By - Cy);
  
  // Check if projection of point falls within line segment
  if (dot_product_1 >= 0 && dot_product_2 >= 0) {
    // Calculate distance using formula for point-to-line distance
    double distance = abs((Ax - Bx)*(Cy - By) - (Cx - Bx)*(Ay - By)) / sqrt(pow(Bx - Cx, 2) + pow(By - Cy, 2));
    return distance;
  } else {
    // Calculate distance to closest endpoint of line segment
    double distance_to_B = sqrt(pow(Ax - Bx, 2) + pow(Ay - By, 2));
    double distance_to_C = sqrt(pow(Ax - Cx, 2) + pow(Ay - Cy, 2));
    double distance = std::min(distance_to_B, distance_to_C);
    return distance;
  }
}*/

double point_to_line_segment_distance(double point[2], double segment[2][2]) {
  // Extract the coordinates of the point and line segment
  double Ax = point[0];
  double Ay = point[1];
  double Bx = segment[0][0];
  double By = segment[0][1];
  double Cx = segment[1][0];
  double Cy = segment[1][1];

  // Calculate the length of the line segment
  double segment_length = sqrt(pow(Bx - Cx, 2) + pow(By - Cy, 2));

  // Calculate the dot products
  double dot_product_1 = (Ax - Bx)*(Cx - Bx) + (Ay - By)*(Cy - By);
  double dot_product_2 = (Ax - Cx)*(Bx - Cx) + (Ay - Cy)*(By - Cy);

  // Check if projection of point falls within line segment
  if (dot_product_1 >= 0 && dot_product_2 >= 0) {
    // Calculate distance using formula for point-to-line distance
    double distance = fabs((Ax - Bx)*(Cy - By) - (Cx - Bx)*(Ay - By)) / sqrt(pow(Bx - Cx, 2) + pow(By - Cy, 2));
    //std::cout << "Projection point falls within line segment..." << std::endl;
    return distance;
  } else {
    // Calculate distance to closest endpoint of line segment
    double distance_to_B = sqrt(pow(Ax - Bx, 2) + pow(Ay - By, 2));
    double distance_to_C = sqrt(pow(Ax - Cx, 2) + pow(Ay - Cy, 2));
    double distance = std::min(distance_to_B, distance_to_C);
    //std::cout << "Projection point falls WITHOUT line segment..." << std::endl;
    return distance;
  }
}


int main(int argc, char * argv[])
{
	HGeometryTree<DIM> h_tree;
	h_tree.readEasyMesh("D");
	//h_tree.readEasyMesh("Coarse");
	IrregularMesh<DIM> irregular_mesh(h_tree);
	irregular_mesh.globalRefine(3);

	do {
		irregular_mesh.semiregularize();
		irregular_mesh.regularize(false);
		RegularMesh<DIM>& regular_mesh = irregular_mesh.regularMesh();
		regular_mesh.writeOpenDXData("D.dx");
		//regular_mesh.writeOpenDXData("Coarse.dx");
		std::cout << "Press ENTER to continue or CTRL+C to stop ..." << std::flush;
		getchar();

		Indicator<DIM> indicator(regular_mesh);
		//AFEPack::Point<DIM> c0(0.495, 0.5);
		//AFEPack::Point<DIM> c1(0.505, 0.5);
	        double Jie1[2][2] = {{0.1, 0.9}, {0.4, 0.9}};
	        double Jie2[2][2] = {{0.2, 0.86}, {0.2, 0.94}};
	        double Jie3[2][2] = {{0.3, 0.86}, {0.3, 0.94}};
	        double Jie4[2][2] = {{0.15, 0.82}, {0.35, 0.82}};
	        double Jie5[2][2] = {{0.25, 0.6}, {0.26, 0.82}};
	        double Jie6[2][2] = {{0.35, 0.65}, {0.35, 0.82}};
	        double Jie7[2][2] = {{0.35, 0.65}, {0.33, 0.67}};

		double Ri1[2][2] = {{0.65,0.9},{0.85,0.9}};
		double Ri2[2][2] = {{0.65,0.6},{0.85,0.6}};
		double Ri3[2][2] = {{0.65,0.75},{0.8,0.75}};
		double Ri4[2][2] = {{0.65,0.9},{0.65,0.6}};
		double Ri5[2][2] = {{0.85,0.9},{0.85,0.6}};

		double Kuai1[2][2] = {{0.06,0.2},{0.12,0.3}};
		double Kuai2[2][2] = {{0.15,0.1},{0.15,0.4}};
		double Kuai3[2][2] = {{0.16,0.29},{0.19,0.25}};
		double Kuai4[2][2] = {{0.25,0.31},{0.37,0.31}};
		double Kuai5[2][2] = {{0.37,0.31},{0.35,0.23}};
		double Kuai6[2][2] = {{0.22,0.23},{0.42,0.23}};
		double Kuai7[2][2] = {{0.30,0.38},{0.30,0.23}};
		double Kuai8[2][2] = {{0.30,0.23},{0.23,0.12}};
		double Kuai9[2][2] = {{0.30,0.23},{0.38,0.11}};

		double Le1[2][2] = {{0.67,0.4},{0.8,0.45}};
		double Le2[2][2] = {{0.67,0.4},{0.65,0.3}};
		double Le3[2][2] = {{0.65,0.3},{0.85,0.3}};
		double Le4[2][2] = {{0.75,0.425},{0.75,0.1}};
		double Le5[2][2] = {{0.75,0.1},{0.71,0.14}};
		double Le6[2][2] = {{0.64,0.16},{0.68,0.25}};
		double Le7[2][2] = {{0.80,0.25},{0.85,0.17}};
		
		for (int i = 0;i < regular_mesh.n_geometry(2);i ++) {
		    AFEPack::Point<DIM>& p0 = regular_mesh.point(regular_mesh.geometry(2,i).vertex(0));
		    AFEPack::Point<DIM>& p1 = regular_mesh.point(regular_mesh.geometry(2,i).vertex(1));
		    AFEPack::Point<DIM>& p2 = regular_mesh.point(regular_mesh.geometry(2,i).vertex(2));
		    AFEPack::Point<DIM> p((p0[0] + p1[0] + p2[0])/3., (p0[1] + p1[1] + p2[1])/3.);
		    double area = (p1[0] - p0[0])*(p2[1] - p0[1])
                                - (p2[0] - p0[0])*(p1[1] - p0[1]);
                    double point[2]={p[0],p[1]};
	            //std::cout<< p<<" Dist:" << point_to_line_segment_distance(point, segment)<<std::endl;	    
		    double min_value = 1.0;
		    if(p[0] <=0.5 && 0.5<=p[1])
		    {
		//	    std::cout<< p << " is in the district 1 ... " << std::endl;
		      double dist1 = point_to_line_segment_distance(point, Jie1);
		      double dist2 = point_to_line_segment_distance(point, Jie2);
		      double dist3 = point_to_line_segment_distance(point, Jie3);
		      double dist4 = point_to_line_segment_distance(point, Jie4);
		      double dist5 = point_to_line_segment_distance(point, Jie5);
		      double dist6 = point_to_line_segment_distance(point, Jie6);
		      double dist7 = point_to_line_segment_distance(point, Jie7);
		      double values[7] = {dist1,dist2,dist3,dist4,dist5,dist6,dist7};
		    min_value = *std::min_element(values, values + 7);
		    }
		    else if(0.5 < p[0] && 0.5<=p[1])
		    {
		//	    std::cout<< p << " is in the district 2 ... " << std::endl;
		      double dist1 = point_to_line_segment_distance(point, Ri1);
		      double dist2 = point_to_line_segment_distance(point, Ri2);
		      double dist3 = point_to_line_segment_distance(point, Ri3);
		      double dist4 = point_to_line_segment_distance(point, Ri4);
		      double dist5 = point_to_line_segment_distance(point, Ri5);
		      double values[5] = {dist1,dist2,dist3,dist4,dist5};
		      min_value = *std::min_element(values, values + 5);
		    }
		    else if( p[0] <= 0.5 && p[1] <= 0.5)
		    {
		//	    std::cout<< p << " is in the district 3 ... " << std::endl;
		      double dist1 = point_to_line_segment_distance(point, Kuai1);
		      double dist2 = point_to_line_segment_distance(point, Kuai2);
		      double dist3 = point_to_line_segment_distance(point, Kuai3);
		      double dist4 = point_to_line_segment_distance(point, Kuai4);
		      double dist5 = point_to_line_segment_distance(point, Kuai5);
		      double dist6 = point_to_line_segment_distance(point, Kuai6);
		      double dist7 = point_to_line_segment_distance(point, Kuai7);
		      double dist8 = point_to_line_segment_distance(point, Kuai8);
		      double dist9 = point_to_line_segment_distance(point, Kuai9);
		      double values[9] = {dist1,dist2,dist3,dist4,dist5,dist6,dist7,dist8,dist9};//,dist4,dist5};
		      min_value = *std::min_element(values, values + 9);
		    }
		    else
		    {
		      double dist1 = point_to_line_segment_distance(point, Le1);
		      double dist2 = point_to_line_segment_distance(point, Le2);
		      double dist3 = point_to_line_segment_distance(point, Le3);
		      double dist4 = point_to_line_segment_distance(point, Le4);
		      double dist5 = point_to_line_segment_distance(point, Le5);
		      double dist6 = point_to_line_segment_distance(point, Le6);
		      double dist7 = point_to_line_segment_distance(point, Le7);
		      double values[7] = {dist1,dist2,dist3,dist4,dist5,dist6,dist7};//,dist4,dist5};
		      min_value = *std::min_element(values, values + 7);
		    }
		   // std::cout << "dist: "<<dist << std::endl;
			if (min_value < 0.01)
			{
			  //std::cout<<"Point is: "<< p << "dist is : "<< dist << std::endl;
			  //std::cout<<"Point is: "<< point[0] << " "<<point[1] << "dist is : "<< dist << std::endl;
			  indicator[i] = area;
			}
		};

		MeshAdaptor<DIM> mesh_adaptor(irregular_mesh);
		mesh_adaptor.convergenceOrder() = 0.;
		mesh_adaptor.refineStep() = 0;
		mesh_adaptor.setIndicator(indicator);
		mesh_adaptor.tolerence() = 2.5e-10;
		mesh_adaptor.adapt();
	} while (1);	
};

#undef DIM

//
// end of file
////////////////////////////////////////////////////////////////////////////////////////////
