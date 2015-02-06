/*
 * Viewer_instructor.cpp
 *
 *  Created on: Aug 26, 2014
 *      Author: dsalinas
 */

#include <utility>
#include "Viewer_instructor.h"
#include "utils/UI_utils.h"
#include "FirstCoordProjector.h"
#include "Projection_from_points.h"


Viewer_instructor::Viewer_instructor(QWidget* parent,
		Viewer* viewer,
		const Complex& mesh
):viewer_(viewer),mesh_(mesh),projector_(new FirstCoordProjector3D<Complex>(mesh)){
	viewer_->set_instructor(this);
}


void
Viewer_instructor::initialize_bounding_box(){
	auto pair_bounding_box = compute_bounding_box_corners();
	std::cout<<"bounding box:\n";
	std::cout<<pair_bounding_box.first<<"  --  ";
	std::cout<<pair_bounding_box.second<<"\n";

	viewer_->set_bounding_box(pair_bounding_box.first,pair_bounding_box.second);
	viewer_->init_scene();
}

std::pair<Viewer_instructor::Point_3,Viewer_instructor::Point_3>
Viewer_instructor::compute_bounding_box_corners(){
	if(mesh_.empty()){
		return std::make_pair(Point_3(-1,-1,-1),Point_3(1,1,1));
	}
	else{
		double x_min = 1e10;
		double y_min = 1e10;
		double z_min = 1e10;
		double x_max = -1e10;
		double y_max = -1e10;
		double z_max = -1e10;
		for( auto vi : mesh_.vertex_range())
		{
			auto pt = proj(vi);
			x_min = (std::min)(x_min,pt.x());
			y_min = (std::min)(y_min,pt.y());
			z_min = (std::min)(z_min,pt.z());

			x_max = (std::max)(x_max,pt.x());
			y_max = (std::max)(y_max,pt.y());
			z_max = (std::max)(z_max,pt.z());

		}
		if(z_min==0&&z_max==0) z_max=1;
		return std::make_pair(
				Point_3(x_min,y_min,z_min),
				Point_3(x_max,y_max,z_max)
		);
	}
}

void
Viewer_instructor::show_entire_scene(){
	viewer_->show_entire_scene();
}

const qglviewer::Camera*
Viewer_instructor::camera() const{
	return viewer_->camera();
}

int
Viewer_instructor::width() const{
	return viewer_->width();
}
int
Viewer_instructor::height() const{
	return viewer_->height();
}

/**
 * to change display parameters
 */
View_parameter&
Viewer_instructor::view_params(){
	return view_params_;
}


void
Viewer_instructor::give_instructions(){
	if(view_params_.relative_light)
		viewer_->set_light_direction();
	else
		viewer_->set_light_direction(view_params_.theta,view_params_.phi);
	viewer_->set_light();

	if (view_params_.edge_mode) draw_edges();
	if (view_params_.triangle_mode)	draw_triangles();
	if (view_params_.vertex_mode) draw_points();

}

void
Viewer_instructor::draw_edges(){
	viewer_->begin_draw_edges(view_params_.size_edges,false);

	for(auto edge : mesh_.edge_range()){
		set_color_edge(edge);
//		const Point& a = mesh_.point(mesh_.first_vertex(edge));
//		const Point& b = mesh_.point(mesh_.second_vertex(edge)) ;
		viewer_->draw_edges(proj(mesh_.first_vertex(edge)),proj(mesh_.second_vertex(edge)));
	}

	viewer_->end_draw_edges();
}

void
Viewer_instructor::draw_triangles(){
	const double size_triangles = 1.0;
	viewer_->begin_draw_triangles(size_triangles,view_params_.light);

	for(const auto& fit : mesh_.triangle_range()) {
		set_color_triangle(fit);
		if(view_params_.triangle_mode){
			auto fit_it  = fit.begin();
			const Point_3& p1 = proj(*fit_it);
			const Point_3& p2 = proj(*(++fit_it));
			const Point_3& p3 = proj(*(++fit_it));
			viewer_->draw_triangles(p1,p2,p3);
		}
	}
	viewer_->end_draw_triangles();
}

void
Viewer_instructor::draw_points(){
	viewer_->begin_draw_points(	view_params_.size_vertices);
	for( auto vi : mesh_.vertex_range())
	{
		viewer_->set_size_point(view_params_.size_vertices);
		set_color_vertex(vi);
		viewer_->draw_points(proj(vi));
	}
	viewer_->end_draw_points();
}


void
Viewer_instructor::draw_edge(const Point&,const Point&){

}

void
Viewer_instructor::draw_point(const Point&){

}


/**
 * set the right color of vertex/edge/triangle considering the view_params choice
 */
void
Viewer_instructor::set_color_vertex(Vertex_handle vh){
	viewer_->set_color(Color(view_params_.light_edges,view_params_.light_edges,view_params_.light_edges));
}

void
Viewer_instructor::set_color_edge(Edge_handle eh)	{
	viewer_->set_color(Color(view_params_.light_edges,view_params_.light_edges,view_params_.light_edges));
}

void
Viewer_instructor::set_color_triangle(const Simplex_handle& triangle){
	viewer_->set_color(Color(view_params_.light_triangles,view_params_.light_triangles,view_params_.light_triangles));
}


void
Viewer_instructor::set_projection_from_first_three_coordinates(){
	projector_.reset(new FirstCoordProjector3D<Complex>(mesh_));
}

void
Viewer_instructor::set_projection_to_specific_points(std::vector<Point_3> points){
	std::cout <<"set new proj\n";
	projector_.reset(new Projection_from_points(points));
}



void
Viewer_instructor::change_projector(Projector3D* new_projector){
	projector_.reset(new_projector);
}


Viewer_instructor::Point_3
Viewer_instructor::proj(Vertex_handle v) const{
	return (*projector_)(v);
}

void
Viewer_instructor::sceneChanged(){
	UIDBG("sceneChanged");
	viewer_->update_GL();
}

void
Viewer_instructor::change_draw_vertices(){
	view_params_.change_vertex_mode();
}
void
Viewer_instructor::change_draw_edges(){
	view_params_.change_edge_mode();
}
void
Viewer_instructor::change_draw_triangles(){
	view_params_.change_triangle_mode();
}
void
Viewer_instructor::change_light(){
	view_params_.light =! view_params_.light ;
}

#include "Viewer_instructor.moc"
