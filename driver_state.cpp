#include "driver_state.h"
#include <cstring>

driver_state::driver_state()
{
}

driver_state::~driver_state()
{
    delete [] image_color;
    delete [] image_depth;
}

// This function should allocate and initialize the arrays that store color and
// depth.  This is not done during the constructor since the width and height
// are not known when this class is constructed.
void initialize_render(driver_state& state, int width, int height)
{
    state.image_width=width;
    state.image_height=height;
    state.image_color=0;
    state.image_depth=0;

    state.image_color = new pixel[width * height];
    for (int i = 0; i < width * height; ++i) {
	state.image_color[i] = make_pixel(0, 0, 0);
    }
}

// This function will be called to render the data that has been stored in this class.
// Valid values of type are:
//   render_type::triangle - Each group of three vertices corresponds to a triangle.
//   render_type::indexed -  Each group of three indices in index_data corresponds
//                           to a triangle.  These numbers are indices into vertex_data.
//   render_type::fan -      The vertices are to be interpreted as a triangle fan.
//   render_type::strip -    The vertices are to be interpreted as a triangle strip.
void render(driver_state& state, render_type type)
{
    switch(type) {
	case render_type::triangle:
		for (int i = 0; i < state.num_vertices; i += 3) { // ith triangle
		    data_geometry** triangle = new data_geometry*[3];
		    for (int j = 0; j < 3; ++j) { // jth vertex
			triangle[j] = new data_geometry;
			triangle[j]->data = new float[MAX_FLOATS_PER_VERTEX];
			for (int k = 0; k < state.floats_per_vertex; ++k) {// kth data
			    triangle[j]->data[k] = state.vertex_data[k + state.floats_per_vertex*(i+j)];
			}
		    }
		    rasterize_triangle(state, (const data_geometry**)triangle);
		}
		break;
	case render_type::indexed: break;
	case render_type::fan: break;
	case render_type::strip: break;
	default: break;
    }
    std::cout<<"TODO: implement rendering for indexed, fan, and strip."<<std::endl;
}


// This function clips a triangle (defined by the three vertices in the "in" array).
// It will be called recursively, once for each clipping face (face=0, 1, ..., 5) to
// clip against each of the clipping faces in turn.  When face=6, clip_triangle should
// simply pass the call on to rasterize_triangle.
void clip_triangle(driver_state& state, const data_geometry* in[3],int face)
{
    if(face==6)
    {
        rasterize_triangle(state, in);
        return;
    }
    std::cout<<"TODO: implement clipping. (The current code passes the triangle through without clipping them.)"<<std::endl;
    clip_triangle(state,in,face+1);
}

// Rasterize the triangle defined by the three vertices in the "in" array.  This
// function is responsible for rasterization, interpolation of data to
// fragments, calling the fragment shader, and z-buffering.
void rasterize_triangle(driver_state& state, const data_geometry* in[3])
{
    data_geometry* out = new data_geometry[3];
    int w = state.image_width;
    int h = state.image_height;
    int pixCor[3][2];
    for (int index = 0; index < 3; ++index) {
	data_vertex v;
	v.data = in[index]->data;
	state.vertex_shader((const data_vertex)v, out[index], state.uniform_data);

	out[index].gl_Position /= out[index].gl_Position[3];
	auto x = out[index].gl_Position[0];
	auto y = out[index].gl_Position[1];
	pixCor[index][0] = (w/2.0)*x + w/2.0 - 0.5;
	pixCor[index][1] = (h/2.0)*y + h/2.0 - 0.5;
	int image_index = pixCor[index][0] + pixCor[index][1] * w;
	state.image_color[image_index] = make_pixel(255,255,255);
    }
    /*
    int ax = (w/2.0)*out[0].gl_Position[0] + (w/2.0) - 0.5;
    int ay = (h/2.0)*out[0].gl_Position[1] + (w/2.0) - 0.5;
    int bx = (w/2.0)*out[1].gl_Position[0] + (w/2.0) - 0.5;
    int by = (h/2.0)*out[1].gl_Position[1] + (w/2.0) - 0.5;
    int cx = (w/2.0)*out[2].gl_Position[0] + (w/2.0) - 0.5;
    int cy = (h/2.0)*out[2].gl_Position[1] + (w/2.0) - 0.5;
    */
    
    int ax = pixCor[0][0]; int ay = pixCor[0][1];
    int bx = pixCor[1][0]; int by = pixCor[1][1];
    int cx = pixCor[2][0]; int cy = pixCor[2][1];
    std::cout << ax << ay << bx << by << cx << cy << std::endl;
    double abc = 0.5 * ((bx*cy - cx*by) - (ax*cy - cx*ay) - (ax*by - bx*ay));

    for (int py = 0; py < h; ++py) {
    for (int px = 0; px < w; ++px) {
	double pbc = 0.5 * (px * (by - cy) + bx * (cy - py) + cx * (py - by));
	double apc = 0.5 * (ax * (py - cy) + px * (cy - ay) + cx * (ay - by));
	double abp = 0.5 * (ax * (by - py) + bx * (py - ay) + px * (ay - by));
	double alpha = pbc/abc;
	double beta = apc/abc;
	double gamma = abp/abc;
	if (alpha >= 0 && beta >= 0 && gamma >= 0) {
	    int index = px + py * w;
	    state.image_color[index] = make_pixel(255,255,255);
	}
    }
    }
    
    std::cout<<"TODO: implement rasterization"<<std::endl;
}

