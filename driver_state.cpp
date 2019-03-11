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
    state.image_depth = new float[width * height];
    for (int i = 0; i < width * height; ++i) {
	state.image_color[i] = make_pixel(0, 0, 0);
	state.image_depth[i] = 1;
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
		{
		for (int i = 0; i < state.num_vertices; i += 3) { // ith triangle
		    data_geometry** triangle = new data_geometry*[3];
		    for (int j = 0; j < 3; ++j) { // jth vertex
			triangle[j] = new data_geometry;
			data_vertex v;
			v.data = new float[MAX_FLOATS_PER_VERTEX];
			triangle[j]->data = new float[MAX_FLOATS_PER_VERTEX];
			for (int k = 0; k < state.floats_per_vertex; ++k) {// kth data
			    v.data[k] = state.vertex_data[k + state.floats_per_vertex*(i+j)];
			    triangle[j]->data[k] = v.data[k];
			}
			state.vertex_shader((const data_vertex)v, *triangle[j], state.uniform_data);
		    }

		    //triangle[0]->gl_Position /= triangle[0]->gl_Position[3];
		    //triangle[1]->gl_Position /= triangle[1]->gl_Position[3];
		    //triangle[2]->gl_Position /= triangle[2]->gl_Position[3];

		    //rasterize_triangle(state, (const data_geometry**)triangle);
		    clip_triangle(state, (const data_geometry**)triangle, 0);
		}
		}
		break;
	case render_type::indexed: 
	{
	    const data_geometry *out[3];
	    data_geometry triangle[3];
	    data_vertex vertices[3];

	    for (int i = 0; i < 3*state.num_triangles; i += 3) {
		for(int j = 0; j < 3; j++) {
		   vertices[j].data = &state.vertex_data[state.index_data[i+j] * state.floats_per_vertex];
		   triangle[j].data = vertices[j].data;
		   state.vertex_shader(vertices[j], triangle[j], state.uniform_data);
		   out[j] = &triangle[j];
		}
		//rasterize_triangle(state, out);
		clip_triangle(state, out, 0);
	    }
	}
	break;
	case render_type::fan: 
	{
	    const data_geometry *out[3];
	    data_geometry triangle[3];
	    data_vertex vertices[3];

	    for (int i = 0; i < state.num_vertices; ++i) {
		for (int j = 0; j < 3; ++j) {
		    int k = i + j;
		    if (j == 0) { k = 0; }
		    vertices[j].data = &state.vertex_data[k*state.floats_per_vertex];
		    triangle[j].data = vertices[j].data;
		    state.vertex_shader(vertices[j], triangle[j], state.uniform_data);
		    out[j] = &triangle[j];
		}
		//rasterize_triangle(state, out);
		clip_triangle(state, out, 0);
	    }
	}
	break;
	case render_type::strip: 
	{
	   const data_geometry *out[3];
	   data_geometry triangle[3];
	   data_vertex vertices[3];

	   for (int i = 0; i < state.num_vertices-2; ++i) {
		for (int j = 0; j < 3; ++j) {
		    vertices[j].data = &state.vertex_data[(i+j)*state.floats_per_vertex];
		    triangle[j].data = vertices[j].data;
		    state.vertex_shader(vertices[j], triangle[j], state.uniform_data);
		    out[j] = &triangle[j];
		}
		//rasterize_triangle(state, out);
		clip_triangle(state, out, 0);
	   }
	}
	break;
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
    //std::cout << "clip\n";
    if(face == 1)
    {
        rasterize_triangle(state, in);
        return;
    }
    //std::cout<<"TODO: implement clipping. (The current code passes the triangle through without clipping them.)"<<std::endl;
    vec4 a = in[0]->gl_Position;
    vec4 b = in[1]->gl_Position;
    vec4 c = in[2]->gl_Position;
    const data_geometry *in2[3] = {in[0], in[1], in[2]};
    data_geometry data1[3];
    data_geometry data2[3];
    float a1, a2, b1, b2;
    vec4 p1, p2;

    if (a[2] < -a[3] && b[2] < -b[3] && c[2] < -c[3]) {
	//std::cout << "not worth\n";
	return;
    } else {
	if (a[2] < -a[3] && b[2] >= -b[3] && c[2] >= -c[3]) {
	    b1 = (-b[3] - b[2]) / (a[2] + a[3] - b[3] - b[2]);
	    b2 = (-a[3] - a[2]) / (c[2] + c[3] - a[3] - a[2]);
	    p1 = b1 * a + (1 - b1) * b;
	    p2 = b2 * c + (1 - b2) * a;

	    data1[0].data = new float[MAX_FLOATS_PER_VERTEX];
	    data1[1] = *in[1];
	    data1[2] = *in[2];

	    for (int i = 0; i < state.floats_per_vertex; ++i) {
		switch (state.interp_rules[i]) {
		    case interp_type::flat:
			data1[0].data[i] = in[0]->data[i];
			break;
		    case interp_type::smooth:
			data1[0].data[i] = b2 * in[2]->data[i] + (1 - b2) * in[0]->data[i];
			break;
		    case interp_type::noperspective:
			a1 = b2 * in[2]->gl_Position[3] / (b2 * in[2]->gl_Position[3] + (1 - b2) * in[0]->gl_Position[3]);
			data1[0].data[i] = a1 * in[2]->data[i] + (1 - a1) * in[0]->data[i];
			break;
		    default:
			break;
		}
	    }
	    data1[0].gl_Position = p2;
	    in2[0] = &data1[0];
	    in2[1] = &data1[1];
	    in2[2] = &data1[2];

	    clip_triangle(state,in2,face+1);

	    data2[0].data = new float[MAX_FLOATS_PER_VERTEX];
	    data2[2] = *in[2];

	    for (int i = 0; i < state.floats_per_vertex; ++i) {
		switch (state.interp_rules[i]) {
		    case interp_type::flat:
			data2[0].data[i] = in[0]->data[i];
			break;
		    case interp_type::smooth:
			data2[0].data[i] = b1 * in[0]->data[i] + (1 - b1) * in[1]->data[i];
			break;
		    case interp_type::noperspective:
			a1 = b1 * in[0]->gl_Position[3] / (b1 * in[0]->gl_Position[3] + (1 - b1) * in[1]->gl_Position[3]);
			data2[0].data[i] = a1 * in[0]->data[i] + (1 - a1) * in[1]->data[i];
			break;
		    default:
			break;
		}
	    }
	    
	    data2[0].gl_Position = p1;
	    in2[0] = &data2[0];
	    in2[1] = &data1[1];
	    in2[2] = &data1[0];
	}
    }

    clip_triangle(state,in2,face+1);
}

// Rasterize the triangle defined by the three vertices in the "in" array.  This
// function is responsible for rasterization, interpolation of data to
// fragments, calling the fragment shader, and z-buffering.
void rasterize_triangle(driver_state& state, const data_geometry* in[3])
{
    //std::cout << "rasterize\n";
    int w = state.image_width;
    int h = state.image_height;
    //int pixCor[3][2];
    /*
    for (int index = 0; index < 3; ++index) {
	auto x = in[index]->gl_Position[0]/in[index]->gl_Position[3];
	auto y = in[index]->gl_Position[1]/in[index]->gl_Position[3];
	pixCor[index][0] = (w/2.0)*x + w/2.0 - 0.5;
	pixCor[index][1] = (h/2.0)*y + h/2.0 - 0.5;
	//int image_index = pixCor[index][0] + pixCor[index][1] * w;
    }
    */
    //float ax = pixCor[0][0]; float ay = pixCor[0][1];
    //float bx = pixCor[1][0]; float by = pixCor[1][1];
    //float cx = pixCor[2][0]; float cy = pixCor[2][1];
   
    float ax = (w/2.0)*(in[0]->gl_Position[0] / in[0]->gl_Position[3]) + (w/2.0) - (0.5);
    float ay = (h/2.0)*(in[0]->gl_Position[1] / in[0]->gl_Position[3]) + (h/2.0) - (0.5);
    float bx = (w/2.0)*(in[1]->gl_Position[0] / in[1]->gl_Position[3]) + (w/2.0) - (0.5);
    float by = (h/2.0)*(in[1]->gl_Position[1] / in[1]->gl_Position[3]) + (h/2.0) - (0.5);
    float cx = (w/2.0)*(in[2]->gl_Position[0] / in[2]->gl_Position[3]) + (w/2.0) - (0.5);
    float cy = (h/2.0)*(in[2]->gl_Position[1] / in[2]->gl_Position[3]) + (h/2.0) - (0.5);

    float minX = std::min(ax, std::min(bx, cx));
    float maxX = std::max(ax, std::max(bx, cx));
    float minY = std::min(ay, std::min(by, cy));
    float maxY = std::max(ay, std::max(by, cy));

    if (minX < 0) { minX = 0; }
    if (minY < 0) { minY = 0; }
    if (maxX > w) { maxX = w; }
    if (maxY > h) { maxY = h; }

    double abc = 0.5 * ((bx*cy - cx*by) - (ax*cy - cx*ay) + (ax*by - bx*ay));

    for (int py = minY; py <= maxY; ++py) {
    for (int px = minX; px <= maxX; ++px) {

	double pbc = 0.5 * ((bx*cy - cx*by) + (by - cy)*px + (cx - bx)*py);
	double apc = 0.5 * ((cx*ay - ax*cy) + (cy - ay)*px + (ax - cx)*py);
	double abp = 0.5 * ((ax*by - bx*ay) + (ay - by)*px + (bx - ax)*py);
	double alpha = pbc/abc;
	double beta = apc/abc;
	double gamma = abp/abc;

	if (alpha >= 0 && beta >= 0 && gamma >= 0) {
	    int index = px + py * w;
	    auto *data = new float[MAX_FLOATS_PER_VERTEX];
	    data_fragment f{data};
	    data_output o;
	    float depth = alpha * in[0]->gl_Position[2]/in[0]->gl_Position[3] + beta * in[1]->gl_Position[2]/in[1]->gl_Position[3] + gamma * in[2]->gl_Position[2]/in[2]->gl_Position[3];

	    if (depth > state.image_depth[index]) {
		continue;
	    }
	    
	    for (int i = 0; i < state.floats_per_vertex; ++i) {
		float k_gour;
		float a_persp, b_persp, g_persp;
		switch(state.interp_rules[i]) {
		    case interp_type::flat:
			f.data[i] = in[0]->data[i];
			break;
		    case interp_type::smooth:
			k_gour = (alpha/in[0]->gl_Position[3] + beta/in[1]->gl_Position[3] + gamma/in[2]->gl_Position[3]);
			a_persp = alpha / k_gour / in[0]->gl_Position[3];
			b_persp = beta / k_gour / in[1]->gl_Position[3];
			g_persp = gamma / k_gour / in[2]->gl_Position[3];
			f.data[i] = a_persp*in[0]->data[i] + b_persp*in[1]->data[i] + g_persp*in[2]->data[i];
			break;
		    case interp_type::noperspective:
			f.data[i] = alpha*in[0]->data[i] + beta*in[1]->data[i] + gamma*in[2]->data[i];
			break;
		    default:
			break;
		}
	    }
	    state.fragment_shader((const data_fragment)f, o, state.uniform_data);
	    o.output_color *= 255;
	    state.image_color[index] = make_pixel(o.output_color[0],o.output_color[1],o.output_color[2]);
	    state.image_depth[index] = depth;
	}
    }
    }
    
}

